#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using NCTSSoS
using ManiDSDP

const MOI = NCTSSoS.MOI

const CHAIN_N = parse(Int, get(ENV, "HEISENBERG_N", "16"))
const RELAX_ORDER = parse(Int, get(ENV, "HEISENBERG_ORDER", "2"))
const MANIDSDP_RANK = parse(Int, get(ENV, "MANIDSDP_RANK", "600"))
const MANIDSDP_MAX_ITER = parse(Int, get(ENV, "MANIDSDP_MAX_ITER", "6"))
const MANIDSDP_TOL = parse(Float64, get(ENV, "MANIDSDP_TOL", "1e-4"))
const MANIDSDP_SIGMA = parse(Float64, get(ENV, "MANIDSDP_SIGMA", "0.12"))
const MANIDSDP_INNER_MAX_ITER = parse(Int, get(ENV, "MANIDSDP_INNER_MAX_ITER", "100"))
const MANIDSDP_INNER_TOL = parse(Float64, get(ENV, "MANIDSDP_INNER_TOL", "1e-5"))
const ED_MAXITER = parse(Int, get(ENV, "ED_MAXITER", "300"))
const ED_TOL = parse(Float64, get(ENV, "ED_TOL", "1e-11"))
const WARMUP_N = parse(Int, get(ENV, "WARMUP_N", "8"))

function heisenberg_periodic_pop(n::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(ComplexF64(0.25) * op[i] * op[mod1(i + 1, n)]
            for op in (σx, σy, σz) for i in 1:n)
    return polyopt(H, registry)
end

function build_manidsdp_problem(n::Integer; order::Integer, rank::Integer)
    optimizer = MOI.OptimizerWithAttributes(
        ManiDSDP.Optimizer,
        MOI.RawOptimizerAttribute("rank") => rank,
        MOI.RawOptimizerAttribute("max_rank") => rank,
        MOI.Silent() => true,
    )
    t_pop = @elapsed pop = heisenberg_periodic_pop(n)
    cfg = SolverConfig(optimizer=optimizer, order=order)
    t_sparsity = @elapsed sparsity = compute_sparsity(pop, cfg)
    t_moment = @elapsed mp = NCTSSoS.moment_relax(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    ext = Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt)
    ext === nothing && error("NCTSSoSManiDSDPExt is not loaded")
    opt = ext._manidsdp_optimizer_from_attributes(optimizer)
    t_problem = @elapsed problem, constant = ext._pauli_moment_problem_to_manidsdp(mp, opt)
    return (; problem, constant, t_pop, t_sparsity, t_moment, t_problem)
end

function solve_manidsdp(problem; rank::Integer=MANIDSDP_RANK)
    rank = min(Int(rank), problem.n)
    GC.gc()
    t = @elapsed result = ManiDSDP.solve(
        problem;
        settings=Settings(max_iter=MANIDSDP_MAX_ITER, tol=MANIDSDP_TOL, verbose=false),
        backend=:cpu,
        rank=rank,
        max_rank=rank,
        sigma=MANIDSDP_SIGMA,
        inner_max_iter=MANIDSDP_INNER_MAX_ITER,
        inner_tol=MANIDSDP_INNER_TOL,
        inner_method=:cg,
        rank_update=false,
        truncate_rank=false,
    )
    return (; result, t, rank)
end

function half_filling_basis(n::Integer)
    n % 2 == 0 || error("half-filling ED requires even n; got $n")
    dim_full = 1 << n
    states = UInt64[]
    sizehint!(states, binomial(n, n ÷ 2))
    for s in UInt64(0):UInt64(dim_full - 1)
        count_ones(s) == n ÷ 2 && push!(states, s)
    end
    index = zeros(Int32, dim_full)
    @inbounds for (i, s) in enumerate(states)
        index[Int(s) + 1] = Int32(i)
    end
    return states, index
end

function heisenberg_sector_mul!(y::Vector{Float64}, x::Vector{Float64}, states::Vector{UInt64}, index::Vector{Int32}, n::Integer)
    fill!(y, 0.0)
    @inbounds for col in eachindex(states)
        s = states[col]
        amp = x[col]
        diag = 0.0
        for site in 0:(n - 1)
            next = site == n - 1 ? 0 : site + 1
            bi = (s >> site) & 0x01
            bj = (s >> next) & 0x01
            if bi == bj
                diag += 0.25
            else
                diag -= 0.25
                flipped = s ⊻ ((UInt64(1) << site) | (UInt64(1) << next))
                row = Int(index[Int(flipped) + 1])
                y[row] += 0.5 * amp
            end
        end
        y[col] += diag * amp
    end
    return y
end

function heisenberg_periodic_ed_lanczos(n::Integer; maxiter::Integer=ED_MAXITER, tol::Float64=ED_TOL, seed::Integer=1234)
    t_basis = @elapsed states, index = half_filling_basis(n)
    dim = length(states)
    rng = MersenneTwister(seed + n)
    q = randn(rng, dim)
    q ./= norm(q)
    qprev = zeros(Float64, dim)
    w = zeros(Float64, dim)
    Q = Matrix{Float64}(undef, dim, maxiter)
    alphas = zeros(Float64, maxiter)
    betas = zeros(Float64, maxiter)
    theta = NaN
    residual = Inf
    used = maxiter
    t_lanczos = @elapsed begin
        for k in 1:maxiter
            @views Q[:, k] .= q
            heisenberg_sector_mul!(w, q, states, index, n)
            if k > 1
                @. w -= betas[k - 1] * qprev
            end
            alpha = dot(q, w)
            alphas[k] = alpha
            @. w -= alpha * q
            for _ in 1:2
                for j in 1:k
                    qj = @view Q[:, j]
                    coeff = dot(qj, w)
                    @. w -= coeff * qj
                end
            end
            beta = norm(w)
            betas[k] = beta
            if k >= 4 && (k % 5 == 0 || beta <= eps(Float64) || k == maxiter)
                T = SymTridiagonal(alphas[1:k], betas[1:max(k - 1, 0)])
                F = eigen(T)
                idx = argmin(F.values)
                theta = F.values[idx]
                residual = beta * abs(F.vectors[end, idx])
                if residual <= tol || beta <= eps(Float64)
                    used = k
                    break
                end
            end
            qprev, q = q, w ./ beta
        end
    end
    return (; value=theta, basis_dim=dim, iterations=used, residual, t_basis, t_lanczos, t_total=t_basis + t_lanczos)
end

function print_manidsdp(label, built, solved)
    objective = solved.result.objective + built.constant
    total = built.t_pop + built.t_sparsity + built.t_moment + built.t_problem + solved.t
    @printf("%s objective=% .12f per_site=% .12f status=%s\n",
        label, objective, objective / CHAIN_N, solved.result.status)
    @printf("%s dims n=%d m=%d nnz=%d rank=%d kkt=%.3e inner_total=%d\n",
        label, built.problem.n, length(built.problem.b), length(built.problem.vals),
        solved.rank, solved.result.dual.kkt.max, solved.result.dual.inner_iterations)
    @printf("%s pop=%7.3f sparsity=%7.3f moment=%7.3f problem=%7.3f solve=%7.3f total=%7.3f\n",
        label, built.t_pop, built.t_sparsity, built.t_moment, built.t_problem, solved.t, total)
    return objective, total
end

function main()
    println("Heisenberg periodic chain: ManiDSDP order-$RELAX_ORDER vs ED")
    println("Julia ", VERSION, " threads=", Threads.nthreads(), " BLAS=", BLAS.get_num_threads(), " syscpu=", Sys.CPU_THREADS)
    println("N=", CHAIN_N, " warmup_N=", WARMUP_N)
    println("ManiDSDP settings: rank=", MANIDSDP_RANK,
            " max_iter=", MANIDSDP_MAX_ITER,
            " tol=", MANIDSDP_TOL,
            " sigma=", MANIDSDP_SIGMA,
            " inner_max_iter=", MANIDSDP_INNER_MAX_ITER,
            " inner_tol=", MANIDSDP_INNER_TOL)
    println("ED settings: half-filling Lanczos maxiter=", ED_MAXITER, " tol=", ED_TOL)
    println("extension loaded=", Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt) !== nothing)

    println("\nwarmup")
    warm_built = build_manidsdp_problem(WARMUP_N; order=RELAX_ORDER, rank=min(MANIDSDP_RANK, 120))
    warm_solved = solve_manidsdp(warm_built.problem; rank=min(MANIDSDP_RANK, warm_built.problem.n))
    warm_ed = heisenberg_periodic_ed_lanczos(WARMUP_N; maxiter=min(80, ED_MAXITER), tol=max(ED_TOL, 1e-10))
    @printf("warmup ManiDSDP N=%d solve=%.3f status=%s obj=%.12f\n", WARMUP_N, warm_solved.t, warm_solved.result.status, warm_solved.result.objective + warm_built.constant)
    @printf("warmup ED       N=%d time=%.3f value=%.12f residual=%.3e iters=%d dim=%d\n", WARMUP_N, warm_ed.t_total, warm_ed.value, warm_ed.residual, warm_ed.iterations, warm_ed.basis_dim)

    println("\nmeasured")
    built = build_manidsdp_problem(CHAIN_N; order=RELAX_ORDER, rank=MANIDSDP_RANK)
    solved = solve_manidsdp(built.problem)
    mani_obj, mani_total = print_manidsdp("ManiDSDP", built, solved)

    ed = heisenberg_periodic_ed_lanczos(CHAIN_N; maxiter=ED_MAXITER, tol=ED_TOL)
    @printf("ED       objective=% .12f per_site=% .12f residual=%.3e iters=%d dim=%d\n",
        ed.value, ed.value / CHAIN_N, ed.residual, ed.iterations, ed.basis_dim)
    @printf("ED       basis=%7.3f lanczos=%7.3f total=%7.3f\n", ed.t_basis, ed.t_lanczos, ed.t_total)

    gap = mani_obj - ed.value
    @printf("\ncomparison: objective gap ManiDSDP-ED=%+.6e per_site=%+.6e\n", gap, gap / CHAIN_N)
    @printf("timing: ManiDSDP solve / ED total = %.3fx, ManiDSDP total / ED total = %.3fx\n",
        solved.t / ed.t_total, mani_total / ed.t_total)
end

main()
