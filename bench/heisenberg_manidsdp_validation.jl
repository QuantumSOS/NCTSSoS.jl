#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using NCTSSoS
using ManiDSDP

const MOI = NCTSSoS.MOI

parse_float_list(raw::AbstractString) = [parse(Float64, strip(x)) for x in split(raw, ',') if !isempty(strip(x))]

const CHAIN_N = parse(Int, get(ENV, "HEISENBERG_N", "20"))
const BASIS_MODE = Symbol(lowercase(get(ENV, "HEISENBERG_BASIS", "order")))
const RELAX_ORDER = parse(Int, get(ENV, "HEISENBERG_ORDER", "2"))
const PAPER_LONG_RANGE_R_ENV = get(ENV, "PAPER_LONG_RANGE_R", "")
const PAPER_INCLUDE_P4 = lowercase(get(ENV, "PAPER_INCLUDE_P4", "false")) in ("1", "true", "yes")
const BUILD_ONLY = lowercase(get(ENV, "BUILD_ONLY", "false")) in ("1", "true", "yes")
const RUN_ED = lowercase(get(ENV, "RUN_ED", "true")) in ("1", "true", "yes")

const MANIDSDP_RANK = parse(Int, get(ENV, "MANIDSDP_RANK", "600"))
const MANIDSDP_MAX_RANK = parse(Int, get(ENV, "MANIDSDP_MAX_RANK", string(MANIDSDP_RANK)))
const MANIDSDP_MAX_ITER = parse(Int, get(ENV, "MANIDSDP_MAX_ITER", "12"))
const MANIDSDP_TOL = parse(Float64, get(ENV, "MANIDSDP_TOL", "1e-5"))
const MANIDSDP_VERBOSE = lowercase(get(ENV, "MANIDSDP_VERBOSE", "true")) in ("1", "true", "yes")
const MANIDSDP_SIGMAS = parse_float_list(get(ENV, "MANIDSDP_SIGMAS", get(ENV, "MANIDSDP_SIGMA", "0.12")))
const MANIDSDP_INNER_MAX_ITER = parse(Int, get(ENV, "MANIDSDP_INNER_MAX_ITER", "120"))
const MANIDSDP_INNER_TOL = parse(Float64, get(ENV, "MANIDSDP_INNER_TOL", "1e-6"))
const MANIDSDP_RANK_UPDATE = lowercase(get(ENV, "MANIDSDP_RANK_UPDATE", "false")) in ("1", "true", "yes")
const MANIDSDP_TRUNCATE_RANK = lowercase(get(ENV, "MANIDSDP_TRUNCATE_RANK", "false")) in ("1", "true", "yes")
const ED_MAXITER = parse(Int, get(ENV, "ED_MAXITER", "700"))
const ED_TOL = parse(Float64, get(ENV, "ED_TOL", "1e-12"))

paper_long_range_r(n::Integer) = isempty(PAPER_LONG_RANGE_R_ENV) ? (n <= 60 ? n ÷ 2 : min(20, n ÷ 2)) : parse(Int, PAPER_LONG_RANGE_R_ENV)

function heisenberg_periodic_pop(n::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(ComplexF64(0.25) * op[i] * op[mod1(i + 1, n)]
            for op in (σx, σy, σz) for i in 1:n)
    return polyopt(H, registry), (σx, σy, σz)
end

function pauli_monomial(ops)
    first_op = first(ops)
    T = eltype(first_op.word)
    word = T[]
    for op in ops
        append!(word, op.word)
    end
    normal_word, _phase = simplify(NCTSSoS.PauliAlgebra, word)
    return NormalMonomial{NCTSSoS.PauliAlgebra,T}(normal_word)
end

function add_site_pattern!(basis::Set, paulis, sites)
    ranges = ntuple(_ -> 1:3, length(sites))
    for types in Iterators.product(ranges...)
        ops = ntuple(k -> paulis[types[k]][sites[k]], length(sites))
        push!(basis, pauli_monomial(ops))
    end
    return basis
end

function paper_1d_heisenberg_basis(paulis, n::Integer; r::Integer=paper_long_range_r(n), include_p4::Bool=PAPER_INCLUDE_P4)
    r < 1 && throw(ArgumentError("long-range cutoff r must be at least 1, got $r"))
    r <= n ÷ 2 || throw(ArgumentError("long-range cutoff r=$r exceeds n÷2=$(n ÷ 2)"))
    basis = Set{typeof(paulis[1][1])}()
    counts = Pair{Symbol,Int}[]

    function mark!(f::Function, label::Symbol)
        before = length(basis)
        f()
        push!(counts, label => (length(basis) - before))
        return nothing
    end

    mark!(:P0_identity) do
        push!(basis, one(paulis[1][1]))
    end
    mark!(:P1_single_site) do
        for i in 1:n, a in 1:3
            push!(basis, paulis[a][i])
        end
    end
    mark!(:P2_nearest_neighbor) do
        for i in 1:n
            add_site_pattern!(basis, paulis, (i, mod1(i + 1, n)))
        end
    end
    mark!(:P2_long_range) do
        for dist in 2:r, i in 1:n
            add_site_pattern!(basis, paulis, (i, mod1(i + dist, n)))
        end
    end
    mark!(:P3_contiguous) do
        for i in 1:n
            add_site_pattern!(basis, paulis, (i, mod1(i + 1, n), mod1(i + 2, n)))
        end
    end
    if include_p4
        mark!(:P4_contiguous) do
            for i in 1:n
                add_site_pattern!(basis, paulis, (i, mod1(i + 1, n), mod1(i + 2, n), mod1(i + 3, n)))
            end
        end
    else
        push!(counts, :P4_contiguous => 0)
    end

    return sort!(collect(basis)), counts
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
            k > 1 && (@. w -= betas[k - 1] * qprev)
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

function block_summary(block_sizes)
    flat = Int[]
    for clique in block_sizes
        append!(flat, clique)
    end
    isempty(flat) && return "empty"
    return @sprintf("nblocks=%d total_dim=%d max_block=%d sumsq=%d", length(flat), sum(flat), maximum(flat), sum(abs2, flat))
end

function build_problem(n::Integer)
    optimizer = MOI.OptimizerWithAttributes(
        ManiDSDP.Optimizer,
        MOI.RawOptimizerAttribute("rank") => MANIDSDP_RANK,
        MOI.RawOptimizerAttribute("max_rank") => MANIDSDP_RANK,
        MOI.Silent() => !MANIDSDP_VERBOSE,
    )

    t_pop = @elapsed pop, paulis = heisenberg_periodic_pop(n)
    basis_label = ""
    counts = Pair{Symbol,Int}[]
    cfg = if BASIS_MODE === :order
        basis_label = "order-$RELAX_ORDER"
        SolverConfig(optimizer=optimizer, order=RELAX_ORDER, cs_algo=NoElimination(), ts_algo=NoElimination())
    elseif BASIS_MODE === :paper
        r = paper_long_range_r(n)
        t_basis = @elapsed basis, counts = paper_1d_heisenberg_basis(paulis, n; r=r, include_p4=PAPER_INCLUDE_P4)
        basis_label = "paper-r$r-P4-$(PAPER_INCLUDE_P4)"
        println(@sprintf("basis_build %.3fs", t_basis))
        SolverConfig(optimizer=optimizer, moment_basis=basis, cs_algo=NoElimination(), ts_algo=NoElimination())
    else
        throw(ArgumentError("HEISENBERG_BASIS must be order or paper, got $BASIS_MODE"))
    end

    t_sparsity = @elapsed sparsity = compute_sparsity(pop, cfg)
    block_sizes = NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)
    t_moment = @elapsed mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    ext = Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt)
    ext === nothing && error("NCTSSoSManiDSDPExt is not loaded")
    opt = ext._manidsdp_optimizer_from_attributes(optimizer)
    t_problem = @elapsed problem, constant = ext._pauli_moment_problem_to_manidsdp(mp, opt)
    return (; problem, constant, basis_label, counts, block_sizes,
            n_unique=mp.n_unique_moment_matrix_elements,
            t_pop, t_sparsity, t_moment, t_problem)
end

function solve_problem(problem; sigma::Float64)
    rank = min(MANIDSDP_RANK, problem.n)
    max_rank = min(MANIDSDP_MAX_RANK, problem.n)
    rank <= max_rank || throw(ArgumentError("MANIDSDP_RANK=$rank exceeds MANIDSDP_MAX_RANK=$max_rank"))
    GC.gc()
    t = @elapsed result = ManiDSDP.solve(
        problem;
        settings=Settings(max_iter=MANIDSDP_MAX_ITER, tol=MANIDSDP_TOL, verbose=MANIDSDP_VERBOSE),
        backend=:cpu,
        rank=rank,
        max_rank=max_rank,
        sigma=sigma,
        inner_max_iter=MANIDSDP_INNER_MAX_ITER,
        inner_tol=MANIDSDP_INNER_TOL,
        inner_method=:cg,
        rank_update=MANIDSDP_RANK_UPDATE,
        truncate_rank=MANIDSDP_TRUNCATE_RANK,
    )
    return (; result, t, rank, max_rank, sigma)
end

function main()
    println("Heisenberg XXX periodic validation: ManiDSDP vs ED")
    println("Julia ", VERSION, " threads=", Threads.nthreads(), " BLAS=", BLAS.get_num_threads(), " syscpu=", Sys.CPU_THREADS)
    println("N=", CHAIN_N, " basis=", BASIS_MODE, " order=", RELAX_ORDER,
            " paper_r=", paper_long_range_r(CHAIN_N), " include_p4=", PAPER_INCLUDE_P4,
            " build_only=", BUILD_ONLY)
    println("ManiDSDP rank=", MANIDSDP_RANK, " max_rank=", MANIDSDP_MAX_RANK,
            " max_iter=", MANIDSDP_MAX_ITER,
            " tol=", MANIDSDP_TOL, " sigmas=", MANIDSDP_SIGMAS,
            " inner_max_iter=", MANIDSDP_INNER_MAX_ITER, " inner_tol=", MANIDSDP_INNER_TOL,
            " rank_update=", MANIDSDP_RANK_UPDATE, " truncate_rank=", MANIDSDP_TRUNCATE_RANK,
            " verbose=", MANIDSDP_VERBOSE)
    println("extension loaded=", Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt) !== nothing)

    ed = nothing
    if RUN_ED
        println("\nED")
        ed = heisenberg_periodic_ed_lanczos(CHAIN_N; maxiter=ED_MAXITER, tol=ED_TOL)
        @printf("ED objective=% .12f per_site=% .12f residual=%.3e iters=%d dim=%d total=%.3f\n",
                ed.value, ed.value / CHAIN_N, ed.residual, ed.iterations, ed.basis_dim, ed.t_total)
        flush(stdout)
    end

    println("\nbuild")
    built = build_problem(CHAIN_N)
    println("basis_label=", built.basis_label)
    isempty(built.counts) || println("basis_counts=", join(["$(k):$(v)" for (k, v) in built.counts], ","))
    println("moment ", block_summary(built.block_sizes), " nuniq=", built.n_unique)
    println("ManiDSDP dims n=", built.problem.n, " m=", length(built.problem.b),
            " nnz=", length(built.problem.vals), " rank=", min(MANIDSDP_RANK, built.problem.n))
    @printf("build_times pop=%.3f sparsity=%.3f moment=%.3f problem=%.3f total=%.3f\n",
            built.t_pop, built.t_sparsity, built.t_moment, built.t_problem,
            built.t_pop + built.t_sparsity + built.t_moment + built.t_problem)
    flush(stdout)

    BUILD_ONLY && return nothing

    for sigma in MANIDSDP_SIGMAS
        println("\nsolve sigma=", sigma)
        solved = solve_problem(built.problem; sigma=sigma)
        objective = solved.result.objective + built.constant
        @printf("ManiDSDP sigma=%.6g objective=% .12f per_site=% .12f status=%s iterations=%d solve=%.3f rank=%d max_rank=%d kkt=%.3e inner_total=%d\n",
                sigma, objective, objective / CHAIN_N, solved.result.status, solved.result.iterations, solved.t, solved.rank, solved.max_rank,
                solved.result.dual.kkt.max, solved.result.dual.inner_iterations)
        if ed !== nothing
            gap = objective - ed.value
            @printf("comparison sigma=%.6g gap ManiDSDP-ED=%+.6e per_site=%+.6e\n", sigma, gap, gap / CHAIN_N)
        end
        flush(stdout)
    end
end

main()
