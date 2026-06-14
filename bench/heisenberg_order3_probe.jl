#!/usr/bin/env julia

using LinearAlgebra
using Printf
using NCTSSoS
using ManiDSDP

const MOI = NCTSSoS.MOI
const CHAIN_N = parse(Int, get(ENV, "HEISENBERG_N", "16"))
const RELAX_ORDER = parse(Int, get(ENV, "HEISENBERG_ORDER", "3"))
const MANIDSDP_RANK = parse(Int, get(ENV, "MANIDSDP_RANK", "600"))
const BUILD_ONLY = lowercase(get(ENV, "BUILD_ONLY", "false")) in ("1", "true", "yes")
const DO_WARMUP = lowercase(get(ENV, "DO_WARMUP", "true")) in ("1", "true", "yes")
const TS_ALGO = Symbol(lowercase(get(ENV, "TS_ALGO", "none")))

function ts_algorithm(name::Symbol)
    name === :none && return NoElimination()
    name === :mmd && return MMD()
    name === :mf && return MF()
    error("unsupported TS_ALGO=$name; expected none, mmd, or mf")
end

function block_summary(block_sizes)
    flat = Int[]
    for clique in block_sizes
        append!(flat, clique)
    end
    isempty(flat) && return "empty"
    return @sprintf("nblocks=%d total_dim=%d max_block=%d sumsq=%d", length(flat), sum(flat), maximum(flat), sum(abs2, flat))
end

function heisenberg_periodic_pop(n::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(ComplexF64(0.25) * op[i] * op[mod1(i + 1, n)]
            for op in (σx, σy, σz) for i in 1:n)
    return polyopt(H, registry)
end

function build_problem(n::Integer, order::Integer, rank::Integer)
    optimizer = MOI.OptimizerWithAttributes(
        ManiDSDP.Optimizer,
        MOI.RawOptimizerAttribute("rank") => rank,
        MOI.RawOptimizerAttribute("max_rank") => rank,
        MOI.Silent() => true,
    )
    t_pop = @elapsed pop = heisenberg_periodic_pop(n)
    cfg = SolverConfig(optimizer=optimizer, order=order, ts_algo=ts_algorithm(TS_ALGO))
    t_sparsity = @elapsed sparsity = compute_sparsity(pop, cfg)
    block_sizes = NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)
    t_moment = @elapsed mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    ext = Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt)
    ext === nothing && error("NCTSSoSManiDSDPExt is not loaded")
    opt = ext._manidsdp_optimizer_from_attributes(optimizer)
    t_problem = @elapsed problem, constant = ext._pauli_moment_problem_to_manidsdp(mp, opt)
    return (; problem, constant, block_sizes, n_unique=mp.n_unique_moment_matrix_elements,
            t_pop, t_sparsity, t_moment, t_problem)
end

function solve_problem(problem)
    rank = min(MANIDSDP_RANK, problem.n)
    GC.gc()
    t = @elapsed result = ManiDSDP.solve(
        problem;
        settings=Settings(max_iter=parse(Int, get(ENV, "MANIDSDP_MAX_ITER", "2")),
                          tol=parse(Float64, get(ENV, "MANIDSDP_TOL", "1e-4")),
                          verbose=false),
        backend=:cpu,
        rank=rank,
        max_rank=rank,
        sigma=parse(Float64, get(ENV, "MANIDSDP_SIGMA", "0.12")),
        inner_max_iter=parse(Int, get(ENV, "MANIDSDP_INNER_MAX_ITER", "20")),
        inner_tol=parse(Float64, get(ENV, "MANIDSDP_INNER_TOL", "1e-5")),
        inner_method=:cg,
        rank_update=false,
        truncate_rank=false,
    )
    return (; result, t, rank)
end

function report(label, built)
    p = built.problem
    @printf("%s N=%d order=%d ts=%s %s nuniq=%d\n",
            label, CHAIN_N, RELAX_ORDER, String(TS_ALGO), block_summary(built.block_sizes), built.n_unique)
    @printf("%s ManiDSDP dims n=%d m=%d nnz=%d rank=%d\n",
            label, p.n, length(p.b), length(p.vals), min(MANIDSDP_RANK, p.n))
    @printf("%s build pop=%.3f sparsity=%.3f moment=%.3f problem=%.3f total=%.3f\n",
            label, built.t_pop, built.t_sparsity, built.t_moment, built.t_problem,
            built.t_pop + built.t_sparsity + built.t_moment + built.t_problem)
end

println("Heisenberg order-3 ManiDSDP probe")
println("Julia ", VERSION, " threads=", Threads.nthreads(), " BLAS=", BLAS.get_num_threads(), " syscpu=", Sys.CPU_THREADS)
println("N=", CHAIN_N, " order=", RELAX_ORDER, " ts=", TS_ALGO, " rank=", MANIDSDP_RANK, " build_only=", BUILD_ONLY)
println("extension loaded=", Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt) !== nothing)

if DO_WARMUP
    println("\nwarmup N=8 order=2")
    warm = build_problem(8, 2, min(MANIDSDP_RANK, 120))
    report("warmup", warm)
end

println("\nmeasured build")
built = build_problem(CHAIN_N, RELAX_ORDER, MANIDSDP_RANK)
report("measured", built)

if !BUILD_ONLY
    println("\nmeasured solve")
    solved = solve_problem(built.problem)
    objective = solved.result.objective + built.constant
    @printf("ManiDSDP objective=% .12f per_site=% .12f status=%s solve=%.3f kkt=%.3e inner_total=%d\n",
            objective, objective / CHAIN_N, solved.result.status, solved.t,
            solved.result.dual.kkt.max, solved.result.dual.inner_iterations)
end
