#!/usr/bin/env julia

using LinearAlgebra
using Printf

using JuMP
using NCTSSoS
using ManiDSDP
import MosekTools

const MOI = NCTSSoS.MOI

const CHAIN_N = parse(Int, get(ENV, "HEISENBERG_N", "6"))
const RELAX_ORDER = parse(Int, get(ENV, "HEISENBERG_ORDER", "2"))
const MANIDSDP_RANK = parse(Int, get(ENV, "MANIDSDP_RANK", "120"))
const MANIDSDP_MAX_ITER = parse(Int, get(ENV, "MANIDSDP_MAX_ITER", "6"))
const MANIDSDP_TOL = parse(Float64, get(ENV, "MANIDSDP_TOL", "1e-4"))
const MANIDSDP_SIGMA = parse(Float64, get(ENV, "MANIDSDP_SIGMA", "0.3"))
const MANIDSDP_INNER_MAX_ITER = parse(Int, get(ENV, "MANIDSDP_INNER_MAX_ITER", "50"))
const MANIDSDP_INNER_TOL = parse(Float64, get(ENV, "MANIDSDP_INNER_TOL", "1e-5"))

function heisenberg_periodic_pop(n::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n)]
            for op in (σx, σy, σz) for i in 1:n)
    return polyopt(H, registry)
end

function manidsdp_optimizer()
    return MOI.OptimizerWithAttributes(
        ManiDSDP.Optimizer,
        MOI.RawOptimizerAttribute("rank") => MANIDSDP_RANK,
        MOI.RawOptimizerAttribute("max_rank") => MANIDSDP_RANK,
        MOI.RawOptimizerAttribute("backend") => :cpu,
        MOI.RawOptimizerAttribute("max_iter") => MANIDSDP_MAX_ITER,
        MOI.RawOptimizerAttribute("tol") => MANIDSDP_TOL,
        MOI.RawOptimizerAttribute("sigma") => MANIDSDP_SIGMA,
        MOI.RawOptimizerAttribute("inner_max_iter") => MANIDSDP_INNER_MAX_ITER,
        MOI.RawOptimizerAttribute("inner_tol") => MANIDSDP_INNER_TOL,
        MOI.RawOptimizerAttribute("inner_method") => :cg,
        MOI.RawOptimizerAttribute("rank_update") => false,
        MOI.RawOptimizerAttribute("truncate_rank") => false,
        MOI.Silent() => true,
    )
end

function mosek_optimizer()
    return optimizer_with_attributes(
        MosekTools.Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, min(Threads.nthreads(), div(Sys.CPU_THREADS, 2))),
        "MSK_IPAR_LOG" => 0,
    )
end

function run_case(label::AbstractString, optimizer)
    GC.gc()
    t_total = @elapsed begin
        t_pop = @elapsed pop = heisenberg_periodic_pop(CHAIN_N)
        cfg = SolverConfig(optimizer=optimizer, order=RELAX_ORDER)
        t_sparsity = @elapsed sparsity = compute_sparsity(pop, cfg)
        t_moment = @elapsed mp = NCTSSoS.moment_relax(
            pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        t_solve = @elapsed solved = NCTSSoS.solve_sdp(
            mp,
            cfg.optimizer;
            dualize=true,
            formulation=:moment_variables,
            representation=:real,
            orphan_policy=:error,
        )
    end
    @printf("%-14s objective=% .12f per_site=% .12f status=%s nuniq=%d\n",
            label, solved.objective, solved.objective / CHAIN_N, solved.status, solved.n_unique_elements)
    @printf("%-14s pop=%7.4f sparsity=%7.4f moment=%7.4f solve=%7.4f total=%7.4f\n",
            label, t_pop, t_sparsity, t_moment, t_solve, t_total)
    return (; objective=Float64(real(solved.objective)), solve=t_solve, total=t_total)
end

function main()
    println("Heisenberg periodic chain N=", CHAIN_N, ", Pauli order=", RELAX_ORDER)
    println("Julia ", VERSION, " threads=", Threads.nthreads(), " syscpu=", Sys.CPU_THREADS)
    println("ManiDSDP settings: rank=", MANIDSDP_RANK,
            " max_iter=", MANIDSDP_MAX_ITER,
            " tol=", MANIDSDP_TOL,
            " sigma=", MANIDSDP_SIGMA,
            " inner_max_iter=", MANIDSDP_INNER_MAX_ITER,
            " inner_tol=", MANIDSDP_INNER_TOL)
    println("NCTSSoSManiDSDPExt loaded=", Base.get_extension(NCTSSoS, :NCTSSoSManiDSDPExt) !== nothing)
    println("\nwarmup")
    run_case("ManiDSDP", manidsdp_optimizer())
    run_case("MosekTools", mosek_optimizer())
    println("\nmeasured")
    mani = run_case("ManiDSDP", manidsdp_optimizer())
    mosek = run_case("MosekTools", mosek_optimizer())
    @printf("\nspeedup: solve %.3fx, total %.3fx, objective diff %.3e\n",
            mosek.solve / mani.solve, mosek.total / mani.total, mani.objective - mosek.objective)
    return nothing
end

main()
