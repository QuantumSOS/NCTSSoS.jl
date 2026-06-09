#!/usr/bin/env julia

# Benchmark periodic XXX Heisenberg chains with Pauli Clifford site symmetries.
#
# Usage:
#   julia --project demos/heisenberg_clifford_mosek_bench.jl 6 8 10
#
# Environment:
#   MOSEK_THREADS=64                 number of Mosek threads
#   HEISENBERG_ORDER=2               relaxation/moment-basis degree
#   HEISENBERG_GLOBAL_SPIN_CLIFFORD=0 add global single-qubit Clifford spin rotations
#   HEISENBERG_TIMELIMIT_SEC=0       Mosek time limit; 0 means unset

using Dates
using JuMP
using MosekTools
using NCTSSoS
using Printf

const ORDER = parse(Int, get(ENV, "HEISENBERG_ORDER", "2"))
const MOSEK_THREADS = parse(Int, get(ENV, "MOSEK_THREADS", string(max(1, Sys.CPU_THREADS ÷ 2))))
const MOSEK_TIMELIMIT = parse(Float64, get(ENV, "HEISENBERG_TIMELIMIT_SEC", "0"))
const GLOBAL_SPIN_CLIFFORD = get(ENV, "HEISENBERG_GLOBAL_SPIN_CLIFFORD", "0") in ("1", "true", "yes")

function site_permutation_clifford(vars, perm::AbstractVector{<:Integer})
    σx, σy, σz = vars
    nsites = length(perm)
    T = eltype(σx).parameters[2]
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()

    for site in 1:nsites
        image_site = Int(perm[site])
        image_site == site && continue
        images[σx[site]] = (1, σx[image_site])
        images[σy[site]] = (1, σy[image_site])
        images[σz[site]] = (1, σz[image_site])
    end

    return CliffordSymmetry(images; nqubits=nsites)
end

function global_h_clifford(vars)
    σx, σy, σz = vars
    nsites = length(σx)
    T = eltype(σx).parameters[2]
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:nsites
        images[σx[site]] = (1, σz[site])
        images[σy[site]] = (-1, σy[site])
        images[σz[site]] = (1, σx[site])
    end
    return CliffordSymmetry(images; nqubits=nsites)
end

function global_s_clifford(vars)
    σx, σy, σz = vars
    nsites = length(σx)
    T = eltype(σx).parameters[2]
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:nsites
        images[σx[site]] = (1, σy[site])
        images[σy[site]] = (-1, σx[site])
    end
    return CliffordSymmetry(images; nqubits=nsites)
end

function heisenberg_symmetry(vars, nsites::Int)
    rotation = site_permutation_clifford(vars, [mod1(site + 1, nsites) for site in 1:nsites])
    reflection = site_permutation_clifford(vars, [nsites + 1 - site for site in 1:nsites])
    generators = CliffordSymmetry[rotation, reflection]
    if GLOBAL_SPIN_CLIFFORD
        push!(generators, global_h_clifford(vars), global_s_clifford(vars))
    end
    return SymmetrySpec(generators)
end

function mosek_optimizer()
    attrs = Pair{String,Any}[
        "MSK_IPAR_NUM_THREADS" => MOSEK_THREADS,
        "MSK_IPAR_LOG" => 0,
    ]
    if MOSEK_TIMELIMIT > 0
        push!(attrs, "MSK_DPAR_OPTIMIZER_MAX_TIME" => MOSEK_TIMELIMIT)
    end
    return optimizer_with_attributes(Mosek.Optimizer, attrs...)
end

function run_case(nsites::Int; order::Int=ORDER)
    started_at = now()
    println("\n=== Heisenberg XXX periodic chain: N=$nsites, order=$order ===")
    println("started_at=$started_at julia_threads=$(Threads.nthreads()) mosek_threads=$MOSEK_THREADS global_spin_clifford=$GLOBAL_SPIN_CLIFFORD")
    flush(stdout)

    registry, vars = create_pauli_variables(1:nsites)
    σx, σy, σz = vars
    hamiltonian = sum(
        ComplexF64(1 / 4) * op[site] * op[mod1(site + 1, nsites)]
        for op in (σx, σy, σz), site in 1:nsites
    )
    pop = polyopt(hamiltonian, registry)
    basis = get_ncbasis(registry, order)

    config = SolverConfig(
        optimizer=mosek_optimizer(),
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=heisenberg_symmetry(vars, nsites),
    )

    elapsed = @elapsed result = cs_nctssos(pop, config)
    status = JuMP.termination_status(result.model)
    primal = JuMP.primal_status(result.model)
    dual = JuMP.dual_status(result.model)
    blocks = result.symmetry.psd_block_sizes

    println("status=$status primal=$primal dual=$dual")
    @printf("objective=%.16g\n", result.objective)
    @printf("energy_per_site=%.16g\n", result.objective / nsites)
    println("basis_size=$(length(basis))")
    println("group_order=$(result.symmetry.group_order)")
    println("invariant_moment_count=$(result.symmetry.invariant_moment_count)")
    println("psd_block_sizes=$(blocks)")
    println("psd_block_count=$(length(blocks)) max_psd_block=$(maximum(blocks)) sum_psd_block_squares=$(sum(abs2, blocks))")
    @printf("elapsed_sec=%.3f\n", elapsed)
    println("finished_at=$(now())")
    flush(stdout)

    return result
end

function main(args)
    ns = isempty(args) ? [6, 8, 10] : parse.(Int, args)
    println("host=$(gethostname()) sys_cpu_threads=$(Sys.CPU_THREADS) order=$ORDER mosek_threads=$MOSEK_THREADS global_spin_clifford=$GLOBAL_SPIN_CLIFFORD")
    for nsites in ns
        try
            run_case(nsites)
        catch err
            println("case_failed N=$nsites error=$(typeof(err)): $err")
            showerror(stdout, err, catch_backtrace())
            println()
            flush(stdout)
            rethrow()
        end
    end
end

main(ARGS)
