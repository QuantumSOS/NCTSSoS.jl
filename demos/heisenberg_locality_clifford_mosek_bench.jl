#!/usr/bin/env julia

# Benchmark periodic XXX Heisenberg chains with a locality-restricted Pauli basis
# and optional Clifford/Wedderburn symmetry reduction.
#
# Default run is intentionally conservative: N=30, r=2, contiguous clusters up
# to degree 4, spatial dihedral symmetry only. Dense/no-symmetry is skipped when
# the basis is too large unless HEISENBERG_MAX_NOSYM_BASIS is raised.
#
# Usage:
#   MOSEK_THREADS=64 julia --project demos/heisenberg_locality_clifford_mosek_bench.jl 30
#
# Useful environment:
#   HEISENBERG_LOCAL_R=2                 cyclic two-body distance cutoff
#   HEISENBERG_CLUSTER_MAX_DEG=4         add contiguous clusters of degree 3..this value; use 2 for degree<=2 only
#   HEISENBERG_RUN_SYMMETRY=1            run the symmetry-reduced SDP
#   HEISENBERG_RUN_NOSYM=1               run the dense/no-symmetry SDP for comparison
#   HEISENBERG_MAX_NOSYM_BASIS=1200      skip no-symmetry solve above this basis size
#   HEISENBERG_SPIN_CLIFFORD=none        none | cyclic | full
#   HEISENBERG_TIMELIMIT_SEC=0           Mosek time limit; 0 means unset
#   MOSEK_THREADS=64                     number of Mosek threads

using Dates
using JuMP
using MosekTools
using NCTSSoS
using Printf

const LOCAL_R = parse(Int, get(ENV, "HEISENBERG_LOCAL_R", "2"))
const CLUSTER_MAX_DEG = parse(Int, get(ENV, "HEISENBERG_CLUSTER_MAX_DEG", "4"))
const RUN_SYMMETRY = get(ENV, "HEISENBERG_RUN_SYMMETRY", "1") in ("1", "true", "yes")
const RUN_NOSYM = get(ENV, "HEISENBERG_RUN_NOSYM", "1") in ("1", "true", "yes")
const MAX_NOSYM_BASIS = parse(Int, get(ENV, "HEISENBERG_MAX_NOSYM_BASIS", "1200"))
const SPIN_CLIFFORD = lowercase(get(ENV, "HEISENBERG_SPIN_CLIFFORD", "none"))
const MOSEK_THREADS = parse(Int, get(ENV, "MOSEK_THREADS", string(max(1, Sys.CPU_THREADS ÷ 2))))
const MOSEK_TIMELIMIT = parse(Float64, get(ENV, "HEISENBERG_TIMELIMIT_SEC", "0"))
const DUALIZE = get(ENV, "HEISENBERG_DUALIZE", "1") in ("1", "true", "yes")

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

function global_axis_cycle_clifford(vars)
    # 120-degree spin rotation: X -> Y -> Z -> X. This gives a small global
    # spin symmetry (order 3) and keeps D_N × C_3 at order 6N.
    σx, σy, σz = vars
    nsites = length(σx)
    T = eltype(σx).parameters[2]
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:nsites
        images[σx[site]] = (1, σy[site])
        images[σy[site]] = (1, σz[site])
        images[σz[site]] = (1, σx[site])
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

    if SPIN_CLIFFORD == "cyclic"
        push!(generators, global_axis_cycle_clifford(vars))
    elseif SPIN_CLIFFORD == "full"
        # Full single-site Clifford action on the spin axes. For N=30 this is
        # usually D_30 × O ≈ 1440 elements, not the cheap 180-element case.
        push!(generators, global_h_clifford(vars), global_s_clifford(vars))
    elseif SPIN_CLIFFORD != "none"
        throw(ArgumentError("HEISENBERG_SPIN_CLIFFORD must be none, cyclic, or full; got '$SPIN_CLIFFORD'."))
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

cyclic_distance(i::Int, j::Int, nsites::Int) = min(abs(i - j), nsites - abs(i - j))

function monomial_from_letters(letters)
    isempty(letters) && throw(ArgumentError("internal error: empty non-identity monomial"))
    poly = reduce(*, letters)
    poly_terms = collect(terms(poly))
    length(poly_terms) == 1 || throw(ArgumentError("expected a single Pauli monomial, got $(length(poly_terms)) terms"))
    coef, mono = only(poly_terms)
    isone(coef) || throw(ArgumentError("locality basis should contain unit-coefficient monomials; got coefficient $coef for $mono"))
    return mono
end

function add_cluster_monomials!(basis, paulis, sites::AbstractVector{Int})
    n = length(sites)
    choices = fill(1, n)
    while true
        push!(basis, monomial_from_letters([paulis[choice][site] for (choice, site) in zip(choices, sites)]))

        pos = n
        while pos >= 1 && choices[pos] == 3
            choices[pos] = 1
            pos -= 1
        end
        pos == 0 && break
        choices[pos] += 1
    end
    return basis
end

function locality_restricted_basis(vars, nsites::Int; r::Int=LOCAL_R, max_cluster_deg::Int=CLUSTER_MAX_DEG)
    r >= 1 || throw(ArgumentError("HEISENBERG_LOCAL_R must be positive; got $r."))
    max_cluster_deg >= 2 || throw(ArgumentError("HEISENBERG_CLUSTER_MAX_DEG must be at least 2; got $max_cluster_deg."))
    max_cluster_deg <= nsites || throw(ArgumentError("cluster max degree $max_cluster_deg exceeds N=$nsites."))

    σx, σy, σz = vars
    T = eltype(σx).parameters[2]
    basis = NormalMonomial{PauliAlgebra,T}[one(NormalMonomial{PauliAlgebra,T})]

    # Degree 1: all single-site Pauli letters.
    for op in (σx, σy, σz), site in 1:nsites
        push!(basis, op[site])
    end

    # Degree 2: all distinct site pairs up to cyclic distance r.
    for i in 1:nsites-1, j in i+1:nsites
        cyclic_distance(i, j, nsites) <= r || continue
        for op_i in (σx, σy, σz), op_j in (σx, σy, σz)
            push!(basis, monomial_from_letters((op_i[i], op_j[j])))
        end
    end

    # Degree 3/4/...: contiguous cyclic clusters, as in the paper's 1D basis.
    # The Set/unique pass removes the expected wrap-around duplicates for k=N.
    if max_cluster_deg >= 3
        for k in 3:max_cluster_deg
            for start in 1:nsites
                sites = [mod1(start + offset, nsites) for offset in 0:k-1]
                add_cluster_monomials!(basis, (σx, σy, σz), sites)
            end
        end
    end

    sort!(unique!(basis))
    return basis
end

function block_summary(blocks::AbstractVector{<:Integer})
    isempty(blocks) && return "[]"
    sorted = sort(collect(blocks); rev=true)
    shown = join(sorted[1:min(end, 20)], ",")
    suffix = length(sorted) > 20 ? ",..." : ""
    return "count=$(length(blocks)) max=$(maximum(blocks)) sumsq=$(sum(abs2, blocks)) top_desc=[$shown$suffix]"
end

function peak_rss_mib()
    return Sys.maxrss() / 2.0^20
end

function build_problem(nsites::Int)
    registry, vars = create_pauli_variables(1:nsites)
    σx, σy, σz = vars
    hamiltonian = sum(
        ComplexF64(1 / 4) * op[site] * op[mod1(site + 1, nsites)]
        for op in (σx, σy, σz), site in 1:nsites
    )
    pop = polyopt(hamiltonian, registry)
    basis = locality_restricted_basis(vars, nsites)
    return registry, vars, pop, basis
end

function solve_mode(pop, vars, basis, nsites::Int, mode::Symbol)
    symmetry_spec = mode == :symmetry ? heisenberg_symmetry(vars, nsites) : nothing
    config = SolverConfig(
        optimizer=mosek_optimizer(),
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry_spec,
    )

    @printf("mode=%s phase=compute_sparsity start=%s peak_rss_mib=%.1f\n", mode, now(), peak_rss_mib())
    flush(stdout)
    sparsity_time = @elapsed sparsity = compute_sparsity(pop, config)

    @printf("mode=%s phase=moment_relax start=%s peak_rss_mib=%.1f\n", mode, now(), peak_rss_mib())
    flush(stdout)
    symmetry_report = nothing
    if mode == :symmetry
        relax_time = @elapsed begin
            moment_problem, symmetry_report = NCTSSoS.moment_relax_symmetric(
                pop,
                sparsity.corr_sparsity,
                sparsity.cliques_term_sparsities,
                symmetry_spec,
            )
        end
    else
        relax_time = @elapsed moment_problem = NCTSSoS.moment_relax(
            pop,
            sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities,
        )
    end

    @printf("mode=%s phase=solve_sdp start=%s peak_rss_mib=%.1f\n", mode, now(), peak_rss_mib())
    flush(stdout)
    solve_time = @elapsed raw_result = NCTSSoS.solve_sdp(moment_problem, config.optimizer; dualize=DUALIZE)

    status = JuMP.termination_status(raw_result.model)
    primal = JuMP.primal_status(raw_result.model)
    dual = JuMP.dual_status(raw_result.model)
    pre_blocks = NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)

    println("mode=$mode status=$status primal=$primal dual=$dual")
    @printf("mode=%s objective=%.16g\n", mode, raw_result.objective)
    @printf("mode=%s energy_per_site=%.16g\n", mode, raw_result.objective / nsites)
    println("mode=$mode pre_symmetry_moment_matrix_sizes=$(pre_blocks)")
    println("mode=$mode unique_moment_elements=$(raw_result.n_unique_elements)")
    if isnothing(symmetry_report)
        flat_blocks = isempty(pre_blocks) ? Int[] : reduce(vcat, pre_blocks)
        println("mode=$mode psd_blocks_no_symmetry $(block_summary(flat_blocks))")
    else
        println("mode=$mode group_order=$(symmetry_report.group_order)")
        println("mode=$mode invariant_moment_count=$(symmetry_report.invariant_moment_count)")
        println("mode=$mode psd_blocks_symmetry $(block_summary(symmetry_report.psd_block_sizes))")
        println("mode=$mode psd_block_sizes=$(symmetry_report.psd_block_sizes)")
    end
    @printf("mode=%s sparsity_sec=%.3f moment_relax_sec=%.3f construction_sec=%.3f solve_sec=%.3f total_sec=%.3f peak_rss_mib=%.1f\n",
        mode, sparsity_time, relax_time, sparsity_time + relax_time, solve_time,
        sparsity_time + relax_time + solve_time, peak_rss_mib())
    flush(stdout)

    return raw_result, symmetry_report
end

function run_case(nsites::Int)
    println("\n=== Heisenberg XXX periodic chain locality benchmark: N=$nsites ===")
    println("started_at=$(now()) host=$(gethostname()) julia=$(VERSION) julia_threads=$(Threads.nthreads()) sys_cpu_threads=$(Sys.CPU_THREADS)")
    println("local_r=$LOCAL_R cluster_max_deg=$CLUSTER_MAX_DEG spin_clifford=$SPIN_CLIFFORD dualize=$DUALIZE mosek_threads=$MOSEK_THREADS mosek_timelimit_sec=$MOSEK_TIMELIMIT")
    if nsites == 30
        println("paper_reference_N30_per_site=-0.444151 paper_reference_basis='r=15 plus contiguous degree 3 and 4 clusters'")
    end
    flush(stdout)

    setup_time = @elapsed registry, vars, pop, basis = build_problem(nsites)
    println("basis_size=$(length(basis)) setup_sec=$(round(setup_time; digits=3)) peak_rss_mib=$(round(peak_rss_mib(); digits=1))")
    println("basis_counts identity=1 degree1=$(3 * nsites) degree2_radius=$LOCAL_R cluster_max_deg=$CLUSTER_MAX_DEG")
    flush(stdout)

    if RUN_NOSYM
        if length(basis) <= MAX_NOSYM_BASIS
            try
                solve_mode(pop, vars, basis, nsites, :nosym)
            catch err
                println("mode=nosym failed error=$(typeof(err)): $err")
                showerror(stdout, err, catch_backtrace())
                println()
                flush(stdout)
            end
        else
            println("mode=nosym skipped reason='basis_size $(length(basis)) exceeds HEISENBERG_MAX_NOSYM_BASIS=$MAX_NOSYM_BASIS'")
            flush(stdout)
        end
    end

    if RUN_SYMMETRY
        try
            solve_mode(pop, vars, basis, nsites, :symmetry)
        catch err
            println("mode=symmetry failed error=$(typeof(err)): $err")
            showerror(stdout, err, catch_backtrace())
            println()
            flush(stdout)
        end
    end

    println("finished_at=$(now()) N=$nsites peak_rss_mib=$(round(peak_rss_mib(); digits=1))")
    flush(stdout)
end

function main(args)
    ns = isempty(args) ? [30] : parse.(Int, args)
    for nsites in ns
        run_case(nsites)
    end
end

main(ARGS)
