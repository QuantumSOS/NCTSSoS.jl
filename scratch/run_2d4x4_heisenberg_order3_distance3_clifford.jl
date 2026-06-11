#!/usr/bin/env julia

# 4x4 periodic spin-1/2 Heisenberg model at a custom order-3 Pauli basis:
# include only Pauli strings whose occupied lattice sites have pairwise torus
# Manhattan distance <= DISTANCE_CUTOFF. Then apply Clifford symmetry reduction.

using Dates
using Printf
using JuMP
using NCTSSoS
import Clarabel

function _import_optional_package(name::Symbol)
    try
        @eval import $(name)
        return true
    catch
        return isdefined(Main, name)
    end
end

const Lx = parse(Int, get(ENV, "LX", "4"))
const Ly = parse(Int, get(ENV, "LY", "4"))
const N = Lx * Ly
const ORDER = parse(Int, get(ENV, "ORDER", "3"))
const DISTANCE_CUTOFF = parse(Int, get(ENV, "DISTANCE_CUTOFF", "3"))
const SOLVE = get(ENV, "SOLVE", "1") != "0"
const SYMMETRY_MODE = lowercase(get(ENV, "SYMMETRY_MODE", "full")) # lattice, spin, full
const OFFBLOCK_CHECK = Symbol(get(ENV, "OFFBLOCK_CHECK", "off"))
const PRINT_FULL_BLOCKS = get(ENV, "PRINT_FULL_BLOCKS", "0") == "1"
const FAST_DENSE_SPARSITY = get(ENV, "FAST_DENSE_SPARSITY", "1") != "0"
const DUALIZE = get(ENV, "DUALIZE", "1") != "0"
const FORMULATION = Symbol(get(ENV, "FORMULATION", "moment_variables"))
const REPRESENTATION = Symbol(get(ENV, "REPRESENTATION", "real"))

const LI = LinearIndices((1:Lx, 1:Ly))
site(i, j) = LI[CartesianIndex(mod1(i, Lx), mod1(j, Ly))]
coord(s) = Tuple(CartesianIndices((1:Lx, 1:Ly))[s])

function stage(msg)
    println("[", Dates.now(), "] ", msg)
    flush(stdout)
    return nothing
end

function solver_optimizer()
    solver = lowercase(get(ENV, "SOLVER", _import_optional_package(:MosekTools) ? "mosek" : "clarabel"))
    if solver == "mosek"
        _import_optional_package(:MosekTools) || throw(ArgumentError("SOLVER=mosek requested but MosekTools is unavailable."))
        return Base.invokelatest(getproperty, Main, :MosekTools).Optimizer
    elseif solver == "clarabel"
        return optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    elseif solver == "cosmo"
        _import_optional_package(:COSMO) || throw(ArgumentError("SOLVER=cosmo requested but COSMO is unavailable."))
        return optimizer_with_attributes(
            Base.invokelatest(getproperty, Main, :COSMO).Optimizer,
            "verbose" => false,
            "eps_abs" => 1e-7,
            "eps_rel" => 1e-7,
            "eps_prim_inf" => 1e-6,
            "eps_dual_inf" => 1e-6,
            "max_iter" => 200_000,
        )
    else
        throw(ArgumentError("SOLVER must be mosek, clarabel, or cosmo; got $solver"))
    end
end

function square_lattice_bonds()
    bonds = Tuple{Int,Int}[]
    sizehint!(bonds, 2N)
    for i in 1:Lx, j in 1:Ly
        u = site(i, j)
        push!(bonds, (u, site(i + 1, j)))
        push!(bonds, (u, site(i, j + 1)))
    end
    return bonds
end

const BONDS = square_lattice_bonds()

function heisenberg_pauli_hamiltonian()
    registry, (σx, σy, σz) = create_pauli_variables(1:N)
    H = sum(
        ComplexF64(1 / 4) * op[u] * op[v]
        for (u, v) in BONDS for op in (σx, σy, σz)
    )
    return registry, (σx, σy, σz), H
end

function torus_distance(a::Int, b::Int)
    ax, ay = coord(a)
    bx, by = coord(b)
    dx = abs(ax - bx)
    dy = abs(ay - by)
    return min(dx, Lx - dx) + min(dy, Ly - dy)
end

function site_diameter(sites::AbstractVector{Int})
    diam = 0
    for j in eachindex(sites), i in firstindex(sites):(j - 1)
        diam = max(diam, torus_distance(sites[i], sites[j]))
    end
    return diam
end

function visit_combinations!(f, n::Int, k::Int, start::Int=1, current::Vector{Int}=Int[])
    if length(current) == k
        f(current)
        return nothing
    end
    remaining = k - length(current)
    for x in start:(n - remaining + 1)
        push!(current, x)
        visit_combinations!(f, n, k, x + 1, current)
        pop!(current)
    end
    return nothing
end

pauli_var_index(::Type{T}, site_index::Int, axis::Int) where {T<:Integer} = T(3 * (site_index - 1) + axis)

function push_axis_assignments!(basis::Vector{NormalMonomial{PauliAlgebra,T}}, sites::Vector{Int}) where {T<:Integer}
    k = length(sites)
    for code0 in 0:(3^k - 1)
        code = code0
        word = T[]
        sizehint!(word, k)
        for site_index in sites
            axis = mod(code, 3) + 1
            code = div(code, 3)
            push!(word, pauli_var_index(T, site_index, axis))
        end
        push!(basis, NormalMonomial{PauliAlgebra,T}(word))
    end
    return nothing
end

function distance_limited_pauli_basis(registry, order::Int, cutoff::Int)
    T = eltype(indices(registry))
    M = NormalMonomial{PauliAlgebra,T}
    basis = M[one(M)]
    excluded_site_sets = 0
    included_site_sets = 1

    for k in 1:min(order, N)
        visit_combinations!(N, k) do sites
            if site_diameter(sites) <= cutoff
                included_site_sets += 1
                push_axis_assignments!(basis, sites)
            else
                excluded_site_sets += 1
            end
        end
    end

    sort!(unique!(basis))
    return basis, included_site_sets, excluded_site_sets
end

function site_permutation_clifford(σx, σy, σz, perm::Vector{Int})
    M = typeof(σx[1])
    images = Dict{M,Tuple{Int,M}}()
    for i in 1:N
        j = perm[i]
        images[σx[i]] = (1, σx[j])
        images[σy[i]] = (1, σy[j])
        images[σz[i]] = (1, σz[j])
    end
    return CliffordSymmetry(images; nqubits=N)
end

function lattice_clifford_generators(σx, σy, σz)
    tx = [site(coord(s)[1] + 1, coord(s)[2]) for s in 1:N]
    ty = [site(coord(s)[1], coord(s)[2] + 1) for s in 1:N]
    rot90 = [site(coord(s)[2], Lx + 1 - coord(s)[1]) for s in 1:N]
    reflx = [site(Lx + 1 - coord(s)[1], coord(s)[2]) for s in 1:N]
    return CliffordSymmetry[
        site_permutation_clifford(σx, σy, σz, tx),
        site_permutation_clifford(σx, σy, σz, ty),
        site_permutation_clifford(σx, σy, σz, rot90),
        site_permutation_clifford(σx, σy, σz, reflx),
    ]
end

function global_spin_clifford_generators(σx, σy, σz)
    M = typeof(σx[1])

    h_images = Dict{M,Tuple{Int,M}}()
    s_images = Dict{M,Tuple{Int,M}}()
    for i in 1:N
        h_images[σx[i]] = (1, σz[i])
        h_images[σy[i]] = (-1, σy[i])
        h_images[σz[i]] = (1, σx[i])

        s_images[σx[i]] = (1, σy[i])
        s_images[σy[i]] = (-1, σx[i])
        # σz is fixed and may be omitted.
    end

    return CliffordSymmetry[
        CliffordSymmetry(h_images; nqubits=N),
        CliffordSymmetry(s_images; nqubits=N),
    ]
end

function symmetry_spec(σx, σy, σz)
    lattice = lattice_clifford_generators(σx, σy, σz)
    spin = global_spin_clifford_generators(σx, σy, σz)
    generators = if SYMMETRY_MODE == "lattice"
        lattice
    elseif SYMMETRY_MODE == "spin"
        spin
    elseif SYMMETRY_MODE == "full"
        [lattice; spin]
    else
        throw(ArgumentError("SYMMETRY_MODE must be lattice, spin, or full; got $SYMMETRY_MODE"))
    end
    return SymmetrySpec(generators; offblock_check=OFFBLOCK_CHECK), length(lattice), length(spin)
end

function print_block_summary(block_sizes)
    if isempty(block_sizes)
        println("  symmetry PSD block count = 0")
        return nothing
    end
    sorted_sizes = sort(collect(block_sizes); rev=true)
    println("  symmetry PSD block count = ", length(sorted_sizes))
    println("  symmetry PSD largest 20  = ", sorted_sizes[1:min(end, 20)])
    println("  symmetry PSD max/median  = ", maximum(sorted_sizes), " / ", sorted_sizes[cld(length(sorted_sizes), 2)])
    println("  symmetry PSD sum         = ", sum(sorted_sizes))
    PRINT_FULL_BLOCKS && println("  symmetry PSD blocks      = ", repr(block_sizes))
    return nothing
end

function fast_dense_custom_sparsity(pop, registry, basis::Vector{M}) where {M}
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) && isempty(pop.moment_eq_constraints) ||
        throw(ArgumentError("FAST_DENSE_SPARSITY is only valid for this unconstrained scratch experiment."))

    basis_set = Set(basis)
    missing_objective_terms = M[m for m in monomials(pop.objective) if !(m in basis_set)]
    isempty(missing_objective_terms) || throw(ArgumentError(
        "Custom basis must contain every Heisenberg objective monomial so identity/objective moments are generated; missing $(length(missing_objective_terms)) terms."
    ))

    T = eltype(indices(registry))
    P = typeof(pop.objective)
    corr = NCTSSoS.CorrelativeSparsity{PauliAlgebra,T,P,M,Nothing}(
        [collect(indices(registry))],
        registry,
        P[],
        [Int[]],
        Int[],
        [basis],
        [Vector{M}[]],
    )
    term_sparsities = Vector{NCTSSoS.TermSparsity{M}}[
        [NCTSSoS.TermSparsity{M}(M[], [basis])],
    ]
    return SparsityResult(corr, [M[]], term_sparsities)
end

function main()
    println("Started: ", Dates.now())
    println("Model: ", Lx, "x", Ly, " periodic Heisenberg, N=", N, ", bonds=", length(BONDS))
    println("Basis: Pauli strings degree <= ", ORDER, ", pairwise torus Manhattan site distance <= ", DISTANCE_CUTOFF)
    println("Symmetry mode: ", SYMMETRY_MODE, ", offblock_check=", OFFBLOCK_CHECK, ", solve=", SOLVE)
    println("Fast dense sparsity: ", FAST_DENSE_SPARSITY)
    println("Solve lowering: dualize=", DUALIZE, ", formulation=", FORMULATION, ", representation=", REPRESENTATION)
    flush(stdout)

    stage("build Hamiltonian")
    registry, (σx, σy, σz), H = heisenberg_pauli_hamiltonian()
    pop = polyopt(H, registry)
    println("  Pauli terms = ", length(terms(H)))
    println("  variables   = ", length(indices(registry)))

    stage("build custom basis")
    basis_time = @elapsed basis, included_site_sets, excluded_site_sets = distance_limited_pauli_basis(registry, ORDER, DISTANCE_CUTOFF)
    dense_formula = sum(binomial(N, k) * 3^k for k in 0:min(ORDER, N))
    println("  custom basis seconds       = ", round(basis_time; digits=3))
    println("  custom basis length        = ", length(basis))
    println("  dense Pauli order-$ORDER count = ", dense_formula)
    println("  kept/excluded site sets    = ", included_site_sets, " / ", excluded_site_sets)
    println("  degree histogram           = ", Dict(d => count(==(d), degree.(basis)) for d in 0:ORDER))
    flush(stdout)

    spec, n_lattice_gens, n_spin_gens = symmetry_spec(σx, σy, σz)
    println("  Clifford generators        = ", length(spec.clifford_generators), " (lattice ", n_lattice_gens, ", spin ", n_spin_gens, ")")

    optimizer = solver_optimizer()
    config = SolverConfig(
        optimizer=optimizer,
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=spec,
    )

    if FAST_DENSE_SPARSITY
        stage("fast dense custom sparsity start")
        sparsity_time = @elapsed sparsity = fast_dense_custom_sparsity(pop, registry, basis)
        stage("fast dense custom sparsity done seconds=$(round(sparsity_time; digits=3))")
    else
        stage("compute_sparsity start")
        sparsity_time = @elapsed sparsity = compute_sparsity(pop, config)
        stage("compute_sparsity done seconds=$(round(sparsity_time; digits=3))")
    end
    println("  clique count               = ", length(sparsity.corr_sparsity.cliques))
    println("  clique variable sizes      = ", length.(sparsity.corr_sparsity.cliques))
    println("  pre-symmetry matrix sizes  = ", NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities))
    flush(stdout)

    stage("symmetry support check start")
    support_time = @elapsed NCTSSoS._check_symmetry_mvp_support(pop, config, sparsity)
    stage("symmetry support check done seconds=$(round(support_time; digits=3))")

    stage("moment_relax_symmetric start")
    relax_time = @elapsed moment_problem, report = NCTSSoS.moment_relax_symmetric(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
        spec,
    )
    stage("moment_relax_symmetric done seconds=$(round(relax_time; digits=3))")
    println("  symmetry group order       = ", report.group_order)
    println("  invariant moments          = ", report.invariant_moment_count)
    println("  reduced total basis length = ", length(moment_problem.total_basis))
    println("  SDP constraints            = ", length(moment_problem.constraints))
    print_block_summary(report.psd_block_sizes)
    flush(stdout)

    if SOLVE
        stage("solve_sdp start")
        solve_time = @elapsed sdp_result = NCTSSoS.solve_sdp(
            moment_problem,
            config.optimizer;
            dualize=DUALIZE,
            formulation=FORMULATION,
            representation=REPRESENTATION,
        )
        stage("solve_sdp done seconds=$(round(solve_time; digits=3))")
        println("  solver status              = ", sdp_result.status)
        @printf("  objective lower bound      = %.12f\n", sdp_result.objective)
        @printf("  objective per site         = %.12f\n", sdp_result.objective / N)
        println("  unique SDP elements        = ", sdp_result.n_unique_elements)
    else
        println("Skipping solve because SOLVE=0.")
    end

    println("Finished: ", Dates.now())
    return nothing
end

main()
