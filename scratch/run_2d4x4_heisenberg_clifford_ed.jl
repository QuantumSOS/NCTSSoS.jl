#!/usr/bin/env julia

# 4x4 periodic spin-1/2 Heisenberg model:
#   H = sum_<ij> S_i · S_j = 1/4 sum_<ij>,a σᵃ_i σᵃ_j
# Compare an order-1 NCTSSoS lower bound using manual Clifford symmetry against
# exact diagonalization in the Sz=0 sector.

using LinearAlgebra
using Random
using Printf
using Dates
using JuMP
using NCTSSoS
import Clarabel

const HAS_MOSEK = try
    @eval import MosekTools
    true
catch err
    @warn "MosekTools unavailable; SOLVER=mosek will fail unless the wrapper installs it" exception=(err, catch_backtrace())
    false
end

const HAS_COSMO = try
    @eval import COSMO
    true
catch err
    @warn "COSMO unavailable; falling back to Clarabel" exception=(err, catch_backtrace())
    false
end

function solver_optimizer()
    solver = lowercase(get(ENV, "SOLVER", HAS_COSMO ? "cosmo" : "clarabel"))
    if solver == "cosmo"
        HAS_COSMO || throw(ArgumentError("SOLVER=cosmo requested but COSMO is unavailable."))
        return optimizer_with_attributes(
            COSMO.Optimizer,
            "verbose" => false,
            "eps_abs" => 1e-7,
            "eps_rel" => 1e-7,
            "eps_prim_inf" => 1e-6,
            "eps_dual_inf" => 1e-6,
            "max_iter" => 200_000,
            "rho" => 1.0,
            "adaptive_rho" => true,
            "alpha" => 1.0,
            "scaling" => 10,
        )
    elseif solver in ("cosmo_cg", "cosmo-cg")
        HAS_COSMO || throw(ArgumentError("SOLVER=cosmo_cg requested but COSMO is unavailable."))
        isdefined(COSMO, :CGIndirectKKTSolver) || throw(ArgumentError(
            "SOLVER=cosmo_cg requires COSMO.CGIndirectKKTSolver. Import/install IterativeSolvers and LinearMaps before loading COSMO."
        ))
        return optimizer_with_attributes(
            COSMO.Optimizer,
            "verbose" => true,
            "eps_abs" => 1e-6,
            "eps_rel" => 1e-6,
            "eps_prim_inf" => 1e-6,
            "eps_dual_inf" => 1e-6,
            "max_iter" => parse(Int, get(ENV, "COSMO_MAX_ITER", "50000")),
            "rho" => 1.0,
            "adaptive_rho" => true,
            "alpha" => 1.0,
            "scaling" => 10,
            "kkt_solver" => COSMO.CGIndirectKKTSolver,
        )
    elseif solver == "clarabel"
        return optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    elseif solver == "mosek"
        HAS_MOSEK || throw(ArgumentError("SOLVER=mosek requested but MosekTools is unavailable."))
        return MosekTools.Optimizer
    else
        throw(ArgumentError("SOLVER must be cosmo, clarabel, or mosek; got $solver"))
    end
end

const Lx = 4
const Ly = 4
const N = Lx * Ly
const LI = LinearIndices((1:Lx, 1:Ly))
const RUN_DENSE = get(ENV, "RUN_DENSE", "1") != "0"
const SYMMETRY_MODE = get(ENV, "SYMMETRY_MODE", "lattice")  # none, tx, translations, lattice, spin, full
const RELAX_ORDER = parse(Int, get(ENV, "RELAX_ORDER", "1"))
const OFFBLOCK_CHECK = Symbol(get(ENV, "OFFBLOCK_CHECK", "randomized"))
const DUALIZE = lowercase(get(ENV, "DUALIZE", "1")) ∉ ("0", "false", "no")
const FORMULATION = Symbol(get(ENV, "FORMULATION", "moment_variables"))
const REPRESENTATION = Symbol(get(ENV, "REPRESENTATION", "real"))
const ORPHAN_POLICY = Symbol(get(ENV, "ORPHAN_POLICY", "error"))
const SILENT_SOLVER = lowercase(get(ENV, "SILENT_SOLVER", "1")) ∉ ("0", "false", "no")

site(i, j) = LI[CartesianIndex(mod1(i, Lx), mod1(j, Ly))]
coord(s) = Tuple(CartesianIndices((1:Lx, 1:Ly))[s])

function square_lattice_bonds()
    bonds = Tuple{Int,Int}[]
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
    basis1 = [one(σx[1]); σx; σy; σz]
    return registry, (σx, σy, σz), H, basis1
end

function site_permutation_clifford(σx, σy, σz, perm::Vector{Int}; name="site permutation")
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
        site_permutation_clifford(σx, σy, σz, tx; name="tx"),
        site_permutation_clifford(σx, σy, σz, ty; name="ty"),
        site_permutation_clifford(σx, σy, σz, rot90; name="rot90"),
        site_permutation_clifford(σx, σy, σz, reflx; name="reflx"),
    ]
end

function global_spin_clifford_generators(σx, σy, σz)
    M = typeof(σx[1])

    h_images = Dict{M,Tuple{Int,M}}()
    for i in 1:N
        h_images[σx[i]] = (1, σz[i])
        h_images[σy[i]] = (-1, σy[i])
        h_images[σz[i]] = (1, σx[i])
    end

    s_images = Dict{M,Tuple{Int,M}}()
    for i in 1:N
        s_images[σx[i]] = (1, σy[i])
        s_images[σy[i]] = (-1, σx[i])
        # σz is fixed and may be omitted.
    end

    return CliffordSymmetry[
        CliffordSymmetry(h_images; nqubits=N),
        CliffordSymmetry(s_images; nqubits=N),
    ]
end

function bit_at(x::UInt32, pos::Int)
    return (x >> (pos - 1)) & UInt32(1)
end

function sz0_basis(n::Int)
    target = n ÷ 2
    states = UInt32[]
    sizehint!(states, binomial(n, target))
    for s in UInt32(0):(UInt32(1) << n) - UInt32(1)
        count_ones(s) == target && push!(states, s)
    end
    index = zeros(Int32, 1 << n)
    for (k, s) in enumerate(states)
        index[Int(s) + 1] = Int32(k)
    end
    return states, index
end

function heisenberg_matvec_sz0!(y::Vector{Float64}, x::Vector{Float64}, states::Vector{UInt32}, index::Vector{Int32})
    fill!(y, 0.0)
    @inbounds for col in eachindex(states)
        state = states[col]
        amp = x[col]
        diag = 0.0
        for (i, j) in BONDS
            bi = bit_at(state, i)
            bj = bit_at(state, j)
            if bi == bj
                diag += 0.25
            else
                diag -= 0.25
                flipped = state ⊻ (UInt32(1) << (i - 1)) ⊻ (UInt32(1) << (j - 1))
                row = Int(index[Int(flipped) + 1])
                y[row] += 0.5 * amp
            end
        end
        y[col] += diag * amp
    end
    return y
end

function lanczos_lowest(matvec!, dim::Int; maxiter::Int=260, tol::Float64=1e-11, seed::Int=4)
    Random.seed!(seed)
    Q = Matrix{Float64}(undef, dim, maxiter)
    q = randn(dim)
    q ./= norm(q)
    qprev = zeros(dim)
    w = zeros(dim)
    alphas = Float64[]
    betas = Float64[]
    prev_ritz = Inf
    best = Inf

    for k in 1:maxiter
        Q[:, k] .= q
        matvec!(w, q)
        if k > 1
            @. w -= betas[end] * qprev
        end
        α = dot(q, w)
        push!(alphas, α)
        @. w -= α * q

        # Full reorthogonalization. Wasteful? A little. Wrong ED is worse.
        for j in 1:k
            c = dot(view(Q, :, j), w)
            @. w -= c * Q[:, j]
        end

        β = norm(w)
        T = SymTridiagonal(alphas, betas)
        F = eigen(T)
        idx = argmin(F.values)
        ritz = F.values[idx]
        residual = β * abs(F.vectors[end, idx])
        best = ritz

        @printf("ED Lanczos iter %3d: E = %.12f, Δ = %.3e, residual = %.3e\n",
                k, ritz, abs(ritz - prev_ritz), residual)
        if k > 10 && residual < tol && abs(ritz - prev_ritz) < 10tol
            return (energy=ritz, iterations=k, residual=residual)
        end
        k == maxiter && return (energy=best, iterations=k, residual=residual)

        push!(betas, β)
        qprev, q = q, w ./ β
        prev_ritz = ritz
    end
end

function run_ed()
    println("\n== Exact diagonalization: Sz=0 matrix-free Lanczos ==")
    states, index = sz0_basis(N)
    println("Hilbert dimension in Sz=0 sector: ", length(states), " / ", 1 << N)
    y = zeros(Float64, length(states))
    matvec! = (out, vec) -> heisenberg_matvec_sz0!(out, vec, states, index)
    elapsed = @elapsed result = lanczos_lowest(matvec!, length(states))
    return merge(result, (; time=elapsed, dimension=length(states)))
end

function stage(msg)
    println("[", Dates.now(), "] ", msg)
    flush(stdout)
end

function run_sdp_case(label, pop, basis, optimizer; symmetry=nothing, order::Int=RELAX_ORDER)
    config = if isnothing(basis)
        SolverConfig(
            optimizer=optimizer,
            order=order,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=symmetry,
        )
    else
        SolverConfig(
            optimizer=optimizer,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=symmetry,
        )
    end

    report = nothing
    moment_sizes = Vector{Vector{Int}}()
    sdp_result = nothing
    timings = Dict{Symbol,Float64}()

    elapsed = @elapsed begin
        stage("compute_sparsity start")
        timings[:compute_sparsity] = @elapsed sparsity = compute_sparsity(pop, config)
        moment_sizes = NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)
        println("  compute_sparsity seconds = ", timings[:compute_sparsity])
        println("  pre-symmetry moment matrix sizes = ", moment_sizes)
        println("  clique count = ", length(sparsity.corr_sparsity.cliques))
        println("  clique variable sizes = ", length.(sparsity.corr_sparsity.cliques))
        flush(stdout)

        stage("symmetry support check start")
        timings[:support_check] = @elapsed NCTSSoS._check_symmetry_mvp_support(pop, config, sparsity)
        println("  support_check seconds = ", timings[:support_check])
        flush(stdout)

        if isnothing(symmetry)
            stage("moment_relax start")
            timings[:moment_relax] = @elapsed moment_problem = NCTSSoS.moment_relax(
                pop,
                sparsity.corr_sparsity,
                sparsity.cliques_term_sparsities,
            )
            println("  moment_relax seconds = ", timings[:moment_relax])
            flush(stdout)
        else
            stage("moment_relax_symmetric start")
            timings[:moment_relax_symmetric] = @elapsed begin
                moment_problem, report = NCTSSoS.moment_relax_symmetric(
                    pop,
                    sparsity.corr_sparsity,
                    sparsity.cliques_term_sparsities,
                    symmetry,
                )
            end
            println("  moment_relax_symmetric seconds = ", timings[:moment_relax_symmetric])
            println("  symmetry group order = ", report.group_order)
            println("  symmetry PSD blocks = ", report.psd_block_sizes)
            println("  invariant moments = ", report.invariant_moment_count)
            flush(stdout)
        end

        stage("solve_sdp start")
        println("  dualize = ", DUALIZE)
        println("  formulation = ", FORMULATION)
        println("  representation = ", REPRESENTATION)
        println("  orphan_policy = ", ORPHAN_POLICY)
        flush(stdout)
        timings[:solve_sdp] = @elapsed sdp_result = if DUALIZE
            NCTSSoS.solve_sdp(
                moment_problem,
                config.optimizer;
                dualize=true,
            )
        else
            direct_result = NCTSSoS.solve_moment_problem(
                moment_problem,
                config.optimizer;
                silent=SILENT_SOLVER,
                formulation=FORMULATION,
                representation=REPRESENTATION,
                orphan_policy=ORPHAN_POLICY,
            )
            (
                objective=direct_result.objective,
                model=direct_result.model,
                n_unique_elements=direct_result.n_unique_elements,
                status=NCTSSoS._check_solver_status(direct_result.model),
            )
        end
        println("  solve_sdp seconds = ", timings[:solve_sdp])
        println("  solver status = ", sdp_result.status)
        flush(stdout)
    end

    return (
        label=label,
        objective=sdp_result.objective,
        time=elapsed,
        stage_timings=timings,
        moment_matrix_sizes=moment_sizes,
        unique_moments=sdp_result.n_unique_elements,
        symmetry_report=report,
        group_order=isnothing(report) ? missing : report.group_order,
        psd_block_sizes=isnothing(report) ? missing : report.psd_block_sizes,
        generators=isnothing(symmetry) ? 0 : length(symmetry.clifford_generators),
    )
end

function print_sdp_result(result, ed_energy)
    println("\n== ", result.label, " ==")
    @printf("objective lower bound = %.12f (per site %.12f)\n", result.objective, result.objective / N)
    @printf("wall time             = %.3f s\n", result.time)
    println("pre-symmetry matrix  = ", result.moment_matrix_sizes)
    println("unique moments       = ", result.unique_moments)
    if result.symmetry_report !== nothing
        println("generators           = ", result.generators)
        println("group order          = ", result.group_order)
        println("PSD blocks           = ", result.psd_block_sizes)
        println("invariant moments    = ", result.symmetry_report.invariant_moment_count)
    end
    @printf("ED - bound           = %.12f\n", ed_energy - result.objective)
end

function main()
    println("Started: ", Dates.now())
    println("Model: 4x4 periodic Heisenberg, N=", N, ", bonds=", length(BONDS), ", relaxation order=", RELAX_ORDER, ", offblock_check=", OFFBLOCK_CHECK)
    println("SDP options: dualize=", DUALIZE, ", formulation=", FORMULATION, ", representation=", REPRESENTATION, ", orphan_policy=", ORPHAN_POLICY, ", silent_solver=", SILENT_SOLVER)

    ed = run_ed()
    @printf("\nED ground energy      = %.12f\n", ed.energy)
    @printf("ED per site           = %.12f\n", ed.energy / N)
    println("ED iterations         = ", ed.iterations)
    @printf("ED residual estimate  = %.3e\n", ed.residual)
    @printf("ED wall time          = %.3f s\n", ed.time)

    println("\n== Build NCTSSoS Pauli Hamiltonian ==")
    registry, (σx, σy, σz), H, basis1 = heisenberg_pauli_hamiltonian()
    pop = polyopt(H, registry)
    println("Pauli terms           = ", length(terms(H)))
    println("order-1 basis size    = ", length(basis1))
    sdp_basis = RELAX_ORDER == 1 ? basis1 : nothing
    isnothing(sdp_basis) && println("using automatic dense order-", RELAX_ORDER, " moment basis")

    optimizer = solver_optimizer()

    results = Any[]
    if RUN_DENSE
        println("\nAbout to solve dense order-", RELAX_ORDER, " baseline...")
        dense = run_sdp_case("NCTSSoS order-$RELAX_ORDER dense baseline", pop, sdp_basis, optimizer; order=RELAX_ORDER)
        print_sdp_result(dense, ed.energy)
        push!(results, dense)
    else
        println("\nSkipping dense baseline because RUN_DENSE=0.")
    end

    if SYMMETRY_MODE == "none"
        println("\nSkipping symmetry solve because SYMMETRY_MODE=none.")
        println("\n== Comparison ==")
        println("case,bound,bound_per_site,ed_minus_bound,time_seconds,group_order,generators,psd_blocks,unique_moments")
        for r in results
            @printf("%s,%.12f,%.12f,%.12f,%.3f,%s,%d,%s,%d\n",
                    r.label,
                    r.objective,
                    r.objective / N,
                    ed.energy - r.objective,
                    r.time,
                    string(r.group_order),
                    r.generators,
                    repr(r.psd_block_sizes),
                    r.unique_moments)
        end
        println("\nFinished: ", Dates.now())
        return nothing
    end

    lattice_gens = lattice_clifford_generators(σx, σy, σz)
    spin_gens = global_spin_clifford_generators(σx, σy, σz)
    selected_gens, label = if SYMMETRY_MODE == "tx"
        (lattice_gens[1:1], "NCTSSoS order-$RELAX_ORDER Clifford symmetry (x-translation subgroup)")
    elseif SYMMETRY_MODE == "translations"
        (lattice_gens[1:2], "NCTSSoS order-$RELAX_ORDER Clifford symmetry (2D translation subgroup)")
    elseif SYMMETRY_MODE == "lattice"
        (lattice_gens, "NCTSSoS order-$RELAX_ORDER Clifford symmetry (4x4 space-group site permutations)")
    elseif SYMMETRY_MODE == "spin"
        (spin_gens, "NCTSSoS order-$RELAX_ORDER Clifford symmetry (global spin Clifford only)")
    elseif SYMMETRY_MODE == "full"
        ([lattice_gens; spin_gens], "NCTSSoS order-$RELAX_ORDER Clifford symmetry (space group + global spin Clifford)")
    else
        throw(ArgumentError("SYMMETRY_MODE must be none, tx, translations, lattice, spin, or full; got $SYMMETRY_MODE"))
    end

    println("\nAbout to solve symmetry mode: ", SYMMETRY_MODE, " with ", length(selected_gens), " generators...")
    flush(stdout)
    spec = SymmetrySpec(selected_gens; offblock_check=OFFBLOCK_CHECK)
    symmetric = run_sdp_case(label, pop, sdp_basis, optimizer; symmetry=spec, order=RELAX_ORDER)
    print_sdp_result(symmetric, ed.energy)
    push!(results, symmetric)

    println("\n== Comparison ==")
    println("case,bound,bound_per_site,ed_minus_bound,time_seconds,group_order,generators,psd_blocks,unique_moments")
    for r in results
        @printf("%s,%.12f,%.12f,%.12f,%.3f,%s,%d,%s,%d\n",
                r.label,
                r.objective,
                r.objective / N,
                ed.energy - r.objective,
                r.time,
                string(r.group_order),
                r.generators,
                repr(r.psd_block_sizes),
                r.unique_moments)
    end

    println("\nFinished: ", Dates.now())
end

main()
