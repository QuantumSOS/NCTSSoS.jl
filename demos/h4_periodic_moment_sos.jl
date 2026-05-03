#!/usr/bin/env julia
# H4/Nk=2 periodic 1D+PQG V2RDM in NCTSSoS.jl language.
#
# This script formulates exactly the SDP we want:
#   minimize the explicit H4 active-space Hamiltonian over one fermionic CAR
#   moment functional, with explicit 1D, D, Q, and G Hermitian PSD blocks.
#
# No SDP solver is called here. The useful artifact is the symbolic
# NCTSSoS MomentProblem.
#
# Usage from the repository root:
#   julia --project=. demos/h4_periodic_moment_sos.jl
#   julia --project=. demos/h4_periodic_moment_sos.jl --integrals=output/h4_chain_nk2_figure1_integrals_drop1e-12.txt
#
# Caveat: the repository's default text asset is the older H4 dump. For the
# Figure-1-matching 1-D PySCF Hamiltonian, use a text dump converted from
# /Users/exaclior/QuantumSOS/replicate/data/h4_chain_nk2_figure1.npz.
#
# The integrals text format is the PySCF dump used in the replication notes:
#   h1e k p q Re Im
#   eri k1 k2 k3 k4 p r q s Re Im
# with ERIs in chemist order (p r | q s). This script normalizes them by
# 1/Nk and 1/Nk^2 before building the per-cell Hamiltonian.

using NCTSSoS
using Printf

const SPINS = (:up, :dn)
const DEFAULT_INTEGRALS = normpath(joinpath(@__DIR__, "..", "test", "data", "assets", "h4_chain_nk2_integrals.txt"))

struct Options
    integrals_path::String
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
end

function parse_options(argv)
    integrals_path = DEFAULT_INTEGRALS
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = false
    singlet_s2 = false

    for arg in argv
        if startswith(arg, "--integrals=")
            integrals_path = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--blocking=")
            value = split(arg, "=", limit = 2)[2]
            if value in ("momentum", "k")
                blocking = :momentum
            elseif value == "spin"
                blocking = :spin
            elseif value == "none"
                blocking = :none
            else
                throw(ArgumentError("unknown --blocking=$value; expected momentum, spin, or none"))
            end
        elseif arg == "--include-1d"
            include_one_d = true
        elseif arg == "--no-1d"
            include_one_d = false
        elseif arg in ("--paper-spin", "--spin-singlet")
            spin_resolved_trace = true
            singlet_s2 = true
        elseif arg == "--spin-resolved-trace"
            spin_resolved_trace = true
        elseif arg == "--no-spin-resolved-trace"
            spin_resolved_trace = false
        elseif arg == "--singlet-s2"
            singlet_s2 = true
        elseif arg == "--no-singlet-s2"
            singlet_s2 = false
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. demos/h4_periodic_moment_sos.jl [options]\n")
            println("Options:")
            println("  --integrals=PATH              text integral dump; default test/data/assets/h4_chain_nk2_integrals.txt")
            println("  --blocking=momentum|spin|none 1D/PQG block grouping; default momentum")
            println("  --include-1d                  add an extra ¹D PSD block (off by default; matches")
            println("                                replicate/scripts/nctssos_pqg_regen default)")
            println("  --no-1d                       explicitly drop the ¹D PSD block (default)")
            println("  --paper-spin, --spin-singlet  add spin-resolved ²D traces and singlet S²")
            println("                                constraint; combine with --include-1d for the")
            println("                                Figure-1 comparison SDP")
            println("  --spin-resolved-trace         replace total TrD with TrDαα/TrDαβ/TrDββ")
            println("  --no-spin-resolved-trace      use one total TrD constraint (default)")
            println("  --singlet-s2                  add Σᵢⱼ D(iα,jβ;jα,iβ) = N/2")
            println("  --no-singlet-s2               do not add singlet S² constraint (default)")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    return Options(integrals_path, blocking, include_one_d, spin_resolved_trace, singlet_s2)
end

bytes_to_gib(bytes::Integer) = bytes / 1024.0^3
rss_string() = @sprintf("%.3f GiB", bytes_to_gib(Sys.maxrss()))

@inline complex_entry(fields, re_idx::Int, im_idx::Int) =
    complex(parse(Float64, fields[re_idx]), parse(Float64, fields[im_idx]))

function load_integrals_txt(path::AbstractString; norb::Int)
    h1e = Dict{Int,Matrix{ComplexF64}}()
    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()

    for (line_no, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        (isempty(line) || startswith(line, '#')) && continue

        fields = split(line)
        tag = first(fields)
        if tag == "h1e"
            length(fields) == 6 || error("malformed h1e line $line_no in $(repr(path))")
            k = parse(Int, fields[2])
            p = parse(Int, fields[3]) + 1
            q = parse(Int, fields[4]) + 1
            value = complex_entry(fields, 5, 6)
            block = get!(h1e, k) do
                zeros(ComplexF64, norb, norb)
            end
            block[p, q] = value
        elseif tag == "eri"
            length(fields) == 11 || error("malformed eri line $line_no in $(repr(path))")
            k1 = parse(Int, fields[2])
            k2 = parse(Int, fields[3])
            k3 = parse(Int, fields[4])
            k4 = parse(Int, fields[5])

            # Text dump stores chemist order (p k1, r k3 | q k2, s k4).
            p = parse(Int, fields[6]) + 1
            r = parse(Int, fields[7]) + 1
            q = parse(Int, fields[8]) + 1
            s = parse(Int, fields[9]) + 1
            value = complex_entry(fields, 10, 11)

            block = get!(eri, (k1, k2, k3, k4)) do
                zeros(ComplexF64, norb, norb, norb, norb)
            end
            block[p, r, q, s] = value
        else
            error("unexpected integral tag $(repr(tag)) at line $line_no in $(repr(path))")
        end
    end

    return h1e, eri
end

function active_hf_energy(h1e, eri; nk::Int, nelec_per_cell::Int)
    nocc = nelec_per_cell ÷ 2

    e1 = 0.0 + 0.0im
    for k in 0:(nk - 1), i in 1:nocc
        e1 += 2.0 * h1e[k][i, i] / nk
    end

    e2 = 0.0 + 0.0im
    for ((k1, k2, k3, k4), Vraw) in eri
        V = Vraw / nk^2
        for i in 1:nocc, j in 1:nocc
            # Closed-shell HF contractions are k-diagonal. Umklapp blocks do
            # not enter the single-determinant density contraction.
            if k1 == k3 && k2 == k4
                e2 += 2.0 * V[i, i, j, j]
            end
            if k1 == k4 && k2 == k3
                e2 -= V[i, j, j, i]
            end
        end
    end

    return real(e1 + e2)
end

function only_monomial(poly)
    mono = monomials(poly)
    length(mono) == 1 || error("expected a single monomial, got $(length(mono)) terms")
    return only(mono)
end

function build_h4_nk2_hamiltonian(h1e, eri; nk::Int, norb::Int)
    registry,
    ((c_up_k0, c_up_k0_dag),
     (c_dn_k0, c_dn_k0_dag),
     (c_up_k1, c_up_k1_dag),
     (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
        ("c_up_k0", 1:norb),
        ("c_dn_k0", 1:norb),
        ("c_up_k1", 1:norb),
        ("c_dn_k1", 1:norb),
    ])

    ann = Dict(
        (0, :up) => c_up_k0,
        (0, :dn) => c_dn_k0,
        (1, :up) => c_up_k1,
        (1, :dn) => c_dn_k1,
    )
    dag = Dict(
        (0, :up) => c_up_k0_dag,
        (0, :dn) => c_dn_k0_dag,
        (1, :up) => c_up_k1_dag,
        (1, :dn) => c_dn_k1_dag,
    )

    ham = (0.0 + 0.0im) * (c_up_k0_dag[1] * c_up_k0[1])

    # One-electron part. Raw h1e is per supercell; divide by Nk for per-cell energy.
    for k in 0:(nk - 1)
        hk = h1e[k] / nk
        for p in 1:norb, q in 1:norb, spin in SPINS
            coeff = hk[p, q]
            iszero(coeff) && continue
            ham += coeff * dag[(k, spin)][p] * ann[(k, spin)][q]
        end
    end

    # Two-electron part. Raw ERI is chemist order (p r | q s), so physicist
    # <p q | r s> is V[p, r, q, s]. Operator order is a†p a†q a_s a_r.
    spin_channels = ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))
    for key in sort!(collect(keys(eri)))
        k1, k2, k3, k4 = key
        V = eri[key] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff = V[p, r, q, s]
            iszero(coeff) && continue
            for (σ, τ) in spin_channels
                ham += 0.5 * coeff *
                    dag[(k1, σ)][p] * dag[(k2, τ)][q] *
                    ann[(k4, τ)][s] * ann[(k3, σ)][r]
            end
        end
    end

    # Remove numerical anti-Hermitian crumbs from integral noise.
    ham = 0.5 * (ham + adjoint(ham))

    vars = (; ann, dag,
        c_up_k0, c_up_k0_dag, c_dn_k0, c_dn_k0_dag,
        c_up_k1, c_up_k1_dag, c_dn_k1, c_dn_k1_dag)
    return registry, vars, ham
end

function build_spin_orbitals(vars; nk::Int, norb::Int)
    spin_orbitals = NamedTuple[]
    for k in 0:(nk - 1), spin in SPINS, orb in 1:norb
        push!(spin_orbitals, (; a = vars.ann[(k, spin)][orb],
                              adag = vars.dag[(k, spin)][orb],
                              k,
                              spin,
                              orb))
    end
    return spin_orbitals
end

function basis_label(m::NormalMonomial; nk::Int, norb::Int)
    ΔN = 0
    K = 0
    two_Sz = 0

    for idx in m.word
        raw = Int(idx)
        mode = abs(raw)
        family, _ = divrem(mode - 1, norb)
        k = family ÷ 2
        spin_sign = iseven(family) ? 1 : -1
        particle_sign = raw < 0 ? 1 : -1  # creator -> +1, annihilator -> -1

        ΔN += particle_sign
        K = mod(K + particle_sign * k, nk)
        two_Sz += particle_sign * spin_sign
    end

    return (; ΔN, K, two_Sz)
end

function block_key(mono::NormalMonomial, blocking::Symbol; nk::Int, norb::Int)
    label = basis_label(mono; nk, norb)
    return blocking == :momentum ? (label.K,) : (label.K, label.two_Sz)
end

function grouped_blocks(basis::Vector{NM}, blocking::Symbol; nk::Int, norb::Int) where {NM<:NormalMonomial}
    if blocking == :none
        return Vector{NM}[basis]
    end

    buckets = Dict{Tuple,Vector{NM}}()
    for mono in basis
        key = block_key(mono, blocking; nk, norb)
        push!(get!(buckets, key, NM[]), mono)
    end

    return [buckets[key] for key in sort!(collect(keys(buckets)))]
end

function grouped_blocks_with_identity(
    basis::Vector{NM},
    identity::NM,
    blocking::Symbol;
    nk::Int,
    norb::Int,
) where {NM<:NormalMonomial}
    if blocking == :none
        return Vector{NM}[NM[identity; basis]]
    end

    buckets = Dict{Tuple,Vector{NM}}()
    for mono in basis
        key = block_key(mono, blocking; nk, norb)
        push!(get!(buckets, key, NM[]), mono)
    end

    identity_key = blocking == :momentum ? (0,) : (0, 0)
    pushfirst!(get!(buckets, identity_key, NM[]), identity)
    return [buckets[key] for key in sort!(collect(keys(buckets)))]
end

function build_rdm_bases(spin_orbitals)
    NM = typeof(only_monomial(spin_orbitals[1].adag * spin_orbitals[1].a))
    M = length(spin_orbitals)

    oneD_basis = NM[]
    twoD_basis = NM[]
    twoQ_basis = NM[]
    twoG_basis = NM[]

    for p in 1:M
        push!(oneD_basis, only_monomial(spin_orbitals[p].a))
    end

    for p in 1:(M - 1), q in (p + 1):M
        push!(twoD_basis, only_monomial(spin_orbitals[p].a * spin_orbitals[q].a))
        push!(twoQ_basis, only_monomial(spin_orbitals[p].adag * spin_orbitals[q].adag))
    end

    for p in 1:M, q in 1:M
        push!(twoG_basis, only_monomial(spin_orbitals[p].adag * spin_orbitals[q].a))
    end

    return (; oneD_basis, twoD_basis, twoQ_basis, twoG_basis)
end

function total_electron_constraint(vars, ham; norb::Int, total_electrons::Int)
    n_up = 1.0 * sum(
        vars.c_up_k0_dag[i] * vars.c_up_k0[i] + vars.c_up_k1_dag[i] * vars.c_up_k1[i]
        for i in 1:norb
    )
    n_dn = 1.0 * sum(
        vars.c_dn_k0_dag[i] * vars.c_dn_k0[i] + vars.c_dn_k1_dag[i] * vars.c_dn_k1[i]
        for i in 1:norb
    )

    return n_up + n_dn - Float64(total_electrons) * one(ham)
end

function trace_constraints(spin_orbitals, ham; total_electrons::Int, include_total_d::Bool = true)
    n_modes = length(spin_orbitals)
    n_holes = n_modes - total_electrons

    trD = zero(ham)
    trQ = zero(ham)
    trG = zero(ham)

    for p in 1:(n_modes - 1), q in (p + 1):n_modes
        drow = spin_orbitals[p].a * spin_orbitals[q].a
        qrow = spin_orbitals[p].adag * spin_orbitals[q].adag
        trD += adjoint(drow) * drow
        trQ += adjoint(qrow) * qrow
    end

    for p in 1:n_modes, q in 1:n_modes
        grow = spin_orbitals[p].adag * spin_orbitals[q].a
        trG += adjoint(grow) * grow
    end

    constraints = typeof(ham)[]
    if include_total_d
        push!(constraints,
            trD - Float64(total_electrons * (total_electrons - 1) ÷ 2) * one(ham))
    end
    push!(constraints,
        trQ - Float64(n_holes * (n_holes - 1) ÷ 2) * one(ham),
        trG - Float64(total_electrons * (n_holes + 1)) * one(ham),
    )
    return constraints
end

function spin_resolved_d_trace_constraints(spin_orbitals, ham; total_electrons::Int)
    iseven(total_electrons) ||
        throw(ArgumentError("spin-resolved trace constraints assume an M_s = 0 sector"))

    n_alpha = total_electrons ÷ 2
    n_beta = total_electrons ÷ 2
    tr_aa = zero(ham)
    tr_ab = zero(ham)
    tr_bb = zero(ham)

    for p in 1:(length(spin_orbitals) - 1), q in (p + 1):length(spin_orbitals)
        drow = spin_orbitals[p].a * spin_orbitals[q].a
        term = adjoint(drow) * drow
        spin_p = spin_orbitals[p].spin
        spin_q = spin_orbitals[q].spin
        if spin_p == :up && spin_q == :up
            tr_aa += term
        elseif spin_p == :dn && spin_q == :dn
            tr_bb += term
        else
            tr_ab += term
        end
    end

    return typeof(ham)[
        tr_aa - Float64(n_alpha * (n_alpha - 1) ÷ 2) * one(ham),
        tr_ab - Float64(n_alpha * n_beta) * one(ham),
        tr_bb - Float64(n_beta * (n_beta - 1) ÷ 2) * one(ham),
    ]
end

function singlet_s2_constraint(vars, ham; nk::Int, norb::Int, total_electrons::Int)
    exchange_sum = zero(ham)
    for ki in 0:(nk - 1), pi in 1:norb, kj in 0:(nk - 1), pj in 1:norb
        exchange_sum +=
            vars.dag[(ki, :up)][pi] * vars.dag[(kj, :dn)][pj] *
            vars.ann[(ki, :dn)][pi] * vars.ann[(kj, :up)][pj]
    end
    exchange_sum = 0.5 * (exchange_sum + adjoint(exchange_sum))
    return exchange_sum - (0.5 * total_electrons) * one(ham)
end

function single_full_system_clique(registry, moment_support_basis::Vector{NM}, objective::P; nk::Int, norb::Int) where {NM<:NormalMonomial,P}
    T = eltype(one(NM).word)
    clique_modes = T.(1:(2 * nk * norb))

    # NCTSSoS.moment_relax wants a CorrelativeSparsity object. We are not
    # computing or using correlative sparsity here: one clique covers the full
    # active space, and the Schouten 1D/PQG blocks below carry the structure.
    return NCTSSoS.CorrelativeSparsity{FermionicAlgebra,T,P,NM,Nothing}(
        [clique_modes],
        registry,
        P[],
        [Int[]],
        Int[],
        [moment_support_basis],
        [Vector{Vector{NM}}()],
    )
end

function build_h4_pqg_moment_problem(options::Options)
    nk = 2
    norb = 8
    nelec_per_cell = 4
    total_electrons = nk * nelec_per_cell

    h1e, eri = load_integrals_txt(options.integrals_path; norb)
    registry, vars, h4_ham = build_h4_nk2_hamiltonian(h1e, eri; nk, norb)
    objective = h4_ham
    spin_orbitals = build_spin_orbitals(vars; nk, norb)

    meq_constraints = typeof(h4_ham)[
        total_electron_constraint(vars, h4_ham; norb, total_electrons),
    ]
    append!(meq_constraints, trace_constraints(spin_orbitals, h4_ham;
        total_electrons,
        include_total_d = !options.spin_resolved_trace,
    ))
    if options.spin_resolved_trace
        append!(meq_constraints,
            spin_resolved_d_trace_constraints(spin_orbitals, h4_ham; total_electrons))
    end
    if options.singlet_s2
        push!(meq_constraints,
            singlet_s2_constraint(vars, h4_ham; nk, norb, total_electrons))
    end

    pop = polyopt(objective, registry; moment_eq_constraints = meq_constraints)

    row_bases = build_rdm_bases(spin_orbitals)
    NM = eltype(row_bases.oneD_basis)
    identity = one(NM)

    oneD_blocks = options.include_one_d ?
        grouped_blocks(row_bases.oneD_basis, options.blocking; nk, norb) :
        Vector{NM}[]
    twoD_blocks = grouped_blocks(row_bases.twoD_basis, options.blocking; nk, norb)
    twoQ_blocks = grouped_blocks(row_bases.twoQ_basis, options.blocking; nk, norb)
    twoG_blocks = grouped_blocks_with_identity(row_bases.twoG_basis, identity, options.blocking; nk, norb)

    moment_support_basis = options.include_one_d ? NM[
        identity;
        row_bases.oneD_basis;
        row_bases.twoD_basis;
        row_bases.twoQ_basis;
        row_bases.twoG_basis;
    ] : NM[
        identity;
        row_bases.twoD_basis;
        row_bases.twoQ_basis;
        row_bases.twoG_basis;
    ]

    schouten_blocks = Vector{NM}[]
    append!(schouten_blocks, oneD_blocks)
    append!(schouten_blocks, twoD_blocks)
    append!(schouten_blocks, twoQ_blocks)
    append!(schouten_blocks, twoG_blocks)

    full_clique = single_full_system_clique(registry, moment_support_basis, objective; nk, norb)
    term_sparsity = NCTSSoS.TermSparsity{NM}(copy(moment_support_basis), schouten_blocks)
    moment_problem = NCTSSoS.moment_relax(pop, full_clique, [[term_sparsity]])
    hf_active = active_hf_energy(h1e, eri; nk, nelec_per_cell)

    return (; nk, norb, nelec_per_cell, total_electrons, h1e, eri,
             registry, vars, h4_ham, objective, pop, meq_constraints,
             spin_orbitals, row_bases, oneD_blocks, twoD_blocks, twoQ_blocks,
             twoG_blocks, moment_support_basis, schouten_blocks, term_sparsity,
             moment_problem, hf_active)
end

real_lift_triangle_rows(n::Integer) = (2n) * (2n + 1) ÷ 2

function constraint_stats(mp)
    hpsd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone == :HPSD]
    psd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone == :PSD]
    zero_sizes = [size(mat) for (cone, mat) in mp.constraints if cone == :Zero]
    total_canonical_moments = length(NCTSSoS._sorted_symmetric_basis(mp.total_basis))
    return (; hpsd_sizes, psd_sizes, zero_sizes,
             total_canonical_moments,
             direct_real_moment_variables = 2 * total_canonical_moments,
             real_lift_psd_rows = sum(real_lift_triangle_rows(n) for n in hpsd_sizes),
             complex_zero_entries = sum((prod(size) for size in zero_sizes); init = 0),
             direct_real_zero_rows = 2 * sum((prod(size) for size in zero_sizes); init = 0))
end

function print_summary(data, options::Options, build_seconds::Real)
    stats = constraint_stats(data.moment_problem)

    formulation_label = options.include_one_d ?
        "H4/Nk=2 periodic 1D+PQG V2RDM as an NCTSSoS MomentProblem" :
        "H4/Nk=2 periodic PQG V2RDM (no ¹D) as an NCTSSoS MomentProblem"
    println("== ", formulation_label, " ==")
    @printf("%-48s %s\n", "integrals", options.integrals_path)
    @printf("%-48s %s\n", "blocking", string(options.blocking))
    @printf("%-48s %s\n", "extra ¹D PSD block",
            options.include_one_d ? "included" : "disabled (paper PQG default)")
    @printf("%-48s %s\n", "objective", "H4 active-space Hamiltonian")
    @printf("%-48s %s\n", "particle constraint", "N=8 via moment_eq_constraints")
    trace_label = options.spin_resolved_trace ?
        "TrDαα=6, TrDαβ=16, TrDββ=6, TrQ=276, TrG=200" :
        "TrD=28, TrQ=276, TrG=200"
    @printf("%-48s %s\n", "trace constraints", trace_label)
    @printf("%-48s %s\n", "singlet S² constraint",
            options.singlet_s2 ? "Σᵢⱼ D(iα,jβ;jα,iβ)=4 included" : "disabled")
    @printf("%-48s %s\n", "correlative sparsity", "none; single full-system clique")
    println()

    @printf("%-48s %d\n", "k-points", data.nk)
    @printf("%-48s %d\n", "active spatial orbitals per k", data.norb)
    @printf("%-48s %d\n", "spin-orbital modes", length(data.spin_orbitals))
    @printf("%-48s %d\n", "active electrons total", data.total_electrons)
    @printf("%-48s %.12f Ha\n", "HF active energy reconstructed", data.hf_active)
    @printf("%-48s %d\n", "Hamiltonian monomials", length(monomials(data.h4_ham)))
    println()

    println("-- Moment support basis --")
    @printf("%-48s %d\n", "identity", 1)
    @printf("%-48s %s\n", "¹D rows a_p",
            options.include_one_d ? string(length(data.row_bases.oneD_basis)) : "not included")
    @printf("%-48s %d\n", "²D rows a_p a_q, p<q", length(data.row_bases.twoD_basis))
    @printf("%-48s %d\n", "²Q rows a†_p a†_q, p<q", length(data.row_bases.twoQ_basis))
    @printf("%-48s %d\n", "²G rows a†_p a_q", length(data.row_bases.twoG_basis))
    @printf("%-48s %d\n", "total support rows", length(data.moment_support_basis))
    println()

    println("-- Schouten HPSD blocks --")
    @printf("%-48s %s\n", "¹D",
            options.include_one_d ? string(length.(data.oneD_blocks)) : "not included")
    @printf("%-48s %s\n", "²D", string(length.(data.twoD_blocks)))
    @printf("%-48s %s\n", "²Q", string(length.(data.twoQ_blocks)))
    @printf("%-48s %s\n", "²G", string(length.(data.twoG_blocks)))
    @printf("%-48s %s\n", "identity", "inside ²G zero-momentum sector")
    @printf("%-48s %d\n", "total HPSD blocks", length(stats.hpsd_sizes))
    @printf("%-48s %d\n", "largest HPSD block", maximum(stats.hpsd_sizes))
    println()

    println("-- Symbolic moment problem --")
    @printf("%-48s %d\n", "unique moment-matrix monomials", data.moment_problem.n_unique_moment_matrix_elements)
    @printf("%-48s %d\n", "total canonical monomials incl. constraints", stats.total_canonical_moments)
    @printf("%-48s %d\n", "direct real moment variables", stats.direct_real_moment_variables)
    @printf("%-48s %d\n", "real-lift PSD scalar rows", stats.real_lift_psd_rows)
    @printf("%-48s %d\n", "moment_eq polynomials", length(data.meq_constraints))
    @printf("%-48s %d\n", "Zero matrix constraints", length(stats.zero_sizes))
    @printf("%-48s %d\n", "Zero complex scalar entries", stats.complex_zero_entries)
    @printf("%-48s %d\n", "direct real zero/equality rows", stats.direct_real_zero_rows)
    @printf("%-48s %d\n", "plus solve-time normalization rows", 2)
    @printf("%-48s %d\n", "direct real rows incl. normalization", stats.real_lift_psd_rows + stats.direct_real_zero_rows + 2)
    println()

    @printf("%-48s %.3f s\n", "build walltime", build_seconds)
    @printf("%-48s %s\n", "max RSS after build", rss_string())
    println()
    println("No optimizer was attached and optimize! was not called.")
    println("Artifact: symbolic 1D+PQG MomentProblem over one shared fermionic CAR moment map.")
    return nothing
end

function main(argv = ARGS)
    options = parse_options(argv)
    @printf("%-48s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h4_pqg_moment_problem(options)
    print_summary(data, options, build_seconds)
    return data
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
