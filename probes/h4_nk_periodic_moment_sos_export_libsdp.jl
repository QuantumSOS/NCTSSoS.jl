#!/usr/bin/env julia
# Build a Julia/NCTSSoS H4 periodic PQG MomentProblem for a chosen Nk and
# export it as a libsdp/BPSDP complex SDPA-sparse `.dat-c` file.
#
# This is intentionally a benchmark probe, not a chemistry asset generator.  If
# --integrals is omitted, it uses a deterministic Julia mock Hamiltonian with
# H4-like electron count and k-momentum structure.  No Python, no PySCF.

using LinearAlgebra
using NCTSSoS
using Printf

include(joinpath(@__DIR__, "..", "demos", "SDPALibsdpExport.jl"))
using .SDPALibsdpExport: export_libsdp

const SPINS = (:up, :dn)

struct Options
    nk::Int
    norb::Int
    nelec_per_cell::Int
    integrals_path::Union{Nothing,String}
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
    outdir::String
    basename::String
    sense::Symbol
    orphans_per_block::Int
end

function parse_options(argv)
    nk = 2
    norb = 4
    nelec_per_cell = 4
    integrals_path = nothing
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = true
    singlet_s2 = true
    outdir = normpath(joinpath(@__DIR__, "..", "output", "h4_nk_sweep"))
    basename = ""
    sense = :min
    orphans_per_block = 32

    for arg in argv
        if startswith(arg, "--nk=")
            nk = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--norb=")
            norb = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--nelec-per-cell=")
            nelec_per_cell = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--integrals=")
            value = String(split(arg, "=", limit = 2)[2])
            integrals_path = isempty(value) || value == "mock" ? nothing : value
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
        elseif startswith(arg, "--outdir=")
            outdir = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--basename=")
            basename = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--sense=")
            sense = Symbol(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--orphans-per-block=")
            orphans_per_block = parse(Int, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. probes/h4_nk_periodic_moment_sos_export_libsdp.jl [options]\n")
            println("Core options:")
            println("  --nk=N                       k-point count; default 2")
            println("  --norb=N                     active spatial orbitals per k; default 4 for the mock sweep")
            println("  --nelec-per-cell=N           active electrons per cell; default 4")
            println("  --integrals=PATH|mock        optional text integral dump; default mock")
            println("  --blocking=momentum|spin|none  PQG block grouping; default momentum")
            println("  --include-1d | --no-1d       include/drop 1D PSD blocks; default drop")
            println("  --spin-resolved-trace | --no-spin-resolved-trace")
            println("  --singlet-s2 | --no-singlet-s2")
            println("Export options:")
            println("  --outdir=PATH                output directory; default output/h4_nk_sweep")
            println("  --basename=STEM              output stem; default h4_nk<N>_mock_pqg")
            println("  --orphans-per-block=N        default 32")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    nk >= 1 || throw(ArgumentError("--nk must be positive"))
    norb >= 1 || throw(ArgumentError("--norb must be positive"))
    nelec_per_cell >= 1 || throw(ArgumentError("--nelec-per-cell must be positive"))
    sense in (:min, :max) || throw(ArgumentError("--sense must be min or max"))
    orphans_per_block >= 1 || throw(ArgumentError("--orphans-per-block must be >= 1"))
    if isempty(basename)
        source = integrals_path === nothing ? "mock" : "asset"
        basename = "h4_nk$(nk)_$(source)_pqg"
    end

    return Options(nk, norb, nelec_per_cell, integrals_path, blocking,
        include_one_d, spin_resolved_trace, singlet_s2, outdir, basename,
        sense, orphans_per_block)
end

bytes_to_gib(bytes::Integer) = bytes / 1024.0^3
rss_string() = @sprintf("%.3f GiB", bytes_to_gib(Sys.maxrss()))

@inline complex_entry(fields, re_idx::Int, im_idx::Int) =
    complex(parse(Float64, fields[re_idx]), parse(Float64, fields[im_idx]))

function load_integrals_txt(path::AbstractString; nk::Int, norb::Int)
    h1e = Dict(k => zeros(ComplexF64, norb, norb) for k in 0:(nk - 1))
    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()

    for (line_no, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        (isempty(line) || startswith(line, '#')) && continue

        fields = split(line)
        tag = first(fields)
        if tag == "h1e"
            length(fields) == 6 || error("malformed h1e line $line_no in $(repr(path))")
            k = parse(Int, fields[2])
            0 <= k < nk || error("h1e line $line_no has k=$k outside 0:$(nk - 1)")
            p = parse(Int, fields[3]) + 1
            q = parse(Int, fields[4]) + 1
            1 <= p <= norb && 1 <= q <= norb || error("h1e line $line_no has orbital outside 0:$(norb - 1)")
            h1e[k][p, q] = complex_entry(fields, 5, 6)
        elseif tag == "eri"
            length(fields) == 11 || error("malformed eri line $line_no in $(repr(path))")
            k1 = parse(Int, fields[2])
            k2 = parse(Int, fields[3])
            k3 = parse(Int, fields[4])
            k4 = parse(Int, fields[5])
            all(k -> 0 <= k < nk, (k1, k2, k3, k4)) ||
                error("eri line $line_no has k outside 0:$(nk - 1)")
            p = parse(Int, fields[6]) + 1
            r = parse(Int, fields[7]) + 1
            q = parse(Int, fields[8]) + 1
            s = parse(Int, fields[9]) + 1
            all(i -> 1 <= i <= norb, (p, r, q, s)) ||
                error("eri line $line_no has orbital outside 0:$(norb - 1)")
            block = get!(eri, (k1, k2, k3, k4)) do
                zeros(ComplexF64, norb, norb, norb, norb)
            end
            block[p, r, q, s] = complex_entry(fields, 10, 11)
        else
            error("unexpected integral tag $(repr(tag)) at line $line_no in $(repr(path))")
        end
    end

    return h1e, eri
end

function mock_h4_integrals(; nk::Int, norb::Int)
    # Deterministic, sparse-ish, momentum-conserving H4-like mock integrals.
    # Values are unnormalized; the Hamiltonian builder divides by Nk and Nk^2.
    h1e = Dict{Int,Matrix{ComplexF64}}()
    for k in 0:(nk - 1)
        θ = 2π * k / nk
        H = zeros(ComplexF64, norb, norb)
        for p in 1:norb
            H[p, p] = -1.30 - 0.11p + 0.035cos(θ + 0.37p)
        end
        for p in 1:(norb - 1)
            amp = -0.055 / (1 + 0.2p)
            val = amp * cis(θ * (p % 2 == 0 ? -1 : 1))
            H[p, p + 1] = val
            H[p + 1, p] = conj(val)
        end
        h1e[k] = H
    end

    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()
    for k1 in 0:(nk - 1), k2 in 0:(nk - 1), k3 in 0:(nk - 1)
        k4 = mod(k1 + k2 - k3, nk)
        V = zeros(ComplexF64, norb, norb, norb, norb)
        transfer = min(mod(k1 - k3, nk), mod(k3 - k1, nk))
        qscale = 1.0 / (1.0 + transfer)
        for p in 1:norb, q in 1:norb
            direct = 0.18 * qscale / (1.0 + abs(p - q))
            exchange = 0.045 * qscale / (1.0 + abs(p - q))
            # Chemist order V[p,r,q,s].
            V[p, p, q, q] += direct
            V[p, q, q, p] -= exchange
        end
        if any(!iszero, V)
            eri[(k1, k2, k3, k4)] = V
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

function create_periodic_fermions(; nk::Int, norb::Int)
    families = Tuple{String,UnitRange{Int}}[]
    for k in 0:(nk - 1)
        push!(families, ("c_up_k$(k)", 1:norb))
        push!(families, ("c_dn_k$(k)", 1:norb))
    end
    registry, groups = create_fermionic_variables(families)

    ann = Dict{Tuple{Int,Symbol},Vector}()
    dag = Dict{Tuple{Int,Symbol},Vector}()
    for k in 0:(nk - 1)
        up_ann, up_dag = groups[2k + 1]
        dn_ann, dn_dag = groups[2k + 2]
        ann[(k, :up)] = up_ann
        dag[(k, :up)] = up_dag
        ann[(k, :dn)] = dn_ann
        dag[(k, :dn)] = dn_dag
    end
    return registry, (; ann, dag)
end

function build_h4_hamiltonian(h1e, eri; nk::Int, norb::Int)
    registry, vars = create_periodic_fermions(; nk, norb)
    seed = vars.dag[(0, :up)][1] * vars.ann[(0, :up)][1]
    ham = (0.0 + 0.0im) * seed

    for k in 0:(nk - 1)
        hk = h1e[k] / nk
        for p in 1:norb, q in 1:norb, spin in SPINS
            coeff = hk[p, q]
            iszero(coeff) && continue
            ham += coeff * vars.dag[(k, spin)][p] * vars.ann[(k, spin)][q]
        end
    end

    spin_channels = ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))
    for key in sort!(collect(keys(eri)))
        k1, k2, k3, k4 = key
        V = eri[key] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff = V[p, r, q, s]
            iszero(coeff) && continue
            for (σ, τ) in spin_channels
                ham += 0.5 * coeff *
                    vars.dag[(k1, σ)][p] * vars.dag[(k2, τ)][q] *
                    vars.ann[(k4, τ)][s] * vars.ann[(k3, σ)][r]
            end
        end
    end

    ham = 0.5 * (ham + adjoint(ham))
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
        particle_sign = raw < 0 ? 1 : -1
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

function grouped_blocks_with_identity(basis::Vector{NM}, identity::NM, blocking::Symbol; nk::Int, norb::Int) where {NM<:NormalMonomial}
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

function total_electron_constraint(vars, ham; nk::Int, norb::Int, total_electrons::Int)
    number = zero(ham)
    for k in 0:(nk - 1), i in 1:norb
        number += vars.dag[(k, :up)][i] * vars.ann[(k, :up)][i]
        number += vars.dag[(k, :dn)][i] * vars.ann[(k, :dn)][i]
    end
    return number - Float64(total_electrons) * one(ham)
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
        push!(constraints, trD - Float64(total_electrons * (total_electrons - 1) ÷ 2) * one(ham))
    end
    push!(constraints,
        trQ - Float64(n_holes * (n_holes - 1) ÷ 2) * one(ham),
        trG - Float64(total_electrons * (n_holes + 1)) * one(ham),
    )
    return constraints
end

function spin_resolved_d_trace_constraints(spin_orbitals, ham; total_electrons::Int)
    iseven(total_electrons) || throw(ArgumentError("spin-resolved traces assume M_s = 0"))
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
    nk = options.nk
    norb = options.norb
    nelec_per_cell = options.nelec_per_cell
    total_electrons = nk * nelec_per_cell

    h1e, eri = options.integrals_path === nothing ?
        mock_h4_integrals(; nk, norb) :
        load_integrals_txt(options.integrals_path; nk, norb)

    registry, vars, h4_ham = build_h4_hamiltonian(h1e, eri; nk, norb)
    objective = h4_ham
    spin_orbitals = build_spin_orbitals(vars; nk, norb)

    # The explicit N̂ - N localizing rows are redundant with TrD/TrG plus
    # the identity-augmented ²G block, and they create orphan moments in the
    # BPSDP lowering. Keep the physical electron count through TrD/TrQ/TrG,
    # spin-resolved D traces, and optional singlet S² instead.
    meq_constraints = typeof(h4_ham)[]
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
    zero_sizes = [size(mat) for (cone, mat) in mp.constraints if cone == :Zero]
    total_canonical_moments = length(NCTSSoS._sorted_symmetric_basis(mp.total_basis))
    return (; hpsd_sizes, zero_sizes,
             total_canonical_moments,
             direct_real_moment_variables = 2 * total_canonical_moments,
             real_lift_psd_rows = sum(real_lift_triangle_rows(n) for n in hpsd_sizes),
             complex_zero_entries = sum((prod(size) for size in zero_sizes); init = 0),
             direct_real_zero_rows = 2 * sum((prod(size) for size in zero_sizes); init = 0))
end

function print_summary(data, options::Options, build_seconds::Real)
    stats = constraint_stats(data.moment_problem)
    println("== H4/Nk=$(options.nk) periodic PQG V2RDM as an NCTSSoS MomentProblem ==")
    @printf("%-44s %s\n", "integrals", options.integrals_path === nothing ? "Julia deterministic mock" : options.integrals_path)
    @printf("%-44s %s\n", "blocking", string(options.blocking))
    @printf("%-44s %s\n", "extra ¹D PSD block", options.include_one_d ? "included" : "disabled")
    @printf("%-44s %s\n", "particle constraint", "dropped; implied by TrD/TrG + identity-augmented ²G")
    @printf("%-44s %s\n", "spin-resolved D traces", string(options.spin_resolved_trace))
    @printf("%-44s %s\n", "singlet S² constraint", string(options.singlet_s2))
    println()
    @printf("%-44s %d\n", "k-points", data.nk)
    @printf("%-44s %d\n", "active spatial orbitals per k", data.norb)
    @printf("%-44s %d\n", "spin-orbital modes", length(data.spin_orbitals))
    @printf("%-44s %d\n", "active electrons total", data.total_electrons)
    @printf("%-44s %.12f Ha\n", "mock/active HF energy", data.hf_active)
    @printf("%-44s %d\n", "Hamiltonian monomials", length(monomials(data.h4_ham)))
    println()
    println("-- Moment support basis --")
    @printf("%-44s %s\n", "¹D rows a_p", options.include_one_d ? string(length(data.row_bases.oneD_basis)) : "not included")
    @printf("%-44s %d\n", "²D rows a_p a_q", length(data.row_bases.twoD_basis))
    @printf("%-44s %d\n", "²Q rows a†_p a†_q", length(data.row_bases.twoQ_basis))
    @printf("%-44s %d\n", "²G rows a†_p a_q", length(data.row_bases.twoG_basis))
    @printf("%-44s %d\n", "total support rows", length(data.moment_support_basis))
    println()
    println("-- Schouten HPSD blocks --")
    @printf("%-44s %s\n", "¹D", options.include_one_d ? string(length.(data.oneD_blocks)) : "not included")
    @printf("%-44s %s\n", "²D", string(length.(data.twoD_blocks)))
    @printf("%-44s %s\n", "²Q", string(length.(data.twoQ_blocks)))
    @printf("%-44s %s\n", "²G", string(length.(data.twoG_blocks)))
    @printf("%-44s %d\n", "total HPSD blocks", length(stats.hpsd_sizes))
    @printf("%-44s %d\n", "largest HPSD block", maximum(stats.hpsd_sizes))
    println()
    println("-- Symbolic moment problem --")
    @printf("%-44s %d\n", "unique moment-matrix monomials", data.moment_problem.n_unique_moment_matrix_elements)
    @printf("%-44s %d\n", "total canonical monomials", stats.total_canonical_moments)
    @printf("%-44s %d\n", "direct real moment variables", stats.direct_real_moment_variables)
    @printf("%-44s %d\n", "real-lift PSD scalar rows", stats.real_lift_psd_rows)
    @printf("%-44s %d\n", "moment_eq polynomials", length(data.meq_constraints))
    @printf("%-44s %d\n", "Zero matrix constraints", length(stats.zero_sizes))
    @printf("%-44s %d\n", "Zero complex scalar entries", stats.complex_zero_entries)
    @printf("%-44s %.3f s\n", "build walltime", build_seconds)
    @printf("%-44s %s\n", "max RSS after build", rss_string())
    println()
    return nothing
end

function main(argv = ARGS)
    options = parse_options(argv)
    @printf("%-44s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h4_pqg_moment_problem(options)
    print_summary(data, options, build_seconds)

    @printf("%-44s %s\n", "writing libsdp export to", options.outdir)
    @printf("%-44s %s\n", "basename", options.basename)
    @printf("%-44s %s\n", "sense", string(options.sense))
    @printf("%-44s %d\n", "orphans per aux block", options.orphans_per_block)
    flush(stdout)

    export_seconds = @elapsed summary = export_libsdp(data.moment_problem;
        outdir = options.outdir,
        basename = options.basename,
        sense = options.sense,
        orphans_per_block = options.orphans_per_block,
    )

    dat_path = joinpath(options.outdir, options.basename * ".dat-c")
    meta_path = joinpath(options.outdir, options.basename * "_meta.txt")
    println()
    println("== libsdp export summary ==")
    @printf("%-44s %s\n", "SDPA file", dat_path)
    @printf("%-44s %s\n", "metadata", meta_path)
    @printf("%-44s %d (HPSD=%d, aux=%d)\n", "PSD blocks", summary.n_blocks, summary.n_hpsd_blocks, summary.n_aux_blocks)
    @printf("%-44s %d / %d\n", "min/max block dim", minimum(summary.block_dims), maximum(summary.block_dims))
    @printf("%-44s %d\n", "unique canonical moments", summary.n_unique_moments)
    @printf("%-44s %d\n", "equality constraints", summary.n_constraints)
    @printf("%-44s %d\n", "objective sparse entries", summary.n_objective_entries)
    @printf("%-44s %d\n", "constraint sparse entries", summary.n_constraint_entries)
    @printf("%-44s %.3f s\n", "export walltime", export_seconds)
    @printf("%-44s %s\n", "max RSS after export", rss_string())
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
