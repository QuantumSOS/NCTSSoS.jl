#!/usr/bin/env julia
raw"""
H2/Nk=2 periodic PQG V2RDM Phase-2 demo.

Builds the H2/Nk=2 active-space MomentProblem, dumps per-symmetry-block
constraint operators for A*A' diagnostics, solves the SDP with BPSDP.jl, and
writes solution-side moment-rank diagnostics.

Usage from the repository root:

    easy-ssh run 'cd /home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark
    tmpdir=$(mktemp -d)
    julia --startup-file=no --project="$tmpdir" -e "using Pkg; \
      Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\"); \
      Pkg.develop(path=pwd()); \
      Pkg.add([\"JuMP\",\"NPZ\",\"JSON3\",\"TOML\",\"Printf\",\"LinearAlgebra\",\"SparseArrays\"]); \
      Pkg.instantiate()"
    julia --startup-file=no --project="$tmpdir" demos/h2_periodic_nk2_moment_sos.jl \
        --output-dir=output/phase2/h2_nk2'

The integral text format is the PySCF dump used in the replication notes:

    h1e k p q Re Im
    eri k1 k2 k3 k4 p r q s Re Im

with ERIs in chemist order (p r | q s). This script normalizes one-electron
integrals by 1/Nk and ERIs by 1/Nk^2 before building the per-cell Hamiltonian.
"""

using NCTSSoS
using JuMP
const MOI = JuMP.MOI
using LinearAlgebra
using SparseArrays
using Printf
using JSON3
using NPZ
using TOML

try
    @eval import BPSDP
catch err
    error("BPSDP.jl is mandatory for Phase 2. Run this on HAI with Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\") as documented in PHASE2_PLAN.md §12.6. Original load error: $(err)")
end

const SPINS = (:up, :dn)
const DEFAULT_INTEGRALS = normpath(joinpath(
    @__DIR__, "..", "test", "data", "assets", "h2_chain_nk2_active_2e4o_integrals.txt"))
const DEFAULT_REFERENCE = normpath(joinpath(
    @__DIR__, "..", "test", "data", "assets", "h2_chain_nk2_reference.toml"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "output", "phase2", "h2_nk2"))
const DEFAULT_RANK_TOLS = [1e-10, 1e-12, 1e-14]

struct Options
    integrals_path::String
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
    output_dir::String
    rank_tols::Vector{Float64}
    diagnostic_block::String
    bpsdp_max_iter::Int
    bpsdp_cg_max_iter::Int
    bpsdp_mu_update_frequency::Int
    bpsdp_penalty::Float64
    bpsdp_cg_tol::Float64
    bpsdp_obj_tol::Float64
    bpsdp_err_tol::Float64
    bpsdp_print_level::Int
end

function parse_rank_tols(value::AbstractString)
    toks = split(value, ",")
    isempty(toks) && throw(ArgumentError("--rank-tol requires at least one tolerance"))
    return [parse(Float64, strip(tok)) for tok in toks]
end

function parse_options(argv)
    integrals_path = DEFAULT_INTEGRALS
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = false
    singlet_s2 = false
    output_dir = DEFAULT_OUTPUT_DIR
    rank_tols = copy(DEFAULT_RANK_TOLS)
    diagnostic_block = ""
    bpsdp_max_iter = 5_000
    bpsdp_cg_max_iter = 100
    bpsdp_mu_update_frequency = 25
    bpsdp_penalty = 0.1
    bpsdp_cg_tol = 1e-12
    bpsdp_obj_tol = 1e-8
    bpsdp_err_tol = 1e-8
    bpsdp_print_level = 1

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
        elseif startswith(arg, "--output-dir=")
            output_dir = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--rank-tol=")
            rank_tols = parse_rank_tols(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--diagnostic-block=")
            diagnostic_block = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--bpsdp-max-iter=")
            bpsdp_max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-cg-max-iter=")
            bpsdp_cg_max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-mu-update-frequency=")
            bpsdp_mu_update_frequency = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-penalty=")
            bpsdp_penalty = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-cg-tol=")
            bpsdp_cg_tol = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-obj-tol=")
            bpsdp_obj_tol = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-err-tol=")
            bpsdp_err_tol = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-print-level=")
            bpsdp_print_level = parse(Int, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=<tmp-with-BPSDP> demos/h2_periodic_nk2_moment_sos.jl [options]\n")
            println("Options:")
            println("  --integrals=PATH                  text integral dump; default test/data/assets/h2_chain_nk2_active_2e4o_integrals.txt")
            println("  --blocking=momentum|spin|none     1D/PQG block grouping; default momentum")
            println("  --include-1d / --no-1d            include/drop extra ¹D PSD block; default --no-1d")
            println("  --paper-spin, --spin-singlet      spin-resolved ²D traces plus singlet S²")
            println("  --spin-resolved-trace / --no-spin-resolved-trace")
            println("  --singlet-s2 / --no-singlet-s2")
            println("  --output-dir=PATH                 default output/phase2/h2_nk2")
            println("  --rank-tol=t1,t2,t3               default 1e-10,1e-12,1e-14")
            println("  --diagnostic-block=LABEL          default largest dumped block")
            println("  --bpsdp-max-iter=N                default 5000")
            println("  --bpsdp-cg-max-iter=N             default 100")
            println("  --bpsdp-mu-update-frequency=N     default 25")
            println("  --bpsdp-penalty=rho               default 0.1")
            println("  --bpsdp-cg-tol=eps                default 1e-12")
            println("  --bpsdp-obj-tol=eps               default 1e-8")
            println("  --bpsdp-err-tol=eps               default 1e-8")
            println("  --bpsdp-print-level=N             default 1")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    return Options(integrals_path, blocking, include_one_d,
        spin_resolved_trace, singlet_s2, output_dir, rank_tols, diagnostic_block,
        bpsdp_max_iter, bpsdp_cg_max_iter, bpsdp_mu_update_frequency,
        bpsdp_penalty, bpsdp_cg_tol, bpsdp_obj_tol, bpsdp_err_tol,
        bpsdp_print_level)
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
            block = get!(h1e, k) do
                zeros(ComplexF64, norb, norb)
            end
            block[p, q] = complex_entry(fields, 5, 6)
        elseif tag == "eri"
            length(fields) == 11 || error("malformed eri line $line_no in $(repr(path))")
            k1 = parse(Int, fields[2])
            k2 = parse(Int, fields[3])
            k3 = parse(Int, fields[4])
            k4 = parse(Int, fields[5])
            p = parse(Int, fields[6]) + 1
            r = parse(Int, fields[7]) + 1
            q = parse(Int, fields[8]) + 1
            s = parse(Int, fields[9]) + 1
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

function build_h2_nk2_hamiltonian(h1e, eri; nk::Int, norb::Int)
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

    for k in 0:(nk - 1)
        hk = h1e[k] / nk
        for p in 1:norb, q in 1:norb, spin in SPINS
            coeff = hk[p, q]
            iszero(coeff) && continue
            ham += coeff * dag[(k, spin)][p] * ann[(k, spin)][q]
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
                    dag[(k1, σ)][p] * dag[(k2, τ)][q] *
                    ann[(k4, τ)][s] * ann[(k3, σ)][r]
            end
        end
    end

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

sanitize_label(s::AbstractString) = replace(s, "-" => "m", "+" => "p", "/" => "_")

function key_suffix(key::Tuple, blocking::Symbol)
    blocking == :none && return "all"
    if blocking == :momentum
        return "K$(key[1])"
    end
    return sanitize_label("K$(key[1])_twoSz$(key[2])")
end

function grouped_named_blocks(basis::Vector{NM}, prefix::AbstractString, blocking::Symbol; nk::Int, norb::Int) where {NM<:NormalMonomial}
    if blocking == :none
        return [(label = string(prefix, "_all"), block = basis, key = ())]
    end

    buckets = Dict{Tuple,Vector{NM}}()
    for mono in basis
        key = block_key(mono, blocking; nk, norb)
        push!(get!(buckets, key, NM[]), mono)
    end

    out = NamedTuple[]
    for key in sort!(collect(keys(buckets)))
        push!(out, (label = string(prefix, "_", key_suffix(key, blocking)), block = buckets[key], key = key))
    end
    return out
end

function grouped_named_blocks_with_identity(
    basis::Vector{NM},
    identity::NM,
    prefix::AbstractString,
    blocking::Symbol;
    nk::Int,
    norb::Int,
) where {NM<:NormalMonomial}
    if blocking == :none
        return [(label = string(prefix, "_all"), block = NM[identity; basis], key = ())]
    end

    buckets = Dict{Tuple,Vector{NM}}()
    for mono in basis
        key = block_key(mono, blocking; nk, norb)
        push!(get!(buckets, key, NM[]), mono)
    end
    identity_key = blocking == :momentum ? (0,) : (0, 0)
    pushfirst!(get!(buckets, identity_key, NM[]), identity)

    out = NamedTuple[]
    for key in sort!(collect(keys(buckets)))
        push!(out, (label = string(prefix, "_", key_suffix(key, blocking)), block = buckets[key], key = key))
    end
    return out
end

function grouped_blocks(basis::Vector{NM}, blocking::Symbol; nk::Int, norb::Int) where {NM<:NormalMonomial}
    return [entry.block for entry in grouped_named_blocks(basis, "block", blocking; nk, norb)]
end

function grouped_blocks_with_identity(
    basis::Vector{NM}, identity::NM, blocking::Symbol; nk::Int, norb::Int,
) where {NM<:NormalMonomial}
    return [entry.block for entry in grouped_named_blocks_with_identity(basis, identity, "block", blocking; nk, norb)]
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

function build_h2_pqg_moment_problem(options::Options)
    nk = 2
    norb = 4
    nelec_per_cell = 2
    total_electrons = nk * nelec_per_cell

    h1e, eri = load_integrals_txt(options.integrals_path; norb)
    registry, vars, h2_ham = build_h2_nk2_hamiltonian(h1e, eri; nk, norb)
    objective = h2_ham
    spin_orbitals = build_spin_orbitals(vars; nk, norb)

    meq_constraints = typeof(h2_ham)[
        total_electron_constraint(vars, h2_ham; norb, total_electrons),
    ]
    append!(meq_constraints, trace_constraints(spin_orbitals, h2_ham;
        total_electrons,
        include_total_d = !options.spin_resolved_trace,
    ))
    if options.spin_resolved_trace
        append!(meq_constraints,
            spin_resolved_d_trace_constraints(spin_orbitals, h2_ham; total_electrons))
    end
    if options.singlet_s2
        push!(meq_constraints,
            singlet_s2_constraint(vars, h2_ham; nk, norb, total_electrons))
    end

    pop = polyopt(objective, registry; moment_eq_constraints = meq_constraints)
    row_bases = build_rdm_bases(spin_orbitals)
    NM = eltype(row_bases.oneD_basis)
    identity = one(NM)

    oneD_named = options.include_one_d ?
        grouped_named_blocks(row_bases.oneD_basis, "oneD", options.blocking; nk, norb) :
        NamedTuple[]
    twoD_named = grouped_named_blocks(row_bases.twoD_basis, "twoD", options.blocking; nk, norb)
    twoQ_named = grouped_named_blocks(row_bases.twoQ_basis, "twoQ", options.blocking; nk, norb)
    twoG_named = grouped_named_blocks_with_identity(row_bases.twoG_basis, identity, "twoG", options.blocking; nk, norb)

    oneD_blocks = [entry.block for entry in oneD_named]
    twoD_blocks = [entry.block for entry in twoD_named]
    twoQ_blocks = [entry.block for entry in twoQ_named]
    twoG_blocks = [entry.block for entry in twoG_named]
    block_labels = String[entry.label for entry in vcat(oneD_named, twoD_named, twoQ_named, twoG_named)]

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
             registry, vars, h2_ham, objective, pop, meq_constraints,
             spin_orbitals, row_bases, oneD_blocks, twoD_blocks, twoQ_blocks,
             twoG_blocks, block_labels, moment_support_basis, schouten_blocks,
             term_sparsity, moment_problem, hf_active)
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

# -----------------------------------------------------------------------------
# Constraint-operator materialization.
# -----------------------------------------------------------------------------

function moment_basis(mp)
    basis = [symmetric_canon(NCTSSoS.expval(m)) for m in mp.total_basis]
    sort!(unique!(basis))
    return basis, Dict(m => i for (i, m) in enumerate(basis))
end

@inline function add_to_row!(row::Dict{Int,Float64}, col::Int, value::Real)
    v = Float64(value)
    iszero(v) && return row
    row[col] = get(row, col, 0.0) + v
    return row
end

function polynomial_real_imag_rows(poly, basis_to_idx::Dict, n_basis::Int)
    row_re = Dict{Int,Float64}()
    row_im = Dict{Int,Float64}()

    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        canon_mono = symmetric_canon(NCTSSoS.expval(mono))
        idx = get(basis_to_idx, canon_mono, 0)
        iszero(idx) && continue
        c = ComplexF64(coef)
        add_to_row!(row_re, idx, real(c))
        add_to_row!(row_re, n_basis + idx, -imag(c))
        add_to_row!(row_im, idx, imag(c))
        add_to_row!(row_im, n_basis + idx, real(c))
    end

    return row_re, row_im
end

negate_row(row::Dict{Int,Float64}) = Dict(k => -v for (k, v) in row)

function rows_to_sparse(rows::Vector{Dict{Int,Float64}}, ncols::Int)
    I = Int[]
    J = Int[]
    V = Float64[]
    for (r, row) in enumerate(rows)
        for (c, v) in row
            abs(v) <= 0.0 && continue
            push!(I, r)
            push!(J, c)
            push!(V, v)
        end
    end
    return sparse(I, J, V, length(rows), ncols)
end

function build_zero_operator(mp, basis_to_idx::Dict, n_basis::Int)
    rows = Dict{Int,Float64}[]
    rhs = Float64[]
    labels = String[]
    ncols = 2 * n_basis

    one_sym = symmetric_canon(NCTSSoS.expval(one(eltype(mp.total_basis))))
    idx_one = get(basis_to_idx, one_sym, 0)
    iszero(idx_one) && error("identity moment missing from basis")

    row = Dict{Int,Float64}(); add_to_row!(row, idx_one, 1.0)
    push!(rows, row); push!(rhs, 1.0); push!(labels, "normalization:Re(<1>)=1")
    row = Dict{Int,Float64}(); add_to_row!(row, n_basis + idx_one, 1.0)
    push!(rows, row); push!(rhs, 0.0); push!(labels, "normalization:Im(<1>)=0")

    zero_counter = 0
    for (cone, mat) in mp.constraints
        cone == :Zero || continue
        zero_counter += 1
        for j in axes(mat, 2), i in axes(mat, 1)
            re_row, im_row = polynomial_real_imag_rows(mat[i, j], basis_to_idx, n_basis)
            push!(rows, re_row); push!(rhs, 0.0); push!(labels, "Zero[$zero_counter]($i,$j):Re")
            push!(rows, im_row); push!(rhs, 0.0); push!(labels, "Zero[$zero_counter]($i,$j):Im")
        end
    end

    return rows_to_sparse(rows, ncols), rhs, labels
end

function build_hpsd_operator(mat, basis_to_idx::Dict, n_basis::Int)
    d = size(mat, 1)
    rows = Dict{Int,Float64}[]
    labels = String[]
    for col in 1:(2d), row in 1:col
        if row <= d && col <= d
            re_row, _ = polynomial_real_imag_rows(mat[row, col], basis_to_idx, n_basis)
            push!(rows, re_row)
            push!(labels, "HPSD($row,$col):Re")
        elseif row <= d && col > d
            _, im_row = polynomial_real_imag_rows(mat[row, col - d], basis_to_idx, n_basis)
            push!(rows, negate_row(im_row))
            push!(labels, "HPSD($row,$(col - d)):minus_Im")
        elseif row > d && col > d
            re_row, _ = polynomial_real_imag_rows(mat[row - d, col - d], basis_to_idx, n_basis)
            push!(rows, re_row)
            push!(labels, "HPSD($(row - d),$(col - d)):Re_copy")
        else
            # Not reached for upper-triangle vectorization, but keep the branch
            # explicit so the real lift is readable.
            _, im_row = polynomial_real_imag_rows(mat[row - d, col], basis_to_idx, n_basis)
            push!(rows, im_row)
            push!(labels, "HPSD($(row - d),$col):Im")
        end
    end
    return rows_to_sparse(rows, 2 * n_basis), labels
end

function row_nonzero_mask(A::SparseMatrixCSC)
    mask = falses(size(A, 1))
    I, _, _ = findnz(A)
    mask[I] .= true
    return mask
end

function active_columns(A::SparseMatrixCSC)
    _, J, _ = findnz(A)
    cols = sort!(unique(J))
    return isempty(cols) ? Int[] : cols
end

function structural_redundancy(A::SparseMatrixCSC)
    M = Matrix(A)
    signatures = Dict{String,Int}()
    zero_rows = 0
    for r in axes(M, 1)
        parts = String[]
        for c in axes(M, 2)
            v = M[r, c]
            abs(v) <= 1e-13 && continue
            push!(parts, string(c, ":", @sprintf("%.12e", v)))
        end
        if isempty(parts)
            zero_rows += 1
        else
            sig = join(parts, ",")
            signatures[sig] = get(signatures, sig, 0) + 1
        end
    end
    nonzero_rows = size(A, 1) - zero_rows
    duplicate_rows = nonzero_rows - length(signatures)
    return (; zero_rows, duplicate_rows, predicted_kernel_dim = zero_rows + duplicate_rows)
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, obj)
        write(io, '\n')
    end
    return path
end

function write_sparse_npz(path::AbstractString, A::SparseMatrixCSC)
    mkpath(dirname(path))
    I, J, V = findnz(A)
    NPZ.npzwrite(path, Dict(
        "row" => Int64.(I .- 1),
        "col" => Int64.(J .- 1),
        "data" => Float64.(V),
        "shape" => Int64[size(A, 1), size(A, 2)],
        "index_base" => Int64[0],
    ))
    return path
end

function finite_or_string(x::Real)
    xf = Float64(x)
    return isfinite(xf) ? xf : string(xf)
end

function tol_key(tol::Real)
    return @sprintf("%.0e", tol)
end

function spectrum_diagnostic(A::SparseMatrixCSC, rank_tols::Vector{Float64})
    m = size(A, 1)
    n = size(A, 2)
    if m == 0
        return Dict(
            "rows" => 0, "cols" => n, "sigma1" => 0.0,
            "rank_by_tol" => Dict{String,Any}(),
            "singular_values" => Float64[],
        ), zeros(Float64, 0, 0)
    end

    G = Matrix(A * transpose(A))
    G = 0.5 .* (G .+ transpose(G))
    F = eigen(Symmetric(G))
    vals_asc = max.(F.values, 0.0)
    order = sortperm(vals_asc; rev = true)
    sigma = vals_asc[order]
    vecs_desc = F.vectors[:, order]
    sigma1 = isempty(sigma) ? 0.0 : sigma[1]

    ranks = Dict{String,Any}()
    for tol in rank_tols
        threshold = tol * sigma1
        r = count(>(threshold), sigma)
        sigma_r = r > 0 ? sigma[r] : 0.0
        sigma_next = r < length(sigma) ? sigma[r + 1] : 0.0
        gap = sigma_next > 0 ? sigma_r / sigma_next : Inf
        kappa = sigma_r > 0 ? sigma1 / sigma_r : Inf
        ranks[tol_key(tol)] = Dict(
            "threshold" => threshold,
            "rank" => r,
            "sigma_rank" => finite_or_string(sigma_r),
            "sigma_next" => finite_or_string(sigma_next),
            "kernel_dim" => m - r,
            "kappa_range" => finite_or_string(kappa),
            "spectral_gap" => finite_or_string(gap),
            "estimated_cg_iterations_1e-8" => finite_or_string(isfinite(kappa) ? 9.0 * sqrt(kappa) : Inf),
        )
    end

    default_tol = 1e-12 in rank_tols ? 1e-12 : rank_tols[min(2, length(rank_tols))]
    kernel_cols = findall(<=(default_tol * sigma1), sigma)
    kernel_basis = isempty(kernel_cols) ? zeros(Float64, m, 0) : vecs_desc[:, kernel_cols]

    return Dict(
        "rows" => m,
        "cols" => n,
        "sigma1" => finite_or_string(sigma1),
        "rank_by_tol" => ranks,
        "singular_values" => sigma,
        "threshold_mode" => "relative_to_sigma1",
    ), kernel_basis
end

function sparsity_diagnostic(A::SparseMatrixCSC)
    G = sparse(A * transpose(A))
    m = size(G, 1)
    nnz_G = nnz(G)
    density = m == 0 ? 0.0 : nnz_G / (m * m)
    if m == 0
        return Dict(
            "aat_rows" => 0, "aat_nnz" => 0, "aat_density" => 0.0,
            "amd_nnz_L" => 0, "amd_memory_bytes" => 0,
            "metis_available" => false,
        )
    end

    delta = max(1.0, maximum(abs, G.nzval; init = 0.0)) * 1e-12
    shifted = G + delta * sparse(I, m, m)
    amd_nnz_L = try
        F = cholesky(Symmetric(shifted); check = false)
        nnz(sparse(F.L))
    catch err
        @warn "AMD Cholesky fill measurement failed" err
        -1
    end

    return Dict(
        "aat_rows" => m,
        "aat_nnz" => nnz_G,
        "aat_density" => density,
        "regularization_delta" => delta,
        "amd_nnz_L" => amd_nnz_L,
        "amd_memory_bytes" => amd_nnz_L < 0 ? -1 : 8 * amd_nnz_L,
        "metis_available" => false,
        "metis_note" => "METIS not required/available in the Phase-2 demo environment; AMD measured with CHOLMOD.",
    )
end

function write_diagnostics!(block_dir::AbstractString, A_full::SparseMatrixCSC, A_eq::SparseMatrixCSC, rank_tols::Vector{Float64})
    spec_full, ker_full = spectrum_diagnostic(A_full, rank_tols)
    spec_eq, ker_eq = spectrum_diagnostic(A_eq, rank_tols)
    write_json(joinpath(block_dir, "spectrum.json"), Dict("A_full" => spec_full, "A_eq" => spec_eq))
    write_json(joinpath(block_dir, "sparsity.json"), Dict(
        "A_full" => sparsity_diagnostic(A_full),
        "A_eq" => sparsity_diagnostic(A_eq),
    ))
    NPZ.npzwrite(joinpath(block_dir, "kernel_basis.npz"), Dict("A_full" => ker_full, "A_eq" => ker_eq))
    return nothing
end

function dump_block_operators(data, options::Options)
    mp = data.moment_problem
    basis, basis_to_idx = moment_basis(mp)
    n_basis = length(basis)
    ncols = 2 * n_basis
    zero_A, zero_rhs, zero_labels = build_zero_operator(mp, basis_to_idx, n_basis)

    block_root = joinpath(options.output_dir, "blocks")
    mkpath(block_root)
    hpsd_mats = [mat for (cone, mat) in mp.constraints if cone == :HPSD]
    length(hpsd_mats) == length(data.block_labels) ||
        error("HPSD block count $(length(hpsd_mats)) does not match labels $(length(data.block_labels))")

    artifacts = NamedTuple[]
    for (block_index, mat) in enumerate(hpsd_mats)
        label = data.block_labels[block_index]
        block_dir = joinpath(block_root, label)
        mkpath(block_dir)

        hpsd_A_global, hpsd_labels = build_hpsd_operator(mat, basis_to_idx, n_basis)
        cols = active_columns(hpsd_A_global)
        isempty(cols) && (cols = collect(1:ncols))
        hpsd_A = hpsd_A_global[:, cols]
        eq_restricted = zero_A[:, cols]
        eq_keep = row_nonzero_mask(eq_restricted)
        A_eq = eq_restricted[eq_keep, :]
        eq_labels = zero_labels[eq_keep]
        A_full = vcat(A_eq, hpsd_A)

        write_sparse_npz(joinpath(block_dir, "A_full.npz"), A_full)
        write_sparse_npz(joinpath(block_dir, "A_eq.npz"), A_eq)

        redundancy_full = structural_redundancy(A_full)
        redundancy_eq = structural_redundancy(A_eq)
        labels_json = Dict(
            "block_label" => label,
            "block_index" => block_index,
            "block_dim_complex" => size(mat, 1),
            "blocking" => string(options.blocking),
            "columns_global_1based" => cols,
            "columns" => [sprint(show, basis[((c - 1) % n_basis) + 1]) * (c <= n_basis ? ":Re" : ":Im") for c in cols],
            "rows" => Dict(
                "A_eq" => eq_labels,
                "A_hpsd" => hpsd_labels,
                "A_full" => vcat(eq_labels, hpsd_labels),
            ),
            "predicted_redundancy" => Dict(
                "A_full" => Dict(
                    "zero_rows" => redundancy_full.zero_rows,
                    "duplicate_rows" => redundancy_full.duplicate_rows,
                    "predicted_kernel_dim" => redundancy_full.predicted_kernel_dim,
                    "note" => "Construction-level lower bound: zero rows plus exact duplicate rows after local column restriction. Excess kernel is empirical minus this count.",
                ),
                "A_eq" => Dict(
                    "zero_rows" => redundancy_eq.zero_rows,
                    "duplicate_rows" => redundancy_eq.duplicate_rows,
                    "predicted_kernel_dim" => redundancy_eq.predicted_kernel_dim,
                    "note" => "Construction-level lower bound: zero rows plus exact duplicate rows after local column restriction.",
                ),
            ),
            "matrix_shapes" => Dict(
                "A_full" => [size(A_full, 1), size(A_full, 2)],
                "A_eq" => [size(A_eq, 1), size(A_eq, 2)],
            ),
        )
        write_json(joinpath(block_dir, "labels.json"), labels_json)

        push!(artifacts, (label = label, dir = block_dir, rows_full = size(A_full, 1), rows_eq = size(A_eq, 1), cols = size(A_full, 2)))
    end

    return (; basis, basis_to_idx, zero_A, zero_rhs, zero_labels, artifacts)
end

# -----------------------------------------------------------------------------
# JuMP/BPSDP solve and solution diagnostics.
# -----------------------------------------------------------------------------

function build_direct_jump_model(mp, block_labels::Vector{String})
    C = real(eltype(coefficients(mp.objective)))
    model = JuMP.GenericModel{C}()

    basis, basis_to_idx = moment_basis(mp)
    n_basis = length(basis)
    @variable(model, y_re[1:n_basis], set_string_name = false)
    @variable(model, y_im[1:n_basis], set_string_name = false)

    one_sym = symmetric_canon(NCTSSoS.expval(one(eltype(mp.total_basis))))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("identity moment missing from basis")

    eq_refs = NamedTuple[]
    c = @constraint(model, y_re[idx_one] == 1)
    push!(eq_refs, (label = "normalization:Re(<1>)=1", ref = c))
    c = @constraint(model, y_im[idx_one] == 0)
    push!(eq_refs, (label = "normalization:Im(<1>)=0", ref = c))

    hpsd_refs = NamedTuple[]
    zero_counter = 0
    hpsd_counter = 0

    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        mat_re = Matrix{Any}(undef, dim, dim)
        mat_im = Matrix{Any}(undef, dim, dim)
        for i in 1:dim, j in 1:dim
            re_expr, im_expr = NCTSSoS._substitute_complex_poly(mat[i, j], basis_to_idx, y_re, y_im)
            mat_re[i, j] = re_expr
            mat_im[i, j] = im_expr
        end

        if cone == :Zero
            zero_counter += 1
            for j in 1:dim, i in 1:dim
                cref = @constraint(model, mat_re[i, j] == 0)
                push!(eq_refs, (label = "Zero[$zero_counter]($i,$j):Re", ref = cref))
                cref = @constraint(model, mat_im[i, j] == 0)
                push!(eq_refs, (label = "Zero[$zero_counter]($i,$j):Im", ref = cref))
            end
        elseif cone == :HPSD
            hpsd_counter += 1
            embedded = [
                [mat_re[i, j] for i in 1:dim, j in 1:dim] [-mat_im[i, j] for i in 1:dim, j in 1:dim]
                [mat_im[i, j] for i in 1:dim, j in 1:dim] [mat_re[i, j] for i in 1:dim, j in 1:dim]
            ]
            cref = @constraint(model, embedded in PSDCone())
            label = hpsd_counter <= length(block_labels) ? block_labels[hpsd_counter] : "HPSD_$hpsd_counter"
            push!(hpsd_refs, (label = label, ref = cref, dim = dim))
        else
            error("unexpected cone $cone")
        end
    end

    obj_re, _ = NCTSSoS._substitute_complex_poly(mp.objective, basis_to_idx, y_re, y_im)
    @objective(model, Min, obj_re)

    return (; model, y_re, y_im, basis, basis_to_idx, eq_refs, hpsd_refs)
end

function progress_monitor_factory(cg_iterations::Vector{Int})
    return function (print_level, outer, inner, primal, dual, mu, perr, derr)
        push!(cg_iterations, Int(inner))
        if print_level > 0 && outer % print_level == 0
            @printf("BPSDP %6d %6d primal=% .10e dual=% .10e mu=% .3e perr=% .3e derr=% .3e\n",
                outer, inner, primal, dual, mu, perr, derr)
            flush(stdout)
        end
        return nothing
    end
end

function bpsdp_optimizer_factory(options::Options, cg_iterations::Vector{Int})
    monitor = progress_monitor_factory(cg_iterations)
    return () -> BPSDP.Optimizer(
        max_iter = options.bpsdp_max_iter,
        cg_max_iter = options.bpsdp_cg_max_iter,
        mu_update_frequency = options.bpsdp_mu_update_frequency,
        penalty_parameter = options.bpsdp_penalty,
        cg_convergence = options.bpsdp_cg_tol,
        dynamic_cg_convergence = false,
        sdp_objective_convergence = options.bpsdp_obj_tol,
        sdp_error_convergence = options.bpsdp_err_tol,
        guess_type = :zero,
        print_level = options.bpsdp_print_level,
        progress_monitor = monitor,
        dependent_rows = :drop,
    )
end

function raw_optimizer(model::JuMP.Model)
    try
        return JuMP.unsafe_backend(model)
    catch
        return nothing
    end
end

has_field(obj, name::Symbol) = name in fieldnames(typeof(obj))
field_or_nothing(obj, name::Symbol) = has_field(obj, name) ? getfield(obj, name) : nothing
function first_field_or_nothing(obj, names::Symbol...)
    for name in names
        has_field(obj, name) && return getfield(obj, name)
    end
    return nothing
end

finite_or_nothing(x) = x === nothing ? nothing : finite_or_string(x)

function bpsdp_state_summary(raw)
    state = raw === nothing ? nothing : field_or_nothing(raw, :state)
    state === nothing && return nothing
    summary = Dict{String,Any}(
        "outer_iterations" => first_field_or_nothing(state, :outer_iterations, :iterations),
        "inner_iterations" => field_or_nothing(state, :inner_iterations),
        "primal_error" => finite_or_nothing(field_or_nothing(state, :primal_error)),
        "dual_error" => finite_or_nothing(field_or_nothing(state, :dual_error)),
        "objective_primal" => finite_or_nothing(first_field_or_nothing(state, :objective_primal, :primal_objective)),
        "objective_dual" => finite_or_nothing(first_field_or_nothing(state, :objective_dual, :dual_objective)),
        "objective_gap" => finite_or_nothing(field_or_nothing(state, :objective_gap)),
        "termination_reason" => string(first_field_or_nothing(state, :termination_reason, :status)),
    )
    return summary
end

function poly_value(poly, monomap)
    acc = 0.0 + 0.0im
    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        key = symmetric_canon(NCTSSoS.expval(mono))
        acc += ComplexF64(coef) * get(monomap, key, 0.0 + 0.0im)
    end
    return acc
end

function evaluate_hpsd_matrix(mat, monomap)
    Y = Matrix{ComplexF64}(undef, size(mat, 1), size(mat, 2))
    for j in axes(mat, 2), i in axes(mat, 1)
        Y[i, j] = poly_value(mat[i, j], monomap)
    end
    return 0.5 .* (Y .+ adjoint(Y))
end

function eig_summary(vals::Vector{Float64}, rank_tols::Vector{Float64})
    vals_desc = sort(vals; rev = true)
    scale = isempty(vals_desc) ? 0.0 : max(abs(vals_desc[1]), 0.0)
    ranks = Dict{String,Any}()
    for tol in rank_tols
        threshold = tol * scale
        r = count(>(threshold), vals_desc)
        smallest = r > 0 ? vals_desc[r] : 0.0
        ranks[tol_key(tol)] = Dict(
            "threshold" => threshold,
            "rank" => r,
            "smallest_nonzero_eigenvalue" => finite_or_string(smallest),
        )
    end
    return Dict("eigenvalues" => vals_desc, "rank_by_tol" => ranks)
end

function write_solution_side_diagnostics(data, options::Options, jump_data, materialized, monomap, cg_iterations::Vector{Int}, solve_seconds::Real)
    model = jump_data.model
    basis = materialized.basis
    y_vec = vcat(
        [real(get(monomap, b, 0.0 + 0.0im)) for b in basis],
        [imag(get(monomap, b, 0.0 + 0.0im)) for b in basis],
    )

    residual = materialized.zero_A * y_vec - materialized.zero_rhs
    residual_linf = isempty(residual) ? 0.0 : norm(residual, Inf)
    residual_l2 = norm(residual)

    refs = isfile(DEFAULT_REFERENCE) ? TOML.parsefile(DEFAULT_REFERENCE) : Dict{String,Any}()
    constant_shift = get(get(refs, "active_hamiltonian", Dict{String,Any}()), "constant_shift_per_cell_Ha", NaN)
    energies = get(refs, "energies", Dict{String,Any}())
    objective = objective_value(model)

    raw = raw_optimizer(model)
    status = termination_status(model)

    solve_dir = joinpath(options.output_dir, "solve")
    mkpath(solve_dir)

    objective_json = Dict{String,Any}(
        "termination_status" => string(status),
        "objective_active_Ha" => objective,
        "objective_total_shifted_Ha" => isfinite(constant_shift) ? objective + constant_shift : NaN,
        "hf_active_reconstructed_Ha" => data.hf_active,
        "reference_total_Ha" => energies,
        "constant_shift_per_cell_Ha" => constant_shift,
        "primal_feasibility_linf" => residual_linf,
        "primal_feasibility_l2" => residual_l2,
        "solve_wall_seconds_measured" => solve_seconds,
        "cg_iterations_per_outer" => cg_iterations,
        "cg_iterations_total_from_monitor" => sum(cg_iterations),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "jump_lowering" => Dict(
            "formulation" => "psd_blocks",
            "representation" => "complex",
            "orphan_policy" => "free_variables",
        ),
        "bpsdp_options" => Dict(
            "max_iter" => options.bpsdp_max_iter,
            "cg_max_iter" => options.bpsdp_cg_max_iter,
            "mu_update_frequency" => options.bpsdp_mu_update_frequency,
            "penalty_parameter" => options.bpsdp_penalty,
            "cg_convergence" => options.bpsdp_cg_tol,
            "sdp_objective_convergence" => options.bpsdp_obj_tol,
            "sdp_error_convergence" => options.bpsdp_err_tol,
            "guess_type" => "zero",
            "dependent_rows" => "drop",
        ),
    )
    state_summary = bpsdp_state_summary(raw)
    state_summary !== nothing && (objective_json["bpsdp_state"] = state_summary)
    write_json(joinpath(solve_dir, "objective.json"), objective_json)

    write_json(joinpath(solve_dir, "multipliers.json"), Dict(
        "equality_multipliers" => Any[],
        "note" => "Unavailable from the PSD-block lowering path; equality labels remain in blocks/*/labels.json.",
    ))

    hpsd_mats = [mat for (cone, mat) in data.moment_problem.constraints if cone == :HPSD]
    for (idx, mat) in enumerate(hpsd_mats)
        label = data.block_labels[idx]
        block_dir = joinpath(options.output_dir, "blocks", label)
        Y = evaluate_hpsd_matrix(mat, monomap)
        vals = eigvals(Hermitian(Y))
        write_json(joinpath(block_dir, "Y_eigvals.json"), eig_summary(real.(vals), options.rank_tols))
    end

    return (; monomap, objective_json)
end

function write_bpsdp_failure!(options::Options, model, cg_iterations::Vector{Int}, solve_seconds::Real; optimize_error = nothing)
    solve_dir = joinpath(options.output_dir, "solve")
    mkpath(solve_dir)
    raw = raw_optimizer(model)
    objective_json = Dict{String,Any}(
        "termination_status" => try string(termination_status(model)) catch err string("unavailable: ", err) end,
        "raw_status" => try string(MOI.get(model, MOI.RawStatusString())) catch err string("unavailable: ", err) end,
        "optimize_error" => optimize_error,
        "solve_wall_seconds_measured" => solve_seconds,
        "cg_iterations_per_outer" => cg_iterations,
        "cg_iterations_total_from_monitor" => sum(cg_iterations),
        "jump_lowering" => Dict(
            "formulation" => "psd_blocks",
            "representation" => "complex",
            "orphan_policy" => "free_variables",
        ),
        "bpsdp_options" => Dict(
            "max_iter" => options.bpsdp_max_iter,
            "cg_max_iter" => options.bpsdp_cg_max_iter,
            "mu_update_frequency" => options.bpsdp_mu_update_frequency,
            "penalty_parameter" => options.bpsdp_penalty,
            "cg_convergence" => options.bpsdp_cg_tol,
            "sdp_objective_convergence" => options.bpsdp_obj_tol,
            "sdp_error_convergence" => options.bpsdp_err_tol,
            "guess_type" => "zero",
            "dependent_rows" => "drop",
        ),
    )
    state_summary = bpsdp_state_summary(raw)
    state_summary !== nothing && (objective_json["bpsdp_state"] = state_summary)
    if raw !== nothing
        objective_json["bpsdp_problem"] = Dict(
            "n_blocks" => has_field(raw, :block_dims) ? length(getfield(raw, :block_dims)) : nothing,
            "largest_block_dim" => has_field(raw, :block_dims) && !isempty(getfield(raw, :block_dims)) ? maximum(getfield(raw, :block_dims)) : nothing,
            "dual_rows" => has_field(raw, :b) ? length(getfield(raw, :b)) : nothing,
            "primal_dim" => has_field(raw, :c) ? length(getfield(raw, :c)) : nothing,
            "A_nnz" => has_field(raw, :A) ? nnz(getfield(raw, :A)) : nothing,
        )
    end
    write_json(joinpath(solve_dir, "objective.json"), objective_json)
    write_json(joinpath(solve_dir, "multipliers.json"), Dict(
        "equality_multipliers" => Any[],
        "note" => "BPSDP did not return OPTIMAL; no reliable dual multipliers available.",
    ))
    return objective_json
end

function solve_with_bpsdp!(data, options::Options, materialized)
    model, extract_monomap = build_jump_model(
        data.moment_problem;
        formulation = :psd_blocks,
        representation = :complex,
        orphan_policy = :free_variables,
    )
    jump_data = (; model, extract_monomap)
    cg_iterations = Int[]
    set_optimizer(model, bpsdp_optimizer_factory(options, cg_iterations))

    optimize_error = nothing
    solve_seconds = @elapsed begin
        try
            optimize!(model)
        catch err
            optimize_error = sprint(showerror, err)
        end
    end
    status = try
        termination_status(model)
    catch
        MOI.OTHER_ERROR
    end
    if optimize_error !== nothing || status != MOI.OPTIMAL
        write_bpsdp_failure!(options, model, cg_iterations, solve_seconds; optimize_error)
        optimize_error === nothing || error("BPSDP optimize! threw: $optimize_error")
        error("BPSDP did not return OPTIMAL; got $status")
    end

    monomap = extract_monomap()
    solution = write_solution_side_diagnostics(data, options, jump_data, materialized, monomap, cg_iterations, solve_seconds)
    return (; jump_data, cg_iterations, solve_seconds, solution)
end

# -----------------------------------------------------------------------------
# Printing.
# -----------------------------------------------------------------------------

function print_summary(data, options::Options, build_seconds::Real, materialized, diagnostic_label::String)
    stats = constraint_stats(data.moment_problem)
    n_modes = length(data.spin_orbitals)
    n_holes = n_modes - data.total_electrons

    formulation_label = options.include_one_d ?
        "H2/Nk=2 periodic 1D+PQG V2RDM Phase-2 build+solve" :
        "H2/Nk=2 periodic PQG V2RDM Phase-2 build+solve (no ¹D)"
    println("== ", formulation_label, " ==")
    @printf("%-48s %s\n", "integrals", options.integrals_path)
    @printf("%-48s %s\n", "blocking", string(options.blocking))
    @printf("%-48s %s\n", "extra ¹D PSD block", options.include_one_d ? "included" : "disabled")
    @printf("%-48s %s\n", "objective", "H2/Nk=2 active-space Hamiltonian")
    @printf("%-48s N=%d via moment_eq_constraints\n", "particle constraint", data.total_electrons)
    if options.spin_resolved_trace
        @printf("%-48s TrDαα=%d, TrDαβ=%d, TrDββ=%d, TrQ=%d, TrG=%d\n", "trace constraints",
            (data.total_electrons ÷ 2) * (data.total_electrons ÷ 2 - 1) ÷ 2,
            (data.total_electrons ÷ 2)^2,
            (data.total_electrons ÷ 2) * (data.total_electrons ÷ 2 - 1) ÷ 2,
            n_holes * (n_holes - 1) ÷ 2,
            data.total_electrons * (n_holes + 1))
    else
        @printf("%-48s TrD=%d, TrQ=%d, TrG=%d\n", "trace constraints",
            data.total_electrons * (data.total_electrons - 1) ÷ 2,
            n_holes * (n_holes - 1) ÷ 2,
            data.total_electrons * (n_holes + 1))
    end
    @printf("%-48s %s\n", "singlet S² constraint", options.singlet_s2 ? "included" : "disabled")
    @printf("%-48s %s\n", "correlative sparsity", "none; single full-system clique")
    println()

    @printf("%-48s %d\n", "k-points", data.nk)
    @printf("%-48s %d\n", "active spatial orbitals per k", data.norb)
    @printf("%-48s %d\n", "spin-orbital modes", n_modes)
    @printf("%-48s %d\n", "active electrons total", data.total_electrons)
    @printf("%-48s %.12f Ha\n", "HF active energy reconstructed", data.hf_active)
    @printf("%-48s %d\n", "Hamiltonian monomials", length(monomials(data.h2_ham)))
    println()

    println("-- Moment support basis --")
    @printf("%-48s %d\n", "identity", 1)
    @printf("%-48s %s\n", "¹D rows a_p", options.include_one_d ? string(length(data.row_bases.oneD_basis)) : "not included")
    @printf("%-48s %d\n", "²D rows a_p a_q, p<q", length(data.row_bases.twoD_basis))
    @printf("%-48s %d\n", "²Q rows a†_p a†_q, p<q", length(data.row_bases.twoQ_basis))
    @printf("%-48s %d\n", "²G rows a†_p a_q", length(data.row_bases.twoG_basis))
    @printf("%-48s %d\n", "total support rows", length(data.moment_support_basis))
    println()

    println("-- Schouten HPSD blocks --")
    @printf("%-48s %s\n", "¹D", options.include_one_d ? string(length.(data.oneD_blocks)) : "not included")
    @printf("%-48s %s\n", "²D", string(length.(data.twoD_blocks)))
    @printf("%-48s %s\n", "²Q", string(length.(data.twoQ_blocks)))
    @printf("%-48s %s\n", "²G", string(length.(data.twoG_blocks)))
    @printf("%-48s %d\n", "total HPSD blocks", length(stats.hpsd_sizes))
    @printf("%-48s %d\n", "largest HPSD block", maximum(stats.hpsd_sizes))
    @printf("%-48s %s\n", "block labels", join(data.block_labels, ", "))
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
    println()

    println("-- Phase-2 block artifacts --")
    @printf("%-48s %s\n", "output dir", options.output_dir)
    @printf("%-48s %d\n", "dumped blocks", length(materialized.artifacts))
    @printf("%-48s %s\n", "inline diagnostic block", diagnostic_label)
    println()

    @printf("%-48s %.3f s\n", "build walltime", build_seconds)
    @printf("%-48s %s\n", "max RSS after build/dump", rss_string())
    println()
    println("BPSDP solve is mandatory; optimize! will be called after this summary.")
    println("Artifact: output/phase2/h2_nk2 blocks + solve diagnostics.")
    return nothing
end

function choose_diagnostic_block(materialized, requested::AbstractString)
    if !isempty(requested)
        any(a.label == requested for a in materialized.artifacts) ||
            error("--diagnostic-block=$requested was not among dumped blocks: $(join((a.label for a in materialized.artifacts), ", "))")
        return requested
    end
    isempty(materialized.artifacts) && error("no block artifacts were dumped")
    idx = argmax([a.rows_full for a in materialized.artifacts])
    return materialized.artifacts[idx].label
end

function main(argv = ARGS)
    options = parse_options(argv)
    for child in ("blocks", "solve", "plots")
        rm(joinpath(options.output_dir, child); force = true, recursive = true)
    end
    mkpath(options.output_dir)
    mkpath(joinpath(options.output_dir, "solve"))
    mkpath(joinpath(options.output_dir, "plots"))

    @printf("%-48s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h2_pqg_moment_problem(options)
    materialized = dump_block_operators(data, options)
    diagnostic_label = choose_diagnostic_block(materialized, options.diagnostic_block)
    diagnostic = only(filter(a -> a.label == diagnostic_label, materialized.artifacts))
    A_full = let d = NPZ.npzread(joinpath(diagnostic.dir, "A_full.npz"))
        sparse(Int.(d["row"]) .+ 1, Int.(d["col"]) .+ 1, Float64.(d["data"]), Int(d["shape"][1]), Int(d["shape"][2]))
    end
    A_eq = let d = NPZ.npzread(joinpath(diagnostic.dir, "A_eq.npz"))
        sparse(Int.(d["row"]) .+ 1, Int.(d["col"]) .+ 1, Float64.(d["data"]), Int(d["shape"][1]), Int(d["shape"][2]))
    end
    write_diagnostics!(diagnostic.dir, A_full, A_eq, options.rank_tols)

    print_summary(data, options, build_seconds, materialized, diagnostic_label)

    solve_data = solve_with_bpsdp!(data, options, materialized)
    @printf("%-48s %.3f s\n", "BPSDP solve walltime", solve_data.solve_seconds)
    @printf("%-48s %s\n", "termination", string(termination_status(solve_data.jump_data.model)))
    @printf("%-48s %.12f Ha\n", "objective", objective_value(solve_data.jump_data.model))
    println("Wrote solve diagnostics under ", joinpath(options.output_dir, "solve"))
    return (; data, materialized, solve_data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
