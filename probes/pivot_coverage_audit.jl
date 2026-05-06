#!/usr/bin/env julia
# Pivot-coverage audit for MomentProblem linear-cache design.
#
# This is a probe, not a solver. It builds representative MomentProblems,
# runs the existing lowering-time pivot discovery, and records how many
# canonical moment keys would become free/orphan keys in the planned cache.

using Dates
using Printf
using NCTSSoS

module PeriodicPQG
include(joinpath(@__DIR__, "..", "demos", "h4_periodic_moment_sos.jl"))
end

const DEFAULT_OUTPUT = normpath(joinpath(@__DIR__, "..", "output", "phase2", "pivot_coverage_audit.md"))
const H2_INTEGRALS = normpath(joinpath(@__DIR__, "..", "test", "data", "assets", "h2_chain_nk2_active_2e4o_integrals.txt"))

struct AuditOptions
    output::String
    only::Set{String}
end

function parse_options(argv)
    output = DEFAULT_OUTPUT
    only = Set{String}()
    for arg in argv
        if startswith(arg, "--output=")
            output = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--only=")
            only = Set(strip.(split(split(arg, "=", limit = 2)[2], ",")))
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. probes/pivot_coverage_audit.jl [--output=PATH] [--only=case1,case2]")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    return AuditOptions(output, only)
end

function build_from_pop(pop; order::Int, cs_algo=NoElimination(), ts_algo=NoElimination())
    config = SolverConfig(
        optimizer = nothing,
        order = order,
        cs_algo = cs_algo,
        ts_algo = ts_algo,
    )
    sparsity = compute_sparsity(pop, config)
    return NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
end

function build_periodic_pqg_mp(; norb::Int, nelec_per_cell::Int, integrals_path::AbstractString,
        blocking::Symbol=:spin, include_one_d::Bool=false,
        spin_resolved_trace::Bool=true, singlet_s2::Bool=true)
    nk = 2
    total_electrons = nk * nelec_per_cell
    h1e, eri = PeriodicPQG.load_integrals_txt(integrals_path; norb)
    registry, vars, ham = PeriodicPQG.build_h4_nk2_hamiltonian(h1e, eri; nk, norb)
    spin_orbitals = PeriodicPQG.build_spin_orbitals(vars; nk, norb)

    meq_constraints = typeof(ham)[
        PeriodicPQG.total_electron_constraint(vars, ham; norb, total_electrons),
    ]
    append!(meq_constraints, PeriodicPQG.trace_constraints(spin_orbitals, ham;
        total_electrons,
        include_total_d = !spin_resolved_trace,
    ))
    spin_resolved_trace && append!(meq_constraints,
        PeriodicPQG.spin_resolved_d_trace_constraints(spin_orbitals, ham; total_electrons))
    singlet_s2 && push!(meq_constraints,
        PeriodicPQG.singlet_s2_constraint(vars, ham; nk, norb, total_electrons))

    pop = polyopt(ham, registry; moment_eq_constraints = meq_constraints)
    row_bases = PeriodicPQG.build_rdm_bases(spin_orbitals)
    NM = eltype(row_bases.oneD_basis)
    identity = one(NM)

    oneD_blocks = include_one_d ?
        PeriodicPQG.grouped_blocks(row_bases.oneD_basis, blocking; nk, norb) :
        Vector{NM}[]
    twoD_blocks = PeriodicPQG.grouped_blocks(row_bases.twoD_basis, blocking; nk, norb)
    twoQ_blocks = PeriodicPQG.grouped_blocks(row_bases.twoQ_basis, blocking; nk, norb)
    twoG_blocks = PeriodicPQG.grouped_blocks_with_identity(row_bases.twoG_basis, identity, blocking; nk, norb)

    moment_support_basis = include_one_d ? NM[
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

    full_clique = PeriodicPQG.single_full_system_clique(registry, moment_support_basis, ham; nk, norb)
    term_sparsity = NCTSSoS.TermSparsity{NM}(copy(moment_support_basis), schouten_blocks)
    return NCTSSoS.moment_relax(pop, full_clique, [[term_sparsity]])
end

function build_h2_chain_nk2()
    return build_periodic_pqg_mp(
        norb = 4,
        nelec_per_cell = 2,
        integrals_path = H2_INTEGRALS,
        blocking = :spin,
        include_one_d = false,
        spin_resolved_trace = true,
        singlet_s2 = true,
    )
end

function build_h4_chain_nk2_proxy_small()
    # Same 2-k, 4-orbital active shape as the H2 asset, but half-filled like H4.
    # This is deliberately a pivot-coverage proxy, not a physical H4 energy model.
    return build_periodic_pqg_mp(
        norb = 4,
        nelec_per_cell = 4,
        integrals_path = H2_INTEGRALS,
        blocking = :spin,
        include_one_d = false,
        spin_resolved_trace = true,
        singlet_s2 = true,
    )
end

function build_pauli_3qubit_ground_state()
    registry, (x, y, z) = create_pauli_variables(1:3)
    ham = sum(0.25 * op[i] * op[mod1(i + 1, 3)] for op in (x, y, z) for i in 1:3)
    return build_from_pop(polyopt(ham, registry); order = 2)
end

function build_fermionic_4site_hubbard()
    nsites = 4
    registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
        ("c_up", 1:nsites),
        ("c_dn", 1:nsites),
    ])
    bonds = [(i, mod1(i + 1, nsites)) for i in 1:nsites]
    hopping = -1.0 * sum(
        c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
        c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
        for (i, j) in bonds
    )
    interaction = 4.0 * sum(
        (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
        for i in 1:nsites
    )
    ham = hopping + interaction
    n_up_total = 1.0 * sum(c_up_dag[i] * c_up[i] for i in 1:nsites)
    n_dn_total = 1.0 * sum(c_dn_dag[i] * c_dn[i] for i in 1:nsites)
    pop = polyopt(ham, registry; moment_eq_constraints = [
        n_up_total - 2.0 * one(ham),
        n_dn_total - 2.0 * one(ham),
    ])
    return build_from_pop(pop; order = 2)
end

function build_bosonic_2mode_truncated()
    registry, (b, b_dag) = create_bosonic_variables(1:2)
    n = [b_dag[i] * b[i] for i in 1:2]
    hopping = -1.0 * (b_dag[1] * b[2] + b_dag[2] * b[1])
    interaction = sum(n[i] * n[i] - n[i] for i in 1:2)
    ham = hopping + interaction
    pop = polyopt(ham, registry; moment_eq_constraints = [sum(n) - 2.0 * one(ham)])
    return build_from_pop(pop; order = 2)
end

function build_cs_failure_E10()
    nvars = 10
    registry, (x,) = create_noncommutative_variables([("X", 1:nvars)])
    poly_type = typeof(1.0 * x[1])
    objective = zero(poly_type)
    for i in 1:nvars
        a = x[i]
        b = x[mod1(i + 1, nvars)]
        c = x[mod1(i + 3, nvars)]
        term = a + 0.35 * b * c - 0.2 * c * b
        objective += term' * term
    end
    constraints = [1.0 - x[i]^2 for i in 1:nvars]
    pop = polyopt(objective, registry; ineq_constraints = constraints)
    return build_from_pop(pop; order = 2, cs_algo = MF(), ts_algo = MMD())
end

const CASES = [
    ("h2_chain_nk2", build_h2_chain_nk2, "primary H2/Nk=2 PQG BPSDP target; spin blocks + paper spin constraints"),
    ("h4_chain_nk2_proxy_small", build_h4_chain_nk2_proxy_small, "2-k, 4-orbital half-filled PQG proxy for H4; structural audit only"),
    ("pauli_3qubit_ground_state", build_pauli_3qubit_ground_state, "3-qubit periodic Heisenberg, order 2"),
    ("fermionic_4site_hubbard", build_fermionic_4site_hubbard, "4-site periodic Hubbard, half-filled N_up=N_dn=2, order 2"),
    ("bosonic_2mode_truncated", build_bosonic_2mode_truncated, "2-mode Bose-Hubbard with N=2 moment-equality truncation, order 2"),
    ("cs_failure_E10", build_cs_failure_E10, "deterministic 10-variable NC polyball stressor; historical E10 coefficients are unavailable"),
]

function polynomial_keys(poly)
    keys = Any[]
    seen = Set{Any}()
    simplified = simplify(poly)
    for (coef, mono) in simplified
        iszero(coef) && continue
        key = symmetric_canon(NCTSSoS.expval(mono))
        if !(key in seen)
            push!(seen, key)
            push!(keys, key)
        end
    end
    return keys
end

function key_presence(mp)
    presence = Dict{Any,Set{Symbol}}()
    function mark!(key, bucket)
        push!(get!(presence, key, Set{Symbol}()), bucket)
        return nothing
    end
    for key in polynomial_keys(mp.objective)
        mark!(key, :objective)
    end
    for (cone, mat) in mp.constraints
        bucket = cone in (:PSD, :HPSD) ? :psd_or_hpsd : cone == :Zero ? :zero : Symbol(lowercase(String(cone)))
        for entry in mat, key in polynomial_keys(entry)
            mark!(key, bucket)
        end
    end
    return presence
end

function orphan_category(n_orphans::Int, n_moments::Int)
    n_orphans == 0 && return "No orphans"
    ratio = n_moments == 0 ? 0.0 : n_orphans / n_moments
    return ratio < 0.001 ? "Few orphans (< 0.1%)" : "Many orphans"
end

function orphan_reason(presence::Set{Symbol})
    if :psd_or_hpsd in presence
        return "PSD/HPSD non-pivot entry"
    elseif :zero in presence && :objective in presence
        return "zero+objective only"
    elseif :zero in presence
        return "zero constraints only"
    elseif :objective in presence
        return "objective only"
    else
        return "not found by occurrence scan"
    end
end

function reason_counts(mp, orphans)
    presence = key_presence(mp)
    counts = Dict{String,Int}()
    for key in orphans
        reason = orphan_reason(get(presence, key, Set{Symbol}()))
        counts[reason] = get(counts, reason, 0) + 1
    end
    return counts
end

key_string(key) = sprint(show, key)

function sorted_key_sample(keys; limit::Int=12)
    strs = sort!(key_string.(keys))
    return strs[1:min(length(strs), limit)]
end

function audit_case(name::String, builder, note::String)
    GC.gc()
    build_seconds = @elapsed mp = builder()
    pivots = NCTSSoS._discover_pivots_unchecked(mp)
    orphans = NCTSSoS.orphan_keys(mp, pivots)
    all_keys = NCTSSoS._all_moment_keys(mp)
    psd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone in (:PSD, :HPSD)]
    zero_constraints = count(c -> c[1] == :Zero, mp.constraints)
    n_moments = length(all_keys)
    n_orphans = length(orphans)
    pct = n_moments == 0 ? 0.0 : 100.0 * n_orphans / n_moments
    return (; name, note, build_seconds, n_moments, n_pivots = length(pivots),
        n_orphans, orphan_pct = pct, category = orphan_category(n_orphans, n_moments),
        n_constraints = length(mp.constraints), n_psd_blocks = length(psd_sizes),
        max_psd_block = isempty(psd_sizes) ? 0 : maximum(psd_sizes),
        zero_constraints, reason_counts = reason_counts(mp, orphans),
        sample_orphans = sorted_key_sample(orphans))
end

function git_commit()
    env_commit = get(ENV, "GIT_COMMIT", "")
    !isempty(env_commit) && return env_commit
    try
        value = read(pipeline(`git rev-parse --short HEAD`; stderr = devnull), String)
        value = strip(value)
        return isempty(value) ? "unknown" : value
    catch
        return "unknown"
    end
end

function write_report(path::AbstractString, records)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Pivot-coverage audit")
        println(io)
        println(io, "Generated: ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
        println(io, "Git commit: `", git_commit(), "`")
        println(io)
        println(io, "Purpose: measure whether existing runtime pivot discovery covers every canonical moment key, or whether planned `MomentLinearData.free_keys` must carry orphan moments.")
        println(io)
        println(io, "| case | build s | moments | pivots | orphans | orphan % | category | PSD/HPSD blocks | max block | zero constraints |")
        println(io, "|---|---:|---:|---:|---:|---:|---|---:|---:|---:|")
        for r in records
            @printf(io, "| `%s` | %.3f | %d | %d | %d | %.4f | %s | %d | %d | %d |\n",
                r.name, r.build_seconds, r.n_moments, r.n_pivots, r.n_orphans,
                r.orphan_pct, r.category, r.n_psd_blocks, r.max_psd_block, r.zero_constraints)
        end
        println(io)
        println(io, "## Details")
        for r in records
            println(io)
            println(io, "### `", r.name, "`")
            println(io)
            println(io, r.note)
            println(io)
            println(io, "- constraints: ", r.n_constraints)
            println(io, "- orphan category: ", r.category)
            println(io, "- orphan reason counts:")
            if isempty(r.reason_counts)
                println(io, "  - none")
            else
                for key in sort!(collect(keys(r.reason_counts)))
                    println(io, "  - ", key, ": ", r.reason_counts[key])
                end
            end
            if !isempty(r.sample_orphans)
                println(io, "- sample orphan keys:")
                for key in r.sample_orphans
                    println(io, "  - `", replace(key, "`" => "\\`"), "`")
                end
            end
        end
        println(io)
        println(io, "## Gate decision")
        if all(r -> r.n_orphans == 0, records)
            println(io, "All audited cases have no orphans. Production lowering can treat `free_keys` as empty for these workloads.")
        elseif all(r -> r.orphan_pct < 0.1, records)
            println(io, "Every audited case is below the 0.1% orphan threshold. `:free_variables` remains the right default orphan policy; `:aux_psd_free` should stay compatibility-only.")
        else
            println(io, "At least one audited case has many orphans. Keep explicit `free_keys`; do not pretend pivot coverage is universal.")
        end
    end
    return path
end

function main(argv = ARGS)
    options = parse_options(argv)
    records = Any[]
    for (name, builder, note) in CASES
        (!isempty(options.only) && !(name in options.only)) && continue
        println("auditing ", name, " ...")
        push!(records, audit_case(name, builder, note))
    end
    isempty(records) && error("no cases selected")
    output = write_report(options.output, records)
    println("wrote ", output)
    return records
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
