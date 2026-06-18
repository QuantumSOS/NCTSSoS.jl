#!/usr/bin/env julia

# No-solver benchmark for the Pauli U(1) charge + spatial Clifford + SU(2)
# singlet symbolic-preparation path.
#
# Usage:
#   NCTS_PERF_NS=8,12,16 julia --project --startup-file=no perf/pauli_charge_singlet_prep.jl
#
# This intentionally stops before optimization. It times:
#   1. compute_sparsity
#   2. moment_relax_symmetric
#   3. sos_dualize / JuMP model construction

using Dates
using Printf
using JuMP
using NCTSSoS

function _fmt_bytes(bytes::Integer)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value = Float64(bytes)
    unit = 1
    while value >= 1024 && unit < length(units)
        value /= 1024
        unit += 1
    end
    return @sprintf("%.2f %s", value, units[unit])
end

function _timed(label::AbstractString, f)
    GC.gc(true)
    stats = @timed f()
    @printf("| `%s` | %.6f | %s | %.6f |\n", label, stats.time, _fmt_bytes(stats.bytes), stats.gctime)
    flush(stdout)
    return stats.value
end

function _build_case(N::Integer; order::Integer=2, offblock_check::Symbol=:off, check_invariance::Bool=true)
    registry, (σx, σy, σz) = create_pauli_variables(1:N)
    H = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
        for op in (σx, σy, σz) for i in 1:N
    )
    pop = polyopt(H, registry)

    translation = pauli_site_permutation([2:N; 1])
    reflection = pauli_site_permutation(reverse(1:N))
    symmetry = SymmetrySpec(
        [translation, reflection];
        pauli_charge=PauliChargeSectorSpec(nqubits=Int(N)),
        pauli_singlet=PauliSingletConstraintSpec(nqubits=Int(N)),
        offblock_check=offblock_check,
        check_invariance=check_invariance,
    )
    cfg = SolverConfig(
        optimizer=nothing,
        order=Int(order),
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )
    return (; N=Int(N), pop, symmetry, cfg)
end

_block_variable_count(sizes) = sum(n * (n + 1) ÷ 2 for n in sizes)

function _run_one(N::Integer; offblock_check::Symbol=:off, check_invariance::Bool=true)
    println("\n## N = $N")
    println("\n| phase | wall time (s) | allocated | GC time (s) |")
    println("|:--|--:|--:|--:|")

    case = _timed("setup Heisenberg case", () -> _build_case(N; offblock_check, check_invariance))
    sparsity = _timed("compute_sparsity", () -> compute_sparsity(case.pop, case.cfg))

    basis = only(sparsity.corr_sparsity.clq_mom_mtx_bases)
    dense_dim = length(basis)
    println("\n- dense half-basis dimension: `$dense_dim`")
    println("- dense PSD scalar variables estimate: `$(dense_dim * (dense_dim + 1) ÷ 2)`")

    mp, report = _timed("moment_relax_symmetric", () -> NCTSSoS.moment_relax_symmetric(
        case.pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
        case.symmetry,
    ))
    sos = _timed("sos_dualize (JuMP model build, no solver)", () -> NCTSSoS.sos_dualize(mp))

    sizes = report.psd_block_sizes
    println("\n#### Symmetry report")
    println("\n- group order: `$(report.group_order)`")
    println("- PSD blocks: `$(length(sizes))`")
    println("- largest PSD block: `$(maximum(sizes))`")
    println("- reduced PSD scalar variables: `$(_block_variable_count(sizes))`")
    println("- invariant moment count: `$(report.invariant_moment_count)`")
    println("- symbolic constraints: `$(length(mp.constraints))`")
    println("- linear moments: `$(length(mp.linear.moments))`")
    println("- zero scalar constraints in linear cache: `$(length(mp.linear.zero_constraints))`")
    println("- JuMP variables: `$(JuMP.num_variables(sos.model))`")
end

function main()
    ns = parse.(Int, split(get(ENV, "NCTS_PERF_NS", "8,12,16"), ","))
    offblock_check = Symbol(get(ENV, "NCTS_PERF_OFFBLOCK_CHECK", "off"))
    check_invariance = parse(Bool, get(ENV, "NCTS_PERF_CHECK_INVARIANCE", "true"))

    println("# Pauli charge/spatial/singlet no-solver benchmark")
    println("\n- generated: `$(Dates.now())`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- offblock_check: `$offblock_check`")
    println("- solver calls: none")

    # Warm the compiler on the smallest nontrivial ring; keep reported N timings useful.
    println("\n## Warmup N = 4")
    println("\n| phase | wall time (s) | allocated | GC time (s) |")
    println("|:--|--:|--:|--:|")
    warm = _timed("setup Heisenberg case", () -> _build_case(4; offblock_check, check_invariance))
    warm_sparsity = _timed("compute_sparsity", () -> compute_sparsity(warm.pop, warm.cfg))
    warm_mp, _ = _timed("moment_relax_symmetric", () -> NCTSSoS.moment_relax_symmetric(
        warm.pop,
        warm_sparsity.corr_sparsity,
        warm_sparsity.cliques_term_sparsities,
        warm.symmetry,
    ))
    _timed("sos_dualize (JuMP model build, no solver)", () -> NCTSSoS.sos_dualize(warm_mp))

    for N in ns
        _run_one(N; offblock_check, check_invariance)
    end
end

main()
