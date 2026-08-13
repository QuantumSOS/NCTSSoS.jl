#!/usr/bin/env julia

# Phase-boundary memory profile of the translation-invariant (DFT) Pauli-chain
# pipeline: symbolic moment relaxation -> sos_dualize (JuMP model build) ->
# optional symbolic-data drop -> Mosek solve.
#
# For each N it records, at every phase boundary:
#   - wall time and Julia @allocated within the phase
#   - VmRSS after a full GC (live memory) and VmHWM (process peak, monotonic)
#   - Base.summarysize of the main retained structures
#
# The "drop symbolic" phase releases the MomentProblem (mp.constraints +
# mp.linear caches) after the JuMP model has been built, quantifying how much
# of the peak is avoidable retention versus JuMP/Mosek working memory.
#
# Usage:  julia perf/ti_memory_profile.jl [N1 N2 ...]   (default: 16 24 32)
# Env:    NCTS_MOSEK_THREADS (default cpu/4), NCTS_MOSEK_LOG (default 0),
#         NCTS_PROFILE_ORDER (default 4), NCTS_HWM_BUDGET_GIB (default 24;
#         escalation to the next N stops once VmHWM exceeds this).

using Pkg
const REPO = dirname(@__DIR__)
Pkg.activate(joinpath(homedir(), ".nctssos_ti_mosek_env"))
Pkg.develop(path=REPO)
Pkg.add(["MosekTools", "JuMP"])
Pkg.instantiate()

using NCTSSoS, JuMP, MosekTools, Printf, Dates
using NCTSSoS: pauli_translation_invariant_moment_relaxation, sos_dualize,
    heisenberg_chain_hamiltonian, create_pauli_variables, polyopt

# ---------------------------------------------------------------------------
# /proc helpers
# ---------------------------------------------------------------------------

function _proc_status_kib(field::AbstractString)
    for line in eachline("/proc/self/status")
        if startswith(line, field * ":")
            return parse(Int, split(line)[2])  # kB
        end
    end
    return -1
end

vmrss_gib() = _proc_status_kib("VmRSS") / (1024^2)
vmhwm_gib() = _proc_status_kib("VmHWM") / (1024^2)

gib(bytes) = bytes / (1024^3)

function full_gc()
    GC.gc(true)
    GC.gc(true)
end

mutable struct PhaseLog
    rows::Vector{NamedTuple}
end
PhaseLog() = PhaseLog(NamedTuple[])

function record!(
    log::PhaseLog,
    n::Int,
    phase::String,
    wall::Float64,
    alloc_bytes::Int;
    extra...,
)
    full_gc()
    row = (; n, phase, wall_s=round(wall; digits=2),
        alloc_gib=round(gib(alloc_bytes); digits=3),
        rss_gib=round(vmrss_gib(); digits=3),
        hwm_gib=round(vmhwm_gib(); digits=3),
        extra...)
    push!(log.rows, row)
    @printf(
        "PHASE N=%-4d %-16s wall=%8.2fs alloc=%8.3fGiB rss=%7.3fGiB hwm=%7.3fGiB",
        n, phase, row.wall_s, row.alloc_gib, row.rss_gib, row.hwm_gib,
    )
    for (k, v) in pairs((; extra...))
        print("  $k=$v")
    end
    println()
    flush(stdout)
    return row
end

function mosek_optimizer()
    threads = parse(Int, get(ENV, "NCTS_MOSEK_THREADS", string(max(1, div(Sys.CPU_THREADS, 4)))))
    return optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_LOG" => parse(Int, get(ENV, "NCTS_MOSEK_LOG", "0")),
        "MSK_IPAR_NUM_THREADS" => threads,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8,
    )
end

# ---------------------------------------------------------------------------
# Profiled run
# ---------------------------------------------------------------------------

function profile_instance!(log::PhaseLog, n::Int, order::Int)
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    full_gc()

    # Phase 1: symbolic moment relaxation (mp.constraints + mp.linear caches)
    local mp, report
    a0 = Base.gc_bytes()
    t0 = time()
    mp, report = pauli_translation_invariant_moment_relaxation(pop, ops, order)
    wall = time() - t0
    alloc = Base.gc_bytes() - a0
    sz_constraints = Base.summarysize(mp.constraints)
    sz_linear = Base.summarysize(mp.linear)
    sz_psd_lin = Base.summarysize(mp.linear.psd_blocks_lin)
    sz_mp = Base.summarysize(mp)
    record!(log, n, "build_mp", wall, alloc;
        n_blocks=length(report.psd_block_sizes),
        max_block=maximum(report.psd_block_sizes),
        n_moments=report.n_unique_moment_matrix_elements,
        constraints_gib=round(gib(sz_constraints); digits=3),
        linear_gib=round(gib(sz_linear); digits=3),
        psd_lin_gib=round(gib(sz_psd_lin); digits=3),
        mp_gib=round(gib(sz_mp); digits=3))

    # Phase 2: SOS dualization (JuMP model construction)
    a0 = Base.gc_bytes()
    t0 = time()
    sos = sos_dualize(mp)
    wall = time() - t0
    alloc = Base.gc_bytes() - a0
    sz_model = Base.summarysize(sos.model)
    record!(log, n, "sos_dualize", wall, alloc;
        model_gib=round(gib(sz_model); digits=3))

    # Phase 3: drop symbolic data — measures avoidable retention through solve
    full_gc()
    rss_before_drop = vmrss_gib()
    mp = nothing
    report = nothing
    pop = nothing
    registry = nothing
    ops = nothing
    full_gc()
    rss_after_drop = vmrss_gib()
    record!(log, n, "drop_symbolic", 0.0, 0;
        rss_drop_gib=round(rss_before_drop - rss_after_drop; digits=3))

    # Phase 4: attach Mosek and solve
    a0 = Base.gc_bytes()
    t0 = time()
    set_optimizer(sos.model, mosek_optimizer())
    optimize!(sos.model)
    wall = time() - t0
    alloc = Base.gc_bytes() - a0
    obj = objective_value(sos.model)
    record!(log, n, "mosek_solve", wall, alloc;
        status=string(termination_status(sos.model)),
        bound=round(obj; digits=10),
        per_site=round(obj / n; digits=10))

    # Phase 5: release the JuMP model, keep only the solver-independent answer
    sos = nothing
    record!(log, n, "drop_model", 0.0, 0)
    return obj
end

function main()
    order = parse(Int, get(ENV, "NCTS_PROFILE_ORDER", "4"))
    budget = parse(Float64, get(ENV, "NCTS_HWM_BUDGET_GIB", "24"))
    ns = isempty(ARGS) ? [16, 24, 32] : parse.(Int, ARGS)

    println("ti_memory_profile order=$order ns=$ns budget=$(budget)GiB ",
        "julia=$(VERSION) threads=$(Sys.CPU_THREADS) date=$(now())")
    log = PhaseLog()
    for n in ns
        if vmhwm_gib() > budget
            println("STOP: VmHWM $(round(vmhwm_gib(); digits=2)) GiB exceeds budget $(budget) GiB; skipping N>=$(n)")
            break
        end
        profile_instance!(log, n, order)
    end

    println("\n=== summary (per phase) ===")
    for row in log.rows
        println(row)
    end
end

main()
