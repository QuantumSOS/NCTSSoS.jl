#!/usr/bin/env julia

# No-solver small-N comparison between the generic Pauli charge/spatial/sign
# symmetry-preparation path and the specialized translation DFT path.
#
# Usage:
#   NCTS_COMPARE_NS=8,10,12 julia --project=. --startup-file=no perf/pauli_translation_compare.jl
#
# Safety:
#   By default this refuses N > 13.  Large comparisons need
#   NCTS_COMPARE_ALLOW_LARGE=true plus explicit wall/RSS estimates and safe
#   load/memory telemetry; otherwise downsize or delete the candidate.

Base.include(@__MODULE__, joinpath(@__DIR__, "shared_load_guard.jl"))

function _translation_compare_preimport_large_run_pressure_guard()
    _ncts_load_guard_parse_bool("NCTS_COMPARE_ALLOW_LARGE", false) ||
        return nothing
    ns_label = strip(get(ENV, "NCTS_COMPARE_NS", "8,10,12"))
    _ncts_check_large_run_pressure_guard(
        env_prefix="NCTS_COMPARE",
        label="translation comparison Ns=$ns_label",
    )
    return nothing
end

_translation_compare_preimport_large_run_pressure_guard()

using Dates
using Printf
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

function _timed(f)
    GC.gc(true)
    stats = @timed f()
    return stats.value, stats
end

function _stage_timed(stages::Vector{Pair{Symbol,Float64}}, stage::Symbol, f)
    value, stats = _timed(f)
    push!(stages, stage => stats.time)
    return value
end

function _parse_bool_env(name::AbstractString, default::Bool)
    raw = lowercase(strip(get(ENV, name, string(default))))
    raw in ("true", "1", "yes", "y") && return true
    raw in ("false", "0", "no", "n") && return false
    throw(ArgumentError("$name must be a boolean value, got $(repr(raw))."))
end

function _parse_int_list_env(name::AbstractString, default::AbstractString)
    values = Int[]
    for raw in split(get(ENV, name, default), ",")
        item = strip(raw)
        isempty(item) && continue
        push!(values, parse(Int, item))
    end
    isempty(values) && throw(ArgumentError("$name did not contain any integer sizes."))
    return values
end

function _check_size_guard(ns::AbstractVector{Int})
    allow_large = _parse_bool_env("NCTS_COMPARE_ALLOW_LARGE", false)
    max_n = parse(Int, get(ENV, "NCTS_COMPARE_MAX_N", "13"))
    if allow_large
        _ncts_check_large_run_pressure_guard(
            env_prefix="NCTS_COMPARE",
            label="translation comparison Ns=$(join(ns, ","))",
        )
        return nothing
    end
    for n in ns
        n <= max_n || throw(ArgumentError(
            "Refusing N=$n because NCTS_COMPARE_ALLOW_LARGE=false and " *
            "NCTS_COMPARE_MAX_N=$max_n.  Set NCTS_COMPARE_ALLOW_LARGE=true " *
            "only for intentional structural runs."
        ))
    end
    return nothing
end

_block_variable_count(sizes) = sum(n * (n + 1) ÷ 2 for n in sizes)

function _generic_charge_spatial_sign_summary(
    n::Integer;
    order::Integer,
    reflection_symmetry::Bool,
)
    stages = Pair{Symbol,Float64}[]
    registry, basis = _stage_timed(stages, :basis, () -> begin
        registry = first(create_pauli_variables(1:Int(n)))
        basis = pauli_contiguous_chain_basis(registry, Int(order))
        return registry, basis
    end)
    T = eltype(basis[1].word)
    generators = _stage_timed(stages, :generators, () -> begin
        translation = pauli_site_permutation([2:Int(n); 1])
        sign = pauli_sign_symmetry(Int(n); integer_type=T)
        if reflection_symmetry
            reflection = pauli_site_permutation(reverse(1:Int(n)))
            return [translation, reflection, sign]
        end
        return [translation, sign]
    end)
    group = _stage_timed(stages, :group_closure, () ->
        CliffordSymmetryGroup(generators; nqubits=Int(n), integer_type=T)
    )
    charge_groups = _stage_timed(stages, :charge_wedderburn_blocks, () ->
        NCTSSoS._pauli_charge_transform_groups(
            basis,
            PauliChargeSectorSpec(nqubits=Int(n), max_degree=Int(order)),
            group,
        )
    )
    block_sizes = sort!(Int[size(block.row_basis, 1) for group in charge_groups for block in group])
    return (
        basis_size=length(basis),
        group_order=length(group),
        block_count=length(block_sizes),
        max_block=maximum(block_sizes),
        scalar_variables=_block_variable_count(block_sizes),
        largest_blocks=block_sizes[max(1, length(block_sizes) - 9):end],
        stage_times=stages,
    )
end

function _specialized_translation_summary(
    n::Integer;
    order::Integer,
    reflection_symmetry::Bool,
)
    registry, ops = create_pauli_variables(1:Int(n))
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    _, report = pauli_translation_invariant_moment_relaxation(
        pop,
        ops,
        Int(order);
        reflection_symmetry,
    )
    metrics = translation_report_metrics(report)
    return (
        basis_size=metrics.basis_size,
        orbit_basis_size=metrics.orbit_basis_size,
        axis_orbit_closed=metrics.axis_orbit_closed,
        axis_orbit_basis_size=metrics.axis_orbit_basis_size,
        axis_reduction_ratio=metrics.axis_reduction_ratio,
        axis_orbit_size_histogram=metrics.axis_orbit_size_histogram,
        block_count=metrics.n_blocks,
        logical_max_block=metrics.logical_max_block,
        psd_max_block=metrics.psd_max_block,
        psd_scalar_variables=metrics.psd_symmetric_entries,
        product_cache_hit_rate=metrics.product_cache_hit_rate,
        construction_stage_time_seconds=metrics.construction_stage_time_seconds,
        largest_logical_blocks=last(sort(report.logical_block_sizes), min(10, length(report.logical_block_sizes))),
        largest_psd_blocks=last(sort(report.psd_block_sizes), min(10, length(report.psd_block_sizes))),
    )
end

function _print_stage_table(title::AbstractString, stages)
    println("\n#### $title")
    println("\n| stage | wall time (s) |")
    println("|:--|--:|")
    for (stage, seconds) in stages
        @printf("| `%s` | %.9f |\n", string(stage), seconds)
    end
end

function _sorted_stage_pairs(stage_times::AbstractDict)
    return [stage => stage_times[stage] for stage in sort!(collect(keys(stage_times)); by=string)]
end

function _print_timing_row(path, stats, summary)
    @printf(
        "| `%s` | %.6f | %s | %.6f | %d | %d | %d |\n",
        path,
        stats.time,
        _fmt_bytes(stats.bytes),
        stats.gctime,
        summary.block_count,
        get(summary, :max_block, get(summary, :logical_max_block, 0)),
        get(summary, :scalar_variables, get(summary, :psd_scalar_variables, 0)),
    )
end

function _run_one(n::Integer; order::Integer, reflection_symmetry::Bool)
    println("\n## N = $n, order = $order")
    println("\n| path | wall time (s) | allocated | GC time (s) | blocks | max logical block | symmetric scalar vars |")
    println("|:--|--:|--:|--:|--:|--:|--:|")

    generic, generic_stats = _timed(() -> _generic_charge_spatial_sign_summary(
        n;
        order,
        reflection_symmetry,
    ))
    specialized, specialized_stats = _timed(() -> _specialized_translation_summary(
        n;
        order,
        reflection_symmetry,
    ))
    _print_timing_row("generic charge/spatial/sign", generic_stats, generic)
    _print_timing_row("specialized translation DFT", specialized_stats, specialized)

    println("\n#### Notes")
    println("\n- generic group order: `$(generic.group_order)`")
    _print_stage_table("Generic stage attribution", generic.stage_times)
    _print_stage_table(
        "Specialized translation stage attribution",
        _sorted_stage_pairs(specialized.construction_stage_time_seconds),
    )
    println("- generic largest blocks: `$(generic.largest_blocks)`")
    println("- specialized orbit basis size: `$(specialized.orbit_basis_size)`")
    println("- specialized axis-orbit closed: `$(specialized.axis_orbit_closed)`")
    println("- specialized axis-orbit basis size: `$(specialized.axis_orbit_basis_size)`")
    println("- specialized axis-orbit reduction ratio: `$(@sprintf("%.3f", specialized.axis_reduction_ratio))x`")
    println("- specialized axis-orbit size histogram: `$(specialized.axis_orbit_size_histogram)`")
    println("- specialized solver-facing max block: `$(specialized.psd_max_block)`")
    println("- specialized largest logical blocks: `$(specialized.largest_logical_blocks)`")
    println("- specialized largest solver-facing blocks: `$(specialized.largest_psd_blocks)`")
    println("- specialized product-cache hit rate: `$(@sprintf("%.6f", specialized.product_cache_hit_rate))`")
    @printf(
        "- wall-time ratio, generic/specialized: `%.3fx`\n",
        generic_stats.time / max(specialized_stats.time, eps()),
    )
end

function main()
    ns = _parse_int_list_env("NCTS_COMPARE_NS", "8")
    _check_size_guard(ns)

    order = parse(Int, get(ENV, "NCTS_COMPARE_ORDER", "4"))
    reflection_symmetry = _parse_bool_env("NCTS_COMPARE_REFLECTION", true)

    println("# Pauli translation small-N structural comparison")
    println("\n- generated: `$(Dates.now())`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- solver calls: none")
    println("- reflection_symmetry: `$reflection_symmetry`")
    println("\nThis is a structural timing comparison.  The generic baseline uses charge")
    println("sector splitting plus a finite spatial/sign group; the specialized path")
    println("uses the Pauli-chain translation DFT backend.  Compare timings and block")
    println("geometry, not bit-for-bit model equivalence.")

    warm_n = min(4, minimum(ns))
    warm_order = min(order, max(0, warm_n - 1))
    println("\n## Warmup N = $warm_n")
    _timed(() -> _generic_charge_spatial_sign_summary(
        warm_n;
        order=warm_order,
        reflection_symmetry,
    ))
    _timed(() -> _specialized_translation_summary(
        warm_n;
        order=warm_order,
        reflection_symmetry,
    ))

    for n in ns
        _run_one(n; order, reflection_symmetry)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
