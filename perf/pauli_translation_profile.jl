#!/usr/bin/env julia

# No-solver benchmark for the specialized translation-invariant Pauli-chain
# construction path.
#
# Usage:
#   NCTS_TRANSLATION_NS=8,12 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_DIRECT_LINEAR=true NCTS_TRANSLATION_NS=8,12 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_WARMUP=false NCTS_TRANSLATION_NS=8 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_REPEATS=3 NCTS_TRANSLATION_NS=8 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_SU2=true NCTS_TRANSLATION_MOMENT_EQ_H2=true NCTS_TRANSLATION_NS=6 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_SIGN=false NCTS_TRANSLATION_SU2=true NCTS_TRANSLATION_SINGLET_EQUALITIES=true NCTS_TRANSLATION_NS=6 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_SIGN=false NCTS_TRANSLATION_SU2=true NCTS_TRANSLATION_RDM_K=2 NCTS_TRANSLATION_RDM_DECOMPOSITION=su2 NCTS_TRANSLATION_NS=6 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_SIGN=false NCTS_TRANSLATION_SU2=true NCTS_TRANSLATION_PSO_WIDTH=1 NCTS_TRANSLATION_NS=6 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_DIRECT_LINEAR=true NCTS_TRANSLATION_SIGN=false NCTS_TRANSLATION_SU2=true NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true NCTS_TRANSLATION_SU2_MOMENT_QUOTIENT=true NCTS_TRANSLATION_RDM_DECOMPOSITION=su2 NCTS_TRANSLATION_RDM_SUPPORT=extend NCTS_TRANSLATION_NS=10 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_TARGET_ONLY=true NCTS_TRANSLATION_AXIS_EQUALITIES=true NCTS_TRANSLATION_NS=8 NCTS_TRANSLATION_ORDER=2 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_TARGET_ONLY=true NCTS_TRANSLATION_MOMENT_EQ_H2=true NCTS_TRANSLATION_NS=8 NCTS_TRANSLATION_ORDER=2 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_TARGET_ONLY=true NCTS_TRANSLATION_NS=12 NCTS_TRANSLATION_ORDER=4 julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#   NCTS_TRANSLATION_ALLOW_LARGE=true NCTS_TRANSLATION_TARGET_ONLY=true NCTS_TRANSLATION_NS=100 NCTS_TRANSLATION_ORDER=4 NCTS_TRANSLATION_ESTIMATED_WALL_SECONDS=<estimate-before-launch> NCTS_TRANSLATION_ESTIMATED_RSS_GIB=<estimate-before-launch> julia --project=. --startup-file=no perf/pauli_translation_profile.jl
#
# Safety:
#   By default this refuses N > 13.  For intentional structural benchmarks,
#   set NCTS_TRANSLATION_ALLOW_LARGE=true only with explicit wall/RSS
#   estimates and safe load/memory telemetry.  Target-only mode is checked too;
#   use the same pressure gate for large analytic shape reports.

Base.include(@__MODULE__, joinpath(@__DIR__, "shared_load_guard.jl"))

function _translation_profile_preimport_large_run_pressure_guard()
    _ncts_load_guard_parse_bool("NCTS_TRANSLATION_ALLOW_LARGE", false) ||
        return nothing
    ns_label = strip(get(ENV, "NCTS_TRANSLATION_NS", "8,12"))
    _ncts_check_large_run_pressure_guard(
        env_prefix="NCTS_TRANSLATION",
        label="translation profile Ns=$ns_label",
    )
    return nothing
end

_translation_profile_preimport_large_run_pressure_guard()

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

function _translation_profile_model_size_gate_reason(;
    status::AbstractString,
    estimated_rss_bytes::Integer,
    mem_available_bytes,
    max_rss_fraction::Real,
)
    status == "ok" && return ""
    status == "blocked_missing_memory_telemetry" &&
        return "available-memory telemetry is unavailable"
    status == "blocked_insufficient_memory" && return (
        "estimated dense-Schur RSS proxy $(_fmt_bytes(estimated_rss_bytes)) " *
        "exceeds available memory $(_fmt_bytes(mem_available_bytes))"
    )
    status == "blocked_insufficient_memory_headroom" && return (
        "estimated dense-Schur RSS proxy $(_fmt_bytes(estimated_rss_bytes)) " *
        "would consume more than $(Float64(max_rss_fraction) * 100)% of " *
        "available memory $(_fmt_bytes(mem_available_bytes))"
    )
    return "model-size gate returned $status"
end

function _translation_profile_model_size_gate(estimated_rss_bytes::Integer)
    mem_available_bytes = _ncts_mem_available_bytes()
    max_rss_fraction = _ncts_load_guard_parse_float(
        "NCTS_TRANSLATION_MODEL_SIZE_MAX_RSS_FRACTION",
        _ncts_load_guard_parse_float("NCTS_LOAD_GUARD_MAX_RSS_FRACTION", 0.8),
    )
    status = _ncts_memory_pressure_status(
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    reason = _translation_profile_model_size_gate_reason(
        status=status,
        estimated_rss_bytes=estimated_rss_bytes,
        mem_available_bytes=mem_available_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    return (
        status=status,
        reason=reason,
        estimated_rss_bytes=estimated_rss_bytes,
        mem_available_bytes=mem_available_bytes,
        max_rss_fraction=Float64(max_rss_fraction),
    )
end

function _translation_profile_model_size_gate_status(estimated_rss_bytes::Integer)
    return _translation_profile_model_size_gate(estimated_rss_bytes).status
end

function _proc_status_kb(field::AbstractString)
    path = "/proc/self/status"
    isfile(path) || return nothing
    prefix = field * ":"
    for line in eachline(path)
        startswith(line, prefix) || continue
        parts = split(line)
        length(parts) >= 2 || return nothing
        return parse(Int, parts[2])
    end
    return nothing
end

function _rss_string()
    rss = _proc_status_kb("VmRSS")
    hwm = _proc_status_kb("VmHWM")
    isnothing(rss) && isnothing(hwm) && return "n/a"
    rss_s = isnothing(rss) ? "n/a" : _fmt_bytes(rss * 1024)
    hwm_s = isnothing(hwm) ? "n/a" : _fmt_bytes(hwm * 1024)
    return "rss=$rss_s, hwm=$hwm_s"
end

function _timed(label::AbstractString, f)
    GC.gc(true)
    stats = @timed f()
    @printf(
        "| `%s` | %.6f | %s | %.6f | %s |\n",
        label,
        stats.time,
        _fmt_bytes(stats.bytes),
        stats.gctime,
        _rss_string(),
    )
    flush(stdout)
    return stats.value, stats
end

function _parse_bool_env(name::AbstractString, default::Bool)
    raw = lowercase(strip(get(ENV, name, string(default))))
    raw in ("true", "1", "yes", "y") && return true
    raw in ("false", "0", "no", "n") && return false
    throw(ArgumentError("$name must be a boolean value, got $(repr(raw))."))
end

function _parse_optional_int_env(name::AbstractString)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    lowercase(raw) in ("nothing", "none", "false", "0") && return nothing
    return parse(Int, raw)
end

function _parse_positive_int_env(name::AbstractString, default::Integer)
    value = parse(Int, get(ENV, name, string(default)))
    value >= 1 || throw(ArgumentError("$name must be >= 1, got $value."))
    return value
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

function _parse_int_tuple2_env(name::AbstractString, default::AbstractString)
    values = _parse_int_list_env(name, default)
    length(values) == 2 || throw(ArgumentError("$name must contain exactly two integers, got $values."))
    return (values[1], values[2])
end

function _check_size_guard(ns::AbstractVector{Int})
    allow_large = _parse_bool_env("NCTS_TRANSLATION_ALLOW_LARGE", false)
    max_n = parse(Int, get(ENV, "NCTS_TRANSLATION_MAX_N", "13"))
    if allow_large
        _ncts_check_large_run_pressure_guard(
            env_prefix="NCTS_TRANSLATION",
            label="translation profile Ns=$(join(ns, ","))",
        )
        return nothing
    end
    for n in ns
        n <= max_n || throw(ArgumentError(
            "Refusing N=$n because NCTS_TRANSLATION_ALLOW_LARGE=false and " *
            "NCTS_TRANSLATION_MAX_N=$max_n.  Set NCTS_TRANSLATION_ALLOW_LARGE=true " *
            "for intentional large structural runs."
        ))
    end
    return nothing
end

function _check_base_su2_extend_rdm_profile_guard(;
    target_only::Bool,
    direct_linear::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool=false,
    su2_symmetry::Bool,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    real_moment_matrix::Bool,
    reflection_symmetry::Bool,
)
    su2_moment_quotient && !base_su2_extend_rdm && throw(ArgumentError(
        "NCTS_TRANSLATION_SU2_MOMENT_QUOTIENT=true requires " *
        "NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true."
    ))
    su2_moment_quotient && !(target_only || direct_linear) && throw(ArgumentError(
        "NCTS_TRANSLATION_SU2_MOMENT_QUOTIENT=true requires " *
        "NCTS_TRANSLATION_DIRECT_LINEAR=true for constructed profiles."
    ))
    base_su2_extend_rdm || return nothing
    (target_only || direct_linear) || throw(ArgumentError(
        "NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true requires " *
        "NCTS_TRANSLATION_DIRECT_LINEAR=true for constructed profiles or " *
        "NCTS_TRANSLATION_TARGET_ONLY=true for structural reports."
    ))
    (
        su2_symmetry &&
        contiguous_rdm_decomposition == :su2 &&
        contiguous_rdm_support == :extend
    ) || throw(ArgumentError(
        "NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true requires " *
        "NCTS_TRANSLATION_SU2=true, NCTS_TRANSLATION_RDM_DECOMPOSITION=su2, " *
        "and NCTS_TRANSLATION_RDM_SUPPORT=extend."
    ))
    target_only && !real_moment_matrix && throw(ArgumentError(
        "NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true target-only reports require " *
        "NCTS_TRANSLATION_REAL=true."
    ))
    !target_only && !real_moment_matrix && reflection_symmetry && throw(ArgumentError(
        "NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true with " *
        "NCTS_TRANSLATION_REFLECTION=true requires NCTS_TRANSLATION_REAL=true."
    ))
    return nothing
end

function _check_singlet_channel_profile_guard(;
    singlet_channel_equalities::Bool,
    su2_symmetry::Bool,
)
    singlet_channel_equalities || return nothing
    su2_symmetry || throw(ArgumentError(
        "NCTS_TRANSLATION_SINGLET_EQUALITIES=true requires " *
        "NCTS_TRANSLATION_SU2=true."
    ))
    return nothing
end

function _heisenberg_case(n::Integer; moment_eq_h2::Bool=false)
    registry, ops = create_pauli_variables(1:Int(n))
    hamiltonian = heisenberg_chain_hamiltonian(ops)
    pop = moment_eq_h2 ?
        polyopt(hamiltonian, registry; moment_eq_constraints=[hamiltonian * hamiltonian]) :
        polyopt(hamiltonian, registry)
    return (; n=Int(n), registry, ops, pop)
end

function _print_stage_times(metrics)
    println("\n| construction stage | internal time (s) |")
    println("|:--|--:|")
    for stage in sort!(collect(keys(metrics.construction_stage_time_seconds)); by=string)
        @printf("| `%s` | %.9f |\n", string(stage), metrics.construction_stage_time_seconds[stage])
    end
end

function _print_histogram(
    label::AbstractString,
    histogram;
    key_label::AbstractString="block size",
)
    println("\n#### $label")
    println("\n| $key_label | count |")
    println("|:--|--:|")
    for (key, count) in histogram
        println("| `$key` | $count |")
    end
end

function _histogram_pairs(values)
    counts = Dict{Int,Int}()
    for value in values
        key = Int(value)
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by=first)
end

function _value_histogram_pairs(values)
    counts = Dict{Any,Int}()
    for value in values
        counts[value] = get(counts, value, 0) + 1
    end
    return sort!(collect(counts); by=pair -> string(first(pair)))
end

function _profile_repeat_shape(result)
    return (
        linear_moments=result.linear_moments,
        zero_constraints=result.zero_constraints,
        psd_blocks=result.psd_blocks,
        solver_facing_max_block=result.solver_facing_max_block,
        contiguous_rdm_blocks=result.contiguous_rdm_blocks,
        psd_state_opt_blocks=result.psd_state_opt_blocks,
    )
end

function _format_seconds_range(values)
    collected = collect(values)
    lo = minimum(collected)
    hi = maximum(collected)
    lo == hi && return @sprintf("%.9f", lo)
    return @sprintf("%.9f .. %.9f", lo, hi)
end

function _format_bytes_range(values)
    collected = collect(values)
    lo = minimum(collected)
    hi = maximum(collected)
    lo == hi && return _fmt_bytes(lo)
    return "$(_fmt_bytes(lo)) .. $(_fmt_bytes(hi))"
end

function _profile_repeat_stage_rows(results; limit::Integer=6)
    limit >= 1 || return NamedTuple[]
    rows = NamedTuple[]
    ns = sort!(unique(Int(result.n) for result in results))
    for n in ns
        group = [result for result in results if Int(result.n) == n]
        length(group) > 1 || continue

        stage_values = Dict{Symbol,Vector{Float64}}()
        for result in group
            :construction_stage_time_seconds in propertynames(result) || continue
            for (stage, seconds) in result.construction_stage_time_seconds
                stage_key = stage isa Symbol ? stage : Symbol(stage)
                push!(get!(stage_values, stage_key, Float64[]), Float64(seconds))
            end
        end
        isempty(stage_values) && continue

        ordered_stages = sort!(
            collect(stage_values);
            by=pair -> (-maximum(last(pair)), string(first(pair))),
        )
        for (index, (stage, values)) in enumerate(ordered_stages)
            index > limit && break
            push!(rows, (
                n=n,
                stage=stage,
                repeats=length(values),
                seconds_range=_format_seconds_range(values),
            ))
        end
    end
    return rows
end

function _print_repeat_stage_summary(results, io::IO=stdout; limit::Integer=6)
    rows = _profile_repeat_stage_rows(results; limit)
    isempty(rows) && return nothing

    println(io, "\n#### Repeat Stage Summary")
    println(io, "\n| N | construction stage | repeats | seconds range |")
    println(io, "|--:|:--|--:|--:|")
    for row in rows
        println(io, "| $(row.n) | $(row.stage) | $(row.repeats) | $(row.seconds_range) |")
    end
    return nothing
end

function _print_repeat_summary(results, io::IO=stdout; stage_limit::Integer=6)
    isempty(results) && return nothing
    ns = sort!(unique(Int(result.n) for result in results))
    any(count(result -> Int(result.n) == n, results) > 1 for n in ns) || return nothing

    println(io, "\n#### Repeat Summary")
    println(io, "\n| N | repeats | construction time range (s) | outer wall range (s) | allocation range | shape stable | shape |")
    println(io, "|--:|--:|--:|--:|:--|:--|:--|")
    for n in ns
        group = [result for result in results if Int(result.n) == n]
        length(group) > 1 || continue
        shapes = [_profile_repeat_shape(result) for result in group]
        shape_stable = length(unique(shapes)) == 1
        shape = first(shapes)
        shape_text = "linear=$(shape.linear_moments), zero=$(shape.zero_constraints), " *
            "psd=$(shape.psd_blocks), max=$(shape.solver_facing_max_block), " *
            "rdm=$(shape.contiguous_rdm_blocks), pso=$(shape.psd_state_opt_blocks)"
        println(
            io,
            "| $n | $(length(group)) | $(_format_seconds_range(result.construction_time_seconds for result in group)) | " *
            "$(_format_seconds_range(result.outer_relaxation_wall_time_seconds for result in group)) | " *
            "$(_format_bytes_range(result.allocated_bytes for result in group)) | $shape_stable | $shape_text |"
        )
    end
    _print_repeat_stage_summary(results, io; limit=stage_limit)
    return nothing
end

function _print_coefficient_domain_histograms(
    label_prefix::AbstractString,
    coefficient_domains,
    exact_coefficient_domains,
)
    isempty(coefficient_domains) || _print_histogram(
        "$label_prefix coefficient-domain histogram",
        _value_histogram_pairs(coefficient_domains);
        key_label="coefficient domain",
    )
    isempty(exact_coefficient_domains) || _print_histogram(
        "$label_prefix exact coefficient-domain histogram",
        _value_histogram_pairs(exact_coefficient_domains);
        key_label="exact coefficient domain",
    )
    return nothing
end

function _print_qmbcertify_rdm_references(references)
    isempty(references) && return nothing
    println("\n#### QMBCertify RDM reference blocks")
    println("\n| k | blocks | max block | total side | symmetric entries |")
    println("|--:|:--|--:|--:|--:|")
    for ref in references
        println(
            "| $(ref.k) | `$(ref.block_sizes)` | $(ref.max_block) | " *
            "$(ref.total_block_side) | $(ref.symmetric_entries) |",
        )
    end
    return nothing
end

function _print_qmbcertify_rdm_reference_aggregate(target)
    isempty(target.qmbcertify_rdm_reference_block_sizes) && return nothing
    println("\n#### QMBCertify RDM reference aggregate")
    println("\n| metric | value |")
    println("|:--|:--|")
    println("| blocks | `$(target.qmbcertify_rdm_reference_block_sizes)` |")
    println("| block count | $(target.qmbcertify_rdm_reference_n_blocks) |")
    println("| max block | $(target.qmbcertify_rdm_reference_max_block) |")
    println("| total side | $(target.qmbcertify_rdm_reference_total_block_side) |")
    println("| dense entries | $(target.qmbcertify_rdm_reference_dense_entries) |")
    println("| symmetric entries | $(target.qmbcertify_rdm_reference_symmetric_entries) |")
    println("| dense bytes | $(_fmt_bytes(target.qmbcertify_rdm_reference_dense_bytes)) |")
    println("| symmetric bytes | $(_fmt_bytes(target.qmbcertify_rdm_reference_symmetric_bytes)) |")
    println("| requires construction | $(target.qmbcertify_rdm_reference_requires_construction) |")
    return nothing
end

function _print_qmbcertify_base_reference(target)
    isnothing(target.qmbcertify_base_reference_extra) && return nothing
    println("\n#### QMBCertify base reference")
    println("\n| metric | value |")
    println("|:--|:--|")
    println("| extra | $(target.qmbcertify_base_reference_extra) |")
    println("| active construction | $(target.qmbcertify_base_reference_active) |")
    println("| requires construction | $(target.qmbcertify_base_reference_requires_construction) |")
    println("| family count by parity | `$(target.qmbcertify_base_reference_family_count_by_parity)` |")
    println("| blocks | `$(target.qmbcertify_base_reference_block_sizes)` |")
    println("| block count | $(target.qmbcertify_base_reference_n_blocks) |")
    println("| max block | $(target.qmbcertify_base_reference_max_block) |")
    println("| total side | $(target.qmbcertify_base_reference_total_block_side) |")
    println("| dense entries | $(target.qmbcertify_base_reference_dense_entries) |")
    println("| symmetric entries | $(target.qmbcertify_base_reference_symmetric_entries) |")
    println("| dense bytes | $(_fmt_bytes(target.qmbcertify_base_reference_dense_bytes)) |")
    println("| symmetric bytes | $(_fmt_bytes(target.qmbcertify_base_reference_symmetric_bytes)) |")
    println(
        "| support nonzero rows | " *
        "$(target.qmbcertify_base_reference_support_nonzero_row_count) |",
    )
    println(
        "| support zero rows | " *
        "$(target.qmbcertify_base_reference_support_zero_row_count) |",
    )
    println(
        "| support diagonal nonzero rows | " *
        "$(target.qmbcertify_base_reference_support_diagonal_nonzero_row_count) |",
    )
    println(
        "| support off-diagonal nonzero rows | " *
        "$(target.qmbcertify_base_reference_support_offdiagonal_nonzero_row_count) |",
    )
    println(
        "| support unique count | " *
        "$(target.qmbcertify_base_reference_support_unique_count) |",
    )
    println(
        "| support diagonal unique count | " *
        "$(target.qmbcertify_base_reference_support_diagonal_unique_count) |",
    )
    println(
        "| support off-diagonal unique count | " *
        "$(target.qmbcertify_base_reference_support_offdiagonal_unique_count) |",
    )
    _print_histogram(
        "QMBCertify base reference block histogram",
        target.qmbcertify_base_reference_block_histogram,
    )
    _print_histogram(
        "QMBCertify base reference support word-length histogram",
        target.qmbcertify_base_reference_support_word_length_histogram;
        key_label="word length",
    )
    return nothing
end

function _print_su2_accounting_records(records)
    isempty(records) && return nothing
    if length(records) > 20
        println(
            "\nSU(2) accounting record table omitted: $(length(records)) sector records.",
        )
        return nothing
    end

    println("\n#### SU(2) accounting records")
    println("\n| label | full dense | active dense | reduced dense | off-block | magnetic-copy |")
    println("|:--|--:|--:|--:|--:|--:|")
    for record in records
        println(
            "| `$(record.label)` | $(record.su2_full_dense_entries) | " *
            "$(record.su2_active_dense_entries) | " *
            "$(record.su2_reduced_dense_entries) | " *
            "$(record.offblock_entry_count) | $(record.copy_entry_count) |",
        )
    end
    return nothing
end

function _print_su2_structural_target(
    n::Integer,
    order::Integer;
    real_moment_matrix::Bool,
    reflection_symmetry::Bool,
)
    2 * Int(order) < Int(n) || return nothing
    orbit_target = pauli_su2_translation_orbit_structural_targets(
        n,
        order;
        real_moment_matrix,
        reflection_symmetry,
    )
    println("\n#### SU(2) translation-orbit structural target")
    println("\n| metric | value |")
    println("|:--|--:|")
    println("| construction performed | false |")
    println("| translation combined | $(orbit_target.translation_combined) |")
    println("| reflection combined | $(orbit_target.reflection_combined) |")
    println("| orbit basis size | $(orbit_target.orbit_basis_size) |")
    println("| momentum sectors | $(orbit_target.momentum_sector_count) |")
    println("| PSD blocks | $(orbit_target.n_blocks) |")
    println("| logical max block | $(orbit_target.logical_max_block) |")
    println("| solver-facing max block | $(orbit_target.psd_max_block) |")
    println("| logical total block side | $(orbit_target.logical_total_block_side) |")
    println("| solver-facing total block side | $(orbit_target.psd_total_block_side) |")
    println("| solver-facing symmetric entries | $(orbit_target.psd_symmetric_entries) |")
    println("| solver-facing symmetric bytes | $(_fmt_bytes(orbit_target.psd_symmetric_bytes)) |")
    println("| model-size gate status | `$(orbit_target.estimated_model_size_gate_status)` |")
    println("| model-size gate reason | $(orbit_target.estimated_model_size_gate_reason) |")
    println("| singlet channels | $(orbit_target.singlet_channel_count) |")
    println("| singlet channel support counts | `$(orbit_target.singlet_channel_support_counts)` |")
    println("| singlet-channel equality candidate rows | $(orbit_target.singlet_channel_equality_row_count) |")
    println("| SU(2) full dense entries | $(orbit_target.su2_full_dense_entries) |")
    println("| SU(2) active dense entries | $(orbit_target.su2_active_dense_entries) |")
    println("| SU(2) reduced dense entries | $(orbit_target.su2_reduced_dense_entries) |")
    println("| SU(2) off-block entries | $(orbit_target.offblock_entry_count) |")
    println("| SU(2) magnetic-copy entries | $(orbit_target.copy_entry_count) |")
    println("| SU(2) off-block zero-row budget | $(orbit_target.wigner_offblock_zero_row_budget) |")
    println("| SU(2) magnetic-copy zero-row budget | $(orbit_target.wigner_magnetic_copy_zero_row_budget) |")
    println("| SU(2) Wigner zero-row budget | $(orbit_target.wigner_zero_row_budget) |")
    println("| SU(2) accounting records | $(orbit_target.su2_accounting_record_count) |")
    println("| solve supported | $(orbit_target.solve_supported) |")
    println("| solve blocker | `$(orbit_target.solve_blocker)` |")
    println("| solve blocker reason | $(orbit_target.solve_blocker_reason) |")
    println("| solve unsupported block features | `$(orbit_target.solve_unsupported_block_features)` |")
    println("| solve unsupported zero-row features | `$(orbit_target.solve_unsupported_zero_features)` |")
    _print_su2_accounting_records(orbit_target.su2_accounting_records)
    _print_histogram(
        "SU(2) translation-orbit logical block histogram",
        orbit_target.logical_block_histogram,
    )
    _print_histogram(
        "SU(2) translation-orbit solver-facing block histogram",
        orbit_target.psd_block_histogram,
    )
    _print_coefficient_domain_histograms(
        "SU(2) structural",
        orbit_target.block_coefficient_domains,
        orbit_target.block_exact_coefficient_domains,
    )

    target = pauli_su2_contiguous_chain_structural_targets(n, order)
    println("\n#### SU(2) full-basis contiguous-chain structural target")
    println("\n| metric | value |")
    println("|:--|--:|")
    println("| construction performed | false |")
    println("| translation combined | $(target.translation_combined) |")
    println("| basis size | $(target.basis_size) |")
    println("| support counts | `$(target.support_counts)` |")
    println("| reduced PSD blocks | $(length(target.reduced_block_sizes)) |")
    println("| reduced max block | $(target.reduced_max_block) |")
    println("| reduced total block side | $(target.reduced_total_block_side) |")
    println("| full dense entries | $(target.full_dense_entries) |")
    println("| reduced dense entries | $(target.reduced_dense_entries) |")
    println("| reduced symmetric entries | $(target.reduced_symmetric_entries) |")
    println("| reduced symmetric bytes | $(_fmt_bytes(target.reduced_symmetric_bytes)) |")
    println("| model-size gate status | `$(target.estimated_model_size_gate_status)` |")
    println("| model-size gate reason | $(target.estimated_model_size_gate_reason) |")
    println("| singlet channels | $(target.singlet_channel_count) |")
    println("| singlet channel support counts | `$(target.singlet_channel_support_counts)` |")
    println("| singlet-channel equality candidate rows | $(target.singlet_channel_equality_row_count) |")
    println("| off-block entries | $(target.offblock_entry_count) |")
    println("| magnetic-copy entries | $(target.copy_entry_count) |")
    println("| solve supported | $(target.solve_supported) |")
    println("| solve blocker | `$(target.solve_blocker)` |")
    println("| solve blocker reason | $(target.solve_blocker_reason) |")
    println("| solve unsupported block features | `$(target.solve_unsupported_block_features)` |")
    println("| solve unsupported zero-row features | `$(target.solve_unsupported_zero_features)` |")
    _print_histogram(
        "SU(2) reduced block histogram",
        _histogram_pairs(target.reduced_block_sizes);
        key_label="block size",
    )
    _print_coefficient_domain_histograms(
        "SU(2) full-basis structural",
        target.block_coefficient_domains,
        target.block_exact_coefficient_domains,
    )
    return nothing
end

function _build_translation_profile(
    case,
    order::Integer;
    qmbcertify_base_construct::Bool,
    direct_linear::Bool,
    sign_symmetry::Bool,
    reflection_symmetry::Bool,
    axis_rotation_symmetry::Bool,
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    real_moment_matrix::Bool,
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool,
    singlet_channel_equalities::Bool,
    moment_eq_h2::Bool,
    linear_state_opt_width,
    psd_state_opt_width,
    qmbcertify_base_extra,
    qmbcertify_base_three_type::Tuple{Int,Int},
)
    if qmbcertify_base_construct
        if !isnothing(contiguous_rdm_k)
            contiguous_rdm_decomposition == :qmbcertify || throw(ArgumentError(
                "NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true with RDM requires " *
                "NCTS_TRANSLATION_RDM_DECOMPOSITION=qmbcertify."
            ))
            contiguous_rdm_support == :extend || throw(ArgumentError(
                "NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true with RDM requires " *
                "NCTS_TRANSLATION_RDM_SUPPORT=extend."
            ))
            real_moment_matrix || throw(ArgumentError(
                "NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true with RDM requires " *
                "NCTS_TRANSLATION_REAL=true."
            ))
        end
        return NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
            case.pop,
            case.ops,
            order;
            extra=something(qmbcertify_base_extra, 0),
            three_type=qmbcertify_base_three_type,
            real_moment_matrix,
            contiguous_rdm_k,
            contiguous_rdm_decomposition,
            contiguous_rdm_support,
            linear_state_opt_width,
            psd_state_opt_width,
        )
    end

    kwargs = (
        sign_symmetry=sign_symmetry,
        reflection_symmetry=reflection_symmetry,
        axis_rotation_symmetry=axis_rotation_symmetry,
        axis_rotation_equalities=axis_rotation_equalities,
        axis_rotation_quotient=axis_rotation_quotient,
        real_moment_matrix=real_moment_matrix,
        contiguous_rdm_k=contiguous_rdm_k,
        contiguous_rdm_decomposition=contiguous_rdm_decomposition,
        contiguous_rdm_support=contiguous_rdm_support,
        u1_symmetry=u1_symmetry,
        su2_symmetry=su2_symmetry,
        singlet_channel_equalities=singlet_channel_equalities,
        linear_state_opt_width=linear_state_opt_width,
        psd_state_opt_width=psd_state_opt_width,
    )
    if direct_linear
        return NCTSSoS._pauli_translation_base_linear_relaxation(
            case.pop,
            case.ops,
            order;
            kwargs...,
            base_su2_extend_rdm,
            su2_moment_quotient,
        )
    end
    return pauli_translation_invariant_moment_relaxation(
        case.pop,
        case.ops,
        order;
        kwargs...,
    )
end

function _print_target_only(
    n::Integer;
    order::Integer,
    sign_symmetry::Bool,
    reflection_symmetry::Bool,
    axis_rotation_symmetry::Bool,
    axis_rotation_equalities::Bool,
    real_moment_matrix::Bool,
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool,
    singlet_channel_equalities::Bool,
    moment_eq_h2::Bool,
    linear_state_opt_width,
    psd_state_opt_width,
    qmbcertify_base_construct::Bool,
    qmbcertify_base_extra,
    qmbcertify_base_three_type::Tuple{Int,Int},
)
    target = pauli_translation_structural_targets(
        n,
        order;
        sign_symmetry,
        reflection_symmetry,
        real_moment_matrix,
        axis_rotation_symmetry,
        axis_rotation_equalities,
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        contiguous_rdm_support,
        u1_symmetry,
        su2_symmetry,
        base_su2_extend_rdm,
        moment_eq_h2,
        linear_state_opt_width,
        psd_state_opt_width,
        qmbcertify_base_construct,
        qmbcertify_base_extra,
        qmbcertify_base_three_type,
    )

    println("\n## N = $n")
    println("\n#### Analytic structural target")
    println("\n| metric | value |")
    println("|:--|--:|")
    println("| construction performed | false |")
    println("| basis size | $(target.basis_size) |")
    println("| orbit basis size | $(target.orbit_basis_size) |")
    println("| axis-orbit basis size | $(target.axis_orbit_basis_size) |")
    println("| axis-orbit reduction ratio | $(@sprintf("%.6f", target.axis_reduction_ratio)) |")
    println("| axis-rotation PSD split | $(target.axis_rotation_symmetry) |")
    println("| momentum sectors | $(target.momentum_sector_count) |")
    println("| PSD blocks | $(target.n_blocks) |")
    println("| logical max block | $(target.logical_max_block) |")
    println("| solver-facing max block | $(target.psd_max_block) |")
    println("| logical total block side | $(target.logical_total_block_side) |")
    println("| solver-facing total block side | $(target.psd_total_block_side) |")
    println("| solver-facing symmetric entries | $(target.psd_symmetric_entries) |")
    println("| solver-facing symmetric bytes | $(_fmt_bytes(target.psd_symmetric_bytes)) |")
    println("| model-size gate status | `$(target.estimated_model_size_gate_status)` |")
    println("| model-size gate reason | $(target.estimated_model_size_gate_reason) |")
    println("| product-cache hits | $(target.product_cache_hits) |")
    println("| product-cache misses | $(target.product_cache_misses) |")
    println("| product-cache lookups | $(target.product_cache_lookups) |")
    println("| product-cache entries | $(target.product_cache_entries) |")
    println("| product-cache hit rate | $(@sprintf("%.6f", target.product_cache_hit_rate)) |")
    println("| Hamiltonian width | $(target.hamiltonian_width) |")
    println("| max closed contiguous RDM k | $(target.max_contiguous_rdm_k) |")
    println("| max linear state-opt width | $(target.max_linear_state_opt_width) |")
    println("| max PSD state-opt width | $(target.max_psd_state_opt_width) |")
    println("| contiguous RDM k | `$(target.contiguous_rdm_k)` |")
    println("| contiguous RDM decomposition | `$(target.contiguous_rdm_decomposition)` |")
    println("| contiguous RDM support | `$(target.contiguous_rdm_support)` |")
    println("| base SU(2) extend-RDM target | $(target.base_su2_extend_rdm) |")
    println("| contiguous RDM PSD blocks | $(length(target.rdm_psd_block_sizes)) |")
    println("| contiguous RDM zero rows | $(target.contiguous_rdm_zero_row_count) |")
    println("| SU(2) base full dense entries | $(target.su2_base_full_dense_entries) |")
    println("| SU(2) base active dense entries | $(target.su2_base_active_dense_entries) |")
    println("| SU(2) base reduced dense entries | $(target.su2_base_reduced_dense_entries) |")
    println("| SU(2) base off-block entries | $(target.su2_base_offblock_entry_count) |")
    println("| SU(2) base magnetic-copy entries | $(target.su2_base_copy_entry_count) |")
    println("| SU(2) base off-block zero-row budget | $(target.su2_base_offblock_zero_row_budget) |")
    println("| SU(2) base magnetic-copy zero-row budget | $(target.su2_base_magnetic_copy_zero_row_budget) |")
    println("| SU(2) base Wigner zero-row budget | $(target.su2_base_zero_row_budget) |")
    println("| SU(2) base singlet channels | $(target.su2_base_singlet_channel_count) |")
    println("| SU(2) base singlet channel support counts | `$(target.su2_base_singlet_channel_support_counts)` |")
    println("| SU(2) base singlet-channel equality candidate rows | $(target.su2_base_singlet_channel_equality_row_count) |")
    println("| axis-rotation equality rows | $(target.axis_rotation_equality_row_count) |")
    println("| axis-rotation raw moment keys | $(target.axis_rotation_raw_moment_key_count) |")
    println("| axis-rotation moment classes | $(target.axis_rotation_moment_class_count) |")
    println("| axis-rotation quotient moment keys | $(target.axis_rotation_quotient_moment_key_count) |")
    println("| axis-rotation forced-zero moment classes | $(target.axis_rotation_forced_zero_moment_class_count) |")
    println("| axis-rotation moment quotient ratio | $(@sprintf("%.6f", target.axis_rotation_moment_quotient_reduction_ratio)) |")
    println("| moment-equality rows | $(target.moment_equality_row_count) |")
    println("| known add-on zero rows | $(target.add_on_zero_row_count) |")
    println("| linear state-opt width | $(target.linear_state_opt_width) |")
    println("| linear state-opt rows | $(target.linear_state_opt_row_count) |")
    println("| linear state-opt candidate count | $(target.linear_state_opt_candidate_count) |")
    println("| PSD state-opt width | $(target.psd_state_opt_width) |")
    println("| PSD state-opt candidate count | $(target.psd_state_opt_candidate_count) |")
    println("| PSD state-opt PSD blocks | $(length(target.psd_state_opt_psd_block_sizes)) |")
    println("| solve supported | $(target.solve_supported) |")
    println("| solve blocker | `$(target.solve_blocker)` |")
    println("| solve blocker reason | $(target.solve_blocker_reason) |")
    println("| solve unsupported block features | `$(target.solve_unsupported_block_features)` |")
    println("| solve unsupported zero-row features | `$(target.solve_unsupported_zero_features)` |")

    _print_histogram("Logical block histogram", target.logical_block_histogram)
    _print_histogram("Solver-facing block histogram", target.psd_block_histogram)
    _print_histogram(
        "Logical feature block histogram",
        target.logical_block_feature_histogram;
        key_label="feature / decomposition / size",
    )
    _print_histogram(
        "Solver-facing feature block histogram",
        target.psd_block_feature_histogram;
        key_label="feature / decomposition / size",
    )
    isempty(target.known_zero_constraint_feature_histogram) || _print_histogram(
        "Known zero-row feature histogram",
        target.known_zero_constraint_feature_histogram;
        key_label="feature / decomposition / reason",
    )
    _print_coefficient_domain_histograms(
        "Structural target",
        target.block_coefficient_domains,
        target.block_exact_coefficient_domains,
    )
    _print_histogram("Axis orbit size histogram", target.axis_orbit_size_histogram; key_label="axis orbit size")
    isempty(target.rdm_logical_block_sizes) || _print_histogram(
        "Contiguous RDM logical block histogram",
        _histogram_pairs(target.rdm_logical_block_sizes),
    )
    isempty(target.rdm_psd_block_sizes) || _print_histogram(
        "Contiguous RDM solver-facing block histogram",
        _histogram_pairs(target.rdm_psd_block_sizes),
    )
    isempty(target.psd_state_opt_logical_block_sizes) || _print_histogram(
        "PSD state-opt logical block histogram",
        _histogram_pairs(target.psd_state_opt_logical_block_sizes),
    )
    isempty(target.psd_state_opt_psd_block_sizes) || _print_histogram(
        "PSD state-opt solver-facing block histogram",
        _histogram_pairs(target.psd_state_opt_psd_block_sizes),
    )
    _print_qmbcertify_rdm_reference_aggregate(target)
    _print_qmbcertify_rdm_references(target.qmbcertify_rdm_reference_blocks)
    _print_qmbcertify_base_reference(target)
    _print_su2_accounting_records(target.su2_base_accounting_records)
    _print_su2_structural_target(
        n,
        order;
        real_moment_matrix,
        reflection_symmetry,
    )
end

function _run_one(
    n::Integer;
    order::Integer,
    qmbcertify_base_construct::Bool,
    direct_linear::Bool,
    sign_symmetry::Bool,
    reflection_symmetry::Bool,
    axis_rotation_symmetry::Bool,
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    real_moment_matrix::Bool,
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool,
    singlet_channel_equalities::Bool,
    moment_eq_h2::Bool,
    linear_state_opt_width,
    psd_state_opt_width,
    qmbcertify_base_extra,
    qmbcertify_base_three_type::Tuple{Int,Int},
    run_index::Integer=1,
    repeats::Integer=1,
)
    heading = repeats == 1 ? "N = $n" : "N = $n (repeat $run_index/$repeats)"
    println("\n## $heading")
    println("\n| phase | wall time (s) | allocated | GC time (s) | process memory |")
    println("|:--|--:|--:|--:|:--|")

    case, _ = _timed("setup Heisenberg case", () -> _heisenberg_case(n; moment_eq_h2))
    relaxation_label = qmbcertify_base_construct ?
        "_pauli_qmbcertify_chain_base_linear_relaxation" :
        direct_linear ?
            "_pauli_translation_base_linear_relaxation" :
            "pauli_translation_invariant_moment_relaxation"
    data_report, relax_stats = _timed(
        relaxation_label,
        () -> _build_translation_profile(
            case,
            order;
            qmbcertify_base_construct,
            direct_linear,
            sign_symmetry,
            reflection_symmetry,
            axis_rotation_symmetry,
            axis_rotation_equalities,
            axis_rotation_quotient,
            real_moment_matrix,
            contiguous_rdm_k,
            contiguous_rdm_decomposition,
            contiguous_rdm_support,
            u1_symmetry,
            su2_symmetry,
            base_su2_extend_rdm,
            su2_moment_quotient,
            singlet_channel_equalities,
            moment_eq_h2,
            linear_state_opt_width,
            psd_state_opt_width,
            qmbcertify_base_extra,
            qmbcertify_base_three_type,
        ),
    )
    data, report = data_report
    linear = (direct_linear || qmbcertify_base_construct) ? data : data.linear
    metrics = translation_report_metrics(report)
    sos_dual_model_size_gate =
        _translation_profile_model_size_gate(metrics.estimated_sos_dual_dense_schur_bytes)

    println("\n#### Summary")
    println("\n| metric | value |")
    println("|:--|--:|")
    println("| basis size | $(metrics.basis_size) |")
    println("| orbit basis size | $(metrics.orbit_basis_size) |")
    println("| axis-orbit closed | $(metrics.axis_orbit_closed) |")
    println("| axis-orbit basis size | $(metrics.axis_orbit_basis_size) |")
    println("| axis-orbit reduction ratio | $(@sprintf("%.6f", metrics.axis_reduction_ratio)) |")
    println("| sign symmetry | $(metrics.sign_symmetry) |")
    println("| reflection symmetry | $(metrics.reflection_symmetry) |")
    println("| axis-rotation symmetry | $(metrics.axis_rotation_symmetry) |")
    println("| axis-rotation equalities | $(metrics.axis_rotation_equalities) |")
    println("| axis-rotation quotient | $(metrics.axis_rotation_quotient) |")
    println("| full discrete profile | $(metrics.full_discrete_profile) |")
    println("| missing discrete symmetries | `$(metrics.missing_discrete_symmetries)` |")
    println("| SU(2) moment symmetry | $(metrics.su2_moment_symmetry) |")
    println("| SU(2) reflection moment symmetry | $(metrics.su2_reflection_moment_symmetry) |")
    println("| contiguous RDM | $(metrics.contiguous_rdm) |")
    println("| U(1) RDM symmetry | $(metrics.u1_rdm_symmetry) |")
    println("| SU(2) RDM symmetry | $(metrics.su2_rdm_symmetry) |")
    println("| linear state-opt | $(metrics.linear_state_opt) |")
    println("| PSD state-opt | $(metrics.psd_state_opt) |")
    println("| momentum sectors | $(metrics.momentum_sector_count) |")
    println("| PSD blocks | $(metrics.n_blocks) |")
    println("| logical max block | $(metrics.logical_max_block) |")
    println("| solver-facing max block | $(metrics.psd_max_block) |")
    println("| moment count | $(metrics.moment_count) |")
    println("| linear moments | $(metrics.linear_moment_count) |")
    println("| SU(2) moment quotient | $(metrics.su2_moment_quotient) |")
    println("| SU(2) raw moments | $(metrics.su2_moment_raw_count) |")
    println("| SU(2) quotient moments | $(metrics.su2_moment_quotient_count) |")
    println("| SU(2) quotient ratio | $(@sprintf("%.6f", metrics.su2_moment_quotient_reduction_ratio)) |")
    println("| SU(2) quotient support orbits | $(metrics.su2_moment_support_orbit_count) |")
    println("| SU(2) quotient support/channel counts | `$(metrics.su2_moment_singlet_channel_support_counts)` |")
    println("| SU(2) quotient pivot residual | $(metrics.su2_moment_max_pivot_residual) |")
    println("| SU(2) quotient invariant residual | $(metrics.su2_moment_max_invariant_residual) |")
    println("| SU(2) quotient reconstruction residual | $(metrics.su2_moment_max_reconstruction_residual) |")
    println("| SU(2) quotient max condition | $(metrics.su2_moment_max_condition) |")
    println("| SU(2) quotient eliminated zero rows | $(metrics.su2_moment_eliminated_zero_row_count) |")
    println("| SOS-dual equality upper bound | $(metrics.estimated_sos_dual_scalar_equalities_upper_bound) |")
    println("| SOS-dual dense-Schur bytes | $(_fmt_bytes(metrics.estimated_sos_dual_dense_schur_bytes)) |")
    println("| SOS-dual model-size gate status | `$(sos_dual_model_size_gate.status)` |")
    println("| SOS-dual model-size gate reason | $(isempty(sos_dual_model_size_gate.reason) ? "N/A" : sos_dual_model_size_gate.reason) |")
    println("| SOS-dual model-size gate estimated RSS bytes | $(sos_dual_model_size_gate.estimated_rss_bytes) |")
    println("| SOS-dual model-size gate MemAvailable bytes | $(isnothing(sos_dual_model_size_gate.mem_available_bytes) ? "unknown" : sos_dual_model_size_gate.mem_available_bytes) |")
    println("| SOS-dual model-size gate max RSS fraction | $(sos_dual_model_size_gate.max_rss_fraction) |")
    println("| zero constraints | $(metrics.zero_constraint_count) |")
    println("| contiguous RDM blocks | $(metrics.contiguous_rdm_block_count) |")
    println("| contiguous RDM zero rows | $(metrics.contiguous_rdm_zero_row_count) |")
    println("| scalar equality rows | $(metrics.scalar_equality_row_count) |")
    println("| axis-rotation equality rows | $(metrics.axis_rotation_equality_row_count) |")
    println("| axis-rotation raw moment keys | $(metrics.axis_rotation_raw_moment_key_count) |")
    println("| axis-rotation moment classes | $(metrics.axis_rotation_moment_class_count) |")
    println("| axis-rotation quotient moment keys | $(metrics.axis_rotation_quotient_moment_key_count) |")
    println("| axis-rotation forced-zero moment classes | $(metrics.axis_rotation_forced_zero_moment_class_count) |")
    println("| axis-rotation moment quotient ratio | $(@sprintf("%.6f", metrics.axis_rotation_moment_quotient_reduction_ratio)) |")
    println("| moment-equality rows | $(metrics.moment_equality_row_count) |")
    println("| linear state-opt rows | $(metrics.linear_state_opt_row_count) |")
    println("| SU(2) singlet-channel equality rows | $(metrics.su2_singlet_channel_equality_row_count) |")
    println("| SU(2) base zero rows | $(metrics.su2_base_zero_row_count) |")
    println("| SU(2) spin off-block rows | $(metrics.su2_base_spin_offblock_row_count) |")
    println("| SU(2) magnetic off-diagonal rows | $(metrics.su2_base_magnetic_offdiag_row_count) |")
    println("| SU(2) magnetic-copy rows | $(metrics.su2_base_magnetic_copy_row_count) |")
    println("| PSD state-opt blocks | $(metrics.psd_state_opt_block_count) |")
    println("| solve supported | $(metrics.solve_supported) |")
    println("| solve blocker | `$(isnothing(metrics.solve_blocker) ? "nothing" : metrics.solve_blocker)` |")
    if !metrics.solve_supported
        println("| solve blocker reason | $(metrics.solve_blocker_reason) |")
        println("| solve unsupported block features | `$(metrics.solve_unsupported_block_features)` |")
        println("| solve unsupported zero-row features | `$(metrics.solve_unsupported_zero_features)` |")
    end
    println("| product-cache hits | $(metrics.product_cache_hits) |")
    println("| product-cache misses | $(metrics.product_cache_misses) |")
    println("| product-cache lookups | $(metrics.product_cache_lookups) |")
    println("| product-cache entries | $(metrics.product_cache_entries) |")
    println("| product-cache hit rate | $(@sprintf("%.6f", metrics.product_cache_hit_rate)) |")
    println("| report construction time (s) | $(@sprintf("%.9f", metrics.construction_time_seconds)) |")
    println("| outer relaxation wall time (s) | $(@sprintf("%.9f", relax_stats.time)) |")
    println("| QMBCertify base construct | $qmbcertify_base_construct |")
    println("| direct linear | $direct_linear |")
    println("| symbolic constraints | $((direct_linear || qmbcertify_base_construct) ? "n/a" : string(length(data.constraints))) |")
    println("| linear PSD blocks | $(length(linear.psd_blocks_lin)) |")
    println("| linear zero constraints | $(length(linear.zero_constraints)) |")

    _print_su2_structural_target(
        n,
        order;
        real_moment_matrix,
        reflection_symmetry,
    )
    _print_stage_times(metrics)
    _print_histogram("Logical block histogram", metrics.logical_block_histogram)
    _print_histogram("Solver-facing block histogram", metrics.psd_block_histogram)
    _print_histogram(
        "Logical feature block histogram",
        metrics.logical_block_feature_histogram;
        key_label="feature / decomposition / size",
    )
    _print_histogram(
        "Solver-facing feature block histogram",
        metrics.psd_block_feature_histogram;
        key_label="feature / decomposition / size",
    )
    isempty(metrics.zero_constraint_feature_histogram) || _print_histogram(
        "Zero constraint feature histogram",
        metrics.zero_constraint_feature_histogram;
        key_label="feature / decomposition / reason",
    )
    _print_coefficient_domain_histograms(
        "Block",
        metrics.block_coefficient_domains,
        metrics.block_exact_coefficient_domains,
    )
    isempty(metrics.contiguous_rdm_logical_block_sizes) || _print_histogram(
        "Constructed contiguous RDM logical block histogram",
        _histogram_pairs(metrics.contiguous_rdm_logical_block_sizes),
    )
    isempty(metrics.contiguous_rdm_psd_block_sizes) || _print_histogram(
        "Constructed contiguous RDM solver-facing block histogram",
        _histogram_pairs(metrics.contiguous_rdm_psd_block_sizes),
    )
    isempty(metrics.psd_state_opt_logical_block_sizes) || _print_histogram(
        "Constructed PSD state-opt logical block histogram",
        _histogram_pairs(metrics.psd_state_opt_logical_block_sizes),
    )
    isempty(metrics.psd_state_opt_psd_block_sizes) || _print_histogram(
        "Constructed PSD state-opt solver-facing block histogram",
        _histogram_pairs(metrics.psd_state_opt_psd_block_sizes),
    )
    return (
        n=Int(n),
        run_index=Int(run_index),
        construction_time_seconds=metrics.construction_time_seconds,
        outer_relaxation_wall_time_seconds=relax_stats.time,
        allocated_bytes=Int(relax_stats.bytes),
        gc_time_seconds=relax_stats.gctime,
        linear_moments=metrics.linear_moment_count,
        zero_constraints=metrics.zero_constraint_count,
        psd_blocks=metrics.n_blocks,
        solver_facing_max_block=metrics.psd_max_block,
        contiguous_rdm_blocks=metrics.contiguous_rdm_block_count,
        psd_state_opt_blocks=metrics.psd_state_opt_block_count,
        construction_stage_time_seconds=Dict{Symbol,Float64}(
            (stage isa Symbol ? stage : Symbol(stage)) => Float64(seconds)
            for (stage, seconds) in metrics.construction_stage_time_seconds
        ),
    )
end

function main()
    ns = _parse_int_list_env("NCTS_TRANSLATION_NS", "8,12")

    order = parse(Int, get(ENV, "NCTS_TRANSLATION_ORDER", "2"))
    target_only = _parse_bool_env("NCTS_TRANSLATION_TARGET_ONLY", false)
    warmup = _parse_bool_env("NCTS_TRANSLATION_WARMUP", true)
    repeats = _parse_positive_int_env("NCTS_TRANSLATION_REPEATS", 1)
    repeat_stage_limit = _parse_positive_int_env("NCTS_TRANSLATION_REPEAT_STAGE_LIMIT", 6)
    qmbcertify_base_construct =
        _parse_bool_env("NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT", false)
    direct_linear = _parse_bool_env("NCTS_TRANSLATION_DIRECT_LINEAR", false)
    sign_symmetry = _parse_bool_env("NCTS_TRANSLATION_SIGN", true)
    reflection_symmetry = _parse_bool_env("NCTS_TRANSLATION_REFLECTION", false)
    axis_rotation_symmetry = _parse_bool_env("NCTS_TRANSLATION_AXIS_ROTATION", false)
    axis_rotation_equalities = _parse_bool_env("NCTS_TRANSLATION_AXIS_EQUALITIES", false)
    axis_rotation_quotient = _parse_bool_env("NCTS_TRANSLATION_AXIS_QUOTIENT", false)
    real_moment_matrix = _parse_bool_env("NCTS_TRANSLATION_REAL", true)
    contiguous_rdm_k = _parse_optional_int_env("NCTS_TRANSLATION_RDM_K")
    default_contiguous_rdm_decomposition =
        qmbcertify_base_construct ? :qmbcertify : :full
    default_contiguous_rdm_support = qmbcertify_base_construct ? :extend : :closed
    contiguous_rdm_decomposition = Symbol(get(
        ENV,
        "NCTS_TRANSLATION_RDM_DECOMPOSITION",
        string(default_contiguous_rdm_decomposition),
    ))
    contiguous_rdm_support = Symbol(get(
        ENV,
        "NCTS_TRANSLATION_RDM_SUPPORT",
        string(default_contiguous_rdm_support),
    ))
    u1_symmetry = _parse_bool_env("NCTS_TRANSLATION_U1", false)
    su2_symmetry = _parse_bool_env("NCTS_TRANSLATION_SU2", false)
    base_su2_extend_rdm = _parse_bool_env("NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM", false)
    su2_moment_quotient = _parse_bool_env(
        "NCTS_TRANSLATION_SU2_MOMENT_QUOTIENT",
        base_su2_extend_rdm,
    )
    singlet_channel_equalities = _parse_bool_env("NCTS_TRANSLATION_SINGLET_EQUALITIES", false)
    moment_eq_h2 = _parse_bool_env("NCTS_TRANSLATION_MOMENT_EQ_H2", false)
    linear_state_opt_width = _parse_optional_int_env("NCTS_TRANSLATION_LSO_WIDTH")
    psd_state_opt_width = _parse_optional_int_env("NCTS_TRANSLATION_PSO_WIDTH")
    qmbcertify_base_extra =
        _parse_optional_int_env("NCTS_TRANSLATION_QMBCERTIFY_BASE_EXTRA")
    qmbcertify_base_three_type =
        _parse_int_tuple2_env("NCTS_TRANSLATION_QMBCERTIFY_BASE_THREE_TYPE", "1,1")

    _check_size_guard(ns)
    _check_base_su2_extend_rdm_profile_guard(
        target_only=target_only,
        direct_linear=direct_linear,
        base_su2_extend_rdm=base_su2_extend_rdm,
        su2_moment_quotient=su2_moment_quotient,
        su2_symmetry=su2_symmetry,
        contiguous_rdm_decomposition=contiguous_rdm_decomposition,
        contiguous_rdm_support=contiguous_rdm_support,
        real_moment_matrix=real_moment_matrix,
        reflection_symmetry=reflection_symmetry,
    )
    _check_singlet_channel_profile_guard(
        singlet_channel_equalities=singlet_channel_equalities,
        su2_symmetry=su2_symmetry,
    )
    println("# Pauli translation-invariant no-solver profile")
    println("\n- generated: `$(Dates.now())`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- solver calls: none")
    println("- order: `$order`")
    println("- target_only: `$target_only`")
    println("- warmup: `$warmup`")
    println("- repeats: `$repeats`")
    println("- repeat_stage_limit: `$repeat_stage_limit`")
    println("- qmbcertify_base_construct: `$qmbcertify_base_construct`")
    println("- direct_linear: `$direct_linear`")
    println("- sign_symmetry: `$sign_symmetry`")
    println("- reflection_symmetry: `$reflection_symmetry`")
    println("- axis_rotation_symmetry: `$axis_rotation_symmetry`")
    println("- axis_rotation_equalities: `$axis_rotation_equalities`")
    println("- axis_rotation_quotient: `$axis_rotation_quotient`")
    println("- real_moment_matrix: `$real_moment_matrix`")
    println("- contiguous_rdm_k: `$(isnothing(contiguous_rdm_k) ? "nothing" : contiguous_rdm_k)`")
    println("- contiguous_rdm_decomposition: `$contiguous_rdm_decomposition`")
    println("- contiguous_rdm_support: `$contiguous_rdm_support`")
    println("- u1_symmetry: `$u1_symmetry`")
    println("- su2_symmetry: `$su2_symmetry`")
    println("- base_su2_extend_rdm: `$base_su2_extend_rdm`")
    println("- su2_moment_quotient: `$su2_moment_quotient`")
    println("- singlet_channel_equalities: `$singlet_channel_equalities`")
    println("- moment_eq_h2: `$moment_eq_h2`")
    println("- linear_state_opt_width: `$(isnothing(linear_state_opt_width) ? "nothing" : linear_state_opt_width)`")
    println("- psd_state_opt_width: `$(isnothing(psd_state_opt_width) ? "nothing" : psd_state_opt_width)`")
    println("- qmbcertify_base_extra: `$(isnothing(qmbcertify_base_extra) ? "nothing" : qmbcertify_base_extra)`")
    println("- qmbcertify_base_three_type: `$qmbcertify_base_three_type`")

    if target_only
        for n in ns
            _print_target_only(
                n;
                order,
                sign_symmetry,
                reflection_symmetry,
                axis_rotation_symmetry,
                axis_rotation_equalities,
                real_moment_matrix,
                contiguous_rdm_k,
                contiguous_rdm_decomposition,
                contiguous_rdm_support,
                u1_symmetry,
                su2_symmetry,
                base_su2_extend_rdm,
                su2_moment_quotient,
                singlet_channel_equalities,
                moment_eq_h2,
                linear_state_opt_width,
                psd_state_opt_width,
                qmbcertify_base_construct,
                qmbcertify_base_extra,
                qmbcertify_base_three_type,
            )
        end
        return nothing
    end

    if warmup
        warm_n = minimum(ns)
        println("\n## Warmup N = $warm_n")
        println("\n| phase | wall time (s) | allocated | GC time (s) | process memory |")
        println("|:--|--:|--:|--:|:--|")
        warm_case, _ = _timed("setup Heisenberg case", () -> _heisenberg_case(warm_n; moment_eq_h2))
        warm_relaxation_label = qmbcertify_base_construct ?
            "_pauli_qmbcertify_chain_base_linear_relaxation" :
            direct_linear ?
                "_pauli_translation_base_linear_relaxation" :
                "pauli_translation_invariant_moment_relaxation"
        _timed(
            warm_relaxation_label,
            () -> _build_translation_profile(
                warm_case,
                order;
                qmbcertify_base_construct,
                direct_linear,
                sign_symmetry,
                reflection_symmetry,
                axis_rotation_symmetry,
                axis_rotation_equalities,
                axis_rotation_quotient,
                real_moment_matrix,
                contiguous_rdm_k,
                contiguous_rdm_decomposition,
                contiguous_rdm_support,
                u1_symmetry,
                su2_symmetry,
                base_su2_extend_rdm,
                su2_moment_quotient,
                singlet_channel_equalities,
                moment_eq_h2,
                linear_state_opt_width,
                psd_state_opt_width,
                qmbcertify_base_extra,
                qmbcertify_base_three_type,
            ),
        )
    end

    repeat_results = NamedTuple[]
    for n in ns, run_index in 1:repeats
        push!(repeat_results, _run_one(
            n;
            order,
            qmbcertify_base_construct,
            direct_linear,
            sign_symmetry,
            reflection_symmetry,
            axis_rotation_symmetry,
            axis_rotation_equalities,
            axis_rotation_quotient,
            real_moment_matrix,
            contiguous_rdm_k,
            contiguous_rdm_decomposition,
            contiguous_rdm_support,
            u1_symmetry,
            su2_symmetry,
            base_su2_extend_rdm,
            su2_moment_quotient,
            singlet_channel_equalities,
            moment_eq_h2,
            linear_state_opt_width,
            psd_state_opt_width,
            qmbcertify_base_extra,
            qmbcertify_base_three_type,
            run_index=run_index,
            repeats=repeats,
        ))
    end
    _print_repeat_summary(repeat_results; stage_limit=repeat_stage_limit)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
