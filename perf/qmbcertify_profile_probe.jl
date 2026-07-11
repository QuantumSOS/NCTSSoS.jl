#!/usr/bin/env julia

# Probe QMBCertify profile variants without emitting reviewed fixture rows.
#
# Example:
#   NCTS_QMB_PROFILE=A1 \
#   NCTS_QMB_PROBE_LABEL=A1_rdm_only_probe \
#   NCTS_QMB_ARG_OVERRIDES='lso = false' \
#   julia --project=. --startup-file=no perf/qmbcertify_profile_probe.jl
#
# Safety:
#   Large L probes use the same NCTS_QMB_ALLOW_LARGE pressure gate as reviewed
#   reference runs.  Set explicit wall/RSS estimates and require safe
#   load/memory telemetry before launching a real profile solve.

Base.include(@__MODULE__, joinpath(@__DIR__, "qmbcertify_reference_runs.jl"))

function _probe_overrides()
    raw = strip(get(ENV, "NCTS_QMB_ARG_OVERRIDES", ""))
    path = strip(get(ENV, "NCTS_QMB_ARG_OVERRIDES_FILE", ""))
    !isempty(raw) && !isempty(path) && throw(ArgumentError(
        "Set only one of NCTS_QMB_ARG_OVERRIDES or NCTS_QMB_ARG_OVERRIDES_FILE."
    ))
    !isempty(path) && return TOML.parsefile(path)
    !isempty(raw) && return TOML.parse(raw)
    return Dict{String,Any}()
end

function _probe_args(profile::AbstractString)
    args = copy(_profile_args(profile))
    overrides = _probe_overrides()
    for (key, value) in overrides
        args[String(key)] = value
    end
    return args, overrides
end

function _probe_label(profile::AbstractString)
    raw = strip(get(ENV, "NCTS_QMB_PROBE_LABEL", ""))
    return isempty(raw) ? "$(profile)_probe" : raw
end

function _probe_case_id(profile::AbstractString, len::Integer)
    raw = strip(get(ENV, "NCTS_QMB_CASE_ID", ""))
    return isempty(raw) ? "$(profile)_L$(len)" : raw
end

_profile_probe_lengths() = _requested_lengths()

_probe_accepted_statuses(profile::AbstractString) =
    _reviewed_solve_status_policy(profile)

function _check_large_profile_probe_variant(
    profile::AbstractString,
    lengths::AbstractVector{Int},
    overrides,
)
    !isempty(overrides) && return nothing
    fixture = _profile_fixture(profile)
    for len in lengths
        case = _pending_case(fixture, profile, len)
        case === nothing && continue
        get(case, "run_queue_status", "queued") == "deleted_until_evidence_gate" ||
            continue
        throw(ArgumentError(
            "Refusing deleted QMBCertify profile probe $(profile)_L$(len) " *
            "without a declared variant.  Set NCTS_QMB_ARG_OVERRIDES or " *
            "NCTS_QMB_ARG_OVERRIDES_FILE before probing a deleted large row."
        ))
    end
    return nothing
end

function _toml_string_array(values)
    return "[" * join((_toml_string(value) for value in values), ", ") * "]"
end

function _probe_arg_summary(args)
    isempty(args) && return "none"
    parts = ["$(key)=$(repr(value))" for (key, value) in args]
    sort!(parts)
    return join(parts, ", ")
end

function _print_profile_probe_rows(rows, qmb_path::AbstractString, env_id::AbstractString)
    isempty(rows) && return nothing
    println("\n## Profile Probe Rows")
    for row in rows
        println("\n[[profile_probes]]")
        println("id = \"$(row.profile)_L$(row.length)_$(_fixture_date_id())\"")
        println("case_id = \"$(row.case_id)\"")
        println("profile = \"$(row.profile)\"")
        println("parent_profile = \"$(row.parent_profile)\"")
        println("accepted_termination_statuses = $(_toml_string_array(row.accepted_statuses))")
        println("length = $(row.length)")
        println("status = \"probe_run\"")
        println("environment = \"$env_id\"")
        println("command_args_summary = $(_toml_string(_probe_arg_summary(row.command_args)))")
        println("overrides_summary = $(_toml_string(_probe_arg_summary(row.overrides)))")
        println("objective = $(row.objective)")
        println("objective_per_site = $(row.objective_per_site)")
        println("objective_total_estimate = $(row.objective_total_estimate)")
        _print_toml_field("moment_count", row.moment_count)
        println("psd_block_histogram = $(repr([[first(pair), last(pair)] for pair in row.psd_block_histogram]))")
        println("psd_max_block = $(row.psd_max_block)")
        println("psd_scalar_variables = $(row.psd_scalar_variables)")
        _print_toml_field("build_time_seconds", row.build_time_seconds)
        _print_toml_field("solve_time_seconds", row.solve_time_seconds)
        _print_toml_field("termination_status", _toml_string(row.termination_status))
        _print_toml_field("solution_status", _toml_string(row.solution_status))
        _print_toml_field("estimated_wall_seconds", _qmb_estimated_wall_seconds())
        _print_toml_field("estimated_rss_gib", _qmb_estimated_rss_gib())
        println("total_wall_seconds = $(row.total_wall_seconds)")
        _print_toml_field("peak_rss_bytes", row.peak_rss_bytes)
        println("notes = $(_toml_string("Probe row only; confirm the profile before promoting it to reviewed cases."))")
    end
    _print_environment_fixture_block(qmb_path; table="profile_probe_environments.$env_id")
end

function _print_profile_failed_attempt_rows(rows, qmb_path::AbstractString, env_id::AbstractString)
    isempty(rows) && return nothing
    println("\n## Failed Attempt Rows")
    for row in rows
        timings = row.timings
        println("\n[[failed_attempts]]")
        println("id = \"$(row.profile)_L$(row.length)_$(_fixture_date_id())\"")
        println("case_id = \"$(row.case_id)\"")
        println("profile = \"$(row.profile)\"")
        println("parent_profile = \"$(row.parent_profile)\"")
        println("accepted_termination_statuses = $(_toml_string_array(row.accepted_statuses))")
        println("length = $(row.length)")
        println("attempted_at = \"$(Dates.today())\"")
        println("status = \"rejected_nonoptimal_solve\"")
        println("reason = \"Guarded harness rejects non-optimal reference rows.\"")
        println("environment = \"$env_id\"")
        println("command_args_summary = $(_toml_string(_probe_arg_summary(row.command_args)))")
        println("overrides_summary = $(_toml_string(_probe_arg_summary(row.overrides)))")
        _print_toml_field("termination_status", _toml_optional_string(timings.termination_status))
        _print_toml_field("solution_status", _toml_optional_string(timings.solution_status))
        _print_toml_field("objective_reported", row.objective_reported)
        _print_toml_field("estimated_wall_seconds", _qmb_estimated_wall_seconds())
        _print_toml_field("estimated_rss_gib", _qmb_estimated_rss_gib())
        _print_toml_field("reported_sdp_n", timings.max_block)
        _print_toml_field("reported_sdp_m", timings.moment_count)
        _print_toml_field("basis_time_seconds", timings.basis_time)
        _print_toml_field("block_diagonalization_time_seconds", timings.block_time)
        _print_toml_field("linear_state_opt_time_seconds", timings.lso_time)
        _print_toml_field("psd_state_opt_time_seconds", timings.pso_time)
        _print_toml_field("rdm_time_seconds", timings.rdm_time)
        _print_toml_field("solve_time_seconds", timings.solve_time)
        println("notes = $(_toml_string("Profile probe failed; confirm the parent profile before promoting any related row to reviewed cases."))")
    end
    _print_environment_fixture_block(qmb_path; table="failed_attempt_environments.$env_id")
end

function probe_main()
    profile = String(get(ENV, "NCTS_QMB_PROFILE", "A1"))
    label = _probe_label(profile)
    args, overrides = _probe_args(profile)
    accepted_statuses = _probe_accepted_statuses(profile)
    source = _profile_source(profile)
    qmb_path = String(get(ENV, "NCTS_QMBCERTIFY_PATH", _default_qmbcertify_path(source)))
    env_id = _environment_id()
    bootstrap_only = _parse_bool_env("NCTS_QMB_BOOTSTRAP_ONLY", false)
    requested_lengths = _profile_probe_lengths()
    _check_size_guard(requested_lengths; pressure_guard=!bootstrap_only)
    bootstrap_only ||
        _check_large_profile_probe_variant(profile, requested_lengths, overrides)

    println("# QMBCertify profile probe")
    println("\n- generated: `$(Dates.now())`")
    println("- base profile: `$profile`")
    println("- probe label: `$label`")
    println("- L list: `$(join(requested_lengths, ","))`")
    println("- overrides: `$(isempty(overrides) ? "none" : overrides)`")
    println("- accepted termination statuses: `$(join(accepted_statuses, ", "))`")
    println("- QMBCertify path: `$qmb_path`")
    println("- environment id: `$env_id`")
    println("- bootstrap_only: `$bootstrap_only`")

    if bootstrap_only
        println("\nBootstrap complete; no profile probe solve requested.")
        _print_environment_fixture_block(qmb_path; table="profile_probe_environments.$env_id")
        return nothing
    end

    qmb = _load_qmbcertify(qmb_path, source)
    rows = NamedTuple[]
    failed_attempts = NamedTuple[]
    failed = false
    for len in requested_lengths
        println("\n## L = $len")
        try
            push!(rows, merge(
                _run_one(len, label, args, qmb; status_policy_profile=profile),
                (;
                    parent_profile=profile,
                    case_id=_probe_case_id(profile, len),
                    command_args=copy(args),
                    overrides=copy(overrides),
                    accepted_statuses=accepted_statuses,
                ),
            ))
        catch err
            failed = true
            if err isa NonOptimalReferenceSolve
                push!(failed_attempts, (;
                    profile=label,
                    parent_profile=profile,
                    case_id=_probe_case_id(profile, len),
                    length=len,
                    objective_reported=err.objective,
                    timings=err.timings,
                    command_args=copy(args),
                    overrides=copy(overrides),
                    accepted_statuses=accepted_statuses,
                ))
            end
            println("\n## L = $len failed")
            println("\n- error type: `$(typeof(err))`")
            println("- error: `$(sprint(showerror, err))`")
            Base.show_backtrace(stdout, catch_backtrace())
            println()
            break
        end
        flush(stdout)
    end
    _print_profile_probe_rows(rows, qmb_path, env_id)
    _print_profile_failed_attempt_rows(failed_attempts, qmb_path, env_id)
    failed && exit(1)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    probe_main()
end
