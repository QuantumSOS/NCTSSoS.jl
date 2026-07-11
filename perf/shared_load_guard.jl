function _ncts_load_guard_parse_bool(name::AbstractString, default::Bool)
    raw = lowercase(strip(get(ENV, name, string(default))))
    raw in ("true", "1", "yes", "y") && return true
    raw in ("false", "0", "no", "n") && return false
    throw(ArgumentError("$name must be a boolean value, got $(repr(raw))."))
end

function _ncts_load_guard_parse_float(name::AbstractString, default::Real)
    return parse(Float64, get(ENV, name, string(default)))
end

function _ncts_load_guard_parse_int(name::AbstractString, default::Integer)
    return parse(Int, get(ENV, name, string(default)))
end

function _ncts_load_guard_optional_float(name::AbstractString)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    lowercase(raw) in ("nothing", "none", "false", "0") && return nothing
    return parse(Float64, raw)
end

function _ncts_cpu_max_quota_count(raw::AbstractString)
    parts = split(strip(raw))
    length(parts) >= 2 || return nothing
    parts[1] == "max" && return nothing
    quota = parse(Float64, parts[1])
    period = parse(Float64, parts[2])
    quota > 0 && period > 0 || return nothing
    return quota / period
end

function _ncts_cgroup_v2_cpu_count(path::AbstractString="/sys/fs/cgroup/cpu.max")
    isfile(path) || return nothing
    return _ncts_cpu_max_quota_count(read(path, String))
end

function _ncts_nproc_cpu_count()
    try
        raw = strip(readchomp(`nproc`))
        isempty(raw) && return nothing
        return parse(Float64, raw)
    catch
        return nothing
    end
end

function _ncts_effective_cpu_count()
    cgroup_count = _ncts_cgroup_v2_cpu_count()
    !isnothing(cgroup_count) && return cgroup_count
    nproc_count = _ncts_nproc_cpu_count()
    !isnothing(nproc_count) && return nproc_count
    return Float64(Sys.CPU_THREADS)
end

function _ncts_mem_available_bytes(path::AbstractString="/proc/meminfo")
    isfile(path) || return nothing
    for line in eachline(path)
        startswith(line, "MemAvailable:") || continue
        parts = split(line)
        length(parts) >= 2 || return nothing
        return parse(Int, parts[2]) * 1024
    end
    return nothing
end

function _ncts_loadavg_triplet(path::AbstractString="/proc/loadavg")
    isfile(path) || return nothing
    raw = strip(read(path, String))
    isempty(raw) && return nothing
    parts = split(raw)
    length(parts) >= 3 || return nothing
    return (
        parse(Float64, parts[1]),
        parse(Float64, parts[2]),
        parse(Float64, parts[3]),
    )
end

function _ncts_large_run_pressure_status(;
    load1,
    load5,
    load15,
    cpu_count::Real,
    max_load_per_cpu::Real,
)
    all(isnothing, (load1, load5, load15)) && return "blocked_missing_loadavg"
    threshold = Float64(cpu_count) * Float64(max_load_per_cpu)
    for load in (load1, load5, load15)
        isnothing(load) && continue
        load > threshold && return "blocked_overloaded_remote"
    end
    return "ok"
end

function _ncts_memory_pressure_status(;
    mem_available_bytes,
    estimated_rss_bytes,
    max_rss_fraction::Real=1.0,
)
    isnothing(estimated_rss_bytes) && return "ok"
    isnothing(mem_available_bytes) && return "blocked_missing_memory_telemetry"
    estimated_rss_bytes > mem_available_bytes && return "blocked_insufficient_memory"
    estimated_rss_bytes <= Float64(max_rss_fraction) * mem_available_bytes ||
        return "blocked_insufficient_memory_headroom"
    return "ok"
end

function _ncts_large_run_estimate_status(;
    estimated_wall_seconds,
    max_wall_seconds,
    estimated_rss_bytes,
)
    wall_status = _ncts_wall_estimate_status(
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
    )
    wall_status == "ok" || return wall_status
    rss_status = _ncts_rss_estimate_status(estimated_rss_bytes=estimated_rss_bytes)
    rss_status == "ok" || return rss_status
    return "ok"
end

function _ncts_wall_estimate_status(; estimated_wall_seconds, max_wall_seconds)
    isnothing(estimated_wall_seconds) && return "blocked_missing_wall_estimate"
    estimated_wall_seconds > 0 || return "blocked_invalid_wall_estimate"
    estimated_wall_seconds <= max_wall_seconds ||
        return "blocked_excessive_wall_estimate"
    return "ok"
end

function _ncts_rss_estimate_status(; estimated_rss_bytes)
    isnothing(estimated_rss_bytes) && return "blocked_missing_rss_estimate"
    estimated_rss_bytes > 0 || return "blocked_invalid_rss_estimate"
    return "ok"
end

function _ncts_load_guard_allow_overcommitted(env_prefix::AbstractString)
    global_allow = _ncts_load_guard_parse_bool(
        "NCTS_LOAD_GUARD_ALLOW_OVERCOMMITTED",
        false,
    )
    global_allow && return true
    return _ncts_load_guard_parse_bool(
        "$(env_prefix)_ALLOW_OVERCOMMITTED",
        false,
    )
end

function _ncts_print_large_run_pressure_gate(;
    status::AbstractString,
    reason::AbstractString,
    blocked_status::AbstractString,
    load_status::AbstractString,
    load_reason::AbstractString,
    estimate_status::AbstractString,
    estimate_reason::AbstractString,
    wall_estimate_status::AbstractString,
    wall_estimate_reason::AbstractString,
    rss_estimate_status::AbstractString,
    rss_estimate_reason::AbstractString,
    memory_status::AbstractString,
    memory_reason::AbstractString,
    label::AbstractString,
    load1,
    load5,
    load15,
    cpu_count::Real,
    max_load_per_cpu::Real,
    estimated_wall_seconds,
    max_wall_seconds::Real,
    mem_available_bytes,
    estimated_rss_bytes,
    max_rss_fraction::Real,
    override::Bool,
)
    println("large_run_pressure_gate_status=$status")
    println("large_run_pressure_gate_blocked_status=$blocked_status")
    println("large_run_pressure_gate_load_status=$load_status")
    println("large_run_pressure_gate_estimate_status=$estimate_status")
    println("large_run_pressure_gate_wall_estimate_status=$wall_estimate_status")
    println("large_run_pressure_gate_rss_estimate_status=$rss_estimate_status")
    println("large_run_pressure_gate_memory_status=$memory_status")
    println("large_run_pressure_gate_reason=$reason")
    println("large_run_pressure_gate_load_reason=$load_reason")
    println("large_run_pressure_gate_estimate_reason=$estimate_reason")
    println("large_run_pressure_gate_wall_estimate_reason=$wall_estimate_reason")
    println("large_run_pressure_gate_rss_estimate_reason=$rss_estimate_reason")
    println("large_run_pressure_gate_memory_reason=$memory_reason")
    println("large_run_pressure_gate_label=$label")
    println("large_run_pressure_gate_load1=$(isnothing(load1) ? "unknown" : load1)")
    println("large_run_pressure_gate_load5=$(isnothing(load5) ? "unknown" : load5)")
    println("large_run_pressure_gate_load15=$(isnothing(load15) ? "unknown" : load15)")
    println("large_run_pressure_gate_cpus=$cpu_count")
    println("large_run_pressure_gate_max_load_per_cpu=$max_load_per_cpu")
    println("large_run_pressure_gate_estimated_wall_seconds=$(isnothing(estimated_wall_seconds) ? "unknown" : estimated_wall_seconds)")
    println("large_run_pressure_gate_max_wall_seconds=$max_wall_seconds")
    println("large_run_pressure_gate_mem_available_bytes=$(isnothing(mem_available_bytes) ? "unknown" : mem_available_bytes)")
    println("large_run_pressure_gate_estimated_rss_bytes=$(isnothing(estimated_rss_bytes) ? "unknown" : estimated_rss_bytes)")
    println("large_run_pressure_gate_max_rss_fraction=$max_rss_fraction")
    println("large_run_pressure_gate_override=$override")
    return nothing
end

function _ncts_load_guard_gib_label(bytes)
    isnothing(bytes) && return "unknown"
    return "$(round(Float64(bytes) / 1024.0^3; digits=3)) GiB"
end

function _ncts_large_run_pressure_failure_detail(;
    status::AbstractString,
    env_prefix::AbstractString,
    load1,
    load5,
    load15,
    cpu_count::Real,
    max_load_per_cpu::Real,
    estimated_wall_seconds,
    max_wall_seconds::Real,
    mem_available_bytes,
    estimated_rss_bytes,
    max_rss_fraction::Real,
)
    if status == "blocked_overloaded_remote"
        return "load averages $(load1), $(load5), $(load15) exceed " *
               "$(cpu_count) * $(max_load_per_cpu)"
    elseif status == "blocked_missing_loadavg"
        return "load-average telemetry is unavailable"
    elseif status == "blocked_missing_wall_estimate"
        return "no wall-time estimate was provided; set " *
               "$(env_prefix)_ESTIMATED_WALL_SECONDS or " *
               "NCTS_LOAD_GUARD_ESTIMATED_WALL_SECONDS"
    elseif status == "blocked_invalid_wall_estimate"
        return "wall-time estimate $(estimated_wall_seconds) is not positive"
    elseif status == "blocked_excessive_wall_estimate"
        return "wall-time estimate $(estimated_wall_seconds) exceeds " *
               "the configured maximum $(max_wall_seconds); lower the " *
               "run size or set $(env_prefix)_MAX_WALL_SECONDS / " *
               "NCTS_LOAD_GUARD_MAX_WALL_SECONDS only with fresh evidence"
    elseif status == "blocked_missing_rss_estimate"
        return "no RSS estimate was provided; set " *
               "$(env_prefix)_ESTIMATED_RSS_GIB or " *
               "NCTS_LOAD_GUARD_ESTIMATED_RSS_GIB"
    elseif status == "blocked_invalid_rss_estimate"
        return "RSS estimate $(_ncts_load_guard_gib_label(estimated_rss_bytes)) " *
               "is not positive"
    elseif status == "blocked_missing_memory_telemetry"
        return "available-memory telemetry is unavailable"
    elseif status == "blocked_insufficient_memory"
        return "estimated RSS $(_ncts_load_guard_gib_label(estimated_rss_bytes)) " *
               "exceeds available memory " *
               "$(_ncts_load_guard_gib_label(mem_available_bytes))"
    elseif status == "blocked_insufficient_memory_headroom"
        return "estimated RSS $(_ncts_load_guard_gib_label(estimated_rss_bytes)) " *
               "would consume more than $(max_rss_fraction * 100)% of " *
               "available memory " *
               "$(_ncts_load_guard_gib_label(mem_available_bytes)); lower " *
               "the run size or set $(env_prefix)_MAX_RSS_FRACTION / " *
               "NCTS_LOAD_GUARD_MAX_RSS_FRACTION only with fresh evidence"
    end
    return "status $status"
end

function _ncts_large_run_pressure_status_detail(;
    status::AbstractString,
    env_prefix::AbstractString,
    load1,
    load5,
    load15,
    cpu_count::Real,
    max_load_per_cpu::Real,
    estimated_wall_seconds,
    max_wall_seconds::Real,
    mem_available_bytes,
    estimated_rss_bytes,
    max_rss_fraction::Real,
)
    status == "ok" && return ""
    return _ncts_large_run_pressure_failure_detail(
        status=status,
        env_prefix=env_prefix,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
end

function _ncts_check_large_run_pressure_guard(;
    env_prefix::AbstractString,
    label::AbstractString,
)
    loadavg = _ncts_loadavg_triplet()
    load1, load5, load15 = isnothing(loadavg) ? (nothing, nothing, nothing) : loadavg
    cpu_count = _ncts_load_guard_parse_float(
        "$(env_prefix)_LOAD_GUARD_CPUS",
        _ncts_effective_cpu_count(),
    )
    max_load_per_cpu = _ncts_load_guard_parse_float(
        "$(env_prefix)_MAX_LOAD_PER_CPU",
        1.0,
    )
    load_status = _ncts_large_run_pressure_status(
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
    )
    estimated_wall_seconds = _ncts_load_guard_optional_float(
        "$(env_prefix)_ESTIMATED_WALL_SECONDS",
    )
    isnothing(estimated_wall_seconds) && (estimated_wall_seconds =
        _ncts_load_guard_optional_float("NCTS_LOAD_GUARD_ESTIMATED_WALL_SECONDS"))
    max_wall_seconds = _ncts_load_guard_parse_float(
        "$(env_prefix)_MAX_WALL_SECONDS",
        _ncts_load_guard_parse_float("NCTS_LOAD_GUARD_MAX_WALL_SECONDS", 1800.0),
    )
    mem_available_bytes = _ncts_mem_available_bytes()
    estimated_rss_gib = _ncts_load_guard_optional_float(
        "$(env_prefix)_ESTIMATED_RSS_GIB",
    )
    isnothing(estimated_rss_gib) && (estimated_rss_gib =
        _ncts_load_guard_optional_float("NCTS_LOAD_GUARD_ESTIMATED_RSS_GIB"))
    estimated_rss_bytes = isnothing(estimated_rss_gib) ?
        nothing :
        estimated_rss_gib * 1024.0^3
    max_rss_fraction = _ncts_load_guard_parse_float(
        "$(env_prefix)_MAX_RSS_FRACTION",
        _ncts_load_guard_parse_float("NCTS_LOAD_GUARD_MAX_RSS_FRACTION", 0.8),
    )
    wall_estimate_status = _ncts_wall_estimate_status(
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
    )
    rss_estimate_status = _ncts_rss_estimate_status(
        estimated_rss_bytes=estimated_rss_bytes,
    )
    estimate_status = wall_estimate_status == "ok" ? rss_estimate_status : wall_estimate_status
    memory_status = _ncts_memory_pressure_status(
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )

    status = load_status
    status == "ok" && estimate_status != "ok" && (status = estimate_status)
    status == "ok" && memory_status != "ok" && (status = memory_status)

    override = _ncts_load_guard_allow_overcommitted(env_prefix)
    load_reason = _ncts_large_run_pressure_status_detail(
        status=load_status,
        env_prefix=env_prefix,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    wall_estimate_reason = _ncts_large_run_pressure_status_detail(
        status=wall_estimate_status,
        env_prefix=env_prefix,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    rss_estimate_reason = _ncts_large_run_pressure_status_detail(
        status=rss_estimate_status,
        env_prefix=env_prefix,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    estimate_reason = wall_estimate_status == "ok" ?
        rss_estimate_reason :
        wall_estimate_reason
    memory_reason = _ncts_large_run_pressure_status_detail(
        status=memory_status,
        env_prefix=env_prefix,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    detail = _ncts_large_run_pressure_status_detail(
        status=status,
        env_prefix=env_prefix,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    printed_status = override ? "override_overcommitted_remote" : status
    printed_reason = override && status != "ok" ?
        "override accepted despite $status: $detail" :
        detail
    _ncts_print_large_run_pressure_gate(
        status=printed_status,
        reason=printed_reason,
        blocked_status=status,
        load_status=load_status,
        load_reason=load_reason,
        estimate_status=estimate_status,
        estimate_reason=estimate_reason,
        wall_estimate_status=wall_estimate_status,
        wall_estimate_reason=wall_estimate_reason,
        rss_estimate_status=rss_estimate_status,
        rss_estimate_reason=rss_estimate_reason,
        memory_status=memory_status,
        memory_reason=memory_reason,
        label=label,
        load1=load1,
        load5=load5,
        load15=load15,
        cpu_count=cpu_count,
        max_load_per_cpu=max_load_per_cpu,
        estimated_wall_seconds=estimated_wall_seconds,
        max_wall_seconds=max_wall_seconds,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
        override=override,
    )

    override && return nothing
    status == "ok" && return nothing
    throw(ArgumentError(
        "Refusing $label because the large-run pressure gate returned " *
        "$status: $detail.  Wait for the remote load/memory risk to " *
        "drop or remove this run from the live plan; use " *
        "NCTS_LOAD_GUARD_ALLOW_OVERCOMMITTED=true only on isolated hardware."
    ))
end
