#!/usr/bin/env julia

# Small guarded QMBCertify reference-run harness for the 1D Heisenberg parity
# fixtures.  This is for reviewed reference data collection, not broad scaling.
#
# Usage:
#   NCTS_QMB_PROFILE=A0 NCTS_QMB_NS=10 julia --project=. --startup-file=no perf/qmbcertify_reference_runs.jl
#   NCTS_QMBCERTIFY_PATH=/path/to/QMBCertify NCTS_QMB_PROFILE=A1 NCTS_QMB_NS=10 julia --project=. --startup-file=no perf/qmbcertify_reference_runs.jl
#
# Safety:
#   By default this refuses L > 13.  The required fixture list includes L=20
#   and L=30, but those should only be run deliberately with explicit wall/RSS
#   estimates and safe load/memory telemetry.

Base.include(@__MODULE__, joinpath(@__DIR__, "shared_load_guard.jl"))

function _qmb_reference_preimport_large_run_pressure_guard()
    _ncts_load_guard_parse_bool("NCTS_QMB_BOOTSTRAP_ONLY", false) &&
        return nothing
    _ncts_load_guard_parse_bool("NCTS_QMB_ALLOW_LARGE", false) ||
        return nothing
    lengths_label = strip(get(ENV, "NCTS_QMB_NS", "10"))
    _ncts_check_large_run_pressure_guard(
        env_prefix="NCTS_QMB",
        label="QMBCertify reference Ls=$lengths_label",
    )
    return nothing
end

_qmb_reference_preimport_large_run_pressure_guard()

using Dates
using LinearAlgebra
using Pkg
using Printf
using TOML

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

function _check_size_guard(lengths::AbstractVector{Int}; pressure_guard::Bool=true)
    allow_large = _parse_bool_env("NCTS_QMB_ALLOW_LARGE", false)
    max_l = parse(Int, get(ENV, "NCTS_QMB_MAX_L", "13"))
    if allow_large
        pressure_guard && _ncts_check_large_run_pressure_guard(
            env_prefix="NCTS_QMB",
            label="QMBCertify reference Ls=$(join(lengths, ","))",
        )
        return nothing
    end
    for len in lengths
        len <= max_l || throw(ArgumentError(
            "Refusing L=$len because NCTS_QMB_ALLOW_LARGE=false and " *
            "NCTS_QMB_MAX_L=$max_l.  Set NCTS_QMB_ALLOW_LARGE=true " *
            "only for intentional reviewed QMBCertify reference runs."
        ))
    end
    return nothing
end

_requested_lengths() = _parse_int_list_env("NCTS_QMB_NS", "10")

function _qmb_estimated_wall_seconds()
    value = _ncts_load_guard_optional_float("NCTS_QMB_ESTIMATED_WALL_SECONDS")
    isnothing(value) &&
        (value = _ncts_load_guard_optional_float("NCTS_LOAD_GUARD_ESTIMATED_WALL_SECONDS"))
    return value
end

function _qmb_estimated_rss_gib()
    value = _ncts_load_guard_optional_float("NCTS_QMB_ESTIMATED_RSS_GIB")
    isnothing(value) &&
        (value = _ncts_load_guard_optional_float("NCTS_LOAD_GUARD_ESTIMATED_RSS_GIB"))
    return value
end

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

function _rss_bytes()
    rss = _proc_status_kb("VmRSS")
    return isnothing(rss) ? nothing : rss * 1024
end

function _hwm_bytes()
    hwm = _proc_status_kb("VmHWM")
    return isnothing(hwm) ? nothing : hwm * 1024
end

function _mem_total_bytes()
    path = "/proc/meminfo"
    isfile(path) || return nothing
    for line in eachline(path)
        startswith(line, "MemTotal:") || continue
        parts = split(line)
        length(parts) >= 2 || return nothing
        return parse(Int, parts[2]) * 1024
    end
    return nothing
end

function _bytes_string(bytes)
    isnothing(bytes) && return "n/a"
    return _fmt_bytes(bytes)
end

function _repo_commit()
    env_commit = get(ENV, "NCTS_NCTSSOS_COMMIT", "")
    !isempty(strip(env_commit)) && return strip(env_commit)

    repo = normpath(joinpath(@__DIR__, ".."))
    try
        return readchomp(pipeline(`git -C $repo rev-parse HEAD`; stderr=devnull))
    catch
        return nothing
    end
end

function _allow_missing_reference_commit()
    return _parse_bool_env("NCTS_QMB_ALLOW_MISSING_NCTSSOS_COMMIT", false)
end

function _reviewed_reference_commit(commit; allow_missing::Bool)
    commit !== nothing && return commit
    allow_missing && return nothing
    throw(ArgumentError(
        "Reviewed QMBCertify reference runs require NCTSSoS commit provenance. " *
        "Set NCTS_NCTSSOS_COMMIT to the local repository commit when running " *
        "through easy-ssh, or set NCTS_QMB_ALLOW_MISSING_NCTSSOS_COMMIT=true " *
        "only for unreviewed probes."
    ))
end

function _package_version(name::AbstractString)
    for dependency in values(Pkg.dependencies())
        dependency.name == name || continue
        return isnothing(dependency.version) ? nothing : string(dependency.version)
    end
    return nothing
end

_toml_string(value) = "\"" * replace(string(value), "\\" => "\\\\", "\"" => "\\\"") * "\""
_toml_optional_string(value) = (ismissing(value) || isnothing(value)) ? missing : _toml_string(value)

function _fixture_date_id()
    return replace(string(Dates.today()), "-" => "_")
end

function _environment_id()
    raw = strip(get(ENV, "NCTS_QMB_ENVIRONMENT_ID", ""))
    !isempty(raw) && return raw
    return "qmb_" * _fixture_date_id()
end

function _fixture_for_profile(profile::AbstractString)
    profile == "A0" && return joinpath(@__DIR__, "..", "test", "data", "expectations", "heisenberg_qmbcertify_base.toml")
    profile in ("A1", "A2") && return joinpath(@__DIR__, "..", "test", "data", "expectations", "heisenberg_qmbcertify_rdm.toml")
    throw(ArgumentError("Unsupported NCTS_QMB_PROFILE=$(repr(profile)); expected A0, A1, or A2."))
end

function _profile_fixture(profile::AbstractString)
    return TOML.parsefile(_fixture_for_profile(profile))
end

function _profile_args(profile::AbstractString)
    fixture = _profile_fixture(profile)
    definitions = fixture["profile_definitions"]
    haskey(definitions, profile) || throw(ArgumentError("Fixture has no profile_definitions.$profile."))
    return definitions[profile]["command_args"]
end

function _pending_case(fixture, profile::AbstractString, len::Integer)
    for case in get(fixture["required_runs"], "pending_cases", Any[])
        String(case["profile"]) == profile || continue
        Int(case["length"]) == Int(len) || continue
        return case
    end
    return nothing
end

function _profile_probe_by_id(fixture, evidence_id::AbstractString)
    for probe in get(fixture, "profile_probes", Any[])
        String(probe["id"]) == evidence_id && return probe
    end
    return nothing
end

function _check_positive_probe_number(probe, field::AbstractString, evidence_id::AbstractString)
    haskey(probe, field) || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) is " *
        "missing required probe field $(repr(field)).",
    ))
    value = probe[field]
    value isa Real && value > 0 || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) has " *
        "$(field)=$(repr(value)); expected a positive number.",
    ))
    return nothing
end

function _check_accepted_profile_variant_restore_evidence(
    fixture,
    profile::AbstractString,
    evidence_id::AbstractString,
)
    probe = _profile_probe_by_id(fixture, evidence_id)
    probe !== nothing || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) was not " *
        "found in the fixture profile_probes table.",
    ))

    String(get(probe, "status", "")) == "probe_run" || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) does not " *
        "name a probe_run row.",
    ))

    parent_profile = String(get(probe, "parent_profile", ""))
    parent_profile == profile || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) belongs " *
        "to parent_profile=$(repr(parent_profile)), not profile=$(repr(profile)).",
    ))

    overrides_summary = strip(String(get(probe, "overrides_summary", "")))
    !isempty(overrides_summary) && overrides_summary != "none" || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) does not " *
        "name an exploratory profile variant; overrides_summary is " *
        "$(repr(overrides_summary)).",
    ))

    _check_positive_probe_number(probe, "total_wall_seconds", evidence_id)
    _check_positive_probe_number(probe, "peak_rss_bytes", evidence_id)
    _check_positive_probe_number(probe, "estimated_wall_seconds", evidence_id)
    _check_positive_probe_number(probe, "estimated_rss_gib", evidence_id)

    termination = String(get(probe, "termination_status", ""))
    accepted = Set(String.(get(
        probe,
        "accepted_termination_statuses",
        _reviewed_solve_status_policy(profile),
    )))
    termination in accepted || throw(ArgumentError(
        "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID=$(repr(evidence_id)) has " *
        "termination_status=$(repr(termination)), outside accepted statuses " *
        "$(collect(accepted)).",
    ))
    return nothing
end

function _check_deleted_reference_case_restore_gate(
    profile::AbstractString,
    lengths::AbstractVector{Int},
)
    fixture = _profile_fixture(profile)
    evidence_id = strip(get(ENV, "NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID", ""))
    evidence_checked = false
    for len in lengths
        case = _pending_case(fixture, profile, len)
        case === nothing && continue
        get(case, "run_queue_status", "queued") == "deleted_until_evidence_gate" ||
            continue
        if !isempty(evidence_id)
            if !evidence_checked
                _check_accepted_profile_variant_restore_evidence(
                    fixture,
                    profile,
                    evidence_id,
                )
                evidence_checked = true
            end
            continue
        end
        restore_gate = String(get(case, "restore_gate", "fresh_wall_rss_load_estimate"))
        throw(ArgumentError(
            "Refusing deleted QMBCertify reviewed-reference row $(profile)_L$(len). " *
            "Fixture restore_gate=$restore_gate requires explicit accepted profile " *
            "variant evidence.  Set NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID only " *
            "after an exploratory perf/qmbcertify_profile_probe.jl row has an " *
            "accepted termination status and the run also has fresh wall/RSS/load estimates."
        ))
    end
    return nothing
end

function _profile_source(profile::AbstractString)
    fixture = _profile_fixture(profile)
    return fixture["source"]
end

function _default_qmbcertify_path(source)
    local_path = String(get(source, "local_path", "/Users/exaclior/projects/QMBCertify"))
    isdir(local_path) && return local_path
    return normpath(joinpath(pwd(), "..", "QMBCertify"))
end

function _heisenberg_input(args)
    get(args, "supp", nothing) == "[UInt16[1;4]]" || throw(ArgumentError(
        "This harness currently supports only supp=[UInt16[1;4]]."
    ))
    get(args, "coe", nothing) == "[3/4]" || throw(ArgumentError(
        "This harness currently supports only coe=[3/4]."
    ))
    return Vector{UInt16}[UInt16[1, 4]], [3 / 4]
end

function _lol_value(raw, len::Int)
    raw isa Integer && return Int(raw)
    raw isa AbstractString || throw(ArgumentError("Unsupported lol value $(repr(raw))."))
    raw == "L" && return len
    raw == "Int(L/2)" && return Int(len / 2)
    return parse(Int, raw)
end

function _rdm_value(raw)
    raw isa Bool && return raw
    raw isa Integer && return Int(raw)
    raw isa AbstractString && lowercase(raw) == "false" && return false
    raw isa AbstractString && return parse(Int, raw)
    throw(ArgumentError("Unsupported rdm value $(repr(raw))."))
end

function _https_remote(remote::AbstractString)
    if startswith(remote, "git@github.com:")
        repo = replace(chop(remote; head=length("git@github.com:"), tail=0), ".git" => "")
        return "https://github.com/$repo.git"
    end
    return remote
end

function _checkout_commit(path::AbstractString)
    isdir(joinpath(path, ".git")) || return nothing
    try
        return readchomp(`git -C $path rev-parse HEAD`)
    catch
        return nothing
    end
end

function _bootstrap_qmbcertify!(path::AbstractString, source)
    isdir(path) && return nothing
    _parse_bool_env("NCTS_QMB_BOOTSTRAP", false) || throw(ArgumentError(
        "QMBCertify checkout not found at $(repr(path)).  Set NCTS_QMBCERTIFY_PATH, " *
        "or set NCTS_QMB_BOOTSTRAP=true to clone the pinned reference checkout."
    ))

    remote = String(get(source, "remote", "git@github.com:wangjie212/QMBCertify.git"))
    url = String(get(ENV, "NCTS_QMBCERTIFY_URL", _https_remote(remote)))
    commit = String(source["commit"])
    mkpath(dirname(path))
    println("- bootstrapping QMBCertify from `$url`")
    run(`git clone --quiet $url $path`)
    run(`git -C $path checkout --quiet $commit`)
    return nothing
end

function _verify_qmbcertify_commit(path::AbstractString, source)
    expected = String(source["commit"])
    actual = _checkout_commit(path)
    actual === nothing && return nothing
    actual == expected || throw(ArgumentError(
        "QMBCertify checkout at $(repr(path)) is at commit $actual, expected $expected."
    ))
    return nothing
end

function _activate_qmbcertify_environment(path::AbstractString)
    _parse_bool_env("NCTS_QMB_ACTIVATE", true) || return nothing
    Pkg.activate(path; io=devnull)
    instantiate = _parse_bool_env("NCTS_QMB_INSTANTIATE", true)
    instantiate && Pkg.instantiate(; io=devnull)
    return nothing
end

function _load_qmbcertify(path::AbstractString, source)
    _bootstrap_qmbcertify!(path, source)
    _verify_qmbcertify_commit(path, source)
    _activate_qmbcertify_environment(path)
    pushfirst!(LOAD_PATH, path)
    return Base.require(Base.PkgId(Base.UUID("71df4310-1531-11eb-256d-190d0a52716f"), "QMBCertify"))
end

function _extract_seconds(log::AbstractString, pattern::Regex)
    match_result = match(pattern, log)
    match_result === nothing && return missing
    return parse(Float64, match_result.captures[1])
end

function _qmb_log_timings(log::AbstractString)
    basis_time = _extract_seconds(log, r"Obtained the monomial basis in ([0-9.eE+-]+) seconds\.")
    block_time = _extract_seconds(log, r"Finished block-diagonalization in ([0-9.eE+-]+) seconds\.")
    lso_time = _extract_seconds(log, r"Added linear state optimality constraints in ([0-9.eE+-]+) seconds\.")
    pso_time = _extract_seconds(log, r"Added PSD state optimality constraints in ([0-9.eE+-]+) seconds\.")
    rdm_time = _extract_seconds(log, r"Added the rdm constraint in ([0-9.eE+-]+) seconds\.")
    solve_time = _extract_seconds(log, r"SDP solving time: ([0-9.eE+-]+) seconds\.")
    sdp_size = match(r"SDP size: n = ([0-9]+), m = ([0-9]+)", log)
    max_block = sdp_size === nothing ? missing : parse(Int, sdp_size.captures[1])
    moment_count = sdp_size === nothing ? missing : parse(Int, sdp_size.captures[2])
    termination = match(r"termination status: ([A-Z_]+)", log)
    solution = match(r"solution status: ([A-Z_]+)", log)
    termination_status = termination === nothing ?
        (occursin(r"optimum = [-+0-9.eE]+", log) ? "OPTIMAL" : missing) :
        termination.captures[1]
    solution_status = solution === nothing ? missing : solution.captures[1]
    return (;
        basis_time,
        block_time,
        lso_time,
        pso_time,
        rdm_time,
        solve_time,
        max_block,
        moment_count,
        termination_status,
        solution_status,
    )
end

_maybe_zero(::Missing) = 0.0
_maybe_zero(x) = x

struct NonOptimalReferenceSolve <: Exception
    objective
    timings
end

function Base.showerror(io::IO, err::NonOptimalReferenceSolve)
    print(
        io,
        "QMBCertify reference solve did not terminate optimally: " *
        "termination_status=$(repr(err.timings.termination_status)), " *
        "solution_status=$(repr(err.timings.solution_status)).",
    )
end

_reviewed_solve_status_policy() = ("OPTIMAL", "ALMOST_OPTIMAL")

function _reviewed_solve_status_policy(profile::AbstractString)
    fixture = _profile_fixture(profile)
    policy = fixture["reviewed_reference_status_policy"]
    statuses = String.(policy["accepted_termination_statuses"])
    isempty(statuses) && throw(ArgumentError(
        "reviewed_reference_status_policy.accepted_termination_statuses is empty for profile $profile.",
    ))
    return Tuple(statuses)
end

function _check_reviewed_solve_status!(
    timings,
    objective;
    accepted_statuses=_reviewed_solve_status_policy(),
)
    termination = timings.termination_status
    !ismissing(termination) && termination in accepted_statuses &&
        return nothing
    throw(NonOptimalReferenceSolve(objective, timings))
end

function _histogram_pairs(values::AbstractVector{<:Integer})
    counts = Dict{Int,Int}()
    for value in values
        key = Int(value)
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by=first)
end

function _rdm_block_sizes(rdm)
    rdm == 8 && return [72, 64, 56]
    rdm == 9 && return [136, 120]
    rdm == 10 && return [256, 272, 240]
    return Int[]
end

function _qmb_base_block_sizes(basis, len::Int)
    sizes = Int[]
    for parity in 1:2
        k = length(basis[parity]) ÷ len
        for sector in 1:(Int(len / 2) + 1)
            if parity == 1 && sector == 1
                push!(sizes, k + 1)
            elseif sector == 1 || sector == Int(len / 2) + 1
                push!(sizes, k)
            else
                push!(sizes, 2k)
            end
        end
    end
    return sizes
end

function _qmb_pso_block_sizes(basis, len::Int, pso::Int)
    pso <= 0 && return Int[]
    sizes = Int[]
    for parity in 1:2
        pso_basis_count = count(word -> length(word) <= pso, basis[parity])
        k = pso_basis_count ÷ len
        for sector in 1:(Int(len / 2) + 1)
            push!(sizes, (sector == 1 || sector == Int(len / 2) + 1) ? k : 2k)
        end
    end
    return sizes
end

function _block_variable_count(sizes::AbstractVector{<:Integer})
    return sum(Int(size) * (Int(size) + 1) ÷ 2 for size in sizes; init=0)
end

function _run_one(
    len::Int,
    profile::AbstractString,
    args,
    qmb::Module;
    status_policy_profile::AbstractString=profile,
)
    supp, coe = _heisenberg_input(args)
    d = Int(args["d"])
    pso = Int(get(args, "pso", 0))
    rdm = _rdm_value(get(args, "rdm", false))
    lso = get(args, "lso", true)
    lol = _lol_value(get(args, "lol", "Int(L/2)"), len)
    extra = Int(get(args, "extra", 0))
    lattice = String(get(args, "lattice", "chain"))
    correlation = Bool(get(args, "correlation", false))
    quiet = Bool(get(args, "quiet", false))
    tolerances = get(args, "mosek_tolerances", [1.0e-8, 1.0e-8, 1.0e-8])
    threads = parse(Int, get(ENV, "NCTS_QMB_MOSEK_THREADS", string(Int(get(args, "mosek_threads", 0)))))

    result = nothing
    log_path = tempname()
    run_stats = @timed begin
        open(log_path, "w") do log_io
            redirect_stdout(log_io) do
                result = Base.invokelatest(
                    qmb.GSB,
                    supp,
                    coe,
                    len,
                    d;
                    lso,
                    lol,
                    pso,
                    rdm,
                    extra,
                    lattice,
                    correlation,
                    QUIET=quiet,
                    mosek_setting=Base.invokelatest(
                        qmb.mosek_para,
                        Float64(tolerances[1]),
                        Float64(tolerances[2]),
                        Float64(tolerances[3]),
                        threads,
                    ),
                )
            end
        end
    end
    log = read(log_path, String)
    rm(log_path; force=true)
    print(log)
    opt, data = result
    timings = _qmb_log_timings(log)
    _check_reviewed_solve_status!(
        timings,
        opt;
        accepted_statuses=_reviewed_solve_status_policy(status_policy_profile),
    )
    base_blocks = _qmb_base_block_sizes(data.basis, len)
    pso_blocks = _qmb_pso_block_sizes(data.basis, len, pso)
    rdm_blocks = _rdm_block_sizes(rdm)
    psd_block_sizes = vcat(base_blocks, pso_blocks, rdm_blocks)
    build_time = _maybe_zero(timings.basis_time) +
        _maybe_zero(timings.block_time) +
        _maybe_zero(timings.lso_time) +
        _maybe_zero(timings.pso_time) +
        _maybe_zero(timings.rdm_time)

    return (
        profile=profile,
        length=len,
        objective=opt,
        objective_per_site=opt,
        objective_total_estimate=opt * len,
        total_wall_seconds=run_stats.time,
        build_time_seconds=ismissing(timings.basis_time) || ismissing(timings.block_time) ? missing : build_time,
        solve_time_seconds=timings.solve_time,
        termination_status=timings.termination_status,
        solution_status=timings.solution_status,
        moment_count=timings.moment_count,
        reported_max_block=timings.max_block,
        psd_block_count=length(psd_block_sizes),
        psd_max_block=isempty(psd_block_sizes) ? 0 : maximum(psd_block_sizes),
        psd_block_histogram=_histogram_pairs(psd_block_sizes),
        psd_scalar_variables=_block_variable_count(psd_block_sizes),
        peak_rss_bytes=_hwm_bytes(),
        current_rss_bytes=_rss_bytes(),
    )
end

function _print_toml_field(name::AbstractString, value)
    if ismissing(value) || isnothing(value)
        println("# $name = missing")
    else
        println("$name = $value")
    end
end

function _print_environment_fixture_block(qmb_path::AbstractString; table::AbstractString="environment")
    println("\n## Environment Fixture Block")
    println("\n[$table]")
    _print_toml_field("nctssos_commit", isnothing(_repo_commit()) ? nothing : _toml_string(_repo_commit()))
    println("julia_version = $(_toml_string(VERSION))")
    _print_toml_field("mosek_version", isnothing(_package_version("Mosek")) ? nothing : _toml_string(_package_version("Mosek")))
    println("cpu_model = $(_toml_string(Sys.cpu_info()[1].model))")
    _print_toml_field("ram_bytes", _mem_total_bytes())
    println("thread_count = $(Threads.nthreads())")
    println("blas_vendor = $(_toml_string(BLAS.vendor()))")
    println("qmbcertify_path = $(_toml_string(qmb_path))")
end

function _print_failed_attempt_rows(rows, qmb_path::AbstractString, env_id::AbstractString)
    isempty(rows) && return nothing
    println("\n## Failed Attempt Rows")
    for row in rows
        timings = row.timings
        println("\n[[failed_attempts]]")
        println("id = \"$(row.profile)_L$(row.length)_$(_fixture_date_id())\"")
        println("case_id = \"$(row.profile)_L$(row.length)\"")
        println("profile = \"$(row.profile)\"")
        println("length = $(row.length)")
        println("attempted_at = \"$(Dates.today())\"")
        println("status = \"rejected_nonoptimal_solve\"")
        println("reason = \"Guarded harness rejects non-optimal reference rows.\"")
        println("environment = \"$env_id\"")
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
        println(
            "notes = " *
            _toml_string(
                "Generated by perf/qmbcertify_reference_runs.jl; confirm the intended parity profile before moving this case into reviewed cases.",
            ),
        )
    end
    _print_environment_fixture_block(qmb_path; table="failed_attempt_environments.$env_id")
end

function _print_rollup(rows, qmb_path::AbstractString, env_id::AbstractString; failed_attempts=NamedTuple[])
    println("\n## Rollup")
    println("\n| profile | L | objective/site | total wall (s) | build (s) | solve (s) | moments | PSD blocks | max PSD block | PSD vars | peak RSS |")
    println("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for row in rows
        @printf(
            "| %s | %d | %.12g | %.3f | %s | %s | %s | %d | %d | %d | %s |\n",
            row.profile,
            row.length,
            row.objective_per_site,
            row.total_wall_seconds,
            ismissing(row.build_time_seconds) ? "missing" : @sprintf("%.3f", row.build_time_seconds),
            ismissing(row.solve_time_seconds) ? "missing" : @sprintf("%.3f", row.solve_time_seconds),
            ismissing(row.moment_count) ? "missing" : string(row.moment_count),
            row.psd_block_count,
            row.psd_max_block,
            row.psd_scalar_variables,
            _bytes_string(row.peak_rss_bytes),
        )
    end

    if !isempty(rows)
        _print_environment_fixture_block(qmb_path; table="reference_environments.$env_id")
    end

    println("\n## Fixture Rows")
    for row in rows
        println("\n[[cases]]")
        println("id = \"$(row.profile)_L$(row.length)\"")
        println("profile = \"$(row.profile)\"")
        println("length = $(row.length)")
        println("status = \"reviewed_run\"")
        println("environment = \"$env_id\"")
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
    end
    _print_failed_attempt_rows(failed_attempts, qmb_path, env_id)
end

function main()
    profile = String(get(ENV, "NCTS_QMB_PROFILE", "A0"))
    env_id = _environment_id()
    args = _profile_args(profile)
    source = _profile_source(profile)
    qmb_path = String(get(ENV, "NCTS_QMBCERTIFY_PATH", _default_qmbcertify_path(source)))
    bootstrap_only = _parse_bool_env("NCTS_QMB_BOOTSTRAP_ONLY", false)
    requested_lengths = _requested_lengths()
    _check_size_guard(requested_lengths; pressure_guard=!bootstrap_only)
    bootstrap_only || _check_deleted_reference_case_restore_gate(profile, requested_lengths)

    println("# QMBCertify guarded reference runs")
    println("\n- generated: `$(Dates.now())`")
    println("- profile: `$profile`")
    println("- L list: `$(join(requested_lengths, ","))`")
    println("- QMBCertify path: `$qmb_path`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- initial RSS: `$(_bytes_string(_rss_bytes()))`")
    println("- environment id: `$env_id`")
    println("- bootstrap_only: `$bootstrap_only`")

    if !bootstrap_only
        _reviewed_reference_commit(_repo_commit(); allow_missing=_allow_missing_reference_commit())
    end

    qmb = _load_qmbcertify(qmb_path, source)
    if bootstrap_only
        println("\nBootstrap complete; no reference solve requested.")
        _print_environment_fixture_block(qmb_path; table="reference_environments.$env_id")
        return nothing
    end

    rows = NamedTuple[]
    failed_attempts = NamedTuple[]
    failed = false
    for len in requested_lengths
        println("\n## L = $len")
        try
            push!(rows, _run_one(len, profile, args, qmb))
        catch err
            failed = true
            if err isa NonOptimalReferenceSolve
                push!(failed_attempts, (;
                    profile,
                    length=len,
                    objective_reported=err.objective,
                    timings=err.timings,
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
    _print_rollup(rows, qmb_path, env_id; failed_attempts)
    failed && exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
