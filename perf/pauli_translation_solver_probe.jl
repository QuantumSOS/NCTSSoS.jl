#!/usr/bin/env julia

# Small-N probe for Pauli translation solver lowering choices.
#
# Default mode constructs the direct-linear relaxation and lowers it to a JuMP
# model without solving.  Set NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true
# to construct the QMBCertify source-family base instead.  Set
# NCTS_SOLVER_PROBE_LOWER_MODEL=false to stop after construction.  Set
# NCTS_SOLVER_PROBE_DUALIZE=true to build the SOS dual model instead.  Set
# NCTS_SOLVER_PROBE_SOLVE=true to optimize.
# Set NCTS_SOLVER_PROBE_EMIT_FIXTURE=true to print a copyable TOML probe row.
# The script refuses N > 13 unless NCTS_SOLVER_PROBE_ALLOW_LARGE=true and the
# shared pressure gate sees explicit wall/RSS estimates plus safe load/memory
# telemetry.
#
# QMBCertify source-base construction-only smoke:
#
#   NCTS_SOLVER_PROBE_ALLOW_LARGE=true \
#   NCTS_SOLVER_PROBE_N=20 NCTS_SOLVER_PROBE_ORDER=4 \
#   NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true \
#   NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_EXTRA=9 \
#   NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true \
#   NCTS_SOLVER_PROBE_LOWER_MODEL=false \
#   NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=<estimate-before-launch> \
#   NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=<estimate-before-launch> \
#   julia --project=. --startup-file=no perf/pauli_translation_solver_probe.jl
#
# Finite-axis A2 dual/Mosek smoke, run from the docs environment so MosekTools
# is available:
#
#   NCTS_SOLVER_PROBE_N=8 NCTS_SOLVER_PROBE_ORDER=4 \
#   NCTS_SOLVER_PROBE_DUALIZE=true NCTS_SOLVER_PROBE_OPTIMIZER=mosek \
#   NCTS_SOLVER_PROBE_MOSEK_THREADS=1 \
#   NCTS_SOLVER_PROBE_SOLVE=true NCTS_SOLVER_PROBE_REFLECTION=true \
#   NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true \
#   NCTS_SOLVER_PROBE_AXIS_SYMMETRY=true \
#   NCTS_SOLVER_PROBE_AXIS_EQUALITIES=true \
#   NCTS_SOLVER_PROBE_AXIS_QUOTIENT=true \
#   NCTS_SOLVER_PROBE_RDM_K=8 \
#   NCTS_SOLVER_PROBE_RDM_DECOMPOSITION=qmbcertify \
#   NCTS_SOLVER_PROBE_RDM_SUPPORT=extend \
#   NCTS_SOLVER_PROBE_LSO_WIDTH=7 \
#   NCTS_SOLVER_PROBE_LSO_MODE=qmbcertify \
#   NCTS_SOLVER_PROBE_PSO_WIDTH=3 \
#   julia --project=docs --startup-file=no -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); include("perf/pauli_translation_solver_probe.jl"); main()'
#
# Add NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true and
# NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT=true to probe the exact invariant
# moment-coordinate quotient.  Large runs still require explicit wall/RSS
# estimates and the shared pressure gate.

Base.include(@__MODULE__, joinpath(@__DIR__, "shared_load_guard.jl"))

function _solver_probe_preimport_large_run_pressure_guard()
    _ncts_load_guard_parse_bool("NCTS_SOLVER_PROBE_ALLOW_LARGE", false) ||
        return nothing
    n_label = strip(get(ENV, "NCTS_SOLVER_PROBE_N", "6"))
    _ncts_check_large_run_pressure_guard(
        env_prefix="NCTS_SOLVER_PROBE",
        label="solver probe N=$n_label",
    )
    return nothing
end

_solver_probe_preimport_large_run_pressure_guard()

using Dates
using LinearAlgebra
using Pkg
using Printf
using JuMP
using NCTSSoS
import Clarabel

const _SOLVER_PROBE_HARNESS = "perf/pauli_translation_solver_probe.jl"

function _solver_probe_comparison_scope(n_sites::Integer)
    scale = n_sites >= 20 ? "large-N" : "small-N"
    return "$scale NCTSSoS backend evidence; not a reviewed QMBCertify parity run"
end

if lowercase(strip(get(ENV, "NCTS_SOLVER_PROBE_OPTIMIZER", "clarabel"))) == "mosek"
    @eval using MosekTools
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

function _rss_string()
    rss = _proc_status_kb("VmRSS")
    hwm = _proc_status_kb("VmHWM")
    isnothing(rss) && isnothing(hwm) && return "n/a"
    rss_s = isnothing(rss) ? "n/a" : _fmt_bytes(rss * 1024)
    hwm_s = isnothing(hwm) ? "n/a" : _fmt_bytes(hwm * 1024)
    return "rss=$rss_s, hwm=$hwm_s"
end

function _parse_bool_env(name::AbstractString, default::Bool)
    raw = lowercase(strip(get(ENV, name, string(default))))
    raw in ("true", "1", "yes", "y") && return true
    raw in ("false", "0", "no", "n") && return false
    throw(ArgumentError("$name must be a boolean value, got $(repr(raw))."))
end

function _parse_lower_model_env()
    if haskey(ENV, "NCTS_SOLVER_PROBE_LOWER_MODEL")
        lower_model = _parse_bool_env("NCTS_SOLVER_PROBE_LOWER_MODEL", true)
        if haskey(ENV, "NCTS_SOLVER_PROBE_LOWER")
            legacy_lower_model = _parse_bool_env("NCTS_SOLVER_PROBE_LOWER", true)
            lower_model == legacy_lower_model || throw(ArgumentError(
                "NCTS_SOLVER_PROBE_LOWER_MODEL and NCTS_SOLVER_PROBE_LOWER disagree."
            ))
        end
        return lower_model
    end
    return _parse_bool_env("NCTS_SOLVER_PROBE_LOWER", true)
end

function _parse_optional_int_env(name::AbstractString)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    lowercase(raw) in ("nothing", "none", "false", "0") && return nothing
    return parse(Int, raw)
end

function _parse_optional_string_env(name::AbstractString)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    lowercase(raw) in ("nothing", "none", "false", "0") && return nothing
    return raw
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
    length(values) == 2 || throw(ArgumentError(
        "$name must contain exactly two integers, got $values."
    ))
    return (values[1], values[2])
end

function _parse_symbol_env(name::AbstractString, default::Symbol)
    raw = strip(get(ENV, name, String(default)))
    isempty(raw) && return default
    return Symbol(raw)
end

function _parse_optional_float_env(name::AbstractString)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    lowercase(raw) in ("nothing", "none", "false", "0") && return nothing
    return parse(Float64, raw)
end

function _push_optional_mosek_float_attribute!(
    attrs::Vector{Pair{String,Any}},
    env_name::AbstractString,
    attr_name::AbstractString,
)
    value = _parse_optional_float_env(env_name)
    value === nothing && return attrs
    push!(attrs, attr_name => value)
    return attrs
end

function _check_size_guard(n::Integer)
    allow_large = _parse_bool_env("NCTS_SOLVER_PROBE_ALLOW_LARGE", false)
    max_n = parse(Int, get(ENV, "NCTS_SOLVER_PROBE_MAX_N", "13"))
    if allow_large
        _ncts_check_large_run_pressure_guard(
            env_prefix="NCTS_SOLVER_PROBE",
            label="solver probe N=$n",
        )
        return nothing
    end
    n <= max_n || throw(ArgumentError(
        "Refusing N=$n because NCTS_SOLVER_PROBE_ALLOW_LARGE=false and " *
        "NCTS_SOLVER_PROBE_MAX_N=$max_n."
    ))
    return nothing
end

function _check_base_su2_extend_rdm_solver_probe_guard(;
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool=false,
    su2_symmetry::Bool,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    real_moment_matrix::Bool,
    reflection_symmetry::Bool,
)
    su2_moment_quotient && !base_su2_extend_rdm && throw(ArgumentError(
        "NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT=true requires " *
        "NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true."
    ))
    base_su2_extend_rdm || return nothing
    (
        su2_symmetry &&
        contiguous_rdm_decomposition == :su2 &&
        contiguous_rdm_support == :extend
    ) || throw(ArgumentError(
        "NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true requires " *
        "NCTS_SOLVER_PROBE_SU2=true, NCTS_SOLVER_PROBE_RDM_DECOMPOSITION=su2, " *
        "and NCTS_SOLVER_PROBE_RDM_SUPPORT=extend."
    ))
    real_moment_matrix || !reflection_symmetry || throw(ArgumentError(
        "NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true with " *
        "NCTS_SOLVER_PROBE_REFLECTION=true requires " *
        "NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true."
    ))
    return nothing
end

function _check_deleted_primal_su2_solved_probe_guard(;
    n::Integer,
    solve::Bool,
    dualize::Bool,
    formulation::Symbol,
    su2_symmetry::Bool,
    su2_moment_quotient::Bool=false,
    contiguous_rdm_decomposition::Symbol,
    linear_state_opt_width,
    psd_state_opt_width,
)
    solve || return nothing
    dualize && return nothing
    formulation == :moment_variables || return nothing
    su2_symmetry || return nothing
    su2_moment_quotient && return nothing
    contiguous_rdm_decomposition == :su2 || return nothing
    linear_state_opt_width == 7 || return nothing
    psd_state_opt_width == 3 || return nothing

    min_n = parse(Int, get(ENV, "NCTS_SOLVER_PROBE_DELETED_PRIMAL_SU2_SOLVE_MIN_N", "14"))
    Int(n) >= min_n || return nothing
    _parse_bool_env("NCTS_SOLVER_PROBE_ALLOW_DELETED_PRIMAL_SU2_SOLVE", false) &&
        return nothing

    throw(ArgumentError(
        "Refusing solved primal moment-variable SU(2) RDM/LSO/PSO probe at N=$n. " *
        "The guarded N=12 solve timed out at 20:02.64 with 26,226,188 KiB RSS, " *
        "so N >= $min_n solved probes are deleted under the current formulation. " *
        "Use SOLVE=false for model sizing, reduce the rung, or set " *
        "NCTS_SOLVER_PROBE_ALLOW_DELETED_PRIMAL_SU2_SOLVE=true only after a " *
        "formulation change plus fresh wall/RSS/load evidence."
    ))
end

function _check_qmbcertify_base_solver_probe_guard(;
    qmbcertify_base_construct::Bool,
    sign_symmetry::Bool,
    reflection_symmetry::Bool,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool=false,
    axis_rotation_symmetry::Bool,
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    real_moment_matrix::Bool,
    linear_state_opt_width,
    linear_state_opt_mode::Symbol,
    psd_state_opt_width,
)
    qmbcertify_base_construct || return nothing
    unsupported = String[]
    sign_symmetry || push!(unsupported, "NCTS_SOLVER_PROBE_SIGN=false")
    reflection_symmetry && push!(unsupported, "NCTS_SOLVER_PROBE_REFLECTION=true")
    u1_symmetry && push!(unsupported, "NCTS_SOLVER_PROBE_U1=true")
    su2_symmetry && push!(unsupported, "NCTS_SOLVER_PROBE_SU2=true")
    base_su2_extend_rdm &&
        push!(unsupported, "NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true")
    su2_moment_quotient &&
        push!(unsupported, "NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT=true")
    axis_rotation_symmetry &&
        push!(unsupported, "NCTS_SOLVER_PROBE_AXIS_SYMMETRY=true")
    axis_rotation_equalities &&
        push!(unsupported, "NCTS_SOLVER_PROBE_AXIS_EQUALITIES=true")
    axis_rotation_quotient &&
        push!(unsupported, "NCTS_SOLVER_PROBE_AXIS_QUOTIENT=true")
    isempty(unsupported) || throw(ArgumentError(
        "NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true does not support " *
        join(unsupported, ", ") * "."
    ))

    if !isnothing(contiguous_rdm_k) && contiguous_rdm_k > 0
        contiguous_rdm_decomposition == :qmbcertify || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true with RDM requires " *
            "NCTS_SOLVER_PROBE_RDM_DECOMPOSITION=qmbcertify."
        ))
        contiguous_rdm_support == :extend || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true with RDM requires " *
            "NCTS_SOLVER_PROBE_RDM_SUPPORT=extend."
        ))
        real_moment_matrix || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true with RDM requires " *
            "NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true."
        ))
    end
    if !isnothing(linear_state_opt_width) && linear_state_opt_width > 0
        linear_state_opt_mode == :qmbcertify || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true with LSO requires " *
            "NCTS_SOLVER_PROBE_LSO_MODE=qmbcertify."
        ))
    end
    if !isnothing(psd_state_opt_width) && psd_state_opt_width > 0
        real_moment_matrix || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true with PSO requires " *
            "NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true."
        ))
    end
    return nothing
end

function _timed(label::AbstractString, f)
    GC.gc(true)
    stats = @timed f()
    @printf(
        "%s_seconds=%.9f %s_allocated=%s %s_gc_seconds=%.9f %s\n",
        label,
        stats.time,
        label,
        _fmt_bytes(stats.bytes),
        label,
        stats.gctime,
        _rss_string(),
    )
    flush(stdout)
    return stats.value, stats
end

function _solver_probe_optional_model_attribute(f)
    try
        return f()
    catch
        return nothing
    end
end

function _solver_probe_solve_result(model::JuMP.Model)
    result_count = JuMP.result_count(model)
    result_count > 0 || throw(ArgumentError(
        "Solver returned no result after optimize!.",
    ))
    return (
        status=JuMP.termination_status(model),
        primal_status=JuMP.primal_status(model),
        dual_status=JuMP.dual_status(model),
        raw_status=JuMP.raw_status(model),
        result_count,
        objective=JuMP.objective_value(model),
        objective_bound=_solver_probe_optional_model_attribute(
            () -> JuMP.objective_bound(model),
        ),
        relative_gap=_solver_probe_optional_model_attribute(
            () -> JuMP.MOI.get(JuMP.backend(model), JuMP.MOI.RelativeGap()),
        ),
    )
end

function _optimizer(name::Symbol)
    if name == :clarabel
        return optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    elseif name == :mosek
        isdefined(@__MODULE__, :Mosek) || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_OPTIMIZER=mosek requires MosekTools in the active Julia project."
        ))
        mosek_threads = parse(Int, get(ENV, "NCTS_SOLVER_PROBE_MOSEK_THREADS", "1"))
        mosek_threads > 0 || throw(ArgumentError(
            "NCTS_SOLVER_PROBE_MOSEK_THREADS must be positive, got $mosek_threads."
        ))
        attrs = Pair{String,Any}[
            "QUIET" => _parse_bool_env("NCTS_SOLVER_PROBE_MOSEK_QUIET", true),
            "MSK_IPAR_NUM_THREADS" => mosek_threads,
        ]
        _push_optional_mosek_float_attribute!(
            attrs,
            "NCTS_SOLVER_PROBE_MOSEK_TOL_PFEAS",
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS",
        )
        _push_optional_mosek_float_attribute!(
            attrs,
            "NCTS_SOLVER_PROBE_MOSEK_TOL_DFEAS",
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS",
        )
        _push_optional_mosek_float_attribute!(
            attrs,
            "NCTS_SOLVER_PROBE_MOSEK_TOL_REL_GAP",
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP",
        )
        return optimizer_with_attributes(Mosek.Optimizer, attrs...)
    elseif name == :cosmo
        @eval import COSMO
        return optimizer_with_attributes(COSMO.Optimizer, "verbose" => false)
    end
    throw(ArgumentError("Unsupported NCTS_SOLVER_PROBE_OPTIMIZER=$(repr(name))."))
end

function _block_histogram(sizes)
    counts = Dict{Int,Int}()
    for size in sizes
        counts[Int(size)] = get(counts, Int(size), 0) + 1
    end
    return sort!(collect(counts); by=first)
end

function _max_native_hermitian_residual(psd_block_diagnostics)
    isempty(psd_block_diagnostics) && return 0.0
    return maximum(diag -> Float64(diag.native_hermitian_residual), psd_block_diagnostics)
end

function _min_native_eigenvalue(psd_block_diagnostics)
    isempty(psd_block_diagnostics) && return nothing
    return minimum(diag -> Float64(diag.native_min_eigenvalue), psd_block_diagnostics)
end

function _max_zero_dual_abs_value(zero_dual_diagnostics)
    isempty(zero_dual_diagnostics) && return 0.0
    return maximum(diag -> Float64(diag.abs_value), zero_dual_diagnostics)
end

function _constraint_type_counts(model::JuMP.Model)
    rows = String[]
    for (F, S) in JuMP.list_of_constraint_types(model)
        n = JuMP.num_constraints(model, F, S)
        n == 0 && continue
        push!(rows, "$(F) in $(S) => $n")
    end
    return rows
end

function _scalar_equality_count(model::JuMP.Model)
    return sum(
        JuMP.num_constraints(model, F, S)
        for (F, S) in JuMP.list_of_constraint_types(model)
        if S <: JuMP.MOI.EqualTo;
        init=0,
    )
end

function _package_version(name::AbstractString)
    for dependency in values(Pkg.dependencies())
        dependency.name == name || continue
        return isnothing(dependency.version) ? nothing : string(dependency.version)
    end
    return nothing
end

function _mosek_library_version(mosek_module)
    major, minor, revision = mosek_module.getversion()
    return string(VersionNumber(Int(major), Int(minor), Int(revision)))
end

function _mosek_library_version()
    isdefined(@__MODULE__, :Mosek) || return _package_version("Mosek")
    try
        return _mosek_library_version(getfield(@__MODULE__, :Mosek))
    catch
        return _package_version("Mosek")
    end
end

function _active_project_name()
    project = Base.active_project()
    project === nothing && return "unknown"
    return basename(dirname(project))
end

function _repo_commit()
    env_commit = strip(get(ENV, "NCTS_NCTSSOS_COMMIT", ""))
    isempty(env_commit) || return env_commit

    repo = normpath(joinpath(@__DIR__, ".."))
    try
        return readchomp(pipeline(`git -C $repo rev-parse HEAD`; stderr=devnull))
    catch
        return nothing
    end
end

function _optimizer_display_name(name::Symbol)
    name == :mosek && return "Mosek"
    name == :clarabel && return "Clarabel"
    name == :cosmo && return "COSMO"
    return String(name)
end

function _infer_probe_profile(
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    linear_state_opt_width,
    psd_state_opt_width,
)
    has_rdm = !isnothing(contiguous_rdm_k) && contiguous_rdm_k > 0
    has_lso = !isnothing(linear_state_opt_width) && linear_state_opt_width > 0
    has_pso = !isnothing(psd_state_opt_width) && psd_state_opt_width > 0
    !has_rdm && !has_lso && !has_pso && return "A0"
    if has_rdm && contiguous_rdm_decomposition == :qmbcertify && has_lso
        has_pso && return "A2"
        return "A1"
    end
    return "custom"
end

function _probe_profile(
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    linear_state_opt_width,
    psd_state_opt_width,
)
    raw = strip(get(ENV, "NCTS_SOLVER_PROBE_PROFILE", ""))
    isempty(raw) || return raw
    return _infer_probe_profile(
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        linear_state_opt_width,
        psd_state_opt_width,
    )
end

function _toml_escape(value::AbstractString)
    escaped = replace(value, "\\" => "\\\\", "\"" => "\\\"")
    return "\"" * escaped * "\""
end

_toml_value(value::AbstractString) = _toml_escape(value)
_toml_value(value::Symbol) = _toml_escape(String(value))
_toml_value(value::Bool) = string(value)
_toml_value(value::Integer) = string(value)
_toml_value(value::AbstractFloat) = string(value)

function _toml_value(values::AbstractVector)
    return "[" * join((_toml_value(value) for value in values), ", ") * "]"
end

_abs_float(value) = value === nothing ? nothing : Float64(abs(value))

function _toml_table_key(value::AbstractString)
    occursin(r"^[A-Za-z0-9_-]+$", value) && return value
    return _toml_escape(value)
end

function _print_toml_field(name::AbstractString, value)
    value === nothing && return nothing
    println("$name = $(_toml_value(value))")
    return nothing
end

function _probe_fixture_id(
    profile::AbstractString,
    n::Integer,
    model_mode::AbstractString,
    optimizer_name::Symbol,
)
    raw_date = strip(get(ENV, "NCTS_SOLVER_PROBE_ATTEMPT_DATE", ""))
    date = isempty(raw_date) ?
        Dates.format(Dates.now(Dates.UTC), "yyyy_mm_dd") :
        raw_date
    optimizer = String(optimizer_name)
    return "NCTSSOS_$(profile)_L$(n)_source_like_$(model_mode)_$(optimizer)_$(date)"
end

function _emit_solver_probe_environment_fixture(environment::AbstractString)
    println()
    println("# nctssos_solver_probe_environment_fixture_begin")
    println("[nctssos_solver_probe_environments.$(_toml_table_key(environment))]")
    _print_toml_field("nctssos_commit", _repo_commit())
    _print_toml_field("julia_version", string(VERSION))
    _print_toml_field("mosek_version", _mosek_library_version())
    _print_toml_field(
        "mosek_tol_pfeas",
        _parse_optional_float_env("NCTS_SOLVER_PROBE_MOSEK_TOL_PFEAS"),
    )
    _print_toml_field(
        "mosek_tol_dfeas",
        _parse_optional_float_env("NCTS_SOLVER_PROBE_MOSEK_TOL_DFEAS"),
    )
    _print_toml_field(
        "mosek_tol_rel_gap",
        _parse_optional_float_env("NCTS_SOLVER_PROBE_MOSEK_TOL_REL_GAP"),
    )
    _print_toml_field("cpu_model", Sys.cpu_info()[1].model)
    _print_toml_field("ram_bytes", _mem_total_bytes())
    _print_toml_field("thread_count", Threads.nthreads())
    _print_toml_field("blas_vendor", string(BLAS.vendor()))
    _print_toml_field("project", _active_project_name())
    _print_toml_field("probe_harness", _SOLVER_PROBE_HARNESS)
    println("# nctssos_solver_probe_environment_fixture_end")
    return nothing
end

function _probe_execution_state(lower_stats, dualize_stats, solve_result)
    model_built = lower_stats !== nothing || dualize_stats !== nothing
    solved = solve_result !== nothing
    return (construction_only=!model_built, model_built=model_built, solved=solved)
end

function _solver_probe_qmbcertify_reference()
    probe_id = _parse_optional_string_env("NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID")
    objective_total_estimate = _parse_optional_float_env(
        "NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE",
    )
    isnothing(probe_id) == isnothing(objective_total_estimate) || throw(ArgumentError(
        "NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID and " *
        "NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE must be provided together."
    ))
    return (
        probe_id=probe_id,
        objective_total_estimate=objective_total_estimate,
    )
end

function _check_source_base_solved_probe_qmbcertify_reference(
    qmbcertify_base_construct::Bool,
    solve_result,
    qmbcertify_reference,
)
    qmbcertify_base_construct || return nothing
    solve_result === nothing && return nothing
    !isnothing(qmbcertify_reference.probe_id) &&
        !isnothing(qmbcertify_reference.objective_total_estimate) &&
        return nothing
    throw(ArgumentError(
        "Solved source-base QMBCertify solver probes must set " *
        "NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID and " *
        "NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE before " *
        "emitting a solved fixture row."
    ))
end

_triangle_count(n::Integer) = Int(n) * (Int(n) + 1) ÷ 2

function _solver_probe_dense_schur_bytes(n_equalities::Integer)
    n = Int(n_equalities)
    n >= 0 || throw(ArgumentError("n_equalities must be nonnegative, got $n."))
    return 8 * n^2
end

function _solver_probe_n_aux_blocks(n_orphans::Integer; orphans_per_block::Integer=32)
    n = Int(n_orphans)
    n == 0 && return 0
    return cld(n, Int(orphans_per_block))
end

function _normalize_solver_probe_hermitian_representation(representation::Symbol)
    representation in (:real_lift, :lift, :lifted, :hermitian_lift) &&
        return :hermitian_lift
    representation in (:native, :hermitian) && return :hermitian
    throw(ArgumentError(
        "Unsupported hermitian_representation $(repr(representation)); expected :real_lift or :native."
    ))
end

function _normalize_solver_probe_coefficient_scaling(scaling::Symbol)
    scaling in (:none, :max_abs) && return scaling
    throw(ArgumentError(
        "Unsupported SOS coefficient scaling $(repr(scaling)); expected :none or :max_abs."
    ))
end

function _normalize_solver_probe_coefficient_scaling_floor(scaling_floor::Real)
    floor = Float64(scaling_floor)
    (isfinite(floor) && floor >= 0) || throw(ArgumentError(
        "SOS coefficient scaling floor must be finite and nonnegative, got $scaling_floor."
    ))
    return floor
end

function _psd_cone_scalar_variables(size::Integer, cone::Symbol)
    n = Int(size)
    cone == :PSD && return _triangle_count(n)
    cone in (:HPSD, :AuxHPSD) && return n^2
    throw(ArgumentError("Unsupported PSD cone $(repr(cone))."))
end

function _dual_psd_cone_scalar_variables(
    size::Integer,
    cone::Symbol,
    hermitian_representation::Symbol,
)
    n = Int(size)
    cone == :PSD && return _triangle_count(n)
    cone == :HPSD || throw(ArgumentError(
        "SOS dual model-size estimate expected :PSD or :HPSD cone, got $(repr(cone))."
    ))
    representation =
        _normalize_solver_probe_hermitian_representation(hermitian_representation)
    representation == :hermitian_lift && return _triangle_count(2 * n)
    return n^2
end

function _primal_psd_block_binding_upper_bound(size::Integer, cone::Symbol)
    n = Int(size)
    cone == :PSD && return _triangle_count(n)
    cone in (:HPSD, :AuxHPSD) && return 2 * n^2
    throw(ArgumentError("Unsupported PSD cone $(repr(cone))."))
end

function _solver_probe_model_size_estimate(;
    block_sizes,
    block_cones,
    n_moments::Integer,
    n_free_keys::Integer,
    n_zero_constraints::Integer,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
    dualize::Bool,
    sos_hermitian_representation::Symbol,
)
    length(block_sizes) == length(block_cones) || throw(ArgumentError(
        "block_sizes and block_cones must have the same length."
    ))
    n_moments = Int(n_moments)
    n_free_keys = Int(n_free_keys)
    n_zero_constraints = Int(n_zero_constraints)
    model_mode = dualize ? "sos_dual" : "primal"

    if dualize
        psd_variables = sum(
            _dual_psd_cone_scalar_variables(size, cone, sos_hermitian_representation)
            for (size, cone) in zip(block_sizes, block_cones);
            init=0,
        )
        has_hermitian = any(cone -> cone == :HPSD, block_cones)
        scalar_equalities = has_hermitian ?
            max(2 * n_moments - 1, 0) :
            max(n_moments - 1, 0)
        return (
            model_mode=model_mode,
            model_variables=psd_variables + n_zero_constraints,
            psd_cone_scalar_variables=psd_variables,
            scalar_equalities_upper_bound=scalar_equalities,
            dense_schur_bytes=_solver_probe_dense_schur_bytes(scalar_equalities),
            zero_dual_variables=n_zero_constraints,
            free_orphan_variables=0,
            aux_orphan_blocks=0,
            lowering_would_error=false,
        )
    end

    formulation in (:moment_variables, :psd_blocks) || throw(ArgumentError(
        "Unsupported formulation $(repr(formulation)); expected :moment_variables or :psd_blocks."
    ))
    if formulation == :moment_variables
        has_hermitian = any(cone -> cone in (:HPSD, :AuxHPSD), block_cones)
        moment_variables =
            (representation == :complex || has_hermitian) ? 2 * n_moments : n_moments
        psd_cone_scalars = sum(
            cone == :PSD ? _triangle_count(size) : _triangle_count(2 * Int(size))
            for (size, cone) in zip(block_sizes, block_cones);
            init=0,
        )
        identity_equalities = (representation == :complex || has_hermitian) ? 2 : 1
        scalar_equalities = identity_equalities + n_zero_constraints
        return (
            model_mode=model_mode,
            model_variables=moment_variables,
            psd_cone_scalar_variables=psd_cone_scalars,
            scalar_equalities_upper_bound=scalar_equalities,
            dense_schur_bytes=_solver_probe_dense_schur_bytes(scalar_equalities),
            zero_dual_variables=0,
            free_orphan_variables=0,
            aux_orphan_blocks=0,
            lowering_would_error=false,
        )
    end

    psd_variables = sum(
        _psd_cone_scalar_variables(size, cone)
        for (size, cone) in zip(block_sizes, block_cones);
        init=0,
    )
    aux_orphan_blocks = 0
    free_orphan_variables = 0
    aux_orphan_variables = 0
    lowering_would_error = n_free_keys > 0 && orphan_policy == :error
    if n_free_keys > 0 && orphan_policy == :free_variables
        free_orphan_variables =
            representation == :complex ? 2 * n_free_keys : n_free_keys
    elseif n_free_keys > 0 && orphan_policy == :aux_psd_free
        aux_orphan_blocks = _solver_probe_n_aux_blocks(n_free_keys)
        for block_idx in 1:aux_orphan_blocks
            n_in_block = min(32, n_free_keys - 32 * (block_idx - 1))
            aux_size = n_in_block + 1
            aux_orphan_variables += representation == :complex ?
                aux_size^2 :
                _triangle_count(aux_size)
        end
    end
    binding_equalities = sum(
        _primal_psd_block_binding_upper_bound(size, cone)
        for (size, cone) in zip(block_sizes, block_cones);
        init=0,
    )
    identity_equalities = representation == :complex ? 2 : 1
    scalar_equalities = binding_equalities + identity_equalities + n_zero_constraints
    return (
        model_mode=model_mode,
        model_variables=psd_variables + free_orphan_variables + aux_orphan_variables,
        psd_cone_scalar_variables=psd_variables + aux_orphan_variables,
        scalar_equalities_upper_bound=scalar_equalities,
        dense_schur_bytes=_solver_probe_dense_schur_bytes(scalar_equalities),
        zero_dual_variables=0,
        free_orphan_variables=free_orphan_variables,
        aux_orphan_blocks=aux_orphan_blocks,
        lowering_would_error=lowering_would_error,
    )
end

function _solver_probe_model_size_estimate(
    linear;
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
    dualize::Bool,
    sos_hermitian_representation::Symbol,
)
    return _solver_probe_model_size_estimate(
        block_sizes=getfield.(linear.psd_blocks_lin, :size),
        block_cones=[block.meta.cone for block in linear.psd_blocks_lin],
        n_moments=length(linear.moments),
        n_free_keys=length(linear.free_keys),
        n_zero_constraints=length(linear.zero_constraints),
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
        dualize=dualize,
        sos_hermitian_representation=sos_hermitian_representation,
    )
end

function _solver_probe_risk_gate_status(estimate)
    estimate.lowering_would_error && return "blocked_lowering_orphan_policy"
    return "ok"
end

function _solver_probe_model_size_gate_reason(;
    status::AbstractString,
    mem_available_bytes,
    estimated_rss_bytes,
    max_rss_fraction::Real,
)
    status == "ok" && return ""
    status == "blocked_missing_memory_telemetry" &&
        return "available-memory telemetry is unavailable"
    status == "blocked_insufficient_memory" && return (
        "estimated dense-Schur RSS proxy " *
        "$(_ncts_load_guard_gib_label(estimated_rss_bytes)) exceeds " *
        "available memory $(_ncts_load_guard_gib_label(mem_available_bytes))"
    )
    status == "blocked_insufficient_memory_headroom" && return (
        "estimated dense-Schur RSS proxy " *
        "$(_ncts_load_guard_gib_label(estimated_rss_bytes)) would consume " *
        "more than $(Float64(max_rss_fraction) * 100)% of available memory " *
        "$(_ncts_load_guard_gib_label(mem_available_bytes))"
    )
    return "model-size gate returned $status"
end

function _solver_probe_model_size_gate(
    estimate;
    mem_available_bytes=_ncts_mem_available_bytes(),
    max_rss_fraction::Real=_ncts_load_guard_parse_float(
        "NCTS_SOLVER_PROBE_MAX_RSS_FRACTION",
        _ncts_load_guard_parse_float("NCTS_LOAD_GUARD_MAX_RSS_FRACTION", 0.8),
    ),
)
    estimated_rss_bytes = estimate.dense_schur_bytes
    status = _ncts_memory_pressure_status(
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
        max_rss_fraction=max_rss_fraction,
    )
    reason = _solver_probe_model_size_gate_reason(
        status=status,
        mem_available_bytes=mem_available_bytes,
        estimated_rss_bytes=estimated_rss_bytes,
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

function _solver_probe_model_size_gate_status(estimate; kwargs...)
    return _solver_probe_model_size_gate(estimate; kwargs...).status
end

function _solver_probe_enforce_model_size_gate()
    default = _parse_bool_env("NCTS_SOLVER_PROBE_ALLOW_LARGE", false)
    return _parse_bool_env("NCTS_SOLVER_PROBE_ENFORCE_MODEL_SIZE_GATE", default)
end

function _print_model_size_estimate(estimate; model_size_gate)
    println("estimated_model_mode=$(estimate.model_mode)")
    println("estimated_model_variables=$(estimate.model_variables)")
    println("estimated_psd_cone_scalar_variables=$(estimate.psd_cone_scalar_variables)")
    println(
        "estimated_scalar_equalities_upper_bound=$(estimate.scalar_equalities_upper_bound)",
    )
    println("estimated_dense_schur_bytes=$(estimate.dense_schur_bytes)")
    println("estimated_model_size_gate_status=$(model_size_gate.status)")
    println("estimated_model_size_gate_reason=$(model_size_gate.reason)")
    println(
        "estimated_model_size_gate_estimated_rss_bytes=$(model_size_gate.estimated_rss_bytes)",
    )
    println(
        "estimated_model_size_gate_mem_available_bytes=$(isnothing(model_size_gate.mem_available_bytes) ? "unknown" : model_size_gate.mem_available_bytes)",
    )
    println(
        "estimated_model_size_gate_max_rss_fraction=$(model_size_gate.max_rss_fraction)",
    )
    println("estimated_zero_dual_variables=$(estimate.zero_dual_variables)")
    println("estimated_free_orphan_variables=$(estimate.free_orphan_variables)")
    println("estimated_aux_orphan_blocks=$(estimate.aux_orphan_blocks)")
    println("estimated_lowering_would_error=$(estimate.lowering_would_error)")
    println("estimated_risk_gate_status=$(_solver_probe_risk_gate_status(estimate))")
    return nothing
end

function _construction_stage_timing_fields(stage_times)
    fields = Pair{String,Float64}[]
    for (stage, seconds) in stage_times
        push!(fields, "construction_stage_$(stage)_seconds" => Float64(seconds))
    end
    return sort!(fields; by=first)
end

function _emit_solver_probe_fixture_row(;
    n::Integer,
    order::Integer,
    model_mode::AbstractString,
    formulation::Symbol,
    representation::Symbol,
    optimizer_name::Symbol,
    solve_support,
    reflection_symmetry::Bool,
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool,
    axis_rotation_symmetry::Bool,
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    linear_state_opt_width,
    linear_state_opt_mode::Symbol,
    psd_state_opt_width,
    qmbcertify_base_construct::Bool,
    qmbcertify_base_extra,
    qmbcertify_base_three_type::Tuple{Int,Int},
    model_size_estimate,
    model_size_gate,
    construct_stats,
    lower_stats,
    dualize_stats,
    solve_result,
    solve_stats,
    sos_certificate_residual,
    sos_certificate_diagnostics,
    sos_hermitian_representation::Symbol,
    sos_coefficient_scaling::Symbol,
    sos_coefficient_scaling_floor::Float64,
    script_wall_seconds::Float64,
    report_metrics,
    linear,
    dual_variables,
    dual_zero_duals,
    model_variables,
    scalar_equalities,
)
    profile = _probe_profile(
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        linear_state_opt_width,
        psd_state_opt_width,
    )
    environment = strip(get(ENV, "NCTS_SOLVER_PROBE_ENVIRONMENT", "manual"))
    remote_wall_seconds = _parse_optional_float_env("NCTS_SOLVER_PROBE_REMOTE_WALL_SECONDS")
    peak_rss_kb = _proc_status_kb("VmHWM")
    rdm_k = isnothing(contiguous_rdm_k) ? 0 : Int(contiguous_rdm_k)
    rdm_decomposition = iszero(rdm_k) ? "none" : String(contiguous_rdm_decomposition)
    lso_width = isnothing(linear_state_opt_width) ? 0 : Int(linear_state_opt_width)
    pso_width = isnothing(psd_state_opt_width) ? 0 : Int(psd_state_opt_width)
    execution_state = _probe_execution_state(lower_stats, dualize_stats, solve_result)
    qmbcertify_reference = _solver_probe_qmbcertify_reference()
    _check_source_base_solved_probe_qmbcertify_reference(
        qmbcertify_base_construct,
        solve_result,
        qmbcertify_reference,
    )
    objective_minus_qmbcertify_total_estimate =
        solve_result === nothing || isnothing(qmbcertify_reference.objective_total_estimate) ?
        nothing :
        solve_result.objective - qmbcertify_reference.objective_total_estimate

    println()
    println("# nctssos_source_like_solver_probe_fixture_begin")
    println("[[nctssos_source_like_solver_probes]]")
    _print_toml_field("id", _probe_fixture_id(profile, n, model_mode, optimizer_name))
    _print_toml_field("profile", profile)
    _print_toml_field("length", Int(n))
    _print_toml_field("status", "nctssos_solver_probe")
    _print_toml_field("comparison_scope", _solver_probe_comparison_scope(n))
    _print_toml_field("environment", environment)
    _print_toml_field("probe_harness", _SOLVER_PROBE_HARNESS)
    _print_toml_field("order", Int(order))
    _print_toml_field("model_mode", model_mode)
    _print_toml_field("formulation", String(formulation))
    _print_toml_field("representation", String(representation))
    _print_toml_field("sos_hermitian_representation", String(sos_hermitian_representation))
    _print_toml_field("sos_coefficient_scaling", String(sos_coefficient_scaling))
    _print_toml_field("sos_coefficient_scaling_floor", sos_coefficient_scaling_floor)
    _print_toml_field("optimizer", _optimizer_display_name(optimizer_name))
    _print_toml_field("reflection_symmetry", reflection_symmetry)
    _print_toml_field("sign_symmetry", sign_symmetry)
    _print_toml_field("real_moment_matrix", real_moment_matrix)
    _print_toml_field("u1_symmetry", u1_symmetry)
    _print_toml_field("su2_symmetry", su2_symmetry)
    _print_toml_field("base_su2_extend_rdm", base_su2_extend_rdm)
    _print_toml_field("su2_moment_quotient", su2_moment_quotient)
    _print_toml_field("su2_moment_raw_count", report_metrics.su2_moment_raw_count)
    _print_toml_field(
        "su2_moment_quotient_count",
        report_metrics.su2_moment_quotient_count,
    )
    _print_toml_field(
        "su2_moment_quotient_reduction_ratio",
        report_metrics.su2_moment_quotient_reduction_ratio,
    )
    _print_toml_field(
        "su2_moment_support_orbit_count",
        report_metrics.su2_moment_support_orbit_count,
    )
    _print_toml_field(
        "su2_moment_singlet_channel_support_counts",
        [Int[first(pair), last(pair)] for pair in
            report_metrics.su2_moment_singlet_channel_support_counts],
    )
    _print_toml_field(
        "su2_moment_max_pivot_residual",
        report_metrics.su2_moment_max_pivot_residual,
    )
    _print_toml_field(
        "su2_moment_max_invariant_residual",
        report_metrics.su2_moment_max_invariant_residual,
    )
    _print_toml_field(
        "su2_moment_max_reconstruction_residual",
        report_metrics.su2_moment_max_reconstruction_residual,
    )
    _print_toml_field(
        "su2_moment_max_condition",
        report_metrics.su2_moment_max_condition,
    )
    _print_toml_field(
        "su2_moment_eliminated_zero_row_count",
        report_metrics.su2_moment_eliminated_zero_row_count,
    )
    _print_toml_field("axis_rotation_symmetry", axis_rotation_symmetry)
    _print_toml_field("axis_rotation_equalities", axis_rotation_equalities)
    _print_toml_field("axis_rotation_quotient", axis_rotation_quotient)
    _print_toml_field("qmbcertify_base_construct", qmbcertify_base_construct)
    _print_toml_field(
        "qmbcertify_base_extra",
        isnothing(qmbcertify_base_extra) ? 0 : Int(qmbcertify_base_extra),
    )
    _print_toml_field("qmbcertify_base_three_type", collect(qmbcertify_base_three_type))
    _print_toml_field("report_su2_moment_symmetry", report_metrics.su2_moment_symmetry)
    _print_toml_field("report_su2_rdm_symmetry", report_metrics.su2_rdm_symmetry)
    _print_toml_field("report_u1_rdm_symmetry", report_metrics.u1_rdm_symmetry)
    _print_toml_field(
        "report_axis_rotation_symmetry",
        report_metrics.axis_rotation_symmetry,
    )
    _print_toml_field("report_axis_quotient", report_metrics.axis_rotation_quotient)
    _print_toml_field("su2_base_zero_row_count", report_metrics.su2_base_zero_row_count)
    _print_toml_field(
        "su2_base_spin_offblock_row_count",
        report_metrics.su2_base_spin_offblock_row_count,
    )
    _print_toml_field(
        "su2_base_magnetic_offdiag_row_count",
        report_metrics.su2_base_magnetic_offdiag_row_count,
    )
    _print_toml_field(
        "su2_base_magnetic_copy_row_count",
        report_metrics.su2_base_magnetic_copy_row_count,
    )
    _print_toml_field("solve_supported", solve_support.supported)
    _print_toml_field(
        "solve_blocker",
        solve_support.blocker === nothing ? "none" : String(solve_support.blocker),
    )
    _print_toml_field("solve_support_reason", solve_support.reason)
    _print_toml_field(
        "solve_unsupported_block_features",
        Symbol.(solve_support.unsupported_block_features),
    )
    _print_toml_field(
        "solve_unsupported_zero_features",
        Symbol.(solve_support.unsupported_zero_features),
    )
    _print_toml_field("estimated_model_mode", model_size_estimate.model_mode)
    _print_toml_field("estimated_model_variables", model_size_estimate.model_variables)
    _print_toml_field(
        "estimated_psd_cone_scalar_variables",
        model_size_estimate.psd_cone_scalar_variables,
    )
    _print_toml_field(
        "estimated_scalar_equalities_upper_bound",
        model_size_estimate.scalar_equalities_upper_bound,
    )
    _print_toml_field(
        "estimated_dense_schur_bytes",
        model_size_estimate.dense_schur_bytes,
    )
    _print_toml_field("estimated_model_size_gate_status", model_size_gate.status)
    _print_toml_field("estimated_model_size_gate_reason", model_size_gate.reason)
    _print_toml_field(
        "estimated_model_size_gate_estimated_rss_bytes",
        model_size_gate.estimated_rss_bytes,
    )
    _print_toml_field(
        "estimated_model_size_gate_mem_available_bytes",
        model_size_gate.mem_available_bytes,
    )
    _print_toml_field(
        "estimated_model_size_gate_max_rss_fraction",
        model_size_gate.max_rss_fraction,
    )
    _print_toml_field(
        "estimated_zero_dual_variables",
        model_size_estimate.zero_dual_variables,
    )
    _print_toml_field(
        "estimated_free_orphan_variables",
        model_size_estimate.free_orphan_variables,
    )
    _print_toml_field(
        "estimated_aux_orphan_blocks",
        model_size_estimate.aux_orphan_blocks,
    )
    _print_toml_field(
        "estimated_lowering_would_error",
        model_size_estimate.lowering_would_error,
    )
    _print_toml_field(
        "estimated_risk_gate_status",
        _solver_probe_risk_gate_status(model_size_estimate),
    )
    if solve_result !== nothing
        _print_toml_field("termination_status", string(solve_result.status))
        _print_toml_field("primal_status", string(solve_result.primal_status))
        _print_toml_field("dual_status", string(solve_result.dual_status))
        _print_toml_field("raw_status", solve_result.raw_status)
        _print_toml_field("result_count", solve_result.result_count)
        _print_toml_field("objective", solve_result.objective)
        _print_toml_field("objective_bound", solve_result.objective_bound)
        _print_toml_field("relative_gap", solve_result.relative_gap)
    else
        _print_toml_field("termination_status", "not_solved")
    end
    _print_toml_field("qmbcertify_probe_id", qmbcertify_reference.probe_id)
    _print_toml_field(
        "qmbcertify_objective_total_estimate",
        qmbcertify_reference.objective_total_estimate,
    )
    _print_toml_field(
        "objective_minus_qmbcertify_total_estimate",
        objective_minus_qmbcertify_total_estimate,
    )
    _print_toml_field("construction_only", execution_state.construction_only)
    _print_toml_field("model_built", execution_state.model_built)
    _print_toml_field("solved", execution_state.solved)
    if sos_certificate_residual !== nothing
        _print_toml_field(
            "sos_certificate_moment_count",
            sos_certificate_residual.moment_count,
        )
        _print_toml_field(
            "sos_certificate_residual_source",
            "sos_dual_certificate_residual",
        )
        _print_toml_field(
            "sos_certificate_max_abs_residual",
            Float64(sos_certificate_residual.max_abs_residual),
        )
        _print_toml_field(
            "sos_certificate_identity_abs_residual",
            _abs_float(sos_certificate_residual.identity_residual),
        )
        _print_toml_field(
            "sos_certificate_worst_abs_residual",
            _abs_float(sos_certificate_residual.max_residual_value),
        )
    end
    if sos_certificate_diagnostics !== nothing
        _print_toml_field(
            "sos_certificate_diagnostics_source",
            "sos_dual_certificate_diagnostics",
        )
        _print_toml_field(
            "sos_certificate_psd_block_count",
            sos_certificate_diagnostics.psd_block_count,
        )
        _print_toml_field(
            "sos_certificate_zero_dual_count",
            sos_certificate_diagnostics.zero_dual_count,
        )
        _print_toml_field(
            "sos_certificate_max_native_hermitian_residual",
            _max_native_hermitian_residual(sos_certificate_diagnostics.psd_blocks),
        )
        _print_toml_field(
            "sos_certificate_min_native_eigenvalue",
            _min_native_eigenvalue(sos_certificate_diagnostics.psd_blocks),
        )
        _print_toml_field(
            "sos_certificate_max_zero_dual_abs_value",
            _max_zero_dual_abs_value(sos_certificate_diagnostics.zero_duals),
        )
    end
    _print_toml_field("construction_time_seconds", construct_stats.time)
    for (field, seconds) in _construction_stage_timing_fields(
        report_metrics.construction_stage_time_seconds,
    )
        _print_toml_field(field, seconds)
    end
    _print_toml_field("lowering_time_seconds", lower_stats === nothing ? nothing : lower_stats.time)
    _print_toml_field("dualization_time_seconds", dualize_stats === nothing ? nothing : dualize_stats.time)
    _print_toml_field("solve_time_seconds", solve_stats === nothing ? nothing : solve_stats.time)
    _print_toml_field("script_wall_seconds", script_wall_seconds)
    _print_toml_field("remote_wall_seconds", remote_wall_seconds)
    _print_toml_field("peak_rss_bytes", peak_rss_kb === nothing ? nothing : peak_rss_kb * 1024)
    _print_toml_field("linear_moments", length(linear.moments))
    _print_toml_field("free_keys", length(linear.free_keys))
    _print_toml_field("linear_state_opt_width", lso_width)
    _print_toml_field("linear_state_opt_mode", String(linear_state_opt_mode))
    _print_toml_field("linear_state_opt_row_count", report_metrics.linear_state_opt_row_count)
    _print_toml_field("contiguous_rdm_k", rdm_k)
    _print_toml_field("contiguous_rdm_decomposition", rdm_decomposition)
    _print_toml_field("contiguous_rdm_support", String(contiguous_rdm_support))
    _print_toml_field(
        "qmbcertify_rdm_block_sizes",
        Int.(report_metrics.contiguous_rdm_psd_block_sizes),
    )
    _print_toml_field("psd_state_opt_width", pso_width)
    _print_toml_field(
        "psd_state_opt_block_sizes",
        Int.(report_metrics.psd_state_opt_psd_block_sizes),
    )
    _print_toml_field("psd_cones", length(linear.psd_blocks_lin))
    _print_toml_field("dual_variables", dual_variables)
    _print_toml_field("zero_dual_variables", dual_zero_duals)
    _print_toml_field("model_variables", model_variables)
    _print_toml_field("scalar_equalities", scalar_equalities)
    _print_toml_field(
        "dense_schur_bytes",
        scalar_equalities === nothing ? nothing :
            _solver_probe_dense_schur_bytes(scalar_equalities),
    )
    _print_toml_field("max_block", report_metrics.psd_max_block)
    _print_toml_field("product_cache_hit_rate", report_metrics.product_cache_hit_rate)
    println("# nctssos_source_like_solver_probe_fixture_end")
    _emit_solver_probe_environment_fixture(environment)
    return nothing
end

function _heisenberg_case(n::Integer)
    registry, ops = create_pauli_variables(1:Int(n))
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    return registry, ops, pop
end

function _solver_probe_linear_relaxation(
    pop,
    ops,
    order::Integer;
    qmbcertify_base_construct::Bool,
    qmbcertify_base_extra,
    qmbcertify_base_three_type::Tuple{Int,Int},
    reflection_symmetry::Bool,
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    axis_rotation_symmetry::Bool,
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    contiguous_rdm_k,
    contiguous_rdm_decomposition::Symbol,
    contiguous_rdm_support::Symbol,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool,
    linear_state_opt_width,
    linear_state_opt_mode::Symbol,
    psd_state_opt_width,
)
    if qmbcertify_base_construct
        return NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
            pop,
            ops,
            order;
            extra=something(qmbcertify_base_extra, 0),
            three_type=qmbcertify_base_three_type,
            real_moment_matrix,
            contiguous_rdm_k,
            contiguous_rdm_decomposition,
            contiguous_rdm_support,
            linear_state_opt_width,
            linear_state_opt_mode,
            psd_state_opt_width,
        )
    end
    return NCTSSoS._pauli_translation_base_linear_relaxation(
        pop,
        ops,
        order;
        reflection_symmetry,
        sign_symmetry,
        real_moment_matrix,
        axis_rotation_symmetry,
        axis_rotation_equalities,
        axis_rotation_quotient,
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        contiguous_rdm_support,
        u1_symmetry,
        su2_symmetry,
        base_su2_extend_rdm,
        su2_moment_quotient,
        linear_state_opt_width,
        linear_state_opt_mode,
        psd_state_opt_width,
    )
end

function main()
    probe_start_ns = time_ns()
    n = parse(Int, get(ENV, "NCTS_SOLVER_PROBE_N", "6"))
    _check_size_guard(n)

    order = parse(Int, get(ENV, "NCTS_SOLVER_PROBE_ORDER", "4"))
    formulation = _parse_symbol_env("NCTS_SOLVER_PROBE_FORMULATION", :psd_blocks)
    orphan_policy = _parse_symbol_env("NCTS_SOLVER_PROBE_ORPHAN_POLICY", :error)
    optimizer_name = _parse_symbol_env("NCTS_SOLVER_PROBE_OPTIMIZER", :clarabel)
    dualize = _parse_bool_env("NCTS_SOLVER_PROBE_DUALIZE", false)
    sos_hermitian_representation = _parse_symbol_env("NCTS_SOLVER_PROBE_SOS_HERMITIAN", :real_lift)
    sos_coefficient_scaling = _normalize_solver_probe_coefficient_scaling(
        _parse_symbol_env("NCTS_SOLVER_PROBE_SOS_EQUATION_SCALING", :none),
    )
    sos_coefficient_scaling_floor =
        _normalize_solver_probe_coefficient_scaling_floor(parse(
            Float64,
            get(ENV, "NCTS_SOLVER_PROBE_SOS_EQUATION_SCALING_FLOOR", "0"),
        ))
    solve = _parse_bool_env("NCTS_SOLVER_PROBE_SOLVE", false)
    lower_model = _parse_lower_model_env()
    emit_fixture = _parse_bool_env("NCTS_SOLVER_PROBE_EMIT_FIXTURE", false)
    solve && !lower_model && throw(ArgumentError(
        "NCTS_SOLVER_PROBE_SOLVE=true requires NCTS_SOLVER_PROBE_LOWER_MODEL=true."
    ))

    qmbcertify_base_construct =
        _parse_bool_env("NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT", false)
    qmbcertify_base_extra =
        _parse_optional_int_env("NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_EXTRA")
    qmbcertify_base_three_type =
        _parse_int_tuple2_env("NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_THREE_TYPE", "1,1")
    reflection_symmetry = _parse_bool_env("NCTS_SOLVER_PROBE_REFLECTION", false)
    sign_symmetry = _parse_bool_env("NCTS_SOLVER_PROBE_SIGN", true)
    real_moment_matrix = _parse_bool_env("NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX", false)
    representation = haskey(ENV, "NCTS_SOLVER_PROBE_REPRESENTATION") ?
        _parse_symbol_env("NCTS_SOLVER_PROBE_REPRESENTATION", :complex) :
        (real_moment_matrix ? :real : :complex)
    u1_symmetry = _parse_bool_env("NCTS_SOLVER_PROBE_U1", false)
    su2_symmetry = _parse_bool_env("NCTS_SOLVER_PROBE_SU2", false)
    base_su2_extend_rdm = _parse_bool_env("NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM", false)
    probe_su2_moment_quotient = _parse_bool_env(
        "NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT",
        base_su2_extend_rdm,
    )
    axis_rotation_symmetry = _parse_bool_env("NCTS_SOLVER_PROBE_AXIS_SYMMETRY", false)
    axis_rotation_equalities = _parse_bool_env("NCTS_SOLVER_PROBE_AXIS_EQUALITIES", false)
    axis_rotation_quotient = _parse_bool_env("NCTS_SOLVER_PROBE_AXIS_QUOTIENT", false)
    contiguous_rdm_k = _parse_optional_int_env("NCTS_SOLVER_PROBE_RDM_K")
    contiguous_rdm_decomposition = _parse_symbol_env(
        "NCTS_SOLVER_PROBE_RDM_DECOMPOSITION",
        qmbcertify_base_construct ? :qmbcertify : :full,
    )
    contiguous_rdm_support = _parse_symbol_env(
        "NCTS_SOLVER_PROBE_RDM_SUPPORT",
        qmbcertify_base_construct ? :extend : :closed,
    )
    psd_state_opt_width = _parse_optional_int_env("NCTS_SOLVER_PROBE_PSO_WIDTH")
    linear_state_opt_width = _parse_optional_int_env("NCTS_SOLVER_PROBE_LSO_WIDTH")
    linear_state_opt_mode = _parse_symbol_env(
        "NCTS_SOLVER_PROBE_LSO_MODE",
        qmbcertify_base_construct ? :qmbcertify : :contiguous,
    )

    _check_qmbcertify_base_solver_probe_guard(
        qmbcertify_base_construct=qmbcertify_base_construct,
        sign_symmetry=sign_symmetry,
        reflection_symmetry=reflection_symmetry,
        u1_symmetry=u1_symmetry,
        su2_symmetry=su2_symmetry,
        base_su2_extend_rdm=base_su2_extend_rdm,
        su2_moment_quotient=probe_su2_moment_quotient,
        axis_rotation_symmetry=axis_rotation_symmetry,
        axis_rotation_equalities=axis_rotation_equalities,
        axis_rotation_quotient=axis_rotation_quotient,
        contiguous_rdm_k=contiguous_rdm_k,
        contiguous_rdm_decomposition=contiguous_rdm_decomposition,
        contiguous_rdm_support=contiguous_rdm_support,
        real_moment_matrix=real_moment_matrix,
        linear_state_opt_width=linear_state_opt_width,
        linear_state_opt_mode=linear_state_opt_mode,
        psd_state_opt_width=psd_state_opt_width,
    )
    _check_base_su2_extend_rdm_solver_probe_guard(
        base_su2_extend_rdm=base_su2_extend_rdm,
        su2_moment_quotient=probe_su2_moment_quotient,
        su2_symmetry=su2_symmetry,
        contiguous_rdm_decomposition=contiguous_rdm_decomposition,
        contiguous_rdm_support=contiguous_rdm_support,
        real_moment_matrix=real_moment_matrix,
        reflection_symmetry=reflection_symmetry,
    )
    _check_deleted_primal_su2_solved_probe_guard(
        n=n,
        solve=solve,
        dualize=dualize,
        formulation=formulation,
        su2_symmetry=su2_symmetry,
        su2_moment_quotient=probe_su2_moment_quotient,
        contiguous_rdm_decomposition=contiguous_rdm_decomposition,
        linear_state_opt_width=linear_state_opt_width,
        psd_state_opt_width=psd_state_opt_width,
    )

    println("probe_n=$n")
    println("probe_order=$order")
    println("probe_formulation=$formulation")
    println("probe_representation=$representation")
    println("probe_orphan_policy=$orphan_policy")
    println("probe_optimizer=$optimizer_name")
    println("probe_dualize=$dualize")
    println("probe_sos_hermitian_representation=$sos_hermitian_representation")
    println("probe_sos_coefficient_scaling=$sos_coefficient_scaling")
    println("probe_sos_coefficient_scaling_floor=$sos_coefficient_scaling_floor")
    model_mode = dualize ? "sos_dual" : "primal"
    println("probe_model_mode=$model_mode")
    println("probe_solve=$solve")
    println("probe_lower_model=$lower_model")
    println("probe_emit_fixture=$emit_fixture")
    println("probe_qmbcertify_base_construct=$qmbcertify_base_construct")
    println("probe_qmbcertify_base_extra=$qmbcertify_base_extra")
    println("probe_qmbcertify_base_three_type=$qmbcertify_base_three_type")
    println("probe_reflection_symmetry=$reflection_symmetry")
    println("probe_sign_symmetry=$sign_symmetry")
    println("probe_real_moment_matrix=$real_moment_matrix")
    println("probe_u1_symmetry=$u1_symmetry")
    println("probe_su2_symmetry=$su2_symmetry")
    println("probe_base_su2_extend_rdm=$base_su2_extend_rdm")
    println("probe_su2_moment_quotient=$probe_su2_moment_quotient")
    println("probe_axis_rotation_symmetry=$axis_rotation_symmetry")
    println("probe_axis_rotation_equalities=$axis_rotation_equalities")
    println("probe_axis_rotation_quotient=$axis_rotation_quotient")
    println("probe_contiguous_rdm_k=$contiguous_rdm_k")
    println("probe_contiguous_rdm_decomposition=$contiguous_rdm_decomposition")
    println("probe_contiguous_rdm_support=$contiguous_rdm_support")
    println("probe_psd_state_opt_width=$psd_state_opt_width")
    println("probe_linear_state_opt_width=$linear_state_opt_width")
    println("probe_linear_state_opt_mode=$linear_state_opt_mode")

    _, ops, pop = _heisenberg_case(n)
    linear_and_report, construct_stats = _timed("construct", () ->
        _solver_probe_linear_relaxation(
            pop,
            ops,
            order;
            qmbcertify_base_construct,
            qmbcertify_base_extra,
            qmbcertify_base_three_type,
            reflection_symmetry,
            sign_symmetry,
            real_moment_matrix,
            axis_rotation_symmetry,
            axis_rotation_equalities,
            axis_rotation_quotient,
            contiguous_rdm_k,
            contiguous_rdm_decomposition,
            contiguous_rdm_support,
            u1_symmetry,
            su2_symmetry,
            base_su2_extend_rdm,
            su2_moment_quotient=probe_su2_moment_quotient,
            linear_state_opt_width,
            linear_state_opt_mode,
            psd_state_opt_width,
        )
    )
    linear, report = linear_and_report
    support = translation_solve_support(report)
    metrics = translation_report_metrics(report)

    println("solve_supported=$(support.supported)")
    println("solve_support_reason=$(support.reason)")
    println("solve_support_blocker=$(support.blocker)")
    println("linear_moments=$(length(linear.moments))")
    println("linear_free_keys=$(length(linear.free_keys))")
    println("linear_zero_constraints=$(length(linear.zero_constraints))")
    println("psd_block_count=$(length(linear.psd_blocks_lin))")
    println("psd_block_histogram=$(_block_histogram(getfield.(linear.psd_blocks_lin, :size)))")
    println("report_su2_moment_symmetry=$(metrics.su2_moment_symmetry)")
    println("report_su2_rdm_symmetry=$(metrics.su2_rdm_symmetry)")
    println("report_u1_rdm_symmetry=$(metrics.u1_rdm_symmetry)")
    println("report_axis_rotation_symmetry=$(metrics.axis_rotation_symmetry)")
    println("report_psd_max_block=$(metrics.psd_max_block)")
    println("report_rdm_blocks=$(metrics.contiguous_rdm_psd_block_sizes)")
    println("report_pso_blocks=$(metrics.psd_state_opt_psd_block_sizes)")
    println("report_axis_quotient=$(metrics.axis_rotation_quotient)")
    println("report_base_moment_count=$(metrics.base_moment_count)")
    println("report_linear_moment_count=$(metrics.linear_moment_count)")
    println("report_su2_moment_quotient=$(metrics.su2_moment_quotient)")
    println("report_su2_moment_raw_count=$(metrics.su2_moment_raw_count)")
    println("report_su2_moment_quotient_count=$(metrics.su2_moment_quotient_count)")
    println(
        "report_su2_moment_quotient_reduction_ratio=" *
        "$(metrics.su2_moment_quotient_reduction_ratio)",
    )
    println(
        "report_su2_moment_support_orbit_count=" *
        "$(metrics.su2_moment_support_orbit_count)",
    )
    println(
        "report_su2_moment_singlet_channel_support_counts=" *
        "$(metrics.su2_moment_singlet_channel_support_counts)",
    )
    println("report_su2_moment_max_pivot_residual=$(metrics.su2_moment_max_pivot_residual)")
    println("report_su2_moment_max_invariant_residual=$(metrics.su2_moment_max_invariant_residual)")
    println(
        "report_su2_moment_max_reconstruction_residual=" *
        "$(metrics.su2_moment_max_reconstruction_residual)",
    )
    println("report_su2_moment_max_condition=$(metrics.su2_moment_max_condition)")
    println(
        "report_su2_moment_eliminated_zero_row_count=" *
        "$(metrics.su2_moment_eliminated_zero_row_count)",
    )
    println("report_su2_base_zero_row_count=$(metrics.su2_base_zero_row_count)")
    println("report_su2_base_spin_offblock_row_count=$(metrics.su2_base_spin_offblock_row_count)")
    println("report_su2_base_magnetic_offdiag_row_count=$(metrics.su2_base_magnetic_offdiag_row_count)")
    println("report_su2_base_magnetic_copy_row_count=$(metrics.su2_base_magnetic_copy_row_count)")
    println("report_product_cache_hit_rate=$(metrics.product_cache_hit_rate)")
    for stage in sort!(collect(keys(metrics.construction_stage_time_seconds)); by=string)
        println("report_stage_$(stage)_seconds=$(metrics.construction_stage_time_seconds[stage])")
    end
    model_size_estimate = _solver_probe_model_size_estimate(
        linear;
        formulation,
        representation,
        orphan_policy,
        dualize,
        sos_hermitian_representation,
    )
    model_size_gate = _solver_probe_model_size_gate(model_size_estimate)
    _print_model_size_estimate(
        model_size_estimate;
        model_size_gate,
    )
    if lower_model && _solver_probe_enforce_model_size_gate() &&
            model_size_gate.status != "ok"
        throw(ArgumentError(
            "Refusing solver probe lowering because estimated_model_size_gate_status=" *
            "$(model_size_gate.status) for dense-Schur proxy " *
            "$(model_size_estimate.dense_schur_bytes) bytes.  Lower the run size " *
            "or adjust NCTS_SOLVER_PROBE_MAX_RSS_FRACTION only with fresh evidence."
        ))
    end

    lower_model || begin
        if emit_fixture
            _emit_solver_probe_fixture_row(
                n=n,
                order=order,
                model_mode=model_mode,
                formulation=formulation,
                representation=representation,
                optimizer_name=optimizer_name,
                solve_support=support,
                reflection_symmetry=reflection_symmetry,
                sign_symmetry=sign_symmetry,
                real_moment_matrix=real_moment_matrix,
                u1_symmetry=u1_symmetry,
                su2_symmetry=su2_symmetry,
                base_su2_extend_rdm=base_su2_extend_rdm,
                su2_moment_quotient=probe_su2_moment_quotient,
                axis_rotation_symmetry=axis_rotation_symmetry,
                axis_rotation_equalities=axis_rotation_equalities,
                axis_rotation_quotient=axis_rotation_quotient,
                contiguous_rdm_k=contiguous_rdm_k,
                contiguous_rdm_decomposition=contiguous_rdm_decomposition,
                contiguous_rdm_support=contiguous_rdm_support,
                linear_state_opt_width=linear_state_opt_width,
                linear_state_opt_mode=linear_state_opt_mode,
                psd_state_opt_width=psd_state_opt_width,
                qmbcertify_base_construct=qmbcertify_base_construct,
                qmbcertify_base_extra=qmbcertify_base_extra,
                qmbcertify_base_three_type=qmbcertify_base_three_type,
                model_size_estimate=model_size_estimate,
                model_size_gate=model_size_gate,
                construct_stats=construct_stats,
                lower_stats=nothing,
                dualize_stats=nothing,
                solve_result=nothing,
                solve_stats=nothing,
                sos_certificate_residual=nothing,
                sos_certificate_diagnostics=nothing,
                sos_hermitian_representation=sos_hermitian_representation,
                sos_coefficient_scaling=sos_coefficient_scaling,
                sos_coefficient_scaling_floor=sos_coefficient_scaling_floor,
                script_wall_seconds=(time_ns() - probe_start_ns) / 1e9,
                report_metrics=metrics,
                linear=linear,
                dual_variables=nothing,
                dual_zero_duals=nothing,
                model_variables=nothing,
                scalar_equalities=nothing,
            )
        end
        return nothing
    end

    lower_stats = nothing
    dualize_stats = nothing
    solve_result = nothing
    solve_stats = nothing
    sos_certificate_residual = nothing
    sos_certificate_diagnostics = nothing
    dual_variables = nothing
    dual_zero_duals = nothing
    model_variables = nothing
    sos_problem = nothing

    model = if dualize
        sos_problem, dualize_stats = _timed(
            "dualize",
            () -> NCTSSoS.sos_dualize(
                linear;
                hermitian_representation=sos_hermitian_representation,
                coefficient_scaling=sos_coefficient_scaling,
                coefficient_scaling_floor=sos_coefficient_scaling_floor,
            ),
        )
        dual_model = sos_problem.model
        dual_variables = JuMP.num_variables(dual_model)
        dual_zero_duals = length(sos_problem.zero_duals)
        println("dual_model_variables=$(dual_variables)")
        println("dual_model_constraint_types=$(_constraint_type_counts(dual_model))")
        println("dual_psd_block_count=$(length(sos_problem.psd_dual_blocks))")
        println("dual_zero_dual_count=$(dual_zero_duals)")
        dual_model
    else
        model_extract, lower_stats = _timed("lower", () ->
            build_jump_model(
                linear;
                formulation,
                representation,
                orphan_policy,
            )
        )
        primal_model, _ = model_extract
        model_variables = JuMP.num_variables(primal_model)
        println("model_variables=$(model_variables)")
        println("model_constraint_types=$(_constraint_type_counts(primal_model))")
        primal_model
    end
    println("model_objective_sense=$(JuMP.objective_sense(model))")
    scalar_equalities = _scalar_equality_count(model)
    println("scalar_equalities=$scalar_equalities")
    println(
        "dense_schur_bytes=$(_solver_probe_dense_schur_bytes(scalar_equalities))",
    )

    if solve
        set_optimizer(model, _optimizer(optimizer_name))
        solve_result, solve_stats = _timed("solve", () -> begin
            JuMP.optimize!(model)
            return _solver_probe_solve_result(model)
        end)
        println("solve_status=$(solve_result.status)")
        println("solve_primal_status=$(solve_result.primal_status)")
        println("solve_dual_status=$(solve_result.dual_status)")
        println("solve_raw_status=$(solve_result.raw_status)")
        println("solve_result_count=$(solve_result.result_count)")
        println("solve_objective=$(solve_result.objective)")
        println("solve_objective_bound=$(solve_result.objective_bound)")
        println("solve_relative_gap=$(solve_result.relative_gap)")
        if dualize
            sos_certificate_diagnostics = NCTSSoS.sos_dual_certificate_diagnostics(
                linear,
                sos_problem,
            )
            sos_certificate_residual = sos_certificate_diagnostics.residual
            println("sos_certificate_moment_count=$(sos_certificate_residual.moment_count)")
            println(
                "sos_certificate_max_abs_residual=$(Float64(sos_certificate_residual.max_abs_residual))",
            )
            println(
                "sos_certificate_identity_abs_residual=$(_abs_float(sos_certificate_residual.identity_residual))",
            )
            println(
                "sos_certificate_worst_abs_residual=$(_abs_float(sos_certificate_residual.max_residual_value))",
            )
            println(
                "sos_certificate_psd_block_count=$(sos_certificate_diagnostics.psd_block_count)",
            )
            println(
                "sos_certificate_zero_dual_count=$(sos_certificate_diagnostics.zero_dual_count)",
            )
        end
    end

    if emit_fixture
        _emit_solver_probe_fixture_row(
            n=n,
            order=order,
            model_mode=model_mode,
            formulation=formulation,
            representation=representation,
            optimizer_name=optimizer_name,
            solve_support=support,
            reflection_symmetry=reflection_symmetry,
            sign_symmetry=sign_symmetry,
            real_moment_matrix=real_moment_matrix,
            u1_symmetry=u1_symmetry,
            su2_symmetry=su2_symmetry,
            base_su2_extend_rdm=base_su2_extend_rdm,
            su2_moment_quotient=probe_su2_moment_quotient,
            axis_rotation_symmetry=axis_rotation_symmetry,
            axis_rotation_equalities=axis_rotation_equalities,
            axis_rotation_quotient=axis_rotation_quotient,
            contiguous_rdm_k=contiguous_rdm_k,
            contiguous_rdm_decomposition=contiguous_rdm_decomposition,
            contiguous_rdm_support=contiguous_rdm_support,
            linear_state_opt_width=linear_state_opt_width,
            linear_state_opt_mode=linear_state_opt_mode,
            psd_state_opt_width=psd_state_opt_width,
            qmbcertify_base_construct=qmbcertify_base_construct,
            qmbcertify_base_extra=qmbcertify_base_extra,
            qmbcertify_base_three_type=qmbcertify_base_three_type,
            model_size_estimate=model_size_estimate,
            model_size_gate=model_size_gate,
            construct_stats=construct_stats,
            lower_stats=lower_stats,
            dualize_stats=dualize_stats,
            solve_result=solve_result,
            solve_stats=solve_stats,
            sos_certificate_residual=sos_certificate_residual,
            sos_certificate_diagnostics=sos_certificate_diagnostics,
            sos_hermitian_representation=sos_hermitian_representation,
            sos_coefficient_scaling=sos_coefficient_scaling,
            sos_coefficient_scaling_floor=sos_coefficient_scaling_floor,
            script_wall_seconds=(time_ns() - probe_start_ns) / 1e9,
            report_metrics=metrics,
            linear=linear,
            dual_variables=dual_variables,
            dual_zero_duals=dual_zero_duals,
            model_variables=model_variables,
            scalar_equalities=scalar_equalities,
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
