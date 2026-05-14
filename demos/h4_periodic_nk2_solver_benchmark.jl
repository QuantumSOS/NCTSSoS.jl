#!/usr/bin/env julia
raw"""
Benchmark BPSDP.jl, SCS.jl indirect/CG, and COSMO.jl CGIndirectKKTSolver on
H4/Nk=2 periodic PQG V2RDM.

This reuses `h4_periodic_moment_sos.jl` to build the same symbolic
MomentProblem for each solver.  Each solver gets a fresh JuMP lowering:

    build_jump_model(...;
        formulation = :psd_blocks,
        representation = :complex,
        orphan_policy = :aux_psd_free,
    )

The default target is 1e-5 solver accuracy.  This is a benchmark/probe, not a
production driver; H4 can run for a long time.

Usage from the repository root on HAI:

    tmpdir=$(mktemp -d)
    julia --startup-file=no --project="$tmpdir" -e 'using Pkg; \
        Pkg.develop(path=pwd()); \
        Pkg.develop(path="/home/ubuntu/BPSDP.jl"); \
        Pkg.add(["JuMP", "JSON3", "SCS", "COSMO", "IterativeSolvers", "LinearMaps"]); \
        Pkg.instantiate()'
    julia --startup-file=no --project="$tmpdir" \
        demos/h4_periodic_nk2_solver_benchmark.jl \
        --solvers=bpsdp,scs-indirect,cosmo-cg \
        --output-dir=output/phase2/h4_nk2_solver_benchmark_1e5
"""

using NCTSSoS
using JuMP
const MOI = JuMP.MOI
using LinearAlgebra
using SparseArrays
using Printf
using Dates
using JSON3

import BPSDP
import SCS
import COSMO
import IterativeSolvers
import LinearMaps

function install_cosmo_quiet_processconstraints!()
    isdefined(COSMO, :processconstraints!) || return false
    @eval COSMO begin
        function processconstraints!(optimizer::Optimizer{T}, src::MOI.ModelLike) where {T <: AbstractFloat}
            convex_sets = COSMO.AbstractConvexSet{T}[]
            for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
                processSets!(convex_sets, src.constraints, F, S)
            end

            if isempty(convex_sets)
                push!(convex_sets, COSMO.ZeroSet{T}(0))
            end

            model = optimizer.inner
            model.p.A = -Base.convert(SparseMatrixCSC{T,Int}, src.constraints.coefficients)
            model.p.b = src.constraints.constants.b
            model.p.C = CompositeConvexSet{T}(deepcopy(convex_sets))
            model.p.model_size = [size(model.p.A, 1), size(model.p.A, 2)]
            return nothing
        end
    end
    return true
end

const COSMO_QUIET_PROCESSCONSTRAINTS_PATCHED = install_cosmo_quiet_processconstraints!()

module H4PeriodicBase
include(joinpath(@__DIR__, "h4_periodic_moment_sos.jl"))
end

const DEFAULT_BENCHMARK_OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "output", "phase2", "h4_nk2_solver_benchmark_1e5"))

struct H4BuildOptions
    integrals_path::String
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
    nk::Int
    norb::Int
    nelec_per_cell::Int
end

struct BPSDPOptions
    max_iter::Int
    cg_max_iter::Int
    mu_update_frequency::Int
    penalty::Float64
    cg_tol::Float64
    obj_tol::Float64
    err_tol::Float64
    print_level::Int
end

struct SCSOptions
    max_iters::Int
    eps_abs::Float64
    eps_rel::Float64
    eps_infeas::Float64
    alpha::Float64
    normalize::Int
    adaptive_scale::Int
    scale::Float64
    rho_x::Float64
    acceleration_lookback::Int
    acceleration_interval::Int
    time_limit_secs::Float64
    verbose::Int
end

struct COSMOOptions
    max_iter::Int
    eps_abs::Float64
    eps_rel::Float64
    eps_prim_inf::Float64
    eps_dual_inf::Float64
    rho::Float64
    sigma::Float64
    alpha::Float64
    scaling::Int
    check_termination::Int
    check_infeasibility::Int
    time_limit::Float64
    verbose::Bool
    verbose_timing::Bool
    decompose::Bool
    compact_transformation::Bool
    complete_dual::Bool
    tol_constant::Float64
    tol_exponent::Float64
end

struct BenchmarkOptions
    build::H4BuildOptions
    solvers::Vector{Symbol}
    repeats::Int
    output_dir::String
    bpsdp::BPSDPOptions
    scs::SCSOptions
    cosmo::COSMOOptions
end

bytes_to_gib(bytes::Integer) = bytes / 1024.0^3
rss_string() = @sprintf("%.3f GiB", bytes_to_gib(Sys.maxrss()))

function parse_bool_int(value::AbstractString)
    v = lowercase(strip(value))
    v in ("1", "true", "yes", "on") && return 1
    v in ("0", "false", "no", "off") && return 0
    throw(ArgumentError("expected boolean/int 0|1|true|false, got $(repr(value))"))
end

parse_bool(value::AbstractString) = parse_bool_int(value) == 1

function parse_finite_float(name::AbstractString, value::AbstractString; min::Union{Nothing,Float64}=nothing, max::Union{Nothing,Float64}=nothing, inclusive_min::Bool=true, inclusive_max::Bool=true)
    x = parse(Float64, value)
    isfinite(x) || throw(ArgumentError("$name must be finite, got $(repr(value))"))
    if min !== nothing
        ok = inclusive_min ? x >= min : x > min
        ok || throw(ArgumentError("$name must be $(inclusive_min ? ">=" : ">") $min, got $x"))
    end
    if max !== nothing
        ok = inclusive_max ? x <= max : x < max
        ok || throw(ArgumentError("$name must be $(inclusive_max ? "<=" : "<") $max, got $x"))
    end
    return x
end

parse_positive_float(name::AbstractString, value::AbstractString) =
    parse_finite_float(name, value; min = 0.0, inclusive_min = false)
parse_nonnegative_float(name::AbstractString, value::AbstractString) =
    parse_finite_float(name, value; min = 0.0)

function parse_blocking(value::AbstractString)
    v = lowercase(strip(value))
    v in ("momentum", "k") && return :momentum
    v == "spin" && return :spin
    v == "none" && return :none
    throw(ArgumentError("unknown --blocking=$value; expected momentum, spin, or none"))
end

function parse_solvers(value::AbstractString)
    solvers = Symbol[]
    for raw in split(value, ',')
        token = lowercase(strip(raw))
        isempty(token) && continue
        if token == "bpsdp"
            push!(solvers, :bpsdp)
        elseif token in ("scs", "scs-indirect", "scs_indirect")
            push!(solvers, :scs_indirect)
        elseif token in ("cosmo", "cosmo-cg", "cosmo_cg", "cosmo-indirect", "cosmo_indirect")
            push!(solvers, :cosmo_cg)
        elseif token == "all"
            append!(solvers, (:bpsdp, :scs_indirect, :cosmo_cg))
        else
            throw(ArgumentError("unknown solver $(repr(token)); expected bpsdp, scs-indirect, cosmo-cg, or all"))
        end
    end
    isempty(solvers) && throw(ArgumentError("--solvers must name at least one solver"))
    return unique(solvers)
end

function print_help()
    println("Usage: julia --project=<env> demos/h4_periodic_nk2_solver_benchmark.jl [options]\n")
    println("Build options:")
    println("  --integrals=PATH                  default test/data/assets/h4_chain_nk2_integrals.txt")
    println("  --nk=N                            number of k-points; default 2")
    println("  --norb=N                          active spatial orbitals per k-point; default 8")
    println("  --nelec-per-cell=N                active electrons per unit cell; default 4")
    println("  --blocking=momentum|spin|none     default momentum")
    println("  --include-1d / --no-1d            include/drop extra ¹D PSD block; default --no-1d")
    println("  --paper-spin, --spin-singlet      spin-resolved ²D traces plus singlet S²")
    println("  --spin-resolved-trace / --no-spin-resolved-trace")
    println("  --singlet-s2 / --no-singlet-s2")
    println("  --output-dir=PATH                 default output/phase2/h4_nk2_solver_benchmark_1e5")
    println("\nBenchmark options:")
    println("  --solvers=bpsdp,scs-indirect,cosmo-cg  default all")
    println("  --repeats=N                       fresh model per solver per repeat; default 1")
    println("\nBPSDP options:")
    println("  --bpsdp-max-iter=N                default 50000")
    println("  --bpsdp-cg-max-iter=N             default 10000")
    println("  --bpsdp-mu-update-frequency=N     default 500")
    println("  --bpsdp-penalty=rho               default 0.1")
    println("  --bpsdp-cg-tol=eps                default 1e-10")
    println("  --bpsdp-obj-tol=eps               default 1e-5")
    println("  --bpsdp-err-tol=eps               default 1e-5")
    println("  --bpsdp-print-level=N             default 1")
    println("\nSCS indirect/CG options:")
    println("  --scs-max-iters=N                 default 200000")
    println("  --scs-eps-abs=eps                 default 1e-5")
    println("  --scs-eps-rel=eps                 default 1e-5")
    println("  --scs-time-limit-secs=sec         default 0.0, no limit")
    println("\nCOSMO CGIndirectKKTSolver options:")
    println("  --cosmo-max-iter=N                default 200000")
    println("  --cosmo-eps-abs=eps               default 1e-5")
    println("  --cosmo-eps-rel=eps               default 1e-5")
    println("  --cosmo-time-limit=sec            default 0.0, no limit")
    println("  --cosmo-verbose=true|false        default false; keep false for H4")
    println("  --cosmo-verbose-timing=true|false default false")
    println("  --cosmo-tol-constant=x            default 1.0")
    println("  --cosmo-tol-exponent=x            default 1.5")
    return nothing
end

function parse_benchmark_options(argv)
    integrals_path = H4PeriodicBase.DEFAULT_INTEGRALS
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = false
    singlet_s2 = false
    output_dir = DEFAULT_BENCHMARK_OUTPUT_DIR
    nk = 2
    norb = 8
    nelec_per_cell = 4

    bpsdp_max_iter = 50_000
    bpsdp_cg_max_iter = 10_000
    bpsdp_mu_update_frequency = 500
    bpsdp_penalty = 0.1
    bpsdp_cg_tol = 1e-10
    bpsdp_obj_tol = 1e-5
    bpsdp_err_tol = 1e-5
    bpsdp_print_level = 1

    solvers = [:bpsdp, :scs_indirect, :cosmo_cg]
    repeats = 1

    scs_max_iters = 200_000
    scs_eps_abs = 1e-5
    scs_eps_rel = 1e-5
    scs_eps_infeas = 1e-7
    scs_alpha = 1.5
    scs_normalize = 1
    scs_adaptive_scale = 1
    scs_scale = 0.1
    scs_rho_x = 1e-6
    scs_acceleration_lookback = 10
    scs_acceleration_interval = 10
    scs_time_limit_secs = 0.0
    scs_verbose = 1

    cosmo_max_iter = 200_000
    cosmo_eps_abs = 1e-5
    cosmo_eps_rel = 1e-5
    cosmo_eps_prim_inf = 1e-4
    cosmo_eps_dual_inf = 1e-4
    cosmo_rho = 0.1
    cosmo_sigma = 1e-6
    cosmo_alpha = 1.6
    cosmo_scaling = 10
    cosmo_check_termination = 25
    cosmo_check_infeasibility = 40
    cosmo_time_limit = 0.0
    cosmo_verbose = false
    cosmo_verbose_timing = false
    cosmo_decompose = false
    cosmo_compact_transformation = true
    cosmo_complete_dual = false
    cosmo_tol_constant = 1.0
    cosmo_tol_exponent = 1.5

    for arg in argv
        if startswith(arg, "--integrals=")
            integrals_path = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--nk=")
            nk = parse(Int, split(arg, "=", limit = 2)[2])
            nk >= 1 || throw(ArgumentError("--nk must be positive"))
        elseif startswith(arg, "--norb=")
            norb = parse(Int, split(arg, "=", limit = 2)[2])
            norb >= 1 || throw(ArgumentError("--norb must be positive"))
        elseif startswith(arg, "--nelec-per-cell=")
            nelec_per_cell = parse(Int, split(arg, "=", limit = 2)[2])
            nelec_per_cell >= 1 || throw(ArgumentError("--nelec-per-cell must be positive"))
        elseif startswith(arg, "--blocking=")
            blocking = parse_blocking(split(arg, "=", limit = 2)[2])
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
        elseif startswith(arg, "--solvers=")
            solvers = parse_solvers(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--repeats=")
            repeats = parse(Int, split(arg, "=", limit = 2)[2])
            repeats >= 1 || throw(ArgumentError("--repeats must be positive"))
        elseif startswith(arg, "--bpsdp-max-iter=")
            bpsdp_max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-cg-max-iter=")
            bpsdp_cg_max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-mu-update-frequency=")
            bpsdp_mu_update_frequency = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-penalty=")
            bpsdp_penalty = parse_positive_float("--bpsdp-penalty", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-cg-tol=")
            bpsdp_cg_tol = parse_positive_float("--bpsdp-cg-tol", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-obj-tol=")
            bpsdp_obj_tol = parse_positive_float("--bpsdp-obj-tol", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-err-tol=")
            bpsdp_err_tol = parse_positive_float("--bpsdp-err-tol", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--bpsdp-print-level=")
            bpsdp_print_level = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-max-iters=")
            scs_max_iters = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-eps-abs=")
            scs_eps_abs = parse_positive_float("--scs-eps-abs", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-eps-rel=")
            scs_eps_rel = parse_positive_float("--scs-eps-rel", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-eps-infeas=")
            scs_eps_infeas = parse_positive_float("--scs-eps-infeas", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-alpha=")
            scs_alpha = parse_finite_float("--scs-alpha", split(arg, "=", limit = 2)[2]; min = 0.0, max = 2.0, inclusive_min = false, inclusive_max = false)
        elseif startswith(arg, "--scs-normalize=")
            scs_normalize = parse_bool_int(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-adaptive-scale=")
            scs_adaptive_scale = parse_bool_int(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-scale=")
            scs_scale = parse_positive_float("--scs-scale", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-rho-x=")
            scs_rho_x = parse_positive_float("--scs-rho-x", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-acceleration-lookback=")
            scs_acceleration_lookback = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-acceleration-interval=")
            scs_acceleration_interval = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-time-limit-secs=")
            scs_time_limit_secs = parse_nonnegative_float("--scs-time-limit-secs", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--scs-verbose=")
            scs_verbose = parse_bool_int(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-max-iter=")
            cosmo_max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-eps-abs=")
            cosmo_eps_abs = parse_positive_float("--cosmo-eps-abs", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-eps-rel=")
            cosmo_eps_rel = parse_positive_float("--cosmo-eps-rel", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-eps-prim-inf=")
            cosmo_eps_prim_inf = parse_positive_float("--cosmo-eps-prim-inf", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-eps-dual-inf=")
            cosmo_eps_dual_inf = parse_positive_float("--cosmo-eps-dual-inf", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-rho=")
            cosmo_rho = parse_positive_float("--cosmo-rho", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-sigma=")
            cosmo_sigma = parse_positive_float("--cosmo-sigma", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-alpha=")
            cosmo_alpha = parse_finite_float("--cosmo-alpha", split(arg, "=", limit = 2)[2]; min = 0.0, max = 2.0, inclusive_min = false, inclusive_max = false)
        elseif startswith(arg, "--cosmo-scaling=")
            cosmo_scaling = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-check-termination=")
            cosmo_check_termination = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-check-infeasibility=")
            cosmo_check_infeasibility = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-time-limit=")
            cosmo_time_limit = parse_nonnegative_float("--cosmo-time-limit", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-verbose=")
            cosmo_verbose = parse_bool(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-verbose-timing=")
            cosmo_verbose_timing = parse_bool(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-decompose=")
            cosmo_decompose = parse_bool(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-compact-transformation=")
            cosmo_compact_transformation = parse_bool(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-complete-dual=")
            cosmo_complete_dual = parse_bool(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-tol-constant=")
            cosmo_tol_constant = parse_positive_float("--cosmo-tol-constant", split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cosmo-tol-exponent=")
            cosmo_tol_exponent = parse_positive_float("--cosmo-tol-exponent", split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            print_help()
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    build = H4BuildOptions(integrals_path, blocking, include_one_d, spin_resolved_trace,
        singlet_s2, nk, norb, nelec_per_cell)
    bpsdp = BPSDPOptions(bpsdp_max_iter, bpsdp_cg_max_iter, bpsdp_mu_update_frequency,
        bpsdp_penalty, bpsdp_cg_tol, bpsdp_obj_tol, bpsdp_err_tol, bpsdp_print_level)
    scs = SCSOptions(scs_max_iters, scs_eps_abs, scs_eps_rel, scs_eps_infeas, scs_alpha,
        scs_normalize, scs_adaptive_scale, scs_scale, scs_rho_x,
        scs_acceleration_lookback, scs_acceleration_interval, scs_time_limit_secs, scs_verbose)
    cosmo = COSMOOptions(cosmo_max_iter, cosmo_eps_abs, cosmo_eps_rel, cosmo_eps_prim_inf,
        cosmo_eps_dual_inf, cosmo_rho, cosmo_sigma, cosmo_alpha, cosmo_scaling,
        cosmo_check_termination, cosmo_check_infeasibility, cosmo_time_limit,
        cosmo_verbose, cosmo_verbose_timing, cosmo_decompose, cosmo_compact_transformation,
        cosmo_complete_dual, cosmo_tol_constant, cosmo_tol_exponent)
    return BenchmarkOptions(build, solvers, repeats, output_dir, bpsdp, scs, cosmo)
end

function h4_options(build::H4BuildOptions)
    return H4PeriodicBase.Options(build.integrals_path, build.blocking, build.include_one_d,
        build.spin_resolved_trace, build.singlet_s2,
        build.nk, build.norb, build.nelec_per_cell)
end

has_field(obj, name::Symbol) = name in fieldnames(typeof(obj))
field_or_nothing(obj, name::Symbol) = has_field(obj, name) ? getfield(obj, name) : nothing
function first_field_or_nothing(obj, names::Symbol...)
    for name in names
        has_field(obj, name) && return getfield(obj, name)
    end
    return nothing
end

function finite_or_string(x::Real)
    xf = Float64(x)
    return isfinite(xf) ? xf : string(xf)
end
finite_or_nothing(x) = x === nothing ? nothing : finite_or_string(x)

function json_ready(x)
    x === nothing && return nothing
    x isa Missing && return nothing
    x isa AbstractString && return x
    x isa Symbol && return string(x)
    x isa Bool && return x
    x isa Integer && return x
    x isa AbstractFloat && return isfinite(x) ? x : string(x)
    x isa Real && return finite_or_string(x)
    x isa DateTime && return string(x)
    if x isa Pair
        return string(x.first) => json_ready(x.second)
    elseif x isa NamedTuple
        return Dict(String(k) => json_ready(v) for (k, v) in pairs(x))
    elseif x isa AbstractDict
        return Dict(string(k) => json_ready(v) for (k, v) in x)
    elseif x isa Tuple
        return Any[json_ready(v) for v in x]
    elseif x isa AbstractVector
        return Any[json_ready(v) for v in x]
    elseif x isa AbstractMatrix
        return Any[Any[json_ready(x[i, j]) for j in axes(x, 2)] for i in axes(x, 1)]
    end
    return string(x)
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, json_ready(obj))
        write(io, '\n')
    end
    return path
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

function bpsdp_optimizer_factory(options::BPSDPOptions, cg_iterations::Vector{Int})
    monitor = progress_monitor_factory(cg_iterations)
    return () -> BPSDP.Optimizer(
        max_iter = options.max_iter,
        cg_max_iter = options.cg_max_iter,
        mu_update_frequency = options.mu_update_frequency,
        penalty_parameter = options.penalty,
        cg_convergence = options.cg_tol,
        dynamic_cg_convergence = true,
        sdp_objective_convergence = options.obj_tol,
        sdp_error_convergence = options.err_tol,
        guess_type = :zero,
        print_level = options.print_level,
        progress_monitor = monitor,
        dependent_rows = :keep,
    )
end

function scs_indirect_optimizer_factory(options::SCSOptions)
    SCS.is_available(SCS.IndirectSolver) || error("SCS.IndirectSolver is not available in this SCS.jl build")
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "linear_solver" => SCS.IndirectSolver,
        "max_iters" => options.max_iters,
        "eps_abs" => options.eps_abs,
        "eps_rel" => options.eps_rel,
        "eps_infeas" => options.eps_infeas,
        "alpha" => options.alpha,
        "normalize" => options.normalize,
        "adaptive_scale" => options.adaptive_scale,
        "scale" => options.scale,
        "rho_x" => options.rho_x,
        "acceleration_lookback" => options.acceleration_lookback,
        "acceleration_interval" => options.acceleration_interval,
        "time_limit_secs" => options.time_limit_secs,
        "verbose" => options.verbose,
    )
end

function cosmo_cg_optimizer_factory(options::COSMOOptions)
    # Importing IterativeSolvers and LinearMaps above triggers COSMO's optional
    # indirect-KKT implementation via Requires.jl.
    isdefined(COSMO, :CGIndirectKKTSolver) || error("COSMO.CGIndirectKKTSolver is unavailable; load IterativeSolvers and LinearMaps in the active env")
    kkt = COSMO.with_options(COSMO.CGIndirectKKTSolver,
        tol_constant = options.tol_constant,
        tol_exponent = options.tol_exponent,
    )
    return JuMP.optimizer_with_attributes(
        COSMO.Optimizer,
        "kkt_solver" => kkt,
        "max_iter" => options.max_iter,
        "eps_abs" => options.eps_abs,
        "eps_rel" => options.eps_rel,
        "eps_prim_inf" => options.eps_prim_inf,
        "eps_dual_inf" => options.eps_dual_inf,
        "rho" => options.rho,
        "sigma" => options.sigma,
        "alpha" => options.alpha,
        "scaling" => options.scaling,
        "check_termination" => options.check_termination,
        "check_infeasibility" => options.check_infeasibility,
        "time_limit" => options.time_limit,
        "verbose" => options.verbose,
        "verbose_timing" => options.verbose_timing,
        "decompose" => options.decompose,
        "compact_transformation" => options.compact_transformation,
        "complete_dual" => options.complete_dual,
    )
end

function set_solver!(model, solver::Symbol, options::BenchmarkOptions, cg_iterations::Vector{Int})
    if solver == :bpsdp
        JuMP.set_optimizer(model, bpsdp_optimizer_factory(options.bpsdp, cg_iterations))
    elseif solver == :scs_indirect
        JuMP.set_optimizer(model, scs_indirect_optimizer_factory(options.scs))
    elseif solver == :cosmo_cg
        JuMP.set_optimizer(model, cosmo_cg_optimizer_factory(options.cosmo))
        # COSMO's MOI wrapper defaults to verbose output. On H4 that can mean
        # dumping the entire JuMP/MOI model before the solver loop. That's not
        # logging, that's vandalism.
        JuMP.set_silent(model)
    else
        error("unknown solver $solver")
    end
    return nothing
end

function raw_optimizer(model::JuMP.Model)
    try
        return JuMP.unsafe_backend(model)
    catch
        return nothing
    end
end

function safe_termination_status(model)
    try
        return string(JuMP.termination_status(model))
    catch err
        return string("unavailable: ", truncated_error(err))
    end
end

function safe_raw_status(model)
    try
        return string(MOI.get(model, MOI.RawStatusString()))
    catch err
        return string("unavailable: ", truncated_error(err))
    end
end

function safe_objective_value(model)
    try
        return JuMP.objective_value(model)
    catch err
        return string("unavailable: ", truncated_error(err))
    end
end

function safe_dual_objective_value(model)
    try
        return MOI.get(model, MOI.DualObjectiveValue())
    catch err
        return string("unavailable: ", truncated_error(err))
    end
end

function safe_solve_time(model)
    try
        return finite_or_string(MOI.get(model, MOI.SolveTimeSec()))
    catch err
        return string("unavailable: ", truncated_error(err))
    end
end

struct TruncatedShowError <: Exception end

mutable struct TruncatingIO <: IO
    buffer::IOBuffer
    limit::Int
    written::Int
    truncated::Bool
end

TruncatingIO(limit::Int) = TruncatingIO(IOBuffer(), limit, 0, false)
Base.iswritable(::TruncatingIO) = true
Base.isopen(::TruncatingIO) = true
Base.close(::TruncatingIO) = nothing

function Base.unsafe_write(io::TruncatingIO, p::Ptr{UInt8}, nb::UInt)
    n = Int(nb)
    remaining = max(io.limit - io.written, 0)
    if remaining > 0
        Base.unsafe_write(io.buffer, p, UInt(min(n, remaining)))
    end
    io.written += n
    if n > remaining
        io.truncated = true
        throw(TruncatedShowError())
    end
    return n
end

function truncated_error(err; limit::Int = 2000)
    io = TruncatingIO(limit)
    try
        showerror(IOContext(io, :limit => true, :compact => true), err)
    catch caught
        if caught isa TruncatedShowError
            msg = String(take!(io.buffer))
            return string(msg, "\n... [truncated after ", limit, " bytes; original error type ", typeof(err), "]")
        end
        return string(typeof(err), " (showerror failed with ", typeof(caught), ")")
    end
    return String(take!(io.buffer))
end

function bpsdp_state_summary(raw)
    state = raw === nothing ? nothing : field_or_nothing(raw, :state)
    state === nothing && return nothing
    return Dict{String,Any}(
        "outer_iterations" => first_field_or_nothing(state, :outer_iterations, :iterations),
        "inner_iterations" => field_or_nothing(state, :inner_iterations),
        "primal_error" => finite_or_nothing(field_or_nothing(state, :primal_error)),
        "dual_error" => finite_or_nothing(field_or_nothing(state, :dual_error)),
        "objective_primal" => finite_or_nothing(first_field_or_nothing(state, :objective_primal, :primal_objective)),
        "objective_dual" => finite_or_nothing(first_field_or_nothing(state, :objective_dual, :dual_objective)),
        "objective_gap" => finite_or_nothing(field_or_nothing(state, :objective_gap)),
        "termination_reason" => string(first_field_or_nothing(state, :termination_reason, :status)),
    )
end

function scs_state_summary(raw)
    raw === nothing && return nothing
    sol = field_or_nothing(raw, :sol)
    sol === nothing && return nothing
    return Dict{String,Any}(
        "iterations" => field_or_nothing(sol, :iterations),
        "solve_time_sec_reported" => finite_or_nothing(field_or_nothing(sol, :solve_time_sec)),
        "primal_objective" => finite_or_nothing(field_or_nothing(sol, :objective_value)),
        "dual_objective" => finite_or_nothing(field_or_nothing(sol, :dual_objective_value)),
        "objective_gap" => begin
            po = field_or_nothing(sol, :objective_value)
            do_ = field_or_nothing(sol, :dual_objective_value)
            (po isa Real && do_ isa Real) ? finite_or_string(abs(po - do_)) : nothing
        end,
        "raw_status" => string(field_or_nothing(sol, :raw_status)),
        "ret_val" => field_or_nothing(sol, :ret_val),
    )
end

function result_info_summary(info)
    info === nothing && return nothing
    return Dict{String,Any}(
        "r_prim" => finite_or_nothing(field_or_nothing(info, :r_prim)),
        "r_dual" => finite_or_nothing(field_or_nothing(info, :r_dual)),
        "max_norm_prim" => finite_or_nothing(field_or_nothing(info, :max_norm_prim)),
        "max_norm_dual" => finite_or_nothing(field_or_nothing(info, :max_norm_dual)),
        "rho_updates_count" => begin
            updates = field_or_nothing(info, :rho_updates)
            updates === nothing ? nothing : length(updates)
        end,
    )
end

function result_times_summary(times)
    times === nothing && return nothing
    return Dict{String,Any}(
        "solver_time" => finite_or_nothing(field_or_nothing(times, :solver_time)),
        "setup_time" => finite_or_nothing(field_or_nothing(times, :setup_time)),
        "scaling_time" => finite_or_nothing(field_or_nothing(times, :scaling_time)),
        "graph_time" => finite_or_nothing(field_or_nothing(times, :graph_time)),
        "init_factor_time" => finite_or_nothing(field_or_nothing(times, :init_factor_time)),
        "factor_update_time" => finite_or_nothing(field_or_nothing(times, :factor_update_time)),
        "iter_time" => finite_or_nothing(field_or_nothing(times, :iter_time)),
        "proj_time" => finite_or_nothing(field_or_nothing(times, :proj_time)),
        "post_time" => finite_or_nothing(field_or_nothing(times, :post_time)),
        "update_time" => finite_or_nothing(field_or_nothing(times, :update_time)),
        "accelerate_time" => finite_or_nothing(field_or_nothing(times, :accelerate_time)),
    )
end

function cosmo_kkt_summary(raw)
    inner = raw === nothing ? nothing : field_or_nothing(raw, :inner)
    kkt = inner === nothing ? nothing : field_or_nothing(inner, :kkt_solver)
    indirect = kkt === nothing ? nothing : field_or_nothing(kkt, :indirect_kktsolver)
    indirect === nothing && return nothing
    multiplications = field_or_nothing(indirect, :multiplications)
    return Dict{String,Any}(
        "type" => string(typeof(kkt)),
        "inner_type" => string(typeof(indirect)),
        "linear_solves" => multiplications === nothing ? nothing : length(multiplications),
        "matvecs_total" => multiplications === nothing ? nothing : sum(multiplications),
        "matvecs_max" => multiplications === nothing || isempty(multiplications) ? nothing : maximum(multiplications),
        "matvecs_last" => multiplications === nothing || isempty(multiplications) ? nothing : last(multiplications),
    )
end

function cosmo_state_summary(raw)
    raw === nothing && return nothing
    res = field_or_nothing(raw, :results)
    res === nothing && return nothing
    info = field_or_nothing(res, :info)
    times = field_or_nothing(res, :times)
    return Dict{String,Any}(
        "status" => string(field_or_nothing(res, :status)),
        "iterations" => field_or_nothing(res, :iter),
        "safeguarding_iterations" => field_or_nothing(res, :safeguarding_iter),
        "objective" => finite_or_nothing(field_or_nothing(res, :obj_val)),
        "info" => result_info_summary(info),
        "times" => result_times_summary(times),
        "kkt" => cosmo_kkt_summary(raw),
    )
end

function solver_options_json(solver::Symbol, options::BenchmarkOptions)
    if solver == :bpsdp
        b = options.bpsdp
        return Dict{String,Any}(
            "max_iter" => b.max_iter,
            "cg_max_iter" => b.cg_max_iter,
            "mu_update_frequency" => b.mu_update_frequency,
            "penalty_parameter" => b.penalty,
            "cg_convergence" => b.cg_tol,
            "dynamic_cg_convergence" => true,
            "sdp_objective_convergence" => b.obj_tol,
            "sdp_error_convergence" => b.err_tol,
            "guess_type" => "zero",
            "dependent_rows" => "keep",
        )
    elseif solver == :scs_indirect
        s = options.scs
        return Dict{String,Any}(
            "linear_solver" => "SCS.IndirectSolver (conjugate gradient)",
            "max_iters" => s.max_iters,
            "eps_abs" => s.eps_abs,
            "eps_rel" => s.eps_rel,
            "eps_infeas" => s.eps_infeas,
            "alpha" => s.alpha,
            "normalize" => s.normalize,
            "adaptive_scale" => s.adaptive_scale,
            "scale" => s.scale,
            "rho_x" => s.rho_x,
            "acceleration_lookback" => s.acceleration_lookback,
            "acceleration_interval" => s.acceleration_interval,
            "time_limit_secs" => s.time_limit_secs,
            "verbose" => s.verbose,
        )
    elseif solver == :cosmo_cg
        c = options.cosmo
        return Dict{String,Any}(
            "kkt_solver" => "COSMO.CGIndirectKKTSolver",
            "kkt_tol_constant" => c.tol_constant,
            "kkt_tol_exponent" => c.tol_exponent,
            "max_iter" => c.max_iter,
            "eps_abs" => c.eps_abs,
            "eps_rel" => c.eps_rel,
            "eps_prim_inf" => c.eps_prim_inf,
            "eps_dual_inf" => c.eps_dual_inf,
            "rho" => c.rho,
            "sigma" => c.sigma,
            "alpha" => c.alpha,
            "scaling" => c.scaling,
            "check_termination" => c.check_termination,
            "check_infeasibility" => c.check_infeasibility,
            "time_limit" => c.time_limit,
            "verbose" => c.verbose,
            "verbose_timing" => c.verbose_timing,
            "decompose" => c.decompose,
            "compact_transformation" => c.compact_transformation,
            "complete_dual" => c.complete_dual,
            "cosmo_processconstraints_patch" => COSMO_QUIET_PROCESSCONSTRAINTS_PATCHED,
            "cosmo_processconstraints_patch_note" => "benchmark override removes COSMO.MOI_wrapper processconstraints! debug println(src)",
        )
    end
    error("unknown solver $solver")
end

function solver_label(solver::Symbol)
    solver == :scs_indirect && return "scs_indirect"
    solver == :cosmo_cg && return "cosmo_cg"
    return string(solver)
end

function common_runtime_json()
    return Dict{String,Any}(
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "julia_num_threads_env" => get(ENV, "JULIA_NUM_THREADS", nothing),
        "blas_threads" => BLAS.get_num_threads(),
        "cpu_threads" => Sys.CPU_THREADS,
        "hostname" => gethostname(),
        "max_rss" => rss_string(),
    )
end

function truncate_log_file!(path::AbstractString; max_bytes::Int = 262_144)
    isfile(path) || return Dict{String,Any}(
        "path" => path,
        "exists" => false,
        "truncated" => false,
    )

    original_bytes = filesize(path)
    if original_bytes <= max_bytes
        return Dict{String,Any}(
            "path" => path,
            "exists" => true,
            "bytes" => original_bytes,
            "truncated" => false,
        )
    end

    head_bytes = max_bytes ÷ 2
    tail_bytes = max_bytes - head_bytes
    head = UInt8[]
    tail = UInt8[]
    open(path, "r") do io
        head = read(io, head_bytes)
        seek(io, max(Int(original_bytes) - tail_bytes, 0))
        tail = read(io)
    end
    marker = codeunits("\n... [log truncated from $(original_bytes) to $(max_bytes) bytes by benchmark wrapper] ...\n")
    open(path, "w") do io
        write(io, head)
        write(io, marker)
        write(io, tail)
    end

    return Dict{String,Any}(
        "path" => path,
        "exists" => true,
        "bytes" => filesize(path),
        "original_bytes" => original_bytes,
        "truncated" => true,
    )
end

function optimize_with_redirected_stdio!(model, stdout_path::AbstractString, stderr_path::AbstractString)
    mkpath(dirname(stdout_path))
    mkpath(dirname(stderr_path))

    optimize_error = nothing
    open(stdout_path, "w") do out
        open(stderr_path, "w") do errio
            try
                redirect_stdout(out) do
                    redirect_stderr(errio) do
                        JuMP.optimize!(model)
                    end
                end
            catch err
                optimize_error = truncated_error(err)
            end
        end
    end

    return optimize_error, Dict{String,Any}(
        "mode" => "file",
        "stdout" => truncate_log_file!(stdout_path),
        "stderr" => truncate_log_file!(stderr_path),
    )
end

function optimize_solver!(model, solver::Symbol, options::BenchmarkOptions; repeat::Int)
    if solver == :cosmo_cg
        # COSMO/MOI failures can print the whole JuMP/MOI model. Keep that
        # garbage out of .easy-ssh-log; retain bounded per-run artifacts.
        run_dir = joinpath(options.output_dir, solver_label(solver), "repeat_$(repeat)")
        return optimize_with_redirected_stdio!(
            model,
            joinpath(run_dir, "cosmo_stdout.log"),
            joinpath(run_dir, "cosmo_stderr.log"),
        )
    end

    try
        JuMP.optimize!(model)
        return nothing, nothing
    catch err
        return truncated_error(err), nothing
    end
end

function solve_one(data, solver::Symbol, options::BenchmarkOptions; repeat::Int)
    started_at = Dates.now()
    result = Dict{String,Any}(
        "solver" => string(solver),
        "repeat" => repeat,
        "started_at" => string(started_at),
        "jump_lowering" => Dict(
            "formulation" => "psd_blocks",
            "representation" => "complex",
            "orphan_policy" => "aux_psd_free",
        ),
        "options" => solver_options_json(solver, options),
        "runtime" => common_runtime_json(),
    )

    model = nothing
    lowering_error = nothing
    lowering_seconds = @elapsed begin
        try
            model, _ = build_jump_model(
                data.moment_problem;
                formulation = :psd_blocks,
                representation = :complex,
                orphan_policy = :aux_psd_free,
            )
        catch err
            lowering_error = truncated_error(err)
        end
    end
    result["lowering_wall_seconds"] = lowering_seconds
    result["lowering_error"] = lowering_error
    if lowering_error !== nothing
        result["finished_at"] = string(Dates.now())
        result["runtime_after"] = common_runtime_json()
        return result
    end

    cg_iterations = Int[]
    set_solver_error = nothing
    try
        set_solver!(model, solver, options, cg_iterations)
    catch err
        set_solver_error = truncated_error(err)
    end
    result["set_solver_error"] = set_solver_error
    if set_solver_error !== nothing
        result["finished_at"] = string(Dates.now())
        result["runtime_after"] = common_runtime_json()
        return result
    end

    optimize_error = nothing
    optimize_stdio = nothing
    solve_seconds = @elapsed begin
        optimize_error, optimize_stdio = optimize_solver!(model, solver, options; repeat)
    end

    raw = raw_optimizer(model)
    result["solve_wall_seconds_measured"] = solve_seconds
    result["solve_time_sec_reported"] = safe_solve_time(model)
    result["termination_status"] = safe_termination_status(model)
    result["raw_status"] = safe_raw_status(model)
    result["objective_active_Ha"] = safe_objective_value(model)
    result["dual_objective_active_Ha"] = solver == :cosmo_cg ?
        "skipped: MOI.DualObjectiveValue fallback traverses the full H4 constraint matrix" :
        safe_dual_objective_value(model)
    result["optimize_error"] = optimize_error
    optimize_stdio !== nothing && (result["optimize_stdio"] = optimize_stdio)
    result["finished_at"] = string(Dates.now())
    result["runtime_after"] = common_runtime_json()

    if solver == :bpsdp
        result["cg_iterations_per_outer"] = cg_iterations
        result["cg_iterations_total_from_monitor"] = sum(cg_iterations)
        state = bpsdp_state_summary(raw)
        state !== nothing && (result["bpsdp_state"] = state)
    elseif solver == :scs_indirect
        state = scs_state_summary(raw)
        state !== nothing && (result["scs_state"] = state)
    elseif solver == :cosmo_cg
        state = cosmo_state_summary(raw)
        state !== nothing && (result["cosmo_state"] = state)
    end

    return result
end

function constraint_stats(mp)
    hpsd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone == :HPSD]
    psd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone == :PSD]
    zero_sizes = [size(mat) for (cone, mat) in mp.constraints if cone == :Zero]
    total_canonical_moments = length(NCTSSoS._sorted_symmetric_basis(mp.total_basis))
    return (; hpsd_sizes, psd_sizes, zero_sizes, total_canonical_moments,
        direct_real_moment_variables = 2 * total_canonical_moments,
        real_lift_psd_scalar_rows = sum(H4PeriodicBase.real_lift_triangle_rows(n) for n in hpsd_sizes),
        complex_zero_entries = sum((prod(size) for size in zero_sizes); init = 0),
        direct_real_zero_rows = 2 * sum((prod(size) for size in zero_sizes); init = 0))
end

function benchmark_summary(data, build_seconds::Real, options::BenchmarkOptions)
    stats = constraint_stats(data.moment_problem)
    return Dict{String,Any}(
        "problem" => Dict(
            "name" => "H4/Nk=$(data.nk) periodic PQG V2RDM",
            "particle_constraint" => "explicit total-N row dropped; implied by TrD/TrG + identity-augmented ²G",
            "integrals" => options.build.integrals_path,
            "blocking" => string(options.build.blocking),
            "include_one_d" => options.build.include_one_d,
            "spin_resolved_trace" => options.build.spin_resolved_trace,
            "singlet_s2" => options.build.singlet_s2,
            "nk" => data.nk,
            "norb_per_k" => data.norb,
            "spin_orbital_modes" => length(data.spin_orbitals),
            "total_electrons" => data.total_electrons,
            "hf_active_Ha" => data.hf_active,
            "hamiltonian_monomials" => length(H4PeriodicBase.monomials(data.h4_ham)),
            "moment_eq_polynomials" => length(data.meq_constraints),
            "hpsd_block_sizes" => stats.hpsd_sizes,
            "zero_matrix_constraints" => length(stats.zero_sizes),
            "zero_complex_scalar_entries" => stats.complex_zero_entries,
            "total_canonical_moments" => stats.total_canonical_moments,
            "direct_real_moment_variables" => stats.direct_real_moment_variables,
            "real_lift_psd_scalar_rows" => stats.real_lift_psd_scalar_rows,
        ),
        "build_wall_seconds" => build_seconds,
        "output_dir" => options.output_dir,
        "solvers" => string.(options.solvers),
        "repeats" => options.repeats,
        "runtime" => common_runtime_json(),
    )
end

function validate_output_dir!(path::AbstractString)
    abs_path = abspath(path)
    forbidden = Set([
        abspath(@__DIR__),
        abspath(joinpath(@__DIR__, "..")),
        abspath(pwd()),
        homedir(),
        dirname(homedir()),
        "/",
    ])
    abs_path in forbidden && throw(ArgumentError(
        "refusing dangerous --output-dir=$(repr(path)); choose a dedicated benchmark directory"
    ))
    return nothing
end

function prepare_output_dir!(path::AbstractString, solvers::Vector{Symbol})
    validate_output_dir!(path)
    mkpath(path)
    for child in ("problem.json", "benchmark_results.json", "summary.md")
        rm(joinpath(path, child); force = true, recursive = true)
    end
    for solver in solvers
        rm(joinpath(path, solver_label(solver)); force = true, recursive = true)
    end
    return nothing
end

function result_metric(result, keys...)
    x = result
    for k in keys
        x isa AbstractDict || return nothing
        haskey(x, k) || return nothing
        x = x[k]
    end
    return x
end

function write_markdown_summary(path::AbstractString, summary, results)
    mkpath(dirname(path))
    open(path, "w") do io
        p = summary["problem"]
        println(io, "# ", p["name"], " solver benchmark")
        println(io)
        println(io, "Lowering: `formulation=:psd_blocks`, `representation=:complex`, `orphan_policy=:aux_psd_free`.")
        println(io)
        println(io, "## Problem")
        println(io)
        println(io, "| quantity | value |")
        println(io, "|---|---:|")
        for key in ("spin_orbital_modes", "total_electrons", "hamiltonian_monomials", "moment_eq_polynomials", "total_canonical_moments", "direct_real_moment_variables", "real_lift_psd_scalar_rows")
            println(io, "| `", key, "` | ", p[key], " |")
        end
        println(io, "| HPSD block sizes | `", p["hpsd_block_sizes"], "` |")
        println(io)
        println(io, "## Results")
        println(io)
        println(io, "| solver | status | raw status | lower s | solve s | iter | objective | primal res | dual res | gap | error |")
        println(io, "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for r in results
            solver = r["solver"]
            status = get(r, "termination_status", "")
            raw = get(r, "raw_status", "")
            lower = get(r, "lowering_wall_seconds", "")
            solve = get(r, "solve_wall_seconds_measured", "")
            obj = get(r, "objective_active_Ha", "")
            err = something(get(r, "lowering_error", nothing), get(r, "set_solver_error", nothing), get(r, "optimize_error", nothing), "")
            iter = result_metric(r, "bpsdp_state", "outer_iterations")
            iter === nothing && (iter = result_metric(r, "scs_state", "iterations"))
            iter === nothing && (iter = result_metric(r, "cosmo_state", "iterations"))
            prim = result_metric(r, "bpsdp_state", "primal_error")
            prim === nothing && (prim = result_metric(r, "cosmo_state", "info", "r_prim"))
            dual = result_metric(r, "bpsdp_state", "dual_error")
            dual === nothing && (dual = result_metric(r, "cosmo_state", "info", "r_dual"))
            gap = result_metric(r, "bpsdp_state", "objective_gap")
            gap === nothing && (gap = result_metric(r, "scs_state", "objective_gap"))
            println(io, "| `", solver, "` | `", status, "` | `", raw, "` | ", lower, " | ", solve, " | ", something(iter, ""), " | ", obj, " | ", something(prim, ""), " | ", something(dual, ""), " | ", something(gap, ""), " | ", replace(string(err), "\n" => " "), " |")
        end
    end
    return path
end

function print_result_row(result)
    @printf("%-14s repeat=%d lower=%8.3f s solve=%8s s status=%s objective=%s error=%s\n",
        result["solver"], result["repeat"], result["lowering_wall_seconds"],
        string(get(result, "solve_wall_seconds_measured", "")),
        string(get(result, "termination_status", "")),
        string(get(result, "objective_active_Ha", "")),
        string(something(get(result, "lowering_error", nothing), get(result, "set_solver_error", nothing), get(result, "optimize_error", nothing), "")))
    flush(stdout)
end

function main(argv = ARGS)
    options = parse_benchmark_options(argv)
    prepare_output_dir!(options.output_dir, options.solvers)

    @printf("%-48s %s\n", "max RSS at start", rss_string())
    @printf("%-48s %s\n", "Julia threads", string(Threads.nthreads()))
    @printf("%-48s %s\n", "BLAS threads", string(BLAS.get_num_threads()))
    @printf("%-48s %s\n", "output dir", options.output_dir)

    build_seconds = @elapsed data = H4PeriodicBase.build_h4_pqg_moment_problem(h4_options(options.build))
    summary = benchmark_summary(data, build_seconds, options)
    write_json(joinpath(options.output_dir, "problem.json"), summary)

    println("== H4/Nk=$(data.nk) solver benchmark ==")
    @printf("%-48s %.3f s\n", "symbolic build walltime", build_seconds)
    @printf("%-48s %s\n", "HPSD block sizes", string(summary["problem"]["hpsd_block_sizes"]))
    @printf("%-48s %d\n", "direct real moment variables", summary["problem"]["direct_real_moment_variables"])
    @printf("%-48s %d\n", "real-lift PSD scalar rows", summary["problem"]["real_lift_psd_scalar_rows"])
    println()

    all_results = Any[]
    for repeat in 1:options.repeats, solver in options.solvers
        result = solve_one(data, solver, options; repeat)
        push!(all_results, result)
        solver_dir = joinpath(options.output_dir, solver_label(solver), "repeat_$(repeat)")
        write_json(joinpath(solver_dir, "result.json"), result)
        write_json(joinpath(options.output_dir, "benchmark_results.json"), Dict(
            "summary" => summary,
            "results" => all_results,
        ))
        write_markdown_summary(joinpath(options.output_dir, "summary.md"), summary, all_results)
        print_result_row(result)
    end

    println("Wrote benchmark results under ", options.output_dir)
    return (; data, summary, results = all_results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
