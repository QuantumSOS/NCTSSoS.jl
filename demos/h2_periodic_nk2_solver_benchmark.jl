#!/usr/bin/env julia
raw"""
Benchmark BPSDP.jl and SCS.jl on the same H2/Nk=2 periodic PQG V2RDM SDP.

This reuses `h2_periodic_nk2_moment_sos.jl` to build the exact same
MomentProblem: H2, Nk=2, PQG, no explicit total-N constraint.  Each solver gets
a fresh JuMP lowering of that same symbolic problem.

SCS is forced onto its indirect linear solver, i.e. conjugate gradient:

    "linear_solver" => SCS.IndirectSolver

Usage from the repository root:

    julia --project=demos demos/h2_periodic_nk2_solver_benchmark.jl \
        --output-dir=output/phase2/h2_nk2_solver_benchmark

On a fresh solver box:

    julia --project=demos -e 'using Pkg; Pkg.develop(path="."); \
        Pkg.develop(path="/path/to/BPSDP.jl"); Pkg.instantiate()'

Useful knobs:

    --solvers=bpsdp,scs-indirect
    --repeats=3
    --scs-max-iters=200000 --scs-eps-abs=1e-5 --scs-eps-rel=1e-5
    --bpsdp-max-iter=50000 --bpsdp-cg-tol=1e-10
"""

include(joinpath(@__DIR__, "h2_periodic_nk2_moment_sos.jl"))

import SCS

const DEFAULT_BENCHMARK_OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "output", "phase2", "h2_nk2_solver_benchmark"))

struct BenchmarkOptions
    build_options::Options
    solvers::Vector{Symbol}
    repeats::Int
    output_dir::String
    scs_max_iters::Int
    scs_eps_abs::Float64
    scs_eps_rel::Float64
    scs_eps_infeas::Float64
    scs_alpha::Float64
    scs_normalize::Int
    scs_adaptive_scale::Int
    scs_scale::Float64
    scs_rho_x::Float64
    scs_acceleration_lookback::Int
    scs_acceleration_interval::Int
    scs_time_limit_secs::Float64
    scs_verbose::Int
end

function parse_bool_int(value::AbstractString)
    v = lowercase(strip(value))
    v in ("1", "true", "yes", "on") && return 1
    v in ("0", "false", "no", "off") && return 0
    throw(ArgumentError("expected boolean/int 0|1|true|false, got $(repr(value))"))
end

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

function parse_positive_float(name::AbstractString, value::AbstractString)
    return parse_finite_float(name, value; min = 0.0, inclusive_min = false)
end

function parse_nonnegative_float(name::AbstractString, value::AbstractString)
    return parse_finite_float(name, value; min = 0.0)
end

function validate_rank_tols!(rank_tols::Vector{Float64})
    all(t -> isfinite(t) && t > 0.0, rank_tols) ||
        throw(ArgumentError("rank tolerances must be finite positive values"))
    return rank_tols
end

function parse_solvers(value::AbstractString)
    solvers = Symbol[]
    for raw in split(value, ",")
        token = lowercase(strip(raw))
        isempty(token) && continue
        if token == "bpsdp"
            push!(solvers, :bpsdp)
        elseif token in ("scs", "scs-indirect", "scs_indirect")
            push!(solvers, :scs_indirect)
        elseif token == "both"
            append!(solvers, (:bpsdp, :scs_indirect))
        else
            throw(ArgumentError("unknown solver $(repr(token)); expected bpsdp, scs-indirect, or both"))
        end
    end
    isempty(solvers) && throw(ArgumentError("--solvers must name at least one solver"))
    return unique(solvers)
end

function parse_blocking(value::AbstractString)
    if value in ("momentum", "k")
        return :momentum
    elseif value == "spin"
        return :spin
    elseif value == "none"
        return :none
    else
        throw(ArgumentError("unknown --blocking=$value; expected momentum, spin, or none"))
    end
end

function print_benchmark_help()
    println("Usage: julia --project=demos demos/h2_periodic_nk2_solver_benchmark.jl [options]\n")
    println("Build options:")
    println("  --integrals=PATH                  default test/data/assets/h2_chain_nk2_active_2e4o_integrals.txt")
    println("  --blocking=momentum|spin|none     default momentum")
    println("  --include-1d / --no-1d            include/drop extra ¹D PSD block; default --no-1d")
    println("  --paper-spin, --spin-singlet      spin-resolved ²D traces plus singlet S²")
    println("  --spin-resolved-trace / --no-spin-resolved-trace")
    println("  --singlet-s2 / --no-singlet-s2")
    println("  --output-dir=PATH                 default output/phase2/h2_nk2_solver_benchmark")
    println("\nBenchmark options:")
    println("  --solvers=bpsdp,scs-indirect      default bpsdp,scs-indirect")
    println("  --repeats=N                       fresh model per solver per repeat; default 1")
    println("\nBPSDP options:")
    println("  --bpsdp-max-iter=N                default 50000")
    println("  --bpsdp-cg-max-iter=N             default 10000")
    println("  --bpsdp-mu-update-frequency=N     default 500")
    println("  --bpsdp-penalty=rho               default 0.1")
    println("  --bpsdp-cg-tol=eps                default 1e-10")
    println("  --bpsdp-obj-tol=eps               default 1e-7")
    println("  --bpsdp-err-tol=eps               default 1e-7")
    println("  --bpsdp-print-level=N             default 1")
    println("\nSCS indirect/CG options:")
    println("  --scs-max-iters=N                 default 100000")
    println("  --scs-eps-abs=eps                 default 1e-5")
    println("  --scs-eps-rel=eps                 default 1e-5")
    println("  --scs-eps-infeas=eps              default 1e-7")
    println("  --scs-alpha=a                     default 1.5")
    println("  --scs-normalize=0|1               default 1")
    println("  --scs-adaptive-scale=0|1          default 1")
    println("  --scs-scale=s                     default 0.1")
    println("  --scs-rho-x=rho                   default 1e-6")
    println("  --scs-acceleration-lookback=N     default 10; use 0 to disable")
    println("  --scs-acceleration-interval=N     default 10")
    println("  --scs-time-limit-secs=sec         default 0.0, no limit")
    println("  --scs-verbose=0|1                 default 1")
end

function parse_benchmark_options(argv)
    integrals_path = DEFAULT_INTEGRALS
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = false
    singlet_s2 = false
    output_dir = DEFAULT_BENCHMARK_OUTPUT_DIR
    rank_tols = copy(DEFAULT_RANK_TOLS)
    diagnostic_block = ""

    bpsdp_max_iter = 50_000
    bpsdp_cg_max_iter = 10_000
    bpsdp_mu_update_frequency = 500
    bpsdp_penalty = 0.1
    bpsdp_cg_tol = 1e-10
    bpsdp_obj_tol = 1e-7
    bpsdp_err_tol = 1e-7
    bpsdp_print_level = 1

    solvers = [:bpsdp, :scs_indirect]
    repeats = 1

    scs_max_iters = 100_000
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

    for arg in argv
        if startswith(arg, "--integrals=")
            integrals_path = split(arg, "=", limit = 2)[2]
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
        elseif startswith(arg, "--rank-tol=")
            rank_tols = validate_rank_tols!(parse_rank_tols(split(arg, "=", limit = 2)[2]))
        elseif startswith(arg, "--diagnostic-block=")
            diagnostic_block = split(arg, "=", limit = 2)[2]
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
        elseif arg in ("-h", "--help")
            print_benchmark_help()
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    build_options = Options(integrals_path, blocking, include_one_d,
        spin_resolved_trace, singlet_s2, output_dir, rank_tols, diagnostic_block,
        bpsdp_max_iter, bpsdp_cg_max_iter, bpsdp_mu_update_frequency,
        bpsdp_penalty, bpsdp_cg_tol, bpsdp_obj_tol, bpsdp_err_tol,
        bpsdp_print_level)

    return BenchmarkOptions(build_options, solvers, repeats, output_dir,
        scs_max_iters, scs_eps_abs, scs_eps_rel, scs_eps_infeas, scs_alpha,
        scs_normalize, scs_adaptive_scale, scs_scale, scs_rho_x,
        scs_acceleration_lookback, scs_acceleration_interval,
        scs_time_limit_secs, scs_verbose)
end

function scs_indirect_optimizer_factory(options::BenchmarkOptions)
    SCS.is_available(SCS.IndirectSolver) || error("SCS.IndirectSolver is not available in this SCS.jl build")
    return optimizer_with_attributes(
        SCS.Optimizer,
        "linear_solver" => SCS.IndirectSolver,
        "max_iters" => options.scs_max_iters,
        "eps_abs" => options.scs_eps_abs,
        "eps_rel" => options.scs_eps_rel,
        "eps_infeas" => options.scs_eps_infeas,
        "alpha" => options.scs_alpha,
        "normalize" => options.scs_normalize,
        "adaptive_scale" => options.scs_adaptive_scale,
        "scale" => options.scs_scale,
        "rho_x" => options.scs_rho_x,
        "acceleration_lookback" => options.scs_acceleration_lookback,
        "acceleration_interval" => options.scs_acceleration_interval,
        "time_limit_secs" => options.scs_time_limit_secs,
        "verbose" => options.scs_verbose,
    )
end

function safe_termination_status(model)
    try
        return string(termination_status(model))
    catch err
        return string("unavailable: ", sprint(showerror, err))
    end
end

function safe_raw_status(model)
    try
        return string(MOI.get(model, MOI.RawStatusString()))
    catch err
        return string("unavailable: ", sprint(showerror, err))
    end
end

function safe_objective_value(model)
    try
        return objective_value(model)
    catch err
        return string("unavailable: ", sprint(showerror, err))
    end
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
        "raw_status" => string(field_or_nothing(sol, :raw_status)),
        "ret_val" => field_or_nothing(sol, :ret_val),
    )
end

function solver_options_json(solver::Symbol, options::BenchmarkOptions)
    if solver == :bpsdp
        b = options.build_options
        return Dict{String,Any}(
            "max_iter" => b.bpsdp_max_iter,
            "cg_max_iter" => b.bpsdp_cg_max_iter,
            "mu_update_frequency" => b.bpsdp_mu_update_frequency,
            "penalty_parameter" => b.bpsdp_penalty,
            "cg_convergence" => b.bpsdp_cg_tol,
            "dynamic_cg_convergence" => true,
            "sdp_objective_convergence" => b.bpsdp_obj_tol,
            "sdp_error_convergence" => b.bpsdp_err_tol,
            "guess_type" => "zero",
            "dependent_rows" => "keep",
        )
    elseif solver == :scs_indirect
        return Dict{String,Any}(
            "linear_solver" => "SCS.IndirectSolver (conjugate gradient)",
            "max_iters" => options.scs_max_iters,
            "eps_abs" => options.scs_eps_abs,
            "eps_rel" => options.scs_eps_rel,
            "eps_infeas" => options.scs_eps_infeas,
            "alpha" => options.scs_alpha,
            "normalize" => options.scs_normalize,
            "adaptive_scale" => options.scs_adaptive_scale,
            "scale" => options.scs_scale,
            "rho_x" => options.scs_rho_x,
            "acceleration_lookback" => options.scs_acceleration_lookback,
            "acceleration_interval" => options.scs_acceleration_interval,
            "time_limit_secs" => options.scs_time_limit_secs,
            "verbose" => options.scs_verbose,
        )
    else
        error("unknown solver $solver")
    end
end

function set_solver!(model, solver::Symbol, options::BenchmarkOptions, cg_iterations::Vector{Int})
    if solver == :bpsdp
        set_optimizer(model, bpsdp_optimizer_factory(options.build_options, cg_iterations))
    elseif solver == :scs_indirect
        set_optimizer(model, scs_indirect_optimizer_factory(options))
    else
        error("unknown solver $solver")
    end
    return nothing
end

function solve_one(data, solver::Symbol, options::BenchmarkOptions; repeat::Int)
    model = nothing
    extract_monomap = nothing
    lowering_seconds = @elapsed begin
        model, extract_monomap = build_jump_model(
            data.moment_problem;
            formulation = :psd_blocks,
            representation = :complex,
            orphan_policy = :aux_psd_free,
        )
    end

    cg_iterations = Int[]
    set_solver!(model, solver, options, cg_iterations)

    optimize_error = nothing
    solve_seconds = @elapsed begin
        try
            optimize!(model)
        catch err
            optimize_error = sprint(showerror, err)
        end
    end

    raw = raw_optimizer(model)
    result = Dict{String,Any}(
        "solver" => string(solver),
        "repeat" => repeat,
        "lowering_wall_seconds" => lowering_seconds,
        "solve_wall_seconds_measured" => solve_seconds,
        "termination_status" => safe_termination_status(model),
        "raw_status" => safe_raw_status(model),
        "objective_active_Ha" => safe_objective_value(model),
        "optimize_error" => optimize_error,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "jump_lowering" => Dict(
            "formulation" => "psd_blocks",
            "representation" => "complex",
            "orphan_policy" => "aux_psd_free",
        ),
        "options" => solver_options_json(solver, options),
    )

    if solver == :bpsdp
        result["cg_iterations_per_outer"] = cg_iterations
        result["cg_iterations_total_from_monitor"] = sum(cg_iterations)
        state = bpsdp_state_summary(raw)
        state !== nothing && (result["bpsdp_state"] = state)
    elseif solver == :scs_indirect
        state = scs_state_summary(raw)
        state !== nothing && (result["scs_state"] = state)
    end

    return result
end

function solver_label(solver::Symbol)
    solver == :scs_indirect && return "scs_indirect"
    return string(solver)
end

function benchmark_summary(data, build_seconds::Real, options::BenchmarkOptions)
    stats = constraint_stats(data.moment_problem)
    return Dict{String,Any}(
        "problem" => Dict(
            "name" => "H2/Nk=2 periodic PQG V2RDM",
            "particle_constraint" => "explicit total-N row dropped; implied by TrD/TrG + identity-augmented ²G",
            "integrals" => options.build_options.integrals_path,
            "blocking" => string(options.build_options.blocking),
            "include_one_d" => options.build_options.include_one_d,
            "spin_resolved_trace" => options.build_options.spin_resolved_trace,
            "singlet_s2" => options.build_options.singlet_s2,
            "nk" => data.nk,
            "norb_per_k" => data.norb,
            "spin_orbital_modes" => length(data.spin_orbitals),
            "total_electrons" => data.total_electrons,
            "moment_eq_polynomials" => length(data.meq_constraints),
            "hpsd_block_sizes" => stats.hpsd_sizes,
            "zero_matrix_constraints" => length(stats.zero_sizes),
            "total_canonical_moments" => stats.total_canonical_moments,
            "direct_real_moment_variables" => stats.direct_real_moment_variables,
            "real_lift_psd_scalar_rows" => stats.real_lift_psd_rows,
        ),
        "build_wall_seconds" => build_seconds,
        "output_dir" => options.output_dir,
        "solvers" => string.(options.solvers),
        "repeats" => options.repeats,
    )
end

function print_result_row(result)
    @printf("%-14s repeat=%d lower=%8.3f s solve=%8.3f s status=%s objective=%s\n",
        result["solver"], result["repeat"], result["lowering_wall_seconds"],
        result["solve_wall_seconds_measured"], result["termination_status"],
        string(result["objective_active_Ha"]))
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

function prepare_output_dir!(path::AbstractString)
    validate_output_dir!(path)
    mkpath(path)
    for child in ("problem.json", "benchmark_results.json", "bpsdp", "scs_indirect")
        rm(joinpath(path, child); force = true, recursive = true)
    end
    return nothing
end

function main(argv = ARGS)
    options = parse_benchmark_options(argv)
    prepare_output_dir!(options.output_dir)

    @printf("%-48s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h2_pqg_moment_problem(options.build_options)
    summary = benchmark_summary(data, build_seconds, options)
    write_json(joinpath(options.output_dir, "problem.json"), summary)

    println("== H2/Nk=2 solver benchmark ==")
    @printf("%-48s %.3f s\n", "symbolic build walltime", build_seconds)
    @printf("%-48s %s\n", "SCS linear solver", "SCS.IndirectSolver (conjugate gradient)")
    @printf("%-48s %s\n", "output dir", options.output_dir)
    println()

    all_results = Any[]
    for repeat in 1:options.repeats, solver in options.solvers
        result = solve_one(data, solver, options; repeat)
        push!(all_results, result)
        solver_dir = joinpath(options.output_dir, solver_label(solver), "repeat_$(repeat)")
        write_json(joinpath(solver_dir, "result.json"), result)
        print_result_row(result)
        flush(stdout)
    end

    write_json(joinpath(options.output_dir, "benchmark_results.json"), Dict(
        "summary" => summary,
        "results" => all_results,
    ))
    println("Wrote benchmark results under ", options.output_dir)

    thrown = [r for r in all_results if r["optimize_error"] !== nothing]
    isempty(thrown) || error("$(length(thrown)) optimize! call(s) threw; see benchmark_results.json")

    return (; data, summary, results = all_results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
