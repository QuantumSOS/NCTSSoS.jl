#!/usr/bin/env julia
# Sequential H4/Nk PQG smoke/scale sweep through the live solver path:
# NCTSSoS MomentProblem -> JuMP -> MOI bridges/copy_to -> BPSDP.jl.
#
# This intentionally does not run the Phase-2 A*A' diagnostics and does not
# export libsdp/dat-c files. It measures build/lowering/solve sizes and how far
# BPSDP gets under a fixed iteration/tolerance budget.

using Dates
using JSON3
using JuMP
const MOI = JuMP.MOI
using LinearAlgebra
using Printf
using SparseArrays

try
    @eval import BPSDP
catch err
    error("BPSDP.jl is required. Run on HAI with Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\"). Original error: $(err)")
end

include(joinpath(@__DIR__, "h4_nk_periodic_moment_sos_export_libsdp.jl"))

struct SweepOptions
    nks::Vector{Int}
    norb::Int
    nelec_per_cell::Int
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
    integrals_root::String
    integrals_drop_label::String
    outdir::String
    max_iter::Int
    cg_max_iter::Int
    atol::Float64
    penalty::Float64
    mu_update_frequency::Int
    progress_stride::Int
end

function parse_int_list(value::AbstractString)
    nks = Int[]
    for tok in split(value, ',')
        s = strip(tok)
        isempty(s) && continue
        push!(nks, parse(Int, s))
    end
    isempty(nks) && throw(ArgumentError("--nks must contain at least one integer"))
    return nks
end

function parse_sweep_options(argv)
    nks = [2, 3, 4]
    norb = 8
    nelec_per_cell = 4
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = true
    singlet_s2 = true
    integrals_root = ""
    integrals_drop_label = "1e-12"
    outdir = normpath(joinpath(@__DIR__, "..", "output", "h4_nk_jump_bpsdp_atol1e-4_iter2000"))
    max_iter = 2_000
    cg_max_iter = 2_000
    atol = 1e-4
    penalty = 0.1
    mu_update_frequency = 25
    progress_stride = 100

    for arg in argv
        if startswith(arg, "--nks=")
            nks = parse_int_list(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--norb=")
            norb = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--nelec-per-cell=")
            nelec_per_cell = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--blocking=")
            value = split(arg, "=", limit = 2)[2]
            if value in ("momentum", "k")
                blocking = :momentum
            elseif value == "spin"
                blocking = :spin
            elseif value == "none"
                blocking = :none
            else
                throw(ArgumentError("unknown --blocking=$value; expected momentum, spin, or none"))
            end
        elseif arg == "--include-1d"
            include_one_d = true
        elseif arg == "--no-1d"
            include_one_d = false
        elseif arg in ("--paper-spin", "--spin-singlet")
            spin_resolved_trace = true
            singlet_s2 = true
        elseif arg == "--no-paper-spin"
            spin_resolved_trace = false
            singlet_s2 = false
        elseif arg == "--spin-resolved-trace"
            spin_resolved_trace = true
        elseif arg == "--no-spin-resolved-trace"
            spin_resolved_trace = false
        elseif arg == "--singlet-s2"
            singlet_s2 = true
        elseif arg == "--no-singlet-s2"
            singlet_s2 = false
        elseif startswith(arg, "--integrals-root=")
            integrals_root = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--integrals-drop-label=")
            integrals_drop_label = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--max-iter=")
            max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--cg-max-iter=")
            cg_max_iter = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--atol=")
            atol = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--penalty=")
            penalty = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--mu-update-frequency=")
            mu_update_frequency = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--progress-stride=")
            progress_stride = parse(Int, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=<env-with-BPSDP> probes/h4_nk_jump_bpsdp_sweep.jl [options]")
            println("  --nks=2,3,4                 default 2,3,4")
            println("  --norb=N                    default 8 (H4 [4,8]-style active space per k)")
            println("  --nelec-per-cell=N          default 4")
            println("  --blocking=momentum|spin|none  default momentum")
            println("  --paper-spin / --no-paper-spin  default paper-spin")
            println("  --include-1d / --no-1d      default no-1d")
            println("  --integrals-root=PATH       physical integral tree with nkN/integrals files; default mock integrals")
            println("  --integrals-drop-label=STR  filename label in h4_chain_nkN_integrals_drop<STR>.txt; default 1e-12")
            println("  --outdir=PATH")
            println("  --max-iter=N                BPSDP outer iteration cap; default 2000")
            println("  --cg-max-iter=N             BPSDP CG cap; default 2000")
            println("  --atol=FLOAT                BPSDP objective/error/CG tolerance; default 1e-4")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    all(>=(1), nks) || throw(ArgumentError("all Nk values must be positive"))
    norb >= 1 || throw(ArgumentError("--norb must be positive"))
    nelec_per_cell >= 1 || throw(ArgumentError("--nelec-per-cell must be positive"))
    max_iter >= 0 || throw(ArgumentError("--max-iter must be nonnegative"))
    cg_max_iter >= 1 || throw(ArgumentError("--cg-max-iter must be positive"))
    atol > 0 || throw(ArgumentError("--atol must be positive"))
    progress_stride >= 1 || throw(ArgumentError("--progress-stride must be positive"))
    if !isempty(integrals_root)
        isdir(integrals_root) || throw(ArgumentError("--integrals-root does not exist: $integrals_root"))
    end

    return SweepOptions(nks, norb, nelec_per_cell, blocking, include_one_d,
        spin_resolved_trace, singlet_s2, integrals_root, integrals_drop_label,
        outdir, max_iter, cg_max_iter, atol, penalty, mu_update_frequency,
        progress_stride)
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, obj)
        write(io, '\n')
    end
    return path
end

function append_jsonl(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "a") do io
        JSON3.write(io, obj)
        write(io, '\n')
    end
    return path
end

has_field(obj, name::Symbol) = name in fieldnames(typeof(obj))
field_or_nothing(obj, name::Symbol) = has_field(obj, name) ? getfield(obj, name) : nothing
function first_field_or_nothing(obj, names::Symbol...)
    for name in names
        has_field(obj, name) && return getfield(obj, name)
    end
    return nothing
end

function finite_or_string(x)
    x === nothing && return nothing
    y = try Float64(x) catch; return string(x) end
    return isfinite(y) ? y : string(y)
end

function objective_gap_from(primal, dual)
    (primal === nothing || dual === nothing) && return nothing
    gap = try
        abs(Float64(primal) - Float64(dual))
    catch
        return nothing
    end
    return finite_or_string(gap)
end

function raw_optimizer(model::JuMP.Model)
    try
        return JuMP.unsafe_backend(model)
    catch
        return nothing
    end
end

function bpsdp_state_summary(raw)
    state = raw === nothing ? nothing : field_or_nothing(raw, :state)
    state === nothing && return nothing
    return Dict{String,Any}(
        "outer_iterations" => first_field_or_nothing(state, :outer_iterations, :iterations),
        "inner_iterations" => field_or_nothing(state, :inner_iterations),
        "primal_error" => finite_or_string(field_or_nothing(state, :primal_error)),
        "dual_error" => finite_or_string(field_or_nothing(state, :dual_error)),
        "objective_primal" => finite_or_string(first_field_or_nothing(state, :objective_primal, :primal_objective)),
        "objective_dual" => finite_or_string(first_field_or_nothing(state, :objective_dual, :dual_objective)),
        "objective_gap" => finite_or_string(field_or_nothing(state, :objective_gap)),
        "termination_reason" => string(first_field_or_nothing(state, :termination_reason, :status)),
    )
end

function block_origin_string(block)
    origin = block.meta.origin
    return sprint(show, origin)
end

function moment_problem_size(data)
    L = data.moment_problem.linear
    block_dims = [block.size for block in L.psd_blocks_lin]
    origins = [block_origin_string(block) for block in L.psd_blocks_lin]
    return Dict{String,Any}(
        "n_moments" => length(L.moments),
        "n_free_keys" => length(L.free_keys),
        "n_pivots" => length(L.pivots),
        "n_zero_scalar_rows" => length(L.zero_constraints),
        "objective_terms" => length(L.objective_lin),
        "hpsd_blocks" => length(block_dims),
        "hpsd_block_dims" => block_dims,
        "hpsd_block_origins" => origins,
        "largest_hpsd_block" => isempty(block_dims) ? 0 : maximum(block_dims),
        "sum_hpsd_n2" => sum(n -> n * n, block_dims; init = 0),
        "n_unique_moment_matrix_elements" => data.moment_problem.n_unique_moment_matrix_elements,
        "hamiltonian_monomials" => length(monomials(data.h4_ham)),
        "moment_eq_polynomials" => length(data.meq_constraints),
    )
end

key_string(key) = sprint(show, key)

function sorted_key_sample(keys; limit::Int = 12)
    strs = sort!(key_string.(collect(keys)))
    return strs[1:min(length(strs), limit)]
end

function pivot_orphan_audit(mp)
    L = mp.linear
    n_moments = length(L.moments)
    n_orphans = length(L.free_keys)
    return Dict{String,Any}(
        "n_moments" => n_moments,
        "n_pivots" => length(L.pivots),
        "n_orphans" => n_orphans,
        "orphan_pct" => n_moments == 0 ? 0.0 : 100.0 * n_orphans / n_moments,
        "n_zero_scalar_rows" => length(L.zero_constraints),
        "sample_orphans" => sorted_key_sample(L.free_keys),
    )
end

function bpsdp_problem_size(raw)
    raw === nothing && return nothing
    dims = has_field(raw, :block_dims) ? getfield(raw, :block_dims) : Int[]
    kinds = has_field(raw, :block_kinds) ? getfield(raw, :block_kinds) : Symbol[]
    scalar_real = (!isempty(dims) && !isempty(kinds)) ? count(i -> dims[i] == 1 && kinds[i] == :real, eachindex(dims)) : count(==(1), dims)
    return Dict{String,Any}(
        "n_blocks" => length(dims),
        "block_dims" => collect(dims),
        "block_kinds" => string.(collect(kinds)),
        "scalar_1x1_blocks" => scalar_real,
        "largest_block_dim" => isempty(dims) ? 0 : maximum(dims),
        "primal_dim" => has_field(raw, :c) ? length(getfield(raw, :c)) : nothing,
        "dual_rows" => has_field(raw, :b) ? length(getfield(raw, :b)) : nothing,
        "A_nnz" => has_field(raw, :A) ? nnz(getfield(raw, :A)) : nothing,
        "A_shape" => has_field(raw, :A) ? [size(getfield(raw, :A), 1), size(getfield(raw, :A), 2)] : nothing,
        "variable_map_length" => has_field(raw, :variable_map) ? length(getfield(raw, :variable_map)) : nothing,
        "reals_count" => has_field(raw, :reals_count) ? getfield(raw, :reals_count) : nothing,
    )
end

function optimizer_factory(opts::SweepOptions, histories)
    return () -> begin
        cg_total = Ref(0)
        return BPSDP.Optimizer(
            max_iter = opts.max_iter,
            cg_max_iter = opts.cg_max_iter,
            mu_update_frequency = opts.mu_update_frequency,
            penalty_parameter = opts.penalty,
            cg_convergence = opts.atol,
            dynamic_cg_convergence = true,
            sdp_objective_convergence = opts.atol,
            sdp_error_convergence = opts.atol,
            guess_type = :zero,
            print_level = opts.progress_stride,
            dependent_rows = :keep,
            progress_monitor = function (print_level, outer, inner, primal, dual, mu, perr, derr)
                cg_total[] += Int(inner)
                if outer == 0 || outer % opts.progress_stride == 0
                    gap = objective_gap_from(primal, dual)
                    row = Dict{String,Any}(
                        "time" => string(now()),
                        "outer" => Int(outer),
                        "inner" => Int(inner), # backward-compatible alias: CG iterations this outer step
                        "cg_iterations" => Int(inner),
                        "cg_iterations_total" => cg_total[],
                        "primal" => finite_or_string(primal),
                        "dual" => finite_or_string(dual),
                        "objective_primal" => finite_or_string(primal),
                        "objective_dual" => finite_or_string(dual),
                        "objective_gap" => gap,
                        "mu" => finite_or_string(mu),
                        "primal_error" => finite_or_string(perr),
                        "dual_error" => finite_or_string(derr),
                    )
                    append_jsonl(histories.progress, row)
                    append_jsonl(histories.objective, Dict(
                        "time" => row["time"],
                        "outer" => row["outer"],
                        "objective_primal" => row["objective_primal"],
                        "objective_dual" => row["objective_dual"],
                        "objective_gap" => row["objective_gap"],
                    ))
                    append_jsonl(histories.residuals, Dict(
                        "time" => row["time"],
                        "outer" => row["outer"],
                        "primal_error" => row["primal_error"],
                        "dual_error" => row["dual_error"],
                        "objective_gap" => row["objective_gap"],
                        "mu" => row["mu"],
                    ))
                    append_jsonl(histories.cg, Dict(
                        "time" => row["time"],
                        "outer" => row["outer"],
                        "cg_iterations" => row["cg_iterations"],
                        "cg_iterations_total" => row["cg_iterations_total"],
                    ))
                end
                return nothing
            end,
        )
    end
end

function physical_integrals_path(opts::SweepOptions, nk::Int)
    isempty(opts.integrals_root) && return nothing
    path = joinpath(opts.integrals_root, "nk$(nk)", "integrals",
        "h4_chain_nk$(nk)_integrals_drop$(opts.integrals_drop_label).txt")
    isfile(path) || throw(ArgumentError("missing physical integrals for Nk=$nk: $path"))
    return path
end

function physical_meta_path(opts::SweepOptions, nk::Int)
    isempty(opts.integrals_root) && return nothing
    path = joinpath(opts.integrals_root, "nk$(nk)", "integrals", "h4_chain_nk$(nk)_meta.json")
    return isfile(path) ? path : nothing
end

function read_physical_meta(opts::SweepOptions, nk::Int)
    path = physical_meta_path(opts, nk)
    path === nothing && return nothing
    return JSON3.read(read(path, String), Dict{String,Any})
end

function energy_shift_from_meta(meta)
    meta === nothing && return nothing
    haskey(meta, "energy_shift") || return nothing
    return Float64(meta["energy_shift"])
end

function adjusted_objectives(result)
    shift = get(result, "energy_shift", nothing)
    state = get(result, "bpsdp_state", nothing)
    (shift === nothing || state === nothing) && return nothing
    primal = get(state, "objective_primal", nothing)
    dual = get(state, "objective_dual", nothing)
    return Dict{String,Any}(
        "energy_shift" => shift,
        "adjusted_primal" => primal === nothing ? nothing : Float64(primal) + Float64(shift),
        "adjusted_dual" => dual === nothing ? nothing : Float64(dual) + Float64(shift),
    )
end

function build_options_for_nk(opts::SweepOptions, nk::Int, nk_dir::AbstractString)
    return Options(
        nk,
        opts.norb,
        opts.nelec_per_cell,
        physical_integrals_path(opts, nk),
        opts.blocking,
        opts.include_one_d,
        opts.spin_resolved_trace,
        opts.singlet_s2,
        nk_dir,
        "h4_nk$(nk)_jump_bpsdp",
        :min,
        32,
    )
end

function run_one(opts::SweepOptions, nk::Int)
    nk_dir = joinpath(opts.outdir, "nk$(nk)")
    mkpath(nk_dir)
    result_path = joinpath(nk_dir, "result.json")
    histories = (
        progress = joinpath(nk_dir, "progress.jsonl"),
        objective = joinpath(nk_dir, "objective_history.jsonl"),
        residuals = joinpath(nk_dir, "residual_gap_history.jsonl"),
        cg = joinpath(nk_dir, "cg_history.jsonl"),
    )

    result = Dict{String,Any}(
        "nk" => nk,
        "stage" => "starting",
        "started_at" => string(now()),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "rss_start" => rss_string(),
        "formulation" => Dict(
            "route" => "NCTSSoS MomentProblem -> JuMP build_jump_model(:psd_blocks) -> MOI/BPSDP",
            "integrals" => "deterministic mock_h4_integrals",
            "norb" => opts.norb,
            "nelec_per_cell" => opts.nelec_per_cell,
            "blocking" => string(opts.blocking),
            "include_one_d" => opts.include_one_d,
            "spin_resolved_trace" => opts.spin_resolved_trace,
            "singlet_s2" => opts.singlet_s2,
            "explicit_total_particle_constraint" => false,
            "particle_constraint_note" => "N̂-N dropped; electron count retained through TrD/TrQ/TrG, spin-resolved traces, and S²",
        ),
        "bpsdp_options" => Dict(
            "max_iter" => opts.max_iter,
            "cg_max_iter" => opts.cg_max_iter,
            "atol_for_cg_objective_error" => opts.atol,
            "penalty_parameter" => opts.penalty,
            "mu_update_frequency" => opts.mu_update_frequency,
            "dependent_rows" => "keep",
        ),
    )
    write_json(result_path, result)

    try
        bopts = build_options_for_nk(opts, nk, nk_dir)
        meta = read_physical_meta(opts, nk)
        formulation = result["formulation"]
        formulation["integrals"] = bopts.integrals_path === nothing ? "deterministic mock_h4_integrals" : bopts.integrals_path
        meta_path = physical_meta_path(opts, nk)
        meta_path === nothing || (formulation["integrals_meta"] = meta_path)
        if meta !== nothing
            result["physical_meta"] = meta
            shift = energy_shift_from_meta(meta)
            shift === nothing || (result["energy_shift"] = shift)
            haskey(meta, "ehf") && (result["hf_total_reference"] = Float64(meta["ehf"]))
        end

        result["stage"] = "building_moment_problem"
        write_json(result_path, result)
        @printf("== H4/Nk=%d build ==\n", nk); flush(stdout)
        result["build_seconds"] = @elapsed data = build_h4_pqg_moment_problem(bopts)
        result["hf_active"] = data.hf_active
        result["moment_problem"] = moment_problem_size(data)
        audit = pivot_orphan_audit(data.moment_problem)
        result["pivot_orphan_audit"] = audit
        write_json(joinpath(nk_dir, "pivot_orphan_audit.json"), audit)
        result["rss_after_build"] = rss_string()
        result["stage"] = "building_jump_model"
        write_json(result_path, result)

        @printf("== H4/Nk=%d JuMP lowering ==\n", nk); flush(stdout)
        result["jump_build_seconds"] = @elapsed begin
            model, extract_monomap = build_jump_model(
                data.moment_problem;
                formulation = :psd_blocks,
                representation = :complex,
                orphan_policy = :aux_psd_free,
            )
        end
        result["jump"] = Dict{String,Any}(
            "num_variables" => JuMP.num_variables(model),
            "num_constraints_total" => JuMP.num_constraints(model; count_variable_in_set_constraints = true),
        )
        result["rss_after_jump_build"] = rss_string()
        result["stage"] = "optimizing_bpsdp"
        write_json(result_path, result)

        set_optimizer(model, optimizer_factory(opts, histories))
        optimize_error = nothing
        @printf("== H4/Nk=%d BPSDP optimize max_iter=%d atol=%.1e ==\n", nk, opts.max_iter, opts.atol)
        flush(stdout)
        result["optimize_seconds"] = @elapsed begin
            try
                optimize!(model)
            catch err
                optimize_error = sprint(showerror, err, catch_backtrace())
            end
        end
        result["optimize_error"] = optimize_error
        result["termination_status"] = try string(termination_status(model)) catch err string("unavailable: ", err) end
        result["raw_status"] = try string(MOI.get(model, MOI.RawStatusString())) catch err string("unavailable: ", err) end
        result["objective_value"] = try finite_or_string(objective_value(model)) catch err string("unavailable: ", err) end

        raw = raw_optimizer(model)
        dimensions = bpsdp_problem_size(raw)
        result["bpsdp_problem"] = dimensions
        dimensions !== nothing && write_json(joinpath(nk_dir, "bpsdp_dimensions.json"), dimensions)
        state_summary = bpsdp_state_summary(raw)
        state_summary !== nothing && (result["bpsdp_state"] = state_summary)
        result["rss_after_optimize"] = rss_string()
        result["stage"] = "finished"
    catch err
        result["stage"] = "error"
        result["error"] = sprint(showerror, err, catch_backtrace())
        result["rss_after_error"] = rss_string()
    finally
        result["finished_at"] = string(now())
        adjusted = adjusted_objectives(result)
        adjusted === nothing || (result["adjusted_objectives"] = adjusted)
        status_summary = Dict{String,Any}(
            "nk" => nk,
            "stage" => get(result, "stage", "unknown"),
            "termination_status" => get(result, "termination_status", nothing),
            "raw_status" => get(result, "raw_status", nothing),
            "optimize_error" => get(result, "optimize_error", nothing),
            "started_at" => get(result, "started_at", nothing),
            "finished_at" => get(result, "finished_at", nothing),
            "result_path" => result_path,
            "progress_path" => histories.progress,
        )
        objective_summary = Dict{String,Any}(
            "nk" => nk,
            "objective_value" => get(result, "objective_value", nothing),
            "termination_status" => get(result, "termination_status", nothing),
            "raw_status" => get(result, "raw_status", nothing),
            "energy_shift" => get(result, "energy_shift", nothing),
            "adjusted_objectives" => get(result, "adjusted_objectives", nothing),
            "bpsdp_state" => get(result, "bpsdp_state", nothing),
        )
        write_json(joinpath(nk_dir, "status.json"), status_summary)
        write_json(joinpath(nk_dir, "objective.json"), objective_summary)
        write_json(result_path, result)
        append_jsonl(joinpath(opts.outdir, "summary.jsonl"), result)
        @printf("== H4/Nk=%d done: %s ==\n", nk, get(result, "stage", "unknown"))
        @printf("result: %s\n", result_path)
        flush(stdout)
    end

    return result
end

function write_run_metadata(opts::SweepOptions)
    mkpath(opts.outdir)
    return write_json(joinpath(opts.outdir, "run_metadata.json"), Dict{String,Any}(
        "started_at" => string(now()),
        "argv" => ARGS,
        "nks" => opts.nks,
        "norb" => opts.norb,
        "nelec_per_cell" => opts.nelec_per_cell,
        "blocking" => string(opts.blocking),
        "include_one_d" => opts.include_one_d,
        "spin_resolved_trace" => opts.spin_resolved_trace,
        "singlet_s2" => opts.singlet_s2,
        "integrals_root" => opts.integrals_root,
        "integrals_drop_label" => opts.integrals_drop_label,
        "outdir" => opts.outdir,
        "max_iter" => opts.max_iter,
        "cg_max_iter" => opts.cg_max_iter,
        "atol" => opts.atol,
        "note" => isempty(opts.integrals_root) ?
            "Sequential no-AAstar solver-path sweep. Explicit N̂-N constraint is dropped; deterministic mock_h4_integrals are used." :
            "Sequential no-AAstar solver-path sweep. Explicit N̂-N constraint is dropped; physical PySCF integrals and per-Nk energy shifts are loaded from integrals_root.",
    ))
end

function main(argv = ARGS)
    opts = parse_sweep_options(argv)
    write_run_metadata(opts)
    println("== H4 Nk JuMP/MOI/BPSDP sequential sweep ==")
    println("output: ", opts.outdir)
    println("Nk: ", join(opts.nks, ","), " norb=", opts.norb, " nelec/cell=", opts.nelec_per_cell)
    println("integrals: ", isempty(opts.integrals_root) ? "deterministic mock_h4_integrals" : opts.integrals_root)
    println("BPSDP max_iter=", opts.max_iter, " cg_max_iter=", opts.cg_max_iter, " atol=", opts.atol)
    println("No A*A' diagnostics. No libsdp/dat-c export.")
    flush(stdout)
    return [run_one(opts, nk) for nk in opts.nks]
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
