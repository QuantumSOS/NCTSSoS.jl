#!/usr/bin/env julia
raw"""
Aggregate H4 sector-pruning run artifacts into a markdown summary.

Usage from the repository root:

    julia --project=<env-with-JSON3> probes/h4_nk_pruning_diagnostic.jl \
        --input-dir=output/h4_pruning --output=output/h4_pruning/summary.md

The input directories are produced by
`demos/h4_periodic_moment_sos_sector_pruned.jl` and contain `report.json`,
`cones.json`, and `kernel.json`.
"""

using JSON3
using Printf
using Dates

const DEFAULT_INPUT_DIR = normpath(joinpath(@__DIR__, "..", "output", "h4_pruning"))
const DEFAULT_OUTPUT = normpath(joinpath(DEFAULT_INPUT_DIR, "summary.md"))

struct ProbeOptions
    input_dir::String
    output_path::String
end

function print_help()
    println("Usage: julia --project=<env-with-JSON3> probes/h4_nk_pruning_diagnostic.jl [options]\n")
    println("Options:")
    println("  --input-dir=PATH       default output/h4_pruning")
    println("  --output=PATH          default output/h4_pruning/summary.md")
    return nothing
end

function parse_options(argv)
    input_dir = DEFAULT_INPUT_DIR
    output_path = DEFAULT_OUTPUT
    for arg in argv
        if startswith(arg, "--input-dir=")
            input_dir = split(arg, "="; limit = 2)[2]
        elseif startswith(arg, "--output=")
            output_path = split(arg, "="; limit = 2)[2]
        elseif arg in ("-h", "--help")
            print_help()
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    return ProbeOptions(input_dir, output_path)
end

function read_json(path::AbstractString)
    return JSON3.read(read(path, String), Dict{String,Any})
end

function dig(obj, keys...; default = nothing)
    cur = obj
    for key in keys
        if cur isa AbstractDict && haskey(cur, key)
            cur = cur[key]
        else
            return default
        end
    end
    return cur
end

function as_float(x)
    x isa Real && return Float64(x)
    return nothing
end

function fmt(x)
    x === nothing && return "—"
    x isa Real && return isfinite(Float64(x)) ? @sprintf("%.6g", Float64(x)) : string(Float64(x))
    s = string(x)
    isempty(s) && return "—"
    return replace(s, "|" => "\\|")
end

function fmt_int(x)
    x === nothing && return "—"
    x isa Real && return string(Int(x))
    return fmt(x)
end

function run_iters(run)
    run === nothing && return nothing
    iters = dig(run, "bpsdp_state", "outer_iterations")
    iters !== nothing && return iters
    return length(get(run, "cg_iterations_per_outer", Any[]))
end

function run_wall(run)
    run === nothing && return nothing
    return get(run, "solve_wall_seconds_measured", nothing)
end

function run_status(run)
    run === nothing && return nothing
    return get(run, "termination_status", nothing)
end

function run_objective(run)
    run === nothing && return nothing
    return as_float(get(run, "objective_value", nothing))
end

function discover_run_dirs(input_dir::AbstractString)
    isdir(input_dir) || return String[]
    dirs = String[]
    for name in sort(readdir(input_dir))
        path = joinpath(input_dir, name)
        isdir(path) || continue
        isfile(joinpath(path, "report.json")) || continue
        isfile(joinpath(path, "kernel.json")) || continue
        push!(dirs, path)
    end
    return dirs
end

function row_from_dir(dir::AbstractString)
    report = read_json(joinpath(dir, "report.json"))
    kernel = read_json(joinpath(dir, "kernel.json"))
    runs = get(report, "runs", Dict{String,Any}())
    unpruned = get(runs, "unpruned", nothing)
    pruned = get(runs, "pruned", nothing)
    obj_un = run_objective(unpruned)
    obj_pr = run_objective(pruned)
    obj_diff = (obj_un === nothing || obj_pr === nothing) ? nothing : abs(obj_pr - obj_un)
    eq_diff = dig(report, "equivalence", "abs_diff")
    obj_diff = eq_diff === nothing ? obj_diff : obj_diff

    return Dict{String,Any}(
        "dir" => dir,
        "nk" => dig(report, "moment_problem", "nk"; default = replace(basename(dir), r"[^0-9]" => "")),
        "blocking" => dig(report, "moment_problem", "blocking"; default = "?"),
        "target_sector" => get(report, "target_sector", "?"),
        "n_moments" => get(kernel, "n_moments", dig(report, "moment_problem", "n_linear_moments")),
        "n_orphans" => get(kernel, "n_orphans", dig(report, "moment_problem", "n_free_keys")),
        "n_expected_zero" => get(kernel, "n_expected_zero", nothing),
        "n_unexpected" => get(kernel, "n_unexpected", get(kernel, "n_physical_unexpected", nothing)),
        "iters_unpruned" => run_iters(unpruned),
        "wall_unpruned" => run_wall(unpruned),
        "status_unpruned" => run_status(unpruned),
        "iters_pruned" => run_iters(pruned),
        "wall_pruned" => run_wall(pruned),
        "status_pruned" => run_status(pruned),
        "obj_diff" => obj_diff,
        "equivalence" => get(report, "equivalence", Dict{String,Any}()),
    )
end

function write_summary(path::AbstractString, rows::Vector{Dict{String,Any}}, input_dir::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# H4/Nk sector-pruning diagnostic")
        println(io)
        println(io, "- generated: ", now())
        println(io, "- input: `", input_dir, "`")
        println(io, "- red `n_unexpected` means target-sector orphan keys survived the symmetry classifier; that is the excess-kernel proxy.")
        println(io)

        if isempty(rows)
            println(io, "No run directories with `report.json` + `kernel.json` were found.")
            return
        end

        println(io, "| Nk | blocking | target sector | n_moments | n_orphans | n_expected_zero | n_unexpected | status_unpruned | bpsdp_iters_unpruned | walltime_unpruned | status_pruned | bpsdp_iters_pruned | walltime_pruned | obj_diff | artifact |")
        println(io, "|----|----------|---------------|-----------|-----------|-----------------|--------------|-----------------|----------------------|-------------------|---------------|--------------------|-----------------|----------|----------|")
        for row in rows
            n_unexpected = row["n_unexpected"]
            n_unexpected_cell = (n_unexpected isa Real && n_unexpected > 0) ? "🔴 **$(Int(n_unexpected))**" : fmt_int(n_unexpected)
            artifact = replace(row["dir"], input_dir => "")
            isempty(artifact) && (artifact = basename(row["dir"]))
            artifact = strip(artifact, ['/', '\\'])
            println(io,
                "| ", fmt_int(row["nk"]),
                " | ", fmt(row["blocking"]),
                " | `", fmt(row["target_sector"]), "`",
                " | ", fmt_int(row["n_moments"]),
                " | ", fmt_int(row["n_orphans"]),
                " | ", fmt_int(row["n_expected_zero"]),
                " | ", n_unexpected_cell,
                " | ", fmt(row["status_unpruned"]),
                " | ", fmt_int(row["iters_unpruned"]),
                " | ", fmt(row["wall_unpruned"]),
                " | ", fmt(row["status_pruned"]),
                " | ", fmt_int(row["iters_pruned"]),
                " | ", fmt(row["wall_pruned"]),
                " | ", fmt(row["obj_diff"]),
                " | `", artifact, "` |")
        end
        println(io)

        println(io, "## Notes")
        println(io)
        for row in rows
            eq = row["equivalence"]
            checked = dig(eq, "checked"; default = false)
            if checked == true
                verdict = dig(eq, "passed"; default = false) == true ? "PASS" : "FAIL"
                println(io, "- Nk=", row["nk"], ": objective equivalence ", verdict,
                    " (`diff=", fmt(dig(eq, "abs_diff")), "`, `tol=", fmt(dig(eq, "tolerance")), "`).")
            else
                println(io, "- Nk=", row["nk"], ": objective equivalence not checked (",
                    fmt(dig(eq, "reason"; default = "missing both optimal solves")), ").")
            end
        end
    end
    return path
end

function main(argv = ARGS)
    options = parse_options(argv)
    dirs = discover_run_dirs(options.input_dir)
    rows = [row_from_dir(dir) for dir in dirs]
    sort!(rows; by = row -> (try parse(Int, string(row["nk"])) catch; typemax(Int) end, string(row["blocking"]), string(row["target_sector"])))
    out = write_summary(options.output_path, rows, options.input_dir)
    println("Wrote ", out)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
