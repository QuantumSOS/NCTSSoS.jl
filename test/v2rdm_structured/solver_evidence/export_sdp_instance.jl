#!/usr/bin/env julia
# Export a lowered periodic PQG-V2RDM SDP as a BPSDP-native sparse payload.
#
# This script is a manual HAI artifact, not part of the test suite. It builds the
# NCTSSoS moment data once, lowers to the complex PSD-block JuMP model, copies
# that model through BPSDP's MOI bridge, and serializes the resulting
# BPSDP.Problem data (A, b, c, block layout) for fast standalone solver work.
include(joinpath(@__DIR__, "driver.jl"))
include(joinpath(@__DIR__, "bpsdp_native_io.jl"))

using .BPSDPNativeIO
using Serialization
using SHA
using SparseArrays

Base.@kwdef mutable struct ExportOptions
    system::String = ""
    nk::Int = 0
    norb::Int = 0
    nelec_per_cell::Int = 0
    integrals::String = ""
    meta::String = ""
    output_dir::String = ""
    dependent_rows::Vector{Symbol} = [:keep]
    dependent_row_tol::Float64 = 1e-12
    solver_threads::Int = parse(Int, get(ENV, "SOLVER_THREADS", "8"))
    force::Bool = false
    smoke_load::Bool = true
end

function parse_dependent_rows(s::AbstractString)
    t = lowercase(strip(s))
    t in ("both", "all") && return [:keep, :drop]
    rows = Symbol[]
    for raw in split(t, ',')
        v = strip(raw)
        isempty(v) && continue
        v in ("keep", "drop") || throw(ArgumentError("unknown dependent row mode $(repr(v)); expected keep, drop, or both"))
        push!(rows, Symbol(v))
    end
    isempty(rows) && throw(ArgumentError("--dependent-rows names no modes"))
    return unique(rows)
end

function parse_export_options(argv)
    o = ExportOptions()
    for arg in argv
        val() = split(arg, "=", limit = 2)[2]
        if startswith(arg, "--system="); o.system = lowercase(val())
        elseif startswith(arg, "--nk="); o.nk = parse(Int, val())
        elseif startswith(arg, "--norb="); o.norb = parse(Int, val())
        elseif startswith(arg, "--nelec-per-cell="); o.nelec_per_cell = parse(Int, val())
        elseif startswith(arg, "--integrals="); o.integrals = val()
        elseif startswith(arg, "--meta="); o.meta = val()
        elseif startswith(arg, "--output-dir="); o.output_dir = val()
        elseif startswith(arg, "--dependent-rows="); o.dependent_rows = parse_dependent_rows(val())
        elseif startswith(arg, "--dependent-row-tol="); o.dependent_row_tol = parse(Float64, val())
        elseif startswith(arg, "--solver-threads="); o.solver_threads = parse(Int, val())
        elseif arg == "--force"; o.force = true
        elseif arg == "--no-smoke-load"; o.smoke_load = false
        elseif arg in ("-h", "--help"); print_export_help(); exit(0)
        else; throw(ArgumentError("unknown argument: $arg"))
        end
    end
    o.system in ("h2", "h4") || throw(ArgumentError("--system must be h2 or h4"))
    o.nk > 0 && o.norb > 0 && o.nelec_per_cell > 0 || throw(ArgumentError("--nk, --norb, --nelec-per-cell must be positive"))
    for (name, path) in (("--integrals", o.integrals), ("--meta", o.meta), ("--output-dir", o.output_dir))
        isempty(path) && throw(ArgumentError("$name is required"))
    end
    o.dependent_row_tol >= 0 || throw(ArgumentError("--dependent-row-tol must be nonnegative"))
    return o
end

function print_export_help()
    println("Usage: julia --project=ENV export_sdp_instance.jl --system=h2|h4 --nk=N --norb=N --nelec-per-cell=N --integrals=PATH --meta=PATH --output-dir=DIR [--dependent-rows=keep|drop|both]")
    println()
    println("Writes BPSDP-native .jls payload(s): bpsdp_native_<dependent_rows>.jls")
end

function git_dirty()
    try
        return !isempty(readchomp(pipeline(`git status --porcelain`, stderr = devnull)))
    catch
        return nothing
    end
end

function module_git_rev(mod)
    try
        src = pathof(mod)
        src === nothing && return nothing
        repo = dirname(dirname(src))
        return readchomp(pipeline(`git -C $repo rev-parse --short HEAD`, stderr = devnull))
    catch
        return nothing
    end
end

function module_git_dirty(mod)
    try
        src = pathof(mod)
        src === nothing && return nothing
        repo = dirname(dirname(src))
        return !isempty(readchomp(pipeline(`git -C $repo status --porcelain`, stderr = devnull)))
    catch
        return nothing
    end
end

env_or(name::AbstractString, fallback) = haskey(ENV, name) && !isempty(ENV[name]) ? ENV[name] : fallback

function sha256_file(path::AbstractString)
    return open(path, "r") do io
        bytes2hex(SHA.sha256(io))
    end
end

function native_filename(mode::Symbol, tol::Real)
    if mode == :keep
        return "bpsdp_native_keep.jls"
    elseif mode == :drop
        tol_s = replace(string(tol), "+" => "", "." => "p", "-" => "m")
        return "bpsdp_native_drop_tol$(tol_s).jls"
    else
        throw(ArgumentError("unknown dependent row mode $mode"))
    end
end

function copy_to_bpsdp_problem(model; dependent_rows::Symbol, dependent_row_tol::Real)
    # Do not call MOI.copy_to(BPSDP.Optimizer, JuMP.backend(model)) directly here.
    # The JuMP backend still contains complex affine equalities; JuMP's bridge
    # layer is what splits them into real rows before BPSDP's native MOI copy path.
    JuMP.set_optimizer(model, () -> BPSDP.Optimizer(
        max_iter = 0,
        print_level = 0,
        dependent_rows = dependent_rows,
        dependent_row_tol = dependent_row_tol,
    ))
    JuMP.set_silent(model)
    copy_seconds = @elapsed JuMP.optimize!(model)
    raw = raw_optimizer(model)
    problem = field_or_nothing(raw, :problem)
    problem === nothing && error("BPSDP optimize!(max_iter=0) did not build a native Problem")
    return raw, problem, copy_seconds
end

function base_metadata(opts::ExportOptions, meta, problem_stats, build_seconds, lowering_seconds, model)
    sdp_size = jump_model_size(model)
    return Dict{String,Any}(
        "format" => BPSDPNativeIO.FORMAT_NAME,
        "format_version" => BPSDPNativeIO.FORMAT_VERSION,
        "system" => opts.system,
        "nk" => opts.nk,
        "norb" => opts.norb,
        "nelec_per_cell" => opts.nelec_per_cell,
        "integrals" => abspath(opts.integrals),
        "meta" => abspath(opts.meta),
        "pyscf" => Dict(
            "ehf_Ha" => get(meta, "ehf", nothing),
            "active_hf_Ha" => get(meta, "active_hf", nothing),
            "energy_shift_Ha" => get(meta, "energy_shift", nothing),
            "basis" => get(meta, "basis", nothing),
            "density_fit" => get(meta, "density_fit", nothing),
            "text_h1e_entries" => get(meta, "text_h1e_entries", nothing),
            "text_eri_entries" => get(meta, "text_eri_entries", nothing),
        ),
        "lowering" => Dict(
            "formulation" => "psd_blocks",
            "representation" => "complex",
            "orphan_policy" => "error",
            "wall_seconds" => lowering_seconds,
        ),
        "moment_problem" => problem_stats,
        "build_wall_seconds" => build_seconds,
        "sdp_size" => sdp_size,
        "provenance" => Dict(
            "built_at" => string(Dates.now()),
            "hostname" => gethostname(),
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "git_rev" => env_or("NCTSSOS_GIT_REV", git_rev()),
            "git_dirty" => env_or("NCTSSOS_GIT_DIRTY", git_dirty()),
            "bpsdp_path" => String(pathof(BPSDP)),
            "bpsdp_git_rev" => module_git_rev(BPSDP),
            "bpsdp_git_dirty" => module_git_dirty(BPSDP),
            "runtime" => common_runtime_json(),
        ),
    )
end

function validate_native_payload!(problem, data)
    A = data[:A]
    size(A, 1) == length(data[:b]) || throw(DimensionMismatch("size(A,1) != length(b)"))
    size(A, 2) == length(data[:c]) || throw(DimensionMismatch("size(A,2) != length(c)"))
    expected_primal = sum(Int(n)^2 for n in data[:block_dims])
    expected_primal == length(data[:c]) || throw(DimensionMismatch("sum(block_dims.^2) != length(c)"))
    length(problem.c) == length(data[:c]) || throw(DimensionMismatch("reconstructed problem c length mismatch"))
    return nothing
end

function export_variant!(model, opts::ExportOptions, metadata0, mode::Symbol)
    outpath = joinpath(opts.output_dir, native_filename(mode, opts.dependent_row_tol))
    if isfile(outpath) && !opts.force
        @info "skip existing native export" outpath
        loaded_problem, loaded_data = BPSDPNativeIO.read_native(outpath)
        validate_native_payload!(loaded_problem, loaded_data)
        return Dict{String,Any}(
            "path" => outpath,
            "skipped" => true,
            "bytes" => filesize(outpath),
            "sha256" => sha256_file(outpath),
            "dependent_rows" => string(mode),
            "dependent_row_tol" => opts.dependent_row_tol,
            "moi_canonicalize_seconds" => 0.0,
            "native_write_seconds" => 0.0,
            "native_dimensions" => BPSDPNativeIO.dimensions(loaded_data),
        )
    end

    raw = nothing
    problem = nothing
    copy_seconds = 0.0
    raw, problem, copy_seconds = copy_to_bpsdp_problem(
        model;
        dependent_rows = mode,
        dependent_row_tol = opts.dependent_row_tol,
    )

    native_dims = BPSDPNativeIO.dimensions(problem)
    variant_meta = copy(metadata0)
    variant_meta["native"] = merge(native_dims, Dict{String,Any}(
        "dependent_rows" => string(mode),
        "dependent_row_tol" => opts.dependent_row_tol,
        "moi_canonicalize_seconds" => copy_seconds,
        "objective_sign" => field_or_nothing(raw, :objective_sign),
        "objective_constant" => field_or_nothing(raw, :objective_constant),
    ))

    write_seconds = @elapsed BPSDPNativeIO.write_native(outpath, problem; metadata = json_ready(variant_meta))
    checksum = sha256_file(outpath)

    smoke = nothing
    if opts.smoke_load
        loaded_problem = nothing
        loaded_data = nothing
        load_seconds = @elapsed loaded_problem, loaded_data = BPSDPNativeIO.read_native(outpath)
        validate_native_payload!(loaded_problem, loaded_data)
        smoke = Dict{String,Any}(
            "load_seconds" => load_seconds,
            "validated" => true,
            "dimensions" => BPSDPNativeIO.dimensions(loaded_data),
        )
    end

    return Dict{String,Any}(
        "path" => outpath,
        "bytes" => filesize(outpath),
        "sha256" => checksum,
        "dependent_rows" => string(mode),
        "dependent_row_tol" => opts.dependent_row_tol,
        "moi_canonicalize_seconds" => copy_seconds,
        "native_write_seconds" => write_seconds,
        "smoke_load" => smoke,
        "native_dimensions" => native_dims,
        "objective_sign" => field_or_nothing(raw, :objective_sign),
        "objective_constant" => field_or_nothing(raw, :objective_constant),
    )
end

function write_checksums(path::AbstractString, variants)
    open(path, "w") do io
        for v in variants
            get(v, "skipped", false) && continue
            println(io, v["sha256"], "  ", basename(v["path"]))
        end
    end
end

function main(argv = ARGS)
    opts = parse_export_options(argv)
    validate_output_dir!(opts.output_dir)
    pin_threads!(Options(solver_threads = opts.solver_threads, system = opts.system, nk = opts.nk,
        norb = opts.norb, nelec_per_cell = opts.nelec_per_cell, integrals = opts.integrals,
        meta = opts.meta, output_dir = opts.output_dir))

    meta = read_meta(opts.meta)
    driver_opts = Options(
        system = opts.system,
        nk = opts.nk,
        norb = opts.norb,
        nelec_per_cell = opts.nelec_per_cell,
        integrals = opts.integrals,
        meta = opts.meta,
        output_dir = opts.output_dir,
        solver_threads = opts.solver_threads,
    )

    mp, problem_stats, build_seconds = build_problem(driver_opts)
    n_orphans = get(problem_stats, "n_orphans", missing)
    n_orphans == 0 || error("refusing to export SDP with n_orphans=$n_orphans")

    model = nothing
    lowering_seconds = @elapsed model, _ = build_jump_model(mp; formulation = :psd_blocks, representation = :complex)
    metadata0 = base_metadata(opts, meta, problem_stats, build_seconds, lowering_seconds, model)

    variants = Any[]
    total_export_seconds = @elapsed begin
        for mode in opts.dependent_rows
            @info "export BPSDP-native variant" system=opts.system nk=opts.nk dependent_rows=mode output_dir=opts.output_dir
            push!(variants, export_variant!(model, opts, metadata0, mode))
        end
    end

    summary = copy(metadata0)
    summary["variants"] = variants
    summary["export_wall_seconds"] = total_export_seconds
    summary_path = joinpath(opts.output_dir, "metadata.json")
    write_json(summary_path, summary)
    write_checksums(joinpath(opts.output_dir, "checksums.sha256"), variants)

    @info "export complete" output_dir=opts.output_dir metadata=summary_path variants=length(variants)
    for v in variants
        @printf("%s rows=%s A=%dx%d nnz=%d bytes=%s copy=%.3fs write=%.3fs\n",
            basename(v["path"]), v["dependent_rows"],
            v["native_dimensions"]["dual_dimension"], v["native_dimensions"]["primal_dimension"],
            v["native_dimensions"]["A_nnz"], fmt(v["bytes"]),
            v["moi_canonicalize_seconds"], v["native_write_seconds"])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
