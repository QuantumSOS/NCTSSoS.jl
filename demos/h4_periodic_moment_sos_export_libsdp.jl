#!/usr/bin/env julia
# Export demos/h4_periodic_moment_sos.jl as a libsdp-readable complex
# SDPA-sparse problem.
#
# This builds the symbolic H4/Nk=2 periodic 1D+PQG MomentProblem once, then
# writes:
#
#   <outdir>/<basename>.dat-c     — libsdp complex SDPA-sparse text
#   <outdir>/<basename>_meta.txt  — block/pivot/provenance summary
#
# Usage from the repository root:
#   julia --project=. demos/h4_periodic_moment_sos_export_libsdp.jl
#   julia --project=. demos/h4_periodic_moment_sos_export_libsdp.jl \
#       --integrals=output/h4_chain_nk2_figure1_integrals_drop1e-12.txt \
#       --outdir=output/h4_libsdp --basename=h4_periodic_moment_sos
#
# Formulation flags are forwarded to h4_periodic_moment_sos.jl:
#   --integrals=PATH
#   --blocking=momentum|spin|none
#
# Export flags:
#   --outdir=PATH                Default: output/h4_periodic_moment_sos_libsdp
#   --basename=STEM              Default: h4_periodic_moment_sos
#   --sense=min|max              Default: min
#   --orphans-per-block=N        Default: 32 (set 1 for true 2x2 aux blocks)

using Printf

include(joinpath(@__DIR__, "h4_periodic_moment_sos.jl"))
include(joinpath(@__DIR__, "SDPALibsdpExport.jl"))
using .SDPALibsdpExport: export_libsdp

function _print_export_help()
    println("Usage: julia --project=. demos/h4_periodic_moment_sos_export_libsdp.jl [export options] [formulation options]\n")
    println("Export options:")
    println("  --outdir=PATH                output directory; default output/h4_periodic_moment_sos_libsdp")
    println("  --basename=STEM              output stem; default h4_periodic_moment_sos")
    println("  --sense=min|max              objective sense recorded for libsdp; default min")
    println("  --orphans-per-block=N        orphan moment packing; default 32")
    println()
    println("Formulation options forwarded to h4_periodic_moment_sos.jl:")
    println("  --integrals=PATH")
    println("  --blocking=momentum|spin|none")
    println()
    println("After export, parse/solve with:")
    println("  python demos/solve_libsdp_dat_c.py --problem <outdir>/<basename>.dat-c --no-solve")
    return nothing
end

function _split_export_args(argv)
    formulation_args = String[]
    outdir = normpath(joinpath(@__DIR__, "..", "output", "h4_periodic_moment_sos_libsdp"))
    basename = "h4_periodic_moment_sos"
    sense = :min
    orphans_per_block = 32

    for arg in argv
        if startswith(arg, "--outdir=")
            outdir = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--basename=")
            basename = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--sense=")
            sense = Symbol(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--orphans-per-block=")
            orphans_per_block = parse(Int, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            _print_export_help()
            exit(0)
        else
            push!(formulation_args, arg)
        end
    end

    sense in (:min, :max) || throw(ArgumentError("--sense must be min or max"))
    orphans_per_block >= 1 || throw(ArgumentError("--orphans-per-block must be >= 1"))

    return formulation_args, outdir, basename, sense, orphans_per_block
end

function main(argv = ARGS)
    formulation_args, outdir, basename, sense, orphans_per_block =
        _split_export_args(argv)

    options = parse_options(formulation_args)

    @printf("%-48s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h4_pqg_moment_problem(options)
    print_summary(data, options, build_seconds)
    println()

    @printf("%-48s %s\n", "writing libsdp export to", outdir)
    @printf("%-48s %s\n", "basename", basename)
    @printf("%-48s %s\n", "sense", string(sense))
    @printf("%-48s %d\n", "orphans per aux block", orphans_per_block)

    export_seconds = @elapsed summary = export_libsdp(data.moment_problem;
        outdir,
        basename,
        sense,
        orphans_per_block,
    )

    dat_path = joinpath(outdir, basename * ".dat-c")
    meta_path = joinpath(outdir, basename * "_meta.txt")

    println()
    println("== libsdp export summary ==")
    @printf("%-48s %s\n", "SDPA file", dat_path)
    @printf("%-48s %s\n", "metadata", meta_path)
    @printf("%-48s %d (HPSD=%d, aux=%d)\n",
            "PSD blocks (complex Hermitian)",
            summary.n_blocks, summary.n_hpsd_blocks, summary.n_aux_blocks)
    @printf("%-48s %d / %d\n",
            "min/max block dim",
            minimum(summary.block_dims),
            maximum(summary.block_dims))
    @printf("%-48s %d\n", "unique canonical moments (= pivots)", summary.n_unique_moments)
    @printf("%-48s %d\n", "  HPSD-anchored pivots", summary.n_hpsd_pivots)
    @printf("%-48s %d\n", "  aux-block pivots (orphans rehomed)", summary.n_aux_pivots)
    @printf("%-48s %d\n", "equality constraints (libsdp rows)", summary.n_constraints)
    @printf("%-48s %d\n", "objective sparse entries", summary.n_objective_entries)
    @printf("%-48s %d\n", "constraint sparse entries", summary.n_constraint_entries)
    @printf("%-48s %d\n", "Zero polynomial-matrix components", summary.n_zero_constraint_components)
    @printf("%-48s %d\n", "dropped orphan term occurrences", summary.n_dropped_orphan_terms)
    @printf("%-48s %.3f s\n", "export walltime", export_seconds)
    @printf("%-48s %s\n", "max RSS after export", rss_string())
    println()
    println("Parse-only sanity check:")
    println("    python demos/solve_libsdp_dat_c.py --problem ", dat_path, " --no-solve")
    println("Solve with libsdp complex BPSDP:")
    println("    python demos/solve_libsdp_dat_c.py --problem ", dat_path)

    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
