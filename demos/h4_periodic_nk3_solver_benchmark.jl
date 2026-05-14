#!/usr/bin/env julia
raw"""
Benchmark BPSDP.jl, SCS.jl indirect/CG, and COSMO.jl CGIndirectKKTSolver on
H4/Nk=3 periodic 1D+PQG V2RDM.

This is a thin Nk=3 launcher over `h4_periodic_nk2_solver_benchmark.jl`, whose
builder is parameterized by `--nk`. It keeps the same lowering and solver knobs
used for the matched-thread Nk=2 benchmark.

Usage from the repository root on HAI:

    /home/ubuntu/.juliaup/bin/julia --startup-file=no --project=/tmp/tmp.EOd2eUfCLO \
        demos/h4_periodic_nk3_solver_benchmark.jl \
        --include-1d --paper-spin --blocking=momentum \
        --solvers=cosmo-cg,scs-indirect,bpsdp \
        --repeats=1 \
        --output-dir=output/phase2/h4_figure1_nk3_1d_paperspin_threads8
"""

module H4NkBenchmark
include(joinpath(@__DIR__, "h4_periodic_nk2_solver_benchmark.jl"))
end

const DEFAULT_NK3_INTEGRALS = normpath(joinpath(
    @__DIR__, "..", "output", "h4_chain_nk3_figure1_integrals_drop1e-12.txt"))
const DEFAULT_NK3_OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "output", "phase2", "h4_figure1_nk3_1d_paperspin_solver_benchmark"))

function has_option(argv, name::AbstractString)
    return any(arg -> arg == name || startswith(arg, string(name, "=")), argv)
end

function main(argv = ARGS)
    args = String[]
    has_option(argv, "--nk") || push!(args, "--nk=3")
    has_option(argv, "--norb") || push!(args, "--norb=8")
    has_option(argv, "--nelec-per-cell") || push!(args, "--nelec-per-cell=4")
    has_option(argv, "--integrals") || push!(args, "--integrals=$DEFAULT_NK3_INTEGRALS")
    has_option(argv, "--output-dir") || push!(args, "--output-dir=$DEFAULT_NK3_OUTPUT_DIR")
    append!(args, argv)
    return H4NkBenchmark.main(args)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
