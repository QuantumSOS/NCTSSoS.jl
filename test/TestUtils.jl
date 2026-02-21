# test/TestUtils.jl - Shared Test Infrastructure
# Provides:
#   - SOLVER constant (COSMO default, Mosek with --local)
#   - USE_LOCAL flag
#   - flatten_sizes helper

using JuMP

# Parse --local flag (may already be defined by runtests.jl)
if !@isdefined(USE_LOCAL)
    const USE_LOCAL = "--local" in ARGS
end

# Configure solver
if !@isdefined(SOLVER)
    const SOLVER = if USE_LOCAL
        @info "Using Mosek (--local)"
        using MosekTools
        optimizer_with_attributes(
            Mosek.Optimizer,
            "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
            "MSK_IPAR_LOG" => 0
        )
    else
        @info "Using COSMO (default open-source solver)"
        using COSMO
        optimizer_with_attributes(
            COSMO.Optimizer,
            "verbose" => false,
            "eps_abs" => 1e-7,
            "eps_rel" => 1e-7
        )
    end
end

# Provide a stable solver name for oracle selection in test files.
if !@isdefined(SOLVER_NAME)
    const SOLVER_NAME = USE_LOCAL ? :mosek : :cosmo
end

# Shared Helpers

if !@isdefined(flatten_sizes)
    """
        flatten_sizes(sizes)

    Flatten nested moment matrix sizes for comparison with oracle values.

    # Example
    ```julia
    flatten_sizes([[3, 3], [2]]) == [3, 3, 2]
    ```
    """
    flatten_sizes(sizes) = reduce(vcat, sizes)
end
