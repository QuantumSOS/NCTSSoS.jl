# =============================================================================
# NCTSSoS Test Suite - Shared Setup
# =============================================================================
# This file configures the SDP solver for all tests.
#
# Solver priority:
#   1. --local flag → Mosek (commercial, fast, high precision)
#   2. Default      → COSMO (open-source, reasonable performance)
#
# Usage:
#   Pkg.test("NCTSSoS"; test_args=["--local"])  # Use Mosek
#   Pkg.test("NCTSSoS")                          # Use COSMO
#
# In test files:
#   include("setup.jl")  # or rely on runtests.jl to include it
#   # SOLVER is now available
# =============================================================================

using JuMP

# Check if --local flag was passed (may already be defined by runtests.jl)
if !@isdefined(USE_LOCAL)
    const USE_LOCAL = "--local" in ARGS
end

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
