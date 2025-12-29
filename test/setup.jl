# =============================================================================
# NCTSSoS Test Suite - Shared Setup
# =============================================================================
# This file configures the SDP solver for all tests.
#
# Solver priority:
#   1. LOCAL_TESTING=true → Mosek (commercial, fast, high precision)
#   2. Default → COSMO (open-source, reasonable performance)
#
# Usage in test files:
#   include("setup.jl")  # or rely on runtests.jl to include it
#   # SOLVER is now available
# =============================================================================

using JuMP

# Try to load Mosek if LOCAL_TESTING is set
const LOCAL_TESTING = haskey(ENV, "LOCAL_TESTING")

const SOLVER = if LOCAL_TESTING
    @info "LOCAL_TESTING enabled, using Mosek"
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

# Quick solver for less demanding tests (e.g., polynomial algebra tests)
const QUICK_SOLVER = SOLVER
