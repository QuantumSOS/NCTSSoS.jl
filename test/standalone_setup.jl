# =============================================================================
# NCTSSoS Standalone Test Setup - Always uses Mosek
# =============================================================================
# Use this for running individual test files directly.
# For suite testing, use runtests.jl which loads setup.jl.
#
# Usage:
#   make test-file FILE=test/relaxations/sparsity.jl
#   julia --project -e 'include("test/relaxations/sparsity.jl")'
#
# This file is included by test files when SOLVER is not already defined.
# =============================================================================

using NCTSSoS
using Test
using JuMP

# Always use Mosek for standalone testing (faster, more precise than COSMO)
const USE_LOCAL = true

using MosekTools
const SOLVER = optimizer_with_attributes(
    Mosek.Optimizer,
    "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
    "MSK_IPAR_LOG" => 0
)

@info "Standalone setup: Using Mosek"
