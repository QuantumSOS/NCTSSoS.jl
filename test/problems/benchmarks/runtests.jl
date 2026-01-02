# =============================================================================
# Benchmark Tests Runner
# =============================================================================
# Tests for classical optimization benchmarks:
#   - NCPOP benchmarks (Rosenbrock, etc.)
# =============================================================================

using Test

@testset "Benchmarks" begin
    include("ncpop_benchmarks.jl")
end
