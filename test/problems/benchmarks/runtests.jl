# =============================================================================
# test/problems/benchmarks/runtests.jl
# =============================================================================
# Tests: Classical polynomial optimization benchmarks
# Dependencies: SOLVER
# Requires --local: no
# =============================================================================

using Test

@testset "Benchmarks" begin
    include("ncpop_benchmarks.jl")
end

