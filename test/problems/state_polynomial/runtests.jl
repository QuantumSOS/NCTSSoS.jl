# =============================================================================
# test/problems/state_polynomial/runtests.jl
# =============================================================================
# Tests: State polynomial optimization examples
# Dependencies: SOLVER
# Requires --local: no
# =============================================================================

using Test

@testset "State Polynomial Examples" begin
    include("state_polynomial.jl")
end

