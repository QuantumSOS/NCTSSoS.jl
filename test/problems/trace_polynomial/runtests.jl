# =============================================================================
# test/problems/trace_polynomial/runtests.jl
# =============================================================================
# Tests: Trace polynomial optimization examples + sparsity variants
# Dependencies: SOLVER
# Requires --local: no
# =============================================================================

using Test

@testset "Trace Polynomial Examples" begin
    include("trace_polynomial.jl")
    include("sparsity_variants.jl")
end

