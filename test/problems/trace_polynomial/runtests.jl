# =============================================================================
# Trace Polynomial Tests Runner
# =============================================================================
# Tests for trace polynomial optimization (Section 6.x examples)
# =============================================================================

using Test

@testset "Trace Polynomial" begin
    include("trace_polynomial.jl")
    include("sparsity_variants.jl")
end
