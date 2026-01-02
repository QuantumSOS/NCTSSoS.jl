# =============================================================================
# NC Polynomial Examples Tests Runner
# =============================================================================
# Tests for noncommutative polynomial optimization examples:
#   - Example 1, 2 from NCTSSOS paper
#   - Correlative sparsity examples
# =============================================================================

using Test

@testset "NC Polynomial Examples" begin
    include("nc_examples.jl")
end
