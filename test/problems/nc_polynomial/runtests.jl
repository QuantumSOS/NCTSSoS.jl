# test/problems/nc_polynomial/runtests.jl
# Tests: Noncommutative polynomial optimization examples
#
# Coverage:
# - nc_example1.jl: Example 1 (unconstrained, 3 vars)
# - nc_example2.jl: Example 2 (constrained, 2 vars)
# - nc_correlative.jl: Correlative sparsity with constraints
# - nc_large_scale.jl: CS+TS n=10

using Test

@testset "NC Polynomial Examples" begin
    include("nc_example1.jl")
    include("nc_example2.jl")
    include("nc_correlative.jl")
    include("nc_large_scale.jl")
end
