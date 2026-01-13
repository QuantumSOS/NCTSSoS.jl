# test/problems/trace_polynomial/runtests.jl
# Tests: Trace polynomial optimization examples + sparsity variants

using Test

@testset "Trace Polynomial Examples" begin
    include("trace_polynomial.jl")
    include("sparsity_variants.jl")
end
