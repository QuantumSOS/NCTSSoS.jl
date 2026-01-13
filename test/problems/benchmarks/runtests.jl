# test/problems/benchmarks/runtests.jl
# Tests: Classical polynomial optimization benchmarks

using Test

@testset "Benchmarks" begin
    include("ncpop_benchmarks.jl")
end
