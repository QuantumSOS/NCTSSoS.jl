# test/quality/runtests.jl - Code Quality Checks

using Test

@testset "Quality" begin
    include("Aqua.jl")
    include("ExplicitImports.jl")
    include("Doctest.jl")
end
