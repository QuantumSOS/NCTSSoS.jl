# Aqua.jl code quality checks
using Test

@testset "Aqua" begin
    using Aqua, NCTSSoS
    Aqua.test_all(NCTSSoS)
end
