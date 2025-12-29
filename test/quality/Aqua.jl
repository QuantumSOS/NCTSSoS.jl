# Aqua.jl code quality checks
# Aqua is a test dependency - only loaded through Pkg.test()
using Test

@testset "Aqua" begin
    using Aqua, NCTSSoS
    Aqua.test_all(NCTSSoS)
end
