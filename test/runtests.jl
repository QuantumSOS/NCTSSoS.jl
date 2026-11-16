using NCTSSoS, Test

include("Expectations.jl")
using .TestExpectations: expectations_oracle

@testset "NCTSSoS.jl" begin
    include("expectations_loader.jl")

    # Code quality (Aqua, ExplicitImports, doctests).
    include("quality/runtests.jl")

    # Shared solver and the complete ported suite.
    include("TestUtils.jl")
    include("next/runtests.jl")
end
