using NCTSSoS, Test

@testset "NCTSSoS.jl" begin

    # 1. Polynomial algebra (no solver)
    include("polynomials/runtests.jl")

    # 2. Code quality (Aqua, ExplicitImports, doctests)
    include("quality/runtests.jl")

    # 3. Solver-dependent tests
    include("TestUtils.jl")

    # 4. Relaxation components
    include("relaxations/runtests.jl")

    # 5. Curated problems
    @testset "Problems" begin
        include("problems/bell_inequalities/chsh_simple.jl")
        include("problems/nc_polynomial/nc_example1.jl")
        include("problems/nc_polynomial/nc_example2.jl")
        include("problems/bell_inequalities/chsh_trace.jl")
        include("problems/bell_inequalities/chsh_state.jl")
    end
end
