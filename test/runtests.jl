using NCTSSoS, Test
import NCTSSoS: degree

include("Expectations.jl")
using .TestExpectations: expectations_oracle

@testset "NCTSSoS.jl" begin

    include("expectations_loader.jl")

    # 1. Polynomial algebra (no solver)
    include("polynomials/runtests.jl")

    # 2. Code quality (Aqua, ExplicitImports, doctests)
    include("quality/runtests.jl")

    # 3. Solver-dependent tests
    include("TestUtils.jl")

    # 4. Relaxation components
    include("relaxations/runtests.jl")

    # 5. State polynomials suite
    include("state_poly/runtests.jl")

    # 6. Correlated sparsity suite
    include("correlated_sparsity/runtests.jl")

    # 7. Trace polynomial suite
    include("trace_poly/runtests.jl")

    # 8. Curated problems
    @testset "Problems" begin
        include("problems/bell_inequalities/chsh_simple.jl")
        include("problems/benchmarks/e1_broyden_banded.jl")
        include("problems/benchmarks/e2_chained_singular.jl")
        include("problems/benchmarks/e3_generalized_rosenbrock.jl")
        include("problems/benchmarks/e4_chained_wood.jl")
        include("problems/benchmarks/e5_broyden_tridiagonal.jl")
        include("problems/nc_polynomial/nc_example1.jl")
        include("problems/nc_polynomial/nc_example2.jl")
        include("problems/trace_polynomial/wang_magron_example_5_3.jl")
    end
end
