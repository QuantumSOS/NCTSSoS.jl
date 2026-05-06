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
        include("problems/bell_inequalities/pironio_toy.jl")
        include("problems/benchmarks/e1_broyden_banded.jl")
        include("problems/benchmarks/e2_chained_singular.jl")
        include("problems/benchmarks/e3_generalized_rosenbrock.jl")
        include("problems/benchmarks/e4_chained_wood.jl")
        include("problems/benchmarks/e5_broyden_tridiagonal.jl")
        include("problems/condensed_matter/ising.jl")
        include("problems/condensed_matter/hubbard.jl")
        include("problems/condensed_matter/bose_hubbard.jl")
        include("problems/fermionic/fermionic_chain.jl")
        include("problems/fermionic/xy_model.jl")
        include("problems/fermionic/free_fermion.jl")
        include("problems/condensed_matter/harmonic_oscillator.jl")
        include("problems/nc_polynomial/nc_example1.jl")
        include("problems/nc_polynomial/nc_example2.jl")
        include("problems/state_polynomial/survey_comparative.jl")
        include("problems/trace_polynomial/t2_broyden_tridiagonal_trace.jl")
        include("problems/trace_polynomial/t3_ncsostools_demo_polynomial.jl")
        include("problems/trace_polynomial/t5_constrained_trace_semialgebraic_set.jl")
        include("problems/trace_polynomial/t6_projection_optimization.jl")
    end
end
