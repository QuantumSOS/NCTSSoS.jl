# test/problems/bell_inequalities/pironio_toy.jl
# Tests: E11 Pironio toy problem.
#
# Coverage: dense moment + dense SOS at order 2.
# Expected optimal value: -3/4.
#
# Note on E12:
#   Pironio's "toy with extra constraints" adds generalized state constraints
#       r(X)|φ⟩ = 0
#   and
#       ⟨φ|s(X)|φ⟩ ≥ 0.
#   Those are state-polynomial features, not a clean Bell/eigenvalue benchmark in
#   this file. A surrogate reformulation can be written with the current API, but
#   that would test a state-polynomial encoding choice rather than the original
#   Bell-problem formulation. So E12 is intentionally not included here; it should
#   be covered in the state-polynomial test suite instead.

using Test, NCTSSoS, JuMP

function create_pironio_toy_problem()
    reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
    X1, X2 = vars

    objective = X1 * X2 + X2 * X1
    projector_constraint = X1^2 - X1
    localizing_constraint = -X2^2 + X2 + 0.5
    pop = polyopt(
        objective,
        reg;
        eq_constraints=[projector_constraint],
        ineq_constraints=[localizing_constraint],
    )

    return pop, reg
end

@testset "Pironio Toy Problem (E11)" begin
    oracle = expectations_oracle("expectations/pironio_toy.toml", "Dense_d2")

    @testset "Dense (Moment)" begin
        pop, _ = create_pironio_toy_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        result = cs_nctssos(pop, config; dualize=false)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Dense (SOS)" begin
        pop, _ = create_pironio_toy_problem()
        config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        result = cs_nctssos(pop, config; dualize=true)

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
