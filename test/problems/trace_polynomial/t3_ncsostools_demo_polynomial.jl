# test/problems/trace_polynomial/t3_ncsostools_demo_polynomial.jl
# T3: NCSOStools Demo Polynomial
#
# Source: Bhardwaj, Klep & Magron, "Noncommutative Polynomial
# Optimization" (Example 4.4); benchmark metadata mirrored in
# `test-examples.typ` as T3.
#
# The example polynomial
#
#     t(x, y) = 1 + 2x + x^2 + x y^2 + 2y^2 + y^2 x + y x^2 y + y^4
#
# is a sum of Hermitian squares, so the exact tracial optimum is 0. The
# literature reports tr_Θ^(2)(t) ≈ 6.97e-10, so this regression keeps only the
# dense order-2 tracial relaxation that is actually literature-backed.
#
# Sparse and basis-reduced variants are useful engineering checks, but they are
# not part of the reviewed expectation fixture for this benchmark.

using Test, NCTSSoS

# Expectations in test/data/expectations/t3_ncsostools_demo_polynomial.toml

const T3_OBJECTIVE_ATOL = 1e-6

function t3_ncsostools_demo_problem()
    reg, (vars,) = create_noncommutative_variables([("x", 1:2)])
    x, y = vars

    t = 1.0 + 2.0 * x + x^2 + x * y^2 + 2.0 * y^2 + y^2 * x + y * x^2 * y + y^4
    return polyopt(tr(t) * one(x), reg)
end

@testset "T3 NCSOStools Demo Polynomial" begin
    pop = t3_ncsostools_demo_problem()

    dense_oracle = expectations_oracle("expectations/t3_ncsostools_demo_polynomial.toml", "dense_order2")
    dense_config = SolverConfig(
        optimizer=SOLVER,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )
    dense_result = cs_nctssos(pop, dense_config; dualize=true)

    @test dense_result.objective ≈ dense_oracle.opt atol = T3_OBJECTIVE_ATOL
    @test flatten_sizes(dense_result.moment_matrix_sizes) == dense_oracle.sides
    @test dense_result.n_unique_moment_matrix_elements == dense_oracle.nuniq
end
