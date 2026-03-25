# test/problems/trace_polynomial/t5_constrained_trace_semialgebraic_set.jl
# T5: Constrained Trace on Semialgebraic Set
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Example 5.3); benchmark metadata mirrored in
# `test-examples.typ` as T5.
#
# This is the trace-formulation counterpart of the eigenvalue benchmark E7.
# The constrained problem
#
#     f = 2 - X^2 + X Y^2 X - Y^2
#     S = {4 - X^2 - Y^2, X Y + Y X - 2}
#
# has exact order-2 optimum -1, and Wang--Magron report that the first
# term-sparse step already matches the dense relaxation:
#
#     μ₂(f, S) = μ^(ts)_(2,1)(f, S) = -1.
#
# The regression here is small and literature-backed, so it belongs in CI.

using Test, NCTSSoS

using NCTSSoS:
    MaxEntangled,
    NCStatePolynomial,
    NCStateWord,
    NonCommutativeAlgebra,
    StateWord

const T5_OBJECTIVE_ATOL = 1e-6
const T5_EXPECTATIONS_PATH = "expectations/t5_constrained_trace_semialgebraic_set.toml"

function t5_trace_constraint(poly, state_word_id)
    coeffs = [coef for (coef, _) in poly.terms]
    words = [NCStateWord(state_word_id, mono) for (_, mono) in poly.terms]
    return NCStatePolynomial(coeffs, words)
end

function t5_constrained_trace_problem()
    reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
    X, Y = vars

    f = 2.0 - X^2 + X * Y^2 * X - Y^2
    g_ball = 4.0 - X^2 - Y^2
    g_link = X * Y + Y * X - 2.0

    state_word_id = one(StateWord{MaxEntangled, NonCommutativeAlgebra, eltype(X.word)})
    return polyopt(
        tr(f) * one(typeof(X)),
        reg;
        ineq_constraints=[t5_trace_constraint(g_ball, state_word_id)],
        eq_constraints=[t5_trace_constraint(g_link, state_word_id)],
    )
end

@testset "T5 Constrained Trace on Semialgebraic Set" begin
    pop = t5_constrained_trace_problem()
    dense_oracle = expectations_oracle(T5_EXPECTATIONS_PATH, "dense_order2")
    sparse_oracle = expectations_oracle(T5_EXPECTATIONS_PATH, "ts_mmd_order2")

    dense_config = SolverConfig(
        optimizer=SOLVER,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
    )
    dense_result = cs_nctssos(pop, dense_config; dualize=true)
    @test dense_result.objective ≈ dense_oracle.opt atol = T5_OBJECTIVE_ATOL

    sparse_config = SolverConfig(
        optimizer=SOLVER,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=MMD(),
    )
    sparse_result = cs_nctssos(pop, sparse_config; dualize=true)
    @test sparse_result.objective ≈ sparse_oracle.opt atol = T5_OBJECTIVE_ATOL
end
