# test/problems/trace_polynomial/t6_projection_optimization.jl
# T6: Projection Optimization
#
# Source: Klep, Magron & Volčič, "Optimization Over Trace Polynomials"
# (Section 6.1); benchmark metadata mirrored in `test-examples.typ` as T6.
#
# The trace objective
#
#     tr(X₁ X₂ X₃) + tr(X₁ X₂) tr(X₃)
#
# over three projections has a strictly negative optimum, unlike the commutative
# case where every feasible value is nonnegative. The literature reports the
# order-2 dense relaxation value ≈ -0.0467 and the exact order-3 optimum -1/32.
#
# This example already appears in the public docs as the toy trace-polynomial
# walkthrough (`docs/src/examples/literate/trace_poly.jl`), so the backlog item
# belongs in CI as a curated regression rather than as new docs.

using Test, NCTSSoS

const T6_OBJECTIVE_ATOL = 1e-6
const T6_EXPECTATIONS_PATH = "expectations/t6_projection_optimization.toml"

function t6_projection_optimization_problem()
    reg, (P,) = create_projector_variables([("X", 1:3)])
    objective = (tr(P[1] * P[2] * P[3]) + tr(P[1] * P[2]) * tr(P[3])) * one(typeof(P[1]))
    return polyopt(objective, reg)
end

@testset "T6 Projection Optimization" begin
    pop = t6_projection_optimization_problem()

    @testset "dense order $(order)" for (case_id, order) in [
        ("dense_order2", 2),
        ("dense_order3", 3),
    ]
        oracle = expectations_oracle(T6_EXPECTATIONS_PATH, case_id)
        config = SolverConfig(
            optimizer=SOLVER,
            order=order,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        result = cs_nctssos(pop, config; dualize=true)

        @test result.objective ≈ oracle.opt atol = T6_OBJECTIVE_ATOL
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
