# Survey comparative example tests (S6, S7)
#
# Huber, Magron & Volčič (2024), Section 3.
# One projector polynomial separates the state, trace, and commutative moment
# regimes:
#   - state optimum: -1/16
#   - trace optimum: -1/27
#   - commutative moment regime: nonnegative
#
# We keep the active regression minimal:
#   1. verify explicit optimal witnesses directly;
#   2. solve the dense state and trace relaxations on the same one-site projector
#      model and compare against reviewed expectations.
#
# Important: this benchmark must use a *single* projector family. Splitting the
# two projections across multiple prefix groups would place them on different
# sites, and NCTSSoS would commute them during simplification. The survey example
# needs two genuinely noncommuting projections.

using Test, NCTSSoS
using LinearAlgebra: Diagonal, dot, norm, tr

const SURVEY_COMPARATIVE_EXPECTATIONS_PATH = "expectations/state_survey_comparative.toml"
const S6_STATE_EXACT = -1 / 16
const S7_TRACE_EXACT = -1 / 27

function survey_projector_variables()
    reg, (P,) = create_projector_variables([("P", 1:2)])
    return reg, P[1], P[2]
end

function survey_state_problem()
    reg, X, Y = survey_projector_variables()
    objective = (
        0.5 * ς(X * Y * X * Y) +
        0.5 * ς(Y * X * Y * X) -
        ς(X * Y * X) * ς(Y)
    ) * one(typeof(X))
    return polyopt(objective, reg)
end

function survey_trace_problem()
    reg, X, Y = survey_projector_variables()
    objective = (
        0.5 * NCTSSoS.tr(X * Y * X * Y) +
        0.5 * NCTSSoS.tr(Y * X * Y * X) -
        NCTSSoS.tr(X * Y * X) * NCTSSoS.tr(Y)
    ) * one(typeof(X))
    return polyopt(objective, reg)
end

function survey_state_witness_value()
    X = [1.0 0.0; 0.0 0.0]
    Y = 0.5 .* [1.0 1.0; 1.0 1.0]

    # The survey witness is unique only up to basis relabeling / phase choices.
    # This basis-convention-equivalent vector evaluates the same projector pair to
    # the reported optimum -1/16.
    v = [sqrt(2 + sqrt(2)) / 2, sqrt(2 - sqrt(2)) / 2]

    xyxy = real(dot(v, (X * Y * X * Y) * v))
    yxyx = real(dot(v, (Y * X * Y * X) * v))
    xyx = real(dot(v, (X * Y * X) * v))
    y = real(dot(v, Y * v))

    return X, Y, v, 0.5 * (xyxy + yxyx) - xyx * y
end

function survey_trace_witness_value()
    # Direct-sum witness from the survey proof: one scalar (0, 1) block together
    # with the irreducible 2×2 projector pair at t = 1/3.
    X = Matrix(Diagonal([0.0, 1.0, 0.0]))
    Y = [
        1.0 0.0 0.0
        0.0 (1.0 / 3.0) (sqrt(2.0) / 3.0)
        0.0 (sqrt(2.0) / 3.0) (2.0 / 3.0)
    ]
    τ(M) = real(tr(M)) / size(M, 1)

    value = 0.5 * (τ(X * Y * X * Y) + τ(Y * X * Y * X)) - τ(X * Y * X) * τ(Y)
    return X, Y, value
end

@testset "Survey Comparative Example (S6, S7)" begin
    @testset "Benchmark setup" begin
        _, X, Y = survey_projector_variables()
        @test X * Y != Y * X
    end

    @testset "Explicit witnesses" begin
        Xs, Ys, v, state_value = survey_state_witness_value()
        @test norm(Xs * Xs - Xs) ≤ 1e-12
        @test norm(Ys * Ys - Ys) ≤ 1e-12
        @test norm(v) ≈ 1.0 atol = 1e-12
        @test state_value ≈ S6_STATE_EXACT atol = 1e-12

        Xt, Yt, trace_value = survey_trace_witness_value()
        @test norm(Xt * Xt - Xt) ≤ 1e-12
        @test norm(Yt * Yt - Yt) ≤ 1e-12
        @test trace_value ≈ S7_TRACE_EXACT atol = 1e-12
        @test state_value < trace_value < 0
    end

    @testset "Dense state relaxation (S6)" begin
        oracle = expectations_oracle(SURVEY_COMPARATIVE_EXPECTATIONS_PATH, "S6_state_dense_order3")
        result = cs_nctssos(
            survey_state_problem(),
            SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )

        @test result.objective ≈ oracle.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Dense trace relaxation (S7)" begin
        oracle = expectations_oracle(SURVEY_COMPARATIVE_EXPECTATIONS_PATH, "S7_trace_dense_order3")
        result = cs_nctssos(
            survey_trace_problem(),
            SolverConfig(
                optimizer=SOLVER,
                order=3,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )

        @test result.objective ≈ oracle.opt atol = 1e-5
        @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
        @test result.objective > S6_STATE_EXACT
    end
end
