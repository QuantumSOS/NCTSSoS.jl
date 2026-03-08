using Test, NCTSSoS, JuMP

using NCTSSoS:
    MaxEntangled,
    NCStatePolynomial,
    NCStateWord,
    NonCommutativeAlgebra,
    StateWord

function _wang_magron_expectation_case(id::AbstractString)
    data = TestExpectations.expectations_load("expectations/wang_magron_example_5_3.json")
    case = TestExpectations.expectations_case(data, id)
    haskey(case, "expected") || error("Missing key `expected` for case $(repr(id)).")
    return case["expected"]
end

_json_float(v) = Float64(v)

function _trace_constraint(poly::Polynomial{NonCommutativeAlgebra,T,C}, sw_id::StateWord{MaxEntangled,NonCommutativeAlgebra,T}) where {T<:Integer,C<:Number}
    coeffs = C[coef for (coef, _) in poly.terms]
    words = [NCStateWord(sw_id, mono) for (_, mono) in poly.terms]
    return NCStatePolynomial(coeffs, words)
end

@testset "Wang-Magron Example 5.3 (trace)" begin
    lower_bounds = _wang_magron_expectation_case("example_5_3_lower_bounds")

    reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
    X, Y = vars

    f = 2.0 - X^2 + X * Y^2 * X - Y^2
    g_ball = 4.0 - X^2 - Y^2
    g_link = X * Y + Y * X - 2.0

    sw_id = one(StateWord{MaxEntangled,NonCommutativeAlgebra,eltype(X.word)})
    tpop = polyopt(
        tr(f) * one(typeof(X)),
        reg;
        ineq_constraints=[_trace_constraint(g_ball, sw_id)],
        eq_constraints=[_trace_constraint(g_link, sw_id)]
    )

    @testset "Lower Bounds" begin
        dense_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        dense_result = cs_nctssos(tpop, dense_config; dualize=false)
        @test dense_result.objective ≈ _json_float(lower_bounds["mu_2"]) atol = 1e-6

        sparse_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        sparse_result = cs_nctssos(tpop, sparse_config; dualize=false)
        @test sparse_result.objective ≈ _json_float(lower_bounds["mu_ts_2_1"]) atol = 1e-6
    end
end
