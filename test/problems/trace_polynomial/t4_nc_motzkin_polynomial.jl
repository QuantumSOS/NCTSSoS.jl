# test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl
# T4: NC Motzkin Polynomial (SOHS Failure)
#
# Source: Bhardwaj, Klep & Magron, "Noncommutative Polynomial
# Optimization" (Example 4.3); benchmark metadata mirrored in
# `test-examples.typ` as T4.
#
# The nc Motzkin polynomial
#
#     M_nc(x, y) = x y^4 x + y x^4 y - 3 x y^2 x + 1
#
# is trace-positive, so the true tracial minimum is 0. However, the standard
# denominator-free SOHS hierarchy has value -Inf at every order.
#
# This is a local-only regression, not part of the curated CI suite in
# `test/runtests.jl`. Clarabel/MOI combinations are not stable enough yet to
# make the exact failure manifestation deterministic across the package-test
# environment, but the benchmark remains useful as a manual check of the
# dualized SOHS failure path.
#
# Run manually with:
#   make test-local-one SCRIPT=test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl
# or
#   julia --project test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl
#
# We solve the dualized SOS side with Clarabel without relying on an artificially
# tiny iteration cap.

using Test, NCTSSoS, JuMP, Clarabel

const T4_FAILURE_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
)

function t4_nc_motzkin_problem()
    reg, (vars,) = create_noncommutative_variables([("x", 1:2)])
    x, y = vars

    motzkin = x * y^4 * x + y * x^4 * y - 3.0 * x * y^2 * x + 1.0
    return polyopt(tr(motzkin) * one(typeof(x)), reg)
end

function _t4_is_expected_failure(err::NCTSSoS.SolverStatusError)
    err.termination in (JuMP.MOI.INFEASIBLE, JuMP.MOI.INFEASIBLE_OR_UNBOUNDED, JuMP.MOI.DUAL_INFEASIBLE) &&
        return true

    err.primal == JuMP.MOI.INFEASIBILITY_CERTIFICATE && return true
    err.dual == JuMP.MOI.INFEASIBILITY_CERTIFICATE && return true

    if err.termination in (JuMP.MOI.ITERATION_LIMIT, JuMP.MOI.SLOW_PROGRESS)
        return err.primal ∉ (JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT) &&
               err.dual ∉ (JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT)
    end

    return false
end

function _t4_is_expected_noncertificate(result::NCTSSoS.PolyOptResult)
    status = termination_status(result.model)
    primal = primal_status(result.model)
    dual = dual_status(result.model)

    status in (JuMP.MOI.ITERATION_LIMIT, JuMP.MOI.SLOW_PROGRESS) || return false

    # Clarabel/MOI version combinations can either throw a solver-status error
    # or return an early-stop iterate that is clearly not a meaningful finite
    # certificate for this benchmark. Accept the latter only when it is both
    # non-optimal and obviously bogus relative to the true trace minimum 0.
    return result.objective < -1.0 ||
           primal ∉ (JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT) ||
           dual ∉ (JuMP.MOI.FEASIBLE_POINT, JuMP.MOI.NEARLY_FEASIBLE_POINT)
end

@testset "T4 NC Motzkin Polynomial" begin
    pop = t4_nc_motzkin_problem()

    @test NCTSSoS.compute_relaxation_order(pop, 0) == 3

    # Use the embedded nc-word trace basis to match the literature SOHS hierarchy,
    # not the broader compound state-word basis.
    trace_basis = NCTSSoS._embedded_trace_moment_basis(pop.registry, 3)
    @test NCTSSoS._uses_embedded_operator_basis(trace_basis)

    config = SolverConfig(
        optimizer=T4_FAILURE_SOLVER,
        moment_basis=trace_basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
    )

    # T4 is a statement about failure of the denominator-free SOHS hierarchy,
    # so keep the regression on the dualized SOS path. Depending on the
    # Clarabel/MOI combination, this can surface either as an explicit solver
    # failure or as a clearly bogus early-stop iterate.
    outcome = try
        cs_nctssos(pop, config; dualize=true)
    catch err
        err
    end

    @test outcome isa Union{NCTSSoS.SolverStatusError,NCTSSoS.PolyOptResult}
    if outcome isa NCTSSoS.SolverStatusError
        @test _t4_is_expected_failure(outcome)
    elseif outcome isa NCTSSoS.PolyOptResult
        @test _t4_is_expected_noncertificate(outcome)
    end
end
