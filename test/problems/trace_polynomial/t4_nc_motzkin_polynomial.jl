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
# This is a local-only Mosek regression, not part of the curated CI suite in
# `test/runtests.jl`. It requires a Mosek license.
#
# Run manually with:
#   make test-local-one SCRIPT=test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl
# or
#   julia --project test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl
#
# We solve the dualized SOS side with Mosek without an artificial iteration cap.
# Clarabel 0.11.1 can report ALMOST_OPTIMAL with a finite objective here, despite
# Example 4.3 ruling out every finite SOHS bound:
# https://arxiv.org/pdf/2109.00786
# Do not accept that numerical outcome as a certificate or widen this contract
# to admit it. Mosek instead reports failure with no feasible solution.

using Test, NCTSSoS, JuMP, MosekTools

const T4_FAILURE_SOLVER = optimizer_with_attributes(
    Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,
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
    # so keep the regression on the dualized SOS path. A finite returned
    # objective is not evidence for the literature's infeasibility statement.
    outcome = try
        cs_nctssos(pop, config; dualize=true)
    catch err
        err
    end

    @test outcome isa NCTSSoS.SolverStatusError
    if outcome isa NCTSSoS.SolverStatusError
        @test _t4_is_expected_failure(outcome)
    end
end
