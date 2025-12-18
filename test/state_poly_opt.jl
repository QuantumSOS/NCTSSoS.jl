# State Polynomial Optimization Tests
#
# STATUS: SKIPPED - Solver pipeline doesn't support StatePolyOpt yet
#
# These tests are syntactically migrated to the new API but cannot run because:
# - correlative_sparsity() only accepts Polynomial, not NCStatePolynomial
# - moment_relax() only accepts PolyOpt, not StatePolyOpt
#
# See .claude/tasks/state-trace-migration/plan.md for details.
#
# To enable: Implement StatePolyOpt support in src/sparse.jl and src/moment_solver.jl

using Test

@testset "State Polynomial Optimization" begin
    @test_skip "StatePolyOpt solver support not yet implemented"
end

#=
# Original tests preserved for when StatePolyOpt support is added:

using Test, NCTSSoS, NCTSSoS.FastPolynomials

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

using COSMO
const QUICK_SOLVER = COSMO.Optimizer
using JuMP
using NCTSSoS:
    neat_dot,
    constrain_moment_matrix!,
    substitute_variables,
    NoElimination

using NCTSSoS.FastPolynomials: expval, terms, Arbitrary, get_state_basis, NCStateWord, ς, Monomial

@testset "State Polynomial Opt 7.2.0" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp =
        -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) +
        1.0 * ς(x[2] * y[2])
    # Use typed identity monomial to convert StatePolynomial to NCStatePolynomial
    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 1

    solver_config = SolverConfig(; optimizer=SOLVER, order=d)

    if haskey(ENV, "LOCAL_TESTING")
        result_mom = cs_nctssos(spop, solver_config; dualize=false)
        @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-5)
    end

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)


    @testset "Sparse" begin
        solver_config = SolverConfig(; optimizer=SOLVER, order=d, cs_algo=NoElimination(), ts_algo=MMD())

        result = cs_nctssos(spop, solver_config)

        @test result.objective ≈ -2.8284271321623202 atol = 1e-5
    end
end

@testset "State Polynomial Opt 7.2.1" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 3

    solver_config = SolverConfig(; optimizer=QUICK_SOLVER, order=d)

    if haskey(ENV, "LOCAL_TESTING")
        result_mom = cs_nctssos(spop, solver_config; dualize=false)
        @test isapprox(result_mom.objective, -4.0, atol=1e-4)
    end

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -4.0, atol=1e-4)
end

@testset "State Polynomial Opt 7.2.2" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
    cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
    sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
    spop = polyopt(sp * one(typeof(x[1])), reg)
    solver_config = SolverConfig(; optimizer=SOLVER, order=2)
    result = cs_nctssos(spop, solver_config)
    @test result.objective ≈ -5.0 atol = 1e-2
end

# NOTE: "Constrain Moment matrix" test requires internal API migration
# This test directly accesses internal functions that need additional work
# Skipping for now - can be migrated later if needed
=#
