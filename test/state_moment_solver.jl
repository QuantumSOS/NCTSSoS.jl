# State Moment Solver Tests
# Tests for both direct moment solving (dualize=false) and SOS dualization (dualize=true)

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
    get_state_basis,
    neat_dot,
    NCStateWord,
    NCStatePolynomial,
    moment_relax,
    NoElimination


@testset "State Polynomial Opt 7.2.0 - Direct vs SOS" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp =
        -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) +
        1.0 * ς(x[2] * y[2])
    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 1

    solver_config = SolverConfig(; optimizer = SOLVER, order = d)

    # Test both direct moment solving and SOS dualization
    result_mom = cs_nctssos(spop, solver_config; dualize=false)
    result_sos = cs_nctssos(spop, solver_config; dualize=true)

    @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-4)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-4)
    # Both methods should give similar results
    @test isapprox(result_mom.objective, result_sos.objective, atol=1e-3)
end

# Test 7.2.1: Objectives with squared expectations <A><A>
# This test requires order=3 to get the tight bound of -4.0.
# The complete basis generation now handles compound StateWords correctly.
@testset "State Polynomial Opt 7.2.1" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 3
    solver_config = SolverConfig(; optimizer = QUICK_SOLVER, order = d)

    # This test now passes with the complete basis generation
    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -4.0, atol = 1e-4)
end

# Test 7.2.2 with covariance expression: cov(a,b) = <xy> - <x><y>
# This contains compound StateWords <x><y> which is a known limitation.
# The test works because the optimal solution happens to be achievable with the
# single-expectation constraints.
@testset "State Polynomial Opt 7.2.2 - Covariance" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
    cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
    sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)

    spop = polyopt(sp*one(typeof(x[1])), reg)

    solver_config = SolverConfig(; optimizer = SOLVER, order = 2)

    result = cs_nctssos(spop, solver_config)
    @test result.objective ≈ -5.0 atol = 1e-2
end

# Test direct moment solving for a simple single-expectation problem
@testset "State Polynomial Direct Moment Solve" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    # Simple CHSH-like objective with only single expectations
    sp = -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) + 1.0 * ς(x[2] * y[2])
    spop = polyopt(sp * one(typeof(x[1])), reg)

    solver_config = SolverConfig(; optimizer = SOLVER, order = 1)

    # Test direct moment solving
    result_mom = cs_nctssos(spop, solver_config; dualize=false)
    result_sos = cs_nctssos(spop, solver_config; dualize=true)

    @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-4)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-4)
    @test isapprox(result_mom.objective, result_sos.objective, atol=1e-3)
end
