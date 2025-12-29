# State Polynomial Optimization Tests
# ====================================
# Tests optimization over state polynomials using the ς() operator.
# Includes:
#   - Basic state polynomial tests
#   - Direct moment solving vs SOS dualization (merged from state_moment_solver.jl)

using Test, NCTSSoS
using NCTSSoS:
    get_state_basis,
    neat_dot,
    NCStateWord,
    NCStatePolynomial,
    moment_relax,
    NoElimination,
    expval,
    Arbitrary

# For quick solver tests that don't need high precision
# Use the already-configured SOLVER from setup.jl
const QUICK_SOLVER = SOLVER

# =============================================================================
# Basic State Polynomial Tests
# =============================================================================

@testset "State Polynomial Opt 7.2.0" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    # Test unary minus on StateWord (uses new -sw syntax)
    sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
    # Use typed identity monomial to convert StatePolynomial to NCStatePolynomial
    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 1

    solver_config = SolverConfig(; optimizer=SOLVER, order=d)

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)


    # Term sparsity (MMD) now works correctly for state polynomial optimization
    # after fixing init_activated_supp to include all pairwise basis products
    @testset "Sparse" begin
        solver_config = SolverConfig(; optimizer=SOLVER, order=d, cs_algo=NoElimination(), ts_algo=MMD())

        result = cs_nctssos(spop, solver_config)

        @test result.objective ≈ -2.8284271321623202 atol = 1e-5
    end
end

# Test 7.2.1: Known limitation - objectives with squared expectations <A><A>
# At order=3, the relaxation correctly gives the tight bound of -4.0.
@testset "State Polynomial Opt 7.2.1" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 3

    solver_config = SolverConfig(; optimizer=SOLVER, order=d)

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

# Test 7.2.3: Mixed state polynomial with products of expectations ς(a)*ς(b) and ς(a)²
# From NCTSSOS example with n=4 variables, testing squared expectations and mixed terms
# Expected optimal value: ≈ -3.5114802
@testset "State Polynomial Opt 7.2.3" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    # Coefficients from NCTSSOS: coe = -[1; 1; 1; -1; 1; 1; 1; -1; -1; -1; -1; -1]
    # After negation: [-1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1]
    # supp = [[[2]], [[3]], [[4]], [[1;3]], [[2;3]], [[1;4]], [[2;4]], 
    #         [[1], [3]], [[2], [3]], [[2], [4]], [[1], [1]], [[4], [4]]]
    sp = -1.0 * ς(x[2]) - 1.0 * ς(y[1]) - 1.0 * ς(y[2]) +    # linear terms
         1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[2] * y[1]) -        # mixed monomial terms
         1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[2]) +        # more mixed
         1.0 * ς(x[1]) * ς(y[1]) + 1.0 * ς(x[2]) * ς(y[1]) +  # products of expectations
         1.0 * ς(x[2]) * ς(y[2]) +                             # more products
         1.0 * ς(x[1]) * ς(x[1]) + 1.0 * ς(y[2]) * ς(y[2])    # squared expectations

    spop = polyopt(sp * one(typeof(x[1])), reg)
    solver_config = SolverConfig(; optimizer=SOLVER, order=2)
    result = cs_nctssos(spop, solver_config)
    @test result.objective ≈ -3.5114802 atol = 1e-2

    @testset "Sparse" begin
        solver_config_sparse = SolverConfig(; optimizer=SOLVER, order=2, cs_algo=NoElimination(), ts_algo=MMD())
        result_sparse = cs_nctssos(spop, solver_config_sparse)
        @test result_sparse.objective ≈ -3.5114802 atol = 1e-2
    end
end

# =============================================================================
# Direct Moment Solving vs SOS Dualization (merged from state_moment_solver.jl)
# =============================================================================

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

    @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-5)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)
    # Both methods should give similar results
    @test isapprox(result_mom.objective, result_sos.objective, atol=1e-5)
end

# Test 7.2.1: Objectives with squared expectations <A><A>
# This test requires order=3 to get the tight bound of -4.0.
# The complete basis generation now handles compound StateWords correctly.
@testset "State Polynomial Opt 7.2.1 - Order 3 Tight Bound" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 3
    solver_config = SolverConfig(; optimizer = QUICK_SOLVER, order = d)

    # This test now passes with the complete basis generation
    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -4.0, atol=1e-5)
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

# Test 7.2.3: Mixed state polynomial with squared expectations - Direct vs SOS
# Tests both moment solving and SOS dualization for the 7.2.3 problem
@testset "State Polynomial Opt 7.2.3 - Direct vs SOS" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp = -1.0 * ς(x[2]) - 1.0 * ς(y[1]) - 1.0 * ς(y[2]) +
         1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[2] * y[1]) -
         1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[2]) +
         1.0 * ς(x[1]) * ς(y[1]) + 1.0 * ς(x[2]) * ς(y[1]) +
         1.0 * ς(x[2]) * ς(y[2]) +
         1.0 * ς(x[1]) * ς(x[1]) + 1.0 * ς(y[2]) * ς(y[2])

    spop = polyopt(sp * one(typeof(x[1])), reg)

    solver_config = SolverConfig(; optimizer = SOLVER, order = 2)

    result_mom = cs_nctssos(spop, solver_config; dualize=false)
    result_sos = cs_nctssos(spop, solver_config; dualize=true)

    @test isapprox(result_mom.objective, -3.5114802, atol=1e-2)
    @test isapprox(result_sos.objective, -3.5114802, atol=1e-2)
    @test isapprox(result_mom.objective, result_sos.objective, atol=1e-3)
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

    @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-5)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)
    @test isapprox(result_mom.objective, result_sos.objective, atol=1e-5)
end
