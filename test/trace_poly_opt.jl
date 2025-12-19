# Trace Polynomial Optimization Tests

using Test, NCTSSoS
using NCTSSoS.FastPolynomials: tr, Monomial

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

# Example 6.1 involves compound StateWords from products of traces (tr(xy) * tr(z))
# This requires compound StateWord support which is a known limitation.
# See .claude/tasks/statepolyopt-solver/context.md for details.
@testset "Example 6.1" begin
    reg, (x,) = create_projector_variables([("x", 1:3)])

    p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(typeof(x[1]))

    tpop = polyopt(p, reg)

    solver_config = SolverConfig(; optimizer=SOLVER, order=2)

    result = cs_nctssos(tpop, solver_config)

    @test_skip result.objective ≈ -0.046717378455438933 atol = 1e-6

    if haskey(ENV, "LOCAL_TESTING")
        solver_config = SolverConfig(; optimizer=SOLVER, order=3)

        result = cs_nctssos(tpop, solver_config)

        @test_skip result.objective ≈ -0.03124998978001017 atol = 1e-6
    end
end

# Note: Term sparsity (MaximalElimination) doesn't work correctly with state polynomial optimization
# Using NoElimination for ts_algo instead
@testset "Example 6.2.0" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

    p = -1.0 * tr(x[1] * y[1]) - 1.0 * tr(x[1] * y[2]) - 1.0 * tr(x[2] * y[1]) + 1.0 * tr(x[2] * y[2])

    tpop = polyopt(p * one(typeof(x[1])), reg)

    solver_config = SolverConfig(; optimizer=SOLVER, order=1)

    result = cs_nctssos(tpop, solver_config)

    @test result.objective ≈ -2.8284271157283083 atol = 1e-5
end

# Example 6.2.1 involves squared trace expressions (tr(xy) * tr(xy))
# This is the same compound StateWord limitation as test 7.2.1
@testset "Example 6.2.1" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])

    p = (1.0 * tr(x[1] * y[2]) + tr(x[2] * y[1])) * (1.0 * tr(x[1] * y[2]) + tr(x[2] * y[1])) + (1.0 * tr(x[1] * y[1]) - tr(x[2] * y[2])) * (1.0 * tr(x[1] * y[1]) - tr(x[2] * y[2]))

    tpop = polyopt((-1.0 * p) * one(typeof(x[1])), reg)

    solver_config = SolverConfig(; optimizer=SOLVER, order=2)

    result = cs_nctssos(tpop, solver_config)

    @test_skip result.objective ≈ -4.000000007460838 atol = 1e-5
end

@testset "Example 6.2.2" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
    cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
    p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
    tpop = polyopt(p * one(typeof(x[1])), reg)
    solver_config = SolverConfig(; optimizer=SOLVER, order=2)
    result = cs_nctssos(tpop, solver_config)
    @test result.objective ≈ -5.0 atol = 1e-5
end
