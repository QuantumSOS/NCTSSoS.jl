# State Polynomial Optimization Tests

using Test, NCTSSoS, NCTSSoS.FastPolynomials

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2))
    )
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

using COSMO
const QUICK_SOLVER = COSMO.Optimizer
using JuMP
using NCTSSoS:
    neat_dot,
    NoElimination

using NCTSSoS.FastPolynomials: expval, terms, Arbitrary, get_state_basis, NCStateWord, ς, Monomial

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
# The current basis (NCStateWords with identity StateWord) cannot generate
# compound StateWords like <A><B> through _neat_dot3. These terms are ignored
# in the optimization, leading to incorrect results for this test case.
# See .claude/tasks/statepolyopt-solver/context.md for detailed analysis.
@testset "State Polynomial Opt 7.2.1" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(typeof(x[1])), reg)

    d = 3

    solver_config = SolverConfig(; optimizer=QUICK_SOLVER, order=d)

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
