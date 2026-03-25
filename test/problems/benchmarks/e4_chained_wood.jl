# test/problems/benchmarks/e4_chained_wood.jl
# E4: NC Chained Wood Function
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Table 5); benchmark metadata mirrored in
# `test-examples.typ` as E4.
#
# The literature benchmark is 1 plus a sum of Hermitian squares of quadratic
# and linear residuals. We keep that factorized form so the exact optimum 1 is
# obvious and the Newton-chip basis matches the paper rows. This test checks
# the literature rows n = 20, 40, verifies the dense Newton-chip basis sizes
# and sparse term-sparsity block sizes, and solves both cases with Clarabel
# near the exact optimum 1.
#
# The paper reports small solver-noise deviations around 1 for some sparse
# rows, so we assert near 1 rather than reproducing a specific noisy decimal.

using Test, TOML, NCTSSoS, JuMP, Clarabel

const E4_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
)
const E4_OBJECTIVE_ATOL = 5e-4

function _e4_chained_wood_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "e4_chained_wood.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in E4 expectations TOML.")
    return data["cases"]
end

function chained_wood_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[min(2, n)])
    one_poly = one(poly_type)
    objective = one_poly

    for i in 1:2:(n - 3)
        residual1 = x[i + 1] - x[i]^2
        shift1 = one_poly - x[i]
        residual2 = x[i + 3] - x[i + 2]^2
        shift2 = one_poly - x[i + 2]
        coupling = x[i + 1] + x[i + 3] - 2.0 * one_poly
        difference = x[i + 1] - x[i + 3]

        objective += 100.0 * (adjoint(residual1) * residual1)
        objective += adjoint(shift1) * shift1
        objective += 90.0 * (adjoint(residual2) * residual2)
        objective += adjoint(shift2) * shift2
        objective += 10.0 * (adjoint(coupling) * coupling)
        objective += 0.1 * (adjoint(difference) * difference)
    end

    return polyopt(objective, reg)
end

function _e4_ts_block_sizes(sparsity)
    return reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end)
end

@testset "E4 NC Chained Wood Function" begin
    cases = _e4_chained_wood_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 40]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        pop = chained_wood_problem(n)
        basis = newton_chip_basis(pop, 2)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            monomial_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        sparsity = compute_sparsity(pop, sparsity_config)
        flat_sizes = _e4_ts_block_sizes(sparsity)

        @test length(basis) == expected["dense_mb"]
        @test maximum(flat_sizes) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=E4_SOLVER,
            monomial_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, solve_config; dualize=false)

        # The objective is 1 plus a sum of Hermitian squares, so the exact
        # optimum is 1. Test near that value rather than the paper's
        # solver-specific noise.
        @test result.objective ≈ 1.0 atol = E4_OBJECTIVE_ATOL
    end
end
