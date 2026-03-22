# test/problems/benchmarks/e5_broyden_tridiagonal.jl
# E5: NC Broyden Tridiagonal Function
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Table 6); benchmark metadata mirrored in
# `test-examples.typ` as E5.
#
# The literature benchmark is a sum of Hermitian squares of tridiagonal
# residuals. We keep that factorized form so the exact optimum 0 is obvious
# and the Newton-chip basis matches the published monomial-basis counts.
# This test checks the literature rows n = 20, 60, verifies the dense
# Newton-chip basis sizes and sparse term-sparsity block sizes, and solves
# both cases with Clarabel near the exact optimum 0.

using Test, TOML, NCTSSoS, JuMP, Clarabel

const E5_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
)
const E5_OBJECTIVE_ATOL = 1e-6

function _e5_broyden_tridiagonal_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "e5_broyden_tridiagonal.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in E5 expectations TOML.")
    return data["cases"]
end

function broyden_tridiagonal_problem(n::Int)
    n >= 2 || throw(ArgumentError("Broyden tridiagonal benchmark requires n >= 2."))

    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[2])
    one_poly = one(poly_type)

    first_residual = 3.0 * x[1] - 2.0 * x[1]^2 - 2.0 * x[2] + one_poly
    objective = adjoint(first_residual) * first_residual

    for i in 2:(n - 1)
        residual = 3.0 * x[i] - 2.0 * x[i]^2 - x[i - 1] - 2.0 * x[i + 1] + one_poly
        objective += adjoint(residual) * residual
    end

    last_residual = 3.0 * x[n] - 2.0 * x[n]^2 - x[n - 1] + one_poly
    objective += adjoint(last_residual) * last_residual

    return polyopt(objective, reg)
end

function _e5_ts_block_sizes(sparsity)
    return reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end)
end

@testset "E5 NC Broyden Tridiagonal Function" begin
    cases = _e5_broyden_tridiagonal_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 60]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        pop = broyden_tridiagonal_problem(n)
        basis = newton_chip_basis(pop, 2)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        sparsity = compute_sparsity(pop, sparsity_config)
        flat_sizes = _e5_ts_block_sizes(sparsity)

        @test length(basis) == expected["dense_mb"]
        @test maximum(flat_sizes) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=E5_SOLVER,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, solve_config; dualize=false)

        @test result.objective ≈ expected["sparse_opt"] atol = E5_OBJECTIVE_ATOL
    end
end
