# test/problems/benchmarks/e2_chained_singular.jl
# E2: NC Chained Singular Function
#
# Sources: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Table 3); Klep, Magron & Povh,
# "Sparse noncommutative polynomial optimization" (Table 1); benchmark
# metadata mirrored in `test-examples.typ` as E2.
#
# The literature benchmark is a sum of Hermitian squares of linear/quadratic
# expressions. That matters: replacing the quartic terms by naive NC fourth
# powers changes the Newton-chip basis and no longer matches the published
# monomial-basis counts. This test checks the literature rows n = 20, 60,
# verifies the dense Newton-chip basis sizes and sparse term-sparsity block
# sizes, and solves both cases with Clarabel near the reported optimum 0.
#
# Clarabel is reliable here on the direct moment formulation; the dualized SOS
# path stalls on this benchmark, so we keep the solve path explicit.

using Test, TOML, NCTSSoS, JuMP, Clarabel

const E2_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
)
const E2_OBJECTIVE_ATOL = 1e-5

function _e2_chained_singular_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "e2_chained_singular.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in E2 expectations TOML.")
    return data["cases"]
end

function chained_singular_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[min(2, n)])
    objective = zero(poly_type)

    for i in 1:2:(n - 3)
        linear1 = x[i] + 10.0 * x[i + 1]
        linear2 = x[i + 2] - x[i + 3]
        quadratic1 = x[i + 1]^2 - 4.0 * x[i + 1] * x[i + 2] + 4.0 * x[i + 2]^2
        quadratic2 = x[i]^2 - 20.0 * x[i] * x[i + 3] + 100.0 * x[i + 3]^2

        objective += linear1' * linear1
        objective += 5.0 * (linear2' * linear2)
        objective += quadratic1' * quadratic1
        objective += 10.0 * (quadratic2' * quadratic2)
    end

    return polyopt(objective, reg)
end

function _e2_ts_block_sizes(sparsity)
    return reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end)
end

@testset "E2 NC Chained Singular Function" begin
    cases = _e2_chained_singular_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 60]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        pop = chained_singular_problem(n)
        basis = newton_chip_basis(pop, 2)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        sparsity = compute_sparsity(pop, sparsity_config)
        flat_sizes = _e2_ts_block_sizes(sparsity)

        @test length(basis) == expected["dense_mb"]
        @test maximum(flat_sizes) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=E2_SOLVER,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, solve_config; dualize=false)

        # The literature reports tiny solver-noise deviations around the exact
        # optimum 0, so test near-zero rather than the paper's specific noise.
        @test abs(result.objective) <= E2_OBJECTIVE_ATOL
        @test maximum(reduce(vcat, result.moment_matrix_sizes)) == expected["sparse_mb"]
    end
end
