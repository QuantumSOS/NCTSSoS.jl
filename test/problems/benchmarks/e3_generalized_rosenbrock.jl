# test/problems/benchmarks/e3_generalized_rosenbrock.jl
# E3: NC Generalized Rosenbrock Function
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Table 4); benchmark metadata mirrored in
# `test-examples.typ` as E3.
#
# This benchmark matters because its optimum is 1 rather than 0: the
# Rosenbrock objective is a constant one plus a sum of Hermitian squares.
# We keep a lightweight regression on the literature rows n = 20, 60,
# checking the dense Newton-chip basis sizes, the sparse term-sparsity block
# sizes, and that Clarabel solves both cases near the reported optimum 1.

using Test, TOML, NCTSSoS, JuMP, Clarabel

const E3_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
)
const E3_OBJECTIVE_ATOL = 2e-5

function _e3_generalized_rosenbrock_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "e3_generalized_rosenbrock.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in E3 expectations TOML.")
    return data["cases"]
end

function generalized_rosenbrock_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[min(2, n)])
    one_poly = one(poly_type)
    objective = one_poly

    for i in 2:n
        residual = x[i] - x[i - 1]^2
        shift = one_poly - x[i]
        objective += 100.0 * adjoint(residual) * residual
        objective += adjoint(shift) * shift
    end

    return polyopt(objective, reg)
end

function _e3_ts_block_sizes(sparsity)
    return reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end)
end

@testset "E3 NC Generalized Rosenbrock Function" begin
    cases = _e3_generalized_rosenbrock_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 60]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        pop = generalized_rosenbrock_problem(n)
        basis = newton_chip_basis(pop, 2)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            monomial_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        sparsity = compute_sparsity(pop, sparsity_config)
        flat_sizes = _e3_ts_block_sizes(sparsity)

        @test length(basis) == expected["dense_mb"]
        @test maximum(flat_sizes) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=E3_SOLVER,
            monomial_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, solve_config; dualize=false)

        @test result.objective ≈ expected["sparse_opt"] atol = E3_OBJECTIVE_ATOL
        @test maximum(reduce(vcat, result.moment_matrix_sizes)) == expected["sparse_mb"]
    end
end
