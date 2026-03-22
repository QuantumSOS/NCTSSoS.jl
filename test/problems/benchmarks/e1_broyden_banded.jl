# test/problems/benchmarks/e1_broyden_banded.jl
# E1: NC Broyden banded function
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative Polynomial
# Optimization" (Table 2); benchmark metadata mirrored in
# `test-examples.typ` as E1.
#
# The paper-facing sparse formulation uses a Newton-chip basis. This test keeps
# a lightweight structural check for the literature rows n = 20, 40 and
# solves both sizes with Clarabel, asserting the optimum is near 0.
# The full n = 20:20:200 structural sweep lives in
# `test/correlated_sparsity/benchmark_structure.jl`.

using Test, TOML, NCTSSoS, JuMP, Clarabel

const E1_SOLVER = optimizer_with_attributes(
    Clarabel.Optimizer,
    "verbose" => false,
)
const E1_OBJECTIVE_ATOL = 2e-5

function _e1_broyden_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "e1_broyden_banded.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in E1 expectations TOML.")
    return data["cases"]
end

function broyden_banded_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[min(2, n)])
    objective = zero(poly_type)

    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        neighbor_sum = sum(x[j] + x[j]^2 for j in jset; init=zero(poly_type))
        single_term = 2.0 * x[i] + 5.0 * x[i]^3 + 1.0 - neighbor_sum
        objective += single_term' * single_term
    end

    return polyopt(objective, reg)
end

function _ts_block_sizes(sparsity)
    return reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end)
end

@testset "E1 NC Broyden Banded Function" begin
    cases = _e1_broyden_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 40]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        pop = broyden_banded_problem(n)
        basis = newton_chip_basis(pop, 3)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        sparsity = compute_sparsity(pop, sparsity_config)
        flat_sizes = _ts_block_sizes(sparsity)

        @test length(basis) == expected["dense_mb"]
        @test maximum(flat_sizes) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=E1_SOLVER,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, solve_config; dualize=true)

        @test result.objective ≈ expected["sparse_opt"] atol = E1_OBJECTIVE_ATOL
        @test maximum(reduce(vcat, result.moment_matrix_sizes)) == expected["sparse_mb"]
    end
end
