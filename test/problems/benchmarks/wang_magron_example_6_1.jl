using Test, NCTSSoS, JuMP, COSMO

if !@isdefined(SOLVER)
    const SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-7,
        "eps_rel" => 1e-7,
    )
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const EXAMPLE_SOURCE = joinpath(
    REPO_ROOT, "docs", "src", "examples", "literate", "wang_magron_2021_example_6_1.jl"
)
const EXAMPLE_GENERATED = joinpath(
    REPO_ROOT, "docs", "src", "examples", "generated", "wang_magron_2021_example_6_1.md"
)
const DOCS_MAKE = joinpath(REPO_ROOT, "docs", "make.jl")

const TABLE_7_REFERENCE = Dict(5 => 3.113, 10 => 3.011)

flatten_sizes(sizes) = reduce(vcat, sizes)
max_block_size(result) = maximum(flatten_sizes(result.moment_matrix_sizes))

dense_config(order::Int) = SolverConfig(
    optimizer=SOLVER,
    order=order,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
)

ts_config(order::Int) = SolverConfig(
    optimizer=SOLVER,
    order=order,
    cs_algo=NoElimination(),
    ts_algo=MMD(),
)

cs_config(order::Int) = SolverConfig(
    optimizer=SOLVER,
    order=order,
    cs_algo=MF(),
    ts_algo=NoElimination(),
)

csts_config(order::Int) = SolverConfig(
    optimizer=SOLVER,
    order=order,
    cs_algo=MF(),
    ts_algo=MMD(),
)

function generalized_rosenbrock_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = Float64(n) * one(typeof(x[1]))
    for i in 2:n
        f += 100.0 * x[i - 1]^4 - 200.0 * x[i - 1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
    end
    return polyopt(f, reg)
end

function broyden_banded_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[min(n, 2)])
    f = Float64(n) * one(P)
    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        f += 4.0 * x[i] + 4.0 * x[i]^2 + 10.0 * x[i]^3 + 20.0 * x[i]^4 + 25.0 * x[i]^6
        for j in jset
            f += -2.0 * x[j] - 2.0 * x[j]^2 - 4.0 * x[i] * x[j] - 4.0 * x[i] * x[j]^2
            f += -10.0 * x[i]^3 * x[j] - 10.0 * x[i]^3 * x[j]^2
        end
        for j in jset
            for k in jset
                f += x[j] * x[k] + 2.0 * x[j]^2 * x[k] + x[j]^2 * x[k]^2
            end
        end
    end
    return polyopt(f, reg)
end

function chained_singular_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[2])
    f = zero(P)

    nc_square(a, b) = a^2 + a * b + b * a + b^2
    function nc_fourth_power(a, b)
        ab = a + b
        ab2 = ab * ab
        return ab2 * ab2
    end

    for i in 1:2:(n - 3)
        f += nc_square(x[i], 10.0 * x[i + 1])
        f += 5.0 * nc_square(x[i + 2], -1.0 * x[i + 3])
        f += nc_fourth_power(x[i + 1], -2.0 * x[i + 2])
        f += 10.0 * nc_fourth_power(x[i], -10.0 * x[i + 3])
    end
    return polyopt(f, reg)
end

function chained_wood_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[2])
    f = Float64(21 * n - 41) * one(P)
    for i in 1:2:(n - 3)
        f += -2.0 * x[i] + x[i]^2 + 100.0 * x[i]^4 - 200.0 * x[i]^2 * x[i + 1]
        f += -40.0 * x[i + 1] + 110.1 * x[i + 1]^2 + 19.8 * x[i + 1] * x[i + 3]
        f += -2.0 * x[i + 2] + x[i + 2]^2 + 90.0 * x[i + 2]^4 - 180.0 * x[i + 2]^2 * x[i + 3]
        f += -40.0 * x[i + 3] + 100.1 * x[i + 3]^2
    end
    return polyopt(f, reg)
end

function broyden_tridiagonal_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    P = typeof(x[1] + x[2])
    f = Float64(n) * one(P)

    f += 5.0 * x[1]^2 + 4.0 * x[1]^4 + 4.0 * x[2]^2 - 12.0 * x[1]^3
    f += -12.0 * x[1] * x[2] + 6.0 * x[1] + 8.0 * x[1]^2 * x[2] - 4.0 * x[2]

    for i in 2:(n - 1)
        f += 5.0 * x[i]^2 + 4.0 * x[i]^4 + x[i - 1]^2 + 4.0 * x[i + 1]^2
        f += -12.0 * x[i]^3 - 6.0 * x[i - 1] * x[i] - 12.0 * x[i] * x[i + 1]
        f += 6.0 * x[i] + 4.0 * x[i - 1] * x[i]^2 + 8.0 * x[i]^2 * x[i + 1]
        f += 4.0 * x[i - 1] * x[i + 1] - 2.0 * x[i - 1] - 4.0 * x[i + 1]
    end

    f += 5.0 * x[n]^2 + 4.0 * x[n]^4 + x[n - 1]^2 - 12.0 * x[n]^3
    f += -6.0 * x[n - 1] * x[n] + 6.0 * x[n] + 4.0 * x[n - 1] * x[n]^2 - 2.0 * x[n - 1]
    return polyopt(f, reg)
end

function constrained_broyden_banded_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = sum(1:n) do i
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        g = isempty(jset) ? 0.0 * x[1] : sum(x[j] + x[j]^2 for j in jset)
        (2.0 * x[i] + 5.0 * x[i]^3 + 1 - g)^2
    end
    ineq_constraints = [[1.0 - x[i]^2 for i in 1:n]; [x[i] - 1//3 for i in 1:n]]
    return polyopt(f, reg; ineq_constraints=ineq_constraints)
end

@testset "Wang-Magron 2021 Example 6.1" begin
    @testset "Docs Artifact" begin
        @test isfile(EXAMPLE_SOURCE)
        @test isfile(EXAMPLE_GENERATED)

        if isfile(EXAMPLE_GENERATED)
            generated = read(EXAMPLE_GENERATED, String)
            @test occursin("# Wang-Magron 2021 Example 6.1: Broyden-Banded Table 2 Row", generated)
            @test occursin("## Problem formulation", generated)
            @test occursin("## Run the sparse and dense relaxations", generated)
            @test occursin("## Compare against the published row", generated)
            @test occursin("## Interpretation", generated)
            @test occursin("`n = 20`", generated)
            @test occursin("sparse: `mb = 15`", generated)
            @test occursin("dense: `mb = 61`", generated)
            @test occursin("`|opt| <= 1e-4`", generated)
        end

        makefile = read(DOCS_MAKE, String)
        @test occursin("wang_magron_2021_example_6_1.md", makefile)
    end

    @testset "Table 2 Mini: Broyden Banded" begin
        pop = broyden_banded_problem(4)
        ts_result = cs_nctssos(pop, ts_config(3))
        dense_result = cs_nctssos(pop, dense_config(3))

        @test abs(ts_result.objective) < 1e-3
        @test abs(dense_result.objective) < 1e-3
        @test max_block_size(ts_result) < max_block_size(dense_result)
    end

    @testset "Table 3 Mini: Chained Singular" begin
        pop = chained_singular_problem(8)
        ts_result = cs_nctssos(pop, ts_config(2))
        dense_result = cs_nctssos(pop, dense_config(2))

        @test abs(ts_result.objective) < 0.1
        @test abs(dense_result.objective) < 0.1
        @test max_block_size(ts_result) < max_block_size(dense_result)
    end

    @testset "Table 4 Mini: Generalized Rosenbrock" begin
        pop = generalized_rosenbrock_problem(6)
        ts_result = cs_nctssos(pop, ts_config(2))
        dense_result = cs_nctssos(pop, dense_config(2))

        @test ts_result.objective ≈ 1.0 atol = 2e-3
        @test dense_result.objective ≈ 1.0 atol = 2e-3
        @test max_block_size(ts_result) < max_block_size(dense_result)
    end

    @testset "Table 5 Mini: Chained Wood" begin
        pop = chained_wood_problem(8)
        ts_result = cs_nctssos(pop, ts_config(2))
        dense_result = cs_nctssos(pop, dense_config(2))

        @test abs(ts_result.objective - 1.0) < 0.1
        @test dense_result.objective < 3.5
        @test max_block_size(ts_result) < max_block_size(dense_result)
    end

    @testset "Table 6 Mini: Broyden Tridiagonal" begin
        pop = broyden_tridiagonal_problem(6)
        ts_result = cs_nctssos(pop, ts_config(2))
        dense_result = cs_nctssos(pop, dense_config(2))

        @test abs(ts_result.objective) < 1e-3
        @test abs(dense_result.objective) < 1e-3
        @test max_block_size(ts_result) < max_block_size(dense_result)
    end

    @testset "Table 7 Mini: Broyden Banded Over D" begin
        pop5 = constrained_broyden_banded_problem(5)
        csts5 = cs_nctssos(pop5, csts_config(3))
        cs5 = cs_nctssos(pop5, cs_config(3))
        dense5 = cs_nctssos(pop5, dense_config(3))

        @test csts5.objective ≈ TABLE_7_REFERENCE[5] atol = 5e-3
        @test cs5.objective ≈ TABLE_7_REFERENCE[5] atol = 5e-3
        @test dense5.objective ≈ TABLE_7_REFERENCE[5] atol = 5e-3
        @test max_block_size(csts5) < max_block_size(cs5)

    end
end
