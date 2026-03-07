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

function generalized_rosenbrock_problem(n::Int)
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    f = Float64(n) * one(typeof(x[1]))
    for i in 2:n
        f += 100.0 * x[i - 1]^4 - 200.0 * x[i - 1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
    end
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
            @test occursin("## Meaning", generated)
            @test occursin("Table 4", generated)
            @test occursin("Table 7", generated)
        end

        makefile = read(DOCS_MAKE, String)
        @test occursin("wang_magron_2021_example_6_1.md", makefile)
    end

    @testset "Unconstrained Generalized Rosenbrock" begin
        pop = generalized_rosenbrock_problem(10)

        cs_result = cs_nctssos(pop, SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=NoElimination(),
        ))
        csts_result = cs_nctssos(pop, SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=MMD(),
        ))

        @test cs_result.objective ≈ 1.0 atol = 1e-3
        @test csts_result.objective ≈ 1.0 atol = 1e-3
        @test csts_result.n_unique_moment_matrix_elements < cs_result.n_unique_moment_matrix_elements
    end

    @testset "Constrained Broyden Banded Over D" begin
        pop = constrained_broyden_banded_problem(5)
        result = cs_nctssos(pop, SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=MF(),
            ts_algo=MMD(),
        ))

        @test result.objective ≈ 3.113 atol = 5e-3
    end
end
