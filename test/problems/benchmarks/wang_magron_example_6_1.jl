using Test, NCTSSoS, JuMP, COSMO

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const EXAMPLE_SOURCE = joinpath(
    REPO_ROOT, "docs", "src", "examples", "literate", "wang_magron_2021_example_6_1.jl"
)
const EXAMPLE_GENERATED = joinpath(
    REPO_ROOT, "docs", "src", "examples", "generated", "wang_magron_2021_example_6_1.md"
)
const DOCS_MAKE = joinpath(REPO_ROOT, "docs", "make.jl")
const EXAMPLE_DATA_DIR = joinpath(
    REPO_ROOT, "docs", "src", "examples", "literate", "data"
)
const EXAMPLE_MOSEK_RESULTS = joinpath(
    EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_mosek_n20_results.json"
)

include(EXAMPLE_SOURCE)

if !@isdefined(SOLVER)
    const SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-7,
        "eps_rel" => 1e-7,
    )
end

@testset "Wang-Magron 2021 Example 6.1" begin
    @testset "Docs Artifact" begin
        @test isfile(EXAMPLE_SOURCE)
        @test isfile(EXAMPLE_GENERATED)
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_2.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_3.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_4.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_5.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_6.json"))
        @test isfile(EXAMPLE_MOSEK_RESULTS)

        makefile = read(DOCS_MAKE, String)
        @test occursin("wang_magron_2021_example_6_1.md", makefile)
    end

    @testset "Fixture Metadata And Configs" begin
        dense_row = example_6_1_reference_row(2, :dense)
        sparse_row = example_6_1_reference_row(2, :sparse)
        table_2_pop = build_table_2_problem(20)
        table_4_pop = build_table_4_problem(20)
        dense_cfg = example_6_1_build_config(2, :dense, table_2_pop; optimizer=SOLVER)
        sparse_cfg = example_6_1_build_config(2, :sparse, table_2_pop; optimizer=SOLVER)
        table_4_dense_cfg = example_6_1_build_config(4, :dense, table_4_pop; optimizer=SOLVER)
        table_4_sparse_cfg = example_6_1_build_config(4, :sparse, table_4_pop; optimizer=SOLVER)

        @test dense_row.family == "f_Bb"
        @test dense_row.mode == "dense"
        @test dense_row.n == 20
        @test dense_row.paper == 0.0
        @test dense_row.paper_mb == 61

        @test sparse_row.family == "f_Bb"
        @test sparse_row.mode == "sparse"
        @test sparse_row.n == 20
        @test sparse_row.paper == 0.0
        @test sparse_row.paper_mb == 15

        @test dense_cfg.order == 3
        @test dense_cfg.moment_basis === nothing
        @test dense_cfg.cs_algo isa NoElimination
        @test dense_cfg.ts_algo isa NoElimination

        @test sparse_cfg.order == 3
        @test sparse_cfg.moment_basis === nothing
        @test sparse_cfg.cs_algo isa NoElimination
        @test sparse_cfg.ts_algo isa MMD

        expected_table_4_basis = newton_chip_basis(table_4_pop, 2)
        @test table_4_dense_cfg.moment_basis == expected_table_4_basis
        @test table_4_dense_cfg.cs_algo isa NoElimination
        @test table_4_dense_cfg.ts_algo isa NoElimination
        @test table_4_sparse_cfg.moment_basis == expected_table_4_basis
        @test table_4_sparse_cfg.cs_algo isa MF
        @test table_4_sparse_cfg.ts_algo isa MMD
        @test length(EXAMPLE_6_1_N20_CASE_SPECS) == 10
    end

    @testset "Tracked Mosek n = 20 Rows" begin
        tracked_rows = example_6_1_verify_mosek_n20_results()

        @test length(tracked_rows) == 8
        @test all(row -> row.n == 20, tracked_rows)
        @test all(row -> row.table in 3:6, tracked_rows)
        @test all(row -> row.mb == row.paper_mb, tracked_rows)
        @test any(row -> row.table == 4 && row.mode == "dense", tracked_rows)
        @test !any(row -> row.table == 2, tracked_rows)
    end

    @testset "COSMO Smoke Cases" begin
        io = IOBuffer()
        cases = [
            (runner=test_table_4_sparse_n20, family="f_gR", mode="sparse", mb=3, atol=2e-3),
            (runner=test_table_6_sparse_n20, family="f_Bt", mode="sparse", mb=5, atol=2e-3),
        ]
        reports = NamedTuple[]

        for spec in cases
            report = spec.runner(optimizer=SOLVER, atol=spec.atol, io=io)
            push!(reports, report)
            @test report.row.family == spec.family
            @test report.row.mode == spec.mode
            @test report.row.n == 20
            @test report.row.mb == spec.mb
        end

        output = String(take!(io))
        @test occursin("sparse", output)
        @test occursin("T4", output)
        @test occursin("T6", output)
    end

    @testset "Generated Docs Mention Solved Mosek Rows" begin
        generated = read(EXAMPLE_GENERATED, String)
        @test occursin("## Recorded Mosek results (`n = 20`)", generated)
        @test occursin("wang_magron_2021_example_6_1_mosek_n20_results.json", generated)
        @test occursin("if EXAMPLE_6_1_HAS_MOSEKTOOLS", generated)
        @test occursin("test_all_n20_reference_cases()        # all 10 cases at once", generated)
        @test occursin("| T3 dense   |  4.72e-4  | -1e-4 |", generated)
        @test occursin("| T6 sparse  | -5.96e-9  |  0.0  |", generated)
        @test occursin("Table 2 expectations (`paper = 0`, `mb = 61` dense / `mb = 15` sparse) are", generated)
    end
end
