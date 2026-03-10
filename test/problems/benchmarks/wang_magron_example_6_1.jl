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
const EXAMPLE_DATA_DIR = joinpath(
    REPO_ROOT, "docs", "src", "examples", "literate", "data"
)

include(EXAMPLE_SOURCE)

@testset "Wang-Magron 2021 Example 6.1" begin
    @testset "Docs Artifact" begin
        @test isfile(EXAMPLE_SOURCE)
        @test isfile(EXAMPLE_GENERATED)
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_3.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_4.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_5.json"))
        @test isfile(joinpath(EXAMPLE_DATA_DIR, "wang_magron_2021_example_6_1_table_6.json"))

        makefile = read(DOCS_MAKE, String)
        @test occursin("wang_magron_2021_example_6_1.md", makefile)
    end

    @testset "Fixture-Backed Fast Sparse Cases" begin
        io = IOBuffer()

        table3 = test_table_3_sparse_n20(optimizer=SOLVER, atol=2e-3, io=io)
        @test table3.row.family == "f_cs"
        @test table3.row.n == 20
        @test table3.row.mb == 3

        table4 = test_table_4_sparse_n20(optimizer=SOLVER, atol=2e-3, io=io)
        @test table4.row.family == "f_gR"
        @test table4.row.n == 20
        @test table4.row.mb == 3

        table5 = test_table_5_sparse_n20(optimizer=SOLVER, atol=2e-3, io=io)
        @test table5.row.family == "f_cW"
        @test table5.row.n == 20
        @test table5.row.mb == 3

        table6 = test_table_6_sparse_n20(optimizer=SOLVER, atol=2e-3, io=io)
        @test table6.row.family == "f_Bt"
        @test table6.row.n == 20
        @test table6.row.mb == 5

        all_reports = test_all_fast_cases(optimizer=SOLVER, io=io)
        @test length(all_reports) == 4

        output = String(take!(io))
        @test occursin("T3", output)
        @test occursin("T4", output)
        @test occursin("T5", output)
        @test occursin("T6", output)
    end
end
