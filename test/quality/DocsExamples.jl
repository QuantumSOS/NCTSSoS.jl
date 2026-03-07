using Test

@testset "Docs Examples" begin
    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    docs_make = joinpath(repo_root, "docs", "make.jl")
    literate_path = joinpath(repo_root, "docs", "src", "examples", "literate", "wang_magron_2021_example_6_2.jl")
    generated_path = joinpath(repo_root, "docs", "src", "examples", "generated", "wang_magron_2021_example_6_2.md")

    @test isfile(literate_path)
    @test isfile(generated_path)

    make_text = read(docs_make, String)
    @test occursin("\"Wang-Magron 2021 Example 6.2\"=>\"examples/generated/wang_magron_2021_example_6_2.md\"", make_text)

    if isfile(generated_path)
        generated_text = read(generated_path, String)
        @test occursin("Example 6.2", generated_text)
        @test occursin("Broyden banded", generated_text)
        @test occursin("Broyden tridiagonal", generated_text)
        @test occursin("over D", generated_text)
        @test occursin("Tables 9-12", generated_text) || occursin("Tables 9–12", generated_text)
        @test occursin("## Meaning", generated_text)
    end
end
