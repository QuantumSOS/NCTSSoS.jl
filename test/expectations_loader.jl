using Test

@testset "Expectations JSON Fixtures" begin
    data = TestExpectations.expectations_load("expectations/chsh_simple.json")
    @test data["schema_version"] == 1

    case = TestExpectations.expectations_case(data, "Dense_d1")
    @test haskey(case, "expected")
    @test haskey(case["expected"], "objective")

    oracle = expectations_oracle("expectations/chsh_simple.json", "Dense_d1")
    @test oracle.opt isa Float64
    @test oracle.nuniq isa Int
    @test oracle.sides isa Vector{Int}
end

