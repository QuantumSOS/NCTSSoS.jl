# test/correlated_sparsity/literature_matrix.jl

const CORRELATED_LITERATURE_MATRIX = (
    CorrSparsity_CS_d3 = (
        citation=:wangExploitingTermSparsity2021,
        reference="Correlative + term sparsity NC examples",
        oracle_key="CorrSparsity_CS_d3",
        tolerance=1e-5,
        pinned_literature_value=false,
    ),
    CorrSparsity_TS_d3 = (
        citation=:wangExploitingTermSparsity2021,
        reference="Term-sparsity NC examples",
        oracle_key="CorrSparsity_TS_d3",
        tolerance=1e-5,
        pinned_literature_value=false,
    ),
)

const CORRELATED_LITERATURE_OMISSIONS = (
    E10_random_cubic_polyball = (
        source="Klep, Magron & Povh (2019), Table 3",
        tested=false,
        reason="No public exact coefficients, fixed seed, generator, or reference value for a specific historical E10 instance has been found, so E10 is intentionally omitted from the correlated sparsity regression suite.",
    ),
)

@testset "Literature Matrix" begin
    @test hasproperty(CORRELATED_LITERATURE_MATRIX, :CorrSparsity_CS_d3)
    @test hasproperty(CORRELATED_LITERATURE_MATRIX, :CorrSparsity_TS_d3)
    @test hasproperty(CORRELATED_LITERATURE_OMISSIONS, :E10_random_cubic_polyball)

    @test CORRELATED_LITERATURE_MATRIX.CorrSparsity_CS_d3.citation == :wangExploitingTermSparsity2021
    @test CORRELATED_LITERATURE_MATRIX.CorrSparsity_TS_d3.citation == :wangExploitingTermSparsity2021

    @test CORRELATED_LITERATURE_MATRIX.CorrSparsity_CS_d3.oracle_key == "CorrSparsity_CS_d3"
    @test CORRELATED_LITERATURE_MATRIX.CorrSparsity_TS_d3.oracle_key == "CorrSparsity_TS_d3"

    @test CORRELATED_LITERATURE_OMISSIONS.E10_random_cubic_polyball.source ==
        "Klep, Magron & Povh (2019), Table 3"
    @test !CORRELATED_LITERATURE_OMISSIONS.E10_random_cubic_polyball.tested
    @test occursin(
        "No public exact coefficients",
        CORRELATED_LITERATURE_OMISSIONS.E10_random_cubic_polyball.reason
    )

    # Placeholder policy from issue #248:
    # use @test_broken until literature-pinned numbers are explicitly quoted in refs/docs.
    # TODO(issue-248): replace with non-broken numeric checks once pinned references are added.
    @test_broken CORRELATED_LITERATURE_MATRIX.CorrSparsity_CS_d3.pinned_literature_value
    @test_broken CORRELATED_LITERATURE_MATRIX.CorrSparsity_TS_d3.pinned_literature_value
end
