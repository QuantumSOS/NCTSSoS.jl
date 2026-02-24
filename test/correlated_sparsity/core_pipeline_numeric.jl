# test/correlated_sparsity/core_pipeline_numeric.jl

@testset "Core Pipeline Numeric" begin
    pop = build_nc_correlative_problem()

    # Source mapping:
    # - oracle verifier: test/oracles/nctssos_corr_sparsity.jl
    # - literature context: @wangExploitingTermSparsity2021
    @testset "CS MF (moment, order=3)" begin
        oracle = CORRELATED_PIPELINE_ORACLES.CS_d3
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "CS MF (SOS, order=3)" begin
        oracle = CORRELATED_PIPELINE_ORACLES.CS_d3
        config = SolverConfig(optimizer=SOLVER, order=3, cs_algo=MF())
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 1e-5
    end

    # Source mapping:
    # - oracle verifier: test/oracles/nctssos_corr_sparsity.jl
    # - literature context: @wangExploitingTermSparsity2021
    @testset "TS MMD (moment, order=3 + higher)" begin
        oracle = CORRELATED_PIPELINE_ORACLES.TS_d3
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=false)
        result = cs_nctssos_higher(pop, result, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "TS MMD (SOS, order=3)" begin
        oracle = CORRELATED_PIPELINE_ORACLES.TS_d3
        config = SolverConfig(optimizer=SOLVER, order=3, ts_algo=MMD())
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 1e-4
    end
end
