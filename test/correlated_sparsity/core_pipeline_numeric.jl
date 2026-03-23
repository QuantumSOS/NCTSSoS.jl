# test/correlated_sparsity/core_pipeline_numeric.jl

@testset "Core Pipeline Numeric" begin
    pop = build_nc_correlative_problem()

    # Source mapping:
    # - oracle verifier: test/oracles/nctssos_corr_sparsity.jl
    # - literature context: @wangExploitingTermSparsity2021
    @testset "CS MF (moment, order=3)" begin
        oracle = CORRELATED_PIPELINE_ORACLES.CS_d3
        structure = correlated_structure_case("correlated_pipeline_cs_d3_numeric_structure")
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        result = cs_nctssos(pop, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(structure["sides"])
        @test result.n_unique_moment_matrix_elements ==
            structure["n_unique_moment_matrix_elements"]
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
        structure = correlated_structure_case("correlated_pipeline_ts_d3_numeric_structure")
        config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=NoElimination(),
            ts_algo=MMD()
        )
        result = cs_nctssos(pop, config; dualize=false)
        result = cs_nctssos_higher(pop, result, config; dualize=false)
        @test result.objective ≈ oracle.opt atol = 1e-5
        @test result.n_unique_moment_matrix_elements ==
            structure["n_unique_moment_matrix_elements"]
    end

    @testset "TS MMD (SOS, order=3)" begin
        oracle = CORRELATED_PIPELINE_ORACLES.TS_d3
        config = SolverConfig(optimizer=SOLVER, order=3, ts_algo=MMD())
        result = cs_nctssos(pop, config; dualize=true)
        @test result.objective ≈ oracle.opt atol = 2e-4
    end

    @testset "E8 constrained polyball (dualized SOS hierarchy)" begin
        # Source mapping:
        # - literature context: Klep-Magron-Povh 2019, Examples 5.10 and 5.15
        # - formulation note: the paper example is self-adjoint, so the test
        #   encodes f as f₁ + f₁' + f₂ + f₂'.
        pop = build_e8_polyball_problem()

        sparse_d2_oracle = CORRELATED_PIPELINE_ORACLES.E8_sparse_d2
        sparse_d3_oracle = CORRELATED_PIPELINE_ORACLES.E8_sparse_d3
        dense_d2_oracle = CORRELATED_PIPELINE_ORACLES.E8_dense_d2

        sparse_d2_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        sparse_d2_result = cs_nctssos(pop, sparse_d2_config; dualize=true)
        @test sparse_d2_result.objective ≈ sparse_d2_oracle.opt atol = 1e-6

        sparse_d3_config = SolverConfig(
            optimizer=SOLVER,
            order=3,
            cs_algo=MF(),
            ts_algo=NoElimination()
        )
        sparse_d3_result = cs_nctssos(pop, sparse_d3_config; dualize=true)
        @test sparse_d3_result.objective ≈ sparse_d3_oracle.opt atol = 1e-6

        dense_d2_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        dense_d2_result = cs_nctssos(pop, dense_d2_config; dualize=true)
        @test dense_d2_result.objective ≈ dense_d2_oracle.opt atol = 1e-6

        @test sparse_d2_result.objective < sparse_d3_result.objective
        @test sparse_d3_result.objective ≈ dense_d2_result.objective atol = 1e-6
    end
end
