# =============================================================================
# I_3322 Bell Inequality Tests
# =============================================================================
# Consolidates all I_3322-related polynomial optimization tests:
#   - Standard I_3322 with moment/SOS methods
#   - Sparsity variants (dense, correlative, term)
#
# The I_3322 inequality uses 3 measurements per party with 2 outcomes each.
# Expected optimal value: ≈ -0.2509 (quantum bound)
# =============================================================================

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__DIR__), "..", "standalone_setup.jl"))

# Load oracle values
include(joinpath(dirname(@__DIR__), "..", "oracles", "results", "i3322_oracles.jl"))

# Helper: flatten moment_matrix_sizes for comparison with oracle
flatten_sizes(sizes) = reduce(vcat, sizes)

const I3322_EXPECTED = -0.2508753049688358

@testset "I_3322 Bell Inequality" begin

    # =========================================================================
    # Problem Setup Helper
    # =========================================================================
    function create_i3322_problem()
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) +
            x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) -
            x[1] - 2.0 * y[1] - y[2]
        pop = polyopt(-f, reg)
        return pop, reg
    end

    # =========================================================================
    # SOS Method (requires high precision - Mosek only)
    # =========================================================================
    if USE_LOCAL
        @testset "SOS Method (order=3)" begin
            pop, _ = create_i3322_problem()
            config = SolverConfig(optimizer=SOLVER, order=3)
            result = cs_nctssos(pop, config)
            @test result.objective ≈ I3322_EXPECTED atol = 1e-6
        end
    end

    # =========================================================================
    # Sparsity Variants
    # =========================================================================
    @testset "Sparsity Methods (order=2)" begin
        pop, _ = create_i3322_problem()

        @testset "Dense (No Sparsity)" begin
            oracle = I3322_ORACLES["I3322_Dense_d2"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            # Tolerance relaxed to 1e-2 due to COSMO solver precision limitations
            # Mosek achieves 1e-5, but COSMO requires looser tolerance for this problem
            @test result.objective ≈ oracle.opt atol = 1e-2
            @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end

        @testset "Correlative Sparsity (MF)" begin
            oracle = I3322_ORACLES["I3322_CS_d2"]
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config)
            # Tolerance relaxed to 1e-2 due to COSMO solver precision limitations
            # Mosek achieves 1e-5, but COSMO requires looser tolerance for this problem
            @test result.objective ≈ oracle.opt atol = 1e-2
            @test sort(flatten_sizes(result.moment_matrix_sizes)) == sort(oracle.sides)
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end
    end

    # =========================================================================
    # Correlative Sparsity Algorithm Comparison
    # =========================================================================
    @testset "Correlative Sparsity Algorithms" begin
        pop, _ = create_i3322_problem()

        oracle_dense = I3322_ORACLES["I3322_Dense_d2"]
        oracle_cs = I3322_ORACLES["I3322_CS_d2"]

        dense_config = SolverConfig(
            optimizer=SOLVER,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination()
        )
        dense_result = cs_nctssos(pop, dense_config)

        @testset "NoElimination (dense)" begin
            # Tolerance relaxed to 1e-2 due to COSMO solver precision limitations
            @test dense_result.objective ≈ oracle_dense.opt atol = 1e-2
            @test flatten_sizes(dense_result.moment_matrix_sizes) == oracle_dense.sides
            @test dense_result.n_unique_moment_matrix_elements == oracle_dense.nuniq
        end

        @testset "MF" begin
            mf_config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            mf_result = cs_nctssos(pop, mf_config)
            # Tolerance relaxed to 1e-2 due to COSMO solver precision limitations
            @test mf_result.objective ≈ oracle_cs.opt atol = 1e-2
            @test sort(flatten_sizes(mf_result.moment_matrix_sizes)) == sort(oracle_cs.sides)
            @test mf_result.n_unique_moment_matrix_elements == oracle_cs.nuniq
        end

        @testset "AsIsElimination" begin
            asis_config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=AsIsElimination(),
                ts_algo=NoElimination()
            )
            asis_result = cs_nctssos(pop, asis_config)
            # AsIsElimination can give looser (more negative) bounds than dense
            @test asis_result.objective <= dense_result.objective + 1e-6
        end
    end

    # =========================================================================
    # Combined Sparsity (requires Mosek for precision)
    # =========================================================================
    if USE_LOCAL
        @testset "Combined Sparsity (order=3)" begin
            pop, _ = create_i3322_problem()

            sparsity_configs = [
                ("Dense", NoElimination(), NoElimination(), I3322_EXPECTED),
                ("MF + MMD", MF(), MMD(), -0.9999999892255513),
                ("MF + MaximalElimination", MF(), MaximalElimination(), -0.3507010331201541),
            ]

            for (name, cs_algo, ts_algo, expected) in sparsity_configs
                @testset "$name" begin
                    config = SolverConfig(
                        optimizer=SOLVER,
                        order=3,
                        cs_algo=cs_algo,
                        ts_algo=ts_algo
                    )
                    result = cs_nctssos(pop, config)
                    @test result.objective ≈ expected atol = 1e-5
                end
            end
        end

        # =========================================================================
        # Combined CS+TS (order=4) - Verified against NCTSSOS oracle
        # =========================================================================
        # NOTE: CS+TS does NOT converge to the optimal -0.2509 even at high order.
        # This is similar to CHSH - combined sparsity loses cross-clique constraints.
        # 
        # We only verify the objective value matches NCTSSOS. Block structure differs
        # due to internal representation choices (NCTSSoS.jl uses different indexing).
        # =========================================================================
        @testset "CS+TS (order=4) - KNOWN LOOSE BOUND" begin
            pop, _ = create_i3322_problem()
            oracle = I3322_ORACLES["I3322_CS_TS_d4"]

            config = SolverConfig(
                optimizer=SOLVER,
                order=4,
                cs_algo=MF(),
                ts_algo=MMD()
            )
            result = cs_nctssos(pop, config)

            # Objective is ≈-1.0, NOT -0.2509 (this is expected due to CS+TS limitation)
            @test result.objective ≈ oracle.opt atol = 1e-5
        end
    end
end
