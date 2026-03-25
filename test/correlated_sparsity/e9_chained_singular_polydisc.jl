# test/correlated_sparsity/e9_chained_singular_polydisc.jl
# E9: Chained Singular on NC Polydisc (correlative sparsity)
#
# Source: Klep, Magron & Povh (2019), Table 2.
# This benchmark belongs in the correlative-sparsity suite because the value we
# want to protect is the clique decomposition and the sparse-vs-dense hierarchy
# behavior, not just the final objective on a single end-to-end instance.
#
# The paper reports objectives only to table precision, and these dualized SOS
# values drift slightly across solver/dependency resolutions. In the repository
# environment, the stored E9 objectives reproduce exactly, but under `Pkg.test()`
# re-resolution (JuMP/MOI/CliqueTrees upgrades) we observed objective drift up
# to about 6.2e-3 while the structural oracles (`sides`, `nuniq`, clique data)
# stayed exact. Use 1e-2 so CI remains sensitive to formulation regressions
# without overfitting one exact floating-point dependency graph.
const E9_OBJECTIVE_ATOL = 1e-2

# Dense and sparse dualized SOS relaxations should agree to tight numerical
# tolerance, but Linux/COSMO runs have shown cross-formulation gaps on the order
# of a few 1e-4 even when both sides still match their reviewed objective
# oracles and exact structural invariants. Keep this stricter than the oracle
# tolerance while allowing ordinary platform/solver drift.
const E9_DENSE_SPARSE_ATOL = 1e-3

@testset "E9 chained singular on NC polydisc" begin
    @testset "Structure (order=2)" begin
        expected_n4 = correlated_structure_case("correlative_sparsity_e9_polydisc_n4_order2")
        _, x4, pop4 = build_e9_chained_singular_polydisc_problem(4)

        dense_corr = correlative_sparsity(pop4, 2, NoElimination())
        @test clique_variable_positions(x4, dense_corr.cliques) ==
            expected_n4["no_elimination"]["cliques"]
        @test isempty(dense_corr.global_cons)

        sparse_corr_n4 = correlative_sparsity(pop4, 2, MF())
        @test clique_variable_positions(x4, sparse_corr_n4.cliques) ==
            expected_n4["mf"]["cliques"]
        @test isempty(sparse_corr_n4.global_cons)

        @testset "MF n=$n" for (n, expected_id) in (
            (8, "correlative_sparsity_e9_polydisc_n8_order2"),
            (12, "correlative_sparsity_e9_polydisc_n12_order2"),
        )
            expected = correlated_structure_case(expected_id)
            _, x, pop = build_e9_chained_singular_polydisc_problem(n)
            corr = correlative_sparsity(pop, 2, MF())

            @test clique_variable_positions(x, corr.cliques) == expected["mf"]["cliques"]
            @test isempty(corr.global_cons)
        end
    end

    # Use the dualized SOS path here because it keeps the E9 coverage cheap
    # enough for always-on CI while preserving the benchmark values.
    @testset "Dualized SOS hierarchy (order=2)" begin
        # The request here is intentionally asymmetric: dense is covered only at
        # n=4, while n=8 and n=12 stay sparse-only to keep always-on CI cheap.
        @testset "n=4 dense vs sparse" begin
            _, _, pop = build_e9_chained_singular_polydisc_problem(4)
            dense_oracle = CORRELATED_PIPELINE_ORACLES.E9_dense_n4_d2
            sparse_oracle = CORRELATED_PIPELINE_ORACLES.E9_sparse_n4_d2

            sparse_config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            sparse_result = cs_nctssos(pop, sparse_config; dualize=true)
            @test sparse_result.objective ≈ sparse_oracle.opt atol = E9_OBJECTIVE_ATOL
            @test flatten_sizes(sparse_result.moment_matrix_sizes) == sparse_oracle.sides
            @test sparse_result.n_unique_moment_matrix_elements == sparse_oracle.nuniq

            dense_config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=NoElimination(),
                ts_algo=NoElimination()
            )
            dense_result = cs_nctssos(pop, dense_config; dualize=true)
            @test dense_result.objective ≈ dense_oracle.opt atol = E9_OBJECTIVE_ATOL
            @test flatten_sizes(dense_result.moment_matrix_sizes) == dense_oracle.sides
            @test dense_result.n_unique_moment_matrix_elements == dense_oracle.nuniq

            # Table 2 reports that sparse matches dense for n ≤ 12.
            @test dense_result.objective ≈ sparse_result.objective atol = E9_DENSE_SPARSE_ATOL
        end

        @testset "Sparse MF n=$n" for (n, oracle) in (
            (8, CORRELATED_PIPELINE_ORACLES.E9_sparse_n8_d2),
            (12, CORRELATED_PIPELINE_ORACLES.E9_sparse_n12_d2),
        )
            _, _, pop = build_e9_chained_singular_polydisc_problem(n)
            config = SolverConfig(
                optimizer=SOLVER,
                order=2,
                cs_algo=MF(),
                ts_algo=NoElimination()
            )
            result = cs_nctssos(pop, config; dualize=true)

            @test result.objective ≈ oracle.opt atol = E9_OBJECTIVE_ATOL
            @test flatten_sizes(result.moment_matrix_sizes) == oracle.sides
            @test result.n_unique_moment_matrix_elements == oracle.nuniq
        end
    end
end
