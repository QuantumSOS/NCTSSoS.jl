# test/problems/nc_polynomial/nc_example1_structural.jl
# Structural tests: NC Example 1 - Term sparsity without solving
#
# Purpose: verify term sparsity block structure + unique moment indexing are correct
# without relying on any SDP optimizer.

using Test, NCTSSoS

include(joinpath(@__DIR__, "..", "..", "SparsityAsserts.jl"))

const NC_EXAMPLE1_STRUCT_ORACLES = (
    Dense_d2 = (sides=[13], nuniq=73),
    TS_d2    = (sides=[1, 1, 2, 2, 2, 2, 3, 3, 3, 4], nuniq=21),
)

@testset "NC Example 1 (structural, unconstrained)" begin
    n = 3
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
        2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
        x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
        54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

    pop = polyopt(f, reg)

    @testset "Dense (order=2)" begin
        oracle = NC_EXAMPLE1_STRUCT_ORACLES.Dense_d2
        cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=NoElimination())
        sparsity = compute_sparsity(pop, cfg)
        assert_term_sparsity_alignment(sparsity)

        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        @test flatten_moment_block_sizes(sparsity) == oracle.sides
        @test mp.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Term Sparsity MMD (order=2)" begin
        oracle = NC_EXAMPLE1_STRUCT_ORACLES.TS_d2
        cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=MMD())
        sparsity = compute_sparsity(pop, cfg)
        assert_term_sparsity_alignment(sparsity)

        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        @test sort(flatten_moment_block_sizes(sparsity)) == sort(oracle.sides)
        @test mp.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Dense vs TS monotonicity" begin
        dense_cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=NoElimination())
        ts_cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=MMD())

        dense_sp = compute_sparsity(pop, dense_cfg)
        ts_sp = compute_sparsity(pop, ts_cfg)

        dense_mp = NCTSSoS.moment_relax(pop, dense_sp.corr_sparsity, dense_sp.cliques_term_sparsities)
        ts_mp = NCTSSoS.moment_relax(pop, ts_sp.corr_sparsity, ts_sp.cliques_term_sparsities)

        @test ts_mp.n_unique_moment_matrix_elements ≤ dense_mp.n_unique_moment_matrix_elements
        @test max_moment_block_size(ts_sp) ≤ max_moment_block_size(dense_sp)
    end
end

