# test/problems/nc_polynomial/nc_example2_structural.jl
# Structural tests: NC Example 2 - Correlative + term sparsity without solving

using Test, NCTSSoS

include(joinpath(@__DIR__, "..", "..", "SparsityAsserts.jl"))

@testset "NC Example 2 (structural, constrained)" begin
    n = 2
    reg, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0

    pop = polyopt(f, reg; eq_constraints=[h1], ineq_constraints=[g])

    dense_cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=NoElimination())
    sparse_cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=MF(), ts_algo=MMD())

    dense_sp = compute_sparsity(pop, dense_cfg)
    sparse_sp = compute_sparsity(pop, sparse_cfg)

    assert_term_sparsity_alignment(dense_sp)
    assert_term_sparsity_alignment(sparse_sp)

    dense_mp = NCTSSoS.moment_relax(pop, dense_sp.corr_sparsity, dense_sp.cliques_term_sparsities)
    sparse_mp = NCTSSoS.moment_relax(pop, sparse_sp.corr_sparsity, sparse_sp.cliques_term_sparsities)

    @test sparse_mp.n_unique_moment_matrix_elements ≤ dense_mp.n_unique_moment_matrix_elements
    @test max_moment_block_size(sparse_sp) ≤ max_moment_block_size(dense_sp)
end

