# test/problems/state_polynomial/state_polynomial_structural.jl
# Structural tests: State polynomial term sparsity without solving

using Test, NCTSSoS

include(joinpath(@__DIR__, "..", "..", "SparsityAsserts.jl"))

const STATE_EX_7_2_3_STRUCT_ORACLES = (
    Dense_d2 = (sides=[49], nuniq=233),
    TS_d2    = (nuniq=41,),
)

@testset "State Polynomial Example 7.2.3 (structural)" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    sp = -1.0 * ς(x[2]) - 1.0 * ς(y[1]) - 1.0 * ς(y[2]) +
         1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[2] * y[1]) -
         1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[2]) +
         1.0 * ς(x[1]) * ς(y[1]) + 1.0 * ς(x[2]) * ς(y[1]) +
         1.0 * ς(x[2]) * ς(y[2]) +
         1.0 * ς(x[1]) * ς(x[1]) + 1.0 * ς(y[2]) * ς(y[2])

    pop = polyopt(sp * one(typeof(x[1])), reg)

    dense_cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=NoElimination())
    ts_cfg = SolverConfig(optimizer=nothing, order=2, cs_algo=NoElimination(), ts_algo=MMD())

    dense_sp = compute_sparsity(pop, dense_cfg)
    ts_sp = compute_sparsity(pop, ts_cfg)

    assert_term_sparsity_alignment(dense_sp)
    assert_term_sparsity_alignment(ts_sp)

    @test flatten_moment_block_sizes(dense_sp) == STATE_EX_7_2_3_STRUCT_ORACLES.Dense_d2.sides
    @test max_moment_block_size(ts_sp) < max_moment_block_size(dense_sp)

    dense_mp = NCTSSoS.moment_relax(pop, dense_sp.corr_sparsity, dense_sp.cliques_term_sparsities)
    ts_mp = NCTSSoS.moment_relax(pop, ts_sp.corr_sparsity, ts_sp.cliques_term_sparsities)

    @test dense_mp.n_unique_moment_matrix_elements == STATE_EX_7_2_3_STRUCT_ORACLES.Dense_d2.nuniq
    @test ts_mp.n_unique_moment_matrix_elements == STATE_EX_7_2_3_STRUCT_ORACLES.TS_d2.nuniq
    @test ts_mp.n_unique_moment_matrix_elements < dense_mp.n_unique_moment_matrix_elements
end
