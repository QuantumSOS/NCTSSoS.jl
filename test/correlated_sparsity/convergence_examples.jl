# test/correlated_sparsity/convergence_examples.jl

function build_wang_magron_example_3_4()
    reg, (vars,) = create_noncommutative_variables([("X", 1:3)])
    x, y, z = vars
    objective = x^2 - x * y - y * x + 3.0 * y^2 -
        2.0 * x * y * x + 2.0 * x * y^2 * x -
        y * z - z * y + 6.0 * z^2 + 9.0 * x^2 * y +
        9.0 * z^2 * y - 54.0 * z * y * z +
        142.0 * z * y^2 * z
    return polyopt(objective, reg)
end

function build_wang_magron_example_3_8()
    reg, (vars,) = create_noncommutative_variables([("X", 1:2)])
    x, y = vars
    objective = 2.0 - x^2 + x * y^2 * x - y^2
    g_ball = 4.0 - x^2 - y^2
    g_link = x * y + y * x - 2.0
    return polyopt(objective, reg; ineq_constraints=[g_ball], eq_constraints=[g_link])
end

function term_summary(term_sparsities)
    (
        support_sizes = [length(ts.term_sparse_graph_supp) for ts in term_sparsities],
        block_sizes = [length.(ts.block_bases) for ts in term_sparsities],
    )
end

function next_term_summary(sparsity, clique_idx, ts_algo)
    corr = sparsity.corr_sparsity
    current = sparsity.cliques_term_sparsities[clique_idx]
    activated_supp = reduce(
        NCTSSoS.sorted_union,
        [ts.term_sparse_graph_supp for ts in current]
    )
    cons_idx = corr.clq_cons[clique_idx]
    updated = NCTSSoS.term_sparsities(
        activated_supp,
        corr.cons[cons_idx],
        corr.clq_mom_mtx_bases[clique_idx],
        corr.clq_localizing_mtx_bases[clique_idx],
        ts_algo,
    )
    return term_summary(updated)
end

@testset "Wang-Magron convergence examples" begin
    @testset "Example 3.4 stabilizes structurally after the first sparse step" begin
        pop = build_wang_magron_example_3_4()
        config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())
        sparsity = compute_sparsity(pop, config)

        k1 = term_summary(sparsity.cliques_term_sparsities[1])
        k2 = next_term_summary(sparsity, 1, MMD())

        @test k1 == k2
    end

    @testset "Example 3.8 is stable and exact with no chordal extension" begin
        pop = build_wang_magron_example_3_8()
        config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=NoElimination())

        result_k1 = cs_nctssos(pop, config)
        result_k2 = cs_nctssos_higher(pop, result_k1, config)

        @test term_summary(result_k1.sparsity.cliques_term_sparsities[1]) ==
            term_summary(result_k2.sparsity.cliques_term_sparsities[1])
        @test result_k1.objective ≈ -1.0 atol=1e-6
        @test result_k2.objective ≈ -1.0 atol=1e-6
    end
end
