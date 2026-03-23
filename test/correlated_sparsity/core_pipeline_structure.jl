# test/correlated_sparsity/core_pipeline_structure.jl

@testset "Core Pipeline Structure" begin
    @testset "Correlative Sparsity Decomposition" begin
        @testset "Example 2 (x[1:3], y[1:3])" begin
            expected = correlated_structure_case("correlative_sparsity_example2_order3")
            max_clique_by_algo = expected["max_clique_by_algo"]
            reg, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            objective = nc_bipartite_objective(x, y)
            pop = polyopt(-objective, reg)
            order = 3

            @test maximum(length.(correlative_sparsity(pop, order, NoElimination()).cliques)) ==
                max_clique_by_algo["no_elimination"]
            @test maximum(length.(correlative_sparsity(pop, order, MF()).cliques)) ==
                max_clique_by_algo["mf"]
            @test maximum(length.(correlative_sparsity(pop, order, AsIsElimination()).cliques)) ==
                max_clique_by_algo["as_is_elimination"]
        end

        @testset "Example 1 (n=10 large scale)" begin
            expected = correlated_structure_case("correlative_sparsity_example1_n10_order3")
            max_clique_by_algo = expected["max_clique_by_algo"]
            reg, (x,) = create_noncommutative_variables([("x", 1:10)])
            objective = nc_large_scale_objective(x)
            pop = polyopt(objective, reg)
            order = 3

            @test maximum(length.(correlative_sparsity(pop, order, NoElimination()).cliques)) ==
                max_clique_by_algo["no_elimination"]
            @test maximum(length.(correlative_sparsity(pop, order, MF()).cliques)) ==
                max_clique_by_algo["mf"]
            @test maximum(length.(correlative_sparsity(pop, order, AsIsElimination()).cliques)) ==
                max_clique_by_algo["as_is_elimination"]
        end

        @testset "Constrained n=2" begin
            expected = correlated_structure_case("correlative_sparsity_constrained_n2_order2")
            reg, (x,) = create_noncommutative_variables([("x", 1:2)])
            objective = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
            g = 4.0 - x[1]^2 - x[2]^2
            h1 = x[1] * x[2] + x[2] * x[1] - 2.0
            pop = polyopt(objective, reg; ineq_constraints=[g], eq_constraints=[h1])
            order = 2

            for (algo_key, algo) in (
                ("no_elimination", NoElimination()),
                ("mf", MF()),
                ("as_is_elimination", AsIsElimination()),
            )
                algo_expected = expected[algo_key]
                sparsity = correlative_sparsity(pop, order, algo)
                @test maximum(length.(sparsity.cliques)) == algo_expected["max_clique"]
                @test length.(sparsity.clq_mom_mtx_bases) ==
                    algo_expected["moment_basis_lengths"]
                @test length.(sparsity.clq_localizing_mtx_bases[1]) ==
                    algo_expected["localizing_basis_lengths"]
            end
        end

        @testset "E8 constrained polyball (order=2)" begin
            expected = correlated_structure_case("correlative_sparsity_e8_polyball_order2")
            pop = build_e8_polyball_problem()
            order = 2

            for (algo_key, algo) in (
                ("no_elimination", NoElimination()),
                ("mf", MF()),
                ("as_is_elimination", AsIsElimination()),
            )
                algo_expected = expected[algo_key]
                sparsity = correlative_sparsity(pop, order, algo)

                @test clique_symbol_names(sparsity.registry, sparsity.cliques) ==
                    normalize_cliques(algo_expected["cliques"])
                @test isempty(sparsity.global_cons)
            end
        end
    end
end
