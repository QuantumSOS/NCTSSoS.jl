# test/correlated_sparsity/graph_and_cliques.jl

@testset "Graph and Cliques" begin
    @testset "get_correlative_graph" begin
        @testset "Ring graph (n=4)" begin
            expected = correlated_structure_case("get_correlative_graph_ring_n4")
            expected_neighbors = json_int_vec_vec(expected["adjacency_by_position"])
            n = 4
            reg, (x,) = create_noncommutative_variables([("x", 1:n)])
            x_idx = var_indices(x)
            objective = sum(x[i] * x[mod1(i + 1, n)] for i = 1:n)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            for i in 1:n
                @test sort(adjacency[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
            end
        end

        @testset "Chain graph (n=3)" begin
            expected = correlated_structure_case("get_correlative_graph_chain_n3")
            expected_neighbors = json_int_vec_vec(expected["adjacency_by_position"])
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            x_idx = var_indices(x)
            objective = nc_chain_objective(x)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            for i in eachindex(x_idx)
                @test sort(adjacency[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
            end
        end

        @testset "Dense graph (n=10)" begin
            expected = correlated_structure_case("get_correlative_graph_dense_n10")
            expected_neighbors = json_int_vec_vec(expected["adjacency_by_position"])
            n = 10
            reg, (x,) = create_noncommutative_variables([("x", 1:n)])
            x_idx = var_indices(x)
            objective = nc_large_scale_objective(x)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            for i in 1:n
                @test sort(adjacency[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
            end
        end

        @testset "Chain with constraints (n=3)" begin
            expected = correlated_structure_case("get_correlative_graph_chain_constraints_n3")
            expected_neighbors = json_int_vec_vec(expected["adjacency_by_position"])
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            x_idx = var_indices(x)
            objective = nc_chain_objective(x)
            constraints = vcat([1.0 - x[i]^2 for i = 1:3], [x[i] - 1.0 / 3 for i = 1:3])

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, constraints)
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            for i in eachindex(x_idx)
                @test sort(adjacency[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
            end
        end

        @testset "Bipartite graph (x[1:3], y[1:3])" begin
            expected = correlated_structure_case("get_correlative_graph_bipartite_x3y3")
            expected_neighbors = json_int_vec_vec(expected["adjacency_by_position"])
            reg, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            x_idx = var_indices(x)
            y_idx = var_indices(y)
            objective = nc_bipartite_objective(x, y)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)
            all_idx = vcat(x_idx, y_idx)

            for i in eachindex(all_idx)
                @test sort(adjacency[all_idx[i]]) == sort(all_idx[expected_neighbors[i]])
            end
        end
    end

    @testset "clique_decomp" begin
        @testset "Ring graph (n=4)" begin
            expected = correlated_structure_case("clique_decomp_ring_n4")
            reg, (x,) = create_noncommutative_variables([("x", 1:4)])
            objective = sum(x[i] * x[mod1(i + 1, 4)] for i = 1:4)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) ==
                json_int_vec_vec(expected["no_elimination"])
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                json_int_vec_vec(expected["as_is_elimination"])
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                json_int_vec_vec(expected["mf"])
        end

        @testset "Chain graph (n=3)" begin
            expected = correlated_structure_case("clique_decomp_chain_n3")
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            objective = nc_chain_objective(x)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) ==
                json_int_vec_vec(expected["no_elimination"])
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                json_int_vec_vec(expected["as_is_elimination"])
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                json_int_vec_vec(expected["mf"])
        end

        @testset "Dense graph (n=10)" begin
            expected = correlated_structure_case("clique_decomp_dense_n10")
            reg, (x,) = create_noncommutative_variables([("x", 1:10)])
            objective = nc_large_scale_objective(x)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) ==
                json_int_vec_vec(expected["no_elimination"])
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                json_int_vec_vec(expected["as_is_elimination"])
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                json_int_vec_vec(expected["mf"])
        end

        @testset "Chain with constraints (n=3)" begin
            expected = correlated_structure_case("clique_decomp_chain_constraints_n3")
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            objective = nc_chain_objective(x)
            constraints = vcat([1.0 - x[i]^2 for i = 1:3], [x[i] - 1.0 / 3 for i = 1:3])
            graph, _, _ = get_correlative_graph(reg, objective, constraints)

            @test normalize_cliques(clique_decomp(graph, NoElimination())) ==
                json_int_vec_vec(expected["no_elimination"])
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                json_int_vec_vec(expected["as_is_elimination"])
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                json_int_vec_vec(expected["mf"])
        end

        @testset "Bipartite graph (x[1:3], y[1:3])" begin
            expected = correlated_structure_case("clique_decomp_bipartite_x3y3")
            reg, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            objective = nc_bipartite_objective(x, y)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) ==
                json_int_vec_vec(expected["no_elimination"])
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                json_int_vec_vec(expected["as_is_elimination"])
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                json_int_vec_vec(expected["mf"])
        end

        @testset "Disconnected components (MaximalElimination)" begin
            expected = correlated_structure_case("clique_decomp_disconnected_maximal")
            reg, (x,) = create_noncommutative_variables([("x", 1:4)])
            objective = 1.0 * x[1] * x[2] + 1.0 * x[3] * x[4]
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, MaximalElimination())) ==
                json_int_vec_vec(expected["maximal_elimination"])
        end
    end

    @testset "assign_constraint" begin
        @testset "Constraint partition across cliques" begin
            expected = correlated_structure_case("assign_constraint_partition_across_cliques")
            reg, (x,) = create_noncommutative_variables([("x", 1:4)])
            idx_type = eltype(indices(reg))
            x_idx = var_indices(x)
            cliques = [idx_type[x_idx[1], x_idx[2], x_idx[4]], idx_type[x_idx[2], x_idx[3], x_idx[4]]]
            constraints = [1.0 * x[1] * x[2], 1.0 * x[2] * x[3], 1.0 * x[3] * x[4], 1.0 * x[4] * x[1]]

            clique_constraints, global_constraints = assign_constraint(cliques, constraints, reg)
            @test clique_constraints == json_int_vec_vec(expected["clique_constraints"])
            @test global_constraints == json_int_vec(expected["global_constraints"])
        end

        @testset "Single clique captures all constraints" begin
            expected = correlated_structure_case("assign_constraint_single_clique")
            reg, (x,) = create_noncommutative_variables([("x", 1:2)])
            idx_type = eltype(indices(reg))
            x_idx = var_indices(x)
            g = 4.0 - x[1]^2 - x[2]^2
            h1 = x[1] * x[2] + x[2] * x[1] - 2.0
            h2 = -h1
            constraints = [g, h1, h2]
            cliques = [idx_type[x_idx[1], x_idx[2]]]

            clique_constraints, global_constraints = assign_constraint(cliques, constraints, reg)
            @test clique_constraints == json_int_vec_vec(expected["clique_constraints"])
            @test global_constraints == json_int_vec(expected["global_constraints"])
        end
    end
end
