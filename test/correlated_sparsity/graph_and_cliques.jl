# test/correlated_sparsity/graph_and_cliques.jl

@testset "Graph and Cliques" begin
    @testset "get_correlative_graph" begin
        @testset "Ring graph (n=4)" begin
            n = 4
            reg, (x,) = create_noncommutative_variables([("x", 1:n)])
            x_idx = var_indices(x)
            objective = sum(x[i] * x[mod1(i + 1, n)] for i = 1:n)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            @test sort(adjacency[x_idx[1]]) == sort([x_idx[2], x_idx[4]])
            @test sort(adjacency[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
            @test sort(adjacency[x_idx[3]]) == sort([x_idx[2], x_idx[4]])
            @test sort(adjacency[x_idx[4]]) == sort([x_idx[1], x_idx[3]])
        end

        @testset "Chain graph (n=3)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            x_idx = var_indices(x)
            objective = nc_chain_objective(x)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            @test adjacency[x_idx[1]] == [x_idx[2]]
            @test sort(adjacency[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
            @test adjacency[x_idx[3]] == [x_idx[2]]
        end

        @testset "Dense graph (n=10)" begin
            n = 10
            reg, (x,) = create_noncommutative_variables([("x", 1:n)])
            x_idx = var_indices(x)
            objective = nc_large_scale_objective(x)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)
            expected_neighbors = [
                [2, 3, 4, 5, 6, 7],
                [1, 3, 4, 5, 6, 7, 8],
                [1, 2, 4, 5, 6, 7, 8, 9],
                [1, 2, 3, 5, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 6, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 7, 8, 9, 10],
                [1, 2, 3, 4, 5, 6, 8, 9, 10],
                [2, 3, 4, 5, 6, 7, 9, 10],
                [3, 4, 5, 6, 7, 8, 10],
                [4, 5, 6, 7, 8, 9],
            ]

            for i in 1:n
                @test sort(adjacency[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
            end
        end

        @testset "Chain with constraints (n=3)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            x_idx = var_indices(x)
            objective = nc_chain_objective(x)
            constraints = vcat([1.0 - x[i]^2 for i = 1:3], [x[i] - 1.0 / 3 for i = 1:3])

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, constraints)
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            @test adjacency[x_idx[1]] == [x_idx[2]]
            @test sort(adjacency[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
            @test adjacency[x_idx[3]] == [x_idx[2]]
        end

        @testset "Bipartite graph (x[1:3], y[1:3])" begin
            reg, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            x_idx = var_indices(x)
            y_idx = var_indices(y)
            objective = nc_bipartite_objective(x, y)

            graph, sorted_indices, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            adjacency = graph_adjacency_by_index(graph, sorted_indices)

            @test sort(adjacency[x_idx[1]]) == sort([y_idx[1], y_idx[2], y_idx[3]])
            @test sort(adjacency[x_idx[2]]) == sort([y_idx[1], y_idx[2], y_idx[3]])
            @test sort(adjacency[x_idx[3]]) == sort([y_idx[1], y_idx[2]])
            @test sort(adjacency[y_idx[1]]) == sort([x_idx[1], x_idx[2], x_idx[3]])
            @test sort(adjacency[y_idx[2]]) == sort([x_idx[1], x_idx[2], x_idx[3]])
            @test sort(adjacency[y_idx[3]]) == sort([x_idx[1], x_idx[2]])
        end
    end

    @testset "clique_decomp" begin
        @testset "Ring graph (n=4)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:4)])
            objective = sum(x[i] * x[mod1(i + 1, 4)] for i = 1:4)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) == [collect(1:4)]
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                [[1, 2], [1, 4], [2, 3], [3, 4]]
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                [[1, 2, 4], [2, 3, 4]]
        end

        @testset "Chain graph (n=3)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            objective = nc_chain_objective(x)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) == [collect(1:3)]
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) == [[1, 2], [2, 3]]
            @test normalize_cliques(clique_decomp(graph, MF())) == [[1, 2], [2, 3]]
        end

        @testset "Dense graph (n=10)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:10)])
            objective = nc_large_scale_objective(x)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])
            expected = [
                [1, 2, 3, 4, 5, 6, 7],
                [2, 3, 4, 5, 6, 7, 8],
                [3, 4, 5, 6, 7, 8, 9],
                [4, 5, 6, 7, 8, 9, 10],
            ]

            @test normalize_cliques(clique_decomp(graph, NoElimination())) == [collect(1:10)]
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) == expected
            @test normalize_cliques(clique_decomp(graph, MF())) == expected
        end

        @testset "Chain with constraints (n=3)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:3)])
            objective = nc_chain_objective(x)
            constraints = vcat([1.0 - x[i]^2 for i = 1:3], [x[i] - 1.0 / 3 for i = 1:3])
            graph, _, _ = get_correlative_graph(reg, objective, constraints)

            @test normalize_cliques(clique_decomp(graph, NoElimination())) == [collect(1:3)]
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) == [[1, 2], [2, 3]]
            @test normalize_cliques(clique_decomp(graph, MF())) == [[1, 2], [2, 3]]
        end

        @testset "Bipartite graph (x[1:3], y[1:3])" begin
            reg, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            objective = nc_bipartite_objective(x, y)
            graph, _, _ = get_correlative_graph(reg, objective, typeof(objective)[])

            @test normalize_cliques(clique_decomp(graph, NoElimination())) == [collect(1:6)]
            @test normalize_cliques(clique_decomp(graph, AsIsElimination())) ==
                [[1, 4], [1, 5], [1, 6], [2, 4], [2, 5], [2, 6], [3, 4], [3, 5]]
            @test normalize_cliques(clique_decomp(graph, MF())) ==
                [[1, 2, 4, 5], [1, 2, 6], [3, 4, 5]]
        end
    end

    @testset "assign_constraint" begin
        @testset "Constraint partition across cliques" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:4)])
            idx_type = eltype(indices(reg))
            x_idx = var_indices(x)
            cliques = [idx_type[x_idx[1], x_idx[2], x_idx[4]], idx_type[x_idx[2], x_idx[3], x_idx[4]]]
            constraints = [1.0 * x[1] * x[2], 1.0 * x[2] * x[3], 1.0 * x[3] * x[4], 1.0 * x[4] * x[1]]

            clique_constraints, global_constraints = assign_constraint(cliques, constraints, reg)
            @test clique_constraints == [[1, 4], [2, 3]]
            @test global_constraints == Int[]
        end

        @testset "Single clique captures all constraints" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:2)])
            idx_type = eltype(indices(reg))
            x_idx = var_indices(x)
            g = 4.0 - x[1]^2 - x[2]^2
            h1 = x[1] * x[2] + x[2] * x[1] - 2.0
            h2 = -h1
            constraints = [g, h1, h2]
            cliques = [idx_type[x_idx[1], x_idx[2]]]

            clique_constraints, global_constraints = assign_constraint(cliques, constraints, reg)
            @test clique_constraints == [[1, 2, 3]]
            @test global_constraints == Int[]
        end
    end
end
