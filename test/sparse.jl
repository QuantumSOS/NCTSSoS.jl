using Test, NCTSSoS, NCTSSoS.FastPolynomials
using Graphs, CliqueTrees

using NCTSSoS:
    assign_constraint,
    get_correlative_graph,
    clique_decomp,
    get_term_sparsity_graph,
    term_sparsity_graph_supp,
    correlative_sparsity,
    init_activated_supp

using NCTSSoS.FastPolynomials:
    Monomial, Polynomial, Term,
    NonCommutativeAlgebra, VariableRegistry,
    create_noncommutative_variables, variable_indices,
    monomials, neat_dot

# Helper to create NC polynomials from registry variables
function nc_poly(registry::VariableRegistry{NonCommutativeAlgebra,T}, indices::Vector{T}) where T
    if isempty(indices)
        return Polynomial{NonCommutativeAlgebra,T,Float64}([Term(1.0, Monomial{NonCommutativeAlgebra}(T[]))])
    end
    m = Monomial{NonCommutativeAlgebra}(indices)
    return Polynomial{NonCommutativeAlgebra,T,Float64}([Term(1.0, m)])
end

@testset "Correlative Sparsity without constraints" begin
    @testset "Example 2" begin
        # Create variables: x[1:3], y[1:3]
        registry, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])

        # Build objective: x[1]*(y[1]+y[2]+y[3]) + x[2]*(y[1]+y[2]-y[3]) + x[3]*(y[1]-y[2]) - x[1] - 2*y[1] - y[2]
        f = x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2.0 * y[1] - y[2]

        pop = polyopt(-f, registry)
        order = 3

        @testset "No Elimination" begin
            corr_sparsity = correlative_sparsity(pop, order, NoElimination())
            @test maximum(length.(corr_sparsity.cliques)) == 6
        end

        @testset "MF" begin
            corr_sparsity = correlative_sparsity(pop, order, MF())
            @test maximum(length.(corr_sparsity.cliques)) == 4
        end

        @testset "AsIS" begin
            corr_sparsity = correlative_sparsity(pop, order, AsIsElimination())
            @test maximum(length.(corr_sparsity.cliques)) == 2
        end
    end

    @testset "Example 1" begin
        n = 10
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])

        # Build objective
        P = typeof(x[1] + x[2])  # Get the polynomial type
        f = zero(P)
        for i = 1:n
            jset = max(1, i - 5):min(n, i + 1)
            jset = setdiff(jset, i)
            g = sum(x[j] + x[j]^2 for j in jset)
            f += (2.0 * x[i] + 5.0 * x[i]^3 + 1.0 - g)^2
        end

        pop = polyopt(f, registry)
        order = 3

        @testset "No Elimination" begin
            corr_sparsity = correlative_sparsity(pop, order, NoElimination())
            @test maximum(length.(corr_sparsity.cliques)) == 10
        end

        @testset "MF" begin
            corr_sparsity = correlative_sparsity(pop, order, MF())
            @test maximum(length.(corr_sparsity.cliques)) == 7
        end

        @testset "AsIS" begin
            corr_sparsity = correlative_sparsity(pop, order, AsIsElimination())
            @test maximum(length.(corr_sparsity.cliques)) == 7
        end
    end
end

@testset "Correlative Sparsity with constraints" begin
    n = 2
    registry, (x,) = create_noncommutative_variables([("x", 1:n)])

    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0

    pop = polyopt(f, registry; ineq_constraints=[g], eq_constraints=[h1])
    order = 2

    @testset "No Elimination" begin
        corr_sparsity = correlative_sparsity(pop, order, NoElimination())
        @test maximum(length.(corr_sparsity.cliques)) == 2
        @test length.(corr_sparsity.clq_mom_mtx_bases) == [7]
        @test length.(corr_sparsity.clq_localizing_mtx_bases[1]) == [3, 3]
    end

    @testset "MF" begin
        corr_sparsity = correlative_sparsity(pop, order, MF())
        @test maximum(length.(corr_sparsity.cliques)) == 2
        @test length.(corr_sparsity.clq_mom_mtx_bases) == [7]
        @test length.(corr_sparsity.clq_localizing_mtx_bases[1]) == [3, 3]
    end

    @testset "AsIS" begin
        corr_sparsity = correlative_sparsity(pop, order, AsIsElimination())
        @test maximum(length.(corr_sparsity.cliques)) == 2
        @test length.(corr_sparsity.clq_mom_mtx_bases) == [7]
        @test length.(corr_sparsity.clq_localizing_mtx_bases[1]) == [3, 3]
    end
end

@testset "Correlative Sparsity Components" begin
    @testset "Get Correlative Graph" begin
        # Helper functions for comparing graph structure
        # Maps graph node positions back to variable indices for comparison
        function graph_adjacency_by_index(G, sorted_indices)
            Dict(
                sorted_indices[i] => sort([sorted_indices[j] for j in neighbors(G, i)])
                for i in 1:nv(G)
            )
        end

        @testset "Ring graph (n=4)" begin
            n = 4
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])
            T = eltype(indices(registry))

            # Get actual indices for each variable
            x_idx = [x[i].word[1] for i in 1:n]

            # f = sum(x[i] * x[i+1]) for ring
            f = sum(x[i] * x[mod1(i + 1, n)] for i = 1:n)

            G, sorted_indices, idx_to_node = get_correlative_graph(registry, f, typeof(f)[])
            savegraph("example1.lgz", G)

            adj = graph_adjacency_by_index(G, sorted_indices)
            # x[1] connects to x[2], x[4]; x[2] to x[1], x[3]; etc.
            @test sort(adj[x_idx[1]]) == sort([x_idx[2], x_idx[4]])
            @test sort(adj[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
            @test sort(adj[x_idx[3]]) == sort([x_idx[2], x_idx[4]])
            @test sort(adj[x_idx[4]]) == sort([x_idx[1], x_idx[3]])
        end

        @testset "Chain graph (n=3)" begin
            n = 3
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])
            T = eltype(indices(registry))

            # Get actual indices for each variable
            x_idx = [x[i].word[1] for i in 1:n]

            f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 - 2.0 * x[1] * x[2] * x[1] +
                2.0 * x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
                6.0 * x[3]^2 +
                9.0 * x[2]^2 * x[3] +
                9.0 * x[3] * x[2]^2 - 54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

            G, sorted_indices, idx_to_node = get_correlative_graph(registry, f, typeof(f)[])
            savegraph("example2.lgz", G)

            adj = graph_adjacency_by_index(G, sorted_indices)
            @test adj[x_idx[1]] == [x_idx[2]]
            @test sort(adj[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
            @test adj[x_idx[3]] == [x_idx[2]]
        end

        @testset "Dense graph (n=10)" begin
            n = 10
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])
            T = eltype(indices(registry))

            # Get actual indices for each variable
            x_idx = [x[i].word[1] for i in 1:n]

            P = typeof(x[1] + x[2])  # Get the polynomial type
            f = zero(P)
            for i = 1:n
                jset = max(1, i - 5):min(n, i + 1)
                jset = setdiff(jset, i)
                f += (2.0 * x[i] + 5.0 * x[i]^3 + 1.0)^2
                f -= sum([
                    4.0 * x[i] * x[j] +
                    10.0 * x[i]^3 * x[j] +
                    2.0 * x[j] +
                    4.0 * x[i] * x[j]^2 +
                    10.0 * x[i]^3 * x[j]^2 +
                    2.0 * x[j]^2 for j in jset
                ])
                f += sum([
                    x[j] * x[k] + 2.0 * x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
                ])
            end
            G, sorted_indices, idx_to_node = get_correlative_graph(registry, f, typeof(f)[])
            savegraph("example3.lgz", G)

            # Check expected number of neighbors for each variable
            # x[1] connects to 6 vars, x[2] to 7, x[3] to 8, x[4-7] to 9, x[8] to 8, x[9] to 7, x[10] to 6
            expected_neighbor_counts = [6, 7, 8, 9, 9, 9, 9, 8, 7, 6]
            for i in 1:n
                node = idx_to_node[x_idx[i]]
                @test length(neighbors(G, node)) == expected_neighbor_counts[i]
            end
        end

        @testset "Chain with constraints (n=3)" begin
            n = 3
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])
            T = eltype(indices(registry))

            # Get actual indices for each variable
            x_idx = [x[i].word[1] for i in 1:n]

            f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 - 2.0 * x[1] * x[2] * x[1] +
                2.0 * x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
                6.0 * x[3]^2 +
                9.0 * x[2]^2 * x[3] +
                9.0 * x[3] * x[2]^2 - 54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]

            cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
            G, sorted_indices, idx_to_node = get_correlative_graph(registry, f, cons)
            savegraph("example4.lgz", G)

            # Chain structure: x[1]-x[2]-x[3], so x[1] and x[3] have 1 neighbor, x[2] has 2
            @test length(neighbors(G, idx_to_node[x_idx[1]])) == 1
            @test length(neighbors(G, idx_to_node[x_idx[2]])) == 2
            @test length(neighbors(G, idx_to_node[x_idx[3]])) == 1
        end

        @testset "Bipartite (x[1:3], y[1:3])" begin
            registry, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            T = eltype(indices(registry))

            # Get actual indices for each variable
            x_idx = [x[i].word[1] for i in 1:3]
            y_idx = [y[i].word[1] for i in 1:3]

            f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
                x[3] * (y[1] - y[2]) - x[1] - 2.0 * y[1] - y[2]

            G, sorted_indices, idx_to_node = get_correlative_graph(registry, f, typeof(f)[])
            savegraph("example5.lgz", G)

            # Check connectivity pattern: x[i] connects to y[j]
            # x[1] -> y[1], y[2], y[3] (3 neighbors)
            # x[2] -> y[1], y[2], y[3] (3 neighbors)
            # x[3] -> y[1], y[2] (2 neighbors)
            @test length(neighbors(G, idx_to_node[x_idx[1]])) == 3
            @test length(neighbors(G, idx_to_node[x_idx[2]])) == 3
            @test length(neighbors(G, idx_to_node[x_idx[3]])) == 2
        end

        # Clean up saved graphs
        for f in ["example1.lgz", "example2.lgz", "example3.lgz", "example4.lgz", "example5.lgz"]
            isfile(f) && rm(f)
        end
    end

    @testset "Clique Decomposition" begin
        @testset "Ring graph cliques" begin
            # Recreate example1 graph
            n = 4
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])
            f = sum(x[i] * x[mod1(i + 1, n)] for i = 1:n)
            G, _, _ = get_correlative_graph(registry, f, typeof(f)[])

            @test sort.(clique_decomp(G, NoElimination())) == [collect(1:4)]
            @test length(clique_decomp(G, AsIsElimination())) == 4
            @test length(clique_decomp(G, MF())) == 2
        end

        @testset "Chain graph cliques" begin
            n = 3
            registry, (x,) = create_noncommutative_variables([("x", 1:n)])
            f = x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 - 2.0 * x[1] * x[2] * x[1] +
                2.0 * x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
                6.0 * x[3]^2

            G, _, _ = get_correlative_graph(registry, f, typeof(f)[])

            @test sort.(clique_decomp(G, NoElimination())) == [collect(1:3)]
            @test length(clique_decomp(G, AsIsElimination())) == 2
            @test length(clique_decomp(G, MF())) == 2
        end
    end

    @testset "Assign Constraint" begin
        n = 4
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])
        T = eltype(indices(registry))

        # Get actual variable indices
        x_idx = [x[i].word[1] for i in 1:n]

        # Cliques as index vectors - use actual variable indices
        cliques = [T[x_idx[1], x_idx[2], x_idx[4]], T[x_idx[2], x_idx[3], x_idx[4]]]
        cons = [1.0 * x[1] * x[2], 1.0 * x[2] * x[3], 1.0 * x[3] * x[4], 1.0 * x[4] * x[1]]

        clq_cons, global_cons = assign_constraint(cliques, cons, registry)
        @test clq_cons == [[1, 4], [2, 3]]
        @test global_cons == Int[]

        # Test with single clique containing all variables
        n = 2
        registry2, (x2,) = create_noncommutative_variables([("x", 1:n)])
        T2 = eltype(indices(registry2))

        # Get actual variable indices
        x2_idx = [x2[i].word[1] for i in 1:n]

        g = 4.0 - x2[1]^2 - x2[2]^2
        h1 = x2[1] * x2[2] + x2[2] * x2[1] - 2.0
        h2 = -h1
        cons2 = [g, h1, h2]
        cliques2 = [T2[x2_idx[1], x2_idx[2]]]

        clq_cons2, global_cons2 = assign_constraint(cliques2, cons2, registry2)
        @test clq_cons2 == [[1, 2, 3]]
        @test global_cons2 == Int[]
    end
end

@testset "Term Sparsity" begin
    # NOTE: get_term_sparsity_graph is designed for state polynomial optimization
    # (uses expval which is only defined for NCStateWord). Testing with regular
    # Monomials would fail. These tests verify init_activated_supp which works
    # with regular polynomials.

    @testset "Init Activated Support" begin
        registry, (x,) = create_noncommutative_variables([("x", 1:2)])
        T = eltype(indices(registry))
        M = Monomial{NonCommutativeAlgebra,T}
        P = Polynomial{NonCommutativeAlgebra,T,Float64}

        # Get actual variable indices
        x_idx = [x[i].word[1] for i in 1:2]

        # Simple objective and empty constraints
        obj = x[1] + x[2]
        cons = P[]

        # Simple basis with actual variable indices
        one_mono = one(M)
        basis = M[one_mono, Monomial{NonCommutativeAlgebra}(T[x_idx[1]]), Monomial{NonCommutativeAlgebra}(T[x_idx[2]])]

        supp = init_activated_supp(obj, cons, basis)
        @test !isempty(supp)
        # Should contain monomials from objective and diagonal entries
        @test length(supp) >= 2  # At least identity + objective monomials
    end
end
