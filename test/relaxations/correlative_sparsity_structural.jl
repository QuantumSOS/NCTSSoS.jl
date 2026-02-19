# test/relaxations/correlative_sparsity_structural.jl
# Structural tests for correlative sparsity (solver-free).
#
# Included from the Minimal Suite so CI default always exercises the core
# graph/clique/assignment logic without solving any SDPs.

using Test, NCTSSoS
using Graphs, CliqueTrees

using NCTSSoS: assign_constraint, clique_decomp, get_correlative_graph

function _adjacency_by_index(G::SimpleGraph, sorted_indices)
    Dict(
        sorted_indices[i] => sort([sorted_indices[j] for j in neighbors(G, i)])
        for i in 1:nv(G)
    )
end

@testset "Correlative Sparsity (Structural, Solver-Free)" begin
    @testset "get_correlative_graph: ring (n=4)" begin
        n = 4
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])
        x_idx = [x[i].word[1] for i in 1:n]

        f = sum(x[i] * x[mod1(i + 1, n)] for i in 1:n)
        G, sorted_indices, _ = get_correlative_graph(registry, f, typeof(f)[])
        adj = _adjacency_by_index(G, sorted_indices)

        @test sort(adj[x_idx[1]]) == sort([x_idx[2], x_idx[4]])
        @test sort(adj[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
        @test sort(adj[x_idx[3]]) == sort([x_idx[2], x_idx[4]])
        @test sort(adj[x_idx[4]]) == sort([x_idx[1], x_idx[3]])
    end

    @testset "get_correlative_graph: chain with constraints (n=3)" begin
        n = 3
        registry, (x,) = create_noncommutative_variables([("x", 1:n)])
        x_idx = [x[i].word[1] for i in 1:n]

        f =
            x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 - 2.0 * x[1] * x[2] * x[1] +
            2.0 * x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2

        cons = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
        G, sorted_indices, _ = get_correlative_graph(registry, f, cons)
        adj = _adjacency_by_index(G, sorted_indices)

        @test adj[x_idx[1]] == [x_idx[2]]
        @test sort(adj[x_idx[2]]) == sort([x_idx[1], x_idx[3]])
        @test adj[x_idx[3]] == [x_idx[2]]
    end

    @testset "clique_decomp: elimination variants (ring n=4)" begin
        G = SimpleGraph(4)
        add_edge!(G, 1, 2)
        add_edge!(G, 2, 3)
        add_edge!(G, 3, 4)
        add_edge!(G, 4, 1)

        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:4)]
        @test sort(sort.(clique_decomp(G, AsIsElimination()))) ==
              [[1, 2], [1, 4], [2, 3], [3, 4]]
        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 2, 4], [2, 3, 4]]
    end

    @testset "assign_constraint: clique vs global" begin
        registry, (x,) = create_noncommutative_variables([("x", 1:4)])
        T = eltype(indices(registry))
        x_idx = [x[i].word[1] for i in 1:4]

        cliques = [T[x_idx[1], x_idx[2]], T[x_idx[3], x_idx[4]]]
        cons = [1.0 * x[2] * x[3], 1.0 * x[1] * x[2]]

        clq_cons, global_cons = assign_constraint(cliques, cons, registry)

        @test clq_cons == [[2], Int[]]
        @test global_cons == [1]
    end
end

