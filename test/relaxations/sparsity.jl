# =============================================================================
# Sparsity Component Tests
# =============================================================================
# Tests for sparsity algorithm components:
#   - Correlative graph construction
#   - Clique decomposition
#   - Constraint assignment
#   - Term sparsity initialization
#   - PolyOptResult fields
#
# NOTE: These tests are commented out as they take too long.
#       Sparsity Algorithm Variants tests moved to problems/trace_polynomial/
# =============================================================================

using Test, NCTSSoS
using Graphs, CliqueTrees

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "standalone_setup.jl"))

using NCTSSoS:
    assign_constraint,
    get_correlative_graph,
    clique_decomp,
    get_term_sparsity_graph,
    term_sparsity_graph_supp,
    correlative_sparsity,
    init_activated_supp,
    variable_indices,
    neat_dot

# =============================================================================
# Correlative Sparsity Component Tests
# =============================================================================

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

            expected_neighbors = [[2, 3, 4, 5, 6, 7], [1, 3, 4, 5, 6, 7, 8], [1, 2, 4, 5, 6, 7, 8, 9], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 8, 9, 10], [2, 3, 4, 5, 6, 7, 9, 10], [3, 4, 5, 6, 7, 8, 10], [4, 5, 6, 7, 8, 9]]

            adj = graph_adjacency_by_index(G, sorted_indices)
            for i in 1:n
                @test sort(adj[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
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

            expected_neighbors = [[2], [1, 3], [2]]

            adj = graph_adjacency_by_index(G, sorted_indices)

            for i in 1:n
                @test sort(adj[x_idx[i]]) == sort(x_idx[expected_neighbors[i]])
            end
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

            adj = graph_adjacency_by_index(G, sorted_indices)

            # Check connectivity pattern: x[i] connects to y[j], y[j] connects to x[i]
            # x[1] -> y[1], y[2], y[3]
            @test sort(adj[x_idx[1]]) == sort([y_idx[1], y_idx[2], y_idx[3]])
            # x[2] -> y[1], y[2], y[3]
            @test sort(adj[x_idx[2]]) == sort([y_idx[1], y_idx[2], y_idx[3]])
            # x[3] -> y[1], y[2]
            @test sort(adj[x_idx[3]]) == sort([y_idx[1], y_idx[2]])
            # y[1] -> x[1], x[2], x[3]
            @test sort(adj[y_idx[1]]) == sort([x_idx[1], x_idx[2], x_idx[3]])
            # y[2] -> x[1], x[2], x[3]
            @test sort(adj[y_idx[2]]) == sort([x_idx[1], x_idx[2], x_idx[3]])
            # y[3] -> x[1], x[2]
            @test sort(adj[y_idx[3]]) == sort([x_idx[1], x_idx[2]])
        end
    end

    @testset "Clique Decomposition" begin
        G = loadgraph("example1.lgz")
        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:4)]
        @test sort(sort.(clique_decomp(G, AsIsElimination()))) == [[1, 2], [1, 4], [2, 3], [3, 4]]
        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 2, 4], [2, 3, 4]]
        rm("example1.lgz")

        G = loadgraph("example2.lgz")
        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:3)]
        @test sort(sort.(clique_decomp(G, AsIsElimination()))) == [[1, 2], [2, 3]]
        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 2], [2, 3]]
        rm("example2.lgz")

        G = loadgraph("example3.lgz")
        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:10)]
        @test sort(sort.(clique_decomp(G, AsIsElimination()))) ==  [[1, 2, 3, 4, 5, 6, 7], [2, 3, 4, 5, 6, 7, 8], [3, 4, 5, 6, 7, 8, 9], [4, 5, 6, 7, 8, 9, 10]]
        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 2, 3, 4, 5, 6, 7], [2, 3, 4, 5, 6, 7, 8], [3, 4, 5, 6, 7, 8, 9], [4, 5, 6, 7, 8, 9, 10]]
        rm("example3.lgz")

        G = loadgraph("example4.lgz")
        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:3)]
        @test sort(sort.(clique_decomp(G, AsIsElimination()))) == [[1, 2], [2, 3]]
        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 2], [2, 3]]
        rm("example4.lgz")

        G = loadgraph("example5.lgz")
        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:6)]
        @test sort(sort.(clique_decomp(G, AsIsElimination()))) == [[1, 4], [1, 5], [1, 6], [2, 4], [2, 5], [2, 6], [3, 4], [3, 5]]
        @test sort(sort.(clique_decomp(G, MF()))) == [
            [1, 2, 4, 5], [1, 2, 6], [3, 4, 5]
        ]
        rm("example5.lgz")
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

# =============================================================================
# PolyOptResult Fields Tests
# =============================================================================

@testset "PolyOptResult Fields" begin
    @testset "moment_matrix_sizes and n_unique_moment_matrix_elements" begin
        # Tiny unipotent problem with fully connected moment graph
        # Use objective that includes both u[1] and u[2] to ensure full basis
        reg, (u,) = create_unipotent_variables([("u", 1:2)])

        # Objective includes both u[1] and u[2] to ensure full basis is used
        f = 1.0 * u[1] * u[2]
        pop = polyopt(f, reg)

        # Use order=1 with no sparsity exploitation for predictable structure
        config = SolverConfig(optimizer=SOLVER, order=1,
                             cs_algo=NoElimination(), ts_algo=NoElimination())
        result = cs_nctssos(pop, config)

        # With 2 variables at order=1:
        # Basis is {1, u[1], u[2]} (3 elements)
        # Moment matrix is 3x3, so sizes should be [[3]]
        @test result.moment_matrix_sizes == [[3]]

        # Unique moment elements: from basis[i]† * basis[j]
        # After symmetric_canon, we get unique monomials:
        # 1*1=1, 1*u[1]=u[1], 1*u[2]=u[2], u[1]*1=u[1], u[1]*u[1]=1 (u[1]^2=1),
        # u[1]*u[2]=u[1]u[2], u[2]*1=u[2], u[2]*u[1]=u[2]u[1]->u[1]u[2] (symmetric), u[2]*u[2]=1
        # Unique after canon: 1, u[1], u[2], u[1]u[2]
        @test result.n_unique_moment_matrix_elements == 4
    end

    @testset "minimal single variable" begin
        # Even simpler: single variable
        reg, (u,) = create_unipotent_variables([("u", 1:1)])
        f = 1.0 * u[1]
        pop = polyopt(f, reg)

        config = SolverConfig(optimizer=SOLVER, order=1,
                             cs_algo=NoElimination(), ts_algo=NoElimination())
        result = cs_nctssos(pop, config)

        # With 1 variable at order=1: basis = {1, u[1]}, size = 2
        @test result.moment_matrix_sizes == [[2]]

        # Moment entries: 1*1=1, 1*u[1]=u[1], u[1]*1=u[1], u[1]*u[1]=1
        # Unique after canon: {1, u[1]}
        @test result.n_unique_moment_matrix_elements == 2
    end
end
