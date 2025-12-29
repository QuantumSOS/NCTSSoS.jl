# Sparsity Tests
# ===============
# Includes:
#   - Correlative sparsity component tests (formerly sparse.jl)
#   - Sparsity verification tests (dense vs sparse comparison)

using Test, NCTSSoS
using Graphs, CliqueTrees

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
# Correlative Sparsity Component Tests (merged from sparse.jl)
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
# Sparsity Verification Tests
# =============================================================================

@testset "Bell Inequalities with Sparsity" begin
    @testset "CHSH - Dense vs Sparse" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        
        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(-f, reg)
        
        # Dense (no sparsity)
        config_dense = SolverConfig(optimizer=SOLVER, order=1, 
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # Correlative sparsity only (MF)
        config_cs = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_cs = cs_nctssos(pop, config_cs)
        
        # Term sparsity only (MMD)
        config_ts = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=NoElimination(), ts_algo=MMD())
        result_ts = cs_nctssos(pop, config_ts)
        
        expected = -2.8284
        
        @test result_dense.objective ≈ expected atol=1e-4
        @test result_cs.objective ≈ expected atol=1e-4
        @test result_ts.objective ≈ expected atol=1e-4
        
        @test result_dense.objective ≈ result_cs.objective atol=1e-5
        @test result_dense.objective ≈ result_ts.objective atol=1e-5
    end
    
    @testset "I_3322 - Dense vs Sparse" begin
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + 
            x[2] * (y[1] + y[2] - y[3]) + 
            x[3] * (y[1] - y[2]) - 
            x[1] - 2.0 * y[1] - y[2]
        
        pop = polyopt(-f, reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # With MF correlative sparsity
        config_mf = SolverConfig(optimizer=SOLVER, order=2,
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_mf = cs_nctssos(pop, config_mf)
        
        expected = -0.25
        
        @test result_dense.objective ≈ expected atol=1e-2
        @test result_mf.objective ≈ expected atol=1e-2
        @test result_mf.objective <= result_dense.objective + 1e-2
    end
end

@testset "Trace Polynomials with Sparsity" begin
    @testset "CHSH Trace - Dense vs Sparse" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]
        
        p = -1.0 * tr(x[1] * y[1]) - tr(x[1] * y[2]) - tr(x[2] * y[1]) + tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=1,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(tpop, config_dense)
        
        # With MaximalElimination term sparsity
        config_ts = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=NoElimination(), ts_algo=MaximalElimination())
        result_ts = cs_nctssos(tpop, config_ts)
        
        expected = -2.8284
        
        @test result_dense.objective ≈ expected atol=1e-4
        @test result_ts.objective ≈ expected atol=1e-4
        @test result_dense.objective ≈ result_ts.objective atol=1e-5
    end
    
    @testset "Covariance Trace - Dense vs Sparse" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:6)])
        x = vars[1:3]
        y = vars[4:6]
        
        cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
        p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
        tpop = polyopt(p * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(tpop, config_dense)
        
        # With MF + MMD
        config_sparse = SolverConfig(optimizer=SOLVER, order=2,
                                     cs_algo=MF(), ts_algo=MMD())
        result_sparse = cs_nctssos(tpop, config_sparse)
        
        expected = -5.0
        
        @test result_dense.objective ≈ expected atol=1e-4
        @test result_sparse.objective ≈ expected atol=1e-4
        @test result_dense.objective ≈ result_sparse.objective atol=1e-4
    end
end

@testset "State Polynomials with Sparsity" begin
    @testset "CHSH State - Dense vs Sparse" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        
        sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
        spop = polyopt(sp * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=1,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(spop, config_dense)
        
        # With MMD term sparsity
        config_ts = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=NoElimination(), ts_algo=MMD())
        result_ts = cs_nctssos(spop, config_ts)
        
        expected = -2.8284
        
        @test result_dense.objective ≈ expected atol=1e-4
        @test result_ts.objective ≈ expected atol=1e-4
        @test result_dense.objective ≈ result_ts.objective atol=1e-5
    end
    
    @testset "Covariance State - Dense vs Sparse" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
        
        cov(a, b) = 1.0 * ς(x[a] * y[b]) - ς(x[a]) * ς(y[b])
        sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
        spop = polyopt(sp * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(spop, config_dense)
        
        # With MF + MMD
        config_sparse = SolverConfig(optimizer=SOLVER, order=2,
                                     cs_algo=MF(), ts_algo=MMD())
        result_sparse = cs_nctssos(spop, config_sparse)
        
        expected = -5.0
        
        @test result_dense.objective ≈ expected atol=1e-2
        @test result_sparse.objective ≈ expected atol=1e-2
        @test result_dense.objective ≈ result_sparse.objective atol=1e-2
    end
end

@testset "Constrained POP with Sparsity" begin
    @testset "Ball Constraint - Dense vs Sparse" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        
        f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2 + x[1]*x[2]*x[1]*x[2] + x[2]*x[1]*x[2]*x[1] +
            x[1]^3*x[2] + x[2]*x[1]^3 + x[1]*x[2]^3 + x[2]^3*x[1]
        
        g1 = 1.0 - x[1]^2
        g2 = 1.0 - x[2]^2
        
        pop = polyopt(f, reg; ineq_constraints=[g1, g2])
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # With MMD term sparsity
        config_ts = SolverConfig(optimizer=SOLVER, order=2,
                                 cs_algo=NoElimination(), ts_algo=MMD())
        result_ts = cs_nctssos(pop, config_ts)
        
        # Both should produce negative lower bounds for this minimization problem
        @test result_dense.objective < 0
        @test result_ts.objective < 0
        @test result_dense.objective < 5.0
        @test result_ts.objective < 5.0
    end
    
    @testset "Rosenbrock - Dense vs Sparse" begin
        n = 6
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])
        
        f = Float64(n) * one(typeof(x[1]))
        for i = 2:n
            f = f + 100.0 * x[i-1]^4 - 200.0 * x[i-1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
        end
        
        pop = polyopt(f, reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # With MF correlative sparsity
        config_cs = SolverConfig(optimizer=SOLVER, order=2,
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_cs = cs_nctssos(pop, config_cs)
        
        # With both sparsities
        config_both = SolverConfig(optimizer=SOLVER, order=2,
                                   cs_algo=MF(), ts_algo=MMD())
        result_both = cs_nctssos(pop, config_both)
        
        # All should give similar results (lower bounds)
        @test result_cs.objective <= result_dense.objective + 1e-2
        @test result_both.objective <= result_dense.objective + 1e-2
    end
end

@testset "Sparsity Algorithm Variants" begin
    @testset "Correlative Sparsity Algorithms" begin
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + 
            x[2] * (y[1] + y[2] - y[3]) + 
            x[3] * (y[1] - y[2]) - 
            x[1] - 2.0 * y[1] - y[2]
        
        pop = polyopt(-f, reg)
        
        expected = -0.25
        
        # NoElimination (dense) - should give correct result
        config_dense = SolverConfig(optimizer=SOLVER, order=2, 
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        @test result_dense.objective ≈ expected atol=1e-2
        
        # MF - should give correct result (good clique tree)
        config_mf = SolverConfig(optimizer=SOLVER, order=2, 
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_mf = cs_nctssos(pop, config_mf)
        @test result_mf.objective ≈ expected atol=1e-2
        
        # AsIsElimination can give looser bounds
        config_asis = SolverConfig(optimizer=SOLVER, order=2, 
                                   cs_algo=AsIsElimination(), ts_algo=NoElimination())
        result_asis = cs_nctssos(pop, config_asis)
        @test result_asis.objective < 0
    end
    
    @testset "Term Sparsity Algorithms" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]
        
        p = -1.0 * tr(x[1] * y[1]) - tr(x[1] * y[2]) - tr(x[2] * y[1]) + tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)
        
        expected = -2.8284
        
        for (name, algo) in [
            ("NoElimination", NoElimination()),
            ("MMD", MMD()),
            ("MaximalElimination", MaximalElimination())
        ]
            config = SolverConfig(optimizer=SOLVER, order=1, cs_algo=NoElimination(), ts_algo=algo)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ expected atol=1e-4
        end
    end
end
