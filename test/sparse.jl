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

using NCTSSoS.FastPolynomials: Arbitrary, get_state_basis, NCStateWord

@testset "Correlative Sparsity without constraints" begin
    @testset "Example 2" begin
        @ncpolyvar x[1:3]
        @ncpolyvar y[1:3]
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function
        pop = polyopt(-f)
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
        @ncpolyvar x[1:n]
        f = 0.
        for i = 1:n
            jset = max(1, i - 5):min(n, i + 1)
            jset = setdiff(jset, i)
            g = sum(x[j] + x[j]^2 for j in jset)
            f += (2 * x[i] + 5 * x[i]^3 + 1 - g)^2
        end

        pop = polyopt(f)

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

@testset "Correlative Sparsity with constrains" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    pop = polyopt(f; ineq_constraints=[g], eq_constraints=[h1])
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
        # compared with G.fadjlist printed in NCTSSOS's `clique_decomp` function
        map2var_sorted(xs, idcs) = sort(map(i -> xs[i], idcs))
        var2vars_dict(xs, idcs) = Dict(zip(xs, map2var_sorted.(Ref(xs), idcs)))
        n = 4
        @ncpolyvar x[1:n]
        f = sum(x[i] * x[mod1(i + 1, n)] for i = 1:n)

        G = get_correlative_graph(sort(x), f, typeof(f)[])

        savegraph("example1.lgz", G)

        @test var2vars_dict(sort(x), G.fadjlist) == Dict(
            x[1] => [x[2], x[4]],
            x[2] => [x[1], x[3]],
            x[3] => [x[2], x[4]],
            x[4] => [x[1], x[3]],
        )

        n = 3
        @ncpolyvar x[1:3]
        f =
            x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3x[2]^2 - 2x[1] * x[2] * x[1] +
            2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
            6.0 * x[3]^2 +
            9x[2]^2 * x[3] +
            9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

        G = get_correlative_graph(sort(x), f, typeof(f)[])
        @test var2vars_dict(sort(x), G.fadjlist) == Dict(
            x[1] => [x[2]],
            x[2] => [x[1], x[3]],
            x[3] => [x[2]],
        )

        savegraph("example2.lgz", G)

        n = 10
        @ncpolyvar x[1:n]
        f = 0.0
        for i = 1:n
            jset = max(1, i - 5):min(n, i + 1)
            jset = setdiff(jset, i)
            f += (2x[i] + 5 * x[i]^3 + 1)^2
            f -= sum([
                4x[i] * x[j] +
                10x[i]^3 * x[j] +
                2x[j] +
                4x[i] * x[j]^2 +
                10x[i]^3 * x[j]^2 +
                2x[j]^2 for j in jset
            ])
            f += sum([
                x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
            ])
        end
        G = get_correlative_graph(sort(x), f, typeof(f)[])

        @test var2vars_dict(sort(x), G.fadjlist) == Dict(zip(x, map2var_sorted.(Ref(x), [[2, 3, 4, 5, 6, 7], [1, 3, 4, 5, 6, 7, 8], [1, 2, 4, 5, 6, 7, 8, 9], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 8, 9, 10], [2, 3, 4, 5, 6, 7, 9, 10], [3, 4, 5, 6, 7, 8, 10], [4, 5, 6, 7, 8, 9]])))

        savegraph("example3.lgz", G)

        n = 3
        @ncpolyvar x[1:3]
        f =
            x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0x[2]^2 - 2x[1] * x[2] * x[1] +
            2x[1] * x[2]^2 * x[1] - x[2] * x[3] - x[3] * x[2] +
            6.0 * x[3]^2 +
            9x[2]^2 * x[3] +
            9x[3] * x[2]^2 - 54x[3] * x[2] * x[3] + 142x[3] * x[2]^2 * x[3]

        cons = vcat([1.0 - x[i]^2 for i = 1:n], [x[i] - 1.0 / 3 for i = 1:n])
        G = get_correlative_graph(sort(x), f, cons)
        @test var2vars_dict(sort(x), G.fadjlist) == Dict(zip(x, map2var_sorted.(Ref(x),
            [[2], [1, 3], [2]]
        )))

        savegraph("example4.lgz", G)

        @ncpolyvar x[1:3] y[1:3]
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function

        G = get_correlative_graph(sort(vcat(x, y)), f, typeof(f)[])
        @test var2vars_dict(sort(vcat(x, y)), G.fadjlist) == Dict(
            x[1] => [y[1], y[2], y[3]],
            x[2] => [y[1], y[2], y[3]],
            x[3] => [y[1], y[2]],
            y[1] => [x[1], x[2], x[3]],
            y[2] => [x[1], x[2], x[3]],
            y[3] => [x[1], x[2]],
        )

        savegraph("example5.lgz", G)

    end

    @testset "Clique Decomposition" begin
        G = loadgraph("example1.lgz")

        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:4)]
        @test sort.(clique_decomp(G, AsIsElimination())) == [[3, 4], [2, 3], [1, 4], [1, 2]]
        @test sort(sort.(clique_decomp(G, MF()))) == sort([[0x0002, 0x0003, 0x0004], [0x0001, 0x0002, 0x0004]])

        rm("example1.lgz")

        G = loadgraph("example2.lgz")

        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:3)]
        @test sort.(clique_decomp(G, AsIsElimination())) == [[0x0002, 0x0003], [0x0001, 0x0002]]
        @test sort.(clique_decomp(G, MF())) == [[0x0002, 0x0003], [0x0001, 0x0002]]


        rm("example2.lgz")

        G = loadgraph("example3.lgz")

        @test sort.(clique_decomp(G, NoElimination())) == [[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]]


        @test sort(sort.(clique_decomp(G, AsIsElimination()))) == [[1, 3, 4, 5, 6, 7, 8], [2, 5, 6, 7, 8, 9, 10], [3, 4, 5, 6, 7, 8, 9], [4, 5, 6, 7, 8, 9, 10]]  # raw results from NCTSSOS needs to be processed to match variable order

        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 3, 4, 5, 6, 7, 8], [2, 5, 6, 7, 8, 9, 10], [3, 4, 5, 6, 7, 8, 9], [4, 5, 6, 7, 8, 9, 10]]

        rm("example3.lgz")

        G = loadgraph("example4.lgz")

        @test sort.(clique_decomp(G, NoElimination())) == [[1, 2, 3]]

        @test sort(sort.(clique_decomp(G, AsIsElimination()))) == [[1, 2], [2, 3]]

        @test sort(sort.(clique_decomp(G, MF()))) == [[1, 2], [2, 3]]
        rm("example4.lgz")

        G = loadgraph("example5.lgz")

        @test sort.(clique_decomp(G, NoElimination())) == [collect(1:6)]
        @test sort.(clique_decomp(G, AsIsElimination())) == [[0x0002, 0x0005], [0x0002, 0x0004], [0x0002, 0x0006], [0x0003, 0x0005], [0x0003, 0x0004], [0x0001, 0x0005], [0x0001, 0x0004], [0x0001, 0x0006]]

        @test sort(sort.(clique_decomp(G, MF()))) == sort([[0x0003, 0x0004, 0x0005], [0x0001, 0x0002, 0x0006], [0x0001, 0x0002, 0x0004, 0x0005]])

        rm("example5.lgz")
    end

    @testset "Assign Constraint" begin
        n = 4
        @ncpolyvar x[1:n]

        cliques = [x[[1, 2, 4]], x[[2, 3, 4]]]
        cons = [1.0 * x[1] * x[2], 1.0 * x[2] * x[3], 1.0 * x[3] * x[4], 1.0 * x[4] * x[1]]

        @test assign_constraint(cliques, cons) == ([[1, 4], [2, 3]], Int[])

        n = 2
        @ncpolyvar x[1:2]
        g = 4.0 - x[1]^2 - x[2]^2
        h1 = x[1] * x[2] + x[2] * x[1] - 2.0
        h2 = -h1
        cons = [g, h1, h2]

        cliques = [x]

        @test assign_constraint(cliques, cons) == ([[1, 2, 3]], Int[])
    end
end

@testset "Term Sparsity" begin
    @testset "Term Sparsity Graph Poly Opt" begin
        # Example 10.2
        @ncpolyvar x y
        activated_support = [
            one(x),
            x^2,
            x * y^2 * x,
            y^2,
            x * y * x * y,
            y * x * y * x,
            x^3 * y,
            y * x^3,
            x * y^3,
            y^3 * x,
        ]

        mtx_basis = [one(x), x, y, x^2, y^2, x * y, y * x]

        sa = SimplifyAlgorithm(comm_gps=[[x, y]], is_unipotent=false, is_projective=false)

        G_tsp = get_term_sparsity_graph([one(x)], activated_support, mtx_basis, sa)
        @test G_tsp.fadjlist == [[4, 5], Int[], Int[], [1, 6], [1, 7], [4, 7], [5, 6]]
        @test sort(term_sparsity_graph_supp(G_tsp, mtx_basis, one(1.0 * x * y), sa)) == sort([
            one(x * y),
            x^2,
            y^2,
            x^4,
            y^4,
            y * x^2 * y,
            x * y^2 * x,
            x^3 * y,
            y^3 * x,
            y * x * y * x,
        ])
        @test sort(term_sparsity_graph_supp(G_tsp, mtx_basis, 1.0 - x^2, sa)) == sort([
            one(x * y),
            x^2,
            y^2,
            y * x * y * x,
            y * x^2 * y,
            y^3 * x,
            y^4,
            x * y^2 * x,
            x^2 * y^2,
            x^3 * y,
            x^4,
            y * x^3 * y * x,
            y * x^4 * y,
            y^2 * x^2 * y * x,
            y^2 * x^2 * y^2,
            x * y * x^2 * y * x,
            x^5 * y,
            x^6,
        ])
    end

    @testset "Test Case 7.2.0" begin
        @ncpolyvar x[1:2] y[1:2]
        sp =
            (-1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) +
             1.0 * ς(x[2] * y[2])) * one(Monomial)

        d = 1

        sa = SimplifyAlgorithm(comm_gps=[x, y], is_unipotent=true, is_projective=false)
        basis = get_state_basis(Arbitrary, [x; y], d, sa)

        init_act_supp = [one(NCStateWord{Arbitrary}), ς(x[1]) * ς(x[1]) * one(Monomial), ς(x[2]) * ς(x[2]) * one(Monomial), ς(y[1]) * ς(y[1]) * one(Monomial), ς(y[2]) * ς(y[2]) * one(Monomial), ς(x[1] * y[1]) * one(Monomial), ς(x[1] * y[2]) * one(Monomial), ς(x[2] * y[1]) * one(Monomial), ς(x[2] * y[2]) * one(Monomial)]


        @testset "Initial Activated Support" begin

            @test init_activated_supp(sp, typeof(sp)[], basis, sa) == init_act_supp
        end

        @testset "Get Term Sparsity Graph" begin
            using NCTSSoS: get_term_sparsity_graph

            G = get_term_sparsity_graph([one(NCStateWord{Arbitrary})], init_act_supp, basis, sa)

            @test G.fadjlist == [[], [6], [7], [8], [9], [2, 8, 9], [3, 8, 9], [4, 6, 7], [5, 6, 7]]
        end

        @testset "Iterate Term Sparse Supp" begin
            using NCTSSoS: iterate_term_sparse_supp
            ts = iterate_term_sparse_supp(init_act_supp, 1.0 * one(NCStateWord{Arbitrary}), basis, MMD(), sa)
        end
    end
end
