# test/correlated_sparsity/core_pipeline_structure.jl

@testset "Core Pipeline Structure" begin
    @testset "Correlative Sparsity Decomposition" begin
        @testset "Example 2 (x[1:3], y[1:3])" begin
            reg, (x, y) = create_noncommutative_variables([("x", 1:3), ("y", 1:3)])
            objective = nc_bipartite_objective(x, y)
            pop = polyopt(-objective, reg)
            order = 3

            @test maximum(length.(correlative_sparsity(pop, order, NoElimination()).cliques)) == 6
            @test maximum(length.(correlative_sparsity(pop, order, MF()).cliques)) == 4
            @test maximum(length.(correlative_sparsity(pop, order, AsIsElimination()).cliques)) == 2
        end

        @testset "Example 1 (n=10 large scale)" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:10)])
            objective = nc_large_scale_objective(x)
            pop = polyopt(objective, reg)
            order = 3

            @test maximum(length.(correlative_sparsity(pop, order, NoElimination()).cliques)) == 10
            @test maximum(length.(correlative_sparsity(pop, order, MF()).cliques)) == 7
            @test maximum(length.(correlative_sparsity(pop, order, AsIsElimination()).cliques)) == 7
        end

        @testset "Constrained n=2" begin
            reg, (x,) = create_noncommutative_variables([("x", 1:2)])
            objective = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
            g = 4.0 - x[1]^2 - x[2]^2
            h1 = x[1] * x[2] + x[2] * x[1] - 2.0
            pop = polyopt(objective, reg; ineq_constraints=[g], eq_constraints=[h1])
            order = 2

            for algo in (NoElimination(), MF(), AsIsElimination())
                sparsity = correlative_sparsity(pop, order, algo)
                @test maximum(length.(sparsity.cliques)) == 2
                @test length.(sparsity.clq_mom_mtx_bases) == [7]
                @test length.(sparsity.clq_localizing_mtx_bases[1]) == [3, 3]
            end
        end
    end
end
