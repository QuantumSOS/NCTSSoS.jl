struct TestCounterAction <: Nx.AbstractAction
    shift::Int
    order::Int
end
Base.one(g::TestCounterAction) = TestCounterAction(0,g.order)
Nx.compose(g::TestCounterAction,h::TestCounterAction) = TestCounterAction(mod(g.shift+h.shift,g.order),g.order)

@testset "Shared finite-group projections" begin
    @test length(Nx.OrbitGroup([TestCounterAction(1,65535)])) == 65535
    @test_throws r"65,535" Nx.OrbitGroup([TestCounterAction(1,65536)])
    n = 4
    generators = [test_shift(n), test_reflect(n)]
    @test_throws ArgumentError Nx.OrbitGroup(generators; limit=7)
    for model in (:heisenberg, :ising)
        rows = toy_words(n, 1) # X, Y and Z each carry the 2D D4 irrep.
        objective = ring_objective(n, model)
        dense, dt = dense_program(rows, objective)
        p, table, report = Nx.finite_program(rows, objective, generators)
        @test report.group_order == 8
        @test solve_reference(p)[1] ≈ solve_reference(dense)[1] atol=1e-6
        # Check matrix entries including cross-orbit couplings, not just a
        # symmetric objective that might hide a projection mismatch.
        y = randn(MersenneTwister(15), p.nvars); y[1] = 1
        yd = map(dt.labels) do l
            id, phase = Nx.orbit_moment!(table, l.word)
            id == 0 ? 0. : Nx.scale(phase, y[id])
        end
        D = only(evaluated_blocks(dense, yd))
        @test eigmin(Hermitian(D)) ≈ minimum(eigmin(Hermitian(A)) for A in evaluated_blocks(p,y)) atol=1e-10
        index = Dict(w => i for (i,w) in enumerate(rows))
        for (b, A) in zip(p.blocks, evaluated_blocks(p,y))
            sector = report.sectors[b.label]
            V = zeros(ComplexF64, length(rows), length(sector))
            for (j, row) in enumerate(sector), (c,w) in row
                V[index[w], j] = c
            end
            @test V'D*V ≈ A atol=1e-11
        end
    end
end
