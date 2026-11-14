@testset "Momentum little groups" begin
    for (n,degree) in ((3,1), (4,1), (4,2)), model in (:heisenberg, :ising)
        rows = toy_words(n, degree)
        objective = ring_objective(n, model)
        point = if model == :heisenberg
            [test_reflect(n), test_map(n; axes=(1,2,0)), test_map(n; axes=(1,0,2), signs=(1,-1,1))]
        else
            [test_reflect(n), test_map(n; signs=(1,-1,-1))]
        end
        p, table, report = Nx.lattice_program(rows, objective, test_shift(n), n, point)
        dense, dt = dense_program(rows, objective)
        @test report.point_order == (model == :heisenberg ? 48 : 4)
        @test solve_reference(p)[1] ≈ solve_reference(dense)[1] atol=1e-6
        relations = [Nx.orbit_moment!(table,l.word) for l in dt.labels]
        y = randn(MersenneTwister(617), length(table.labels)); y[1] = 1
        yd = [id == 0 ? 0. : Nx.scale(phase,y[id]) for (id,phase) in relations]
        @test eigmin(Hermitian(only(evaluated_blocks(dense,yd)))) ≈
            minimum(eigmin(Hermitian(A)) for A in evaluated_blocks(p,y)) atol=1e-9
        panti, _, _ = Nx.lattice_program(rows, objective, test_shift(n), n, point;
            moment_generators=(test_conjugate(n),))
        @test all(b.kind == Nx.PSD for b in panti.blocks)
        @test [b.side for b in panti.blocks] == [b.side for b in p.blocks]
        @test solve_reference(panti)[1] ≈ solve_reference(dense)[1] atol=1e-6
    end
end
