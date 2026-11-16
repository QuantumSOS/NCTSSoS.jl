# Direct tests of retained elimination code, not the retired graph builder.
function elimination_contracts()
    # #492 moves the unchanged helper into Next; support both merge orders.
    decompose=isdefined(Next,:clique_decomp) ? Next.clique_decomp : NCTSSoS.clique_decomp
    @testset "Reviewed elimination graphs" begin
        # Edges are the adjacency fixtures from graph_and_cliques.jl. Keeping
        # the graph input explicit separates decomposition from graph creation.
        cases=(
            ("ring_n4",4,[(1,2),(2,3),(3,4),(1,4)]),
            ("chain_n3",3,[(1,2),(2,3)]),
            ("dense_n10",10,[(i,j) for i in 1:10 for j in i+1:10 if j-i<=6]),
            ("chain_constraints_n3",3,[(1,2),(2,3)]),
            ("bipartite_x3y3",6,[(1,4),(1,5),(1,6),(2,4),(2,5),(2,6),(3,4),(3,5)]))
        for (id,n,edges) in cases
            graph=Graphs.SimpleGraph(n)
            for (a,b) in edges
                Graphs.add_edge!(graph,a,b)
            end
            expected=fixture("correlated_structure","clique_decomp_$id")
            for (key,algorithm) in (("no_elimination",NoElimination()),("mf",MF()),("as_is_elimination",AsIsElimination()))
                @test sort(sort.(decompose(graph,algorithm)))==expected[key]
            end
        end
        disconnected=Graphs.SimpleGraph(4)
        Graphs.add_edge!(disconnected,1,2); Graphs.add_edge!(disconnected,3,4)
        @test sort(sort.(decompose(disconnected,MaximalElimination())))==
            fixture("correlated_structure","clique_decomp_disconnected_maximal")["maximal_elimination"]
        expected=fixture("correlated_structure","wang_magron_example_3_maximal_chordal")
        graph=Graphs.cycle_graph(expected["n_vertices"])
        @test sort(sort.(decompose(graph,MaximalElimination())))==expected["maximal_elimination"]
        cliques=decompose(graph,MF())
        @test maximum(length.(cliques))==expected["mf_max_clique_size"]
        @test length(cliques)==expected["mf_n_cliques"]
    end
    @testset "Natural elimination preserves the input graph" begin
        graph=Graphs.cycle_graph(4)
        @test decompose(graph,NaturalElimination())==[[1,2,4],[2,3,4]]
        @test Graphs.ne(graph)==4
        @test decompose(Graphs.SimpleGraph(3),NaturalElimination())==[[1],[2],[3]]
        @test isempty(decompose(Graphs.SimpleGraph(0),NaturalElimination()))
    end
end
