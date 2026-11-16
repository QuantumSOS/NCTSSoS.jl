using Test,NCTSSoS,Mosek,TOML
include("functional-sparse-models.jl")
const M=SparseFunctionalModels
const N=NCTSSoS.Next
const backend=N.MosekBackend(MSK_IPAR_NUM_THREADS=1)
for trace in (false,true)
    model=M.chsh(trace)
    old=compute_sparsity(model.old,SolverConfig(optimizer=nothing,order=1,cs_algo=NoElimination(),ts_algo=MMD()))
    oldrows=only(old.corr_sparsity.clq_mom_mtx_bases)
    W=N.word_type(model.registry,2)
    rows=[N._pack_word(W,last(only(N.terms(model.convert_word(w))))) for w in oldrows]
    objective=[(ComplexF64(c),N._pack_word(W,w)) for (c,w) in N.terms(model.new.objective)]
    block=N.RowBlock(rows,[(1.0+0im,one(W))],N.HPSD,0)
    expected=sort!(length.(only(only(old.cliques_term_sparsities)).block_bases))
    for order in ("old","new"), seed in ("current","legacy")
        ordered=N.RowBlock(order=="old" ? rows : sort(rows),block.polynomial,N.HPSD,0)
        if seed=="current"
            active=N._initial_support([ordered],objective)
        else
            sink=N.WordSupportSink(Set{W}())
            for (c,w) in objective; N.push_term!(sink,c,w); end
            if trace
                for w in rows; N._block_product!(sink,w,w,block.polynomial); end
            end
            active=sink.support
        end
        graph=N._support_graph(ordered,active)
        family=[N.RowBlock(sort(ordered.rows[c]),block.polynomial,N.HPSD,0) for c in NCTSSoS.clique_decomp(graph,MMD())]
        program,table,_=N.family_program(family,objective;real_form=true)
        solution=N.solve_standard(N.moment_standard_form(program),backend)
        sides=sort!([b.side for b in program.blocks])
        println((;trace,order,seed,rows=length(rows),active=length(active),sides,bound=solution.bound))
        @test solution.bound ≈ -2sqrt(2) atol=1e-6
        order=="old" && seed=="legacy" && @test sides==expected
        flush(stdout)
    end
    expected=TOML.parsefile(joinpath(pkgdir(NCTSSoS),"test/data/expectations",trace ? "chsh_trace.toml" : "chsh_state.toml"))["cases"][2]["expected"]["objective"]
    @test expected ≈ -2sqrt(2) atol=1e-6
end
for (name,model) in (("bilocal",M.bilocal()),("t2-n4",M.t2(4)))
    println(name," degree=",N.degree(model.new.objective)," terms=",length(collect(N.terms(model.new.objective))))
end
