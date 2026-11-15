# The same I_nn22 coefficients as perf/next/probe.jl, including its I3322 case.
function letter_bell(a,b)
    n=length(a)
    c=if n==3
        [0 -1 -2 0; 0 1 1 -1; -1 1 1 1; 0 -1 1 0]
    else
        matrix=zeros(Int,n+1,n+1)
        for i in 1:n-1
            matrix[1,i+1]=matrix[i+1,1]=-1
        end
        for i in 1:n,j in 1:n
            matrix[i+1,j+1]=(i==n || j==n) ? (i==j ? -1 : 1) :
                i+j<n ? 0 : i+j==n ? 1 : -1
        end
        matrix
    end
    sum(c[i+1,1]*a[i]+c[1,i+1]*b[i] for i in 1:n)+
        sum(c[i+1,j+1]*a[i]*b[j] for i in 1:n,j in 1:n)
end

@testset "Correlative site-word groups are letter cliques" begin
    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        registry,(a,b)=Nx.site_variables(rule,2,3)
        # Different letters on the same physical site need not share a group.
        problem=Nx.Problem(a[1]*b[1]+a[2]*b[2];constraints=(Nx.Positive(2+a[1]+a[3]),))
        built=Nx.build(problem,Nx.Relaxation(order=2,grouping=Nx.Correlative()))
        groups=[filter(w->Nx.degree(w)==1,block.rows) for block in built.sizes.closure.family if block.source==0]
        letter_key(group)=Tuple(sort!([only(Nx.legacy_letters(w)) for w in group]))
        expected=((a[1],b[1]),(a[2],b[2]),(a[1],a[3]),(b[3],))
        @test Set(letter_key.(groups))==Set(letter_key([only(p.data.words) for p in group]) for group in expected)
        for block in built.sizes.closure.family
            atoms=Set(v for w in block.rows for v in Nx.word_letters(w))
            @test any(group->issubset(atoms,group),groups)
            if block.source==1
                @test any(group->issubset(atoms,group) && all(cw->all(in(group),Nx.word_letters(last(cw))),block.polynomial),groups)
            end
        end
        contiguous=Nx.build(problem,Nx.Relaxation(order=2,grouping=Nx.Correlative(),row_rule=Nx.Contiguous()))
        @test contiguous.program.forms.coefs==built.program.forms.coefs
        @test contiguous.program.forms.ids==built.program.forms.ids
        # PhysicsInformed still describes physical-site regions.
        physical=Nx.build(problem,Nx.Relaxation(order=1,grouping=Nx.PhysicsInformed(level=1)))
        @test [sort!(unique!(reduce(vcat,collect.(Nx.sites.(block.rows));init=Int[])))
            for block in physical.sizes.closure.family if block.source==0]==[[1,2],[1]]
        callback=(W,r,group,order)->Nx.basis_rows(W,r.algebra,group,order)
        selected=Nx.build(problem,Nx.Relaxation(order=2,grouping=Nx.Correlative(),row_rule=callback))
        @test selected.program.forms.ids==built.program.forms.ids
        @test selected.program.forms.coefs==built.program.forms.coefs
        outside=(W,r,group,order)->Nx.basis_rows(W,r.algebra,1:r.nsites,order)
        @test_throws r"on its group" Nx.build(problem,Nx.Relaxation(order=2,grouping=Nx.Correlative(),row_rule=outside))
    end
    n=5
    registry,(a,b)=Nx.site_variables(Nx.Projector(),2,n)
    old_registry,(old_a,old_b)=NCTSSoS.create_projector_variables([("A",1:n),("B",1:n)])
    old_pop=NCTSSoS.polyopt(-letter_bell(old_a,old_b),old_registry)
    old=NCTSSoS.correlative_sparsity(old_pop,5,NCTSSoS.MMD())
    W=Nx.word_type(registry,10)
    cliques=Nx.build(Nx.Problem(-letter_bell(a,b)),Nx.Relaxation(order=1,grouping=Nx.Correlative(algorithm=NCTSSoS.MMD())))
    groups=[[Nx.word_from_legacy(W,Nx.legacy_letters(w)) for w in block.rows if !isone(w)]
        for block in cliques.sizes.closure.family]
    # Registry letters need not use the UInt32 encoding of legacy_letters.
    old_ids=Dict(Nx._pack_word(W,only(p.data.words))=>Int(only(q.word))
        for (p,q) in zip([a;b],[old_a;old_b]))
    key(group)=sort!([old_ids[w] for w in group])
    @test sort(key.(groups))==sort([Int.(clique) for clique in old.cliques])
    # Test row generation at the acceptance order without constructing its SDP.
    bases=[Nx.basis_rows(W,registry.algebra,group,5) for group in groups]
    by_clique=Dict(Tuple(Int.(clique))=>basis for (clique,basis) in zip(old.cliques,old.clq_mom_mtx_bases))
    for (group,basis) in zip(groups,bases)
        @test Set(Tuple(old_ids[v] for v in Nx.word_letters(w)) for w in basis)==
            Set(Tuple(Int.(w.word)) for w in by_clique[Tuple(key(group))])
    end
end

@testset "I3322 correlative and term-sparse bound equals current package" begin
    registry,(a,b)=Nx.site_variables(Nx.Projector(),2,3)
    problem=Nx.Problem(-letter_bell(a,b))
    description=Nx.Relaxation(order=3,grouping=Nx.Correlative(algorithm=NCTSSoS.MMD()),split=Nx.TermSplit())
    result=Nx.solve(problem,description;solver=Nx.JuMPBackend(SOLVER))
    old_registry,(a,b)=NCTSSoS.create_projector_variables([("A",1:3),("B",1:3)])
    old=NCTSSoS.cs_nctssos(NCTSSoS.polyopt(-letter_bell(a,b),old_registry),
        NCTSSoS.SolverConfig(optimizer=SOLVER,order=3,cs_algo=NCTSSoS.MMD(),ts_algo=NCTSSoS.MaximalElimination());dualize=false)
    @test result.status==MOI.OPTIMAL
    @test result.bound≈old.objective atol=1e-6
end
