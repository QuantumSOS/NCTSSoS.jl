@testset "Named letter actions and flat result sizes" begin
    r=Nx.Registry(Nx.SiteAlgebra(Nx.Unipotent(),2),2,
        [:AliceX=>(1,1),:AliceZ=>(1,2),:BobX=>(2,1),:BobZ=>(2,2)])
    vars=Nx.variables(r); a,c,b,d=vars[:AliceX],vars[:AliceZ],vars[:BobX],vars[:BobZ]
    g=Nx.LetterPermutation(a=>d,c=>b,b=>c,d=>a)
    @test sprint(show,g)=="LetterPermutation(AliceX => BobZ, AliceZ => BobX, BobX => AliceZ, BobZ => AliceX)"
    @test occursin(sprint(show,g),sprint(show,Nx.Problem(a*d+c*b;symmetry=(g,))))
    @test Nx.LetterPermutation((-(-a))=>d,c=>b,b=>c,d=>a)==g
    @test_throws r"AliceX.*BobX.*bijective" Nx.LetterPermutation(a=>b)
    @test_throws r"partition" Nx.LetterPermutation(a=>b,b=>a)
    @test_throws r"repeats" Nx.LetterPermutation(a=>c,a=>c)
    for bad in (one(a),-a,2a,a+c,a*c)
        @test_throws r"coefficient-one" Nx.LetterPermutation(bad=>a)
    end
    @test_throws r"operator-letter" Nx.LetterPermutation(1=>2)
    @test_throws r"operator-letter" Nx.LetterPermutation([])
    _,(other,)=Nx.site_variables(Nx.Free(),1,2)
    @test_throws r"one registry" Nx.LetterPermutation(a=>other[1])
    @test_throws r"incompatible algebra" Nx.Problem(a;symmetry=(Nx.LetterPermutation(other[1]=>other[2],other[2]=>other[1]),))
    _,(settings,)=Nx.site_variables(Nx.Unipotent(),1,3)
    @test_throws r"settings outside" Nx.Problem(a;symmetry=(Nx.LetterPermutation(settings[1]=>settings[3],settings[3]=>settings[1]),))
    for factory in (Nx.pauli_variables,Nx.fermionic_variables,Nx.bosonic_variables)
        _,ops=factory(1)
        @test_throws r"SitePermutation.*CliffordMap" Nx.LetterPermutation(ops[1][1]=>ops[1][1])
    end
    for functional in (Nx.RealState(),Nx.RealTrace())
        _,(x,)=Nx.site_variables(Nx.Free(),1,2;functional)
        symbol=functional isa Nx.RealState ? Nx.state_symbol : Nx.trace_symbol
        @test_throws r"operator letters" Nx.LetterPermutation(symbol(x[1])=>x[2])
    end
    # Both routes expose the same flat names; the existing nested field stays.
    example=Nx.solve(Nx.Problem(a*d+c*b;symmetry=(g,)),Nx.Relaxation(order=1);solver=Nx.JuMPBackend(SOLVER))
    @test example.bound ≈ -2 atol=1e-6
    _,(x,_,_)=Nx.pauli_variables(2)
    problem=Nx.Problem(-x[1]-x[2]); backend=Nx.JuMPBackend(SOLVER)
    results=[Nx.solve(problem,relaxation;solver=backend) for relaxation in
        (Nx.Relaxation(order=1),Nx.Relaxation(order=1,grouping=Nx.Correlative()))]
    keys=(:rows,:moments,:real_variables,:blocks,:largest_block,:localizing_rows,:closure)
    for result in results
        @test result.bound ≈ -2 atol=1e-6
        @test all(k->k in propertynames(result.sizes),keys)
        values=NamedTuple{keys}(map(k->getproperty(result.sizes,k),keys))
        @test all(>(0),(values.rows,values.moments,values.real_variables,values.blocks,values.largest_block))
        @test isempty(values.localizing_rows)
    end
    @test results[1].sizes.closure===nothing
    @test results[2].sizes.closure!==nothing
end
