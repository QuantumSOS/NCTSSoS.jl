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

@testset "Signed CHSH symmetry and input guards" begin
    for F in (nothing,Nx.RealState(),Nx.RealTrace())
        _,(p,)=Nx.site_variables(Nx.Projector(),1,1;functional=F)
        @test_throws r"negative Projector.*-P is not a projector" Nx.LetterPermutation(p[1]=>-p[1])
    end
    _,(A,B)=Nx.site_variables(Nx.Unipotent(),2,2)
    g=Nx.LetterPermutation(A[1]=>A[2],A[2]=>A[1],B[2]=>-B[2])
    @test occursin("s2_2 => -s2_2",sprint(show,g))
    @test !isone(g) && Nx.compose(g,g)==one(g)
    for bad in (2A[1],im*A[1],zero(A[1]))
        @test_throws r"target coefficients" Nx.LetterPermutation(A[1]=>bad)
    end
    @test_throws r"bijective" Nx.LetterPermutation(A[1]=>-A[2])
    _,(extra,)=Nx.site_variables(Nx.Unipotent(),1,3)
    @test_throws r"settings outside" Nx.Problem(A[1];symmetry=(Nx.LetterPermutation(extra[3]=>-extra[3]),))
    H=A[1]*B[1]+A[1]*B[2]+A[2]*B[1]-A[2]*B[2]
    _,(otherA,otherB)=Nx.site_variables(Nx.Unipotent(),2,2)
    foreign=Nx.LetterPermutation(otherA[1]=>otherA[2],otherA[2]=>otherA[1],otherB[2]=>-otherB[2])
    @test_throws r"different registry.*this problem's variables" Nx.Problem(H;symmetry=(foreign,))
    @test_throws r"different registry" Nx.Problem(H;symmetry=Nx.Symmetry(;moment_only=(one(foreign),)))
    for same_registry in (one(g),inv(g),Nx.compose(g,g))
        @test Nx.Problem(H;symmetry=(same_registry,)).objective===H
    end
    positional=Nx.LetterPermutation(eltype(H.data.words),Nx.SitePermutation([2,1]))
    @test Nx.Problem(H;symmetry=(positional,)).objective===H
    solver=Nx.JuMPBackend(SOLVER); relaxation=Nx.Relaxation(order=1)
    plain=Nx.solve(Nx.Problem(H),relaxation;solver)
    reduced=Nx.solve(Nx.Problem(H;symmetry=(g,)),relaxation;solver)
    @test plain.bound ≈ -2sqrt(2) atol=1e-6
    @test reduced.bound ≈ plain.bound atol=1e-6
    @test reduced.sizes.real_variables < plain.sizes.real_variables
    @test Nx.expectation(reduced,H) ≈ reduced.bound atol=1e-6
    println("Signed CHSH: bounds=",(plain.bound,reduced.bound)," variables=",(plain.sizes.real_variables,reduced.sizes.real_variables))

    swap=Nx.SitePermutation([2,1])
    letters=Nx.LetterPermutation(A[1]=>B[1],B[1]=>A[1],A[2]=>B[2],B[2]=>A[2])
    for (mixed,explicit) in (((g,swap),(g,letters)),((swap,g),(letters,g)))
        built=Nx.build(Nx.Problem(-H;symmetry=mixed),relaxation)
        reference=Nx.build(Nx.Problem(-H;symmetry=explicit),relaxation)
        f,r=Nx.moment_standard_form(built.program),Nx.moment_standard_form(reference.program)
        @test (f.A,f.b,f.c,f.offset,f.cones)==(r.A,r.b,r.c,r.offset,r.cones)
        @test built.sizes.real_variables==1
        @test [c.side for c in built.program.blocks]==[1,1,1]
        @test Nx.solve(Nx.Problem(-H;symmetry=mixed),relaxation;solver).bound ≈ -2sqrt(2) atol=1e-6
    end
end

@testset "Pairs of public symmetry actions" begin
    relaxation=Nx.Relaxation(order=1)
    site=Nx.SitePermutation([2,1])
    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        _,(a,b)=Nx.site_variables(rule,2,2)
        letter=Nx.LetterPermutation(a[1]=>a[2],a[2]=>a[1])
        swap=Nx.LetterPermutation(a[1]=>b[1],b[1]=>a[1],a[2]=>b[2],b[2]=>a[2])
        h=sum(x^2 for x in (a...,b...))
        for (mixed,explicit) in (((site,letter),(swap,letter)),((letter,site),(letter,swap)))
            f=Nx.moment_standard_form(Nx.build(Nx.Problem(h;symmetry=mixed),relaxation).program)
            r=Nx.moment_standard_form(Nx.build(Nx.Problem(h;symmetry=explicit),relaxation).program)
            @test (f.A,f.b,f.c,f.offset,f.cones)==(r.A,r.b,r.c,r.offset,r.cones)
        end
    end
    _,(x,y,z)=Nx.pauli_variables(2)
    word(p)=last(only(Nx.terms(p)))
    clifford=Nx.CliffordMap(word.([z[1],z[2],x[1],x[2]]))
    anti=Nx.AntiunitaryPauli(one(typeof(word(x[1]))))
    h=-(x[1]*x[2]+z[1]*z[2])
    for pair in ((site,clifford),(site,anti),(clifford,anti)), actions in (pair,reverse(pair))
        result=Nx.solve(Nx.Problem(h;symmetry=actions),relaxation;solver=Nx.JuMPBackend(SOLVER))
        @test result.bound ≈ -2 atol=1e-6
    end
    _,(a,b)=Nx.site_variables(Nx.Unipotent(),2,2)
    letter=Nx.LetterPermutation(a[1]=>a[2],a[2]=>a[1])
    for other in (clifford,anti), pair in ((letter,other),(other,letter)), objective in (h,a[1]+a[2])
        @test_throws r"incompatible algebra.*construct its images with this problem's variables" Nx.build(Nx.Problem(objective;symmetry=pair),relaxation)
    end
    for functional in (Nx.RealState(),Nx.RealTrace())
        _,(a,b)=Nx.site_variables(Nx.Unipotent(),2,2;functional)
        letter=Nx.LetterPermutation(a[1]=>a[2],a[2]=>a[1],b[2]=>-b[2])
        h=-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
        for pair in ((site,letter),(letter,site))
            @test_throws r"State/trace.*remove symmetry.*functional block routing" Nx.build(Nx.Problem(h;symmetry=pair),relaxation)
        end
    end
end
