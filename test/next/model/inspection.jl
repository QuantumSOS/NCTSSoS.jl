module InspectionChecks
using Test
using NCTSSoS
import NCTSSoS as N

word(p)=last(only(N.terms(p)))
names(r,w)=Tuple(N.operator_name(r,l) for l in N.word_letters(w))
function consume(p)
    result=0.0
    for (c,w) in N.terms(p)
        result+=real(c)
        for l in N.word_letters(w), site in N.sites(l)
            result+=site+N.degree(l)
        end
    end
    return result
end
bytes(p)=(consume(p); @allocated consume(p))

@testset "Public word and polynomial inspection" begin
    r,(a,ad)=fermionic_variables(2)
    p=ad[2]*a[1]+2a[2]*ad[2]
    @test length(N.terms(p))==3
    @test Dict(names(r,w)=>c for (c,w) in N.terms(p))==
        Dict(()=>2,(:ad2,:a2)=>-2,(:ad2,:a1)=>1)
    @test N.degree(p)==2
    @test all(iseven(N.degree(w)) for (_,w) in N.terms(p))
    @test isodd(N.degree(word(a[1])))
    @test collect(N.sites(word(ad[2]*a[2])))==[2]
    @test N.variables(r)[:ad2]==ad[2]
    @test sprint(show,word(ad[2]*a[1]);context=:registry=>r)=="ad2 a1"
    @test occursin("operator(site=2",sprint(show,word(ad[2]*a[1])))
    @test bytes(p)==0
    @test_throws ArgumentError N.operator_name(r,word(one(a[1])))
    @test_throws ArgumentError N.operator_name(r,word(ad[2]*a[1]))
    rb,(b,bd)=bosonic_variables(1)
    @test Dict(names(rb,w)=>c for (c,w) in N.terms(b[1]*bd[1]))==Dict(()=>1,(:bd1,:b1)=>1)
    @test_throws ArgumentError N.operator_name(r,word(b[1]))
    for rule in (Free(),Projector(),Unipotent())
        rs,(x,y)=site_variables(rule,2,2)
        w=word(x[2]*x[1]*y[1])
        @test names(rs,w)==(:s1_2,:s1_1,:s2_1)
        @test sprint(show,w;context=:registry=>rs)=="s1_2 s1_1 s2_1"
        @test bytes(x[2]*x[1]*y[1])==0
    end
    for (functional,symbol,mark) in ((N.RealState(),state_symbol,"⟨"),(N.RealTrace(),trace_symbol,"tr("))
        rf,(x,)=site_variables(Free(),1,2;functional)
        w=word(symbol(x[1]*x[2])*x[1])
        @test N.degree(w)==3
        @test_throws ArgumentError N.word_letters(w)
        @test sort!([N.operator_name(rf,l) for l in N.support_letters(w)])==[:s1_1,:s1_1,:s1_2]
        @test occursin(mark,sprint(show,w;context=:registry=>rf))
        @test occursin(mark,sprint(show,w))
        @test_throws ArgumentError N.operator_name(rf,word(symbol(x[1])))
        rows=collect(N.basis(rf,1))
        @test length(rows)==5 # identity, two operators, two scalar symbols
        @test all(first(only(N.terms(p)))==1 for p in rows)
        @test_throws ArgumentError N.basis(rf,1000)
    end
    rp,(x,y,z)=pauli_variables(1)
    sparse=N.Registry(Pauli(),1,[:spin=>(1,0)])
    @test N.operator_name(sparse,word(x[1]))==:spin
    @test_throws ArgumentError N.operator_name(sparse,word(z[1]))
    @test bytes(x[1]+y[1]+z[1])==0
    @test_throws ArgumentError N.basis(rb,256)
    rs,_=site_variables(Free(),1,2)
    @test_throws ArgumentError N.basis(rs,1000)
end

@testset "Public basis and closed block rows" begin
    for (r,ops,count) in ((pauli_variables(2)...,16),(fermionic_variables(2)...,11),
            (bosonic_variables(1)...,6),(site_variables(Free(),1,2)...,7),
            (site_variables(Unipotent(),1,2)...,5),(site_variables(Projector(),1,2)...,5))
        rows=collect(N.basis(r,2))
        @test length(rows)==count
        @test all(length(N.terms(p))==1 && first(only(N.terms(p)))==1 for p in rows)
        @test all(N.degree(p)<=2 for p in rows)
        @test isone(word(first(rows)))
        objective=first(first(ops)); objective=objective+objective'
        plain=N.build(Problem(objective),Relaxation(order=2))
        custom=N.build(Problem(objective),Relaxation(rows=rows))
        # Compare complete standard forms, independently of display or packed width.
        a,b=N.moment_standard_form(plain.program),N.moment_standard_form(custom.program)
        @test all(getproperty(a,f)==getproperty(b,f) for f in (:A,:b,:c,:offset,:cones))
        @test [names(r,w) for w in only(N.block_rows(plain))]==[names(r,word(p)) for p in rows]
        @test length(N.basis(r,0))==1
        @test_throws ArgumentError N.basis(r,-1)
    end
    r,(x,)=site_variables(Free(),1,2)
    problem=Problem(x[1];constraints=(Positive(1-x[1]^2),Zero(x[2]),
        ExpectationZero(x[1]),ExpectationPositive(1+x[1])))
    for relaxation in (Relaxation(order=2),Relaxation(rows=collect(N.basis(r,2))))
        built=N.build(problem,relaxation)
        @test length.(N.block_rows(built))==[7,3,3]
        @test N.block_sides(built)==[7,3] # no ZERO or scalar NONNEGATIVE blocks
    end
    for (d,nrows) in ((1,7),(2,7),(3,3),(4,3),(5,1))
        built=N.build(Problem(x[1];constraints=(Positive(1-x[1]^d),)),Relaxation(order=3))
        @test length.(N.block_rows(built))==N.block_sides(built)==[15,nrows]
    end
    r,(x,y,z)=pauli_variables(2)
    h=x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    @test N.block_sides(N.build(Problem(h),Relaxation(order=1)))==[7]
    swap=N.build(Problem(h;symmetry=(SitePermutation([2,1]),)),Relaxation(order=1))
    @test sort(N.block_sides(swap))==[3,4]
    @test length(N.block_rows(swap))==1
    problem=Problem(-z[1]-z[2];symmetry=(SitePermutation([2,1]),))
    built=N.build(problem,Relaxation(order=1,grouping=Correlative()))
    rows=collect(N.block_rows(built))
    @test length(rows)==2
    @test Set(Tuple(names(r,w) for w in block) for block in rows)==
        Set([((),(:Z1,),(:X1,),(:Y1,)),((),(:Z2,),(:X2,),(:Y2,))])
    @test all(rows[i]===built.sizes.closure.family[i].rows for i in eachindex(rows))
    r,(x,y,z)=pauli_variables(4)
    problem=Problem(sum(z[i]*z[mod1(i+1,4)] for i in 1:4);symmetry=(SitePermutation([2,3,4,1]),))
    built=N.build(problem,Relaxation(order=1,grouping=Correlative()))
    @test built.sizes.closure.added_blocks>0
    @test all(rows===block.rows for (rows,block) in zip(N.block_rows(built),built.sizes.closure.family))
end

@testset "Short user displays" begin
    r,ops=pauli_variables(4)
    @test sprint(show,(r,ops);context=:limit=>true)==
        "(Registry(Pauli, 4 sites, 12 variables), ([X1, X2, X3, X4], [Y1, Y2, Y3, Y4], [Z1, Z2, Z3, Z4]))"
    @test sprint(show,MIME"text/plain"(),ops[1])=="[X1, X2, X3, X4]"
    @test occursin("Names: X1, Y1, Z1",sprint(show,MIME"text/plain"(),r))
    r,(x,y,z)=pauli_variables(20)
    @test occursin("… (20 entries)",sprint(show,x))
    @test occursin("… (20 terms)",sprint(show,sum(x)))
    @test !occursin("…",sprint(show,sum(x);context=:limit=>false))
    @test occursin("X20",sprint(show,x;context=:limit=>false))
    @test occursin("…",sprint(show,MIME"text/plain"(),r))
    @test occursin("Z20",sprint(show,MIME"text/plain"(),r;context=:limit=>false))
    @test sprint(show,zero(x[1]))=="0"
    @test N.degree(zero(x[1]))==0
end

function check_cache(backend)
    @testset "Inspection leaves prepared solves unchanged" begin
        r,(x,y,z)=pauli_variables(2)
        for symmetry in ((),(SitePermutation([2,1]),))
            built=N.build(Problem(-z[1]-z[2];symmetry),Relaxation(order=1))
            handle=N.prepare(built;solver=backend)
            before=N.moment_standard_form(built.program)
            collect(N.block_rows(built))
            @test N.prepare(built;solver=backend)===handle
            after=N.moment_standard_form(built.program)
            @test all(getproperty(before,f)==getproperty(after,f) for f in (:A,:b,:c,:offset,:cones))
        end
    end
end
end
