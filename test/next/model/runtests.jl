# Tests of the model layer of NCTSSoS.Next.

include("mixed.jl")
MixedModelChecks.check(Nx.JuMPBackend(SOLVER))

@testset "Public imports in fresh processes" begin
    project=dirname(Base.active_project())
    standalone=raw"""
        using Test
        using NCTSSoS.Next
        @testset "Every exported binding is usable" begin
            for name in names(Next)
                Base.isexported(Next,name) || continue
                @test Core.eval(Main,name) === getfield(Next,name)
            end
            _,parts=Next.mixed_variables(spin=(Pauli(),1),matter=(Fermionic(),1))
            @test parts.spin[1][1]*parts.matter[1][1]==parts.matter[1][1]*parts.spin[1][1]
            registry,(X,Y,Z)=pauli_variables(4)
            H=sum(X[i]*X[mod1(i+1,4)]+Y[i]*Y[mod1(i+1,4)]+Z[i]*Z[mod1(i+1,4)] for i in 1:4)
            translation=SitePermutation([2,3,4,1])
            reflection=SitePermutation([4,3,2,1])
            problem=Problem(H;symmetry=(translation,reflection,spin_rotations(registry)...))
            @test occursin("Symmetry: 2 block generators",sprint(show,problem))
            @test occursin("3 continuous generators",sprint(show,problem))
            @test !occursin("FamilySymmetry",sprint(show,problem))
            @test inspect(problem,Relaxation(order=1)).sizes.program.rows>0
        end
        import NCTSSoS
        @testset "Qualified legacy API" begin
            for name in (:Symmetry,:variables,:newton_chip_basis)
                @test Base.ispublic(Next,name)
                @test !Base.isexported(Next,name)
                @test getfield(NCTSSoS,name) !== getfield(Next,name)
            end
            for name in (:RealState,:RealTrace,:build,:BuiltRelaxation,:prepare,:prepare_standard,:solve_standard,:task_sizes,
                    :terms,:word_letters,:support_letters,:sites,:degree,:operator_name,:basis,:block_rows,:block_sides)
                @test Base.ispublic(Next,name)
                @test !Base.isexported(Next,name)
            end
        end
        """
    combined=raw"""
        using Test
        using NCTSSoS
        using NCTSSoS.Next
        @testset "No clashing exports" begin
            for name in (:MF,:MMD,:NoElimination,:AsIsElimination,:MaximalElimination,:NaturalElimination)
                @test Base.isexported(NCTSSoS.Next,name)
                @test Core.eval(Main,name) === getfield(NCTSSoS,name) === getfield(NCTSSoS.Next,name)
            end
            for name in (:gns,:evaluate,:bosonic_variables,:Bosonic,:state_symbol,:trace_symbol,:ExpectationPositive,:ExpectationZero)
                @test Base.isexported(NCTSSoS.Next,name)
                @test Core.eval(Main,name) === getfield(NCTSSoS.Next,name)
            end
            clashes=sort!([name for name in intersect(names(NCTSSoS),names(NCTSSoS.Next))
                if Base.isexported(NCTSSoS,name) && Base.isexported(NCTSSoS.Next,name) &&
                    getfield(NCTSSoS,name) !== getfield(NCTSSoS.Next,name)])
            @test isempty(clashes)
            @test Symmetry === NCTSSoS.Symmetry
            @test variables === NCTSSoS.variables
        end
        """
    for script in (standalone,combined)
        @test success(pipeline(`$(Base.julia_cmd()) --startup-file=no --project=$project -e $script`;stdout,stderr))
    end
end

include("front.jl")
FrontDoorChecks.check(Nx.JuMPBackend(SOLVER))
include("functional.jl")
FunctionalModelChecks.check(Nx.JuMPBackend(SOLVER))
include("newton_chip.jl")
NewtonChipChecks.check(Nx.JuMPBackend(SOLVER))
include("letter_permutation.jl")
include("inspection.jl")
InspectionChecks.check_cache(Nx.JuMPBackend(SOLVER))
include("staged.jl")
StagedSolveChecks.check(Nx.JuMPBackend(SOLVER))
include("expectation_constraints.jl")
ExpectationConstraintChecks.check(Nx.JuMPBackend(SOLVER))
include("gns.jl")
include("bosonic.jl")
# Order-three bosonic Hankel matrices need tighter solver residuals to compare
# the two lowerings within the unchanged 1e-6 bound tolerance.
let solver=optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-10,"eps_rel"=>1e-10,"max_iter"=>100_000)
    BosonModelChecks.check(Nx.JuMPBackend(solver),solver)
end
include("detection.jl")
isdefined(@__MODULE__,:check_expectations) || include("expectations.jl")
@testset "Reading expectations through JuMP" begin
    check_expectations(Nx.JuMPBackend(SOLVER),SOLVER)
end

@testset "Registry display and minimal relaxation storage" begin
    r = Nx.Registry(Nx.Pauli(),3,[:Alice=>(1,0),:Bob=>(2,0),:Carol=>(3,0)])
    v = Nx.variables(r)
    p = v[:Alice]*v[:Bob]*v[:Carol]
    @test sprint(show,p)=="Alice Bob Carol"
    @test sprint(show,p.data.words[1],r)=="Alice Bob Carol"
    @test sprint(show,last(only(Nx.terms(p)));context=:registry=>r)=="Alice Bob Carol"
    @test_throws r"Alice Bob Carol.*degree 2" Nx.build(Nx.Problem(p),Nx.Relaxation(order=1))
    mixed,parts=Nx.mixed_variables(spin=(Nx.Pauli(),1),sign=(Nx.SiteAlgebra(Nx.Unipotent(),1),1))
    X,Y,Z=parts.spin; u=parts.sign[1][1]
    phase_product=X[1]*Y[1]*u
    coefficient,w=only(Nx.terms(phase_product))
    @test signbit(real(coefficient))
    @test sprint(show,phase_product)=="(0.0 + 1.0im) spin.Z1 sign.s1_1"
    @test sprint(show,phase_product*Z[1]*u)=="0.0 + 1.0im"
    @test sprint(show,w;context=:registry=>mixed)=="spin.Z1 sign.s1_1"
    @test isequal(first(only(Nx.terms(phase_product))),coefficient)
    r,(a,ad) = Nx.fermionic_variables(2)
    @test sprint(show,ad[2]*a[1])=="ad2 a1"
    r,(x,y,z) = Nx.pauli_variables(130)
    @test sprint(show,x[1]*y[64]*z[130])=="X1 Y64 Z130"
    r,(a,ad) = Nx.fermionic_variables(130)
    @test sprint(show,ad[130]*ad[64]*a[64]*a[130])=="ad130 ad64 a64 a130"
    r,ops = Nx.site_variables(Nx.Free(),1,2)
    x,y = ops[1]
    built = Nx.build(Nx.Problem(x*y+y*x;constraints=(Nx.Positive(1-x^2-y^2),)),Nx.Relaxation(order=2))
    @test eltype(built.table.keys)===Nx.SiteWord{Nx.Free,UInt16,3,1}
    @test occursin("s1_1",sprint(show,built.valid_for))
    @test occursin("⪰ 0",sprint(show,built.valid_for))
    @test occursin("rows",sprint(show,built.sizes))
    r,(x,y) = Nx.site_variables(Nx.Unipotent(),2,2)
    bell = Nx.Problem(x[1]*y[1]+x[2]*y[2])
    @test eltype(Nx.build(bell,Nx.Relaxation(order=1)).table.keys)===Nx.SiteWord{Nx.Unipotent,UInt8,4,2}
    @test eltype(Nx.build(bell,Nx.Relaxation(rows=[one(x[1]),x[1],x[2],y[1],y[2]])).table.keys)===Nx.SiteWord{Nx.Unipotent,UInt8,4,2}
end

@testset "Model arithmetic and checked inputs" begin
    r,(x,y,z) = Nx.pauli_variables(1:3)
    @test x[1]*y[1]==im*z[1]
    @test x[1]^2==one(x[1])
    @test (x[1]*y[1])' == -im*z[1]
    @test 2+x[1]-2==x[1]
    @test iszero(x[1]-x[1])
    @test ishermitian(x[1]*x[2]+0.3*z[3])
    @test_throws ArgumentError Nx.Problem(x[1]*y[1])
    @test_throws ArgumentError Inf*x[1]
    @test_throws ArgumentError x[1]^(-1)
    r2,(x2,_,_) = Nx.pauli_variables(3)
    @test_throws ArgumentError x[1]+x2[1]
    @test_throws ArgumentError Nx.Problem(x[1];constraints=[Nx.Positive(x2[1])])
    @test_throws ArgumentError Nx.Problem(x[1];constraints=[true])
    @test_throws ArgumentError Nx.Relaxation(order=-1)
    @test_throws ArgumentError Nx.Relaxation(order=1,rows=[one(x[1])])
    @test_throws ArgumentError Nx.Relaxation(rows=[x[1]])
    @test_throws ArgumentError Nx.Relaxation(rows=[one(x[1]),x[1],x[1]])
    @test_throws ArgumentError Nx.Relaxation(rows=[one(x[1]),2*x[1]])
    @test_throws ArgumentError Nx.build(Nx.Problem(x[1]),Nx.Relaxation(rows=[one(x2[1]),x2[1]]))
    @test_throws ArgumentError Nx.build(Nx.Problem(x[1]*x[2]*x[3]),Nx.Relaxation(order=1))
    @test_throws r"using SCS.*SCSBackend.*JuMPBackend" Nx.solve(Nx.Problem(x[1]),Nx.Relaxation(order=1))
    custom = [x[1],one(x[1])]
    built = Nx.build(Nx.Problem(x[1]),Nx.Relaxation(rows=custom))
    @test built.sizes.rows==2 && custom==[x[1],one(x[1])]
    r,(a,ad) = Nx.fermionic_variables(3)
    @test ad[2]==a[2]'
    @test a[2]*ad[2]+ad[2]*a[2]==one(a[2])
    @test iszero(a[1]*a[3]+a[3]*a[1])
    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        r,ops = Nx.site_variables(rule,2,2)
        @test ops[1][1]*ops[2][2]==ops[2][2]*ops[1][1]
        @test ops[1][1]*ops[1][2] != ops[1][2]*ops[1][1]
        if rule isa Nx.Projector
            @test ops[1][1]^2==ops[1][1]
        elseif rule isa Nx.Unipotent
            @test ops[1][1]^2==one(ops[1][1])
        end
    end
end

@testset "One-call model solve equals frozen legacy bounds" begin
    backend = Nx.JuMPBackend(SOLVER)
    old = LEGACY_EXPECTED["heisenberg-open4-order1"]["objective"]
    r,(X,Y,Z) = Nx.pauli_variables(4)
    H = sum(X[i]*X[i+1]+Y[i]*Y[i+1]+Z[i]*Z[i+1] for i in 1:3)
    result = Nx.solve(Nx.Problem(H),Nx.Relaxation(order=1);solver=backend)
    @test result.bound ≈ old atol=2e-5
    @test result.status==JuMP.MOI.OPTIMAL
    @test result.bound_kind==:numerical_optimum_not_certificate
    @test occursin("numerical_optimum_not_certificate",sprint(show,result))
    @test occursin("status: OPTIMAL",sprint(show,result))
    @test result.sizes==Nx.inspect(result.problem,result.relaxation).sizes
    @test isempty(result.valid_for.constraints)
    @info "One call: Heisenberg4" old new=result.bound sizes=result.sizes

    # Two spatial sites, two spin orbitals each: an actual two-site Hubbard model.
    old = LEGACY_EXPECTED["hubbard2-order2"]["objective"]
    r,(a,ad) = Nx.fermionic_variables(4)
    H = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    result = Nx.solve(Nx.Problem(H),Nx.Relaxation(order=2);solver=backend)
    @test result.bound ≈ old atol=2e-5
    @test result.bound ≈ 1-sqrt(5) atol=2e-5
    @info "One call: Hubbard2" old new=result.bound sizes=result.sizes

    old = LEGACY_EXPECTED["chsh-order1"]["objective"]
    r,(x,y) = Nx.site_variables(Nx.Unipotent(),2,2)
    H = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    built = Nx.build(Nx.Problem(H),Nx.Relaxation(order=1))
    @test only(built.program.blocks).kind==Nx.PSD
    @test only(built.program.blocks).side==5
    @test built.program.nvars==11
    result = Nx.solve(Nx.Problem(H),Nx.Relaxation(order=1);solver=backend)
    @test result.bound ≈ old atol=2e-5
    @test result.bound ≈ -2sqrt(2) atol=2e-5
    @info "One call: CHSH" old new=result.bound sizes=result.sizes
end

@testset "Value hashing of polynomials and sizes" begin
    registries=[Nx.pauli_variables(2)[1],Nx.fermionic_variables(2)[1],Nx.bosonic_variables(2)[1],
        (Nx.site_variables(rule,1,2)[1] for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent()))...,
        (Nx.site_variables(Nx.Free(),1,2;functional=f)[1] for f in (Nx.RealState(),Nx.RealTrace()))...,
        Nx.mixed_variables(spin=(Nx.Pauli(),1),matter=(Nx.Fermionic(),1))[1]]
    for r in registries
        p=first(values(Nx.variables(r)))
        q=adjoint(adjoint(p))
        @test p==q && isequal(p,q)
        @test hash(p)==hash(q)
        @test length(unique([p,q,p+p,2p]))==2
        @test q in Set([p])
        @test Dict(p=>7)[q]==7
        W=typeof(last(only(Nx.terms(p))))
        # Construct canonical terms with distinct coefficient representations.
        for coefficients in ((1,1.0,1//1,ComplexF32(1),ComplexF64(1,-0.0)),
                (complex(0.0,1.0),complex(-0.0,1.0)))
            polys=map(coefficients) do c
                data=Nx.PolynomialBuffer{W,typeof(c)}(4)
                Nx.push_term!(data,c,last(only(Nx.terms(p))))
                Nx.Polynomial(r,Nx.canonicalize!(data))
            end
            @test all(q->q==first(polys),polys)
            @test all(q->hash(q,UInt(17))==hash(first(polys),UInt(17)),polys)
            @test length(Set(polys))==1
            @test all(q->haskey(Dict(first(polys)=>1),q),polys)
        end
        @test length(unique([zero(p),p-p]))==1
    end
    r,(letters,)=Nx.site_variables(Nx.Free(),1,2)
    a,b=letters
    rows=unique([a,a',a*b,a*b])
    @test length(rows)==2
    built=Nx.build(Nx.Problem(a^2),Nx.Relaxation(rows=[one(a);rows]))
    explicit=Nx.build(Nx.Problem(a^2),Nx.Relaxation(rows=[one(a),a,a*b]))
    f,g=Nx.moment_standard_form(built.program),Nx.moment_standard_form(explicit.program)
    @test (f.A,f.b,f.c,f.offset,f.cones)==(g.A,g.b,g.c,g.offset,g.cones)
    s,t=built.sizes,explicit.sizes
    @test s==t && hash(s)==hash(t)
    @test length(unique([s,t]))==1
    @test t in Set([s]) && Dict(s=>1)[t]==1
    _,(other,)=Nx.site_variables(Nx.Free(),1,2)
    @test a!=first(other) && length(Set([a,first(other)]))==2
end

@testset "Bulk polynomial construction" begin
    r,(x,y,z)=Nx.pauli_variables(3)
    ps=[3x[2],(2+3im)*x[1],-3x[2],x[3]*y[3],-im*z[3],2one(x[1])]
    saved=deepcopy.(collect.(Nx.terms.(ps)))
    expected=(2+3im)*x[1]+2
    @test Nx.polynomial(r,(p for p in ps))==expected
    @test sum(ps)==expected
    @test sum(ps;init=one(x[1]))==expected+1
    @test sum(ps;init=2)==expected+2
    @test collect.(Nx.terms.(ps))==saved
    @test Nx.polynomial(r,())==zero(x[1])
    @test Nx.polynomial(r,(zero(x[1]),))==zero(x[1])
    @test sum(ps[1:0];init=zero(x[1]))==zero(x[1])
    @test_throws r"registry" sum(ps[1:0])
    @test_throws r"polynomials" Nx.polynomial(r,(1,))
    @test_throws r"polynomials" Nx.polynomial(r,(x[1],1))
    _,(other,_,_)=Nx.pauli_variables(3)
    @test_throws r"registr" Nx.polynomial(r,(other[1],))
    @test_throws r"registr" sum([x[1],other[1]])
    @test_throws r"finite" Nx.polynomial(r,(c*x[1] for c in (1,NaN)))
    @test_throws r"finite" Nx.polynomial(r,(floatmax(Float64)*x[1],floatmax(Float64)*x[1]))
    # A one-pass iterator must not be restarted when the coefficient type widens.
    bigterm=big"1.000000000000000000001"*x[2]
    stream=Iterators.Stateful((x[1],bigterm,-x[1]))
    @test Nx.polynomial(r,stream)==bigterm
    @test isempty(stream)
    matrix=[x[1] x[2]; x[3] -x[1]]
    @test sum(matrix)==x[2]+x[3]
    @test sum(matrix;dims=1)==[x[1]+x[3] x[2]-x[1]]
    # Insert, combine and cancel at each boundary without changing the saved prefix.
    prefix=2x[2]; original=collect(Nx.terms(prefix))
    for term in (x[1],3x[2],-2x[2],x[3],zero(x[1]))
        @test prefix+term==Nx.polynomial(r,(prefix,term))
        @test prefix-term==Nx.polynomial(r,(prefix,-term))
        @test collect(Nx.terms(prefix))==original
    end
    fr,(a,ad)=Nx.fermionic_variables(2)
    @test Nx.polynomial(fr,(a[1]*ad[1],ad[1]*a[1],a[2]*a[2]))==one(a[1])
    # Multiple limbs and degree >255 exercise all radix-key fields and byte order.
    wide,(a,ad)=Nx.fermionic_variables(130)
    modes=(1,8,9,63,64,65,129,130)
    ps=vec([(i+2im)*ad[j]*a[k] for (i,(j,k)) in enumerate(Iterators.product(modes,reverse(modes)))])
    occupied=foldl(*,(ad[i]*a[i] for i in 1:130))
    append!(ps,(occupied,one(a[1]),-ps[3],2ps[7]))
    @test Nx.polynomial(wide,ps)==foldl(+,ps)
    @test Nx.polynomial(wide,reverse(ps))==foldl(+,reverse(ps))
    @test Nx.polynomial(wide,(occupied,-occupied))==zero(occupied)
    br,(b,bd)=Nx.bosonic_variables(1)
    @test Nx.polynomial(br,(b[1]*bd[1],-bd[1]*b[1]))==one(b[1])
    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        sr,(s,)=Nx.site_variables(rule,1,2)
        @test Nx.polynomial(sr,(s[1]*s[2],-s[2]*s[1]))==s[1]*s[2]-s[2]*s[1]
    end
    for (functional,symbol) in ((Nx.RealState(),Nx.state_symbol),(Nx.RealTrace(),Nx.trace_symbol))
        sr,(s,)=Nx.site_variables(Nx.Free(),1,2;functional)
        ps=(symbol(s[1]*s[2]),symbol(s[1])*symbol(s[2]),-symbol(s[2]*s[1]))
        @test Nx.polynomial(sr,ps)==foldl(+,ps)
    end
    @test Base.ispublic(Nx,:polynomial) && !Base.isexported(Nx,:polynomial)
end

@testset "Real form ignores unused polynomial capacity" begin
    _,(letters,)=Nx.site_variables(Nx.Free(),1,2)
    x,y=letters
    f=1.0*x*y*x
    g,h=1.0-x^2,1.0-y^2
    # Make the clean baseline deterministic, then contaminate each tail alone.
    for p in (f,g,h)
        Nx.reserve!(p.data,max(length(p.data.words),length(p.data)+1))
        fill!(view(p.data.coefficients,length(p.data)+1:length(p.data.coefficients)),0)
    end
    problem=Nx.Problem(f;constraints=(Nx.Positive(g),Nx.Positive(h)))
    relaxation=Nx.Relaxation(order=3)
    clean=Nx.build(problem,relaxation)
    @test eltype(clean.program.forms.coefs)===Float64
    @test clean.sizes.real_variables==77
    reference=Nx.moment_standard_form(clean.program)
    for p in (f,g)
        used=collect(Nx.terms(p))
        tail=length(p.data)+1:length(p.data.coefficients)
        fill!(view(p.data.coefficients,tail),1im)
        fill!(view(p.data.words,tail),last(only(Nx.terms(x))))
        @test collect(Nx.terms(p))==used
        built=Nx.build(problem,relaxation)
        @test eltype(built.program.forms.coefs)===Float64
        form=Nx.moment_standard_form(built.program)
        @test all(getproperty(form,k)==getproperty(reference,k) for k in (:A,:b,:c,:offset,:cones))
        result=Nx.solve(built;solver=Nx.JuMPBackend(SOLVER))
        @test result.sizes==clean.sizes
        @test result.bound≈-1 atol=2e-5
        fill!(view(p.data.coefficients,tail),0)
    end
    # A non-real live coefficient must still select the Hermitian form.
    complex=Nx.build(Nx.Problem(f+im*(x*y-y*x)),relaxation)
    @test eltype(complex.program.forms.coefs)===ComplexF64
end

@testset "Real site-algebra form and complex fallback" begin
    backend=Nx.JuMPBackend(SOLVER)
    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        _,ops=Nx.site_variables(rule,1,2)
        a,b=ops[1]
        real_build=Nx.build(Nx.Problem(a*b+b*a),Nx.Relaxation(order=1))
        @test only(real_build.program.blocks).kind==Nx.PSD
        complex_build=Nx.build(Nx.Problem(im*(a*b-b*a)),Nx.Relaxation(order=1))
        @test only(complex_build.program.blocks).kind==Nx.HPSD
        for constraint in (Nx.Positive(1+0.25im*(a*b-b*a)),Nx.Zero(im*(a*b-b*a)))
            constrained=Nx.build(Nx.Problem(a+b;constraints=(constraint,)),Nx.Relaxation(order=1))
            @test first(constrained.program.blocks).kind==Nx.HPSD
        end
    end
    _,ops=Nx.site_variables(Nx.Unipotent(),1,2)
    a,b=ops[1]
    @test Nx.solve(Nx.Problem(im*(a*b-b*a)),Nx.Relaxation(order=1);solver=backend).bound ≈ -2 atol=2e-5
    restricted=Nx.Problem(im*(a*b-b*a);constraints=(Nx.Positive(0.5+im*(a*b-b*a)),))
    @test Nx.solve(restricted,Nx.Relaxation(order=1);solver=backend).bound ≈ -0.5 atol=2e-5

    _,ops=Nx.site_variables(Nx.Free(),1,2)
    a,b=ops[1]
    h,g=a*b+b*a,1-a^2-b^2
    problem=Nx.Problem(h;constraints=(Nx.Positive(g),))
    built=Nx.build(problem,Nx.Relaxation(order=2))
    @test all(block.kind==Nx.PSD for block in built.program.blocks)
    # Independent complex reference, with the same rows and localizing polynomial.
    W=Nx.word_type(h.registry,4)
    packed(w)=Nx.word_from_legacy(W,Nx.legacy_letters(w))
    rows=Nx.basis_rows(W,h.registry.algebra,1:1,2)
    reference,table=Nx.plain_program(rows,[(c,packed(w)) for (c,w) in Nx.terms(h)])
    sink=Nx.MomentSink(table,reference.forms)
    local_rows=[one(a),a,b]
    for j in 1:3,i in j:3
        for (c,w) in Nx.terms(local_rows[i]'*g*local_rows[j])
            Nx.push_term!(sink,c,packed(w))
        end
        Nx.close_form!(reference.forms)
    end
    Nx.push_block!(reference,Nx.HPSD,3,1)
    reference.nvars=Nx.nvars(table)
    complex=Nx.solve_standard(Nx.moment_standard_form(reference),backend)
    real=Nx.solve(problem,Nx.Relaxation(order=2);solver=backend)
    @test real.bound ≈ -1 atol=2e-5
    @test real.bound ≈ complex.bound atol=2e-5
end

@testset "Localizers and re-solving an open backend" begin
    r,ops = Nx.site_variables(Nx.Free(),1,1)
    x = ops[1][1]
    p = Nx.Problem(x;constraints=(Nx.Positive(1-x^2),))
    result = Nx.solve(p,Nx.Relaxation(order=2);solver=Nx.JuMPBackend(SOLVER))
    @test result.bound ≈ -1 atol=2e-5
    @test result.sizes.localizing_rows==[2]
    limited = Nx.solve(result,Nx.Positive(x-0.25))
    @test limited.bound ≈ 0.25 atol=2e-5
    @test limited.solution.handle === result.solution.handle
    @test length(limited.valid_for.constraints)==2
    fixed = Nx.solve(limited,Nx.Zero(x-0.6))
    @test fixed.bound ≈ 0.6 atol=2e-5
    @test fixed.solution.handle === result.solution.handle
    @test Nx.expectation(fixed,x) ≈ 0.6 atol=2e-5
    @test Nx.expectation(result,x) ≈ -1 atol=2e-5
    failed = Nx.solve(fixed,Nx.Positive(x-0.8))
    @test failed.bound === nothing
    @test failed.bound_kind==:unavailable
    @test failed.status in (JuMP.MOI.INFEASIBLE,JuMP.MOI.INFEASIBLE_OR_UNBOUNDED)
    @test_throws r"INFEASIBLE.*no accepted solution" Nx.expectation(failed,x)
    @test_throws ArgumentError Nx.build(Nx.Problem(x;constraints=(Nx.Positive(1-x^4),)),Nx.Relaxation(order=1))
    # Direct use of both standard forms shares the same backend entry point.
    built = Nx.build(p,Nx.Relaxation(order=2))
    sos = Nx.solve_standard(Nx.sos_standard_form(built.program),Nx.JuMPBackend(SOLVER))
    @test sos.bound ≈ -1 atol=2e-5
end
