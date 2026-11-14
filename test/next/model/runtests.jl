# Tests of the model layer of NCTSSoS.Next.

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
    @test_throws r"Alice Bob Carol.*degree 2" Nx.build(Nx.Problem(p),Nx.Relaxation(order=1))
    r,(a,ad) = Nx.fermionic_variables(2)
    @test sprint(show,ad[2]*a[1])=="ad2 a1"
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

@testset "One-call model solve equals current package" begin
    backend = Nx.JuMPBackend(SOLVER)
    reg,(x,y,z) = create_pauli_variables(1:4)
    oldH = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
    old = cs_nctssos(polyopt(oldH,reg),SolverConfig(optimizer=SOLVER,order=1);dualize=false).objective
    r,(X,Y,Z) = Nx.pauli_variables(4)
    H = sum(X[i]*X[i+1]+Y[i]*Y[i+1]+Z[i]*Z[i+1] for i in 1:3)
    result = Nx.solve(Nx.Problem(H),Nx.Relaxation(order=1);solver=backend)
    @test result.bound ≈ old atol=2e-5
    @test result.status==JuMP.MOI.OPTIMAL
    @test result.bound_kind==:numerical_optimum_not_certificate
    @test occursin("numerical_optimum_not_certificate",sprint(show,result))
    @test occursin("status: OPTIMAL",sprint(show,result))
    @test result.sizes==Nx.inspect(result.problem,result.relaxation)
    @test isempty(result.valid_for.constraints)
    @info "One call: Heisenberg4" old new=result.bound sizes=result.sizes

    # Two spatial sites, two spin orbitals each: an actual two-site Hubbard model.
    reg,(a,ad) = create_fermionic_variables(1:4)
    oldH = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    old = cs_nctssos(polyopt(oldH,reg),SolverConfig(optimizer=SOLVER,order=2);dualize=false).objective
    r,(a,ad) = Nx.fermionic_variables(4)
    H = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    result = Nx.solve(Nx.Problem(H),Nx.Relaxation(order=2);solver=backend)
    @test result.bound ≈ old atol=2e-5
    @test result.bound ≈ 1-sqrt(5) atol=2e-5
    @info "One call: Hubbard2" old new=result.bound sizes=result.sizes

    reg,(x,y) = create_unipotent_variables([("x",1:2),("y",1:2)])
    oldH = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    old = cs_nctssos(polyopt(oldH,reg),SolverConfig(optimizer=SOLVER,order=1);dualize=false).objective
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
