# Tests of the model layer of NCTSSoS.Next.

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
    @test_throws ArgumentError Nx.solve(Nx.Problem(x[1]),Nx.Relaxation(order=1))
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
    result = Nx.solve(Nx.Problem(H),Nx.Relaxation(order=1);solver=backend)
    @test result.bound ≈ old atol=2e-5
    @test result.bound ≈ -2sqrt(2) atol=2e-5
    @info "One call: CHSH" old new=result.bound sizes=result.sizes
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
    failed = Nx.solve(fixed,Nx.Positive(x-0.8))
    @test failed.bound === nothing
    @test failed.bound_kind==:unavailable
    @test failed.status in (JuMP.MOI.INFEASIBLE,JuMP.MOI.INFEASIBLE_OR_UNBOUNDED)
    @test_throws ArgumentError Nx.build(Nx.Problem(x;constraints=(Nx.Positive(1-x^4),)),Nx.Relaxation(order=1))
    # Direct use of both standard forms shares the same backend entry point.
    built = Nx.build(p,Nx.Relaxation(order=2))
    sos = Nx.solve_standard(Nx.sos_standard_form(built.program),Nx.JuMPBackend(SOLVER))
    @test sos.bound ≈ -1 atol=2e-5
end
