# Shared observable contract for the JuMP and native SCS routes.
function check_expectations(backend,reference_optimizer)
    N=NCTSSoS.Next
    _,(a,b)=N.site_variables(N.Unipotent(),2,2)
    correlators=(-a[1]*b[1],-a[1]*b[2],-a[2]*b[1],a[2]*b[2])
    bell=N.solve(N.Problem(sum(correlators)),N.Relaxation(order=1);solver=backend)
    @test sum(N.expectation(bell,p) for p in correlators) ≈ bell.bound atol=2e-5
    @test N.expectation(bell,first(correlators)) isa Real
    @test N.expectation(bell,one(a[1])) ≈ 1 atol=2e-5
    @test N.expectation(bell,zero(a[1])) === 0.0
    @test hasproperty(bell,:bound_kind)
    @test N.expectation(bell,a[1]*a[2]) isa Complex
    @test imag(N.expectation(bell,a[1]*a[2]))==0
    @test N.expectation(bell,a[1]*a[2]) ≈ N.expectation(bell,a[2]*a[1]) atol=2e-5
    # Three letters do not fit the order-1 table's narrow site-word storage.
    @test_throws r"s1_1.*s1_2.*s2_1.*order=2.*custom rows" N.expectation(bell,a[1]*a[2]*b[1])
    _,(other,_)=N.site_variables(N.Unipotent(),2,2)
    @test_throws r"different registry" N.expectation(bell,other[1])

    _,(x,y,z)=N.pauli_variables(4)
    h=sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
    problem=N.Problem(h)
    relaxation=N.Relaxation(order=1)
    result=N.solve(problem,relaxation;solver=backend)
    built=N.build(problem,relaxation)
    model,variables=N.to_jump(built.program)
    set_optimizer(model,reference_optimizer)
    optimize!(model)
    @test termination_status(model)==JuMP.MOI.OPTIMAL
    # Read a variable directly from the independent, expression-based JuMP path.
    id=only(i for (i,label) in enumerate(built.table.labels) if N.legacy_letters(label.word)==[1,4])
    @test N.expectation(result,x[1]*x[2]) ≈ value(variables[id]) atol=2e-5
    @test N.expectation(result,im*h) isa Complex
    @test N.expectation(result,im*h) ≈ im*result.bound atol=2e-5

    # A split can omit a degree-one moment even at order two. Increasing the
    # plain Relaxation order is not the right advice for this result's route.
    ring=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+1.5*z[i]*z[mod1(i+1,4)] for i in 1:4)
    blocks=N.BlockRelaxation(order=2,grouping=N.Correlative(),split=N.TermSplit(),
        symmetry=N.FamilySymmetry([N.SitePermutation([2,3,4,1])]))
    blocked=N.solve(N.Problem(ring),blocks;solver=backend)
    @test_throws r"Moment Z2.*selected blocks.*higher order.*no split.*group.*2" N.expectation(blocked,z[2])

    _,(a,ad)=N.fermionic_variables(2)
    t=ad[1]*a[2]
    result=N.solve(N.Problem(-((1+2im)*t+(1-2im)*t')),N.Relaxation(order=1);solver=backend)
    # The one-particle Hamiltonian has eigenvalues ±sqrt(5); its minimizing
    # eigenvector gives <ad1*a2> = (1-2im)/(2sqrt(5)).
    @test result.bound ≈ -sqrt(5) atol=2e-5
    @test N.expectation(result,t) isa Complex
    @test N.expectation(result,t) ≈ (1-2im)/(2sqrt(5)) atol=2e-5
    @test N.expectation(result,t') ≈ (1+2im)/(2sqrt(5)) atol=2e-5
    @test N.expectation(result,t+t') isa Real
    @test N.expectation(result,t+t') ≈ 1/sqrt(5) atol=2e-5
    @test N.expectation(result,result.problem.objective) ≈ result.bound atol=2e-5

    registry=N.Registry(N.Pauli(),4,[:Alice=>(1,0),:Bob=>(2,0),:Carol=>(3,0),:Dave=>(4,0)])
    ops=N.variables(registry)
    alice,bob=ops[:Alice],ops[:Bob]
    left,right=alice*bob,ops[:Carol]*ops[:Dave]
    custom=N.solve(N.Problem(left*right),N.Relaxation(rows=[one(alice),left,right]);solver=backend)
    @test N.expectation(custom,left*right) ≈ -1 atol=2e-5
    # Suggested complete rows must cover the existing objective too, not only Bob.
    @test_throws r"Moment Bob.*order=2.*custom rows" N.expectation(custom,bob)
end
