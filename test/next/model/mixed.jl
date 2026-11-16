module MixedModelChecks
using Test, NCTSSoS, LinearAlgebra
import NCTSSoS as N

function check(solver)
    @testset "Mixed modeling and joint moments" begin
        @test Base.ispublic(N,:mixed_variables) && !Base.isexported(N,:mixed_variables)
        for specs in ((;), (;spin=(N.Pauli(),0)), (;spin=(N.Pauli,1)),
                (;a=(N.Fermionic(),1),b=(N.Fermionic(),1)),
                (;state=(N.FunctionalAlgebra(N.RealState(),N.SiteAlgebra(N.Free(),2)),1)))
            @test_throws ArgumentError N.mixed_variables(;specs...)
        end
        r,p=N.mixed_variables(spin=(N.Pauli(),1),matter=(N.Fermionic(),2))
        X,Y,Z=p.spin; a,ad=p.matter
        n1,n2=ad[1]*a[1],ad[2]*a[2]
        @test X[1]*a[1]==a[1]*X[1]
        @test a[1]*ad[2]==-ad[2]*a[1]
        @test X[1]*Y[1]==im*Z[1]
        @test a[1]*ad[1]+n1==one(X[1])
        @test N.degree(X[1]*n1)==3
        @test N.variables(r)[Symbol("spin.X1")]==X[1]
        h=Z[1]*(n1-n2)+X[1]*(ad[1]*a[2]+ad[2]*a[1])
        @test occursin("matter.ad",string(h))
        rebuilt=N.polynomial(r,(c*foldl(*,(N.variables(r)[N.operator_name(r,l)] for l in N.word_letters(w));init=one(X[1])) for (c,w) in N.terms(h)))
        @test rebuilt==h
        # Independent Jordan–Wigner matrices; a Bell eigenstate attains -2,
        # whereas product spin/matter states cannot go below -1.
        sx=ComplexF64[0 1;1 0]; sz=ComplexF64[1 0;0 -1]
        id=Matrix{ComplexF64}(I,2,2); annih=ComplexF64[0 1;0 0]
        c1=kron(annih,id); c2=kron(sz,annih)
        matrix=kron(sz,c1'*c1-c2'*c2)+kron(sx,c1'*c2+c2'*c1)
        exact=eigmin(Hermitian(matrix))
        @test exact≈-2 atol=1e-14
        built=N.build(N.Problem(h),N.Relaxation(order=2))
        @test built.sizes.rows==26
        @test N.block_sides(built)==[26]
        @test length(only(N.block_rows(built)))==26
        custom=N.build(N.Problem(h),N.Relaxation(rows=N.basis(r,2)))
        f,g=N.moment_standard_form(built.program),N.moment_standard_form(custom.program)
        @test (f.A,f.b,f.c,f.offset,f.cones)==(g.A,g.b,g.c,g.offset,g.cones)
        prepared=N.prepare(built;solver)
        @test N.prepare(built;solver)===prepared
        result=N.solve(built;solver)
        @test result.bound≈exact atol=2e-5
        @test result.bound < -1.9
        @test N.expectation(result,Z[1]*(n1-n2))≈-1 atol=2e-5
        @test N.expectation(result,X[1]*(ad[1]*a[2]+ad[2]*a[1]))≈-1 atol=2e-5
        @test_throws r"expectation" N.gns(result)
        @test_throws r"remove symmetry" N.Problem(h;symmetry=(N.SitePermutation([1,2,3]),))
        for relaxation in (N.Relaxation(order=2;grouping=N.Correlative()),
                N.Relaxation(order=2;strengthen=N.FermionicConditions(2:3)),
                N.Relaxation(family=[collect(N.basis(r,2))]))
            @test_throws r"plain route" N.build(N.Problem(h),relaxation)
        end
        _,q=N.mixed_variables(spin=(N.Pauli(),1),matter=(N.Fermionic(),1))
        x,_,z=q.spin; b,bd=q.matter; n=bd[1]*b[1]
        energy=n/2+z[1]+2x[1]*n
        exact=eigmin(Hermitian(kron(id,annih'*annih)/2+kron(sz,id)+2kron(sx,annih'*annih)))
        @test exact≈0.5-sqrt(5) atol=1e-14
        @test N.solve(N.Problem(energy),N.Relaxation(order=2);solver).bound≈exact atol=2e-5
        @test_throws r"different registries" energy+h
        _,q=N.mixed_variables(spin=(N.Pauli(),1),osc=(N.Bosonic(),1))
        _,_,z=q.spin; b,bd=q.osc
        energy=2bd[1]*b[1]+z[1]*(b[1]+bd[1])+z[1]/4
        # 2(b+Z/2)'(b+Z/2)+Z/4-1/2 >= -3/4, attained by a coherent state.
        @test N.solve(N.Problem(energy),N.Relaxation(order=2);solver).bound≈-0.75 atol=2e-5
        _,q=N.mixed_variables(free=(N.SiteAlgebra(N.Free(),2),1),sign=(N.SiteAlgebra(N.Unipotent(),1),1))
        u,v=q.free[1]; t=q.sign[1][1]
        problem=N.Problem(u*t;constraints=(N.Positive(1-u^2),N.Zero(v)))
        built=N.build(problem,N.Relaxation(order=2))
        @test all(b->b.kind!=N.HPSD,built.program.blocks)
        @test N.solve(built;solver).bound≈-1 atol=2e-5
    end
end
end
