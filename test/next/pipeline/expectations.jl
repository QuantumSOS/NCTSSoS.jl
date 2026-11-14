@testset "Block-family expectations against exact states" begin
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-10,"eps_rel"=>1e-10,"max_iter"=>100_000))
    _,(x,y,z)=Nx.pauli_variables(4)
    H=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+0.5*z[i]*z[mod1(i+1,4)] for i in 1:4)
    matrix(p)=sum(c*physical_pauli(w,1:4) for (c,w) in Nx.terms(p))
    exact=eigen(Hermitian(matrix(H)))
    @test exact.values[2]-exact.values[1]>1e-3
    ground=exact.vectors[:,1]
    symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,3,4,1]),Nx.SitePermutation([1,4,3,2])])
    for description in (Nx.BlockRelaxation(order=2,symmetry=symmetry),
        Nx.BlockRelaxation(order=1,symmetry=symmetry,strengthen=(Nx.DensityMatrix(1:4),)))
        # The numerical density-sector projection has a 1e-9 residual
        # tolerance; demanding 1e-10 from COSMO stalls on that program.
        density=!isempty(description.strengthen)
        backend=density ? Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,
            "verbose"=>false,"decompose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8)) : solver
        result=Nx.solve(Nx.Problem(H),description;solver=backend)
        @test result.sizes.closure.route==:lattice
        @test result.bound≈exact.values[1] atol=1e-6
        for p in (x[1]*x[2],z[1]*z[2],x[1]*x[3],prod(z),x[1]*y[1]*z[2])
            @test Nx.expectation(result,p)≈dot(ground,matrix(p)*ground) atol=(density ? 5e-6 : 1e-6)
        end
        @test Nx.expectation(result,x[1]*x[2])≈Nx.expectation(result,x[3]*x[4]) atol=1e-12
        @test Nx.expectation(result,H)≈result.bound atol=1e-7
        @test Nx.expectation(result,one(H))==1
    end

    W=typeof(first(x[1].data.words))
    anti=Nx.AntiunitaryPauli(one(W))
    result=Nx.solve(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry(;moment_only=(anti,)));solver)
    before=(length(result.table.labels),length(result.table.cache))
    @test Nx.expectation(result,y[1])==0
    # The degree-three word was never built, but conjugation forces it to zero.
    @test Nx.expectation(result,y[1]*z[2]*z[3])==0
    @test_throws r"Moment Z1 Z2 Z3.*selected blocks.*higher order.*no split.*group.*1, 2, 3" Nx.expectation(result,z[1]*z[2]*z[3])
    @test before==(length(result.table.labels),length(result.table.cache))
end

@testset "Real-form and complex fermionic expectations" begin
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8))
    _,ops=Nx.site_variables(Nx.Unipotent(),3,2)
    H=sum(ops[i][1]*ops[mod1(i+1,3)][1]+0.3ops[i][2] for i in 1:3)
    result=Nx.solve(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,3,1]) ]));solver)
    @test result.table.real_form
    @test result.sizes.closure.route==:block_stabilizers
    @test Nx.expectation(result,H)≈result.bound atol=1e-7
    @test Nx.expectation(result,ops[1][1]*ops[2][1])≈Nx.expectation(result,ops[2][1]*ops[3][1]) atol=1e-12
    @test Nx.expectation(result,ops[1][1]*ops[1][2]) isa Complex
    @test imag(Nx.expectation(result,ops[1][1]*ops[1][2]))==0
    @test_throws r"Moment s1_1 s1_2 s2_1.*selected blocks.*higher order.*no split.*group.*1, 2" Nx.expectation(result,ops[1][1]*ops[1][2]*ops[2][1])

    _,(a,ad)=Nx.fermionic_variables(2)
    t=ad[1]*a[2]
    H=-((1+2im)*t+(1-2im)*t')
    result=Nx.solve(Nx.Problem(H),Nx.BlockRelaxation(order=1);solver)
    @test Nx.expectation(result,t)≈(1-2im)/(2sqrt(5)) atol=1e-6
    @test Nx.expectation(result,t')≈(1+2im)/(2sqrt(5)) atol=1e-6
    @test Nx.expectation(result,H)≈result.bound atol=1e-7
end
