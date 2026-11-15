# This also runs standalone after combination_matrix.jl and TestUtils.jl.
@testset "Antiunitary moment constraints without a real sector basis" begin
    C=CombinationMatrix
    for n in (4,6)
        f=C.fixture(:xxz,n;anisotropy=0.5)
        t=Nx.SitePermutation([2:n;1])
        K=Nx.AntiunitaryPauli(one(eltype(f.H.data.words)))
        flip=Nx.CliffordMap(Nx.PauliConjugation(first(prod(f.operators[1:n]).data.words),false))
        sz=Nx.SpinRotation(f.registry,:z)
        declarations=((t,K),(flip,K),(t,flip),(t,flip,K),(flip,t,K),
            (t,flip,sz),(t,flip,K,sz),(t,K,sz))
        for order in (1,2)
            r=Nx.Relaxation(;order)
            for generators in declarations
                p=Nx.Problem(f.H;symmetry=generators)
                built=Nx.build(p,r)
                checks=C.check_program(f,p,r,built)
                @test checks.state==:pass
                @test checks.spectral_error<1e-7
                @test checks.equality_error<1e-7
            end
            plain=Nx.Problem(f.H;symmetry=(t,))
            p=Nx.Problem(f.H;symmetry=(t,K))
            built=Nx.build(p,r)
            @test any(b->b.kind==Nx.HPSD,built.program.blocks)
            @test any(b->b.kind==Nx.PSD,built.program.blocks)
            # Keeping a complex sector must not discard K's moment relations.
            odd=Nx.moment_terms(built.table,first(f.operators[n+1].data.words))
            even=Nx.moment_terms(built.table,first((f.operators[n+1]*f.operators[n+2]).data.words))
            @test !isnothing(odd) && all(iszero(c) for (_,c) in odd)
            @test !isnothing(even) && any(!iszero(c) for (_,c) in even)
            if n==4 && order==1
                rows=only(built.sizes.closure.family).rows
                program,table,_=Nx.abelian_lattice_program(rows,collect(Nx.terms(f.H)),
                    Nx.TranslationGroup([t],[n]),[one(t)];moment_generators=(K,))
                checks=C.check_program(f,p,r,(;program,table,sizes=built.sizes))
                @test checks.state==:pass
                @test any(b->b.kind==Nx.HPSD,program.blocks)
            end
            solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
                "decompose"=>false,"eps_abs"=>1e-9,"eps_rel"=>1e-9))
            baseline=Nx.solve(plain,r;solver)
            result=Nx.solve(p,r;solver)
            @test !isnothing(result.bound) && !isnothing(baseline.bound)
            @test result.bound≈baseline.bound atol=1e-6 rtol=0
            @test Nx.expectation(result,f.H)≈result.bound atol=1e-8 rtol=0
        end
    end
end
