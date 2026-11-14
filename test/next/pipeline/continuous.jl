# Contract fixture; the production particle-number derivation belongs to core.
struct TestNumberDerivation <: Nx.AbstractDerivation end
function Nx.derive!(sink,::TestNumberDerivation,w::Nx.FermiWord,c::Number=1)
    charge=sum((l<0 ? 1 : -1 for l in Nx.legacy_letters(w));init=0)
    iszero(charge) || Nx.push_term!(sink,im*charge*c,w)
    sink
end

function continuous_true_result(built,problem,description,rho,matrix)
    values=[label.part==2 ? imag(LinearAlgebra.tr(rho*matrix(label.word))) : real(LinearAlgebra.tr(rho*matrix(label.word)))
        for label in built.table.labels]
    solution=Nx.ConicSolution(0.0,MOI.OPTIMAL,"invariant true-state fixture",values[2:end],nothing)
    result=Nx.SolveResult(0.0,MOI.OPTIMAL,built.sizes,built.valid_for,problem,description,solution,built.table,nothing)
    # Verify every old orbit coordinate, not just the objective's coordinates.
    for (id,label) in enumerate(built.table.base.labels)
        expected=label.part==2 ? imag(LinearAlgebra.tr(rho*matrix(label.word))) : real(LinearAlgebra.tr(rho*matrix(label.word)))
        @test sum((c*values[j] for (j,c) in built.table.coordinates[id]);init=0.0)≈expected atol=1e-11
    end
    for block in built.program.blocks
        if block.kind==Nx.ZERO
            @test maximum((abs(sum((built.program.forms.coefs[t]*values[built.program.forms.ids[t]]
                for t in Nx.term_range(built.program.forms,f));init=0.0)) for f in Nx.form_range(block));init=0.0)<1e-10
        end
    end
    psd=deepcopy(built.program)
    filter!(b->b.kind!=Nx.ZERO,psd.blocks)
    @test minimum(extra_values(psd,values,1))>=-1e-10
    result
end

@testset "Continuous moments: Heisenberg and finite symmetries" begin
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-9,"eps_rel"=>1e-9))
    for n in (4,6),finite in (false,true)
        _,(x,y,z)=Nx.pauli_variables(n)
        H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
        ds=[Nx.SpinRotation(axis,(UInt64(2^n-1),)) for axis in 0:2]
        generators=finite ? [Nx.SitePermutation([2:n;1]),Nx.SitePermutation([1;n:-1:2])] : Nx.AbstractAction[]
        plain=Nx.Problem(H;symmetry=Nx.Symmetry(generators...))
        invariant=Nx.Problem(H;symmetry=Nx.Symmetry(generators...;continuous=ds))
        description=Nx.Relaxation(order=1)
        base=Nx.build(plain,description)
        built=Nx.build(invariant,description)
        @test built.program.nvars<base.program.nvars
        @test maximum(b.side for b in built.program.blocks)<maximum(b.side for b in base.program.blocks)
        @test built.sizes.closure.continuous.method==:exact_elimination
        report=built.sizes.closure.continuous
        @test occursin("continuous: $(report.equations) independent generator equations, $(report.initial_variables) → $(report.variables) variables",sprint(show,Nx.inspect(invariant,description)))
        baseline=Nx.solve(plain,description;solver)
        result=Nx.solve(invariant,description;solver)
        @test result.status==baseline.status==MOI.OPTIMAL
        @test result.bound≈baseline.bound atol=1e-6
        @test Nx.expectation(result,H)≈result.bound atol=1e-6
        @test Nx.expectation(result,z[1])==0
        @test Nx.expectation(result,x[1]*x[2])==Nx.expectation(result,z[1]*z[2])
        matrix(w)=physical_pauli(w,1:n)
        h=sum(c*matrix(w) for (c,w) in Nx.terms(H))
        rho=exp(-0.23*h); rho/=LinearAlgebra.tr(rho)
        true_result=continuous_true_result(built,invariant,description,rho,matrix)
        @test Nx.expectation(true_result,x[1]*x[2])≈real(LinearAlgebra.tr(rho*matrix(last(only(Nx.terms(x[1]*x[2])))))) atol=1e-11
        before=(length(built.table.labels),length(built.table.base.cache))
        Nx.expectation(true_result,z[1])
        @test before==(length(built.table.labels),length(built.table.base.cache))
        @info "Continuous Heisenberg" n finite before=base.program.nvars-1 after=built.program.nvars-1 bound=result.bound
    end
end

@testset "Continuous diagonal labels equal general elimination" begin
    registry,(a,ad)=Nx.fermionic_variables(4)
    N=sum(ad[i]*a[i] for i in 1:4)
    H=sum((1+0.2im)*ad[i]*a[mod1(i+1,4)]+(1-0.2im)*ad[mod1(i+1,4)]*a[i]+
        0.4ad[i]*a[i]*ad[mod1(i+1,4)]*a[mod1(i+1,4)] for i in 1:4)
    d=TestNumberDerivation()
    for finite in (false,true)
        generators=finite ? [Nx.SitePermutation([2,3,4,1])] : Nx.AbstractAction[]
        description=Nx.Relaxation(order=1,strengthen=(Nx.FermionicConditions(1:4;conditions=(:P,:Q,:G,:T1,:T2prime)),Nx.Sector(N,2)))
        plain=Nx.Problem(H;symmetry=Nx.Symmetry(generators...))
        invariant=Nx.Problem(H;symmetry=Nx.Symmetry(generators...;continuous=(d,)))
        base=Nx.build(plain,description)
        built=Nx.build(invariant,description)
        general,_,report=Nx._continuous_program(deepcopy(base.program),deepcopy(base.table),(d,);diagonal=false)
        family_blocks=base.sizes.closure.strengthening.first_block-1
        general,_=Nx._continuous_blocks(general;family_blocks)
        @test built.sizes.closure.continuous.method==:diagonal
        @test report.method==:exact_elimination
        @test built.program.nvars<base.program.nvars
        @test built.program.nvars==general.nvars
        @test built.program.forms.ids==general.forms.ids
        @test built.program.forms.offsets==general.forms.offsets
        @test built.program.forms.coefs==general.forms.coefs
        @test built.program.objective_ids==general.objective_ids
        @test built.program.objective_coefs==general.objective_coefs
        production=Nx.build(Nx.Problem(H;symmetry=Nx.Symmetry(generators...;continuous=(Nx.NumberRotation(registry),))),description).program
        @test production.nvars==built.program.nvars
        @test production.forms.ids==built.program.forms.ids
        @test production.forms.offsets==built.program.forms.offsets
        @test production.forms.coefs==built.program.forms.coefs
        @test production.objective_ids==built.program.objective_ids
        @test production.objective_coefs==built.program.objective_coefs
        @test production.blocks==built.program.blocks
        solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
            "decompose"=>false,"eps_abs"=>1e-9,"eps_rel"=>1e-9))
        baseline=Nx.solve(plain,description;solver)
        result=Nx.solve(invariant,description;solver)
        @test result.status==baseline.status==MOI.OPTIMAL
        @test result.bound≈baseline.bound atol=1e-6
        @test Nx.expectation(result,a[1]*a[2]*a[3])==0
        @test Nx.expectation(result,H)≈result.bound atol=1e-6
        matrices=fermi_matrices(4)
        sector=findall(v->abs(v-2)<1e-12,diag(matrices.polynomial(N)))
        rho=zeros(ComplexF64,16,16)
        rho[sector,sector]=exp(-0.31*matrices.polynomial(H)[sector,sector]); rho/=LinearAlgebra.tr(rho)
        true_result=continuous_true_result(built,invariant,description,rho,matrices.word)
        @test Nx.expectation(true_result,ad[1]*a[2])≈LinearAlgebra.tr(rho*matrices.polynomial(ad[1]*a[2])) atol=1e-11
        @info "Continuous fermions" finite before=base.program.nvars-1 after=built.program.nvars-1 bound=result.bound
    end
end

@testset "Continuous guards and family closure" begin
    _,(x,y,z)=Nx.pauli_variables(2)
    H=x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    ds=[Nx.SpinRotation(axis,(UInt64(3),)) for axis in 0:2]
    problem=Nx.Problem(H;symmetry=Nx.Symmetry(continuous=ds))
    @test_throws r"both.*keep it on Problem" Nx.build(problem,Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,1]) ])))
    @test_throws r"derivation 1.*objective.*exactly" Nx.build(Nx.Problem(H+1e-15z[1];symmetry=Nx.Symmetry(continuous=ds)),Nx.Relaxation(order=1))
    @test_throws r"derivation 1.*constraint 1" Nx.build(Nx.Problem(H;constraints=(Nx.Positive(1+z[1]),),symmetry=Nx.Symmetry(continuous=ds)),Nx.Relaxation(order=1))
    @test_throws r"derivation 1.*Sector operator" Nx.build(problem,Nx.Relaxation(order=1,strengthen=Nx.Sector(z[1]+z[2],0)))
    @test_throws r"Custom block 1.*continuous derivation" Nx.build(Nx.Problem(zero(H);symmetry=Nx.Symmetry(continuous=ds)),
        Nx.Relaxation(family=[[one(H),x[1]]]))
    W=typeof(first(x[1].data.words))
    # Conjugation by a local Z changes only one site's X rotation; the global
    # X generator's span is therefore not preserved although the zero H is.
    flip=Nx.PauliConjugation(Nx.pauli_letter(W,1,2),false)
    @test_throws r"Finite generator 1.*add the missing derivations or remove" Nx.build(
        Nx.Problem(zero(H);symmetry=Nx.Symmetry(flip;continuous=(ds[1],))),Nx.Relaxation(order=1))
    # A normalizer need not commute: global X sends the Z generator to -Z.
    global_flip=Nx.PauliConjugation(first((x[1]*x[2]).data.words),false)
    axial=Nx.Problem(H;symmetry=Nx.Symmetry(global_flip;continuous=(ds[3],)))
    normalized=Nx.build(axial,Nx.Relaxation(order=1))
    # Only identity, XX=YY, and ZZ remain; in particular D(Y)=−X is enforced.
    @test normalized.program.nvars==3
    @test isempty(Nx.moment_terms(normalized.table,first(x[1].data.words)))
    built=Nx.build(problem,Nx.Relaxation(order=1,split=Nx.TermSplit()))
    for block in built.sizes.closure.family,w in block.rows,d in ds
        @test all(v->v in block.rows,keys(Nx._derivative(d,w)))
    end
    # Local density and stationarity cones remain valid on a thermal state.
    strengthened=Nx.build(problem,Nx.Relaxation(order=1,strengthen=(Nx.DensityMatrix(1:2),Nx.Stationarity())))
    matrix(w)=physical_pauli(w,1:2)
    rho=exp(-0.2sum(c*matrix(w) for (c,w) in Nx.terms(H))); rho/=LinearAlgebra.tr(rho)
    continuous_true_result(strengthened,problem,Nx.Relaxation(order=1),rho,matrix)
end

@testset "Continuous degree-three invariants and antiunitary actions" begin
    _,(x,y,z)=Nx.pauli_variables(3)
    ds=[Nx.SpinRotation(axis,(UInt64(7),)) for axis in 0:2]
    pair(i,j)=x[i]*x[j]+y[i]*y[j]+z[i]*z[j]
    chirality=x[1]*y[2]*z[3]+y[1]*z[2]*x[3]+z[1]*x[2]*y[3]-
        x[1]*z[2]*y[3]-z[1]*y[2]*x[3]-y[1]*x[2]*z[3]
    W=typeof(first(x[1].data.words))
    matrix(w)=physical_pauli(w,1:3)
    for chiral in (false,true)
        H=pair(1,2)+2pair(2,3)+3pair(1,3)+(chiral ? 0.7chirality : zero(chirality))
        anti=chiral ? () : (Nx.PauliConjugation(one(W),true),)
        problem=Nx.Problem(H;symmetry=Nx.Symmetry(moment_only=anti,continuous=ds))
        description=Nx.Relaxation(order=2)
        built=Nx.build(problem,description)
        rho=exp(-0.17sum(c*matrix(w) for (c,w) in Nx.terms(H))); rho/=LinearAlgebra.tr(rho)
        result=continuous_true_result(built,problem,description,rho,matrix)
        expected=real(LinearAlgebra.tr(rho*sum(c*matrix(w) for (c,w) in Nx.terms(chirality))))
        @test Nx.expectation(result,chirality)≈expected atol=1e-11
        chiral && @test abs(expected)>0.1
        @test Nx.expectation(result,x[1]*y[2]*z[3])≈expected/6 atol=1e-11
        Nx.refresh_coverage!(built.table,built.program)
        @test all(id->Nx.moment_covered(built.table,id),built.program.forms.ids)
    end
end
