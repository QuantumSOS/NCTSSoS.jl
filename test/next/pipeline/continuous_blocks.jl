function coefficient_program(matrices)
    p=Nx.MomentProgram{ComplexF64}(length(matrices))
    n=size(first(matrices),1)
    for j in 1:n,i in j:n
        for (id,A) in enumerate(matrices)
            iszero(A[i,j]) || Nx.push_term!(p.forms,A[i,j],id)
        end
        Nx.close_form!(p.forms)
    end
    Nx.push_block!(p,Nx.HPSD,n)
    Nx.validate(p)
end

function assigned_block(program,block,values)
    A=zeros(ComplexF64,block.side,block.side)
    for j in 1:block.side,i in j:block.side
        f=block.first_form+Nx.tri_offset(block.side,i,j)
        A[i,j]=sum((program.forms.coefs[t]*values[program.forms.ids[t]] for t in Nx.term_range(program.forms,f));init=0.0im)
        A[j,i]=conj(A[i,j])
    end
    Hermitian(A)
end

function check_split_spectra(original,reduced,report)
    rng=MersenneTwister(251)
    for _ in 1:4
        values=randn(rng,original.nvars); values[1]=1
        for record in report.blocks
            input=original.blocks[record.input]
            input.kind in (Nx.PSD,Nx.HPSD) || continue
            expected=eigvals(assigned_block(original,input,values))
            actual=Float64[]
            for (output,multiplicity) in zip(record.outputs,record.multiplicities)
                append!(actual,repeat(eigvals(assigned_block(reduced,reduced.blocks[output],values)),multiplicity))
            end
            @test sort!(actual)≈expected atol=1e-8 rtol=1e-8
        end
    end
end

function continuous_bound(program)
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-9,"eps_rel"=>1e-9))
    solution=Nx.solve_standard(Nx.moment_standard_form(program),solver)
    @test solution.status==MOI.OPTIMAL
    solution.bound
end

@testset "Coefficient certificates reject false copies" begin
    eye=Matrix{ComplexF64}(I,4,4)
    B=[ComplexF64[1 2im;-2im 3],ComplexF64[2 1+im;1-im -1]]
    matrices=[kron(Matrix{ComplexF64}(I,2,2),b) for b in B]
    copies=[[eye[:,1:2],eye[:,3:4]]]
    @test Nx._certify_coefficient_copies(matrices,copies).accepted
    perturbed=deepcopy(copies); perturbed[1][1][1,1]+=1e-3
    @test !Nx._certify_coefficient_copies(matrices,perturbed).accepted
    # The first matrix has two dimension-two eigenspaces. The second gives
    # an identity transport between them, but the third distinguishes copies.
    A0=Matrix{ComplexF64}(Diagonal([1,1,2,2]))
    A1=ComplexF64[0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0]
    A2=Matrix{ComplexF64}(Diagonal([1,-1,1,-1]))
    wrong=[[eye[:,[1,3]],eye[:,[2,4]]]]
    rejected=Nx._certify_coefficient_copies([A0,A1,A2],wrong)
    @test rejected.orthogonality_error==0
    @test !rejected.accepted
    @test rejected.reconstruction_error>0.1

    # Force the deterministic seed to be degenerate although neither
    # coefficient is scalar. Discovery may fail; the output must remain valid.
    a=ComplexF64[0.01 0.02im;-0.02im -0.01]
    accidental=coefficient_program([a,-sin(sqrt(2.0))/sin(2sqrt(2.0))*a])
    for family_blocks in (0,1)
        reduced,report=Nx._continuous_blocks(accidental;family_blocks)
        @test report.kept[:certificate_failed]==1
        @test reduced.forms.coefs==accidental.forms.coefs
        @test reduced.forms.ids==accidental.forms.ids
        check_split_spectra(accidental,reduced,report)
    end
    nested=coefficient_program([cat(A,zeros(ComplexF64,1,1);dims=(1,2))
        for A in (a,-sin(sqrt(2.0))/sin(2sqrt(2.0))*a)])
    reduced,report=Nx._continuous_blocks(nested;family_blocks=0)
    @test only(report.blocks).sides==[2,1]
    @test report.reconstruction_error>1e-9 # The failed component proposal was retained.
    check_split_spectra(nested,reduced,report)

    original=coefficient_program(matrices)
    reduced,report=Nx._continuous_blocks(original)
    @test report.split==1
    @test report.largest_after==2
    @test all(b->b.kind==Nx.HPSD,reduced.blocks)
    check_split_spectra(original,reduced,report)
    for kwargs in ((max_side=1,),(max_bytes=1,))
        kept,reason=Nx._continuous_blocks(original;kwargs...)
        # Exact family components precede the spectral resource caps too.
        @test only(reason.blocks).reason==:support_split
        @test only(reason.blocks).sides==[2,2]
        check_split_spectra(original,kept,reason)
        # Interleaved components test both row mapping and complex conjugation.
        perm=[1,3,2,4]
        disconnected=coefficient_program([A[perm,perm] for A in matrices])
        split,report=Nx._continuous_blocks(disconnected;family_blocks=0,kwargs...)
        @test only(report.blocks).reason==:support_split
        @test only(report.blocks).sides==[2,2]
        @test report.reconstruction_error==report.orthogonality_error==0
        check_split_spectra(disconnected,split,report)
        # Even a very small coupling must prevent an exact support split.
        coupled=deepcopy(matrices)
        coupled[1][1,3]=coupled[1][3,1]=1e-16
        tiny=coefficient_program(coupled)
        kept,report=Nx._continuous_blocks(tiny;family_blocks=0,kwargs...)
        @test report.split==0
        @test kept.forms.coefs==tiny.forms.coefs
        # Diagonal actions never invoke discovery, even on a connected cone.
        kept,report=Nx._continuous_blocks(tiny;diagonal=true,kwargs...)
        @test only(report.blocks).reason==:diagonal
        @test kept.forms.coefs==tiny.forms.coefs
    end
end

@testset "Continuous split: family spectra, field and determinism" begin
    for n in (4,6),finite in (false,true)
        registry,(x,y,z)=Nx.pauli_variables(n)
        H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
        generators=finite ? [Nx.SitePermutation([2:n;1]),Nx.SitePermutation([1;n:-1:2])] : Nx.AbstractAction[]
        ds=Nx.spin_rotations(registry)
        # Never build the six-site order-two dense comparison.
        description=Nx.Relaxation(order=n==4 ? 2 : 1)
        plain=Nx.build(Nx.Problem(H;symmetry=Nx.Symmetry(generators...)),Nx.BlockRelaxation(order=description.order))
        unreduced,_,_=Nx._continuous_program(plain.program,plain.table,ds)
        problem=Nx.Problem(H;symmetry=Nx.Symmetry(generators...;continuous=ds))
        built=Nx.build(problem,description)
        repeated=Nx.build(problem,description)
        @test repeated.program.forms.coefs==built.program.forms.coefs
        @test repeated.program.forms.ids==built.program.forms.ids
        @test repeated.program.forms.offsets==built.program.forms.offsets
        @test repeated.program.nvars==built.program.nvars
        @test repeated.program.objective_ids==built.program.objective_ids
        @test repeated.program.objective_coefs==built.program.objective_coefs
        @test repeated.program.blocks==built.program.blocks
        @test all(b->b.kind==Nx.HPSD,built.program.blocks)
        @test built.sizes.closure.continuous.blocks.split>0
        check_split_spectra(unreduced,built.program,built.sizes.closure.continuous.blocks)
        @test continuous_bound(unreduced)≈continuous_bound(built.program) atol=1e-6
    end

    registry,(a,ad)=Nx.fermionic_variables(3)
    H=ad[1]*a[2]+ad[2]*a[1]+0.3ad[3]*a[3]
    ds=(Nx.NumberRotation(registry;modes=1:2),Nx.NumberRotation(registry;modes=3:3))
    description=Nx.Relaxation(order=1,strengthen=Nx.FermionicConditions(1:3;conditions=(:P,:Q,:G,:T1,:T2prime)))
    plain=Nx.build(Nx.Problem(H),description)
    original,_,_=Nx._continuous_program(plain.program,plain.table,ds)
    built=Nx.build(Nx.Problem(H;symmetry=Nx.Symmetry(continuous=ds)),description)
    report=built.sizes.closure.continuous.blocks
    @test length(report.blocks)==length(original.blocks)
    @test any(r->r.input>1 && length(r.outputs)>1,report.blocks)
    # The +1 and -1 charge blocks remain distinct in the absence of a pairing map.
    @test sort(first(report.blocks).sides)==[1,1,1,2,2]
    @test first(report.blocks).multiplicities==ones(Int,5)
    @test first(report.blocks).reason==:support_split
    @test report.reconstruction_error==report.orthogonality_error==0
    start=length(first(report.blocks).outputs)
    @test built.sizes.closure.strengthening.first_block==start+1
    @test built.sizes.closure.strengthening.added_blocks==length(built.program.blocks)-start
    check_split_spectra(original,built.program,report)
    @test continuous_bound(original)≈continuous_bound(built.program) atol=1e-6
end

@testset "Continuous split with a moment-only antiunitary" begin
    registry,(x,y,z)=Nx.pauli_variables(4)
    H=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+z[i]*z[mod1(i+1,4)] for i in 1:4)
    W=typeof(first(H.data.words))
    anti=Nx.PauliConjugation(one(W),true)
    ds=Nx.spin_rotations(registry)
    plain=Nx.build(Nx.Problem(H;symmetry=Nx.Symmetry(moment_only=(anti,))),Nx.BlockRelaxation(order=1))
    original,_,_=Nx._continuous_program(plain.program,plain.table,ds)
    built=Nx.build(Nx.Problem(H;symmetry=Nx.Symmetry(moment_only=(anti,),continuous=ds)),Nx.Relaxation(order=1))
    @test built.sizes.closure.closed_blocks==plain.sizes.closure.closed_blocks
    @test all(b->b.kind==Nx.PSD,built.program.blocks)
    check_split_spectra(original,built.program,built.sizes.closure.continuous.blocks)
    @test continuous_bound(original)≈continuous_bound(built.program) atol=1e-6
end

@testset "Non-diagonal strengthening refines exact components" begin
    n=6
    registry,(x,y,z)=Nx.pauli_variables(n)
    H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    translation=Nx.SitePermutation([2:n;1])
    ds=Nx.spin_rotations(registry)
    description=Nx.Relaxation(order=1,strengthen=(Nx.DensityMatrix(1:3),Nx.InvariantCone(1:3,"YYY")))
    plain=Nx.build(Nx.Problem(H;symmetry=(translation,)),description)
    original,_,_=Nx._continuous_program(plain.program,plain.table,ds)
    family_blocks=plain.sizes.closure.strengthening.first_block-1
    # A side-one cap prevents any nontrivial spectral split in this control.
    exact,exact_report=Nx._continuous_blocks(original;family_blocks,diagonal=true,max_side=1)
    built=Nx.build(Nx.Problem(H;symmetry=Nx.Symmetry(translation;continuous=ds)),description)
    report=built.sizes.closure.continuous.blocks
    @test built.sizes.closure.continuous.method==:exact_elimination
    @test any(r.input>family_blocks && maximum(r.sides)<maximum(s.sides)
        for (r,s) in zip(report.blocks,exact_report.blocks))
    for (r,s) in zip(report.blocks,exact_report.blocks)
        r.input>family_blocks || continue
        @info "Strengthening cone sides" label=original.blocks[r.input].label unsplit=r.before exact=string(s.sides) spectral=string(r.sides)
    end
    check_split_spectra(original,exact,exact_report)
    check_split_spectra(original,built.program,report)
    @test continuous_bound(original)≈continuous_bound(built.program) atol=1e-6
end

@testset "Continuous grouped and localizing blocks" begin
    for localizer in (false,true)
        n=localizer ? 4 : 6
        registry,(x,y,z)=Nx.pauli_variables(n)
        bond(i,j)=x[i]*x[j]+y[i]*y[j]+z[i]*z[j]
        H=sum(bond(i,mod1(i+1,n)) for i in 1:n)
        ds=Nx.spin_rotations(registry)
        generators=localizer ? Nx.AbstractAction[] : [Nx.SitePermutation([2:n;1]),Nx.SitePermutation([1;n:-1:2])]
        constraints=localizer ? (Nx.Positive(3+bond(1,2)),) : ()
        problem=Nx.Problem(H;constraints,symmetry=Nx.Symmetry(generators...;continuous=ds))
        description=Nx.Relaxation(order=2,grouping=Nx.PhysicsInformed(level=0),split=Nx.TermSplit())
        built=Nx.build(problem,description)
        # Compare the very same closed family, not the pre-closure split.
        p,t,_=Nx.family_program(built.sizes.closure.family,[(ComplexF64(c),w) for (c,w) in Nx.terms(H)];generators)
        unreduced,_,_=Nx._continuous_program(p,t,ds)
        report=built.sizes.closure.continuous.blocks
        @test report.split>0
        check_split_spectra(unreduced,built.program,report)
        @test continuous_bound(unreduced)≈continuous_bound(built.program) atol=1e-6
        if localizer
            @test any(r->unreduced.blocks[r.input].label==1 && maximum(r.sides)<r.before,report.blocks)
        end
    end
end
