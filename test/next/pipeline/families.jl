# Independent reference: select principal submatrices of the model's dense
# moment/localizing matrices. Only tests construct these full matrices.
function family_reference(problem,order,family)
    dense=Nx.build(problem,Nx.Relaxation(order=order))
    p=dense.program
    W=eltype(dense.table.keys)
    registry=problem.objective.registry
    domains=[Nx.basis_rows(W,registry.algebra,1:registry.nsites,order)]
    append!(domains,[Nx.basis_rows(W,registry.algebra,1:registry.nsites,order-cld(Nx.degree(c.polynomial),2)) for c in problem.constraints])
    result=Nx.MomentProgram{ComplexF64}(p.nvars)
    append!(result.objective_ids,p.objective_ids); append!(result.objective_coefs,p.objective_coefs)
    for block in family
        source=block.source+1
        indices=Dict(w=>i for (i,w) in enumerate(domains[source]))
        n=length(block.rows)
        for j in 1:n,i in j:n
            a,b=indices[block.rows[i]],indices[block.rows[j]]
            f=p.blocks[source].first_form+Nx.tri_offset(length(domains[source]),max(a,b),min(a,b))
            for t in Nx.term_range(p.forms,f)
                Nx.push_term!(result.forms,a<b ? conj(p.forms.coefs[t]) : p.forms.coefs[t],p.forms.ids[t])
            end
            Nx.close_form!(result.forms)
        end
        Nx.push_block!(result,block.kind,block.kind==Nx.ZERO ? n*(n+1)÷2 : n,block.source)
    end
    result,dense.table
end

@testset "Block families, closure and one-call model" begin
    n=4; registry,(x,y,z)=Nx.pauli_variables(n)
    H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    problem=Nx.Problem(H)
    W=Nx.word_type(registry,2)
    translation=Nx.SitePermutation([2,3,4,1])
    reflection=Nx.SitePermutation([4,3,2,1])
    rotation=Nx.CliffordMap([Nx.pauli_letter(W,i,a) for a in (1,0) for i in 1:n])
    generators=[Nx._group_action(W,g) for g in (translation,reflection,rotation)]
    for grouping in (nothing,Nx.Correlative(),Nx.Correlative(algorithm=MMD()),Nx.PhysicsInformed(level=1)),
        split in (nothing,Nx.TermSplit(),Nx.PhysicsInformed(level=2))
        description=Nx.BlockRelaxation(order=1,grouping=grouping,split=split,symmetry=Nx.FamilySymmetry(generators))
        built=Nx.build(problem,description)
        closure=built.sizes.closure
        reference,reference_table=family_reference(problem,1,closure.family)
        @test solve_reference(built.program)[1]≈solve_reference(reference)[1] atol=1e-6
        coordinates=[Nx.orbit_coordinates!(built.table,label.word) for label in reference_table.labels]
        values=randn(MersenneTwister(182),length(built.table.labels)); values[1]=1
        old_values=[label.part==2 ? imag(coordinate_value(c,values)) : real(coordinate_value(c,values))
            for (label,c) in zip(reference_table.labels,coordinates)]
        @test minimum(eigmin(Hermitian(A)) for A in evaluated_blocks(built.program,values))≈
            minimum(eigmin(Hermitian(A)) for A in evaluated_blocks(reference,old_values)) atol=1e-9
        @test closure.closed_blocks>=closure.initial_blocks
        @test length(closure.block_orbits)<=closure.closed_blocks
        @test closure.block_generators==generators
        @test isempty(closure.moment_only)
    end
    description=Nx.BlockRelaxation(order=1,grouping=Nx.Correlative(),symmetry=Nx.FamilySymmetry(generators;
        moment_only=(Nx.AntiunitaryPauli(one(W)),)))
    # With only cyclic axis rotations, conjugation exchanges projected
    # sectors. The generic route must report that limit, not double a block.
    @test_throws r"No same-side real basis" Nx.build(problem,description)
    quarter=Nx.CliffordMap([Nx.pauli_letter(W,i,a) for a in (1,2) for i in 1:n])
    description=Nx.BlockRelaxation(order=1,grouping=Nx.Correlative(),symmetry=Nx.FamilySymmetry([generators;quarter];
        moment_only=(Nx.AntiunitaryPauli(one(W)),)))
    built=Nx.build(problem,description)
    @test all(b.kind==Nx.PSD for b in built.program.blocks)
    result=Nx.solve(problem,description;solver=Nx.JuMPBackend(SOLVER))
    @test result.bound≈solve_reference(family_reference(problem,1,built.sizes.closure.family)[1])[1] atol=1e-6
    @test !isempty(result.sizes.closure.moment_only)
    @test occursin("moment-only",sprint(show,Nx.inspect(problem,description)))
    @test isempty(result.valid_for.constraints)

    custom=[[one(x[1]),x[1],x[2]]]
    simple=Nx.Problem(sum(x[i]*x[mod1(i+1,n)] for i in 1:n))
    @test_throws r"not closed" Nx.build(simple,Nx.BlockRelaxation(family=custom,symmetry=Nx.FamilySymmetry([translation])))
    closed=[[one(x[1]),x[i],x[mod1(i+1,n)]] for i in 1:n]
    custom_build=Nx.build(simple,Nx.BlockRelaxation(family=closed,symmetry=Nx.FamilySymmetry([translation])))
    @test custom_build.sizes.closure.added_blocks==0
    @test length(custom_build.sizes.closure.family)==length(closed)
    @test all(length(b.rows)==3 for b in custom_build.sizes.closure.family)
    @test solve_reference(custom_build.program)[1]≈-4 atol=1e-6
    shifted_anti=Nx.compose(Nx._group_action(W,translation),Nx.AntiunitaryPauli(one(W)))
    real_build=Nx.build(simple,Nx.BlockRelaxation(family=closed,
        symmetry=Nx.FamilySymmetry([translation];moment_only=(shifted_anti,))))
    @test all(b.kind==Nx.PSD for b in real_build.program.blocks)
    @test solve_reference(real_build.program)[1]≈-4 atol=1e-6
    @test real_build.sizes.closure.added_blocks==0
    @test_throws r"does not cover" Nx.build(simple,Nx.BlockRelaxation(family=[[one(x[1]),x[1]]]))
    @test_throws ArgumentError Nx.BlockRelaxation(family=closed,split=Nx.TermSplit())
end

@testset "Fermionic families and localizing matrices" begin
    registry,(a,ad)=Nx.fermionic_variables(3)
    H=sum(cis(0.37)*ad[i]*a[mod1(i+1,3)]+cis(-0.37)*ad[mod1(i+1,3)]*a[i] for i in 1:3)
    problem=Nx.Problem(H)
    description=Nx.BlockRelaxation(order=1,grouping=Nx.PhysicsInformed(level=1),split=Nx.TermSplit(),
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,3,1])]))
    built=Nx.build(problem,description)
    reference,_=family_reference(problem,1,built.sizes.closure.family)
    @test solve_reference(built.program)[1]≈solve_reference(reference)[1] atol=1e-6

    registry,(x,y,z)=Nx.pauli_variables(3)
    H=sum(x)
    problem=Nx.Problem(H;constraints=(Nx.Positive(one(H)+0.2sum(z)),Nx.Zero(sum(y))))
    description=Nx.BlockRelaxation(order=2,split=Nx.TermSplit(),
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,3,1])]))
    built=Nx.build(problem,description)
    reference,_=family_reference(problem,2,built.sizes.closure.family)
    @test solve_reference(built.program)[1]≈solve_reference(reference)[1] atol=1e-6
    @test built.valid_for.constraints==problem.constraints
    broken=Nx.Problem(H;constraints=(Nx.Positive(one(H)+z[1]),))
    @test_throws r"constraint 1" Nx.build(broken,description)
end

@testset "Moment-only maps change the comparison reference" begin
    registry,(x,y,z)=Nx.pauli_variables(2)
    W=Nx.word_type(registry,2)
    cnot=Nx.clifford_cnot(W,2,1,2)
    H=x[1]+x[1]*x[2]+z[2]+z[1]*z[2]+0.5x[2]+0.3z[1]
    problem=Nx.Problem(H)
    rows=[[one(H),x[1],z[1],x[2],z[2]]]
    built=Nx.build(problem,Nx.BlockRelaxation(family=rows,symmetry=Nx.FamilySymmetry(;moment_only=(cnot,))))
    reference,table=family_reference(problem,1,built.sizes.closure.family)
    plain=solve_reference(reference)[1]
    for w in Nx.basis_rows(W,Nx.Pauli(),1:2,2)
        phase,image=Nx.act(cnot,w)
        iw=Nx.moment_ids!(table,w)[1]; iv=Nx.moment_ids!(table,image)[1]
        Nx.push_term!(reference.forms,1.,iw); Nx.push_term!(reference.forms,-Nx.scale(phase,1.),iv)
        Nx.close_form!(reference.forms); Nx.push_block!(reference,Nx.ZERO,1)
    end
    reference.nvars=Nx.nvars(table)
    bound=solve_reference(built.program)[1]
    @test bound≈solve_reference(reference)[1] atol=1e-6
    @test bound>plain+0.01
    @test built.sizes.closure.added_blocks==0
    @test isempty(built.sizes.closure.block_generators)
    @test length(built.sizes.closure.moment_only)==1
end

@testset "Selection algorithms, row rules and levels" begin
    registry,(x,y,z)=Nx.pauli_variables(4)
    H=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+z[i]*z[mod1(i+1,4)] for i in 1:4)
    problem=Nx.Problem(H)
    old_registry,old_ops=create_pauli_variables(1:4)
    old_H=heisenberg_chain_hamiltonian(old_ops;periodic=true,coupling=1)
    old_problem=polyopt(old_H,old_registry)
    for algorithm in (MF(),MMD())
        old=NCTSSoS.correlative_sparsity(old_problem,1,algorithm)
        # The old graph has letter vertices; the new graph has site
        # vertices. Their tie-dependent triangulations agree after closure.
        groups=[sort!(unique([(Int(i)-1)÷3+1 for i in clique])) for clique in old.cliques]
        expected=Set(Tuple(sort!([mod1(i+shift,4) for i in group])) for group in groups for shift in 0:3)
        built=Nx.build(problem,Nx.BlockRelaxation(order=1,grouping=Nx.Correlative(;algorithm),
            symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,3,4,1])])))
        actual=Set(Tuple(sort!(unique!(reduce(vcat,[collect(Nx.sites(w)) for w in b.rows];init=Int[])))) for b in built.sizes.closure.family)
        @test actual==expected
    end
    old_basis=get_ncbasis(old_registry,1)
    active=NCTSSoS.init_activated_supp(old_H,typeof(old_H)[],old_basis)
    old_split=NCTSSoS.term_sparsities(active,typeof(old_H)[],old_basis,typeof(old_basis)[],MaximalElimination())
    built=Nx.build(problem,Nx.BlockRelaxation(order=1,split=Nx.TermSplit()))
    expected=Set(Set(Tuple(w.word) for w in block) for block in first(old_split).block_bases)
    actual=Set(Set(Tuple(Nx.legacy_letters(w)) for w in b.rows) for b in built.sizes.closure.family)
    @test actual==expected

    for role in (:grouping,:split)
        previous=nothing
        for level in 0:4
            method=Nx.PhysicsInformed(;level)
            description=role==:grouping ? Nx.BlockRelaxation(order=1,grouping=method) : Nx.BlockRelaxation(order=1,split=method)
            family=Nx.build(problem,description).sizes.closure.family
            if !isnothing(previous)
                @test all(any(issubset(old.rows,new.rows) for new in family) for old in previous)
            end
            previous=family
        end
        @test length(previous)==1
        @test length(only(previous).rows)==13
    end
    @test Nx.next_level(Nx.TermSplit()).level==2
    @test Nx.next_level(Nx.PhysicsInformed()).level==2
    @test_throws ArgumentError Nx.TermSplit(level=0)
    @test_throws ArgumentError Nx.PhysicsInformed(level=-1)
    W=Nx.word_type(registry,4)
    custom_rule=(W,r,group,order)->Nx.basis_rows(W,r.algebra,group,order)
    full=Nx.build(problem,Nx.BlockRelaxation(order=2,grouping=Nx.Correlative()))
    custom=Nx.build(problem,Nx.BlockRelaxation(order=2,grouping=Nx.Correlative(),row_rule=custom_rule))
    @test [b.rows for b in custom.sizes.closure.family]==[b.rows for b in full.sizes.closure.family]
    contiguous=Nx.build(problem,Nx.BlockRelaxation(order=2,row_rule=Nx.Contiguous(periodic=true)))
    @test length(only(contiguous.sizes.closure.family).rows)==49 # 1 + 4*3 + 4*9
    @test solve_reference(contiguous.program)[1]≈solve_reference(family_reference(problem,2,contiguous.sizes.closure.family)[1])[1] atol=1e-6
end

@testset "Automatic real block families and complex fallback" begin
    function hermitian_reference(problem,built;generators=Nx.AbstractAction[])
        W=eltype(built.table.representatives.keys)
        objective=[(c,Nx.word_from_legacy(W,Nx.legacy_letters(w))) for (c,w) in Nx.terms(problem.objective)]
        first(Nx.family_program(built.sizes.closure.family,objective;generators))
    end
    _,(a,b)=Nx.site_variables(Nx.Unipotent(),2,2)
    problem=Nx.Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2]))
    built=Nx.build(problem,Nx.BlockRelaxation(order=1))
    @test only(built.program.blocks).kind==Nx.PSD
    @test only(built.program.blocks).side==5
    @test built.program.nvars==11 # Includes normalization; Hermitian has 13.
    @test eltype(built.program.forms.coefs)==Float64
    reference=hermitian_reference(problem,built)
    @test reference.nvars==13
    bound=solve_reference(built.program)[1]
    @test bound≈-2sqrt(2) atol=1e-6
    @test bound≈solve_reference(reference)[1] atol=1e-6

    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        _,ops=Nx.site_variables(rule,1,2)
        a,b=ops[1]
        real_build=Nx.build(Nx.Problem(a*b+b*a),Nx.BlockRelaxation(order=1))
        @test only(real_build.program.blocks).kind==Nx.PSD
        @test all(label.part!=2 for label in real_build.table.labels)
        complex_problem=Nx.Problem(im*(a*b-b*a))
        complex_build=Nx.build(complex_problem,Nx.BlockRelaxation(order=1))
        @test only(complex_build.program.blocks).kind==Nx.HPSD
        @test any(label.part==2 for label in complex_build.table.labels)
        rule isa Nx.Unipotent && @test solve_reference(complex_build.program)[1]≈-2 atol=1e-6
        for constraint in (Nx.Positive(1+0.25im*(a*b-b*a)),Nx.Zero(im*(a*b-b*a)))
            constrained=Nx.build(Nx.Problem(a+b;constraints=(constraint,)),Nx.BlockRelaxation(order=1))
            @test first(constrained.program.blocks).kind==Nx.HPSD
        end
    end
    _,ops=Nx.site_variables(Nx.Free(),1,2)
    a,b=ops[1]
    problem=Nx.Problem(a*b+b*a;constraints=(Nx.Positive(1-a^2-b^2),))
    built=Nx.build(problem,Nx.BlockRelaxation(order=2))
    @test all(block.kind==Nx.PSD for block in built.program.blocks)
    bound=solve_reference(built.program)[1]
    @test bound≈-1 atol=1e-6
    @test bound≈solve_reference(hermitian_reference(problem,built))[1] atol=1e-6

    # Real finite-group projections must retain the conjugate momentum
    # pair of a three-cycle, rather than discard imaginary Fourier entries.
    _,ops=Nx.site_variables(Nx.Unipotent(),3,2)
    problem=Nx.Problem(sum(ops[i][1]*ops[mod1(i+1,3)][1]+0.3ops[i][2] for i in 1:3))
    for generators in ([Nx.SitePermutation([2,3,1])],
        [Nx.SitePermutation([2,3,1]),Nx.SitePermutation([1,3,2])]),split in (nothing,Nx.TermSplit())
        built=Nx.build(problem,Nx.BlockRelaxation(;order=1,split,symmetry=Nx.FamilySymmetry(generators)))
        @test all(block.kind==Nx.PSD for block in built.program.blocks)
        @test solve_reference(built.program)[1]≈solve_reference(hermitian_reference(problem,built;generators))[1] atol=1e-6
    end
    _,(_,_,z)=Nx.pauli_variables(1)
    @test only(Nx.build(Nx.Problem(z[1]),Nx.BlockRelaxation(order=1)).program.blocks).kind==Nx.HPSD
end
