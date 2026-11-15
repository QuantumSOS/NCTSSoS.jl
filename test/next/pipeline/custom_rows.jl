@testset "Custom starting rows with term splitting" begin
    registry,(x,y,z)=Nx.pauli_variables(2)
    rows=[one(x[1]),x[1],x[2]]
    h=x[1]+2x[2]+0.3x[1]*x[2]
    for symmetry in ((),(Nx.SitePermutation([2,1]),))
        objective=isempty(symmetry) ? h : x[1]+x[2]+0.3x[1]*x[2]
        problem=Nx.Problem(objective;symmetry)
        plain=Nx.build(problem,Nx.Relaxation(;rows))
        split=Nx.build(problem,Nx.Relaxation(;rows,split=Nx.TermSplit(algorithm=MMD(),level=3)))
        @test solve_reference(split.program)[1]≈solve_reference(plain.program)[1] atol=1e-6
        @test collect(Nx.block_rows(split))==collect(Nx.block_rows(plain))
    end
    for options in ((grouping=Nx.Correlative(),),(row_rule=Nx.Contiguous(),),
        (split=Nx.PhysicsInformed(),))
        @test_throws r"rows with split=TermSplit" Nx.Relaxation(;rows,options...)
    end
    @test_throws r"not covered by the starting rows" Nx.build(Nx.Problem(z[2]),
        Nx.Relaxation(;rows,split=Nx.TermSplit()))
    # Symmetry closes algorithm families, but leaves explicit families supplied.
    problem=Nx.Problem(x[1]+x[2];symmetry=(Nx.SitePermutation([2,1]),))
    asymmetric=[one(x[1]),x[1],x[2],z[1]]
    built=Nx.build(problem,Nx.Relaxation(rows=asymmetric,split=Nx.TermSplit()))
    @test any(block->any(w->w==only(Nx.terms(z[2]))[2],block),Nx.block_rows(built))
    @test_throws ArgumentError Nx.build(problem,Nx.Relaxation(family=[asymmetric]))

    # A localizer uses all registry letters, not only letters in custom rows.
    _,(a,)=Nx.site_variables(Nx.Free(),1,2)
    problem=Nx.Problem(a[1]^4;constraints=[Nx.Positive(1-a[2]^2)])
    custom=[one(a[1]),a[1],a[1]^2]
    built=Nx.build(problem,Nx.Relaxation(rows=custom,split=Nx.TermSplit(level=3)))
    localizers=filter(b->b.source==1,built.sizes.closure.family)
    second_letter=Nx.legacy_letters(only(Nx.terms(a[2]))[2])
    @test any(b->second_letter in Nx.legacy_letters.(b.rows),localizers)
    @test solve_reference(built.program)[1]≈0 atol=1e-6
    @test_throws r"exceeds the row order" Nx.build(Nx.Problem(a[1]^2;
        constraints=[Nx.Positive(1-a[2]^6)]),Nx.Relaxation(rows=custom,split=Nx.TermSplit()))
    chip_problem=Nx.Problem((a[1]*a[2])'*(a[1]*a[2])+a[2]^2)
    chips=Nx.newton_chip_basis(chip_problem)
    split=Nx.build(chip_problem,Nx.Relaxation(rows=chips,split=Nx.TermSplit(level=3)))
    plain=Nx.build(chip_problem,Nx.Relaxation(rows=chips))
    @test solve_reference(split.program)[1]≈solve_reference(plain.program)[1] atol=1e-6
end

@testset "Close starting rows before term splitting" begin
    registry,(x,y,z)=Nx.pauli_variables(4)
    h=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+0.6z[i]*z[mod1(i+1,4)] for i in 1:4)
    base=[one(h);x;y;z]
    short=[base;x[1]*x[2];y[1]*y[2];z[1]*z[2]]
    complete=[base;[a[i]*a[mod1(i+1,4)] for a in (x,y,z) for i in 1:4]]
    problem=Nx.Problem(h;symmetry=(Nx.SitePermutation([2,3,4,1]),))
    for level in 1:2
        split=Nx.TermSplit(;algorithm=MMD(),level)
        automatic=Nx.build(problem,Nx.Relaxation(rows=short,split=split))
        manual=Nx.build(problem,Nx.Relaxation(rows=complete,split=split))
        @test collect(Nx.block_rows(automatic))==collect(Nx.block_rows(manual))
        @test solve_reference(automatic.program)[1]≈solve_reference(manual.program)[1] atol=1e-6
    end
    # Closure may supply objective support; moment-only maps must not add rows.
    problem=Nx.Problem(x[1]+x[2]+x[3]+x[4];symmetry=(Nx.SitePermutation([2,3,4,1]),))
    automatic=Nx.build(problem,Nx.Relaxation(rows=[one(h),x[1]],split=Nx.TermSplit()))
    manual=Nx.build(problem,Nx.Relaxation(rows=[one(h);x],split=Nx.TermSplit()))
    @test collect(Nx.block_rows(automatic))==collect(Nx.block_rows(manual))
    moment_only=Nx.Problem(problem.objective;
        symmetry=Nx.Symmetry(;moment_only=(Nx.SitePermutation([2,3,4,1]),)))
    @test_throws r"not covered by the starting rows" Nx.build(moment_only,
        Nx.Relaxation(rows=[one(h),x[1]],split=Nx.TermSplit()))
    problem=Nx.Problem(z[1];symmetry=(Nx.SpinRotation(registry,:z),))
    automatic=Nx.build(problem,Nx.Relaxation(rows=[one(h),x[1]],split=Nx.TermSplit()))
    manual=Nx.build(problem,Nx.Relaxation(rows=[one(h),x[1],y[1]],split=Nx.TermSplit()))
    @test collect(Nx.block_rows(automatic))==collect(Nx.block_rows(manual))
    @test solve_reference(automatic.program)[1]≈-1 atol=1e-6

    # Degree-changing closure must also determine the localizer row order.
    r,(a,b,c)=Nx.pauli_variables(2)
    W=Nx.word_type(r,4)
    h=a[1]+a[1]*a[2]+c[2]+c[1]*c[2]
    problem=Nx.Problem(h;constraints=[Nx.Positive(1+c[1])],
        symmetry=(Nx.clifford_cnot(W,2,1,2),))
    short=[one(h),a[1],c[1],a[2],c[2]]
    complete=[short;a[1]*a[2];c[1]*c[2]]
    automatic=Nx.build(problem,Nx.Relaxation(rows=short,split=Nx.TermSplit()))
    manual=Nx.build(problem,Nx.Relaxation(rows=complete,split=Nx.TermSplit()))
    @test collect(Nx.block_rows(automatic))==collect(Nx.block_rows(manual))
    @test solve_reference(automatic.program)[1]≈solve_reference(manual.program)[1] atol=1e-6
end
