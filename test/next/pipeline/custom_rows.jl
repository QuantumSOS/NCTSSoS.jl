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
