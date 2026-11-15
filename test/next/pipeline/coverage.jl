function check_used_coordinates(program;expected=0)
    covered=union(Set(program.forms.ids),Set(program.objective_ids),Set([1]))
    missing=setdiff(1:program.nvars,covered)
    @test length(missing)==expected
    missing
end

@testset "Discarded sector coordinates remain uncovered" begin
    _,(a,ad)=Nx.fermionic_variables(4)
    N=sum(ad[i]*a[i] for i in 1:4)
    H=-(ad[1]*a[2]+ad[2]*a[1]+ad[3]*a[4]+ad[4]*a[3])
    description=Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,1,4,3])]),strengthen=(Nx.Sector(N,2),))
    built=Nx.build(Nx.Problem(H),description)
    check_used_coordinates(built.program)
    @test built.sizes.closure.unused_coordinates==0
    @test occursin("0 unused coordinates",sprint(show,Nx.inspect(Nx.Problem(H),description)))
    @test length(only(collect(Nx.block_rows(built))))==9
    @test sort(Nx.block_sides(built))==[4,5]
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8))
    result=Nx.solve(Nx.Problem(H),description;solver)
    @test result.status==MOI.OPTIMAL
    @test result.bound≈-2 atol=1e-7
    before=(length(result.table.labels),length(result.table.base.cache),copy(result.table.covered))
    missing=Int[]
    for word in unique(label.word for label in result.table.base.labels)
        pairs=collect(Nx.moment_terms(result.table,word))
        unused=[id for (id,c) in pairs if !iszero(c) && !Nx.moment_covered(result.table,id)]
        isempty(unused) && continue
        append!(missing,unused)
        polynomial=prod(l<0 ? ad[-l] : a[l] for l in Nx.legacy_letters(word))
        @test_throws r"Moment .*selected blocks" Nx.expectation(result,polynomial)
    end
    @test !isempty(missing)
    @test Nx.expectation(result,H)≈result.bound atol=1e-7
    @test Nx.expectation(result,one(H))==1
    @test before==(length(result.table.labels),length(result.table.base.cache),result.table.covered)
    reconstruction=Nx.gns(result)
    @test reconstruction.verdict===:withheld
    @test isnothing(reconstruction.operators)

    # Manual objective edits require an explicit refresh, including removals.
    id=first(missing)
    built.program.nvars=max(built.program.nvars,id)
    push!(built.program.objective_ids,id)
    push!(built.program.objective_coefs,1.0)
    @test !Nx.moment_covered(built.table,id)
    Nx.refresh_coverage!(built.table,built.program)
    @test Nx.moment_covered(built.table,id)
    pop!(built.program.objective_ids); pop!(built.program.objective_coefs)
    Nx.refresh_coverage!(built.table,built.program)
    @test !Nx.moment_covered(built.table,id)
    Nx.push_term!(built.program.forms,1.0,id)
    Nx.close_form!(built.program.forms)
    Nx.push_block!(built.program,Nx.ZERO,1)
    @test !Nx.moment_covered(built.table,id)
    Nx.refresh_coverage!(built.table,built.program)
    @test Nx.moment_covered(built.table,id)
    @test Nx.moment_covered(built.table,1)
    @test !Nx.moment_covered(built.table,0)
end
