function check_used_coordinates(program;expected=0)
    covered=union(Set(program.forms.ids),Set(program.objective_ids),Set([1]))
    missing=setdiff(1:program.nvars,covered)
    @test length(missing)==expected
    missing
end

@testset "Cancelled sector coordinates are reported" begin
    _,(a,ad)=Nx.fermionic_variables(4)
    N=sum(ad[i]*a[i] for i in 1:4)
    H=-(ad[1]*a[2]+ad[2]*a[1]+ad[3]*a[4]+ad[4]*a[3])
    description=Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,1,4,3])]),strengthen=(Nx.Sector(N,2),))
    built=Nx.build(Nx.Problem(H),description)
    check_used_coordinates(built.program;expected=6)
    @test built.sizes.closure.unused_coordinates==6
    @test occursin("6 unused coordinates",sprint(show,Nx.inspect(Nx.Problem(H),description)))
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8))
    result=Nx.solve(Nx.Problem(H),description;solver)
    @test result.status==MOI.OPTIMAL
    @test result.bound≈-2 atol=1e-7
end
