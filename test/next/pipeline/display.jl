@testset "Pipeline option displays" begin
    examples=(
        (Nx.Correlative(),"Correlative()"),
        (Nx.Correlative(algorithm=MMD()),"Correlative(algorithm=MMD())"),
        (Nx.Correlative(algorithm=MMD(delta=2)),"Correlative(algorithm=MMD(delta=2))"),
        (Nx.TermSplit(),"TermSplit()"),
        (Nx.TermSplit(level=2),"TermSplit(level=2)"),
        (Nx.TermSplit(algorithm=MF()),"TermSplit(algorithm=MF())"),
        (Nx.TermSplit(algorithm=MMD(),level=3),"TermSplit(algorithm=MMD(), level=3)"),
        (Nx.PhysicsInformed(),"PhysicsInformed()"),
        (Nx.PhysicsInformed(level=0),"PhysicsInformed(level=0)"),
        (Nx.Contiguous(),"Contiguous()"),
        (Nx.Contiguous(periodic=true),"Contiguous(periodic=true)"),
        (Nx.DensityMatrix(1:4),"DensityMatrix(1:4)"),
        (Nx.DensityMatrix([2,4]),"DensityMatrix([2, 4])"),
        (Nx.DensityMatrix([3]),"DensityMatrix(3:3)"),
        (Nx.InvariantCone(1:3,"YXX"),"InvariantCone(1:3, \"YXX\")"),
        (Nx.Stationarity(),"Stationarity()"),
        (Nx.StateOptimality(),"StateOptimality()"),
        (Nx.FermionicConditions(1:4),"FermionicConditions(1:4; conditions=(:P, :Q, :G))"),
        (Nx.FermionicConditions([1,3];conditions=(:T2prime,)),"FermionicConditions([1, 3]; conditions=(:T2prime,))"),
    )
    for (option,expected) in examples
        @test sprint(show,option)==expected
        @test sprint(show,MIME"text/plain"(),option)==expected
        @test sprint(show,option;context=:compact=>true)==expected
    end
    _,(_,_,z)=Nx.pauli_variables(2)
    sector=Nx.Sector(z[1]+z[2],0)
    @test sprint(show,sector)=="Sector($(sector.operator), 0)"
    @test sprint(show,MIME"text/plain"(),sector)==sprint(show,sector)

    # The public front door and result delegate option printing to this layer.
    problem=Nx.Problem(z[1]+z[2])
    for relaxation in (Nx.Relaxation,Nx.BlockRelaxation)
        description=relaxation(order=1,row_rule=Nx.Contiguous(periodic=true),
            grouping=Nx.Correlative(),split=Nx.TermSplit(level=2),
            strengthen=(Nx.DensityMatrix(1:2),Nx.Stationarity()))
        result=Nx.SolveResult(nothing,:not_solved,nothing,nothing,problem,description,nothing,nothing,nothing)
        for text in (sprint(show,description),sprint(show,result))
            @test occursin("row_rule=Contiguous(periodic=true)",text)
            @test occursin("grouping=Correlative()",text)
            @test occursin("split=TermSplit(level=2)",text)
            @test occursin("strengthen=(DensityMatrix(1:2), Stationarity())",text)
            @test !occursin("CliqueTrees",text)
            @test !occursin("NCTSSoS.",text)
        end
    end
end
