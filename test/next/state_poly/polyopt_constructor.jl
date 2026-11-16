function state_constructor()
    @testset "State polynomial expansions" begin
        _,(x,y) = site_variables(Unipotent(),2,2;functional=NCTSSoS.RealState())
        a,b,c = state_symbol(x[1]*y[2]),state_symbol(x[1]*y[1]),state_symbol(x[2]*y[2])
        # Survey 7.2.1: cross-site commutation makes its first sum 2a.
        p = Problem((state_symbol(x[1]*y[2])+state_symbol(y[2]*x[1]))^2+(b-c)^2)
        @test p.objective == 4a^2+b^2-2b*c+c^2
        @test isempty(p.constraints)
    end
end
