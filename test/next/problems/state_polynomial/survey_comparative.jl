function survey_witnesses()
    @testset "Survey independent projector witnesses" begin
        x = [1.0 0;0 0]
        y = [0.5 0.5;0.5 0.5]
        v = [sqrt(2+sqrt(2))/2,sqrt(2-sqrt(2))/2]
        state(m) = real(dot(v,m*v))
        value = (state(x*y*x*y)+state(y*x*y*x))/2-state(x*y*x)*state(y)
        @test norm(x*x-x)<=1e-12 && norm(y*y-y)<=1e-12
        @test norm(v) ≈ 1 atol=1e-12
        @test value ≈ -1/16 atol=1e-12
        x = Matrix(Diagonal([0.,1.,0.]))
        y = [1. 0 0;0 1/3 sqrt(2)/3;0 sqrt(2)/3 2/3]
        trace(m) = real(tr(m))/3
        value = (trace(x*y*x*y)+trace(y*x*y*x))/2-trace(x*y*x)*trace(y)
        @test norm(x*x-x)<=1e-12 && norm(y*y-y)<=1e-12
        @test value ≈ -1/27 atol=1e-12
    end
end

function survey_comparative(solver,sos_solver)
    data = TOML.parsefile(joinpath(@__DIR__,"../../../data/expectations/state_survey_comparative.toml"))
    @testset "Survey dense state and trace bounds" begin
        for (functional,symbol,id) in ((Next.RealState(),state_symbol,"S6_state_dense_order3"),
                (Next.RealTrace(),trace_symbol,"S7_trace_dense_order3"))
            _,(p,) = site_variables(Projector(),1,2;functional)
            x,y = p
            @test x*y != y*x
            h = 0.5symbol(x*y*x*y)+0.5symbol(y*x*y*x)-symbol(x*y*x)*symbol(y)
            oracle = only(filter(c->c["id"]==id,data["cases"]))["expected"]
            built = Next.build(Problem(h),Relaxation(order=3))
            @test Next.block_sides(built) == oracle["sides"]
            @test built.sizes.moments == oracle["nuniq"]
            if functional isa Next.RealState
                # COSMO reaches its iteration limit with the old S6 options.
                # A raw objective is not a bound. Keep the witness and build
                # checks in CI; run the bound only with the licensed backend.
                @test solve(built;solver=sos_solver).bound ≈ oracle["objective"] atol=1e-5 skip=isnothing(sos_solver)
            else
                @test solve(built;solver).bound ≈ oracle["objective"] atol=1e-5
            end
        end
    end
end
