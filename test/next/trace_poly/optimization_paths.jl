function trace_paths(solver,sos_solver;references=false)
    @testset "Dense trace CHSH sizes and forms" begin
        # One group: the trace formulation does not commute Alice with Bob.
        _,(x,)=site_variables(Unipotent(),1,4;functional=NCTSSoS.RealTrace())
        h=-trace_symbol(x[1]*x[3]+x[1]*x[4]+x[2]*x[3]-x[2]*x[4])
        built=NCTSSoS.build(Problem(h),Relaxation(order=1))
        data=TOML.parsefile(joinpath(@__DIR__,"../../data/expectations/trace_optimization_paths.toml"))
        oracle=only(filter(c->c["id"]=="dense_auto_order",data["cases"]))["expected"]
        @test NCTSSoS.block_sides(built)==oracle["sides"]
        @test built.sizes.moments==oracle["nuniq"]
        # This is also the Dense record in chsh_trace.toml.
        # The cold solve took 28.3 s; keep the cheap size checks routine.
        moment=references ? solve(built;solver) : nothing
        @test moment.bound ≈ oracle["objective"] atol=1e-5 skip=!references
        @test solve(built;solver=sos_solver).bound ≈ oracle["objective"] atol=1e-5 skip=!references || isnothing(sos_solver)
        # Public replacement for substituting a moment absent from the basis.
        if references
            @test_throws ArgumentError expectation(moment,trace_symbol(x[1]*x[2]*x[3]))
        else
            @test_skip expectation(moment,trace_symbol(x[1]*x[2]*x[3]))
        end
    end
    @testset "Trace equality agrees in moment and SOS forms" begin
        _,(x,)=site_variables(Unipotent(),1,2;functional=NCTSSoS.RealTrace())
        h=-trace_symbol(x[1]*x[2])
        problem=Problem(h;constraints=(Zero(trace_symbol(x[1])),))
        # This solve takes 30.0 s when it is the first cold solver case.
        built=references ? NCTSSoS.build(problem,Relaxation(order=1)) : nothing
        moment=references ? solve(built;solver) : nothing
        # The old file gives no independent value. A=B=diag(1,-1) has
        # normalized trace(A)=0 and objective -1; Cauchy–Schwarz gives >=-1.
        @test moment.bound ≈ -1 atol=1e-5 skip=!references
        @test solve(built;solver=sos_solver).bound ≈ moment.bound atol=1e-5 skip=!references || isnothing(sos_solver)
    end
end
