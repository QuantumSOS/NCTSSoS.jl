function nc_example2(solver, sos_solver)
    @testset "NC example 2" begin
        _, (x,) = site_variables(Free(), 1, 2)
        f = 2 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2
        problem = Problem(f; constraints=(Positive(4-x[1]^2-x[2]^2),
            Zero(x[1]*x[2]+x[2]*x[1]-2)))
        opt = fixture("nc_example2", "opt")["objective"]
        @test solve(problem, Relaxation(order=2); solver).bound ≈ opt atol=1e-6
        @test solve(problem, Relaxation(order=2); solver=sos_solver).bound ≈ opt atol=1e-6 skip=isnothing(sos_solver)
        result = solve(problem, Relaxation(order=2, grouping=Correlative(algorithm=MF()),
            split=TermSplit(algorithm=MMD(), level=1)); solver)
        @test result.bound ≈ opt atol=1e-6
        @test solve(problem, Relaxation(order=2, split=TermSplit(algorithm=MMD(), level=1));
            solver=sos_solver).bound ≈ opt atol=1e-6 skip=isnothing(sos_solver)
        higher = solve(problem, Relaxation(order=2, grouping=Correlative(algorithm=MF()),
            split=TermSplit(algorithm=MMD(), level=2)); solver)
        # Test level 2 against the independent optimum, not another package result.
        @test higher.bound ≈ opt atol=1e-4
    end
end
