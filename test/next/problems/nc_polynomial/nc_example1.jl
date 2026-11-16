function nc_example1(solver, sos_solver)
    @testset "NC example 1" begin
        _, (x,) = site_variables(Free(), 1, 3)
        f = x[1]^2 - x[1]*x[2] - x[2]*x[1] + 3x[2]^2 -
            2x[1]*x[2]*x[1] + 2x[1]*x[2]^2*x[1] -
            x[2]*x[3] - x[3]*x[2] + 6x[3]^2 +
            9x[2]^2*x[3] + 9x[3]*x[2]^2 -
            54x[3]*x[2]*x[3] + 142x[3]*x[2]^2*x[3]
        problem = Problem(f)
        for (id, relaxation) in (
            ("Dense_d2", Relaxation(order=2)),
            ("TS_d2", Relaxation(order=2, split=TermSplit(algorithm=MMD(), level=1))),
        )
            @testset "$id" begin
                oracle = fixture("nc_example1", id)
                built = NCTSSoS.build(problem, relaxation)
                result = solve(built; solver)
                @test result.bound ≈ oracle["objective"] atol=1e-4
                structure = fixture("next_problem_pilot", "nc_example1/$id")
                @test sort(NCTSSoS.block_sides(built)) == structure["sides"]
                @test result.sizes.moments == structure["moments"]
            end
        end
        @test solve(problem, Relaxation(order=2); solver=sos_solver).bound ≈
            fixture("nc_example1", "SOS_d2")["objective"] atol=1e-6 skip=isnothing(sos_solver)
    end
end
