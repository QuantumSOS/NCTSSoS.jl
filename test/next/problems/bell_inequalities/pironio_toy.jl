function pironio_toy(solver, sos_solver)
    @testset "Pironio E11" begin
        _, (x,) = site_variables(Free(), 1, 2)
        problem = Problem(x[1]*x[2]+x[2]*x[1];
            constraints=(Zero(x[1]^2-x[1]), Positive(-x[2]^2+x[2]+0.5)))
        built = Next.build(problem, Relaxation(order=2))
        oracle = fixture("pironio_toy", "Dense_d2")
        result = solve(built; solver)
        @test result.bound ≈ oracle["objective"] atol=1e-6
        @test solve(built; solver=sos_solver).bound ≈ oracle["objective"] atol=1e-6 skip=isnothing(sos_solver)
        structure = fixture("next_problem_pilot", "pironio_toy/Dense_d2")
        @test sort(Next.block_sides(built)) == structure["sides"]
        @test result.sizes.moments == structure["moments"]
    end
end
