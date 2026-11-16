function harmonic_oscillator(solver)
    @testset "Harmonic oscillator" begin
        _, (b, bd) = bosonic_variables(1)
        built = Next.build(Problem(bd[1]*b[1]+0.5), Relaxation(order=2))
        result = solve(built; solver)
        @test result.bound ≈ 0.5 atol=1e-6
        @test result.bound ≈ fixture("harmonic_oscillator", "single_mode_order2")["objective"] atol=1e-6
        structure = fixture("next_problem_pilot", "harmonic_oscillator/single_mode_order2")
        @test sort(Next.block_sides(built)) == structure["sides"]
        @test result.sizes.moments == structure["moments"]
    end
end
