function chsh_simple(solver)
    @testset "CHSH" begin
        _, (a, b) = site_variables(Unipotent(), 2, 2)
        h = -(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
        for (id, relaxation) in (
            ("Dense_d1", Relaxation(order=1)),
            ("CS_d1", Relaxation(order=1, grouping=Correlative(algorithm=MF()))),
            ("TS_d1", Relaxation(order=1, split=TermSplit(algorithm=MMD(), level=1))),
        )
            @testset "$id" begin
                oracle = fixture("chsh_simple", id)
                built = Next.build(Problem(h), relaxation)
                result = solve(built; solver)
                @test result.bound ≈ oracle["objective"] atol=1e-6
                structure = fixture("next_problem_pilot", "chsh_simple/$id")
                @test sort(Next.block_sides(built)) == structure["sides"]
                @test result.sizes.moments == structure["moments"]
            end
        end
        symmetry = (
            LetterPermutation(a[1]=>a[2], a[2]=>a[1], b[2]=>-b[2]),
            LetterPermutation(a[2]=>-a[2], b[1]=>b[2], b[2]=>b[1]),
            LetterPermutation(a[1]=>b[1], a[2]=>b[2], b[1]=>a[1], b[2]=>a[2]),
        )
        built = Next.build(Problem(h; symmetry), Relaxation(order=1))
        result = solve(built; solver)
        @test result.bound ≈ fixture("chsh_simple", "Symmetry_d1")["objective"] atol=1e-6
        structure = fixture("next_problem_pilot", "chsh_simple/Symmetry_d1")
        @test sort(Next.block_sides(built)) == structure["sides"]
        @test result.sizes.real_variables == structure["real_variables"]
    end
end
