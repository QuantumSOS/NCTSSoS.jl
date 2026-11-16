function bose_hubbard(solver)
    @testset "Bose-Hubbard number-two sector" begin
        _,(b,bd) = bosonic_variables(2)
        numbers = bd.*b
        h = -(bd[1]*b[2]+bd[2]*b[1])+sum(n^2-n for n in numbers)
        # Old moment_eq_constraints also multiply N-2 by row words. Sector
        # preserves those equations; ExpectationZero alone does not.
        problem = Problem(h)
        matrix = [2. -sqrt(2) 0; -sqrt(2) 0 -sqrt(2);0 -sqrt(2) 2]
        exact = eigmin(Hermitian(matrix))
        @test exact ≈ 1-sqrt(5) atol=1e-12
        oracle = fixture("bose_hubbard","two_site_n2_u2_order2")
        built = NCTSSoS.build(problem,Relaxation(order=2,strengthen=Sector(sum(numbers),2)))
        result = solve(built;solver)
        @test result.bound ≈ oracle["objective"] atol=1e-6
        @test result.bound ≈ exact atol=1e-6
        structure = fixture("next_manybody","bose_hubbard_sector")
        @test NCTSSoS.block_sides(built) == oracle["sides"]
        @test result.sizes.moments == structure["moments"]
        mean = solve(Problem(h;constraints=(ExpectationZero(sum(numbers)-2),)),Relaxation(order=2);solver)
        @test mean.bound <= result.bound+1e-6
    end
end
