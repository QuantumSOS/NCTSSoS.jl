function particle_number(solver)
    @testset "Total and grouped mean particle numbers" begin
        _,(b,bd) = bosonic_variables(3)
        number = sum(bd[i]*b[i] for i in eachindex(b))
        result = solve(Problem(number;constraints=(ExpectationZero(number-2),)),
            Relaxation(order=1);solver)
        @test result.bound ≈ 2 atol=1e-6
        @test expectation(result,number) ≈ 2 atol=1e-6

        _,(a,ad) = fermionic_variables(4)
        up = ad[1]*a[1]+ad[2]*a[2]
        down2 = ad[4]*a[4]
        # Asymmetric groups: the third mode is free, not part of down2.
        # The mean can be fractional; this is not a fixed-number Sector.
        objective = 2up+3down2-ad[3]*a[3]
        result = solve(Problem(objective;constraints=(ExpectationZero(up-1),
            ExpectationZero(down2-0.25))),Relaxation(order=1);solver)
        # Occupy mode 3 and mix mode 4 with probability 1/4: 2+3/4-1.
        @test result.bound ≈ 1.75 atol=1e-6
        @test expectation(result,up) ≈ 1 atol=1e-6
        @test expectation(result,down2) ≈ 0.25 atol=1e-6
        @test isempty(result.sizes.localizing_rows)
    end
end
