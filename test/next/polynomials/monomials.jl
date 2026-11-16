function monomial_powers()
    @testset "Powers of noncommuting words" begin
        _,(x,) = site_variables(Free(),1,2)
        p = x[1]*x[2]
        @test p^0 == one(p)
        @test p^1 == p
        @test p^2 == x[1]*x[2]*x[1]*x[2]
        @test NCTSSoS.degree(p^3) == 6
        @test one(p)^10 == one(p)
        @test p^3 == p^2*p == p*p^2
        _,(a,ad) = fermionic_variables(1)
        n = ad[1]*a[1]
        @test n^2 == n
    end
end
