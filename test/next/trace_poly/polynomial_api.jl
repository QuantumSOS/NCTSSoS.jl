function polynomial_contracts()
    @testset "Public scalar and mixed functional arithmetic" begin
        for (functional,symbol) in ((Next.RealState(),state_symbol),(Next.RealTrace(),trace_symbol))
            registry,(x,) = site_variables(Free(),1,5;functional)
            a,b = symbol(x[1]*x[2]),symbol(x[3]*x[4]*x[5])
            # Different degrees and coefficients expose wrong factor merging.
            @test Next.polynomial(registry,(a,2a,zero(a),b,-b)) == 3a
            @test Next.degree(a+b) == 3
            @test 3a*(2b)-a*b == 5a*b
            mixed = a*x[3]*x[4]
            @test Next.degree(mixed) == 4
            names = Set(Next.operator_name(registry,l) for (_,w) in Next.terms(mixed) for l in Next.support_letters(w))
            @test names == Set((:s1_1,:s1_2,:s1_3,:s1_4))
            # Old D-file accessors become value checks, not concrete types.
            @test Next.degree(a) == 2
            @test Set(Next.operator_name(registry,l) for (_,w) in Next.terms(a)
                for l in Next.support_letters(w)) == Set((:s1_1,:s1_2))
            for p in (a,mixed)
                @test !isone(p) && isone(one(p))
                @test iszero(zero(p)) && iszero(symbol(zero(p)))
                @test collect(Next.terms(2p)) == [(2,last(only(Next.terms(p))))]
            end
            @test a*(2a) == (2a)*a
        end
    end
end
