function fermionic_parity(solver)
    @testset "Fermionic parity and mode swap" begin
        registry,(a,ad)=fermionic_variables(2)
        for nonhermitian in (a[1],ad[1],ad[1]*a[1]+a[1])
            @test_throws ArgumentError Problem(nonhermitian)
        end
        rows=collect(Next.basis(registry,1))
        @test any(p->isodd(Next.degree(p)),rows)
        @test any(p->iseven(Next.degree(p)),rows)
        h=-(ad[1]*a[2]+ad[2]*a[1])
        plain=Next.build(Problem(h),Relaxation(rows=rows))
        symmetric=Next.build(Problem(h;symmetry=(Next.NumberRotation(registry),SitePermutation([2,1]))),Relaxation(rows=rows))
        # The two-mode single-particle energies are -1,+1; fill the negative one.
        result=solve(symmetric;solver)
        @test result.bound ≈ -1 atol=1e-6
        @test maximum(Next.block_sides(symmetric)) < maximum(Next.block_sides(plain))
        @test expectation(result,a[1]+ad[1]) ≈ 0 atol=1e-6
    end
end
