function fermionic_parity(solver)
    @testset "Fermionic parity and mode swap" begin
        registry,(a,ad)=fermionic_variables(2)
        for nonhermitian in (a[1],ad[1],ad[1]*a[1]+a[1])
            @test_throws ArgumentError Problem(nonhermitian)
        end
        rows=collect(NCTSSoS.basis(registry,1))
        @test any(p->isodd(NCTSSoS.degree(p)),rows)
        @test any(p->iseven(NCTSSoS.degree(p)),rows)
        h=-(ad[1]*a[2]+ad[2]*a[1])
        plain=NCTSSoS.build(Problem(h),Relaxation(rows=rows))
        symmetric=NCTSSoS.build(Problem(h;symmetry=(NCTSSoS.NumberRotation(registry),SitePermutation([2,1]))),Relaxation(rows=rows))
        # The two-mode single-particle energies are -1,+1; fill the negative one.
        result=solve(symmetric;solver)
        @test result.bound ≈ -1 atol=1e-6
        @test maximum(NCTSSoS.block_sides(symmetric)) < maximum(NCTSSoS.block_sides(plain))
        @test expectation(result,a[1]+ad[1]) ≈ 0 atol=1e-6
    end
end
