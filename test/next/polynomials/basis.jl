function basis_contracts()
    @testset "Two-site basis spans" begin
        for rule in (Free(),Projector(),Unipotent())
            registry,_ = site_variables(rule,2,2)
            rows = collect(NCTSSoS.basis(registry,4))
            # Convolution of the two single-site word counts through degree 4:
            # Free: sum((k+1)*2^k)=129; reduced: 1+4+8+12+16=41.
            @test length(rows) == (rule isa Free ? 129 : 41)
            @test length(unique(rows)) == length(rows)
            @test all(p->NCTSSoS.degree(p)<=4,rows)
        end
        for (registry,count) in ((fermionic_variables(2)[1],15),(bosonic_variables(2)[1],35))
            # Four CAR letters without repetition: sum(binomial(4,k), k=0:3).
            # Four CCR letters with repetition: binomial(4+3,3).
            rows = collect(NCTSSoS.basis(registry,3))
            @test length(rows) == count
            @test length(unique(rows)) == length(rows)
        end
    end
end
