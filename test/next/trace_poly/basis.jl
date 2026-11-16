function basis_contracts()
    @testset "Trace basis spans" begin
        for (sites,letters,count) in ((1,4,53),(1,2,15),(2,2,49))
            registry,_ = site_variables(Unipotent(),sites,letters;functional=Next.RealTrace())
            rows = collect(Next.basis(registry,2))
            @test length(rows) == count
            @test length(unique(rows)) == count
            @test all(p->Next.degree(p)<=2,rows)
        end
    end
end
