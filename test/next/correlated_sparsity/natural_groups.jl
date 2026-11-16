function natural_groups()
    @testset "Natural-order groups retain constraints and unused letters" begin
        registry,(x,)=site_variables(Free(),1,5)
        f=x[1]*x[4]+x[2]*x[4]
        @test_throws ArgumentError Problem(f)
        built=NCTSSoS.build(Problem((f+f')/2;constraints=(Zero(x[2]+x[3]),)),
            Relaxation(order=1,grouping=Correlative(algorithm=NaturalElimination())))
        position=Dict(NCTSSoS.operator_name(registry,last(only(NCTSSoS.terms(p))))=>i for (i,p) in enumerate(x))
        groups=[sort!(unique!([position[NCTSSoS.operator_name(registry,l)] for w in block for l in NCTSSoS.word_letters(w)]))
            for block in collect(NCTSSoS.block_rows(built))[1:3]]
        @test sort(groups)==[[1,4],[2,3,4],[5]]
        @test length.(collect(NCTSSoS.block_rows(built)))==[3,4,2,1]
    end
    @testset "Natural-order projector Bell groups" begin
        for n in (5,6),order in (2,3)
            registry,(a,b)=site_variables(Projector(),2,n)
            f=sum(a[i]*b[j] for i in 1:n for j in 1:n if i==n || j==n || i+j>=n)
            built=NCTSSoS.build(Problem(f),Relaxation(order=order,grouping=Correlative(algorithm=NaturalElimination())))
            expected=[vcat(i,collect((n+max(1,n-i)):2n)) for i in 1:n]
            # The old registry puts all Alice letters before Bob's; basis()
            # need not use that order. Compare names in this explicit order.
            position=Dict(NCTSSoS.operator_name(registry,last(only(NCTSSoS.terms(p))))=>i for (i,p) in enumerate([a;b]))
            groups=[sort!(unique!([position[NCTSSoS.operator_name(registry,l)] for w in block for l in NCTSSoS.word_letters(w)]))
                for block in NCTSSoS.block_rows(built)]
            @test sort(groups)==expected
            sizes=[sum(k==0 ? 1 : m*(m-1)^(k-1) for k in 0:order)+
                sum(k==0 ? 1 : m*(m-1)^(k-1) for k in 0:order-1) for m in min.(2:n+1,n)]
            @test NCTSSoS.block_sides(built)==sizes
        end
    end
end
