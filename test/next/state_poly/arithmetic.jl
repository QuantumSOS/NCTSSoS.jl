function state_arithmetic()
    @testset "Public state scalar and operator-tail laws" begin
        r,(x,)=site_variables(Free(),1,2;functional=Next.RealState())
        a,b=state_symbol(x[1]),state_symbol(x[2])
        word(p)=last(only(Next.terms(p)))
        terms(p)=Dict(w=>c for (c,w) in Next.terms(p))
        wa,wb,id=word(a),word(b),word(one(a))
        @test iszero(Next.polynomial(r,()))
        @test terms(a+2b)==terms(2b+a)==Dict(wa=>1,wb=>2)
        @test terms(a-2b)==Dict(wa=>1,wb=>-2)
        @test terms(2b-a)==Dict(wa=>-1,wb=>2)
        p=a+2b
        @test terms(p+2b)==terms(2b+p)==Dict(wa=>1,wb=>4)
        @test p-2b==a
        @test terms(a+p)==terms(p+a)==Dict(wa=>2,wb=>2)
        @test p-a==2b && a-p==-2b
        @test a*p==p*a==a^2+2a*b
        # The same scalar laws must hold with a nontrivial operator tail.
        for q in (a,a*x[2])
            wq=word(q)
            @test terms(2+q)==terms(q+2)==Dict(id=>2,wq=>1)
            @test terms(2-q)==Dict(id=>2,wq=>-1)
            @test terms(q-2)==Dict(id=>-2,wq=>1)
            @test terms(2q)==terms(q*2)==terms(q+q)==Dict(wq=>2)
            @test terms(-q)==Dict(wq=>-1)
        end
        @test terms(x[1]*p)==terms(p*x[1])
        @test length(Next.terms(x[1]*p))==2
        @test state_symbol(x[1]*p)==a^2+2a*b
        @test state_symbol(x[1]*x[2])!=a*b
    end
end

function state_chsh(sos_solver)
    @testset "Dense state CHSH fixture" begin
        _,(a,b)=site_variables(Unipotent(),2,2;functional=Next.RealState())
        h=-state_symbol(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
        built=Next.build(Problem(h),Relaxation(order=1))
        data=TOML.parsefile(joinpath(@__DIR__,"../../data/expectations/chsh_state.toml"))
        oracle=only(filter(c->c["id"]=="Dense",data["cases"]))["expected"]
        @test Next.block_sides(built)==oracle["sides"]
        @test built.sizes.moments==oracle["nuniq"]
        @test solve(built;solver=sos_solver).bound ≈ oracle["objective"] atol=1e-5 skip=isnothing(sos_solver)
    end
end
