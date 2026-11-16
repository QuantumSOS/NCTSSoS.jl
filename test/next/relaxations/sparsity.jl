function sparsity(solver)
    data = TOML.parsefile(joinpath(@__DIR__,"../../data/expectations/relaxations_sparsity.toml"))
    @testset "Dense involution fixture contracts" begin
        for n in (2,1)
            _,(u,) = site_variables(Unipotent(),1,n)
            if n==2
                raw = u[1]*u[2]
                @test_throws ArgumentError Problem(raw)
                # The old real form identifies <w> and <w*>. Its objective
                # therefore means the Hermitian part, made explicit here.
                objective = (raw+raw')/2
            else
                objective = u[1]
            end
            oracle = only(filter(c->c["id"]=="unipotent_n$(n)_order1_dense",data["cases"]))["expected"]
            built = Next.build(Problem(objective),Relaxation(order=1))
            result = solve(built;solver)
            # No bound assertion/tolerance existed in the old file. Use the
            # unchanged fixture objective with Julia's default isapprox.
            @test result.bound ≈ oracle["objective"]
            @test Next.block_sides(built) == oracle["sides"]
            @test result.sizes.moments == oracle["nuniq"]
        end
    end
end
