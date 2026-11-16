function free_fermion(solver)
    @testset "Free fermion hopping" begin
        for (n, order, id) in ((2,1,"two_site_hopping_order1"),
            (2,2,"two_site_hopping_order2"), (4,1,"pure_hopping_chain_n4_order1"))
            @testset "$id" begin
                # Independent single-particle spectrum: fill each negative level.
                # This gives the many-body ground energy with no particle restriction.
                levels = eigvals(SymTridiagonal(zeros(n), -ones(n-1)))
                exact = sum(min(0, e) for e in levels)
                if order == 1
                    @test exact ≈ (n == 2 ? -1 : -sqrt(5)) atol=(n == 2 ? 1e-12 : 1e-10)
                end
                _, (a, ad) = fermionic_variables(n)
                H = -sum(ad[i]*a[i+1]+ad[i+1]*a[i] for i in 1:n-1)
                built = Next.build(Problem(H), Relaxation(order=order))
                result = solve(built; solver)
                @test result.bound ≈ exact atol=1e-6
                @test result.bound ≈ fixture("free_fermion", id)["objective"] atol=1e-6
                structure = fixture("next_problem_pilot", "free_fermion/$id")
                @test sort(Next.block_sides(built)) == structure["sides"]
                @test result.sizes.moments == structure["moments"]
            end
        end
    end
end
