# Matrices and the Jordan–Wigner map are independent of the polynomial package.
function ising(solver)
    @testset "Ising and Jordan–Wigner" begin
        n = 3
        J, h = 1.0, 2.0
        eye = Matrix{ComplexF64}(I, 2, 2)
        X, Z = ComplexF64[0 1; 1 0], ComplexF64[1 0; 0 -1]
        localop(A, i) = reduce(kron, [j == i ? A : eye for j in 1:n])
        matrix(A, B, periodic) = -J/4 * sum(localop(A,i)*localop(A,mod1(i+1,n))
            for i in 1:(periodic ? n : n-1)) - h/2 * sum(localop(B,i) for i in 1:n)
        _, (sx, _, sz) = pauli_variables(n)
        for (periodic, id) in ((false,"open_pauli_n3_order2"), (true,"periodic_pauli_n3_order2"))
            oracle = fixture("ising", id)
            @test eigmin(Hermitian(matrix(Z, X, periodic))) ≈ oracle["objective"] atol=1e-12
            H = -J/4 * sum(sz[i]*sz[mod1(i+1,n)] for i in 1:(periodic ? n : n-1)) -
                h/2 * sum(sx)
            built = NCTSSoS.build(Problem(H), Relaxation(order=2))
            result = solve(built; solver)
            @test result.bound ≈ oracle["objective"] atol=1e-6
            structure = fixture("next_problem_pilot", "ising/$id")
            @test sort(NCTSSoS.block_sides(built)) == structure["sides"]
            @test result.sizes.moments == structure["moments"]
        end
        rotated = matrix(X, Z, false)
        exact = eigmin(Hermitian(rotated))
        @test exact ≈ fixture("ising", "open_pauli_n3_order2")["objective"] atol=1e-12
        oracle = fixture("ising", "open_fermionic_n3_order1")
        @test exact ≈ oracle["objective"] atol=1e-12
        r, (a, ad) = fermionic_variables(n)
        H = J/4 * sum((ad[i]-a[i])*(ad[i+1]+a[i+1]) for i in 1:n-1) -
            h/2 * sum(2ad[i]*a[i]-1 for i in 1:n)
        annihilators = [reduce(kron, [j < i ? Z : j == i ? ComplexF64[0 0; 1 0] : eye
            for j in 1:n]) for i in 1:n]
        operators = Dict(NCTSSoS.operator_name(r, last(only(NCTSSoS.terms(p)))) => M
            for (p,M) in zip([a; ad], [annihilators; adjoint.(annihilators)]))
        evaluated = sum(c * foldl(*, (operators[NCTSSoS.operator_name(r,l)] for l in NCTSSoS.word_letters(w));
            init=Matrix{ComplexF64}(I, 2^n, 2^n)) for (c,w) in NCTSSoS.terms(H))
        @test isapprox(evaluated, rotated; atol=1e-12, rtol=0)
        built = NCTSSoS.build(Problem(H), Relaxation(order=1))
        result = solve(built; solver)
        @test result.bound ≈ oracle["objective"] atol=1e-6
        structure = fixture("next_problem_pilot", "ising/open_fermionic_n3_order1")
        @test sort(NCTSSoS.block_sides(built)) == structure["sides"]
        @test result.sizes.moments == structure["moments"]
    end
end
