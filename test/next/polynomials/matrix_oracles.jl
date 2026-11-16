# Evaluate the public term iterator with independent finite matrices. Neither
# the matrix product nor the expected value uses a package simplifier.
function matrix_value(p, registry, operators, dim)
    result = zeros(ComplexF64, dim, dim)
    for (coefficient, word) in Next.terms(p)
        term = Matrix{ComplexF64}(I, dim, dim)
        for letter in Next.word_letters(word)
            term *= operators[Next.operator_name(registry, letter)]
        end
        result .+= coefficient .* term
    end
    return result
end

function matrix_oracles()
    @testset "Independent Pauli and Jordan–Wigner matrices" begin
        eye = Matrix{ComplexF64}(I, 2, 2)
        X, Y, Z = ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1]
        registry, (x,y,z) = pauli_variables(2)
        letters = [x[1],y[1],z[1],x[2],y[2],z[2]]
        matrices = [kron(X,eye),kron(Y,eye),kron(Z,eye),kron(eye,X),kron(eye,Y),kron(eye,Z)]
        operators = Dict(Next.operator_name(registry,last(only(Next.terms(p)))) => M
            for (p,M) in zip(letters,matrices))
        for a in 1:6, b in 1:6, c in 1:6
            p = letters[a]*letters[b]*letters[c]
            @test length(Next.terms(p)) == 1
            @test isapprox(matrix_value(p,registry,operators,4), matrices[a]*matrices[b]*matrices[c]; atol=1e-12,rtol=0)
        end

        registry, (a,ad) = fermionic_variables(3)
        annihilators = [reduce(kron,[j<i ? Z : j==i ? ComplexF64[0 0; 1 0] : eye for j in 1:3]) for i in 1:3]
        ops = Dict([i=>a[i] for i in 1:3]..., [-i=>ad[i] for i in 1:3]...)
        mats = Dict([i=>annihilators[i] for i in 1:3]..., [-i=>Matrix(annihilators[i]') for i in 1:3]...)
        operators = Dict(Next.operator_name(registry,last(only(Next.terms(p)))) => mats[i] for (i,p) in ops)
        # The four original matrix cases, plus the surplus, balanced and full
        # contraction cases from the old simplify tests. Keep their signs/order.
        for word in ([1,1],[-1,-1],[1,-1],[1,-2], [1,-1,1,-1],
            [1,2,1],[-1,-2,-1],[1,1,-1],[-1,1,-1],[1,2,3,1,2],
            [1,2,-1,-2],[1,-1,2,-2],[1,2,3,-3,-2,-1])
            p = foldl(*,(ops[i] for i in word);init=one(a[1]))
            expected = foldl(*,(mats[i] for i in word);init=Matrix{ComplexF64}(I,8,8))
            @test isapprox(matrix_value(p,registry,operators,8),expected;atol=1e-12,rtol=0)
        end
    end
end
