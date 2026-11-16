# Evaluate the public term iterator with independent finite matrices. Neither
# the matrix product nor the expected value uses a package simplifier.
function matrix_value(p, registry, operators, dim)
    result = zeros(ComplexF64, dim, dim)
    for (coefficient, word) in NCTSSoS.terms(p)
        term = Matrix{ComplexF64}(I, dim, dim)
        for letter in NCTSSoS.word_letters(word)
            term *= operators[NCTSSoS.operator_name(registry, letter)]
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
        operators = Dict(NCTSSoS.operator_name(registry,last(only(NCTSSoS.terms(p)))) => M
            for (p,M) in zip(letters,matrices))
        for a in 1:6, b in 1:6, c in 1:6
            p = letters[a]*letters[b]*letters[c]
            @test length(NCTSSoS.terms(p)) == 1
            @test isapprox(matrix_value(p,registry,operators,4), matrices[a]*matrices[b]*matrices[c]; atol=1e-12,rtol=0)
            @test isapprox(matrix_value(p',registry,operators,4), matrix_value(p,registry,operators,4)'; atol=1e-12,rtol=0)
        end

        registry, (a,ad) = fermionic_variables(3)
        annihilators = [reduce(kron,[j<i ? Z : j==i ? ComplexF64[0 0; 1 0] : eye for j in 1:3]) for i in 1:3]
        ops = Dict([i=>a[i] for i in 1:3]..., [-i=>ad[i] for i in 1:3]...)
        mats = Dict([i=>annihilators[i] for i in 1:3]..., [-i=>Matrix(annihilators[i]') for i in 1:3]...)
        operators = Dict(NCTSSoS.operator_name(registry,last(only(NCTSSoS.terms(p)))) => mats[i] for (i,p) in ops)
        # All 1+6+36+216+1296 words, including the empty word.
        for len in 0:4, word in Iterators.product(ntuple(_->(1,2,3,-1,-2,-3),len)...)
            p = foldl(*,(ops[i] for i in word);init=one(a[1]))
            expected = foldl(*,(mats[i] for i in word);init=Matrix{ComplexF64}(I,8,8))
            @test isapprox(matrix_value(p,registry,operators,8),expected;atol=1e-12,rtol=0)
            @test isapprox(matrix_value(p',registry,operators,8),matrix_value(p,registry,operators,8)';atol=1e-12,rtol=0)
        end
    end
end

function site_matrix_oracles()
    @testset "Independent two-site matrix representations" begin
        rng = MersenneTwister(20261116)
        dim = 3
        eye = Matrix{ComplexF64}(I,dim,dim)
        for rule in (Free(),Projector(),Unipotent())
            registry,(a,b) = site_variables(rule,2,2)
            letters = [a;b]
            local_matrices = map(1:4) do i
                raw = randn(rng,ComplexF64,dim,dim)
                rule isa Free && return (raw+raw')/2
                q = Matrix(qr(raw).Q)
                # Alternating ranks avoid identity/zero projectors.
                columns = q[:,1:(isodd(i) ? 1 : 2)]
                projector = columns*columns'
                return rule isa Projector ? projector : 2projector-eye
            end
            matrices = [i<=2 ? kron(m,eye) : kron(eye,m) for (i,m) in enumerate(local_matrices)]
            operators = Dict(NCTSSoS.operator_name(registry,last(only(NCTSSoS.terms(p))))=>m
                for (p,m) in zip(letters,matrices))
            for len in 0:4, word in Iterators.product(ntuple(_->1:4,len)...)
                p = foldl(*,(letters[i] for i in word);init=one(a[1]))
                expected = foldl(*,(matrices[i] for i in word);init=Matrix{ComplexF64}(I,dim^2,dim^2))
                actual = matrix_value(p,registry,operators,dim^2)
                @test isapprox(actual,expected;atol=1e-10,rtol=0)
                @test isapprox(matrix_value(p',registry,operators,dim^2),actual';atol=1e-10,rtol=0)
            end
        end
    end
end
