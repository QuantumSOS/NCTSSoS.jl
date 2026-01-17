using Test
using NCTSSoS
using NCTSSoS: decode_site, decode_operator_id

@testset "Monoid algebra showcase" begin
    @testset "Anatomy (sites + decoding)" begin
        reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])

        @test x[1] isa NormalMonomial
        @test decode_site(x[1].word[1]) == 1
        @test decode_site(y[1].word[1]) == 2

        lhs = x[2] * y[1] * x[1] * y[2]
        w = monomials(lhs)[1].word

        @test decode_site.(w) == [1, 1, 2, 2]
        @test decode_operator_id.(w) == [
            decode_operator_id(x[2].word[1]),
            decode_operator_id(x[1].word[1]),
            decode_operator_id(y[1].word[1]),
            decode_operator_id(y[2].word[1]),
        ]

        @test !isempty(indices(reg))
    end

    @testset "NonCommutativeAlgebra" begin
        reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])

        # Different sites commute (x-site group before y-site group)
        @test x[2] * y[1] * x[1] * y[2] == x[2] * x[1] * y[1] * y[2]

        # Within a site, no relations
        @test x[1] * x[2] != x[2] * x[1]

        # Registry is used (silence unused warning)
        @test !isempty(indices(reg))

        basis_d1 = get_ncbasis(reg, 1)
        @test length(basis_d1) == 1 + length(indices(reg))
    end

    @testset "ProjectorAlgebra" begin
        reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])

        @test monomials(P[1] * P[1]) == [P[1]]
        @test monomials(P[1] * P[1] * P[1]) == [P[1]]

        @test monomials(Q[2] * P[1]) == monomials(P[1] * Q[2])
        @test monomials(Q[1] * P[1] * Q[1] * P[1]) == monomials(P[1] * Q[1])
        @test !isempty(indices(reg))
    end

    @testset "UnipotentAlgebra" begin
        reg, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 1:2)])
        T = eltype(indices(reg))
        ID = one(NormalMonomial{UnipotentAlgebra,T})

        @test monomials(U[1] * U[1]) == [ID]
        @test monomials(U[1] * U[1] * U[1]) == [U[1]]

        @test monomials(V[1] * U[1] * V[1] * U[1]) == [ID]
        @test monomials(U[1] * U[1] * U[2]) == monomials(U[2])

        @test monomials(V[2] * U[1]) == monomials(U[1] * V[2])
        @test !isempty(indices(reg))
    end
end
