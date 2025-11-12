using Test, NCTSSoS.FastPolynomials

@testset "Comparison" begin
    @testset "cmp variables" begin
        @ncpolyvar x y z

        @test cmp(x, x) == 0
        @test cmp(x, y) == -1
        @test cmp(y, x) == 1

        @test isless(x, y) == true
        @test sort([z, y, x]) == [x, y, z]

        @ncpolyvar x[1:10]
        @test cmp(x[1], x[10]) == -1

        @test x[1] == x[1]
        @test x[1] != x[2]

        @test x in [x, y, z]
        @test x ∉ [y, z]
    end

    @testset "cmp monomials" begin
        @ncpolyvar x y z

        mono1 = Monomial([x, y], [1, 2])
        mono2 = Monomial([x, y], [1, 2])
        mono3 = Monomial([x, y], [1, 1])
        mono4 = Monomial([x, y], [1, 2])
        mono5 = Monomial([x, z, y], [1, 1, 1])
        mono6 = Monomial([x], [2])

        # TODO: add more tests

        @test cmp(mono1, mono2) == 0
        @test cmp(mono1, mono3) == 1
        @test cmp(mono1, mono4) == 0
        @test cmp(mono4, mono5) == -1
        @test cmp(mono6, mono3) == -1

        @test isless(mono4, mono5)

        @test mono1 in sort([mono1, mono2, mono3])
        @test mono1 ∉ [mono3, mono3]
    end

    @testset "compare polynomials" begin
        @ncpolyvar x y z
        p1 = Polynomial([1.0, 2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
        p2 = Polynomial([1.0, 2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
        p3 = Polynomial(
            [1.00000001, 2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])]
        )
        @test p1 == p2
        @test !(p1 == p3)
    end
    @testset "Hash Polynomial" begin
        @ncpolyvar x y
        p1 = Polynomial([1.0, 2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
        p2 = Polynomial([1.0, 2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])])
        @test hash(p1) == hash(p2)

        p3 = Polynomial(
            [1.00000001, 2.0], [Monomial([x, y], [1, 2]), Monomial([x, y], [1, 1])]
        )
        p_set = unique!([p1, p2, p1, p3])
        @test length(p_set) == 2
    end
end
