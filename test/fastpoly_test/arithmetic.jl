# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    coefficients, monomials, create_noncommutative_variables

@testset "Arithmetic" begin
    @testset "Monomial Multiplication" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Variable multiplication creates polynomial via Term
        m1 = x[1]
        m2 = x[2]

        result = m1 * m2
        @test result isa Term
        @test degree(result.monomial) == 2
    end

    @testset "Scalar * Monomial -> Polynomial" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])

        # Scalar times monomial gives Term, which can be converted to Polynomial
        t = Term(3.0, m)
        p = Polynomial(t)

        @test p isa Polynomial
        @test coefficients(p) == [3.0]
    end

    @testset "Polynomial Addition" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        p1 = Polynomial(Term(1.0, m1))
        p2 = Polynomial(Term(2.0, m2))

        p_sum = p1 + p2
        @test length(monomials(p_sum)) == 2

        # Adding same monomial
        p3 = Polynomial(Term(3.0, m1))
        p_combined = p1 + p3
        @test length(monomials(p_combined)) == 1
        @test coefficients(p_combined) == [4.0]
    end

    @testset "Polynomial Multiplication" begin
        m1 = Monomial{NonCommutativeAlgebra}(UInt16[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt16[2])

        p1 = Polynomial(Term(2.0, m1))
        p2 = Polynomial(Term(3.0, m2))

        p_prod = p1 * p2
        @test coefficients(p_prod) == [6.0]
        @test degree(p_prod) == 2
    end

    @testset "Monomial-Polynomial Interaction" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        p = Polynomial(Term(2.0, m))

        # Polynomial can be scaled
        p_scaled = 3.0 * p
        @test coefficients(p_scaled) == [6.0]
    end

    @testset "Subtraction" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        p1 = Polynomial(Term(5.0, m1))
        p2 = Polynomial(Term(3.0, m1))

        p_diff = p1 - p2
        @test coefficients(p_diff) == [2.0]

        # Different monomials
        p3 = Polynomial(Term(2.0, m2))
        p_diff2 = p1 - p3
        @test length(monomials(p_diff2)) == 2
    end

    @testset "Scaling Polynomial" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        p = Polynomial([Term(1.0, m), Term(2.0, Monomial{NonCommutativeAlgebra}([2]))])

        p_scaled = 2.0 * p
        @test coefficients(p_scaled) == [2.0, 4.0]

        # Type promotion
        p_scaled_float32 = Float32(2.0) * p
        @test coefficients(p_scaled_float32) == [2.0, 4.0]
    end

    @testset "Polynomial Multiplication Distributivity" begin
        m1 = Monomial{NonCommutativeAlgebra}(UInt16[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt16[2])

        # (a + b) * c = a*c + b*c
        p_ab = Polynomial([Term(1.0, m1), Term(2.0, m2)])
        p_c = Polynomial([Term(3.0, m1)])

        p_prod = p_ab * p_c

        # Should have terms from both products
        @test degree(p_prod) == 2
    end

    @testset "Identity Multiplication" begin
        m = Monomial{NonCommutativeAlgebra}(UInt16[1])
        p = Polynomial(Term(2.0, m))

        p_one = one(typeof(p))

        @test p * p_one == p
        @test p_one * p == p
    end

    @testset "Zero Multiplication" begin
        m = Monomial{NonCommutativeAlgebra}(UInt16[1])
        p = Polynomial(Term(2.0, m))

        p_zero = zero(typeof(p))

        @test iszero(p * p_zero)
        @test iszero(p_zero * p)
    end

    @testset "Type Promotion" begin
        m = Monomial{NonCommutativeAlgebra}([1])

        p_int = Polynomial([Term{Monomial{NonCommutativeAlgebra,Int64},Int}(2, m)])
        p_float = Polynomial([Term(1.5, m)])

        p_sum = p_int + p_float
        @test eltype(coefficients(p_sum)) <: AbstractFloat
    end
end
