using NCTSSoS, Test

@testset "Arithmetic" begin
    @testset "Monomial Multiplication" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Variable multiplication returns a simplified Monomial
        m1 = x[1]
        m2 = x[2]

        result = m1 * m2
        @test result isa Polynomial
        @test degree(result) == 2
    end

    @testset "Scalar * Monomial -> Polynomial" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))

        # Scalar times monomial as a (coef, monomial) pair
        p = Polynomial((3.0, m))

        @test p isa Polynomial
        @test coefficients(p) == [3.0]
    end

    @testset "Polynomial Addition" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        p1 = Polynomial((1.0, m1))
        p2 = Polynomial((2.0, m2))

        p_sum = p1 + p2
        @test length(monomials(p_sum)) == 2

        # Adding same monomial
        p3 = Polynomial((3.0, m1))
        p_combined = p1 + p3
        @test length(monomials(p_combined)) == 1
        @test coefficients(p_combined) == [4.0]
    end

    @testset "Polynomial Multiplication" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        p1 = Polynomial((2.0, m1))
        p2 = Polynomial((3.0, m2))

        p_prod = p1 * p2
        @test coefficients(p_prod) == [6.0]
        @test degree(p_prod) == 2
    end

    @testset "Monomial-Polynomial Interaction" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        p = Polynomial((2.0, m))

        # Polynomial can be scaled
        p_scaled = 3.0 * p
        @test coefficients(p_scaled) == [6.0]
    end

    @testset "Subtraction" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        p1 = Polynomial((5.0, m1))
        p2 = Polynomial((3.0, m1))

        p_diff = p1 - p2
        @test coefficients(p_diff) == [2.0]

        # Different monomials
        p3 = Polynomial((2.0, m2))
        p_diff2 = p1 - p3
        @test length(monomials(p_diff2)) == 2
    end

    @testset "Scaling Polynomial" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p = Polynomial([(1.0, m), (2.0, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2)))])

        p_scaled = 2.0 * p
        @test coefficients(p_scaled) == [2.0, 4.0]

        # Type promotion
        p_scaled_float32 = Float32(2.0) * p
        @test coefficients(p_scaled_float32) == [2.0, 4.0]
    end

    @testset "Polynomial Multiplication Distributivity" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        # (a + b) * c = a*c + b*c
        p_ab = Polynomial([(1.0, m1), (2.0, m2)])
        p_c = Polynomial([(3.0, m1)])

        p_prod = p_ab * p_c

        # Should have terms from both products
        @test degree(p_prod) == 2
    end

    @testset "Identity Multiplication" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p = Polynomial((2.0, m))

        p_one = one(typeof(p))

        @test p * p_one == p
        @test p_one * p == p
    end

    @testset "Zero Multiplication" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p = Polynomial((2.0, m))

        p_zero = zero(typeof(p))

        @test iszero(p * p_zero)
        @test iszero(p_zero * p)
    end

    @testset "Type Promotion" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))

        p_int = Polynomial([(2, m)])
        p_float = Polynomial([(1.5, m)])

        p_sum = p_int + p_float
        @test eltype(coefficients(p_sum)) <: AbstractFloat
    end
end
