# Note: FastPolynomials is loaded by setup.jl
using .FastPolynomials:
    coefficients, monomials, terms, support, create_noncommutative_variables

@testset "Polynomial" begin
    @testset "Creation from Terms" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])  # degree 2
        m2 = Monomial{NonCommutativeAlgebra}([3])     # degree 1
        m3 = Monomial{NonCommutativeAlgebra}([1])     # degree 1

        t1 = Term(1.0, m1)
        t2 = Term(2.0, m2)
        t3 = Term(-3.0, m3)

        p = Polynomial([t1, t2, t3])

        # Polynomial should be sorted by monomial ordering: degree-first, then lexicographic
        # Order: [1] (deg 1), [3] (deg 1), [1,2] (deg 2)
        # But [1] < [3] lexicographically, so: [1], [3], [1,2] -> coeffs: -3.0, 2.0, 1.0
        @test length(terms(p)) == 3
        @test coefficients(p) == [-3.0, 2.0, 1.0]
    end

    @testset "Automatic Deduplication" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2])  # same as m1

        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        # Should combine duplicate monomials
        @test length(terms(p)) == 1
        @test coefficients(p) == [3.0]
    end

    @testset "Zero Coefficient Removal" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        # Zero coefficient term should be removed
        p = Polynomial([Term(0.0, m1), Term(1.0, m2)])
        @test length(terms(p)) == 1
        @test coefficients(p) == [1.0]

        # Cancellation should remove the term entirely
        p2 = Polynomial([Term(2.0, m1), Term(-2.0, m1)])
        @test iszero(p2)
    end

    @testset "Constant Polynomial" begin
        p_const = Polynomial{NonCommutativeAlgebra,Int64,Float64}(5.0)

        @test length(terms(p_const)) == 1
        @test coefficients(p_const) == [5.0]
        @test isone(monomials(p_const)[1])

        # Zero constant
        p_zero = Polynomial{NonCommutativeAlgebra,Int64,Float64}(0.0)
        @test iszero(p_zero)
    end

    @testset "From Monomial" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        p = Polynomial(m)

        @test length(terms(p)) == 1
        @test coefficients(p) == [1.0]
        @test monomials(p)[1] == m
    end

    @testset "From Term" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t = Term(3.0, m)
        p = Polynomial(t)

        @test length(terms(p)) == 1
        @test coefficients(p) == [3.0]
    end

    @testset "Accessors" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3])

        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        @test length(coefficients(p)) == 2
        @test length(monomials(p)) == 2
        @test support(p) == monomials(p)
    end

    @testset "Degree" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3, 4])

        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        @test degree(p) == 3  # max degree

        # Zero polynomial
        p_zero = zero(Polynomial{NonCommutativeAlgebra,Int64,Float64})
        @test degree(p_zero) == -1  # convention for zero polynomial
    end

    @testset "Variables" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3])

        p = Polynomial([Term(1.0, m1), Term(1.0, m2)])

        vars = variables(p)
        @test 1 in vars
        @test 2 in vars
        @test 3 in vars
        @test length(vars) == 3
    end

    @testset "Zero and One" begin
        T = Polynomial{NonCommutativeAlgebra,Int64,Float64}

        p_zero = zero(T)
        @test iszero(p_zero)
        @test isempty(terms(p_zero))

        p_one = one(T)
        @test isone(p_one)
        @test length(terms(p_one)) == 1
        @test isone(monomials(p_one)[1])
        @test coefficients(p_one) == [1.0]
    end

    @testset "Negation" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])  # degree 2
        m2 = Monomial{NonCommutativeAlgebra}([3])    # degree 1
        p = Polynomial([Term(2.0, m), Term(-3.0, m2)])

        p_neg = -p
        # Order: [3] (deg 1) before [1,2] (deg 2), so coeffs: 3.0, -2.0
        @test coefficients(p_neg) == [3.0, -2.0]
        @test monomials(p_neg) == monomials(p)
    end

    @testset "Polynomial Addition" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])
        m3 = Monomial{NonCommutativeAlgebra}([3])

        p1 = Polynomial([Term(1.0, m1), Term(2.0, m2)])
        p2 = Polynomial([Term(3.0, m2), Term(4.0, m3)])

        p_sum = p1 + p2

        @test length(terms(p_sum)) == 3
        # m2 should have coefficient 2 + 3 = 5
        monos = monomials(p_sum)
        coeffs = coefficients(p_sum)
        m2_idx = findfirst(==(m2), monos)
        @test coeffs[m2_idx] == 5.0
    end

    @testset "Polynomial Subtraction" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        p1 = Polynomial([Term(3.0, m1)])
        p2 = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        p_diff = p1 - p2

        # Should have m1 with coeff 2, m2 with coeff -2
        monos = monomials(p_diff)
        coeffs = coefficients(p_diff)
        m1_idx = findfirst(==(m1), monos)
        m2_idx = findfirst(==(m2), monos)
        @test coeffs[m1_idx] == 2.0
        @test coeffs[m2_idx] == -2.0
    end

    @testset "Scalar Multiplication" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        p = Polynomial([Term(2.0, m)])

        p_scaled = 3.0 * p
        @test coefficients(p_scaled) == [6.0]

        p_scaled_right = p * 3.0
        @test coefficients(p_scaled_right) == [6.0]
    end

    @testset "Scalar Division" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        p = Polynomial([Term(6.0, m)])

        p_div = p / 2.0
        @test coefficients(p_div) == [3.0]
    end

    @testset "Polynomial Power" begin
        m = Monomial{NonCommutativeAlgebra}(UInt16[1])
        p = Polynomial(Term(2.0, m))

        p0 = p^0
        @test isone(p0)

        p1 = p^1
        @test p1 == p

        p2 = p^2
        @test degree(p2) == 2
    end

    @testset "Equality" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3])

        p1 = Polynomial([Term(1.0, m1), Term(2.0, m2)])
        p2 = Polynomial([Term(1.0, m1), Term(2.0, m2)])
        p3 = Polynomial([Term(1.0, m1), Term(3.0, m2)])

        @test p1 == p2
        @test p1 != p3
    end

    @testset "Hash" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3])

        p1 = Polynomial([Term(1.0, m1), Term(2.0, m2)])
        p2 = Polynomial([Term(1.0, m1), Term(2.0, m2)])
        p3 = Polynomial([Term(1.0, m1), Term(3.0, m2)])

        @test hash(p1) == hash(p2)
        @test hash(p1) != hash(p3)

        # Unique should work correctly
        p_set = unique([p1, p2, p1, p3])
        @test length(p_set) == 2
    end

    @testset "Copy" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        p = Polynomial([Term(2.0, m)])

        p_copy = copy(p)
        @test p == p_copy
        @test p !== p_copy
    end

    @testset "Type Promotion in Arithmetic" begin
        m = Monomial{NonCommutativeAlgebra}([1])

        p_int = Polynomial([Term{Monomial{NonCommutativeAlgebra,Int64},Int}(2, m)])
        p_float = Polynomial([Term(1.0, m)])

        p_sum = p_int + p_float
        @test eltype(coefficients(p_sum)) == Float64
        @test coefficients(p_sum) == [3.0]
    end
end
