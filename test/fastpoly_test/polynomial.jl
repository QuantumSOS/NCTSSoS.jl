# Note: FastPolynomials is loaded by setup.jl
using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: variable_indices

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

    @testset "variable_indices" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3])

        p = Polynomial([Term(1.0, m1), Term(1.0, m2)])

        var_idxs = variable_indices(p)
        # variable_indices() returns Set{T} of indices
        @test length(var_idxs) == 3
        @test 1 in var_idxs
        @test 2 in var_idxs
        @test 3 in var_idxs
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

    @testset "Display (show)" begin
        # Zero polynomial
        p_zero = zero(Polynomial{NonCommutativeAlgebra,Int64,Float64})
        @test repr(p_zero) == "0"

        # Single term
        m = Monomial{NonCommutativeAlgebra}([1])
        p_single = Polynomial([Term(2.0, m)])
        output = repr(p_single)
        @test contains(output, "2.0")

        # Multiple terms
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])
        p_multi = Polynomial([Term(1.0, m1), Term(3.0, m2)])
        output_multi = repr(p_multi)
        # Tests that both terms appear and + separator exists
        @test contains(output_multi, "+")
        @test contains(output_multi, "[1]")
        @test contains(output_multi, "[2]")
    end

    @testset "Adjoint Operations" begin
        # Test with PauliAlgebra (self-adjoint generators)
        m_pauli = Monomial{PauliAlgebra}([1, 2, 3])
        p_pauli = Polynomial([Term(1.0 + 2.0im, m_pauli)])

        p_adj = adjoint(p_pauli)

        # Coefficient should be conjugated
        @test coefficients(p_adj)[1] == 1.0 - 2.0im

        # Monomial word should be reversed and negated for PauliAlgebra
        @test monomials(p_adj)[1].word == [-3, -2, -1]

        # Test with NonCommutativeAlgebra
        m_nc = Monomial{NonCommutativeAlgebra}([1, 2])
        p_nc = Polynomial([Term(3.0 + 4.0im, m_nc)])
        p_nc_adj = adjoint(p_nc)
        @test coefficients(p_nc_adj)[1] == 3.0 - 4.0im
        # NonCommutativeAlgebra uses negation for adjoint
        @test monomials(p_nc_adj)[1].word == [-2, -1]

        # Zero polynomial adjoint
        p_zero = zero(Polynomial{PauliAlgebra,Int64,ComplexF64})
        @test iszero(adjoint(p_zero))
    end

    @testset "Star (alias for adjoint)" begin
        m = Monomial{PauliAlgebra}([1, 2])
        p = Polynomial([Term(1.0 + 1.0im, m)])

        @test star(p) == adjoint(p)
        @test coefficients(star(p))[1] == 1.0 - 1.0im
    end

    @testset "is_symmetric" begin
        # Identity is symmetric
        p_identity = one(Polynomial{PauliAlgebra,Int64,ComplexF64})
        @test is_symmetric(p_identity)

        # Real scalar is symmetric
        p_real = Polynomial{PauliAlgebra,Int64,ComplexF64}(3.0)
        @test is_symmetric(p_real)

        # Complex coefficient breaks symmetry
        m = Monomial{PauliAlgebra}([1])
        p_complex = Polynomial([Term(1.0 + 1.0im, m)])
        @test !is_symmetric(p_complex)

        # Non-self-adjoint monomial
        m_non_self = Monomial{PauliAlgebra}([1, 2])
        p_non_sym = Polynomial([Term(1.0 + 0.0im, m_non_self)])
        @test !is_symmetric(p_non_sym)

        # Sum of term and its adjoint with proper negation
        # For NonCommutativeAlgebra: adjoint reverses and negates
        # [1,2] -> [-2,-1], so we need both terms
        m_pos = Monomial{NonCommutativeAlgebra}([1, 2])
        m_neg = Monomial{NonCommutativeAlgebra}([-2, -1])
        p_hermitian = Polynomial([Term(1.0 + 0.0im, m_pos), Term(1.0 + 0.0im, m_neg)])
        @test is_symmetric(p_hermitian)
    end

    @testset "variable_indices" begin
        # Multiple variables
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3])
        p = Polynomial([Term(1.0, m1), Term(1.0, m2)])

        var_set = variable_indices(p)
        @test var_set == Set([1, 2, 3])

        # Empty polynomial
        p_zero = zero(Polynomial{NonCommutativeAlgebra,Int64,Float64})
        @test isempty(variable_indices(p_zero))

        # Single variable repeated
        m_single = Monomial{NonCommutativeAlgebra}([1, 1, 1])
        p_single = Polynomial([Term(2.0, m_single)])
        @test variable_indices(p_single) == Set([1])

        # Test with signed types (fermionic/bosonic)
        m_signed = Monomial{FermionicAlgebra}([-1, 2, -3])  # negative = annihilation
        p_signed = Polynomial([Term(1.0, m_signed)])
        var_set_signed = variable_indices(p_signed)
        @test var_set_signed == Set([1, 2, 3])  # abs() applied
    end

    @testset "Error Paths" begin
        # Division by zero
        m = Monomial{NonCommutativeAlgebra}(UInt8[1])
        p = Polynomial([Term(6.0, m)])
        @test_throws DivideError p / 0.0

        # Negative power - actually throws MethodError (no inv method for Polynomial)
        # This tests that negative powers are not supported
        @test_throws MethodError p^(-1)
        @test_throws MethodError p^(-5)

        # Power of 0 and 1 should work
        @test isone(p^0)
        @test p^1 == p
    end

    @testset "Edge Cases: zero and one from instance" begin
        # zero from instance
        m = Monomial{NonCommutativeAlgebra}(UInt8[1, 2])
        p = Polynomial([Term(5.0, m)])
        p_zero_from_inst = zero(p)
        @test iszero(p_zero_from_inst)
        @test typeof(p_zero_from_inst) == typeof(p)

        # one from instance
        p_one_from_inst = one(p)
        @test isone(p_one_from_inst)
        @test typeof(p_one_from_inst) == typeof(p)
    end

    @testset "Scalar + Polynomial" begin
        # Scalar + polynomial (both orderings)
        m = Monomial{NonCommutativeAlgebra}(UInt8[1])
        p = Polynomial([Term(2.0, m)])

        # Polynomial + scalar
        p_plus_scalar = p + 3.0
        @test length(terms(p_plus_scalar)) == 2  # One constant term, one m1 term
        monos = monomials(p_plus_scalar)
        coeffs = coefficients(p_plus_scalar)
        # Identity monomial should appear
        id_idx = findfirst(isone, monos)
        @test !isnothing(id_idx)
        @test coeffs[id_idx] == 3.0

        # Scalar + polynomial
        scalar_plus_p = 3.0 + p
        @test scalar_plus_p == p_plus_scalar
    end

    @testset "Scalar - Polynomial and Polynomial - Scalar" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        p = Polynomial([Term(2.0, m)])

        # Polynomial - scalar
        p_minus_scalar = p - 3.0
        monos = monomials(p_minus_scalar)
        coeffs = coefficients(p_minus_scalar)
        id_idx = findfirst(isone, monos)
        @test coeffs[id_idx] == -3.0

        # Scalar - polynomial
        scalar_minus_p = 5.0 - p
        @test length(terms(scalar_minus_p)) == 2
        # Identity term should have 5.0, m term should have -2.0
        monos2 = monomials(scalar_minus_p)
        coeffs2 = coefficients(scalar_minus_p)
        id_idx2 = findfirst(isone, monos2)
        m_idx = findfirst(==(m), monos2)
        @test coeffs2[id_idx2] == 5.0
        @test coeffs2[m_idx] == -2.0
    end

    @testset "Cross-Algebra Equality (should be false)" begin
        m_pauli = Monomial{PauliAlgebra}([1])
        m_nc = Monomial{NonCommutativeAlgebra}([1])

        p_pauli = Polynomial([Term(1.0, m_pauli)])
        p_nc = Polynomial([Term(1.0, m_nc)])

        # Different algebra types should never be equal
        @test p_pauli != p_nc
    end

    @testset "Polynomial Multiplication (PauliAlgebra with simplification)" begin
        # PauliAlgebra: σ_i^2 = 1, anticommutation
        m1 = Monomial{PauliAlgebra}([1])  # σ_1
        m2 = Monomial{PauliAlgebra}([1])  # σ_1

        p1 = Polynomial([Term(2.0, m1)])
        p2 = Polynomial([Term(3.0, m2)])

        # σ_1 * σ_1 = I (identity)
        p_prod = p1 * p2
        @test length(terms(p_prod)) == 1
        @test coefficients(p_prod)[1] == 6.0  # 2 * 3
        @test isone(monomials(p_prod)[1])  # Identity monomial
    end

    @testset "Polynomial Multiplication with zero polynomial" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        p = Polynomial([Term(2.0, m)])
        p_zero = zero(Polynomial{NonCommutativeAlgebra,Int64,Float64})

        # p * 0 = 0
        @test iszero(p * p_zero)
        # 0 * p = 0
        @test iszero(p_zero * p)
    end

    @testset "Scalar * Monomial (returns Polynomial)" begin
        m = Monomial{PauliAlgebra}([1, 2])

        # Scalar * monomial
        p = 3.0 * m
        @test p isa Polynomial
        @test length(terms(p)) == 1
        @test coefficients(p)[1] == 3.0
        @test monomials(p)[1] == m

        # Monomial * scalar
        p2 = m * 3.0
        @test p2 == p

        # Zero scalar * monomial
        p_zero = 0.0 * m
        @test iszero(p_zero)
    end

    @testset "Polynomial from zero Term" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        t_zero = Term(0.0, m)

        p = Polynomial(t_zero)
        @test iszero(p)
        @test isempty(terms(p))
    end

    @testset "maxdegree alias" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3, 4])
        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        @test maxdegree(p) == degree(p)
        @test maxdegree(p) == 3
    end
end
