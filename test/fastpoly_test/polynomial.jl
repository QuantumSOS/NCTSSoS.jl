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

    @testset "From Monomial - No Simplification" begin
        # Polynomial(m) does NOT simplify the monomial - this is documented behavior
        # For Pauli: σx * σx should simplify to I, but Polynomial(m) should preserve it
        m_reducible = Monomial{PauliAlgebra}([1, 1])  # σx₁ * σx₁ (same Pauli twice)
        p = Polynomial(m_reducible)

        # Should preserve the un-simplified monomial
        @test length(terms(p)) == 1
        @test degree(p) == 2  # Still degree 2, not 0 (identity)
        @test monomials(p)[1] == m_reducible

        # For comparison, multiplication DOES simplify
        σx = Monomial{PauliAlgebra}([1])
        p_simplified = Polynomial(σx) * Polynomial(σx)
        # After simplification, σx₁ * σx₁ = I (identity)
        @test degree(p_simplified) == 0  # Simplified to identity
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
    end

    @testset "Degree" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3, 4])

        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        @test degree(p) == 3  # max degree

        # Zero polynomial
        p_zero = zero(Polynomial{NonCommutativeAlgebra,Int64,Float64})
        @test degree(p_zero) == -Inf  # preserves deg(p*q) = deg(p) + deg(q)
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

    @testset "Polynomial + Monomial Addition" begin
        # Setup: create a polynomial and a different monomial
        m1 = Monomial{NonCommutativeAlgebra}(UInt8[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2])
        m3 = Monomial{NonCommutativeAlgebra}(UInt8[3])

        p = Polynomial([Term(2.0, m1), Term(3.0, m2)])

        # Test Polynomial + Monomial
        result1 = p + m3
        @test result1 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(result1)) == 3  # Original 2 terms + new monomial

        # Test Monomial + Polynomial
        result2 = m3 + p
        @test result2 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test result2 == result1  # Should be commutative

        # Test reduction: summing vector of monomials
        monomials_vec = [m1, m2, m3]
        sum_result = sum(monomials_vec)
        @test sum_result isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(sum_result)) == 3

        # Test with powers: x + x^2 pattern that was failing
        x = Monomial{NonCommutativeAlgebra}(UInt8[5])
        expr_result = x + x^2
        @test expr_result isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(expr_result)) == 2
        # Check the monomials are correct
        monos = monomials(expr_result)
        degrees = [degree(m) for m in monos]
        @test 1 in degrees
        @test 2 in degrees
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

    @testset "Adjoint with Julia syntax" begin
        m = Monomial{PauliAlgebra}([1, 2])
        p = Polynomial([Term(1.0 + 1.0im, m)])

        @test p' == adjoint(p)
        @test coefficients(p')[1] == 1.0 - 1.0im
    end

    @testset "Adjoint for Fermionic Polynomials" begin
        # Fermionic adjoint: reverse order AND negate indices
        # For [1, -2] (a₁ a₂†), adjoint is a₂ a₁† = [2, -1] after reversing [-(-2), -1]
        m = Monomial{FermionicAlgebra}(Int32[1, -2])
        p = Polynomial([Term(2.0, m)])

        p_adj = adjoint(p)

        # Real coefficient stays unchanged
        @test coefficients(p_adj)[1] == 2.0
        # Monomial adjoint: reverse and negate: [1, -2] -> [2, -1]
        @test monomials(p_adj)[1].word == Int32[2, -1]

        # Test involution: adjoint(adjoint(p)) == p
        @test adjoint(adjoint(p)) == p
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
        # Note: negative = creation (a†), positive = annihilation (a)
        m_signed = Monomial{FermionicAlgebra}([-1, 2, -3])  # a₁†, a₂, a₃†
        p_signed = Polynomial([Term(1.0, m_signed)])
        var_set_signed = variable_indices(p_signed)
        @test var_set_signed == Set([1, 2, 3])  # abs() normalizes for sparsity analysis
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

    @testset "Polynomial(m::Monomial) uses default_coeff_type" begin
        # PauliAlgebra should create ComplexF64 coefficients
        m_pauli = Monomial{PauliAlgebra}([1, 2])
        p_pauli = Polynomial(m_pauli)
        @test eltype(coefficients(p_pauli)) == ComplexF64
        @test coefficients(p_pauli)[1] == 1.0 + 0.0im

        # NonCommutativeAlgebra should create Float64 coefficients
        m_nc = Monomial{NonCommutativeAlgebra}(UInt16[1, 2])
        p_nc = Polynomial(m_nc)
        @test eltype(coefficients(p_nc)) == Float64
        @test coefficients(p_nc)[1] == 1.0

        # FermionicAlgebra should create Float64 coefficients
        m_ferm = Monomial{FermionicAlgebra}(Int32[1, -2])
        p_ferm = Polynomial(m_ferm)
        @test eltype(coefficients(p_ferm)) == Float64
    end

    @testset "simplify(p::Polynomial) - direct tests" begin
        # Test 1: Pauli simplification - σx₁ * σx₁ = I
        m_xx = Monomial{PauliAlgebra}([1, 1])
        p_unsimplified = Polynomial(m_xx)
        @test degree(p_unsimplified) == 2  # Before simplify

        p_simplified = simplify(p_unsimplified)
        @test degree(p_simplified) == 0  # After simplify: identity
        @test isone(monomials(p_simplified)[1])
        @test coefficients(p_simplified)[1] == 1.0 + 0.0im

        # Test 2: Pauli simplification with phase - σx₁ * σy₁ = i*σz₁
        m_xy = Monomial{PauliAlgebra}([1, 2])  # σx₁ * σy₁
        p_xy = Polynomial(m_xy)
        p_xy_simplified = simplify(p_xy)
        @test degree(p_xy_simplified) == 1
        @test coefficients(p_xy_simplified)[1] ≈ 0.0 + 1.0im  # Phase is i

        # Test 3: Zero polynomial stays zero
        p_zero = zero(Polynomial{PauliAlgebra,Int64,ComplexF64})
        @test iszero(simplify(p_zero))

        # Test 4: Multiple terms - each simplified independently
        m1 = Monomial{PauliAlgebra}([1, 1])  # σx₁² = I
        m2 = Monomial{PauliAlgebra}([4])     # σx₂ (already simple)
        p_multi = Polynomial([Term(2.0+0im, m1), Term(3.0+0im, m2)])
        p_multi_simplified = simplify(p_multi)
        # Should have identity term (coef 2) and σx₂ term (coef 3)
        @test length(terms(p_multi_simplified)) == 2

        # Test 5: Coefficient type promotion (Float64 → ComplexF64 for Pauli)
        m_pauli = Monomial{PauliAlgebra}([1])
        p_float = Polynomial([Term(2.0, m_pauli)])  # Float64 coefficient
        p_promoted = simplify(p_float)
        @test eltype(coefficients(p_promoted)) == ComplexF64
    end

    @testset "Monomial * Polynomial multiplication" begin
        # m * p
        m = Monomial{NonCommutativeAlgebra}(UInt16[1])
        p = Polynomial([Term(1.0, Monomial{NonCommutativeAlgebra}(UInt16[2])),
                       Term(2.0, Monomial{NonCommutativeAlgebra}(UInt16[3]))])

        result_mp = m * p
        @test result_mp isa Polynomial{NonCommutativeAlgebra,UInt16,Float64}
        @test length(terms(result_mp)) == 2

        # p * m
        result_pm = p * m
        @test result_pm isa Polynomial{NonCommutativeAlgebra,UInt16,Float64}
        @test length(terms(result_pm)) == 2

        # Note: For NonCommutativeAlgebra with site-based encoding,
        # operators on different sites commute after simplification.
        # Both m*p and p*m should give same result when sites differ.
        @test result_mp == result_pm
    end

    @testset "Polynomial * Polynomial (NonCommutativeAlgebra)" begin
        # Test multiplication with UInt16 (required by simplify for NC)
        m1 = Monomial{NonCommutativeAlgebra}(UInt16[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt16[2])

        p1 = Polynomial(Term(2.0, m1))
        p2 = Polynomial(Term(3.0, m2))

        p_prod = p1 * p2
        @test coefficients(p_prod)[1] == 6.0
        @test degree(p_prod) == 2
    end

    @testset "Polynomial - Monomial subtraction" begin
        m1 = Monomial{NonCommutativeAlgebra}(UInt8[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt8[2])
        p = Polynomial([Term(3.0, m1)])

        # p - m
        result1 = p - m2
        @test length(terms(result1)) == 2
        monos = monomials(result1)
        coeffs = coefficients(result1)
        m2_idx = findfirst(==(m2), monos)
        @test coeffs[m2_idx] == -1.0

        # m - p
        result2 = m2 - p
        @test length(terms(result2)) == 2
        monos2 = monomials(result2)
        coeffs2 = coefficients(result2)
        m1_idx = findfirst(==(m1), monos2)
        m2_idx2 = findfirst(==(m2), monos2)
        @test coeffs2[m1_idx] == -3.0
        @test coeffs2[m2_idx2] == 1.0
    end
end
