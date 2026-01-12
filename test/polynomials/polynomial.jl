# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: variable_indices

@testset "Polynomial" begin
    @testset "Creation from Terms" begin
        # NonCommutativeAlgebra requires Unsigned integer types
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])  # degree 2
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[3])     # degree 1
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])     # degree 1

        t1 = (1.0, m1)
        t2 = (2.0, m2)
        t3 = (-3.0, m3)

        p = Polynomial([t1, t2, t3])

        # Polynomial should be sorted by monomial ordering: degree-first, then lexicographic
        # Order: [1] (deg 1), [3] (deg 1), [1,2] (deg 2)
        # But [1] < [3] lexicographically, so: [1], [3], [1,2] -> coeffs: -3.0, 2.0, 1.0
        @test length(terms(p)) == 3
        @test coefficients(p) == [-3.0, 2.0, 1.0]
    end

    @testset "Automatic Deduplication" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])  # same as m1

        p = Polynomial([(1.0, m1), (2.0, m2)])

        # Should combine duplicate monomials
        @test length(terms(p)) == 1
        @test coefficients(p) == [3.0]
    end

    @testset "Zero Coefficient Removal" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[2])

        # Zero coefficient term should be removed
        p = Polynomial([(0.0, m1), (1.0, m2)])
        @test length(terms(p)) == 1
        @test coefficients(p) == [1.0]

        # Cancellation should remove the term entirely
        p2 = Polynomial([(2.0, m1), (-2.0, m1)])
        @test iszero(p2)
    end

    @testset "Constant Polynomial" begin
        # NonCommutativeAlgebra requires Unsigned integer types
        p_const = Polynomial{NonCommutativeAlgebra,UInt,Float64}(5.0)

        @test length(terms(p_const)) == 1
        @test coefficients(p_const) == [5.0]
        @test isone(monomials(p_const)[1])

        # Zero constant
        p_zero = Polynomial{NonCommutativeAlgebra,UInt,Float64}(0.0)
        @test iszero(p_zero)
    end

    @testset "From Monomial" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2, 3])
        p = Polynomial(m)

        @test length(terms(p)) == 1
        @test coefficients(p) == [1.0]
        @test monomials(p)[1] == m
    end

    @testset "From Monomial - Preserves Normal Form" begin
        # Polynomial(m) wraps a NormalMonomial directly without additional simplification
        # NormalMonomial constructor enforces normal form, so the monomial is already canonical
        σx = NormalMonomial{PauliAlgebra,Int}([1])  # Already in normal form
        p = Polynomial(σx)

        # Should preserve the monomial as-is
        @test length(terms(p)) == 1
        @test degree(p) == 1

        # For comparison, multiplication applies algebra rules
        p_squared = p * p
        # After simplification, σx₁ * σx₁ = I (identity)
        @test degree(p_squared) == 0  # Simplified to identity
    end

    @testset "From (coefficient, monomial) pair" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        p = Polynomial((3.0, m))

        @test length(terms(p)) == 1
        @test coefficients(p) == [3.0]
    end

    @testset "Accessors" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[2, 3])

        p = Polynomial([(1.0, m1), (2.0, m2)])

        @test length(coefficients(p)) == 2
        @test length(monomials(p)) == 2
    end

    @testset "Degree" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[2, 3, 4])

        p = Polynomial([(1.0, m1), (2.0, m2)])

        @test degree(p) == 3  # max degree

        # Zero polynomial (NonCommutativeAlgebra requires Unsigned)
        p_zero = zero(Polynomial{NonCommutativeAlgebra,UInt,Float64})
        @test degree(p_zero) == -Inf  # preserves deg(p*q) = deg(p) + deg(q)
    end

    @testset "variable_indices" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[2, 3])

        p = Polynomial([(1.0, m1), (1.0, m2)])

        var_idxs = variable_indices(p)
        # variable_indices() returns Set{T} of indices
        @test length(var_idxs) == 3
        @test UInt(1) in var_idxs
        @test UInt(2) in var_idxs
        @test UInt(3) in var_idxs
    end

    @testset "Zero and One" begin
        # NonCommutativeAlgebra requires Unsigned integer types
        T = Polynomial{NonCommutativeAlgebra,UInt,Float64}

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
        m = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])  # degree 2
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[3])    # degree 1
        p = Polynomial([(2.0, m), (-3.0, m2)])

        p_neg = -p
        # Order: [3] (deg 1) before [1,2] (deg 2), so coeffs: 3.0, -2.0
        @test coefficients(p_neg) == [3.0, -2.0]
        @test monomials(p_neg) == monomials(p)
    end

    @testset "Polynomial Addition" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[2])
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[3])

        p1 = Polynomial([(1.0, m1), (2.0, m2)])
        p2 = Polynomial([(3.0, m2), (4.0, m3)])

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
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[1, 2])
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[3])

        p = Polynomial([(2.0, m1), (3.0, m2)])

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
        x = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[5])
        expr_result = x + x^2
        @test expr_result isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(expr_result)) == 2
        # Check the monomials are correct
        monos = monomials(expr_result)
        degrees = [degree(m) for m in monos]
        @test degrees == [1,2]
    end

    @testset "Polynomial Subtraction" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[2])

        p1 = Polynomial([(3.0, m1)])
        p2 = Polynomial([(1.0, m1), (2.0, m2)])

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
        m = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        p = Polynomial([(2.0, m)])

        p_scaled = 3.0 * p
        @test coefficients(p_scaled) == [6.0]

        p_scaled_right = p * 3.0
        @test coefficients(p_scaled_right) == [6.0]
    end

    @testset "Scalar Division" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1])
        p = Polynomial([(6.0, m)])

        p_div = p / 2.0
        @test coefficients(p_div) == [3.0]
    end

    @testset "Polynomial Power" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(UInt16[1])
        p = Polynomial((2.0, m))

        p0 = p^0
        @test isone(p0)

        p1 = p^1
        @test p1 == p

        p2 = p^2
        @test degree(p2) == 2
    end

    @testset "Equality" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[3])

        p1 = Polynomial([(1.0, m1), (2.0, m2)])
        p2 = Polynomial([(1.0, m1), (2.0, m2)])
        p3 = Polynomial([(1.0, m1), (3.0, m2)])

        @test p1 == p2
        @test p1 != p3
    end

    @testset "Hash" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[3])

        p1 = Polynomial([(1.0, m1), (2.0, m2)])
        p2 = Polynomial([(1.0, m1), (2.0, m2)])
        p3 = Polynomial([(1.0, m1), (3.0, m2)])

        @test hash(p1) == hash(p2)
        @test hash(p1) != hash(p3)

        # Unique should work correctly
        p_set = unique([p1, p2, p1, p3])
        @test length(p_set) == 2
    end

    @testset "Copy" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 2])
        p = Polynomial([(2.0, m)])

        p_copy = copy(p)
        @test p == p_copy
        @test p !== p_copy
    end

    @testset "Type Promotion in Arithmetic" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))

        p_int = Polynomial([(2, m)])
        p_float = Polynomial([(1.0, m)])

        p_sum = p_int + p_float
        @test eltype(coefficients(p_sum)) == Float64
        @test coefficients(p_sum) == [3.0]
    end

    @testset "Display (show)" begin
        # Zero polynomial (NonCommutativeAlgebra requires Unsigned)
        p_zero = zero(Polynomial{NonCommutativeAlgebra,UInt16,Float64})
        @test repr(p_zero) == "0"

        # Single term
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p_single = Polynomial([(2.0, m)])
        output = repr(p_single)
        @test contains(output, "2.0")

        # Multiple terms
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))
        p_multi = Polynomial([(1.0, m1), (3.0, m2)])
        output_multi = repr(p_multi)
        # Tests that both terms appear and + separator exists
        @test contains(output_multi, "+")
        @test contains(output_multi, "0x0011")
        @test contains(output_multi, "0x0021")
    end

    @testset "variable_indices" begin
        # Multiple variables
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3))
        p = Polynomial([(1.0, m1), (1.0, m2)])

        var_set = variable_indices(p)
        @test var_set == Set([nc_idx(1), nc_idx(2), nc_idx(3)])

        # Empty polynomial
        p_zero = zero(Polynomial{NonCommutativeAlgebra,UInt16,Float64})
        @test isempty(variable_indices(p_zero))

        # Single variable repeated
        m_single = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 1, 1))
        p_single = Polynomial([(2.0, m_single)])
        @test variable_indices(p_single) == Set([nc_idx(1)])

        # Test with signed types (fermionic/bosonic)
        # Note: negative = creation (a†), positive = annihilation (a)
        # Use normal-ordered: creators first (sorted by mode descending), then annihilators (sorted by mode ascending)
        m_signed = NormalMonomial{FermionicAlgebra,Int}([-3, -1, 2])  # a₃†, a₁†, a₂ (modes 3,1 then 2)
        p_signed = Polynomial([(1.0, m_signed)])
        var_set_signed = variable_indices(p_signed)
        @test var_set_signed == Set([1, 2, 3])  # abs() normalizes for sparsity analysis
    end

    @testset "Error Paths" begin
        # Division by zero
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        p = Polynomial([(6.0, m)])
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
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1, 2))
        p = Polynomial([(5.0, m)])
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
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        p = Polynomial([(2.0, m)])

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
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p = Polynomial([(2.0, m)])

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
        m_pauli = NormalMonomial{PauliAlgebra,Int}([1])
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))

        p_pauli = Polynomial([(1.0, m_pauli)])
        p_nc = Polynomial([(1.0, m_nc)])

        # Different algebra types should never be equal
        @test p_pauli != p_nc
    end

    @testset "Polynomial Multiplication (PauliAlgebra with simplification)" begin
        # PauliAlgebra: σ_i^2 = 1, anticommutation
        m1 = NormalMonomial{PauliAlgebra,Int}([1])  # σ_1
        m2 = NormalMonomial{PauliAlgebra,Int}([1])  # σ_1

        p1 = Polynomial([(2.0, m1)])
        p2 = Polynomial([(3.0, m2)])

        # σ_1 * σ_1 = I (identity)
        p_prod = p1 * p2
        @test length(terms(p_prod)) == 1
        @test coefficients(p_prod)[1] == 6.0  # 2 * 3
        @test isone(monomials(p_prod)[1])  # Identity monomial
    end

    @testset "Polynomial Multiplication with zero polynomial" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p = Polynomial([(2.0, m)])
        p_zero = zero(typeof(p))

        # p * 0 = 0
        @test iszero(p * p_zero)
        # 0 * p = 0
        @test iszero(p_zero * p)
    end

    @testset "Scalar * Monomial (returns Polynomial)" begin
        # Use indices on different sites: 1,4 = σx on sites 1,2
        m = NormalMonomial{PauliAlgebra,Int}([1, 4])

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

    @testset "Polynomial from zero coefficient pair" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p = Polynomial((0.0, m))
        @test iszero(p)
        @test isempty(terms(p))
    end

    @testset "Polynomial(m::NormalMonomial) uses coeff_type" begin
        # PauliAlgebra should create ComplexF64 coefficients
        # Use indices on different sites: 1,4 = σx on sites 1,2
        m_pauli = NormalMonomial{PauliAlgebra,Int}([1, 4])
        p_pauli = Polynomial(m_pauli)
        @test eltype(coefficients(p_pauli)) == ComplexF64
        @test coefficients(p_pauli)[1] == 1.0 + 0.0im

        # NonCommutativeAlgebra should create Float64 coefficients
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        p_nc = Polynomial(m_nc)
        @test eltype(coefficients(p_nc)) == Float64
        @test coefficients(p_nc)[1] == 1.0

        # FermionicAlgebra should create Float64 coefficients
        # Use normal-ordered: creators first (negative), then annihilators (positive)
        m_ferm = NormalMonomial{FermionicAlgebra,Int32}(Int32[-2, 1])
        p_ferm = Polynomial(m_ferm)
        @test eltype(coefficients(p_ferm)) == Float64
    end

    @testset "Monomial * Polynomial multiplication" begin
        # m * p
        # Use different sites so site-based canonicalization makes multiplication commute.
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1; site=1))
        p = Polynomial([
            (1.0, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2; site=2))),
            (2.0, NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3; site=2))),
        ])

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
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        p1 = Polynomial((2.0, m1))
        p2 = Polynomial((3.0, m2))

        p_prod = p1 * p2
        @test coefficients(p_prod)[1] == 6.0
        @test degree(p_prod) == 2
    end

    @testset "Iteration Protocol" begin
        # Multi-term polynomial iteration
        m1 = NormalMonomial{PauliAlgebra,Int}([1])
        m2 = NormalMonomial{PauliAlgebra,Int}([2])
        p = Polynomial([(1.0 + 0im, m1), (2.0 + 0im, m2)])

        pairs = collect(p)
        @test length(pairs) == 2
        @test pairs[1] == (1.0 + 0im, m1)
        @test pairs[2] == (2.0 + 0im, m2)

        # Iteration with for loop
        coeffs = Float64[]
        monos = NormalMonomial{PauliAlgebra,Int}[]
        for (coef, mono) in p
            push!(coeffs, real(coef))
            push!(monos, mono)
        end
        @test coeffs == [1.0, 2.0]
        @test monos == [m1, m2]

        # Empty polynomial
        p_empty = zero(Polynomial{PauliAlgebra,Int64,ComplexF64})
        @test isempty(collect(p_empty))
        @test iterate(p_empty) === nothing

        # Single-term polynomial
        p_single = Polynomial([(3.5 + 0im, m1)])
        pairs_single = collect(p_single)
        @test length(pairs_single) == 1
        @test pairs_single[1] == (3.5 + 0im, m1)

        # eltype
        @test eltype(typeof(p)) == Tuple{ComplexF64,NormalMonomial{PauliAlgebra,Int64}}

        # NonCommutative algebra (Float64 coefficients)
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        p_nc = Polynomial([(5.0, m_nc)])
        pairs_nc = collect(p_nc)
        @test pairs_nc[1][1] isa Float64
        @test eltype(typeof(p_nc)) == Tuple{Float64,NormalMonomial{NonCommutativeAlgebra,UInt16}}
    end

    @testset "Polynomial - Monomial subtraction" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 2))
        p = Polynomial([(3.0, m1)])

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

    @testset "isless(Polynomial, Polynomial) - Degree-First Ordering" begin
        # isless compares polynomials for sorting:
        # 1. Compare by highest-degree monomial first
        # 2. For equal monomials, compare by coefficient magnitude
        # 3. Fewer terms is "less" when all compared terms are equal

        # Test 1: Different highest-degree monomials
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        p_one_term = Polynomial([(1.0, m1)])
        p_two_terms = Polynomial([(1.0, m1), (2.0, m2)])
        p_three_terms = Polynomial([(1.0, m1), (2.0, m2), (3.0, m3)])

        @test isless(p_one_term, p_two_terms)
        @test isless(p_two_terms, p_three_terms)
        @test isless(p_one_term, p_three_terms)

        # Test 2: Same number of terms - compare by monomial (lexicographic)
        # Monomials: [1] < [2] < [3] (by lexicographic order)
        p_m1 = Polynomial([(1.0, m1)])
        p_m2 = Polynomial([(1.0, m2)])
        p_m3 = Polynomial([(1.0, m3)])

        @test isless(p_m1, p_m2)
        @test isless(p_m2, p_m3)
        @test isless(p_m1, p_m3)

        # Test 3: Same first term, different second term
        p_12 = Polynomial([(1.0, m1), (1.0, m2)])
        p_13 = Polynomial([(1.0, m1), (1.0, m3)])
        @test isless(p_12, p_13)  # [1,2] < [1,3] because m2 < m3

        # Test 4: Equal polynomials - isless returns false
        p_same1 = Polynomial([(2.0, m1)])
        p_same2 = Polynomial([(2.0, m1)])
        @test !isless(p_same1, p_same2)
        @test !isless(p_same2, p_same1)

        # Test 5: Zero polynomial (no terms) should come first
        p_zero = zero(typeof(p_one_term))
        @test isless(p_zero, p_one_term)
        @test !isless(p_one_term, p_zero)

        # Test 6: Graded ordering - degree-first for monomials affects polynomial ordering
        m_deg1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))       # degree 1
        m_deg2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))    # degree 2
        p_deg1 = Polynomial([(1.0, m_deg1)])
        p_deg2 = Polynomial([(1.0, m_deg2)])
        @test isless(p_deg1, p_deg2)  # Monomial with lower degree < higher degree

        # Test 7: Sorting a vector of polynomials
        polynomials = [p_three_terms, p_one_term, p_two_terms]
        sorted_polys = sort(polynomials)
        @test sorted_polys[1] == p_one_term
        @test sorted_polys[2] == p_two_terms
        @test sorted_polys[3] == p_three_terms
    end

    @testset "convert(Polynomial{A,T,C2}, p) - Coefficient Type Conversion" begin
        # Test 1: Int to Float64
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        p_int = Polynomial([(2, m)])
        p_float = convert(Polynomial{NonCommutativeAlgebra,UInt8,Float64}, p_int)

        @test p_float isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test coefficients(p_float)[1] === 2.0
        @test monomials(p_float)[1] == m

        # Test 2: Float64 to ComplexF64
        p_real = Polynomial([(3.0, m)])
        p_complex = convert(Polynomial{NonCommutativeAlgebra,UInt8,ComplexF64}, p_real)

        @test p_complex isa Polynomial{NonCommutativeAlgebra,UInt8,ComplexF64}
        @test coefficients(p_complex)[1] === 3.0 + 0.0im

        # Test 3: Identity conversion (no-op)
        m_u8 = NormalMonomial{NonCommutativeAlgebra,UInt8}(nc_word(UInt8, 1))
        p_float64 = Polynomial([(5.0, m_u8)])
        p_same = convert(Polynomial{NonCommutativeAlgebra,UInt8,Float64}, p_float64)
        @test p_same === p_float64  # Should return the same object

        # Test 4: Multiple terms conversion
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 2))
        p_multi_int = Polynomial([(2, m1), (3, m2)])
        p_multi_float = convert(Polynomial{NonCommutativeAlgebra,UInt8,Float64}, p_multi_int)

        @test length(terms(p_multi_float)) == 2
        @test all(c isa Float64 for c in coefficients(p_multi_float))
        @test coefficients(p_multi_float) == [2.0, 3.0]

        # Test 5: Zero polynomial conversion
        p_zero_int = zero(Polynomial{NonCommutativeAlgebra,UInt8,Int})
        p_zero_float = convert(Polynomial{NonCommutativeAlgebra,UInt8,Float64}, p_zero_int)
        @test iszero(p_zero_float)
        @test p_zero_float isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}

        # Test 6: Pauli algebra conversion (ComplexF64 to other complex types)
        # Use indices on different sites: 1,4 = σx on sites 1,2
        m_pauli = NormalMonomial{PauliAlgebra,Int}([1, 4])
        p_pauli = Polynomial([(1.0 + 2.0im, m_pauli)])
        # Convert to ComplexF32 (narrower precision)
        p_pauli_f32 = convert(Polynomial{PauliAlgebra,Int64,ComplexF32}, p_pauli)
        @test p_pauli_f32 isa Polynomial{PauliAlgebra,Int64,ComplexF32}
        @test coefficients(p_pauli_f32)[1] ≈ ComplexF32(1.0 + 2.0im)
    end

    @testset "coeff_type - Coefficient Type Accessor" begin
        # Test 1: From type
        @test NCTSSoS.coeff_type(Polynomial{NonCommutativeAlgebra,UInt16,Float64}) == Float64
        @test NCTSSoS.coeff_type(Polynomial{PauliAlgebra,Int64,ComplexF64}) == ComplexF64
        @test NCTSSoS.coeff_type(Polynomial{FermionicAlgebra,Int32,Float32}) == Float32

        # Test 2: From instance
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        p_float = Polynomial([(2.0, m)])
        @test NCTSSoS.coeff_type(p_float) == Float64

        m_pauli = NormalMonomial{PauliAlgebra,Int}([1])
        p_complex = Polynomial([(1.0 + 0.0im, m_pauli)])
        @test NCTSSoS.coeff_type(p_complex) == ComplexF64

        # Test 3: Zero polynomial preserves type
        p_zero = zero(Polynomial{NonCommutativeAlgebra,UInt16,Float64})
        @test NCTSSoS.coeff_type(p_zero) == Float64

        # Test 4: Integer coefficients
        m_int = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word( 1))
        p_int = Polynomial([(5, m_int)])
        @test NCTSSoS.coeff_type(p_int) == Int

        # Test 5: Type consistency through operations
        p1 = Polynomial([(1.0, m)])
        p2 = Polynomial([(2.0, m)])
        p_sum = p1 + p2
        @test NCTSSoS.coeff_type(p_sum) == Float64

        # Test 6: Type promotion through operations
        p_int2 = Polynomial([(3, m)])
        p_promoted = p_float + p_int2
        @test NCTSSoS.coeff_type(p_promoted) == Float64  # Int + Float64 → Float64
    end
end
