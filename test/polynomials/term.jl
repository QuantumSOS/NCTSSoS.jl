# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: ComposedMonomial

@testset "Term" begin
    @testset "Construction" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t = Term(2.0, m)

        @test t.coefficient == 2.0
        @test t.monomial == m
        @test t isa Term{Monomial{NonCommutativeAlgebra,Int64},Float64}
    end

    @testset "isone" begin
        m_empty = Monomial{NonCommutativeAlgebra}(Int[])
        t_one = Term(1.0, m_empty)
        @test isone(t_one)

        # Coefficient not 1 -> false
        t_two = Term(2.0, m_empty)
        @test !isone(t_two)

        # Monomial not empty -> false
        m_non_empty = Monomial{NonCommutativeAlgebra}([1])
        t_one_mono = Term(1.0, m_non_empty)
        @test !isone(t_one_mono)

        # Complex coefficient
        t_one_c = Term(1.0 + 0.0im, Monomial{PauliAlgebra}(Int[]))
        @test isone(t_one_c)

        # Complex with imaginary part -> false
        t_not_one = Term(1.0 + 0.1im, Monomial{PauliAlgebra}(Int[]))
        @test !isone(t_not_one)
    end

    @testset "iszero" begin
        m_empty = Monomial{NonCommutativeAlgebra}(Int[])
        m_non_empty = Monomial{NonCommutativeAlgebra}([1])

        # Zero coefficient -> zero
        t_zero = Term(0.0, m_empty)
        @test iszero(t_zero)

        # Zero with non-empty monomial still zero
        t_zero_mono = Term(0.0, m_non_empty)
        @test iszero(t_zero_mono)

        # Non-zero coefficient -> not zero
        t_nonzero = Term(1.0, m_empty)
        @test !iszero(t_nonzero)

        t_nonzero2 = Term(0.001, m_non_empty)
        @test !iszero(t_nonzero2)

        # Complex zero
        t_zero_c = Term(0.0 + 0.0im, Monomial{PauliAlgebra}(Int[]))
        @test iszero(t_zero_c)
    end

    @testset "Equality and Hash" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2])
        m3 = Monomial{NonCommutativeAlgebra}([2, 1])

        t1 = Term(2.0, m1)
        t2 = Term(2.0, m2)
        t3 = Term(2.0, m3)
        t4 = Term(3.0, m1)

        # Same coefficient and monomial
        @test t1 == t2

        # Different monomial
        @test t1 != t3

        # Different coefficient
        @test t1 != t4

        # Cross-algebra comparison should return false
        m_pauli = Monomial{PauliAlgebra}([1, 2])
        t_pauli = Term(2.0, m_pauli)
        @test t1 != t_pauli

        # Different coefficient types but same values
        t_int = Term(2, Monomial{NonCommutativeAlgebra}([1, 2]))
        @test t1 == t_int  # 2.0 == 2

        # Hash contract: equal terms must have equal hashes
        @test hash(t1) == hash(t2)

        # Hash works in collections
        d = Dict{typeof(t1),Int}()
        d[t1] = 1
        @test d[t2] == 1  # t2 == t1, so should find it

        # Set deduplication
        s = Set([t1, t2])
        @test length(s) == 1

        # Different terms have different hashes (with high probability)
        @test hash(t1) != hash(t3)
        @test hash(t1) != hash(t4)
    end

    @testset "one(Type)" begin
        T = Term{Monomial{UnipotentAlgebra,Int},Float64}
        t_one = one(T)
        @test isone(t_one)
        @test t_one.coefficient == 1.0
        @test isempty(t_one.monomial.word)

        # Different coefficient types
        T_complex = Term{Monomial{PauliAlgebra,UInt16},ComplexF64}
        t_one_c = one(T_complex)
        @test isone(t_one_c)
        @test t_one_c.coefficient == (1.0 + 0.0im)
        @test isempty(t_one_c.monomial.word)

        # Different integer types
        T_int32 = Term{Monomial{FermionicAlgebra,Int32},Float64}
        t_one_32 = one(T_int32)
        @test isone(t_one_32)
        @test eltype(t_one_32.monomial.word) == Int32
    end

    @testset "zero(Type)" begin
        T = Term{Monomial{UnipotentAlgebra,Int},Float64}
        t_zero = zero(T)
        @test iszero(t_zero)
        @test t_zero.coefficient == 0.0
        @test isempty(t_zero.monomial.word)

        # Different coefficient types
        T_complex = Term{Monomial{PauliAlgebra,UInt16},ComplexF64}
        t_zero_c = zero(T_complex)
        @test iszero(t_zero_c)
        @test t_zero_c.coefficient == (0.0 + 0.0im)

        # Verify isone(zero(T)) is false
        @test !isone(zero(T))
    end

    @testset "one/zero for ComposedMonomial" begin
        # Create a ComposedMonomial type
        m_pauli = Monomial{PauliAlgebra,UInt16}(UInt16[1, 2])
        m_fermi = Monomial{FermionicAlgebra,Int32}(Int32[1])
        cm = ComposedMonomial((m_pauli, m_fermi))

        # Test one for ComposedMonomial
        cm_one = one(cm)
        @test isone(cm_one)
        @test all(isone, cm_one.components)

        # Test one for Term with ComposedMonomial
        CMType = typeof(cm)
        T_composed = Term{CMType,ComplexF64}
        t_one = one(T_composed)
        @test isone(t_one)
        @test t_one.coefficient == (1.0 + 0.0im)

        # Test zero for Term with ComposedMonomial
        t_zero = zero(T_composed)
        @test iszero(t_zero)
        @test t_zero.coefficient == (0.0 + 0.0im)
    end

    @testset "Scalar Multiplication" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t = Term(2.0, m)

        # Left multiplication
        t_left = 3.0 * t
        @test t_left.coefficient == 6.0
        @test t_left.monomial == m

        # Right multiplication
        t_right = t * 3.0
        @test t_right.coefficient == 6.0
        @test t_right.monomial == m

        # Commutativity
        @test (3.0 * t) == (t * 3.0)

        # Type promotion: Int * Term{..., Float64} -> Float64
        t_int_scaled = 2 * t
        @test t_int_scaled.coefficient == 4.0
        @test t_int_scaled isa Term{Monomial{NonCommutativeAlgebra,Int64},Float64}

        # Type promotion: ComplexF64 * Term{..., Float64} -> ComplexF64
        t_complex_scaled = (1.0 + 0.0im) * t
        @test t_complex_scaled.coefficient == (2.0 + 0.0im)
        @test t_complex_scaled isa Term{Monomial{NonCommutativeAlgebra,Int64},ComplexF64}

        # Multiplication by zero
        t_zeroed = 0.0 * t
        @test iszero(t_zeroed)

        # Multiplication by one
        t_same = 1.0 * t
        @test t_same == t
    end

    @testset "Negation" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t = Term(2.0, m)

        t_neg = -t
        @test t_neg.coefficient == -2.0
        @test t_neg.monomial == m

        # Double negation
        @test -(-t) == t

        # Negation of zero
        t_zero = Term(0.0, m)
        @test iszero(-t_zero)

        # Negation preserves type
        @test typeof(-t) == typeof(t)
    end

    @testset "Display (show)" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        m_empty = Monomial{NonCommutativeAlgebra}(Int[])

        # Zero term
        t_zero = Term(0.0, m)
        @test sprint(show, t_zero) == "0"

        # Identity term (coefficient 1, empty monomial)
        t_one = Term(1.0, m_empty)
        @test sprint(show, t_one) == "1"

        # Empty monomial with coefficient != 1
        t_const = Term(3.5, m_empty)
        @test sprint(show, t_const) == "3.5"

        # Coefficient = 1 with non-empty monomial
        t_coef_one = Term(1.0, m)
        @test sprint(show, t_coef_one) == "[1, 2]"

        # Coefficient = -1
        t_coef_neg_one = Term(-1.0, m)
        @test sprint(show, t_coef_neg_one) == "-[1, 2]"

        # General case
        t_general = Term(2.5, m)
        @test sprint(show, t_general) == "2.5 * [1, 2]"

        # Negative coefficient (not -1)
        t_neg = Term(-2.5, m)
        @test sprint(show, t_neg) == "-2.5 * [1, 2]"

        # Complex coefficient (parenthesized)
        t_complex = Term(1.0 + 2.0im, m)
        @test sprint(show, t_complex) == "(1.0 + 2.0im) * [1, 2]"
    end

    @testset "Iteration Protocol" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t = Term(2.5, m)

        # length (now 1, yields single (coef, mono) tuple)
        @test length(t) == 1

        # New destructuring: yields single tuple
        (coef, mono), = t
        @test coef == 2.5
        @test mono == m

        # Direct field access (alternative to destructuring)
        @test t.coefficient == 2.5
        @test t.monomial == m

        # collect yields vector of tuples
        pairs = collect(t)
        @test length(pairs) == 1
        @test pairs[1] == (2.5, m)

        # Manual iteration: yields (coef, mono) tuple, then nothing
        iter = iterate(t)
        @test iter !== nothing
        @test iter[1] == (2.5, m)

        iter2 = iterate(t, iter[2])
        @test iter2 === nothing

        # eltype
        @test eltype(typeof(t)) == Tuple{Float64,Monomial{NonCommutativeAlgebra,Int64}}
    end

    @testset "Immutability" begin
        # Term is immutable - monomial reference is shared
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t1 = Term(2.0, m)
        t2 = Term(3.0, m)

        # Both terms share the same monomial
        @test t1.monomial === t2.monomial

        # Creating a new term via multiplication shares the monomial
        t3 = 2.0 * t1
        @test t3.monomial === t1.monomial
    end

    @testset "coeff_type" begin
        using NCTSSoS: coeff_type

        # Type-level dispatch (primary usage)
        T_float = Term{Monomial{NonCommutativeAlgebra,Int},Float64}
        @test coeff_type(T_float) === Float64

        T_complex = Term{Monomial{PauliAlgebra,UInt16},ComplexF64}
        @test coeff_type(T_complex) === ComplexF64

        T_int = Term{Monomial{UnipotentAlgebra,Int32},Int}
        @test coeff_type(T_int) === Int

        # Instance-level dispatch
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t_float = Term(2.5, m)
        @test coeff_type(t_float) === Float64

        t_complex = Term(1.0 + 2.0im, Monomial{PauliAlgebra}([1]))
        @test coeff_type(t_complex) === ComplexF64

        # After scalar multiplication (type promotion)
        t_promoted = (1.0 + 0.0im) * t_float
        @test coeff_type(t_promoted) === ComplexF64

        # ComposedMonomial support
        m_pauli = Monomial{PauliAlgebra,UInt16}(UInt16[1, 2])
        m_fermi = Monomial{FermionicAlgebra,Int32}(Int32[1])
        cm = ComposedMonomial((m_pauli, m_fermi))
        CMType = typeof(cm)
        T_composed = Term{CMType,ComplexF64}
        @test coeff_type(T_composed) === ComplexF64

        # Instance with ComposedMonomial
        t_composed = Term(3.0 + 1.0im, cm)
        @test coeff_type(t_composed) === ComplexF64
    end
end
