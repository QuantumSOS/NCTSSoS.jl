# Note: FastPolynomials is loaded by setup.jl
using Test, NCTSSoS.FastPolynomials

@testset "Monomials" begin
    @testset "Creation" begin
        # Create monomials directly with word representation
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 1, 3])
        @test mono1.word == [1, 2, 1, 3]
        @test degree(mono1) == 4

        # Empty monomial (identity)
        mono_identity = Monomial{NonCommutativeAlgebra}(Int[])
        @test isone(mono_identity)
        @test degree(mono_identity) == 0

        # Single element monomial
        mono_single = Monomial{NonCommutativeAlgebra}([5])
        @test mono_single.word == [5]
        @test degree(mono_single) == 1

        # Monomial with zeros should filter them out
        mono_with_zero = Monomial{NonCommutativeAlgebra}([1, 0, 2, 0, 3])
        @test mono_with_zero.word == [1, 2, 3]
        @test degree(mono_with_zero) == 3

        # one(Monomial) should return identity
        mono_one = one(Monomial{NonCommutativeAlgebra,Int64})
        @test isone(mono_one)
    end

    @testset "Degree" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = Monomial{NonCommutativeAlgebra}(Int[])
        mono3 = Monomial{NonCommutativeAlgebra}([1, 1, 1, 1, 1])

        @test degree(mono1) == 3
        @test degree(mono2) == 0
        @test degree(mono3) == 5
    end

    @testset "Hash" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = Monomial{NonCommutativeAlgebra}([1, 3, 3])
        mono3 = Monomial{NonCommutativeAlgebra}([1, 2, 3])

        @test hash(mono1) != hash(mono2)
        @test hash(mono1) == hash(mono3)

        mono_empty = Monomial{NonCommutativeAlgebra}(Int[])
        mono_zero_filtered = Monomial{NonCommutativeAlgebra}([0])
        @test hash(mono_empty) == hash(mono_zero_filtered)
    end

    @testset "Star Operation" begin
        # Star reverses word for self-adjoint algebras (unsigned types)
        mono1 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        mono1_star = star(mono1)
        @test mono1_star.word == [3, 2, 1]

        # Empty monomial star should be empty
        mono_empty = Monomial{NonCommutativeAlgebra}(UInt8[])
        @test isone(star(mono_empty))

        # Single element star
        mono_single = Monomial{NonCommutativeAlgebra}(UInt8[5])
        @test star(mono_single).word == [5]

        # Star is involution: star(star(m)) == m
        mono2 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3, 4])
        @test star(star(mono2)) == mono2
    end

    @testset "Multiplication" begin
        # Note: Multiplication now returns Monomial (word concatenation only)
        # Callers should apply simplify! explicitly if needed
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Same variable multiplication produces longer word (concatenation)
        result = x[1] * x[1]
        @test result isa Monomial
        @test degree(result) == 2

        # Different variables
        result2 = x[1] * x[2]
        @test degree(result2) == 2

        # Identity multiplication
        mono_id = one(Monomial{NonCommutativeAlgebra,UInt8})
        mono_x = Monomial{NonCommutativeAlgebra}(UInt8[1])
        result3 = mono_id * mono_x
        @test result3 == mono_x

        result4 = mono_x * mono_id
        @test result4 == mono_x
    end

    @testset "Comparison" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1])
        mono2 = Monomial{NonCommutativeAlgebra}([1, 2])
        mono3 = Monomial{NonCommutativeAlgebra}([2])

        # Degree-first ordering
        @test isless(mono1, mono2)  # degree 1 < degree 2

        # Same degree: lexicographic
        @test isless(mono1, mono3)  # [1] < [2]

        # Sorting works
        monos = [mono2, mono3, mono1]
        sorted_monos = sort(monos)
        @test sorted_monos == [mono1, mono3, mono2]
    end

    @testset "Equality" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono3 = Monomial{NonCommutativeAlgebra}([1, 3, 2])

        @test mono1 == mono2
        @test mono1 != mono3

        # Different algebra types are never equal
        mono_pauli = Monomial{PauliAlgebra}([1, 2, 3])
        mono_nc = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        @test mono_pauli != mono_nc
    end

    @testset "Algebra Type Preservation" begin
        # Monomials preserve their algebra type
        mono_pauli = Monomial{PauliAlgebra}([1, 2])
        mono_fermi = Monomial{FermionicAlgebra}(Int32[1, 2])
        mono_unipotent = Monomial{UnipotentAlgebra}([1, 2])

        @test mono_pauli isa Monomial{PauliAlgebra}
        @test mono_fermi isa Monomial{FermionicAlgebra}
        @test mono_unipotent isa Monomial{UnipotentAlgebra}
    end

    @testset "adjoint! for Signed Types" begin
        # Test signed type negation (critical for Fermionic/Bosonic algebras)
        word_signed = Int32[1, -2, 3]
        adjoint!(word_signed)
        @test word_signed == Int32[-3, 2, -1]

        # Test with all positive indices
        word_pos = Int64[1, 2, 3]
        adjoint!(word_pos)
        @test word_pos == Int64[-3, -2, -1]

        # Test with all negative indices
        word_neg = Int32[-1, -2, -3]
        adjoint!(word_neg)
        @test word_neg == Int32[3, 2, 1]

        # Test empty word
        word_empty = Int32[]
        adjoint!(word_empty)
        @test word_empty == Int32[]

        # Test single element
        word_single = Int64[5]
        adjoint!(word_single)
        @test word_single == Int64[-5]

        # Test unsigned types (should only reverse, not negate)
        word_unsigned = UInt16[1, 2, 3]
        adjoint!(word_unsigned)
        @test word_unsigned == UInt16[3, 2, 1]
    end

    @testset "star! Direct Tests" begin
        # Verify star! is correctly aliased to adjoint!
        word_unsigned = UInt16[1, 2, 3]
        star!(word_unsigned)
        @test word_unsigned == UInt16[3, 2, 1]

        # Test star! for signed types
        word_signed = Int32[1, -2, 3]
        star!(word_signed)
        @test word_signed == Int32[-3, 2, -1]

        # Verify star! and adjoint! produce identical results
        word1 = Int64[1, 2, 3]
        word2 = Int64[1, 2, 3]
        adjoint!(word1)
        star!(word2)
        @test word1 == word2
    end

    @testset "Adjoint for Fermionic/Signed Types" begin
        # Test FermionicAlgebra adjoint (reverses AND negates)
        m_ferm = Monomial{FermionicAlgebra}(Int32[1, -2, 3])
        m_adj = adjoint(m_ferm)
        @test m_adj.word == Int32[-3, 2, -1]

        # Verify involution property for signed types: adjoint(adjoint(m)) == m
        m_ferm2 = Monomial{FermionicAlgebra}(Int32[1, -2, 3, -4])
        @test adjoint(adjoint(m_ferm2)) == m_ferm2

        # Empty fermionic monomial
        m_empty = Monomial{FermionicAlgebra}(Int32[])
        @test isone(adjoint(m_empty))

        # Single element
        m_single = Monomial{FermionicAlgebra}(Int32[5])
        @test adjoint(m_single).word == Int32[-5]

        # star() should also work for signed types
        @test star(m_ferm) == m_adj
    end

    @testset "one(m::Monomial) Instance Method" begin
        # Test instance method preserves algebra type
        m_pauli = Monomial{PauliAlgebra}(UInt16[1, 2, 3])
        @test one(m_pauli) == Monomial{PauliAlgebra}(UInt16[])
        @test typeof(one(m_pauli)) == Monomial{PauliAlgebra,UInt16}

        # Test with FermionicAlgebra
        m_ferm = Monomial{FermionicAlgebra}(Int32[1, -2])
        @test one(m_ferm) == Monomial{FermionicAlgebra}(Int32[])
        @test typeof(one(m_ferm)) == Monomial{FermionicAlgebra,Int32}

        # Test that one(m) returns identity
        @test isone(one(m_pauli))
        @test isone(one(m_ferm))

        # Verify one(type) and one(instance) return same result
        @test one(Monomial{PauliAlgebra,UInt16}) == one(m_pauli)
    end

    @testset "Monomial Default Constructor" begin
        # Test default constructor uses NonCommutativeAlgebra
        m = Monomial([1, 2, 3])
        @test m isa Monomial{NonCommutativeAlgebra}
        @test m.word == [1, 2, 3]

        # Test with different integer types
        m_int32 = Monomial(Int32[1, 2])
        @test m_int32 isa Monomial{NonCommutativeAlgebra,Int32}

        m_uint16 = Monomial(UInt16[1, 2])
        @test m_uint16 isa Monomial{NonCommutativeAlgebra,UInt16}

        # Test zero filtering in default constructor
        m_zeros = Monomial([1, 0, 2, 0])
        @test m_zeros.word == [1, 2]

        # Empty default constructor
        m_empty = Monomial(Int[])
        @test isone(m_empty)
    end

    @testset "Zero Filtering Edge Cases" begin
        # Leading zeros
        m1 = Monomial{NonCommutativeAlgebra}([0, 0, 1, 2])
        @test m1.word == [1, 2]

        # Trailing zeros
        m2 = Monomial{NonCommutativeAlgebra}([1, 2, 0, 0])
        @test m2.word == [1, 2]

        # All zeros -> identity
        m3 = Monomial{NonCommutativeAlgebra}([0, 0, 0])
        @test isone(m3)

        # Mixed zeros throughout
        m4 = Monomial{NonCommutativeAlgebra}([0, 1, 0, 2, 0, 3, 0])
        @test m4.word == [1, 2, 3]

        # Single zero
        m5 = Monomial{NonCommutativeAlgebra}([0])
        @test isone(m5)
    end

    @testset "Cross-Algebra Type Equality" begin
        # Same word, different algebras should never be equal
        word = [1, 2, 3]
        m_pauli = Monomial{PauliAlgebra}(word)
        m_unipotent = Monomial{UnipotentAlgebra}(word)
        m_projector = Monomial{ProjectorAlgebra}(word)

        @test m_pauli != m_unipotent
        @test m_pauli != m_projector
        @test m_unipotent != m_projector

        # Different algebras with different integer types
        m_fermi = Monomial{FermionicAlgebra}(Int32.(word))
        @test m_pauli != m_fermi
        @test m_unipotent != m_fermi
    end

    @testset "Power Operator" begin
        # Basic power operations
        m = Monomial{NonCommutativeAlgebra}([1, 2])

        # m^0 should be identity
        m0 = m^0
        @test isone(m0)
        @test degree(m0) == 0

        # m^1 should equal m
        m1 = m^1
        @test m1 == m
        @test degree(m1) == 2

        # m^2 should repeat word twice
        m2 = m^2
        @test m2.word == [1, 2, 1, 2]
        @test degree(m2) == 4

        # m^3 should repeat word three times
        m3 = m^3
        @test m3.word == [1, 2, 1, 2, 1, 2]
        @test degree(m3) == 6

        # Power of single variable
        x = Monomial{NonCommutativeAlgebra}([3])
        x5 = x^5
        @test x5.word == [3, 3, 3, 3, 3]
        @test degree(x5) == 5

        # Power of identity stays identity
        id = one(m)
        @test isone(id^10)

        # Different integer types
        m_uint8 = Monomial{PauliAlgebra}(UInt8[1, 3])
        m_sq = m_uint8^2
        @test m_sq.word == UInt8[1, 3, 1, 3]
        @test typeof(m_sq) == Monomial{PauliAlgebra,UInt8}

        m_int32 = Monomial{FermionicAlgebra}(Int32[1, -2])
        m_cubed = m_int32^3
        @test m_cubed.word == Int32[1, -2, 1, -2, 1, -2]

        # Power distribution: m^(a+b) = m^a * m^b
        m_test = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        left = m_test^3
        right = m_test^2 * m_test
        @test left == right

        right2 = m_test * m_test^2
        @test left == right2
    end

    @testset "Monomial addition" begin
        m1 = Monomial{NonCommutativeAlgebra}(UInt8[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2])

        # Test monomial + monomial
        p = m1 + m2
        @test p isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(p)) == 2
        @test coefficients(p) == [1.0, 1.0]

        # Test subtraction
        p2 = m1 - m2
        @test p2 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(p2)) == 2
        coeffs = coefficients(p2)
        @test coeffs[1] == 1.0
        @test coeffs[2] == -1.0

        # Test negation
        t = -m1
        @test t isa Term{Monomial{NonCommutativeAlgebra, UInt8}, Float64}
        @test t.coefficient == -1.0
        @test t.monomial === m1

        # Test scalar addition
        p3 = m1 + 2.0
        @test p3 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(p3)) == 2

        # Test scalar subtraction
        p4 = m1 - 3.0
        @test p4 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}

        p5 = 4.0 - m2
        @test p5 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}

        # Test with powers
        m_sq = m1^2
        p6 = m1 + m_sq
        @test p6 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(p6)) == 2

        # Test that monomial + monomial handles different monomials
        mono_diff = Monomial{NonCommutativeAlgebra}(UInt8[3, 4])
        p7 = m2 + mono_diff
        @test p7 isa Polynomial{NonCommutativeAlgebra, UInt8, Float64}
        @test length(terms(p7)) == 2
    end
end
