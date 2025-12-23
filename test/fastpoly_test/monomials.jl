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

    @testset "Adjoint Operation" begin
        # Adjoint reverses word for self-adjoint algebras (unsigned types)
        mono1 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        mono1_adj = adjoint(mono1)
        @test mono1_adj.word == [3, 2, 1]

        # Empty monomial adjoint should be empty
        mono_empty = Monomial{NonCommutativeAlgebra}(UInt8[])
        @test isone(adjoint(mono_empty))

        # Single element adjoint
        mono_single = Monomial{NonCommutativeAlgebra}(UInt8[5])
        @test adjoint(mono_single).word == [5]

        # Adjoint is involution: adjoint(adjoint(m)) == m
        mono2 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3, 4])
        @test adjoint(adjoint(mono2)) == mono2

        # Julia syntax shorthand
        @test mono1' == adjoint(mono1)

        mono_grouped = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 20, 21, 3])
        # simplification are done separately from adjoint in this version
        @test adjoint(mono_grouped).word == UInt8[3, 21, 20, 2, 1]
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

        # Julia syntax shorthand
        @test m_ferm' == m_adj
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

    @testset "Cross-Algebra Safety" begin
        # Test that hash/equality contract is maintained across algebras
        m_pauli = Monomial{PauliAlgebra}([1, 2, 3])
        m_nc = Monomial{NonCommutativeAlgebra}([1, 2, 3])

        # Same word, different algebra = different hash (ensures Dict safety)
        @test hash(m_pauli) != hash(m_nc)

        # Verify isless throws descriptive error for cross-algebra comparison
        @test_throws ArgumentError isless(m_pauli, m_nc)
        @test_throws ArgumentError m_pauli < m_nc

        # Verify error message is informative
        try
            isless(m_pauli, m_nc)
            @test false  # Should not reach here
        catch e
            @test e isa ArgumentError
            @test occursin("PauliAlgebra", string(e))
            @test occursin("NonCommutativeAlgebra", string(e))
        end
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
        @test p isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p)) == 2
        @test coefficients(p) == [1.0, 1.0]

        # Test subtraction
        p2 = m1 - m2
        @test p2 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p2)) == 2
        coeffs = coefficients(p2)
        @test coeffs[1] == 1.0
        @test coeffs[2] == -1.0

        # Test negation
        t = -m1
        @test t isa Term{Monomial{NonCommutativeAlgebra,UInt8},Float64}
        @test t.coefficient == -1.0
        @test t.monomial === m1

        # Test scalar addition
        p3 = m1 + 2.0
        @test p3 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p3)) == 2

        # Test scalar subtraction
        p4 = m1 - 3.0
        @test p4 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}

        p5 = 4.0 - m2
        @test p5 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}

        # Test with powers
        m_sq = m1^2
        p6 = m1 + m_sq
        @test p6 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p6)) == 2

        # Test that monomial + monomial handles different monomials
        mono_diff = Monomial{NonCommutativeAlgebra}(UInt8[3, 4])
        p7 = m2 + mono_diff
        @test p7 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p7)) == 2
    end

    @testset "Iteration Protocol" begin
        # Monomial iteration yields single (coefficient, monomial) pair
        m = Monomial{PauliAlgebra}([1, 2])
        pairs = collect(m)
        @test length(pairs) == 1
        @test pairs[1] == (1.0 + 0.0im, m)  # ComplexF64 for Pauli

        # Float64 coefficient for non-Pauli algebras
        m_nc = Monomial{NonCommutativeAlgebra}([1, 2])
        pairs_nc = collect(m_nc)
        @test length(pairs_nc) == 1
        @test pairs_nc[1][1] isa Float64
        @test pairs_nc[1][1] == 1.0
        @test pairs_nc[1][2] == m_nc

        # Fermionic algebra also uses Float64
        m_fermi = Monomial{FermionicAlgebra}(Int32[1, 2])
        pairs_fermi = collect(m_fermi)
        @test pairs_fermi[1][1] isa Float64
        @test pairs_fermi[1][1] == 1.0

        # Destructuring pattern
        (coef, mono), = m_nc
        @test coef == 1.0
        @test mono == m_nc

        # eltype
        @test eltype(typeof(m)) == Tuple{ComplexF64,Monomial{PauliAlgebra,Int64}}
        @test eltype(typeof(m_nc)) == Tuple{Float64,Monomial{NonCommutativeAlgebra,Int64}}

        # Manual iteration
        iter = iterate(m)
        @test iter !== nothing
        @test iter[1] == (1.0 + 0.0im, m)
        @test iterate(m, iter[2]) === nothing

        # Identity monomial
        m_id = one(Monomial{PauliAlgebra,Int64})
        pairs_id = collect(m_id)
        @test length(pairs_id) == 1
        @test pairs_id[1][1] == 1.0 + 0.0im
        @test isone(pairs_id[1][2])
    end

    @testset "Display with Registry and Exponents" begin
        # Create a simple test registry
        struct TestRegistry
            idx_to_variables::Dict{UInt8,Symbol}
        end

        reg = TestRegistry(Dict(
            UInt8(1) => :x,
            UInt8(2) => :y,
            UInt8(3) => :z
        ))

        function test_display(m::Monomial, expected_contains::Vector{String};
            expected_not_contains::Vector{String}=String[])
            buf = IOBuffer()
            show(IOContext(buf, :registry => reg), m)
            str = String(take!(buf))

            for expected in expected_contains
                @test occursin(expected, str)
            end
            for not_expected in expected_not_contains
                @test !occursin(not_expected, str)
            end
            return str
        end

        # Test x^3 displays with exponent
        m1 = Monomial{NonCommutativeAlgebra}(UInt8[1, 1, 1])
        str1 = test_display(m1, ["x³"], expected_not_contains=["xxx", "x^3", "^[", ","])

        # Test x*y^2
        m2 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 2])
        str2 = test_display(m2, ["x", "y²"], expected_not_contains=["[2, 2]", "^2", "^3"])

        # Test x^2*y^3
        m3 = Monomial{NonCommutativeAlgebra}(UInt8[1, 1, 2, 2, 2])
        str3 = test_display(m3, ["x²", "y³"], expected_not_contains=["[1, 1]", "[2, 2, 2]"])

        # Test single variable no exponent
        m4 = Monomial{NonCommutativeAlgebra}(UInt8[1])
        str4 = test_display(m4, ["x"], expected_not_contains=["^", "²", "³"])

        # Test three different variables, no repetition
        m5 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        str5 = test_display(m5, ["x", "y", "z"], expected_not_contains=["^", "²", "³"])

        # Test large exponent (10) - should use ^10 syntax
        m6 = Monomial{NonCommutativeAlgebra}(UInt8[fill(1, 10)..., 2])
        str6 = test_display(m6, ["x^10", "y"], expected_not_contains=["xxxxxxxxxx"])

        # Test exponent 4-9 for superscript coverage
        m7 = Monomial{NonCommutativeAlgebra}(UInt8[1, 1, 1, 1])  # 4 x's
        str7 = test_display(m7, ["x⁴"], expected_not_contains=["^4"])

        m8 = Monomial{NonCommutativeAlgebra}(UInt8[fill(1, 9)..., 2])  # 9 x's
        str8 = test_display(m8, ["x⁹", "y"], expected_not_contains=["^9"])

        # Test identity displays as 𝟙 symbol
        m_identity = Monomial{NonCommutativeAlgebra}(UInt8[])
        buf = IOBuffer()
        show(IOContext(buf, :registry => reg), m_identity)
        str_identity = String(take!(buf))
        @test contains(str_identity, "𝟙")
    end
end
