# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: variable_indices, expval

@testset "Monomials" begin
    @testset "Creation" begin
        # Create monomials directly with word representation
        mono1 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 1, 3])
        @test mono1.word == [1, 2, 1, 3]
        @test degree(mono1) == 4

        # Empty monomial (identity)
        mono_identity = NormalMonomial{NonCommutativeAlgebra}(Int[])
        @test isone(mono_identity)
        @test degree(mono_identity) == 0

        # Single element monomial
        mono_single = NormalMonomial{NonCommutativeAlgebra}([5])
        @test mono_single.word == [5]
        @test degree(mono_single) == 1

        # Monomial with zeros should filter them out
        mono_with_zero = NormalMonomial{NonCommutativeAlgebra}([1, 0, 2, 0, 3])
        @test mono_with_zero.word == [1, 2, 3]
        @test degree(mono_with_zero) == 3

        # one(NormalMonomial) should return identity
        mono_one = one(NormalMonomial{NonCommutativeAlgebra,Int64})
        @test isone(mono_one)
    end

    @testset "Degree" begin
        mono1 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = NormalMonomial{NonCommutativeAlgebra}(Int[])
        mono3 = NormalMonomial{NonCommutativeAlgebra}([1, 1, 1, 1, 1])

        @test degree(mono1) == 3
        @test degree(mono2) == 0
        @test degree(mono3) == 5
    end

    @testset "Hash" begin
        mono1 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = NormalMonomial{NonCommutativeAlgebra}([1, 3, 3])
        mono3 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3])

        @test hash(mono1) != hash(mono2)
        @test hash(mono1) == hash(mono3)

        mono_empty = NormalMonomial{NonCommutativeAlgebra}(Int[])
        mono_zero_filtered = NormalMonomial{NonCommutativeAlgebra}([0])
        @test hash(mono_empty) == hash(mono_zero_filtered)
    end

    @testset "Adjoint Operation" begin
        # Adjoint reverses word for self-adjoint encodings (unsigned types).
        # Use indices on the same site so site-based canonicalization does not reorder.
        using NCTSSoS: encode_index
        idx1_s1 = encode_index(UInt8, 1, 1)
        idx2_s1 = encode_index(UInt8, 2, 1)
        idx3_s1 = encode_index(UInt8, 3, 1)
        mono1 = NormalMonomial{NonCommutativeAlgebra}(UInt8[idx1_s1, idx2_s1, idx3_s1])
        mono1_adj = adjoint(mono1)
        @test mono1_adj.word == UInt8[idx3_s1, idx2_s1, idx1_s1]

        # Empty monomial adjoint should be empty
        mono_empty = NormalMonomial{NonCommutativeAlgebra}(UInt8[])
        @test isone(adjoint(mono_empty))

        # Single element adjoint
        mono_single = NormalMonomial{NonCommutativeAlgebra}(UInt8[idx2_s1])
        @test adjoint(mono_single).word == UInt8[idx2_s1]

        # Adjoint is involution: adjoint(adjoint(m)) == m
        # Use indices on same site to avoid sorting interaction
        idx4_s1 = encode_index(UInt8, 4, 1)
        mono2 = NormalMonomial{NonCommutativeAlgebra}(UInt8[idx1_s1, idx2_s1, idx3_s1, idx4_s1])
        @test adjoint(adjoint(mono2)) == mono2

        # Julia syntax shorthand
        @test mono1' == adjoint(mono1)

        # Multi-site monomial using encode_index for predictable behavior
        # Note: Constructor auto-sorts by site (stable sort)
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx2_s1 = encode_index(UInt16, 2, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        # Input [idx2_s1, idx1_s2, idx1_s1] sorts by site to [idx2_s1, idx1_s1, idx1_s2]
        mono_multi = NormalMonomial{NonCommutativeAlgebra}(UInt16[idx2_s1, idx1_s2, idx1_s1])
        # Adjoint reverses and then re-canonicalizes by site
        @test adjoint(mono_multi).word == UInt16[idx1_s1, idx2_s1, idx1_s2]
    end

    @testset "Multiplication" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Same variable multiplication produces longer word (concatenation)
        result = x[1] * x[1]
        @test result isa Monomial
        @test degree(result) == 2

        # Different variables
        result2 = x[1] * x[2]
        @test degree(result2) == 2

        # Identity multiplication
        mono_id = one(x[1])
        @test mono_id * x[1] == x[1]
        @test x[1] * mono_id == x[1]
    end

    @testset "Comparison" begin
        mono1 = NormalMonomial{NonCommutativeAlgebra}([1])
        mono2 = NormalMonomial{NonCommutativeAlgebra}([1, 2])
        mono3 = NormalMonomial{NonCommutativeAlgebra}([2])

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
        mono1 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3])
        mono3 = NormalMonomial{NonCommutativeAlgebra}([1, 3, 2])

        @test mono1 == mono2
        @test mono1 != mono3

        # Different algebra types are never equal
        # Use indices on different Pauli sites: 1,4,7 = σx on sites 1,2,3
        mono_pauli = NormalMonomial{PauliAlgebra}([1, 4, 7])
        mono_nc = NormalMonomial{NonCommutativeAlgebra}([1, 4, 7])
        @test mono_pauli != mono_nc
    end

    @testset "Algebra Type Preservation" begin
        # Monomials preserve their algebra type
        # Use canonical forms: Pauli indices on diff sites, Fermi normal-ordered
        mono_pauli = NormalMonomial{PauliAlgebra}([1, 4])  # sites 1,2
        mono_fermi = NormalMonomial{FermionicAlgebra}(Int32[-1, 1])  # creation before annihilation
        mono_unipotent = NormalMonomial{UnipotentAlgebra}(UInt16[1, 2])

        @test mono_pauli isa NormalMonomial{PauliAlgebra}
        @test mono_fermi isa NormalMonomial{FermionicAlgebra}
        @test mono_unipotent isa NormalMonomial{UnipotentAlgebra}
    end

    @testset "Adjoint for Fermionic/Signed Types" begin
        # Test FermionicAlgebra adjoint (reverses AND negates)
        # Use normal-ordered: creations (negative) first, then annihilations (positive)
        # Sorted by mode within each group
        m_ferm = NormalMonomial{FermionicAlgebra}(Int32[-2, 1, 3])  # a2†, a1, a3
        m_adj = adjoint(m_ferm)
        @test m_adj.word == Int32[-3, -1, 2]

        # Verify involution property for signed types: adjoint(adjoint(m)) == m
        # Note: adjoint reverses AND negates, then the result may not be normal-ordered
        # For testing involution, we use monomials where adjoint(adjoint(m)) preserves form
        # Normal-ordered: creators sorted by mode descending (4 > 2), then annihilators sorted by mode ascending (1 < 3)
        m_ferm2 = NormalMonomial{FermionicAlgebra}(Int32[-4, -2, 1, 3])  # c₄†c₂†a₁a₃ (normal-ordered)
        @test adjoint(adjoint(m_ferm2)) == m_ferm2

        # Empty fermionic monomial
        m_empty = NormalMonomial{FermionicAlgebra}(Int32[])
        @test isone(adjoint(m_empty))

        # Single element (creation)
        m_single = NormalMonomial{FermionicAlgebra}(Int32[-5])  # creation operator
        @test adjoint(m_single).word == Int32[5]  # becomes annihilation

        # Julia syntax shorthand
        @test m_ferm' == m_adj
    end

    @testset "one(m::NormalMonomial) Instance Method" begin
        # Test instance method preserves algebra type
        # Use indices on different Pauli sites: 1,4,7 = σx on sites 1,2,3
        m_pauli = NormalMonomial{PauliAlgebra}(UInt16[1, 4, 7])
        @test one(m_pauli) == NormalMonomial{PauliAlgebra}(UInt16[])
        @test typeof(one(m_pauli)) == NormalMonomial{PauliAlgebra,UInt16}

        # Test with FermionicAlgebra (normal-ordered: creation first)
        m_ferm = NormalMonomial{FermionicAlgebra}(Int32[-2, 1])  # a2†, a1
        @test one(m_ferm) == NormalMonomial{FermionicAlgebra}(Int32[])
        @test typeof(one(m_ferm)) == NormalMonomial{FermionicAlgebra,Int32}

        # Test that one(m) returns identity
        @test isone(one(m_pauli))
        @test isone(one(m_ferm))

        # Verify one(type) and one(instance) return same result
        @test one(NormalMonomial{PauliAlgebra,UInt16}) == one(m_pauli)
    end

    @testset "NormalMonomial Default Constructor" begin
        # Test default constructor uses NonCommutativeAlgebra
        m = NormalMonomial([1, 2, 3])
        @test m isa NormalMonomial{NonCommutativeAlgebra}
        @test m.word == [1, 2, 3]

        # Test with different integer types
        m_int32 = NormalMonomial(Int32[1, 2])
        @test m_int32 isa NormalMonomial{NonCommutativeAlgebra,Int32}

        m_uint16 = NormalMonomial(UInt16[1, 2])
        @test m_uint16 isa NormalMonomial{NonCommutativeAlgebra,UInt16}

        # Test zero filtering in default constructor
        m_zeros = NormalMonomial([1, 0, 2, 0])
        @test m_zeros.word == [1, 2]

        # Empty default constructor
        m_empty = NormalMonomial(Int[])
        @test isone(m_empty)
    end

    @testset "Zero Filtering Edge Cases" begin
        # Leading zeros
        m1 = NormalMonomial{NonCommutativeAlgebra}([0, 0, 1, 2])
        @test m1.word == [1, 2]

        # Trailing zeros
        m2 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 0, 0])
        @test m2.word == [1, 2]

        # All zeros -> identity
        m3 = NormalMonomial{NonCommutativeAlgebra}([0, 0, 0])
        @test isone(m3)

        # Mixed zeros throughout
        m4 = NormalMonomial{NonCommutativeAlgebra}([0, 1, 0, 2, 0, 3, 0])
        @test m4.word == [1, 2, 3]

        # Single zero
        m5 = NormalMonomial{NonCommutativeAlgebra}([0])
        @test isone(m5)
    end

    @testset "Cross-Algebra Type Equality" begin
        # Same word, different algebras should never be equal
        # Use indices on different Pauli sites: 1,4,7 = σx on sites 1,2,3
        word_pauli = [1, 4, 7]
        word_unsigned = UInt16[1, 4, 7]
        m_pauli = NormalMonomial{PauliAlgebra}(word_pauli)
        m_unipotent = NormalMonomial{UnipotentAlgebra}(word_unsigned)
        m_projector = NormalMonomial{ProjectorAlgebra}(word_unsigned)

        @test m_pauli != m_unipotent
        @test m_pauli != m_projector
        @test m_unipotent != m_projector

        # Different algebras with different integer types
        # Use normal-ordered for fermionic: creations first, then annihilations
        m_fermi = NormalMonomial{FermionicAlgebra}(Int32[-4, -1, 7])
        @test m_pauli != m_fermi
        @test m_unipotent != m_fermi
    end

    @testset "Cross-Algebra Safety" begin
        # Test that hash/equality contract is maintained across algebras
        # Use indices on different Pauli sites: 1,4,7 = σx on sites 1,2,3
        m_pauli = NormalMonomial{PauliAlgebra}([1, 4, 7])
        m_nc = NormalMonomial{NonCommutativeAlgebra}([1, 4, 7])

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
        m = NormalMonomial{NonCommutativeAlgebra}([1, 2])

        # m^0 should be identity
        m0 = m^0
        @test isone(m0)
        @test degree(m0) == 0

        # m^1 should equal m
        m1 = m^1
        @test m1.words == m
        @test degree(m1) == 2

        # m^2 should repeat word twice
        m2 = m^2
        @test m2.words.word == [1, 2, 1, 2]
        @test degree(m2) == 4

        # m^3 should repeat word three times
        m3 = m^3
        @test m3.words.word == [1, 2, 1, 2, 1, 2]
        @test degree(m3) == 6

        # Power of single variable
        x = NormalMonomial{NonCommutativeAlgebra}([3])
        x5 = x^5
        @test x5.words.word == [3, 3, 3, 3, 3]
        @test degree(x5) == 5

        # Power of identity stays identity
        id = one(m)
        @test isone(id^10)

        # Different integer types
        # Pauli: power returns a `Monomial` (scalar phase + canonical word)
        # Use indices on different sites: 1,4 = σx on sites 1,2
        m_pauli = NormalMonomial{PauliAlgebra}(UInt8[1, 4])
        pairs_sq = m_pauli^2
        @test length(pairs_sq) == 1
        phase_k, mono_sq = pairs_sq[1]
        @test phase_k == UInt8(0)
        @test isone(mono_sq)

        # Fermionic: power returns a `Monomial` (PBW expansion; integer coefficients)
        m_num = NormalMonomial{FermionicAlgebra}(Int32[-1, 1])  # number operator n₁
        pairs_num_sq = m_num^2
        @test length(pairs_num_sq) == 1
        c_num, mono_num = pairs_num_sq[1]
        @test c_num == 1
        @test mono_num == m_num

        # Power distribution: m^(a+b) = m^a * m^b
        m_test = NormalMonomial{NonCommutativeAlgebra}([1, 2, 3])
        left = m_test^3
        right = m_test^2 * m_test
        @test left == right

        right2 = m_test * m_test^2
        @test left == right2
    end

    @testset "Monomial addition" begin
        m1 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1])
        m2 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 2])

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
        @test t isa Tuple{Float64,NormalMonomial{NonCommutativeAlgebra,UInt8}}
        @test first(t) == -1.0
        @test last(t) === m1

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
        mono_diff = NormalMonomial{NonCommutativeAlgebra}(UInt8[3, 4])
        p7 = m2 + mono_diff
        @test p7 isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p7)) == 2
    end

    @testset "Simplify Result Protocol" begin
        # `simplify(m)` returns `Monomial` (iterable as `(c_internal, mono)` pairs).

        # Use indices on different Pauli sites: 1,4 = σx on sites 1,2
        m = NormalMonomial{PauliAlgebra}([1, 4])
        pairs = simplify(m)
        @test length(pairs) == 1
        c_internal, mono = pairs[1]
        @test c_internal == UInt8(0)  # phase_k = 0 → 1
        @test mono == m
        @test NCTSSoS._coeff_to_number(mono, c_internal) == 1.0 + 0.0im

        # Non-Pauli monoids use unit coefficient encoding `0x01`
        m_nc = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 2])
        pairs_nc = simplify(m_nc)
        @test length(pairs_nc) == 1
        c_nc, mono_nc = pairs_nc[1]
        @test c_nc == UInt8(1)
        @test mono_nc == m_nc
        @test NCTSSoS._coeff_to_number(mono_nc, c_nc) == 1.0

        # Fermionic simplification returns Int coefficients (converted to Float64 for polynomials)
        m_fermi = NormalMonomial{FermionicAlgebra}(Int32[-2, 1])  # creation first (normal-ordered)
        pairs_fermi = simplify(m_fermi)
        @test length(pairs_fermi) == 1
        c_f, mono_f = pairs_fermi[1]
        @test c_f == 1
        @test mono_f == m_fermi
        @test NCTSSoS._coeff_to_number(mono_f, c_f) == 1.0

        # Identity monomial
        m_id = one(NormalMonomial{PauliAlgebra,Int64})
        pairs_id = simplify(m_id)
        @test length(pairs_id) == 1
        c_id, mono_id = pairs_id[1]
        @test c_id == UInt8(0)
        @test isone(mono_id)
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

        function test_display(m::NormalMonomial, expected_contains::Vector{String};
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
        m1 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 1, 1])
        str1 = test_display(m1, ["x³"], expected_not_contains=["xxx", "x^3", "^[", ","])

        # Test x*y^2
        m2 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 2, 2])
        str2 = test_display(m2, ["x", "y²"], expected_not_contains=["[2, 2]", "^2", "^3"])

        # Test x^2*y^3
        m3 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 1, 2, 2, 2])
        str3 = test_display(m3, ["x²", "y³"], expected_not_contains=["[1, 1]", "[2, 2, 2]"])

        # Test single variable no exponent
        m4 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1])
        str4 = test_display(m4, ["x"], expected_not_contains=["^", "²", "³"])

        # Test three different variables, no repetition
        m5 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        str5 = test_display(m5, ["x", "y", "z"], expected_not_contains=["^", "²", "³"])

        # Test large exponent (10) - should use ^10 syntax
        m6 = NormalMonomial{NonCommutativeAlgebra}(UInt8[fill(1, 10)..., 2])
        str6 = test_display(m6, ["x^10", "y"], expected_not_contains=["xxxxxxxxxx"])

        # Test exponent 4-9 for superscript coverage
        m7 = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 1, 1, 1])  # 4 x's
        str7 = test_display(m7, ["x⁴"], expected_not_contains=["^4"])

        m8 = NormalMonomial{NonCommutativeAlgebra}(UInt8[fill(1, 9)..., 2])  # 9 x's
        str8 = test_display(m8, ["x⁹", "y"], expected_not_contains=["^9"])

        # Test identity displays as 𝟙 symbol
        m_identity = NormalMonomial{NonCommutativeAlgebra}(UInt8[])
        buf = IOBuffer()
        show(IOContext(buf, :registry => reg), m_identity)
        str_identity = String(take!(buf))
        @test contains(str_identity, "𝟙")
    end

    @testset "variable_indices" begin
        # Basic case: unique indices (on different Pauli sites)
        m1 = NormalMonomial{PauliAlgebra}([1, 4, 7])  # sites 1,2,3
        @test variable_indices(m1) == Set([1, 4, 7])

        # Repeated indices -> deduplicated in Set
        m2 = NormalMonomial{NonCommutativeAlgebra}([1, 2, 1, 3, 2])
        @test variable_indices(m2) == Set([1, 2, 3])
        @test length(variable_indices(m2)) == 3

        # Identity monomial -> empty set
        m_id = one(NormalMonomial{PauliAlgebra,Int64})
        @test isempty(variable_indices(m_id))

        # Single index
        m_single = NormalMonomial{UnipotentAlgebra}(UInt16[42])
        @test variable_indices(m_single) == Set(UInt16[42])

        # Signed types (Fermionic): abs() normalization
        # Creation (-1) and annihilation (1) refer to same mode
        # Use normal-ordered: creations first (sorted by mode), then annihilators (sorted by mode)
        m_fermi = NormalMonomial{FermionicAlgebra}(Int32[-3, -2, -1, 1, 3])  # creators (desc), then annihilators (asc)
        var_set = variable_indices(m_fermi)
        @test var_set == Set(Int32[1, 2, 3])
        @test length(var_set) == 3

        # Bosonic also uses signed indices (normal-ordered)
        # Use normal-ordered: creators first (sorted by mode), then annihilators (sorted by mode)
        m_bos = NormalMonomial{BosonicAlgebra}(Int32[-2, -1, 1, 2])  # creators (desc), then annihilators (asc)
        @test variable_indices(m_bos) == Set(Int32[1, 2])

        # Different integer types preserve type in Set
        m_uint8 = NormalMonomial{ProjectorAlgebra}(UInt8[1, 2, 3])
        @test eltype(variable_indices(m_uint8)) == UInt8
    end

    @testset "expval (identity for NormalMonomial)" begin
        # expval is identity for regular monomials
        # (exists for API compatibility with NCStateWord)
        m = NormalMonomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        @test expval(m) === m

        # Use indices on different Pauli sites: 1,4 = σx on sites 1,2
        m_pauli = NormalMonomial{PauliAlgebra}([1, 4])
        @test expval(m_pauli) === m_pauli

        # Use normal-ordered: creation first
        m_fermi = NormalMonomial{FermionicAlgebra}(Int32[-2, 1])  # a2†, a1
        @test expval(m_fermi) === m_fermi

        # Identity monomial
        m_id = one(NormalMonomial{PauliAlgebra,Int64})
        @test expval(m_id) === m_id
    end
end
