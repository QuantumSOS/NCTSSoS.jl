# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: variable_indices, expval, encode_index

@testset "Monomials" begin
    # =============================================================================
    # Setup: 4 variables (2 operators × 2 sites) for proper commutativity testing
    # =============================================================================
    nc_s1_op1 = encode_index(UInt, 1, 1)  # operator 1 @ site 1
    nc_s1_op2 = encode_index(UInt, 2, 1)  # operator 2 @ site 1
    nc_s2_op1 = encode_index(UInt, 1, 2)  # operator 1 @ site 2
    nc_s2_op2 = encode_index(UInt, 2, 2)  # operator 2 @ site 2

    @testset "NormalMonomial Properties" begin
        # Data-driven tests: single construction tests multiple properties
        test_cases = [
            # (word, degree, isone, var_count, description)
            (UInt[], 0, true, 0, "identity"),
            ([nc_s1_op1], 1, false, 1, "single operator"),
            ([nc_s1_op1, nc_s1_op2], 2, false, 2, "two ops same site"),
            ([nc_s1_op1, nc_s2_op1], 2, false, 2, "two ops different sites"),
            ([nc_s1_op1, nc_s1_op2, nc_s2_op1, nc_s2_op2], 4, false, 4, "degree-4 all vars"),
        ]

        for (word, exp_degree, exp_isone, exp_var_count, desc) in test_cases
            @testset "$desc" begin
                m = NormalMonomial{NonCommutativeAlgebra,UInt}(word)
                @test degree(m) == exp_degree
                @test isone(m) == exp_isone
                @test length(variable_indices(m)) == exp_var_count
                @test isone(one(m))
                @test typeof(one(m)) == typeof(m)
            end
        end
    end

    @testset "Simplify Site Reordering" begin
        # Degree-4 monomial with interleaved sites and reversed intra-site order
        # Tests: (1) cross-site commutativity, (2) intra-site order preservation
        raw_word = [nc_s2_op2, nc_s1_op1, nc_s2_op1, nc_s1_op2]
        result = simplify(NonCommutativeAlgebra, raw_word)

        # Sites grouped (s1 before s2), intra-site order preserved
        expected = [nc_s1_op1, nc_s1_op2, nc_s2_op2, nc_s2_op1]
        @test result == expected
    end

    @testset "Hash Contract" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2])
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s2_op1])

        @test hash(m1) == hash(m2)  # equal monomials have equal hashes
        @test hash(m1) != hash(m3)  # different monomials have different hashes
    end

    @testset "Equality" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2])
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op2, nc_s1_op1])

        @test m1 == m2
        @test m1 != m3  # order preserved within same site
    end

    @testset "Comparison and Sorting" begin
        m_deg1 = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1])
        m_deg2a = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2])
        m_deg2b = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s2_op1])

        # Degree-first ordering
        @test isless(m_deg1, m_deg2a)

        # Sorting works
        monos = [m_deg2b, m_deg1, m_deg2a]
        sorted = sort(monos)
        @test sorted[1] == m_deg1
        @test degree(sorted[2]) == 2
        @test degree(sorted[3]) == 2
    end

    @testset "Zero Rejection" begin
        @test_throws ArgumentError NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[0, 1, 2])
        @test_throws ArgumentError NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[0])
    end

    @testset "Normal Form Validation" begin
        # Uses 4 variables (2 per site) to test both site-sorting and simplification

        # --- NonCommutativeAlgebra: requires site-sorted ---
        nc16_s1_op1 = encode_index(UInt16, 1, 1)
        nc16_s1_op2 = encode_index(UInt16, 2, 1)
        nc16_s2_op1 = encode_index(UInt16, 1, 2)
        nc16_s2_op2 = encode_index(UInt16, 2, 2)

        @test_throws ArgumentError NormalMonomial{NonCommutativeAlgebra,UInt16}(
            UInt16[nc16_s2_op1, nc16_s1_op1, nc16_s2_op2, nc16_s1_op2]
        )

        # --- PauliAlgebra: requires site-sorted, ≤1 operator per site ---
        pauli_s1_x, pauli_s1_y = UInt16(1), UInt16(2)
        pauli_s2_x = UInt16(4)

        @test_throws ArgumentError NormalMonomial{PauliAlgebra,UInt16}(
            UInt16[pauli_s1_x, pauli_s1_y]  # multiple ops same site
        )
        @test_throws ArgumentError NormalMonomial{PauliAlgebra,UInt16}(
            UInt16[pauli_s2_x, pauli_s1_x]  # non-site-sorted
        )

        # --- FermionicAlgebra: requires normal-ordered, no duplicates ---
        fermi_c1, fermi_c2 = Int32(-1), Int32(-2)
        fermi_a1, fermi_a2 = Int32(1), Int32(2)

        @test_throws ArgumentError NormalMonomial{FermionicAlgebra,Int32}(
            Int32[fermi_a1, fermi_c1]  # annihilator before creator
        )
        @test_throws ArgumentError NormalMonomial{FermionicAlgebra,Int32}(
            Int32[fermi_c1, fermi_c1]  # duplicate (nilpotent)
        )

        # --- BosonicAlgebra: requires normal-ordered (duplicates allowed) ---
        bos_c1, bos_a1 = Int32(-1), Int32(1)

        @test_throws ArgumentError NormalMonomial{BosonicAlgebra,Int32}(
            Int32[bos_a1, bos_c1]  # annihilator before creator
        )

        # --- UnipotentAlgebra: site-sorted, no consecutive duplicates (U² = I) ---
        uni_s1_op1 = encode_index(UInt16, 1, 1)
        uni_s2_op1 = encode_index(UInt16, 1, 2)

        @test_throws ArgumentError NormalMonomial{UnipotentAlgebra,UInt16}(
            UInt16[uni_s1_op1, uni_s1_op1]  # consecutive identical
        )
        @test_throws ArgumentError NormalMonomial{UnipotentAlgebra,UInt16}(
            UInt16[uni_s2_op1, uni_s1_op1]  # non-site-sorted
        )

        # --- ProjectorAlgebra: site-sorted, no consecutive duplicates (P² = P) ---
        proj_s1_op1 = encode_index(UInt16, 1, 1)
        proj_s2_op1 = encode_index(UInt16, 1, 2)

        @test_throws ArgumentError NormalMonomial{ProjectorAlgebra,UInt16}(
            UInt16[proj_s1_op1, proj_s1_op1]  # consecutive identical
        )
        @test_throws ArgumentError NormalMonomial{ProjectorAlgebra,UInt16}(
            UInt16[proj_s2_op1, proj_s1_op1]  # non-site-sorted
        )
    end

    @testset "Cross-Algebra Safety" begin
        # Pauli requires Unsigned types
        m_pauli = NormalMonomial{PauliAlgebra,UInt}(UInt[1, 4, 7])
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[1, 4, 7])

        # Different algebras are never equal
        @test m_pauli != m_nc

        # Different algebra = different hash
        @test hash(m_pauli) != hash(m_nc)

        # Cross-algebra comparison throws
        @test_throws ArgumentError isless(m_pauli, m_nc)

        # Error message is informative
        try
            isless(m_pauli, m_nc)
        catch e
            @test occursin("PauliAlgebra", string(e))
            @test occursin("NonCommutativeAlgebra", string(e))
        end

    end

    @testset "Algebra Type Preservation" begin
        algebras_and_words = [
            (PauliAlgebra, UInt16, UInt16[1, 4]),  # Pauli requires Unsigned
            (FermionicAlgebra, Int32, Int32[-1, 1]),
            (BosonicAlgebra, Int32, Int32[-1, 1]),
            (UnipotentAlgebra, UInt16, UInt16[1, 2]),
            (ProjectorAlgebra, UInt16, UInt16[1, 2]),
        ]

        for (A, T, word) in algebras_and_words
            m = NormalMonomial{A,T}(word)
            @test m isa NormalMonomial{A,T}
            @test one(m) isa NormalMonomial{A,T}
            @test isone(one(m))
        end
    end

    @testset "Multiplication" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Multiplication returns Polynomial
        result = x[1] * x[2]
        @test result isa Polynomial
        @test degree(result) == 2

        # Identity multiplication: result is Polynomial containing the same monomial
        mono_id = one(x[1])
        prod = mono_id * x[1]
        @test prod isa Polynomial
        @test length(terms(prod)) == 1
        @test monomials(prod)[1] == x[1]
    end

    @testset "Power Operator" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2])

        # m^0 returns identity Polynomial
        p0 = m^0
        @test p0 isa Polynomial
        @test length(terms(p0)) == 1
        @test isone(monomials(p0)[1])

        # m^1 returns Polynomial with same monomial
        p1 = m^1
        @test monomials(p1)[1] == m

        # m^2 concatenates word
        p2 = m^2
        expected_word = [nc_s1_op1, nc_s1_op2, nc_s1_op1, nc_s1_op2]
        @test monomials(p2)[1].word == expected_word

        @test degree(m^3) == 6

        # Power of identity stays identity
        p_id = one(m)^10
        @test length(terms(p_id)) == 1
        @test isone(monomials(p_id)[1])

        # Power distribution: m^(a+b) = m^a * m^b
        @test m^3 == m^2 * m
        @test m^3 == m * m^2

        # Pauli power: (σx⊗σx)² = I
        m_pauli = NormalMonomial{PauliAlgebra,UInt8}(UInt8[1, 4])
        p_pauli_sq = m_pauli^2
        @test isone(monomials(p_pauli_sq)[1])

        # Fermionic number operator: n² = n
        m_num = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 1])
        p_num_sq = m_num^2
        @test monomials(p_num_sq)[1] == m_num
    end

    @testset "Addition and Subtraction" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[1])
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[1, 2])

        p_add = m1 + m2
        @test p_add isa Polynomial{NonCommutativeAlgebra,UInt8,Float64}
        @test length(terms(p_add)) == 2

        p_sub = m1 - m2
        @test coefficients(p_sub) == [1.0, -1.0]

        # Negation returns Polynomial with -1 coefficient
        p_neg = -m1
        @test p_neg isa Polynomial
        @test length(terms(p_neg)) == 1
        @test coefficients(p_neg)[1] == -1.0
        @test monomials(p_neg)[1] == m1

        # Scalar arithmetic
        @test (m1 + 2.0) isa Polynomial
        @test (m1 - 3.0) isa Polynomial
        @test (4.0 - m2) isa Polynomial
    end

    @testset "Simplify Word Protocol" begin
        # simplify(AlgebraType, word) returns simplified word/result

        # NonCommutative: returns sorted word
        raw_nc = [nc_s2_op1, nc_s1_op1]  # site2 before site1
        result_nc = simplify(NonCommutativeAlgebra, raw_nc)
        @test result_nc == [nc_s1_op1, nc_s2_op1]  # sorted by site

        # Pauli: returns (word, phase) where phase is accumulated
        pauli_word = UInt16[1, 4]  # σx on sites 1,2
        result_pauli = simplify(PauliAlgebra, pauli_word)
        @test result_pauli[1] == pauli_word  # already canonical
        @test result_pauli[2] == UInt8(0)    # no phase

        # Fermionic: returns Vector{Tuple{Int, word}} for PBW expansion
        fermi_word = Int32[-2, 1]  # creation then annihilation (normal-ordered)
        result_fermi = simplify(FermionicAlgebra, fermi_word)
        @test result_fermi[1][1] == 1  # coefficient
        @test result_fermi[1][2] == fermi_word  # already normal-ordered
    end

    @testset "Display" begin
        struct TestRegistry
            idx_to_variables::Dict{UInt8,Symbol}
        end

        reg = TestRegistry(Dict(UInt8(1) => :x, UInt8(2) => :y, UInt8(3) => :z))

        display_cases = [
            (UInt8[1, 1, 1], ["x³"], ["xxx"]),
            (UInt8[1, 2, 2], ["x", "y²"], String[]),
            (UInt8[1], ["x"], ["²", "³"]),
            (UInt8[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2], ["x^10", "y"], String[]),
        ]

        for (word, must_contain, must_not_contain) in display_cases
            m = NormalMonomial{NonCommutativeAlgebra,UInt8}(word)
            buf = IOBuffer()
            show(IOContext(buf, :registry => reg), m)
            str = String(take!(buf))
            for s in must_contain
                @test occursin(s, str)
            end
            for s in must_not_contain
                @test !occursin(s, str)
            end
        end

        # Identity displays as 𝟙
        m_id = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[])
        buf = IOBuffer()
        show(IOContext(buf, :registry => reg), m_id)
        @test contains(String(take!(buf)), "𝟙")
    end

    @testset "variable_indices" begin
        # Pauli: unique indices (requires Unsigned)
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[1, 4, 7])
        @test variable_indices(m_pauli) == Set(UInt16[1, 4, 7])

        # NonCommutative: repeated indices deduplicated (use site-sorted word)
        # Same-site indices can repeat: [s1_op1, s1_op2, s1_op1, s1_op2]
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt}([nc_s1_op1, nc_s1_op2, nc_s1_op1, nc_s1_op2])
        @test variable_indices(m_nc) == Set([nc_s1_op1, nc_s1_op2])

        # Identity: empty set
        @test isempty(variable_indices(one(NormalMonomial{PauliAlgebra,UInt16})))

        # Fermionic: abs() normalization (creation=-1, annihilation=1 → same mode)
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[-3, -2, -1, 1, 3])
        @test variable_indices(m_fermi) == Set(Int32[1, 2, 3])

        # Type preservation
        m_uint8 = NormalMonomial{ProjectorAlgebra,UInt8}(UInt8[1, 2, 3])
        @test eltype(variable_indices(m_uint8)) == UInt8
    end
end
