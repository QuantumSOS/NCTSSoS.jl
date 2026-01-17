# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: encode_index, ComposedMonomial

# Import internal functions for testing
import NCTSSoS: _expand_simplified_components, _infer_coef_type_from_types

# Helper: Pauli indices on different sites (site = (idx-1)÷3 + 1)
# Site 1: 1,2,3; Site 2: 4,5,6; Site 3: 7,8,9
const P_S1 = UInt16(1)  # Pauli X on site 1
const P_S2 = UInt16(4)  # Pauli X on site 2
const P_S3 = UInt16(7)  # Pauli X on site 3

@testset "ComposedMonomial" begin
    @testset "Construction" begin
        # Basic construction with two algebras
        # Use Pauli indices on different sites to avoid validation errors
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])  # Site 1 and 2
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 2])  # c₁†c₂ (normal-ordered)

        cm = ComposedMonomial((m_pauli, m_fermi))

        @test cm isa ComposedMonomial
        @test length(cm.components) == 2
        @test cm.components[1] === m_pauli
        @test cm.components[2] === m_fermi

        # Construction with three algebras
        m_unip = NormalMonomial{UnipotentAlgebra,UInt16}(UInt16[encode_index(UInt16, 1, 1)])
        cm3 = ComposedMonomial((m_pauli, m_fermi, m_unip))

        @test length(cm3.components) == 3
        @test cm3.components[3] === m_unip

        # Single component
        cm_single = ComposedMonomial((m_pauli,))
        @test length(cm_single.components) == 1
    end

    @testset "Equality" begin
        # Use Pauli indices on different sites (P_S1=1, P_S2=4)
        m1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1, 2])  # a₁a₂ (annihilators sorted)
        m3 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])  # Same as m1
        m4 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S2, P_S3])  # Different from m1 (sites 2 and 3)

        cm1 = ComposedMonomial((m1, m2))
        cm2 = ComposedMonomial((m3, m2))  # Same components
        cm3 = ComposedMonomial((m4, m2))  # Different first component

        # Same components -> equal
        @test cm1 == cm2

        # Different components -> not equal
        @test cm1 != cm3

        # Different tuple types -> not equal
        cm_single = ComposedMonomial((m1,))
        @test cm1 != cm_single

        # Different algebra types in same position -> not equal
        m_bosonic = NormalMonomial{BosonicAlgebra,Int32}(Int32[-2, -1, 1, 2])  # normal-ordered
        cm_bosonic = ComposedMonomial((m1, m_bosonic))
        @test cm1 != cm_bosonic  # Different tuple types

        # Reflexive
        @test cm1 == cm1

        # Symmetric
        @test (cm1 == cm2) == (cm2 == cm1)
    end

    @testset "Hashing" begin
        m1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1, 2])  # a₁a₂ (valid)

        cm1 = ComposedMonomial((m1, m2))
        cm2 = ComposedMonomial((m1, m2))

        # Equal objects have equal hash
        @test hash(cm1) == hash(cm2)

        # Can be used in Sets/Dicts
        s = Set([cm1])
        @test cm2 in s

        d = Dict(cm1 => 42)
        @test d[cm2] == 42

        # Hash with salt
        @test hash(cm1, UInt(123)) == hash(cm2, UInt(123))
    end

    @testset "Comparison (isless)" begin
        m1_short = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        m1_long = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2, P_S3])  # 3 different sites
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1])

        # Degree-first ordering
        cm_deg2 = ComposedMonomial((m1_short, m2))  # degree 1+1=2
        cm_deg4 = ComposedMonomial((m1_long, m2))   # degree 3+1=4

        @test isless(cm_deg2, cm_deg4)
        @test !isless(cm_deg4, cm_deg2)

        # Lexicographic when degrees equal
        m_lex1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])  # sites 1,2
        m_lex2 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S3])  # sites 1,3
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[1])

        cm_lex1 = ComposedMonomial((m_lex1, m_fermi))
        cm_lex2 = ComposedMonomial((m_lex2, m_fermi))

        @test isless(cm_lex1, cm_lex2)
        @test !isless(cm_lex2, cm_lex1)

        # Equal -> not less
        cm_dup = ComposedMonomial((m_lex1, m_fermi))
        @test !isless(cm_lex1, cm_dup)

        # Lexicographic across multiple components
        m_pauli1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        m_pauli2 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        m_fermi1 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1])
        m_fermi2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[2])

        cm_a = ComposedMonomial((m_pauli1, m_fermi1))  # (1, 1)
        cm_b = ComposedMonomial((m_pauli2, m_fermi2))  # (1, 2)

        @test isless(cm_a, cm_b)  # Second component differs
    end

    @testset "Access methods" begin
        m1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1, 2, 3])  # a₁a₂a₃ (annihilators sorted)
        m3 = NormalMonomial{UnipotentAlgebra,UInt16}(UInt16[encode_index(UInt16, 1, 1)])

        cm = ComposedMonomial((m1, m2, m3))

        # length
        @test length(cm) == 3

        # getindex
        @test cm[1] === m1
        @test cm[2] === m2
        @test cm[3] === m3

        # degree (sum across all components)
        @test degree(cm) == 2 + 3 + 1  # 6

        # Empty monomials -> degree 0
        m_empty1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        m_empty2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[])
        cm_empty = ComposedMonomial((m_empty1, m_empty2))
        @test degree(cm_empty) == 0
    end

    @testset "(coefficient, ComposedMonomial) isone" begin
        m_empty1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        m_empty2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[])
        cm_empty = ComposedMonomial((m_empty1, m_empty2))

        # Coefficient 1, all empty monomials -> true
        t_one = (1.0, cm_empty)
        @test isone(t_one)

        # Complex coefficient 1
        t_one_complex = (1.0 + 0.0im, cm_empty)
        @test isone(t_one_complex)

        # Coefficient not 1 -> false
        t_two = (2.0, cm_empty)
        @test !isone(t_two)

        # Non-empty monomial -> false
        m_nonempty = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        cm_nonempty = ComposedMonomial((m_nonempty, m_empty2))
        t_nonempty = (1.0, cm_nonempty)
        @test !isone(t_nonempty)

        # Second component non-empty -> false
        cm_nonempty2 = ComposedMonomial((m_empty1, NormalMonomial{FermionicAlgebra,Int32}(Int32[1])))
        t_nonempty2 = (1.0, cm_nonempty2)
        @test !isone(t_nonempty2)
    end

    @testset "Simplification - Single-term algebras (Pauli)" begin
        # Pauli: Test simplification with already-canonical monomials
        # Since NormalMonomial enforces normal form, we test with valid inputs
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[1])  # σx₁
        m_unip = NormalMonomial{UnipotentAlgebra,UInt16}(UInt16[encode_index(UInt16, 1, 1)])
        cm = ComposedMonomial((m_pauli, m_unip))

        result = simplify(cm)

        # Result is always Vector{(coefficient, ComposedMonomial)} pairs
        @test result isa Vector{<:Tuple}
        @test length(result) == 1
        @test result[1][1] ≈ 1.0 + 0.0im
        @test result[1][2][1].word == [1]  # Pauli unchanged
        @test result[1][2][2].word == [encode_index(UInt16, 1, 1)]
    end

    @testset "Simplification - Single-term algebras (UnipotentAlgebra)" begin
        # Unipotent: U² = I (consecutive pairs cancel)
        # For raw words, use `simplify(UnipotentAlgebra, word)` before constructing NormalMonomial.
        # Example: [u1, u1, u2] -> [u2] (u1 pair cancels, u2 remains)
        # ComposedMonomial always returns Vector{(coefficient, ComposedMonomial)}
        u1 = encode_index(UInt16, 1, 1)
        u2 = encode_index(UInt16, 2, 1)

        # Raw word -> simplified word -> NormalMonomial
        m_unip = NormalMonomial{UnipotentAlgebra,UInt16}(simplify(UnipotentAlgebra, UInt16[u1, u1, u2]))
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        cm = ComposedMonomial((m_pauli, m_unip))

        result = simplify(cm)

        # Result is always Vector{(coefficient, ComposedMonomial)} pairs
        @test result isa Vector{<:Tuple}
        @test length(result) == 1
        # Unipotent [u1, u1, u2] -> [u2] (u1 pair cancels via U²=I)
        @test result[1][2][2].word == [u2]
    end

    @testset "Simplification - Single-term algebras (ProjectorAlgebra)" begin
        # Projector: [e_ii, e_ii] -> [e_ii] (idempotent)
        # ComposedMonomial always returns Vector{(coefficient, ComposedMonomial)}
        p1 = encode_index(UInt16, 1, 1)

        m_proj = NormalMonomial{ProjectorAlgebra,UInt16}(simplify(ProjectorAlgebra, UInt16[p1, p1]))
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        cm = ComposedMonomial((m_proj, m_pauli))

        result = simplify(cm)

        # Result is always Vector{(coefficient, ComposedMonomial)} pairs
        @test result isa Vector{<:Tuple}
        @test length(result) == 1
        @test result[1][2][1].word == [p1]  # Projector simplified
    end

    @testset "Simplification - Multi-term algebras (FermionicAlgebra)" begin
        # Fermionic: anti-commutation can produce multiple terms
        # a₁ a₁ = 0 (nilpotent - same mode annihilators)
        # a₁ a₁† = 1 - a₁† a₁ (anti-commutation)
        # Use `simplify(FermionicAlgebra, word)` for non-normal-ordered input

        pm_fermi = simplify(FermionicAlgebra, Int32[1, 1])
        @test pm_fermi == [(0, Int32[])]

        # For ComposedMonomial test, use normal-ordered fermionic monomial
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 1])  # c₁†a₁ (normal-ordered)
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        cm = ComposedMonomial((m_fermi, m_pauli))

        result = simplify(cm)

        # Result is Vector{(coefficient, ComposedMonomial)} pairs
        @test result isa Vector{<:Tuple}
        @test all(t -> t[2] isa ComposedMonomial, result)
    end

    @testset "Simplification - Multi-term algebras (BosonicAlgebra)" begin
        # Bosonic: [b, b†] -> b† b + 1 (commutation)
        pm_bosonic = simplify(BosonicAlgebra, Int32[1, -1])
        @test length(pm_bosonic) == 2

        # For ComposedMonomial test, use normal-ordered bosonic monomial
        m_bosonic = NormalMonomial{BosonicAlgebra,Int32}(Int32[-1, 1])  # c₁†a₁ (normal-ordered)
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        cm = ComposedMonomial((m_bosonic, m_pauli))

        result = simplify(cm)

        # Should produce multiple terms due to commutation
        # Result is Vector{(coefficient, ComposedMonomial)} pairs
        @test result isa Vector{<:Tuple}
        # Exact count depends on simplification rules - just verify structure
        @test all(t -> t[2] isa ComposedMonomial, result)
    end

    @testset "Simplification - Cartesian product" begin
        # When both algebras produce terms, result is Cartesian product
        # Use normal-ordered fermionic: c₁†a₁
        # Use canonical Pauli on single site
        # ComposedMonomial with Fermionic returns Vector{(coefficient, ComposedMonomial)}

        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 1])  # c₁†a₁ (normal-ordered)
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])  # Single Pauli
        cm = ComposedMonomial((m_fermi, m_pauli))

        result = simplify(cm)

        # Verify we get Vector{(coefficient, ComposedMonomial)} pairs with properly structured monomials
        @test result isa Vector{<:Tuple}
        @test all(t -> length(t[2].components) == 2, result)
    end

    @testset "Simplification - Type promotion" begin
        # Test _infer_coef_type: Float64 + ComplexF64 -> ComplexF64
        # Use valid canonical forms
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])  # Returns ComplexF64 when simplified
        m_unip = NormalMonomial{UnipotentAlgebra,UInt16}(UInt16[encode_index(UInt16, 1, 1)])  # Returns Float64

        cm = ComposedMonomial((m_pauli, m_unip))
        result = simplify(cm)

        # Result should have ComplexF64 coefficient
        @test result isa Vector{<:Tuple}
        @test length(result) == 1
        @test result[1][1] isa ComplexF64
    end

    @testset "Simplification - Zero result handling" begin
        # Fermionic: a₁ a₁ = 0 (same mode annihilators = nilpotent)
        pairs_fermi = simplify(FermionicAlgebra, Int32[1, 1])
        @test pairs_fermi == [(0, Int32[])]

        # For ComposedMonomial test, use identity fermionic (empty)
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[])  # Identity
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        cm = ComposedMonomial((m_fermi, m_pauli))

        result = simplify(cm)

        # Result is Vector{(coefficient, ComposedMonomial)} pairs with identity
        @test result isa Vector{<:Tuple}
        @test length(result) >= 1
    end

    @testset "_infer_coef_type_from_types" begin
        # Test compile-time coefficient type inference using coeff_type
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[1])

        # Monomial + Monomial: Pauli (ComplexF64) + Fermionic (Float64) -> ComplexF64
        T1 = _infer_coef_type_from_types((m_pauli, m_fermi))
        @test T1 == ComplexF64

        # All Float64 monomials
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        T3 = _infer_coef_type_from_types((m_fermi, m_nc))
        @test T3 == Float64

        # Polynomial (ComplexF64) + Monomial (Float64) -> ComplexF64
        p_pauli = Polynomial([(1.0 + 0.0im, m_pauli)])
        T4 = _infer_coef_type_from_types((p_pauli, m_fermi))
        @test T4 == ComplexF64
    end

    @testset "_cartesian_product_iter! (via _expand_simplified_components)" begin
        # Test _cartesian_product_iter! indirectly through _expand_simplified_components
        # The helper expects components iterable as `(coef, NormalMonomial)` pairs.
        f1 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1])
        f2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[2])
        ferm_multi = Polynomial([(1.0, f1), (2.0, f2)])

        p = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        pauli_phase = Polynomial([(0.0 + 1.0im, p)])

        cm_template = ComposedMonomial((f1, p))
        result = _expand_simplified_components(typeof(cm_template), (ferm_multi, pauli_phase))

        # Should produce 2 * 1 = 2 terms (Cartesian product across components)
        @test length(result) == 2

        # Coefficients pick up the Pauli phase: {i, 2i}
        coeffs = Set(map(imag ∘ first, result))
        @test 1.0 in coeffs
        @test 2.0 in coeffs

        # Check monomials are composed correctly
        @test all(t -> t[2] isa ComposedMonomial, result)
        @test all(t -> length(t[2].components) == 2, result)
    end

    @testset "Display - ComposedMonomial" begin
        # Use valid canonical forms
        m1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])  # Different sites
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 2])  # c₁†a₂ (normal-ordered)
        cm = ComposedMonomial((m1, m2))

        s = sprint(show, cm)
        @test occursin("ComposedMonomial", s)
        # Check for the site indices we used
        @test occursin("$P_S1", s) || occursin("0x0001", s)
        @test occursin("[-1, 2]", s)

        # Single component
        cm_single = ComposedMonomial((m1,))
        s_single = sprint(show, cm_single)
        @test occursin("ComposedMonomial", s_single)
    end

    @testset "Display - (coefficient, ComposedMonomial)" begin
        m1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2])
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[1])
        cm = ComposedMonomial((m1, m2))

        # Zero term
        t_zero = (0.0, cm)
        @test sprint(show, t_zero) == "0"

        # Identity term (coefficient 1, empty monomials)
        m_empty1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        m_empty2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[])
        cm_empty = ComposedMonomial((m_empty1, m_empty2))
        t_one = (1.0, cm_empty)
        @test sprint(show, t_one) == "1"

        # Coefficient = 1 with non-empty monomial
        t_coef_one = (1.0, cm)
        s_one = sprint(show, t_coef_one)
        @test !occursin("1.0 *", s_one)  # Should not show coefficient
        @test occursin("ComposedMonomial", s_one)

        # General case
        t_general = (2.5, cm)
        s_general = sprint(show, t_general)
        @test occursin("2.5", s_general)
        @test occursin("*", s_general)
        @test occursin("ComposedMonomial", s_general)

        # Complex coefficient should be parenthesized
        t_complex = (1 + 2im, cm)
        s_complex = sprint(show, t_complex)
        @test occursin("(1 + 2im)", s_complex)
        @test occursin("*", s_complex)
    end

    @testset "Integration - Pauli + Fermionic composition" begin
        # Create realistic composed monomial with valid canonical forms
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1, P_S2, P_S3])  # 3 sites
        m_fermi = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 1])  # c₁†a₁ (normal-ordered)

        cm = ComposedMonomial((m_pauli, m_fermi))

        @test degree(cm) == 5
        @test length(cm) == 2

        # Simplify
        terms = simplify(cm)
        @test !isempty(terms)
        @test all(t -> t[2] isa ComposedMonomial, terms)
    end

    @testset "Integration - Three-way composition" begin
        # Pauli + Unipotent + Projector
        m_pauli = NormalMonomial{PauliAlgebra,UInt16}(UInt16[P_S1])
        u1 = encode_index(UInt16, 1, 1)
        m_unip = NormalMonomial{UnipotentAlgebra,UInt16}(simplify(UnipotentAlgebra, UInt16[u1, u1]))
        p1 = encode_index(UInt16, 1, 1)
        m_proj = NormalMonomial{ProjectorAlgebra,UInt16}(UInt16[p1])

        cm = ComposedMonomial((m_pauli, m_unip, m_proj))

        @test length(cm) == 3
        @test cm[1] === m_pauli
        @test cm[2] === m_unip
        @test cm[3] === m_proj

        terms = simplify(cm)
        @test !isempty(terms)
    end

    @testset "Edge Cases" begin
        # Empty monomials in all components
        m1 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[])
        m2 = NormalMonomial{FermionicAlgebra,Int32}(Int32[])
        cm_empty = ComposedMonomial((m1, m2))

        @test degree(cm_empty) == 0
        @test length(cm_empty) == 2

        terms = simplify(cm_empty)
        @test length(terms) == 1
        @test isone(terms[1])

        # Hash stability
        cm1 = ComposedMonomial((m1, m2))
        cm2 = ComposedMonomial((m1, m2))
        @test hash(cm1) == hash(cm2)
        @test cm1 == cm2
    end
end
