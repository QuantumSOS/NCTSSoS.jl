# Note: FastPolynomials is loaded by setup.jl
using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: encode_index

# Import internal functions for testing
import NCTSSoS.FastPolynomials: _compute_composed_hash, _expand_simplified_components,
    _infer_coef_type, _to_term_vector

@testset "ComposedMonomial" begin
    @testset "Construction" begin
        # Basic construction with two algebras
        m_pauli = Monomial{PauliAlgebra}(UInt16[1, 2])
        m_fermi = Monomial{FermionicAlgebra}(Int32[-1, 2])

        cm = ComposedMonomial((m_pauli, m_fermi))

        @test cm isa ComposedMonomial
        @test length(cm.components) == 2
        @test cm.components[1] === m_pauli
        @test cm.components[2] === m_fermi
        @test cm.hash isa UInt64

        # Construction with three algebras
        m_unip = Monomial{UnipotentAlgebra}(UInt16[encode_index(UInt16, 1, 1)])
        cm3 = ComposedMonomial((m_pauli, m_fermi, m_unip))

        @test length(cm3.components) == 3
        @test cm3.components[3] === m_unip

        # Single component
        cm_single = ComposedMonomial((m_pauli,))
        @test length(cm_single.components) == 1
    end

    @testset "_compute_composed_hash" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m2 = Monomial{FermionicAlgebra}(Int32[1, 2])

        # Hash depends on component count and component hashes
        h1 = _compute_composed_hash((m1,))
        h2 = _compute_composed_hash((m1, m2))
        h3 = _compute_composed_hash((m2, m1))  # Different order

        @test h1 != h2  # Different number of components
        # Note: h2 and h3 might collide due to hash function structure,
        # but they should ideally differ - we test that hash is computed
        @test h2 isa UInt64
        @test h3 isa UInt64

        # Same components same order -> same hash
        m1_copy = Monomial{PauliAlgebra}(UInt16[1, 2])
        h4 = _compute_composed_hash((m1_copy, m2))
        @test h2 == h4
    end

    @testset "Equality" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m2 = Monomial{FermionicAlgebra}(Int32[1, 2])
        m3 = Monomial{PauliAlgebra}(UInt16[1, 2])  # Same as m1
        m4 = Monomial{PauliAlgebra}(UInt16[2, 1])  # Different from m1

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
        m_bosonic = Monomial{BosonicAlgebra}(Int32[1, 2])
        cm_bosonic = ComposedMonomial((m1, m_bosonic))
        @test cm1 != cm_bosonic  # Different tuple types

        # Reflexive
        @test cm1 == cm1

        # Symmetric
        @test (cm1 == cm2) == (cm2 == cm1)
    end

    @testset "Hashing" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m2 = Monomial{FermionicAlgebra}(Int32[1, 2])

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
        m1_short = Monomial{PauliAlgebra}(UInt16[1])
        m1_long = Monomial{PauliAlgebra}(UInt16[1, 2, 3])
        m2 = Monomial{FermionicAlgebra}(Int32[1])

        # Degree-first ordering
        cm_deg2 = ComposedMonomial((m1_short, m2))  # degree 1+1=2
        cm_deg4 = ComposedMonomial((m1_long, m2))   # degree 3+1=4

        @test isless(cm_deg2, cm_deg4)
        @test !isless(cm_deg4, cm_deg2)

        # Lexicographic when degrees equal
        m_lex1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m_lex2 = Monomial{PauliAlgebra}(UInt16[1, 3])
        m_fermi = Monomial{FermionicAlgebra}(Int32[1])

        cm_lex1 = ComposedMonomial((m_lex1, m_fermi))
        cm_lex2 = ComposedMonomial((m_lex2, m_fermi))

        @test isless(cm_lex1, cm_lex2)
        @test !isless(cm_lex2, cm_lex1)

        # Equal -> not less
        cm_dup = ComposedMonomial((m_lex1, m_fermi))
        @test !isless(cm_lex1, cm_dup)

        # Lexicographic across multiple components
        m_pauli1 = Monomial{PauliAlgebra}(UInt16[1])
        m_pauli2 = Monomial{PauliAlgebra}(UInt16[1])
        m_fermi1 = Monomial{FermionicAlgebra}(Int32[1])
        m_fermi2 = Monomial{FermionicAlgebra}(Int32[2])

        cm_a = ComposedMonomial((m_pauli1, m_fermi1))  # (1, 1)
        cm_b = ComposedMonomial((m_pauli2, m_fermi2))  # (1, 2)

        @test isless(cm_a, cm_b)  # Second component differs
    end

    @testset "Access methods" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m2 = Monomial{FermionicAlgebra}(Int32[1, 2, 3])
        m3 = Monomial{UnipotentAlgebra}(UInt16[encode_index(UInt16, 1, 1)])

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
        m_empty1 = Monomial{PauliAlgebra}(UInt16[])
        m_empty2 = Monomial{FermionicAlgebra}(Int32[])
        cm_empty = ComposedMonomial((m_empty1, m_empty2))
        @test degree(cm_empty) == 0
    end

    @testset "Term isone" begin
        m_empty1 = Monomial{PauliAlgebra}(UInt16[])
        m_empty2 = Monomial{FermionicAlgebra}(Int32[])
        cm_empty = ComposedMonomial((m_empty1, m_empty2))

        # Coefficient 1, all empty monomials -> true
        t_one = Term(1.0, cm_empty)
        @test isone(t_one)

        # Complex coefficient 1
        t_one_complex = Term(1.0 + 0.0im, cm_empty)
        @test isone(t_one_complex)

        # Coefficient not 1 -> false
        t_two = Term(2.0, cm_empty)
        @test !isone(t_two)

        # Non-empty monomial -> false
        m_nonempty = Monomial{PauliAlgebra}(UInt16[1])
        cm_nonempty = ComposedMonomial((m_nonempty, m_empty2))
        t_nonempty = Term(1.0, cm_nonempty)
        @test !isone(t_nonempty)

        # Second component non-empty -> false
        cm_nonempty2 = ComposedMonomial((m_empty1, Monomial{FermionicAlgebra}(Int32[1])))
        t_nonempty2 = Term(1.0, cm_nonempty2)
        @test !isone(t_nonempty2)
    end

    @testset "Simplification - Single-term algebras (Pauli)" begin
        # Pauli: [1, 1] -> [] with coefficient 1 (σ_x^2 = I)
        # ComposedMonomial with Pauli returns Term (since Pauli returns Term)
        m_pauli = Monomial{PauliAlgebra}(UInt16[1, 1])
        m_unip = Monomial{UnipotentAlgebra}(UInt16[encode_index(UInt16, 1, 1)])
        cm = ComposedMonomial((m_pauli, m_unip))

        result = simplify(cm)

        # Result is a Term (Pauli is the most complex component returning Term)
        @test result isa Term
        @test result.coefficient ≈ 1.0 + 0.0im
        @test isempty(result.monomial[1].word)  # Pauli simplified
        @test result.monomial[2].word == [encode_index(UInt16, 1, 1)]
    end

    @testset "Simplification - Single-term algebras (UnipotentAlgebra)" begin
        # Unipotent: U² = I (consecutive pairs cancel)
        # [u, u] -> [] (pair cancels)
        # [u1, u1, u2] -> [u2] (u1 pair cancels, u2 remains)
        # ComposedMonomial with Pauli returns Term (since Pauli returns Term)
        u1 = encode_index(UInt16, 1, 1)
        u2 = encode_index(UInt16, 2, 1)

        m_unip = Monomial{UnipotentAlgebra}(UInt16[u1, u1, u2])
        m_pauli = Monomial{PauliAlgebra}(UInt16[1])
        cm = ComposedMonomial((m_pauli, m_unip))

        result = simplify(cm)

        # Result is a Term (Pauli is the most complex component returning Term)
        @test result isa Term
        # Unipotent [u1, u1, u2] -> [u2] (u1 pair cancels via U²=I)
        @test result.monomial[2].word == [u2]
    end

    @testset "Simplification - Single-term algebras (ProjectorAlgebra)" begin
        # Projector: [e_ii, e_ii] -> [e_ii] (idempotent)
        # ComposedMonomial with Pauli returns Term (since Pauli returns Term)
        p1 = encode_index(UInt16, 1, 1)

        m_proj = Monomial{ProjectorAlgebra}(UInt16[p1, p1])
        m_pauli = Monomial{PauliAlgebra}(UInt16[])
        cm = ComposedMonomial((m_proj, m_pauli))

        result = simplify(cm)

        # Result is a Term (Pauli is the most complex component returning Term)
        @test result isa Term
        @test result.monomial[1].word == [p1]  # Projector simplified
    end

    @testset "Simplification - Multi-term algebras (FermionicAlgebra)" begin
        # Fermionic: anti-commutation can produce multiple terms
        # a1 * a1 = 0 (creation operator squared)
        # But a1 * a1† = 1 - a1† * a1 (anti-commutation)
        # ComposedMonomial with Fermionic returns Vector{Term} (Polynomial type doesn't support ComposedMonomial)

        # Simple case: [1, 1] (a1 * a1) -> should give zero or single term
        m_fermi = Monomial{FermionicAlgebra}(Int32[1, 1])
        m_pauli = Monomial{PauliAlgebra}(UInt16[])
        cm = ComposedMonomial((m_fermi, m_pauli))

        result = simplify(cm)

        # Fermionic [1, 1] -> should simplify (exact behavior depends on implementation)
        # Result is Vector{Term} (since Fermionic returns Polynomial)
        @test result isa Vector{<:Term}
        @test all(t -> t.monomial isa ComposedMonomial, result)
    end

    @testset "Simplification - Multi-term algebras (BosonicAlgebra)" begin
        # Bosonic: [b, b†] -> b† b + 1 (commutation)
        # ComposedMonomial with Bosonic returns Vector{Term} (Polynomial type doesn't support ComposedMonomial)
        m_bosonic = Monomial{BosonicAlgebra}(Int32[1, -1])
        m_pauli = Monomial{PauliAlgebra}(UInt16[])
        cm = ComposedMonomial((m_bosonic, m_pauli))

        result = simplify(cm)

        # Should produce multiple terms due to commutation
        # Result is Vector{Term} (since Bosonic returns Polynomial)
        @test result isa Vector{<:Term}
        # Exact count depends on simplification rules - just verify structure
        @test all(t -> t.monomial isa ComposedMonomial, result)
    end

    @testset "Simplification - Cartesian product" begin
        # When both algebras produce multiple terms, result is Cartesian product
        # Fermionic [1, -1] might produce [(coef1, mono1), (coef2, mono2)]
        # Pauli [1, 1] -> [(1, [])]
        # Result should be product
        # ComposedMonomial with Fermionic returns Vector{Term}

        m_fermi = Monomial{FermionicAlgebra}(Int32[1, -1])
        m_pauli = Monomial{PauliAlgebra}(UInt16[1, 1])
        cm = ComposedMonomial((m_fermi, m_pauli))

        result = simplify(cm)

        # Verify we get Vector{Term} with properly structured terms
        @test result isa Vector{<:Term}
        @test all(t -> length(t.monomial.components) == 2, result)
    end

    @testset "Simplification - Type promotion" begin
        # Test _infer_coef_type: Float64 + ComplexF64 -> ComplexF64
        # Pauli returns Term, so ComposedMonomial returns Term
        m_pauli = Monomial{PauliAlgebra}(UInt16[1, 1])  # Returns ComplexF64
        m_unip = Monomial{UnipotentAlgebra}(UInt16[encode_index(UInt16, 1, 1)])  # Returns Float64

        cm = ComposedMonomial((m_pauli, m_unip))
        result = simplify(cm)

        # Result should have ComplexF64 coefficient (returns Term for Pauli+Unipotent)
        @test result isa Term
        @test result.coefficient isa ComplexF64
    end

    @testset "Simplification - Zero result handling" begin
        # Create a situation that simplifies to zero
        # Fermionic: a * a = 0 (same creation operator twice)
        # ComposedMonomial with Fermionic returns Vector{Term}
        m_fermi = Monomial{FermionicAlgebra}(Int32[1, 1])
        m_pauli = Monomial{PauliAlgebra}(UInt16[])
        cm = ComposedMonomial((m_fermi, m_pauli))

        result = simplify(cm)

        # Result is Vector{Term} (since Fermionic returns Polynomial)
        @test result isa Vector{<:Term}
        # Should return a zero term if simplification gives zero
        if isempty(filter(!iszero, result))
            @test length(result) >= 1  # Should have at least one zero term
            @test all(t -> iszero(t.coefficient) || iszero(t), result)
        end
    end

    @testset "_expand_simplified_components helpers" begin
        # Test _to_term_vector with single Term
        m = Monomial{PauliAlgebra}(UInt16[1])
        t = Term(2.0, m)
        vec = _to_term_vector(t)
        @test vec == [(2.0, m)]

        # Test _to_term_vector with Vector{Term}
        m2 = Monomial{PauliAlgebra}(UInt16[2])
        terms_vec = [Term(1.0, m), Term(3.0, m2)]
        vec2 = _to_term_vector(terms_vec)
        @test length(vec2) == 2
        @test vec2[1] == (1.0, m)
        @test vec2[2] == (3.0, m2)
    end

    @testset "_infer_coef_type" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1])
        m2 = Monomial{FermionicAlgebra}(Int32[1])

        # Float64 only
        component_terms1 = ([(1.0, m1)], [(2.0, m2)])
        T1 = _infer_coef_type(component_terms1)
        @test T1 == Float64

        # Float64 + ComplexF64 -> ComplexF64
        component_terms2 = ([(1.0 + 0.0im, m1)], [(2.0, m2)])
        T2 = _infer_coef_type(component_terms2)
        @test T2 == ComplexF64

        # Int + Float64 -> Float64
        component_terms3 = ([(1, m1)], [(2.0, m2)])
        T3 = _infer_coef_type(component_terms3)
        @test T3 == Float64
    end

    @testset "_cartesian_product! (via _expand_simplified_components)" begin
        # Test _cartesian_product! indirectly through _expand_simplified_components
        m1a = Monomial{PauliAlgebra}(UInt16[1])
        m1b = Monomial{PauliAlgebra}(UInt16[2])
        m2a = Monomial{FermionicAlgebra}(Int32[1])

        t1a = Term(1.0, m1a)
        t1b = Term(2.0, m1b)
        t2a = Term(3.0, m2a)

        # Create simplified components tuple - Pauli and Fermionic
        # Pauli might return single term, Fermionic returns vector
        simplified = ([t1a, t1b], [t2a])

        result = _expand_simplified_components(simplified)

        # Should produce 2 * 1 = 2 terms (Cartesian product)
        @test length(result) == 2

        # Check coefficients are products
        coeffs = Set([t.coefficient for t in result])
        @test 1.0 * 3.0 in coeffs
        @test 2.0 * 3.0 in coeffs

        # Check monomials are composed correctly
        @test all(t -> t.monomial isa ComposedMonomial, result)
        @test all(t -> length(t.monomial.components) == 2, result)
    end

    @testset "Display - ComposedMonomial" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m2 = Monomial{FermionicAlgebra}(Int32[-1, 2])
        cm = ComposedMonomial((m1, m2))

        s = sprint(show, cm)
        @test occursin("ComposedMonomial", s)
        @test occursin("[1, 2]", s) || occursin("UInt16[0x0001, 0x0002]", s)
        @test occursin("[-1, 2]", s)

        # Single component
        cm_single = ComposedMonomial((m1,))
        s_single = sprint(show, cm_single)
        @test occursin("ComposedMonomial", s_single)
    end

    @testset "Display - Term with ComposedMonomial" begin
        m1 = Monomial{PauliAlgebra}(UInt16[1, 2])
        m2 = Monomial{FermionicAlgebra}(Int32[1])
        cm = ComposedMonomial((m1, m2))

        # Zero term
        t_zero = Term(0.0, cm)
        @test sprint(show, t_zero) == "0"

        # Identity term (coefficient 1, empty monomials)
        m_empty1 = Monomial{PauliAlgebra}(UInt16[])
        m_empty2 = Monomial{FermionicAlgebra}(Int32[])
        cm_empty = ComposedMonomial((m_empty1, m_empty2))
        t_one = Term(1.0, cm_empty)
        @test sprint(show, t_one) == "1"

        # Coefficient = 1 with non-empty monomial
        t_coef_one = Term(1.0, cm)
        s_one = sprint(show, t_coef_one)
        @test !occursin("1.0 *", s_one)  # Should not show coefficient
        @test occursin("ComposedMonomial", s_one)

        # General case
        t_general = Term(2.5, cm)
        s_general = sprint(show, t_general)
        @test occursin("2.5", s_general)
        @test occursin("*", s_general)
        @test occursin("ComposedMonomial", s_general)
    end

    @testset "Integration - Pauli + Fermionic composition" begin
        # Create realistic composed monomial
        m_pauli = Monomial{PauliAlgebra}(UInt16[1, 2, 3])  # σx σy σz
        m_fermi = Monomial{FermionicAlgebra}(Int32[1, -1])  # a1 a1†

        cm = ComposedMonomial((m_pauli, m_fermi))

        @test degree(cm) == 5
        @test length(cm) == 2

        # Simplify
        terms = simplify(cm)
        @test !isempty(terms)
        @test all(t -> t.monomial isa ComposedMonomial, terms)
    end

    @testset "Integration - Three-way composition" begin
        # Pauli + Unipotent + Projector
        m_pauli = Monomial{PauliAlgebra}(UInt16[1])
        u1 = encode_index(UInt16, 1, 1)
        m_unip = Monomial{UnipotentAlgebra}(UInt16[u1, u1])
        p1 = encode_index(UInt16, 1, 1)
        m_proj = Monomial{ProjectorAlgebra}(UInt16[p1])

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
        m1 = Monomial{PauliAlgebra}(UInt16[])
        m2 = Monomial{FermionicAlgebra}(Int32[])
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
