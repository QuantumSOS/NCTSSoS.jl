using Test, NCTSSoS
# get_ncbasis is exported from NCTSSoS; get_ncbasis_deg, _generate_all_words are internal
using NCTSSoS: get_ncbasis_deg, _generate_all_words

@testset "Basis Generation" begin

    # =========================================================================
    # _generate_all_words Tests
    # =========================================================================

    @testset "_generate_all_words with arbitrary indices" begin
        # Basic test with standard indices
        words = _generate_all_words([1, 2], 2)
        @test length(words) == 4
        @test [1, 1] in words
        @test [1, 2] in words
        @test [2, 1] in words
        @test [2, 2] in words

        # Non-contiguous indices
        words_nc = _generate_all_words([5, 10], 2)
        @test length(words_nc) == 4
        @test [5, 5] in words_nc
        @test [5, 10] in words_nc
        @test [10, 5] in words_nc
        @test [10, 10] in words_nc

        # Degree 0 returns empty word
        words_d0 = _generate_all_words([1, 2, 3], 0)
        @test length(words_d0) == 1
        @test words_d0[1] == []

        # Degree 1
        words_d1 = _generate_all_words([3, 7, 11], 1)
        @test length(words_d1) == 3
        @test [3] in words_d1
        @test [7] in words_d1
        @test [11] in words_d1

        # Empty indices
        words_empty = _generate_all_words(Int[], 2)
        @test length(words_empty) == 0

        # Signed indices (for fermionic/bosonic)
        words_signed = _generate_all_words([-1, 1], 2)
        @test length(words_signed) == 4
        @test [-1, -1] in words_signed
        @test [-1, 1] in words_signed
        @test [1, -1] in words_signed
        @test [1, 1] in words_signed
    end

    # =========================================================================
    # get_ncbasis_deg Tests (Registry-Based)
    # =========================================================================

    @testset "get_ncbasis_deg (NonCommutativeAlgebra)" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Degree 2: returns Vector{Polynomial}
        basis = get_ncbasis_deg(reg, 2)
        @test basis isa Vector{<:Polynomial}
        @test length(basis) == 9  # 3^2 = 9

        # All polynomials should have single degree-2 monomial (NonCommutativeAlgebra doesn't simplify)
        @test all(p -> length(terms(p)) == 1 && degree(monomials(p)[1]) == 2, basis)

        # Degree 0: identity polynomial
        basis_d0 = get_ncbasis_deg(reg, 0)
        @test length(basis_d0) == 1
        @test isone(basis_d0[1])

        # Degree 1
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test length(basis_d1) == 3
        @test all(p -> length(terms(p)) == 1 && degree(monomials(p)[1]) == 1, basis_d1)

        # Negative degree
        basis_neg = get_ncbasis_deg(reg, -1)
        @test isempty(basis_neg)
    end

    @testset "get_ncbasis_deg (PauliAlgebra)" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        # Degree 1: 6 Pauli operators
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test basis_d1 isa Vector{<:Polynomial}
        @test length(basis_d1) == 6

        # Each polynomial should be a single term with degree 1 monomial
        @test all(p -> length(terms(p)) == 1 && degree(monomials(p)[1]) == 1, basis_d1)
    end

    @testset "get_ncbasis_deg (ProjectorAlgebra)" begin
        reg, (P,) = create_projector_variables([("P", 1:3)])

        # Degree 1: 3 projectors
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test basis_d1 isa Vector{<:Polynomial}
        @test length(basis_d1) == 3
    end

    @testset "get_ncbasis_deg (UnipotentAlgebra)" begin
        # UnipotentAlgebra: U^2 = I (consecutive repeats simplify to identity)
        reg, (U,) = create_unipotent_variables([("U", 1:2)])

        # Degree 2: generates 4 words, but U_i * U_i simplifies
        basis = get_ncbasis_deg(reg, 2)
        @test basis isa Vector{<:Polynomial}

        # Should have 4 polynomials (one per input word)
        @test length(basis) == 4
    end

    @testset "get_ncbasis_deg (FermionicAlgebra)" begin
        reg, (a, a_dag) = create_fermionic_variables(1:2)

        # Degree 1: 4 operators (a1, a1†, a2, a2†)
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test basis_d1 isa Vector{<:Polynomial}
        @test length(basis_d1) == 4

        # Degree 2: 16 words, some may produce multi-term polynomials
        basis_d2 = get_ncbasis_deg(reg, 2)
        @test basis_d2 isa Vector{<:Polynomial}
        @test length(basis_d2) == 16  # One polynomial per input word
    end

    @testset "get_ncbasis_deg (BosonicAlgebra)" begin
        reg, (c, c_dag) = create_bosonic_variables(1:2)

        # Degree 1: 4 operators
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test basis_d1 isa Vector{<:Polynomial}
        @test length(basis_d1) == 4
    end

    # =========================================================================
    # get_ncbasis Tests (Registry-Based)
    # =========================================================================

    @testset "get_ncbasis (NonCommutativeAlgebra)" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])

        # Up to degree 2: 1 + 2 + 4 = 7 polynomials
        basis = get_ncbasis(reg, 2)
        @test basis isa Vector{<:Polynomial}
        @test length(basis) == 7

        # Contains identity polynomial
        @test any(isone, basis)

        # Contains degree 1 and 2 polynomials
        @test any(p -> degree(p) == 1, basis)
        @test any(p -> degree(p) == 2, basis)
    end

    @testset "get_ncbasis (PauliAlgebra)" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        # Up to degree 2
        basis = get_ncbasis(reg, 2)
        @test basis isa Vector{<:Polynomial}

        # Contains identity polynomial
        @test any(isone, basis)
    end

    @testset "get_ncbasis (with multi-prefix variables)" begin
        # Test with multiple prefix groups
        reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 3:4)])

        # Should have 4 total variables
        @test length(reg) == 4

        # Degree 1 basis should have 4 polynomials
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test length(basis_d1) == 4

        # Degree 2: 4^2 = 16 polynomials (one per input word)
        basis_d2 = get_ncbasis_deg(reg, 2)
        @test basis_d2 isa Vector{<:Polynomial}
        @test length(basis_d2) == 16
    end

    # =========================================================================
    # Type Consistency Tests
    # =========================================================================

    @testset "Registry-based API type consistency" begin
        # The key benefit: registry indices match monomial indices
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        basis = get_ncbasis_deg(reg, 1)

        # Get indices from registry
        reg_indices = indices(reg)

        # All monomial words should use registry indices
        for poly in basis
            for mono in monomials(poly)
                if !isempty(mono.word)
                    @test all(idx -> idx in reg_indices, mono.word)
                end
            end
        end

        # Index type should match
        T = eltype(keys(reg.idx_to_variables))
        for poly in basis
            for mono in monomials(poly)
                if !isempty(mono.word)
                    @test eltype(mono.word) == T
                end
            end
        end
    end

    # =========================================================================
    # Edge Cases
    # =========================================================================

    @testset "Edge cases" begin
        # Single variable
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        basis = get_ncbasis(reg, 2)
        @test length(basis) == 3  # 1 + 1 + 1

        # Large degree
        reg2, (y,) = create_noncommutative_variables([("y", 1:2)])
        basis_d3 = get_ncbasis_deg(reg2, 3)
        @test length(basis_d3) == 8  # 2^3 = 8
    end
end


# TODO: it does not test for different physical group
