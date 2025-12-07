using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_ncbasis, get_ncbasis_deg, has_consecutive_repeats, _generate_words, _generate_all_words, index_type

@testset "Basis Generation" begin
    @testset "get_ncbasis" begin
        # Test basis generation for non-commutative algebra
        basis_deg2 = get_ncbasis(NonCommutativeAlgebra, 3, 2)

        # Basis up to degree 2 with 3 variables: 1 + 3 + 9 = 13
        @test length(basis_deg2) == 13

        # First element should be identity
        @test isone(basis_deg2[1])

        # Check that basis is sorted
        @test issorted(basis_deg2)

        # Test specific degrees
        deg1_only = get_ncbasis(NonCommutativeAlgebra, 3, 1)
        @test length(deg1_only) == 4  # 1 + 3
    end

    @testset "get_ncbasis with different algebras" begin
        # Pauli algebra
        pauli_basis = get_ncbasis(PauliAlgebra, 2, 2)
        @test all(m -> m isa Monomial{PauliAlgebra}, pauli_basis)
        @test issorted(pauli_basis)

        # Projector algebra
        proj_basis = get_ncbasis(ProjectorAlgebra, 3, 1)
        @test all(m -> m isa Monomial{ProjectorAlgebra}, proj_basis)
    end

    @testset "get_ncbasis_deg (exact degree)" begin
        # Exact degree 2 with 3 variables: 3^2 = 9 monomials
        basis_exact_2 = get_ncbasis_deg(NonCommutativeAlgebra, 3, 2)
        @test length(basis_exact_2) == 9
        @test all(m -> degree(m) == 2, basis_exact_2)

        # Exact degree 0 is just identity
        basis_exact_0 = get_ncbasis_deg(NonCommutativeAlgebra, 3, 0)
        @test length(basis_exact_0) == 1
        @test isone(basis_exact_0[1])

        # Exact degree 1
        basis_exact_1 = get_ncbasis_deg(NonCommutativeAlgebra, 3, 1)
        @test length(basis_exact_1) == 3
        @test all(m -> degree(m) == 1, basis_exact_1)
    end

    @testset "has_consecutive_repeats" begin
        # No consecutive repeats
        @test has_consecutive_repeats([1, 2, 3]) == false
        @test has_consecutive_repeats([1, 2, 1]) == false  # non-consecutive repeat
        @test has_consecutive_repeats([1, 2, 3, 1]) == false

        # Has consecutive repeats
        @test has_consecutive_repeats([1, 1, 2]) == true
        @test has_consecutive_repeats([1, 2, 2]) == true
        @test has_consecutive_repeats([1, 1]) == true
        @test has_consecutive_repeats([2, 2, 2]) == true

        # Edge cases
        @test has_consecutive_repeats([]) == false  # empty
        @test has_consecutive_repeats([1]) == false  # single element
    end

    @testset "_generate_words" begin
        # Basic generation
        words = _generate_words(2, 2, Int)
        @test length(words) == 4  # 2^2 = 4
        @test [1, 1] in words
        @test [1, 2] in words
        @test [2, 1] in words
        @test [2, 2] in words

        # Degree 0
        words_d0 = _generate_words(3, 0, Int)
        @test length(words_d0) == 1
        @test words_d0[1] == []

        # Degree 1
        words_d1 = _generate_words(3, 1, Int)
        @test length(words_d1) == 3
        @test [1] in words_d1
        @test [2] in words_d1
        @test [3] in words_d1

        # No variables (n=0)
        words_n0 = _generate_words(0, 2, Int)
        @test length(words_n0) == 0
        @test words_n0 == []

        # Type parameter
        words_uint16 = _generate_words(2, 1, UInt16)
        @test eltype(words_uint16[1]) == UInt16
        @test UInt16[1] in words_uint16
    end

    @testset "monomials alias" begin
        using NCTSSoS.FastPolynomials: monomials

        # Should produce same result as get_ncbasis_deg
        basis_deg = get_ncbasis_deg(PauliAlgebra, 3, 2)
        basis_monos = monomials(PauliAlgebra, 3, 2)
        @test basis_deg == basis_monos

        # Works with different algebra types
        basis_nc = monomials(NonCommutativeAlgebra, 2, 2)
        @test length(basis_nc) == 4
        @test all(m -> m isa Monomial{NonCommutativeAlgebra}, basis_nc)

        # Type parameter
        basis_uint = monomials(PauliAlgebra, 2, 1; T=UInt16)
        @test all(m -> eltype(m.word) == UInt16, basis_uint)
    end

    @testset "filter_constraint option" begin
        # get_ncbasis_deg with filter_constraint
        # Without filtering: 2^2 = 4 monomials
        basis_unfiltered = get_ncbasis_deg(NonCommutativeAlgebra, 2, 2; filter_constraint=false)
        @test length(basis_unfiltered) == 4

        # With filtering: removes [1,1] and [2,2]
        basis_filtered = get_ncbasis_deg(NonCommutativeAlgebra, 2, 2; filter_constraint=true)
        @test length(basis_filtered) == 2
        @test all(m -> !has_consecutive_repeats(m.word), basis_filtered)

        # get_ncbasis with filter_constraint
        # Without filtering: 1 + 2 + 4 = 7
        full_unfiltered = get_ncbasis(NonCommutativeAlgebra, 2, 2; filter_constraint=false)
        @test length(full_unfiltered) == 7

        # With filtering: 1 + 2 + 2 = 5 (removes [1,1] and [2,2] from degree 2)
        full_filtered = get_ncbasis(NonCommutativeAlgebra, 2, 2; filter_constraint=true)
        @test length(full_filtered) == 5
        @test all(m -> !has_consecutive_repeats(m.word), full_filtered)

        # Degree 3 with filtering
        deg3_filtered = get_ncbasis_deg(NonCommutativeAlgebra, 2, 3; filter_constraint=true)
        @test all(m -> !has_consecutive_repeats(m.word), deg3_filtered)
        # Unfiltered: 2^3 = 8, filtered should have fewer
        @test length(deg3_filtered) < 8
    end

    @testset "Type parameter T" begin
        # UInt16
        basis_uint16 = get_ncbasis_deg(PauliAlgebra, 2, 2; T=UInt16)
        @test all(m -> eltype(m.word) == UInt16, basis_uint16)
        @test length(basis_uint16) == 4

        # Int32
        basis_int32 = get_ncbasis(NonCommutativeAlgebra, 2, 1; T=Int32)
        @test all(m -> isempty(m.word) || eltype(m.word) == Int32, basis_int32)

        # Type preserved through get_ncbasis
        full_uint16 = get_ncbasis(PauliAlgebra, 2, 2; T=UInt16)
        @test all(m -> isempty(m.word) || eltype(m.word) == UInt16, full_uint16)
    end

    @testset "Edge cases: n=0, negative d" begin
        # Negative degree should return empty
        @test get_ncbasis_deg(NonCommutativeAlgebra, 3, -1) == []
        @test get_ncbasis_deg(PauliAlgebra, 2, -5) == []

        # n=0 with d>0 should return empty (no variables)
        @test get_ncbasis_deg(NonCommutativeAlgebra, 0, 1) == []
        @test get_ncbasis_deg(NonCommutativeAlgebra, 0, 2) == []

        # n=0, d=0 should return identity only
        basis_n0_d0 = get_ncbasis_deg(NonCommutativeAlgebra, 0, 0)
        @test length(basis_n0_d0) == 1
        @test isone(basis_n0_d0[1])

        # get_ncbasis with n=0
        full_n0 = get_ncbasis(NonCommutativeAlgebra, 0, 2)
        @test length(full_n0) == 1  # only identity
        @test isone(full_n0[1])
    end

    # =========================================================================
    # Registry-Based API Tests
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

    @testset "Registry-based get_ncbasis_deg (NonCommutativeAlgebra)" begin
        # Create registry
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Degree 2: returns Vector{Term}
        basis = get_ncbasis_deg(reg, 2)
        @test basis isa Vector{<:Term}
        @test length(basis) == 9  # 3^2 = 9

        # All terms should have degree 2 monomials (NonCommutativeAlgebra doesn't simplify)
        @test all(t -> degree(t.monomial) == 2, basis)

        # Degree 0: identity term
        basis_d0 = get_ncbasis_deg(reg, 0)
        @test length(basis_d0) == 1
        @test isone(basis_d0[1].monomial)

        # Degree 1
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test length(basis_d1) == 3
        @test all(t -> degree(t.monomial) == 1, basis_d1)

        # Negative degree
        basis_neg = get_ncbasis_deg(reg, -1)
        @test isempty(basis_neg)
    end

    @testset "Registry-based get_ncbasis (NonCommutativeAlgebra)" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])

        # Up to degree 2: 1 + 2 + 4 = 7 terms
        basis = get_ncbasis(reg, 2)
        @test basis isa Vector{<:Term}
        @test length(basis) == 7

        # Contains identity
        @test any(t -> isone(t.monomial), basis)

        # Contains degree 1 and 2 terms
        @test any(t -> degree(t.monomial) == 1, basis)
        @test any(t -> degree(t.monomial) == 2, basis)
    end

    @testset "Registry-based get_ncbasis_deg (UnipotentAlgebra)" begin
        # UnipotentAlgebra: U^2 = I (consecutive repeats simplify to identity)
        reg, (U,) = create_unipotent_variables([("U", 1:2)])

        # Degree 2: generates 4 words, but U_i * U_i simplifies
        basis = get_ncbasis_deg(reg, 2)
        @test basis isa Vector{<:Term}

        # Should have some terms (exact count depends on simplification)
        @test length(basis) >= 2  # At least [U1,U2] and [U2,U1]
    end

    @testset "Registry-based get_ncbasis with PauliAlgebra" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        # Degree 1: 6 Pauli operators
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test basis_d1 isa Vector{<:Term}
        @test length(basis_d1) == 6

        # All should be degree 1 after simplification
        @test all(t -> degree(t.monomial) == 1, basis_d1)
    end

    @testset "Registry-based get_ncbasis with FermionicAlgebra" begin
        reg, (a, a_dag) = create_fermionic_variables(1:2)

        # Degree 1: 4 operators (a1, a1†, a2, a2†)
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test basis_d1 isa Vector{<:Term}
        @test length(basis_d1) == 4

        # Degree 2: 16 words, but some simplify (e.g., a_i a_i = 0)
        basis_d2 = get_ncbasis_deg(reg, 2)
        @test basis_d2 isa Vector{<:Term}
        # Fermionic simplification may produce mixed-degree terms
    end

    @testset "Registry-based API type consistency" begin
        # The key benefit: registry indices match monomial indices
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        basis = get_ncbasis_deg(reg, 1)

        # Get indices from registry
        reg_indices = indices(reg)

        # All monomial words should use registry indices
        for term in basis
            if !isempty(term.monomial.word)
                @test all(idx -> idx in reg_indices, term.monomial.word)
            end
        end

        # Index type should match
        T = index_type(reg)
        for term in basis
            if !isempty(term.monomial.word)
                @test eltype(term.monomial.word) == T
            end
        end
    end

    @testset "Registry-based API with multi-prefix variables" begin
        # Test with multiple prefix groups
        reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 3:4)])

        # Should have 4 total variables
        @test length(reg) == 4

        # Degree 1 basis should have 4 terms
        basis_d1 = get_ncbasis_deg(reg, 1)
        @test length(basis_d1) == 4

        # Degree 2: 4^2 = 16 words (before simplification)
        basis_d2 = get_ncbasis_deg(reg, 2)
        @test basis_d2 isa Vector{<:Term}
    end
end
