using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_ncbasis, get_ncbasis_deg, has_consecutive_repeats, _generate_words

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
end
