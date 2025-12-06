# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    encode_index,
    decode_operator_id,
    decode_site,
    get_ncbasis,
    get_ncbasis_deg,
    has_consecutive_repeats,
    create_noncommutative_variables,
    star

@testset "Utilities" begin
    @testset "Index Encoding/Decoding" begin
        # Test UInt16 encoding (4-bit site, 12-bit operator)
        idx = encode_index(UInt16, 3, 5)  # operator 3, site 5
        @test decode_operator_id(idx) == 3
        @test decode_site(idx) == 5

        # Test UInt8 encoding (2-bit site, 6-bit operator)
        # UInt8 max sites is 3 (2 bits), so use site 2
        idx8 = encode_index(UInt8, 2, 2)  # operator 2, site 2
        @test decode_operator_id(idx8) == 2
        @test decode_site(idx8) == 2

        # Multiple encodings
        idx1 = encode_index(UInt16, 1, 1)
        idx2 = encode_index(UInt16, 1, 2)
        idx3 = encode_index(UInt16, 2, 1)

        @test idx1 != idx2  # Same operator, different site
        @test idx1 != idx3  # Different operator, same site
    end

    @testset "Basis Generation" begin
        # Test get_ncbasis for different algebra types
        basis_nc = get_ncbasis(NonCommutativeAlgebra, 2, 2)
        @test length(basis_nc) == 7  # 1 + 2 + 4

        basis_pauli = get_ncbasis(PauliAlgebra, 2, 2)
        @test length(basis_pauli) == 7

        # Test single degree
        basis_deg1 = get_ncbasis_deg(NonCommutativeAlgebra, 3, 1)
        @test length(basis_deg1) == 3  # n^d = 3^1 = 3

        basis_deg2 = get_ncbasis_deg(NonCommutativeAlgebra, 3, 2)
        @test length(basis_deg2) == 9  # n^d = 3^2 = 9
    end

    @testset "Filtered Basis" begin
        # Unfiltered: all words including [1,1], [2,2]
        basis_unfiltered = get_ncbasis(UnipotentAlgebra, 2, 2)
        @test length(basis_unfiltered) == 7

        # Filtered: removes consecutive repeats
        basis_filtered = get_ncbasis(UnipotentAlgebra, 2, 2; filter_constraint=true)
        @test length(basis_filtered) == 5  # 1 + 2 + 2 (no [1,1] or [2,2])

        # Check that filtered basis has no consecutive repeats
        for m in basis_filtered
            @test !has_consecutive_repeats(m.word)
        end
    end

    @testset "has_consecutive_repeats" begin
        @test !has_consecutive_repeats([1, 2, 3])
        @test has_consecutive_repeats([1, 1, 2])
        @test has_consecutive_repeats([1, 2, 2])
        @test has_consecutive_repeats([1, 1])
        @test !has_consecutive_repeats([1])
        @test !has_consecutive_repeats(Int[])
    end

    @testset "Basis Ordering" begin
        basis = get_ncbasis(NonCommutativeAlgebra, 2, 3)

        # Verify basis is sorted
        @test issorted(basis)

        # Verify first element is identity
        @test isone(basis[1])

        # Verify degree ordering: deg 0 < deg 1 < deg 2 < deg 3
        degrees = [degree(m) for m in basis]
        @test issorted(degrees)
    end

    @testset "Integer Type Parameter" begin
        # Test with different integer types
        basis_int = get_ncbasis(NonCommutativeAlgebra, 2, 2; T=Int64)
        basis_int32 = get_ncbasis(NonCommutativeAlgebra, 2, 2; T=Int32)
        basis_int16 = get_ncbasis(NonCommutativeAlgebra, 2, 2; T=Int16)

        @test length(basis_int) == length(basis_int32) == length(basis_int16)

        # Check element types
        @test eltype(basis_int[1].word) == Int64
        @test eltype(basis_int32[1].word) == Int32
        @test eltype(basis_int16[1].word) == Int16
    end

    @testset "Empty Basis Edge Cases" begin
        # Degree 0 returns only identity
        basis_deg0 = get_ncbasis(NonCommutativeAlgebra, 3, 0)
        @test length(basis_deg0) == 1
        @test isone(basis_deg0[1])

        # Zero variables, any degree
        basis_0vars = get_ncbasis_deg(NonCommutativeAlgebra, 0, 2)
        @test isempty(basis_0vars)
    end

    @testset "Large Basis Computation" begin
        # Moderate size test to ensure performance is reasonable
        basis_large = get_ncbasis(NonCommutativeAlgebra, 3, 3)
        # 1 + 3 + 9 + 27 = 40
        @test length(basis_large) == 40
    end

    @testset "Variable Registry Integration" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:5)])

        # Variables created should have correct degrees
        for i in 1:5
            @test degree(x[i]) == 1
        end

        # Monomials from registry should multiply correctly
        m = x[1] * x[2]
        @test m isa Term
        @test degree(m.monomial) == 2
    end

    @testset "Star Operation" begin
        # Star operation reverses word for unsigned types (self-adjoint)
        m = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        m_star = star(m)
        @test m_star.word == [3, 2, 1]

        # Star is involution
        @test star(star(m)) == m

        # Empty monomial star
        m_empty = Monomial{NonCommutativeAlgebra}(UInt8[])
        @test isone(star(m_empty))
    end
end
