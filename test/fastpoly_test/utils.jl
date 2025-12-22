# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    encode_index,
    decode_operator_id,
    decode_site,
    get_ncbasis,
    get_ncbasis_deg,
    create_noncommutative_variables,
    create_pauli_variables,
    create_unipotent_variables,
    indices,
    degree,
    Term,
    Polynomial,
    Monomial,
    NonCommutativeAlgebra

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
        # Test get_ncbasis for different algebra types (registry-based API)
        reg_nc, (x,) = create_noncommutative_variables([("x", 1:2)])
        basis_nc = get_ncbasis(reg_nc, 2)
        @test length(basis_nc) == 7  # 1 + 2 + 4

        reg_pauli, _ = create_pauli_variables(1:2)
        basis_pauli = get_ncbasis(reg_pauli, 2)
        @test length(basis_pauli) > 0  # Exact count depends on simplification

        # Test single degree
        reg3, (y,) = create_noncommutative_variables([("y", 1:3)])
        basis_deg1 = get_ncbasis_deg(reg3, 1)
        @test length(basis_deg1) == 3  # n^d = 3^1 = 3

        basis_deg2 = get_ncbasis_deg(reg3, 2)
        @test length(basis_deg2) == 9  # n^d = 3^2 = 9
    end

    @testset "Basis Ordering" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        basis = get_ncbasis(reg, 3)

        # Verify first element is identity (degree 0)
        @test isone(basis[1])

        # Verify degree ordering: deg 0 < deg 1 < deg 2 < deg 3
        degrees = [degree(p) for p in basis]
        @test issorted(degrees)
    end

    @testset "Empty Basis Edge Cases" begin
        # Degree 0 returns only identity
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        basis_deg0 = get_ncbasis(reg, 0)
        @test length(basis_deg0) == 1
        @test isone(basis_deg0[1])
    end

    @testset "Large Basis Computation" begin
        # Moderate size test to ensure performance is reasonable
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        basis_large = get_ncbasis(reg, 3)
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
        # Note: Monomial multiplication now returns Monomial (word concatenation)
        m = x[1] * x[2]
        @test m isa Monomial
        @test degree(m) == 2
    end

    @testset "Adjoint Operation" begin
        # Adjoint operation reverses word for unsigned types (self-adjoint)
        m = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        m_adj = adjoint(m)
        @test m_adj.word == [3, 2, 1]

        # Adjoint is involution
        @test adjoint(adjoint(m)) == m

        # Empty monomial adjoint
        m_empty = Monomial{NonCommutativeAlgebra}(UInt8[])
        @test isone(adjoint(m_empty))

        # Julia syntax shorthand
        @test m' == m_adj
    end
end
