# NCTSSoS is loaded by parent runtests.jl
# Exported: get_ncbasis, create_noncommutative_variables, create_pauli_variables, create_unipotent_variables, indices, degree, Polynomial, NormalMonomial, NonCommutativeAlgebra
# Internal (not exported): encode_index
using NCTSSoS: encode_index

@testset "Utilities" begin
    @testset "Basis Generation" begin
        # Test get_ncbasis for different algebra types (registry-based API)
        reg_nc, (x,) = create_noncommutative_variables([("x", 1:2)])
        basis_nc = get_ncbasis(reg_nc, 2)
        @test length(basis_nc) == 7  # 1 + 2 + 4

        reg_pauli, _ = create_pauli_variables(1:2)
        basis_pauli = get_ncbasis(reg_pauli, 2)
        @test length(basis_pauli) > 0  # Exact count depends on simplification
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
        m = x[1] * x[2]
        @test m isa Polynomial{NonCommutativeAlgebra,<:Integer,<:Number}
        @test degree(m) == 2
    end

    @testset "Adjoint Operation" begin
        # Adjoint operation reverses word for unsigned types (self-adjoint).
        # Use indices on a single site so site-canonicalization does not reorder.
        idx1_s1 = encode_index(UInt8, 1, 1)
        idx2_s1 = encode_index(UInt8, 2, 1)
        idx3_s1 = encode_index(UInt8, 3, 1)
        m = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[idx1_s1, idx2_s1, idx3_s1])
        m_adj = adjoint(m)
        @test m_adj.word == UInt8[idx3_s1, idx2_s1, idx1_s1]

        # Adjoint is involution
        @test adjoint(adjoint(m)) == m

        # Empty monomial adjoint
        m_empty = NormalMonomial{NonCommutativeAlgebra,UInt8}(UInt8[])
        @test isone(adjoint(m_empty))

        # Julia syntax shorthand
        @test m' == m_adj
    end
end
