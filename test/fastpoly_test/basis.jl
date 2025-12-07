using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_ncbasis, get_ncbasis_deg

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
end
