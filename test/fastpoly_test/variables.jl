# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    VariableRegistry,
    create_noncommutative_variables,
    create_pauli_variables,
    create_projector_variables,
    create_unipotent_variables,
    create_fermionic_variables,
    create_bosonic_variables,
    symbols,
    indices,
    get_ncbasis

@testset "Variable Registry" begin
    @testset "Non-Commutative Variable Creation" begin
        # Create variables with single prefix
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        @test reg isa VariableRegistry
        @test length(x) == 3

        # Check that monomials are created correctly
        @test x[1] isa Monomial{NonCommutativeAlgebra}
        @test x[2] isa Monomial{NonCommutativeAlgebra}
        @test x[3] isa Monomial{NonCommutativeAlgebra}

        # Check registry contains symbols
        @test :x₁ in reg
        @test :x₂ in reg
        @test :x₃ in reg
    end

    @testset "Multi-Prefix Variable Creation" begin
        reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 3:4)])

        @test length(x) == 2
        @test length(y) == 2

        @test :x₁ in reg
        @test :x₂ in reg
        @test :y₃ in reg
        @test :y₄ in reg
    end

    @testset "Pauli Variable Creation" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        @test length(σx) == 2
        @test length(σy) == 2
        @test length(σz) == 2

        @test σx[1] isa Monomial{PauliAlgebra}
        @test σy[1] isa Monomial{PauliAlgebra}
        @test σz[1] isa Monomial{PauliAlgebra}

        # Check registry symbols
        @test :σx₁ in reg
        @test :σy₁ in reg
        @test :σz₁ in reg
        @test :σx₂ in reg
        @test :σy₂ in reg
        @test :σz₂ in reg
    end

    @testset "Projector Variable Creation" begin
        reg, (P,) = create_projector_variables([("P", 1:3)])

        @test length(P) == 3
        @test P[1] isa Monomial{ProjectorAlgebra}

        @test :P₁ in reg
        @test :P₂ in reg
        @test :P₃ in reg
    end

    @testset "Unipotent Variable Creation" begin
        reg, (U,) = create_unipotent_variables([("U", 1:3)])

        @test length(U) == 3
        @test U[1] isa Monomial{UnipotentAlgebra}

        @test :U₁ in reg
        @test :U₂ in reg
        @test :U₃ in reg
    end

    @testset "Fermionic Variable Creation" begin
        reg, (a, a_dag) = create_fermionic_variables(1:2)

        @test length(a) == 2
        @test length(a_dag) == 2
        @test a[1] isa Monomial{FermionicAlgebra}
        @test a_dag[1] isa Monomial{FermionicAlgebra}

        # Check both annihilation and creation operators
        @test :a₁ in reg
        @test Symbol("a⁺₁") in reg
        @test :a₂ in reg
        @test Symbol("a⁺₂") in reg
    end

    @testset "Bosonic Variable Creation" begin
        reg, (c, c_dag) = create_bosonic_variables(1:2)

        @test length(c) == 2
        @test length(c_dag) == 2
        @test c[1] isa Monomial{BosonicAlgebra}
        @test c_dag[1] isa Monomial{BosonicAlgebra}

        # Check both annihilation and creation operators
        @test :c₁ in reg
        @test Symbol("c⁺₁") in reg
    end

    @testset "Registry Access" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:5)])

        # Test symbols function
        syms = symbols(reg)
        @test length(syms) == 5
        @test :x₁ in syms
        @test :x₅ in syms

        # Test indices function
        idxs = indices(reg)
        @test length(idxs) == 5
    end

    @testset "Monomial Power" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # x^2 should create a monomial with word [idx, idx]
        x2 = x[1] * x[1]
        @test x2 isa Term
        @test degree(x2.monomial) == 2

        # x^0 should be identity
        x0 = one(x[1])
        @test isone(x0)
    end

    @testset "Get Basis" begin
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
end
