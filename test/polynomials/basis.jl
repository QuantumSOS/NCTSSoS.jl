using Test, NCTSSoS

@testset "Basis Generation" begin
    @testset "NonCommutative" begin
        reg, (x,y) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])
        d = 4
        basis = get_ncbasis(reg, d)
        @test length(basis) == 129
        @test issorted(basis)
        @test length(basis) == length(unique(basis))
        @test all(m -> degree(m) <= d, basis)
    end

    @testset "Projector" begin
        reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])
        d = 4
        basis = get_ncbasis(reg, d)
        @test length(basis) == 41
        @test issorted(basis)
        @test length(basis) == length(unique(basis))
        @test all(m -> degree(m) <= d, basis)
    end

    @testset "Unipotent" begin
        reg, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 1:2)])
        d = 4
        basis = get_ncbasis(reg, d)
        @test length(basis) == 41
        @test issorted(basis)
        @test length(basis) == length(unique(basis))
        @test all(m -> degree(m) <= d, basis)
    end

    @testset "Pauli" begin
        reg, _ = create_pauli_variables(1:2)
        d = 4
        basis = get_ncbasis(reg, d)
        @test length(basis) == 16
        @test issorted(basis)
        @test length(basis) == length(unique(basis))
        @test all(m -> degree(m) <= d, basis)
    end

    @testset "Fermionic" begin
        reg, _ = create_fermionic_variables(1:2)
        d = 3
        basis = get_ncbasis(reg, d)
        @test length(basis) == 15
        @test issorted(basis)
        @test length(basis) == length(unique(basis))
        @test all(m -> degree(m) <= d, basis)
    end

    @testset "Bosonic" begin
        reg, _ = create_bosonic_variables(1:2)
        d = 3
        basis = get_ncbasis(reg, d)
        @test length(basis) == 35
        @test issorted(basis)
        @test length(basis) == length(unique(basis))
        @test all(m -> degree(m) <= d, basis)
    end

    @testset "Newton Chip" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])

        diagonal_objective = 1.0 * x[1]^2 + 2.0 * x[2] * x[1] * x[1] * x[2] + 3.0 * x[2] * x[1]
        pop = polyopt(diagonal_objective, reg)

        x1x2 = only(monomials(x[1] * x[2]))
        expected_order2 = sort([one(x[1]), x[1], x[2], x1x2])

        basis_order2 = newton_chip_basis(pop, 2)
        @test basis_order2 == expected_order2
        @test issorted(basis_order2)
        @test length(basis_order2) == length(unique(basis_order2))
        @test all(m -> degree(m) <= 2, basis_order2)

        basis_order1 = newton_chip_basis(pop, 1)
        @test basis_order1 == sort([one(x[1]), x[1], x[2]])

        nondiagonal_objective = 1.0 * x[1] * x[2] + 1.0 * x[2] * x[1]
        nondiagonal_pop = polyopt(nondiagonal_objective, reg)
        @test newton_chip_basis(nondiagonal_pop, 2) == [one(x[1])]

        reg_with_unused, (x_unused, y_unused) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])
        single_site_pop = polyopt(1.0 * x_unused[1]^2 + 1.0, reg_with_unused)
        @test newton_chip_basis(single_site_pop, 1) == sort([one(x_unused[1]), x_unused[1]])

        reg_ex34, (vars_ex34,) = create_noncommutative_variables([("X", 1:3)])
        X, Y, Z = vars_ex34
        ex34_objective = X^2 - X * Y - Y * X + 3.0 * Y^2 -
            2.0 * X * Y * X + 2.0 * X * Y^2 * X -
            Y * Z - Z * Y + 6.0 * Z^2 + 9.0 * Y^2 * Z +
            9.0 * Z * Y^2 - 54.0 * Z * Y * Z +
            142.0 * Z * Y^2 * Z
        ex34_pop = polyopt(ex34_objective, reg_ex34)

        yx = only(monomials(Y * X))
        yz = only(monomials(Y * Z))
        expected_ex34_basis = [one(X), X, Y, Z, yx, yz]

        @test newton_chip_basis(ex34_pop, 2) == expected_ex34_basis
    end
end
