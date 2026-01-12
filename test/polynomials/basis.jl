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
end
