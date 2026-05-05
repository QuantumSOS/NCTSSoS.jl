using Test, NCTSSoS

@testset "Particle-number constraint helper" begin
    @testset "Bosonic total particle number uses every annihilation mode" begin
        registry, (b, b_dag) = create_bosonic_variables(1:3)

        cons = particle_number_constraint(registry, 2)
        expected = sum(b_dag[i] * b[i] for i in eachindex(b)) - 2.0 * one(cons[1])

        @test cons == [expected]
        @test eltype(cons) == Polynomial{BosonicAlgebra,Int8,Float64}
        @test degree(cons[1]) == 2
    end

    @testset "Fermionic grouped particle numbers match manual number operators" begin
        registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
            ("c_up", 1:2),
            ("c_dn", 1:2),
        ])

        cons = particle_number_constraint(registry, c_up => 1, [c_dn[2]] => 1)
        expected_up = sum(c_up_dag[i] * c_up[i] for i in eachindex(c_up)) - one(cons[1])
        expected_dn2 = c_dn_dag[2] * c_dn[2] - one(cons[2])

        @test cons == [expected_up, expected_dn2]
        @test all(c -> c isa Polynomial{FermionicAlgebra,Int8,Float64}, cons)
    end

    @testset "polyopt promotes real helper constraints for complex objectives" begin
        registry, (a, a_dag) = create_fermionic_variables(1:2)
        ham = ComplexF64(1.0) * sum(a_dag[i] * a[i] for i in eachindex(a))
        helper_constraints = particle_number_constraint(registry, 1)

        pop = polyopt(ham, registry; moment_eq_constraints=helper_constraints)

        @test pop.moment_eq_constraints == helper_constraints
        @test eltype(pop.moment_eq_constraints) == typeof(ham)

        real_ham = sum(a_dag[i] * a[i] for i in eachindex(a))
        bad_constraint = (1.0 + 1.0im) * helper_constraints[1]
        huge_constraint = (big(10.0)^400) * helper_constraints[1]
        @test_throws ArgumentError polyopt(real_ham, registry; moment_eq_constraints=[bad_constraint])
        @test_throws ArgumentError polyopt(real_ham, registry; moment_eq_constraints=[huge_constraint])
    end

    @testset "Validation rejects ambiguous or impossible constraints" begin
        registry, ((c_up, _), (_, c_dn_dag)) = create_fermionic_variables([
            ("c_up", 1:2),
            ("c_dn", 1:2),
        ])
        pauli_registry, _ = create_pauli_variables(1:1)

        bosonic_registry, _ = create_bosonic_variables(1:1)

        @test_throws DomainError particle_number_constraint(registry, c_up => -1)
        @test_throws DomainError particle_number_constraint(registry, c_up => 3)
        @test_throws DomainError particle_number_constraint(bosonic_registry, big(2)^54)
        @test_throws ArgumentError particle_number_constraint(registry, c_dn_dag => 1)
        @test_throws ArgumentError particle_number_constraint(registry, [c_up[1], c_up[1]] => 1)
        @test_throws ArgumentError particle_number_constraint(registry, c_up[1] => 1)
        @test_throws ArgumentError particle_number_constraint(registry, c_up => 1.0)
        @test_throws ArgumentError particle_number_constraint(pauli_registry, 1)
    end
end
