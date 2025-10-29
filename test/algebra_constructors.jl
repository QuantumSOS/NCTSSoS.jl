using Test, NCTSSoS, NCTSSoS.FastPolynomials

@testset "Pauli Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = pauli_algebra(1)

        # Check return type structure
        @test hasfield(typeof(sys), :variables)
        @test hasfield(typeof(sys), :simplify_algo)
        @test hasfield(typeof(sys), :equality_constraints)
        @test hasfield(typeof(sys), :inequality_constraints)
        @test hasfield(typeof(sys), :comm_gps)

        # Check variables
        x, y, z = sys.variables
        @test length(x) == 1
        @test length(y) == 1
        @test length(z) == 1

        # Check constraints (6 per site)
        @test length(sys.equality_constraints) == 6
        @test isempty(sys.inequality_constraints)

        # Check commutation groups
        @test length(sys.comm_gps) == 1
        @test sort(sys.comm_gps[1]) == sort([x[1], y[1], z[1]])
    end

    @testset "Basic Structure N=3" begin
        sys = pauli_algebra(3)
        x, y, z = sys.variables

        @test length(x) == 3
        @test length(y) == 3
        @test length(z) == 3

        # 6 constraints per site × 3 sites = 18 constraints
        @test length(sys.equality_constraints) == 18

        # 3 commutation groups, one per site
        @test length(sys.comm_gps) == 3
        for i in 1:3
            @test sort(sys.comm_gps[i]) == sort([x[i], y[i], z[i]])
        end
    end

    @testset "SimplifyAlgorithm Properties" begin
        sys = pauli_algebra(2)
        sa = sys.simplify_algo

        @test sa.is_unipotent == true
        @test sa.is_projective == false
        @test sa.n_gps == 2
    end

    @testset "Commutation Relations Encoded Correctly" begin
        sys = pauli_algebra(1)
        x, y, z = sys.variables

        # The 6 constraints should be:
        # x*y - im*z, y*x + im*z, y*z - im*x, z*y + im*x, z*x - im*y, x*z + im*y
        constraints = sys.equality_constraints

        @test x[1] * y[1] - im * z[1] in constraints
        @test y[1] * x[1] + im * z[1] in constraints
        @test y[1] * z[1] - im * x[1] in constraints
        @test z[1] * y[1] + im * x[1] in constraints
        @test z[1] * x[1] - im * y[1] in constraints
        @test x[1] * z[1] + im * y[1] in constraints
    end

    @testset "Integration with cpolyopt" begin
        sys = pauli_algebra(2)
        x, y, z = sys.variables

        # Create a simple Hamiltonian
        ham = ComplexF64(0.5) * (x[1] * x[2] + y[1] * y[2] + z[1] * z[2])

        # Should be able to create optimization problem
        pop = cpolyopt(ham;
            eq_constraints=sys.equality_constraints,
            comm_gps=sys.comm_gps,
            is_unipotent=true)

        @test pop.objective == ham
        @test pop.eq_constraints == sys.equality_constraints
        @test pop.is_unipotent == true
        @test pop.is_projective == false
    end

    @testset "Invalid Input" begin
        @test_throws AssertionError pauli_algebra(0)
        @test_throws AssertionError pauli_algebra(-1)
    end
end


@testset "cpolyopt with Algebra Interface" begin
    @testset "Pauli Algebra Interface" begin
        sys = pauli_algebra(2)
        x, y, z = sys.variables

        ham = ComplexF64(0.5) * (x[1] * x[2] + y[1] * y[2])
        pop = cpolyopt(ham, sys)

        @test pop.objective == ham
        @test pop.is_unipotent == true
        @test pop.is_projective == false
        @test length(pop.eq_constraints) == 12
        @test isempty(pop.ineq_constraints)
        @test pop.comm_gps == sys.comm_gps
    end

    @testset "Pauli Algebra with Additional Constraints" begin
        sys = pauli_algebra(2)
        x, y, z = sys.variables
        ham = ComplexF64(0.5) * (x[1] * x[2])

        custom_eq = [x[1] + y[1]]
        pop = cpolyopt(ham, sys; eq_constraints=custom_eq)

        @test length(pop.eq_constraints) == 13
        @test ComplexF64(1.0) * (x[1] + y[1]) in pop.eq_constraints
    end

    @testset "Pauli Algebra with Inequality Constraints" begin
        sys = pauli_algebra(1)
        x, y, z = sys.variables
        ham = ComplexF64(1.0) * x[1]

        ineq = [ComplexF64(1.0) + z[1]]
        pop = cpolyopt(ham, sys; ineq_constraints=ineq)

        @test length(pop.eq_constraints) == 6
        @test length(pop.ineq_constraints) == 1
        @test ComplexF64(1.0) + z[1] in pop.ineq_constraints
    end

    @testset "Interface Equivalence Test" begin
        sys = pauli_algebra(2)
        x, y, z = sys.variables
        ham = ComplexF64(0.5) * (x[1] * y[2])

        pop1 = cpolyopt(ham, sys)
        pop2 = cpolyopt(ham; eq_constraints=sys.equality_constraints,
            comm_gps=sys.comm_gps, is_unipotent=true, is_projective=false)

        @test pop1.objective == pop2.objective
        @test Set(pop1.eq_constraints) == Set(pop2.eq_constraints)
        @test pop1.ineq_constraints == pop2.ineq_constraints
        @test pop1.comm_gps == pop2.comm_gps
        @test pop1.is_unipotent == pop2.is_unipotent
        @test pop1.is_projective == pop2.is_projective
    end
end
