using Test, NCTSSoS, NCTSSoS.FastPolynomials

# Import solver based on environment
if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

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


if haskey(ENV, "LOCAL_TESTING")
    @testset "Integration: XXX Model with Pauli Algebra" begin

        # Test that the new pauli_algebra interface produces the same
        # numerical results as the manual setup for the XXX Heisenberg Model
        T = ComplexF64
        N = 6
        J1 = 1.0

        # Create algebra system using new interface
        sys = pauli_algebra(N)
        x, y, z = sys.variables

        # Construct XXX Heisenberg Hamiltonian
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

        # Create optimization problem using new algebra interface
        pop = cpolyopt(ham, sys)

        # Verify problem structure
        @test pop.objective == ham
        @test pop.is_unipotent == true
        @test pop.is_projective == false
        @test length(pop.eq_constraints) == 6 * N  # 6 constraints per site

        # Solve the optimization problem
        solver_config = SolverConfig(optimizer=SOLVER, order=2)
        res = cs_nctssos(pop, solver_config)

        # Verify against known result from heisenberg.jl test
        @test res.objective / N ≈ -0.467129 atol = 1e-6
    end

    @testset "Integration: J1-J2 Model with Pauli Algebra" begin
        # Test the new pauli_algebra interface with J1-J2 Heisenberg Model
        # This requires higher order relaxation (cs_nctssos_higher)
        T = ComplexF64
        N = 6
        J1 = 1.0
        J2 = 0.2

        # Create algebra system using new interface
        sys = pauli_algebra(N)
        x, y, z = sys.variables

        # Construct J1-J2 Heisenberg Hamiltonian
        # J1: nearest-neighbor interactions, J2: next-nearest-neighbor
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] +
                  T(J2 / 4) * op[i] * op[mod1(i + 2, N)]
                  for op in [x, y, z] for i in 1:N)

        # Create optimization problem using new algebra interface
        pop = cpolyopt(ham, sys)

        # Verify problem structure
        @test pop.objective == ham
        @test pop.is_unipotent == true
        @test pop.is_projective == false
        @test length(pop.eq_constraints) == 6 * N

        # Solve with MMD term sparsity algorithm
        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)
        res = cs_nctssos_higher(pop, res, solver_config)

        # Verify against known result from heisenberg.jl test
        @test res.objective / N ≈ -0.4270083225302217 atol = 1e-6
    end

    @testset "Integration: 1D Transverse Field Ising Model" begin
        # Test the new pauli_algebra interface with Transverse Field Ising Model
        # This model combines ZZ interactions with transverse field X
        N = 3
        J = 1.0
        h = 2.0

        # Create algebra system using new interface
        sys = pauli_algebra(N)
        x, y, z = sys.variables

        # Test both periodic and open boundary conditions
        for (periodic, true_ans) in zip((true, false), (-1.0175918, -1.0104160))
            # Construct Transverse Field Ising Hamiltonian
            # -J/4 * sum(Z_i * Z_{i+1}) - h/2 * sum(X_i)
            ham = sum(-complex(J / 4) * z[i] * z[mod1(i + 1, N)]
                      for i in 1:(periodic ? N : N - 1)) +
                  sum(-h / 2 * x[i] for i in 1:N)

            # Create optimization problem using new algebra interface
            pop = cpolyopt(ham, sys)

            # Verify problem structure
            @test pop.objective == ham
            @test pop.is_unipotent == true
            @test pop.is_projective == false
            @test length(pop.eq_constraints) == 6 * N

            # Solve the optimization problem
            solver_config = SolverConfig(optimizer=SOLVER, order=2)
            res = cs_nctssos(pop, solver_config)

            # Verify against known result from interface.jl test
            @test res.objective / N ≈ true_ans atol = 1e-6
        end
    end
end
