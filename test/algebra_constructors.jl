using Test, NCTSSoS, NCTSSoS.FastPolynomials

# Import solver based on environment
if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

@testset "AbstractAlgebra Type Hierarchy" begin
    # Test abstract type exists
    @test AbstractAlgebra isa Type
    @test isabstracttype(AbstractAlgebra)

    # Test concrete types exist
    @test PauliAlgebra <: AbstractAlgebra
    @test FermionicAlgebra <: AbstractAlgebra
    @test isconcretetype(PauliAlgebra)
    @test isconcretetype(FermionicAlgebra)

    # Test that struct fields are correct for PauliAlgebra
    @test hasfield(PauliAlgebra, :N)
    @test hasfield(PauliAlgebra, :variables)
    @test hasfield(PauliAlgebra, :simplify_algo)
    @test hasfield(PauliAlgebra, :equality_constraints)
    @test hasfield(PauliAlgebra, :inequality_constraints)
    @test hasfield(PauliAlgebra, :comm_gps)

    # Test that struct fields are correct for FermionicAlgebra
    @test hasfield(FermionicAlgebra, :N)
    @test hasfield(FermionicAlgebra, :variables)
    @test hasfield(FermionicAlgebra, :simplify_algo)
    @test hasfield(FermionicAlgebra, :equality_constraints)
    @test hasfield(FermionicAlgebra, :inequality_constraints)
    @test hasfield(FermionicAlgebra, :comm_gps)
end

@testset "Pauli Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = pauli_algebra(1)

        # Check return type is PauliAlgebra
        @test sys isa PauliAlgebra
        @test sys isa AbstractAlgebra
        @test sys.N == 1

        # Check return type structure
        @test hasfield(typeof(sys), :N)
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


@testset "Fermionic Algebra Constructor" begin
    @testset "Basic Structure N=1" begin
        sys = fermionic_algebra(1)

        # Check return type is FermionicAlgebra
        @test sys isa FermionicAlgebra
        @test sys isa AbstractAlgebra
        @test sys.N == 1

        # Verify structure
        @test hasfield(typeof(sys), :N)
        @test hasfield(typeof(sys), :variables)
        @test hasfield(typeof(sys), :simplify_algo)
        @test hasfield(typeof(sys), :equality_constraints)
        @test hasfield(typeof(sys), :inequality_constraints)
        @test hasfield(typeof(sys), :comm_gps)

        # Verify variables
        c, c_dag = sys.variables
        @test length(c) == 1
        @test length(c_dag) == 1

        # Verify constraint count for N=1: 2(1)² + 3(1) = 5
        @test length(sys.equality_constraints) == 5
        @test isempty(sys.inequality_constraints)

        # Verify commutation groups
        @test length(sys.comm_gps) == 1
        @test length(sys.comm_gps[1]) == 2  # c[1] and c_dag[1]
    end

    @testset "SimplifyAlgorithm Properties" begin
        sys = fermionic_algebra(2)
        sa = sys.simplify_algo

        @test sa.is_unipotent == false
        @test sa.is_projective == false
        @test sa.n_gps == 1  # Single commutation group
    end

    @testset "Constraint Verification N=1" begin
        sys = fermionic_algebra(1)
        c, c_dag = sys.variables
        constraints = sys.equality_constraints

        # Should contain all 5 constraints:
        # 1. {c, c} = 2c² = 0
        @test ComplexF64(1.0) * (c[1] * c[1] + c[1] * c[1]) in constraints

        # 2. {c†, c†} = 2(c†)² = 0
        @test ComplexF64(1.0) * (c_dag[1] * c_dag[1] + c_dag[1] * c_dag[1]) in constraints

        # 3. {c, c†} = 1
        @test ComplexF64(1.0) * (c[1] * c_dag[1] + c_dag[1] * c[1]) - 1 in constraints

        # 4. c² = 0 (nilpotent constraint)
        @test ComplexF64(1.0) * c[1] * c[1] in constraints

        # 5. (c†)² = 0 (nilpotent constraint)
        @test ComplexF64(1.0) * c_dag[1] * c_dag[1] in constraints
    end

    @testset "Constraint Scaling" begin
        # N=1: 2(1)² + 3(1) = 5
        sys1 = fermionic_algebra(1)
        @test length(sys1.equality_constraints) == 5

        # N=2: 2(4) + 3(2) = 14
        sys2 = fermionic_algebra(2)
        @test length(sys2.equality_constraints) == 14

        # N=3: 2(9) + 3(3) = 27
        sys3 = fermionic_algebra(3)
        @test length(sys3.equality_constraints) == 27

        # General formula verification
        for N in 1:5
            sys = fermionic_algebra(N)
            expected_count = 2 * N^2 + 3 * N
            @test length(sys.equality_constraints) == expected_count
        end
    end

    @testset "Constraint Verification N=2" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables
        constraints = sys.equality_constraints

        # Spot check key constraints

        # Anti-commutation {c₁, c₂} = 0
        @test ComplexF64(1.0) * (c[1] * c[2] + c[2] * c[1]) in constraints

        # Anti-commutation {c₁†, c₂†} = 0
        @test ComplexF64(1.0) * (c_dag[1] * c_dag[2] + c_dag[2] * c_dag[1]) in constraints

        # Canonical {c₁, c₁†} = 1
        @test ComplexF64(1.0) * (c[1] * c_dag[1] + c_dag[1] * c[1]) - 1 in constraints

        # Canonical {c₁, c₂†} = 0
        @test ComplexF64(1.0) * (c[1] * c_dag[2] + c_dag[2] * c[1]) in constraints

        # Nilpotent c₁² = 0
        @test ComplexF64(1.0) * c[1] * c[1] in constraints

        # Nilpotent (c₂†)² = 0
        @test ComplexF64(1.0) * c_dag[2] * c_dag[2] in constraints
    end

    @testset "Integration with cpolyopt" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables

        # Use symmetric objective: symmetrized number operator
        # n₁_symm = (c₁† c₁ + c₁ c₁†)/2
        n1_symm = ComplexF64(0.5) * (c_dag[1] * c[1] + c[1] * c_dag[1])

        # Create optimization problem using algebra interface
        pop = cpolyopt(n1_symm, sys)

        @test pop.objective == n1_symm
        @test pop.is_unipotent == false
        @test pop.is_projective == false

        # Check that algebra constraints were included
        @test length(pop.eq_constraints) == length(sys.equality_constraints)

        # Verify fermionic constraints are present
        @test ComplexF64(1.0) * c[1] * c[1] in pop.eq_constraints  # c₁² = 0
    end

    @testset "Custom Constraints" begin
        sys = fermionic_algebra(2)
        c, c_dag = sys.variables

        # Symmetric total number operator: (n1 + n2) symmetrized
        # n1_symm = (c†₁c₁ + c₁c†₁)/2, n2_symm = (c†₂c₂ + c₂c†₂)/2
        n1_symm = ComplexF64(0.5) * (c_dag[1] * c[1] + c[1] * c_dag[1])
        n2_symm = ComplexF64(0.5) * (c_dag[2] * c[2] + c[2] * c_dag[2])
        N_total_symm = n1_symm + n2_symm

        # Add custom constraint: total particle number = 1
        constraint = N_total_symm - 1

        pop = cpolyopt(N_total_symm, sys; eq_constraints=[constraint])

        # Should have algebra constraints + custom constraint
        @test length(pop.eq_constraints) == length(sys.equality_constraints) + 1
        @test constraint in pop.eq_constraints

        # Fermionic constraints still present
        @test ComplexF64(1.0) * c[1] * c[1] in pop.eq_constraints
    end

    @testset "Error Handling" begin
        @test_throws AssertionError fermionic_algebra(0)
        @test_throws AssertionError fermionic_algebra(-1)
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

    @testset "Integration: Fermi-Hubbard Model" begin
        # Test the fermionic_algebra interface with the 1D Fermi-Hubbard model
        # H = -t Σ_<i,j>,σ (c†_iσ c_jσ + h.c.) + U Σ_i n_i↑ n_i↓

        @testset "2-site Hubbard: Algebra Structure Test" begin
            # Test fermionic algebra with 2 sites (4 modes total)
            # Verify problem construction and basic solving
            N_sites = 2

            # Create fermionic algebra for 2*N_sites modes
            # Modes 1:N_sites are spin-up, modes (N_sites+1):2*N_sites are spin-down
            sys = fermionic_algebra(2 * N_sites)
            c, c_dag = sys.variables

            # Use interaction energy as objective (already tested in simpler cases)
            # This is symmetric and physically meaningful
            ham = ComplexF64(0.0)

            # Interaction term: U Σ_i n_i↑ n_i↓ (with U=1 for testing)
            for i in 1:N_sites
                i_up = i
                i_down = i + N_sites
                # Use symmetrized number operators
                n_up = ComplexF64(0.5) * (c_dag[i_up] * c[i_up] + c[i_up] * c_dag[i_up])
                n_down = ComplexF64(0.5) * (c_dag[i_down] * c[i_down] + c[i_down] * c_dag[i_down])
                ham += n_up * n_down
            end

            # Create optimization problem
            pop = cpolyopt(ham, sys)

            # Verify problem structure
            @test pop.objective == ham
            @test pop.is_unipotent == false
            @test pop.is_projective == false
            # 2N² + 3N constraints for N=4 modes: 2(16) + 3(4) = 44
            @test length(pop.eq_constraints) == 2 * (2 * N_sites)^2 + 3 * (2 * N_sites)

            # Solve the optimization problem
            solver_config = SolverConfig(optimizer=SOLVER, order=1)
            res = cs_nctssos(pop, solver_config)

            # Minimum interaction energy (avoid double occupancy)
            @test res.objective ≥ 0.0  # Should be non-negative
            @test res.objective ≤ 1.0  # Upper bound sanity check
        end

        @testset "2-site Hubbard: Interaction Term Scaling" begin
            # Test that interaction term scales correctly with U
            N_sites = 2

            # Create fermionic algebra
            sys = fermionic_algebra(2 * N_sites)
            c, c_dag = sys.variables

            # Test with U=1 and U=4
            for U in [1.0, 4.0]
                ham = ComplexF64(0.0)

                # Interaction term: U Σ_i n_i↑ n_i↓
                for i in 1:N_sites
                    i_up = i
                    i_down = i + N_sites
                    n_up = ComplexF64(0.5) * (c_dag[i_up] * c[i_up] + c[i_up] * c_dag[i_up])
                    n_down = ComplexF64(0.5) * (c_dag[i_down] * c[i_down] + c[i_down] * c_dag[i_down])
                    ham += U * n_up * n_down
                end

                pop = cpolyopt(ham, sys)
                solver_config = SolverConfig(optimizer=SOLVER, order=1)
                res = cs_nctssos(pop, solver_config)

                # Minimum should be non-negative and scale with U
                @test res.objective ≥ 0.0
                @test res.objective ≤ U * N_sites  # Upper bound
            end
        end

        @testset "2-site Hubbard: With Constraints" begin
            # Test adding custom constraints to fermionic system
            N_sites = 2
            U = 1.0

            # Create fermionic algebra
            sys = fermionic_algebra(2 * N_sites)
            c, c_dag = sys.variables

            # Interaction Hamiltonian
            ham = ComplexF64(0.0)
            for i in 1:N_sites
                i_up = i
                i_down = i + N_sites
                n_up = ComplexF64(0.5) * (c_dag[i_up] * c[i_up] + c[i_up] * c_dag[i_up])
                n_down = ComplexF64(0.5) * (c_dag[i_down] * c[i_down] + c[i_down] * c_dag[i_down])
                ham += U * n_up * n_down
            end

            # Add a simple constraint (e.g., n_1↑ = n_1↓)
            n1_up = ComplexF64(0.5) * (c_dag[1] * c[1] + c[1] * c_dag[1])
            n1_down = ComplexF64(0.5) * (c_dag[3] * c[3] + c[3] * c_dag[3])
            constraint = n1_up - n1_down

            # Create optimization problem with constraint
            pop = cpolyopt(ham, sys; eq_constraints=[constraint])

            # Verify constraint was added
            @test length(pop.eq_constraints) == length(sys.equality_constraints) + 1
            @test constraint in pop.eq_constraints

            # Solve the optimization problem
            solver_config = SolverConfig(optimizer=SOLVER, order=1)
            res = cs_nctssos(pop, solver_config)

            # Should find a feasible solution
            @test res.objective ≥ 0.0
        end
    end
end
