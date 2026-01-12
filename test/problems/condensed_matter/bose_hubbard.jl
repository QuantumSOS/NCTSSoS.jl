# =============================================================================
# Bose-Hubbard Model Tests
# =============================================================================
# Tests for the Bose-Hubbard model on chain lattices.
#
# The Bose-Hubbard Hamiltonian is:
#   H = -t Σ_{<i,j>} (b†_i b_j + b†_j b_i) + (U/2) Σ_i n̂_i(n̂_i - 1) - μ Σ_i n̂_i
#
# where:
#   <i,j> denotes neighboring sites on the lattice
#   b†_i, b_i are bosonic creation and annihilation operators
#   n̂_i = b†_i b_i is the number operator on site i
#   t = hopping amplitude
#   U = on-site interaction strength
#   μ = chemical potential
#
# The bosonic operators satisfy the canonical commutation relations (CCR):
#   [b_i, b†_j] = δ_{ij}
#   [b_i, b_j] = 0
#   [b†_i, b†_j] = 0
# =============================================================================

using NCTSSoS, Test
using JuMP
using NCTSSoS: simplify, degree  # Disambiguate from JuMP/Graphs

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

@testset "Bose-Hubbard Hamiltonian Construction" begin
    N = 4  # Number of sites

    # Create bosonic variables (b = annihilation, b_dag = creation)
    registry, (b, b_dag) = create_bosonic_variables(1:N)

    @testset "Variable creation" begin
        @test length(b) == N
        @test length(b_dag) == N
        @test b[1] isa NormalMonomial{BosonicAlgebra}
        @test b_dag[1] isa NormalMonomial{BosonicAlgebra}
    end

    @testset "Number operator n_i = b†_i b_i" begin
        # Number operator for site 1
        n1 = b_dag[1] * b[1]
        n1_poly = Polynomial(n1)
        @test n1_poly isa Polynomial{BosonicAlgebra}
        @test degree(n1_poly) == 2

        # Already normal-ordered; should be single-term with unit coefficient
        @test length(n1_poly.terms) == 1
        @test n1_poly.terms[1].coefficient == 1.0
    end

    @testset "Hopping terms b†_i b_j + b†_j b_i" begin
        # Hopping between sites 1 and 2
        hop_12 = b_dag[1] * b[2] + b_dag[2] * b[1]

        @test hop_12 isa Polynomial{BosonicAlgebra}
        @test length(hop_12.terms) == 2
    end

    @testset "Interaction term n_i(n_i - 1) = b†_i b_i b†_i b_i - b†_i b_i" begin
        # For site 1: n₁(n₁ - 1) = (b†₁ b₁)(b†₁ b₁ - 1) = b†₁ b₁ b†₁ b₁ - b†₁ b₁
        n1 = b_dag[1] * b[1]
        n1_squared = n1 * n1  # b†₁ b₁ b†₁ b₁

        # After simplification using CCR: b b† = b† b + 1
        # b†₁ b₁ b†₁ b₁ = b†₁ (b†₁ b₁ + 1) b₁ = b†₁ b†₁ b₁ b₁ + b†₁ b₁
        simplified = Polynomial(n1_squared)

        # Should have 2 terms: (b†)² b² and b† b
        @test length(simplified.terms) == 2

        # Check for the b†₁ b†₁ b₁ b₁ term (degree 4)
        degree4_terms = filter(t -> degree(t.monomial) == 4, simplified.terms)
        @test length(degree4_terms) == 1
        @test degree4_terms[1].coefficient == 1.0

        # Check for the b†₁ b₁ term (degree 2)
        degree2_terms = filter(t -> degree(t.monomial) == 2, simplified.terms)
        @test length(degree2_terms) == 1
        @test degree2_terms[1].coefficient == 1.0
    end

    @testset "Full Hamiltonian construction (OBC chain)" begin
        t_hop = 1.0  # Hopping amplitude
        U = 2.0      # On-site interaction
        μ = 0.5      # Chemical potential

        # Hopping term: -t Σ_{<i,j>} (b†_i b_j + b†_j b_i)
        ham_hop = sum(-t_hop * (b_dag[k] * b[k+1] + b_dag[k+1] * b[k]) for k in 1:N-1)

        # Number operators for each site
        n = [b_dag[k] * b[k] for k in 1:N]

        # Interaction term: (U/2) Σ_i n_i(n_i - 1)
        # n_i(n_i - 1) = n_i² - n_i
        ham_int = sum((U / 2) * (n[k] * n[k] - n[k]) for k in 1:N)

        # Chemical potential term: -μ Σ_i n_i
        ham_chem = sum(-μ * n[k] for k in 1:N)

        # Total Hamiltonian
        ham = ham_hop + ham_int + ham_chem

        @test ham isa Polynomial{BosonicAlgebra}

        # The Hamiltonian should have multiple terms
        @test length(ham.terms) > 0
    end
end

@testset "Bose-Hubbard CCR Verification" begin
    # Verify canonical commutation relations: [b_i, b†_j] = δ_{ij}
    # This is the key algebraic property that the bosonic simplification must satisfy

    registry, (b, b_dag) = create_bosonic_variables(1:3)

    @testset "Same-site commutator: b b† = b† b + 1" begin
        # [b₁, b₁†] = b₁ b₁† - b₁† b₁ = 1
        # So: b₁ b₁† = b₁† b₁ + 1
        lhs = simplify(b[1] * b_dag[1])  # b b†
        rhs = simplify(b_dag[1] * b[1])  # b† b

        # lhs should have 2 terms: b† b + 1
        @test length(lhs.terms) == 2

        # Find the identity term (degree 0) with coefficient 1
        identity_term = filter(t -> degree(t.monomial) == 0, lhs.terms)
        @test length(identity_term) == 1
        @test identity_term[1].coefficient == 1.0

        # Find the b† b term (degree 2) with coefficient 1
        number_term = filter(t -> degree(t.monomial) == 2, lhs.terms)
        @test length(number_term) == 1
        @test number_term[1].coefficient == 1.0

        # rhs should have 1 term: b† b
        @test length(rhs.terms) == 1
        @test rhs.terms[1].coefficient == 1.0
        @test degree(rhs.terms[1].monomial) == 2
    end

    @testset "Different-site commutator: [b_i, b†_j] = 0 for i ≠ j" begin
        # b₁ b₂† = b₂† b₁ (they commute, no delta term)
        lhs = simplify(b[1] * b_dag[2])
        rhs = simplify(b_dag[2] * b[1])

        # Both should be single terms with coefficient 1
        @test length(lhs.terms) == 1
        @test length(rhs.terms) == 1
        @test lhs.terms[1].coefficient == 1.0
        @test rhs.terms[1].coefficient == 1.0

        # Both should have the same normal-ordered form: b₂† b₁
        @test lhs.terms[1].monomial == rhs.terms[1].monomial
    end

    @testset "Annihilators commute: [b_i, b_j] = 0" begin
        # b₁ b₂ = b₂ b₁
        m12 = simplify(b[1] * b[2])
        m21 = simplify(b[2] * b[1])

        @test length(m12.terms) == 1
        @test length(m21.terms) == 1
        @test m12.terms[1].monomial == m21.terms[1].monomial
    end

    @testset "Creators commute: [b†_i, b†_j] = 0" begin
        # b₁† b₂† = b₂† b₁†
        m12 = simplify(b_dag[1] * b_dag[2])
        m21 = simplify(b_dag[2] * b_dag[1])

        @test length(m12.terms) == 1
        @test length(m21.terms) == 1
        @test m12.terms[1].monomial == m21.terms[1].monomial
    end
end

@testset "Bose-Hubbard Ground State (vacuum reference)" begin
    # With large positive chemical potential μ > 2t, the vacuum is the ground state
    # because adding any particle costs more energy than the hopping can provide.
    # In this limit, E₀ = 0 (vacuum energy).

    N = 2
    t_hop = 1.0
    μ = 10.0  # Large chemical potential → vacuum is ground state

    registry, (b, b_dag) = create_bosonic_variables(1:N)

    # Number operators
    n = [b_dag[k] * b[k] for k in 1:N]

    # Hamiltonian: hopping + chemical potential (no interaction)
    # H = -t(b₁†b₂ + b₂†b₁) + μ(n₁ + n₂)
    ham = -t_hop * (b_dag[1] * b[2] + b_dag[2] * b[1]) + μ * (n[1] + n[2])

    # Create optimization problem
    pop = polyopt(ham, registry)

    # Solve with order-1 relaxation
    solver_config = SolverConfig(optimizer=SOLVER, order=1)
    res = cs_nctssos(pop, solver_config)

    # With μ >> t, the vacuum (no particles) is the ground state with E = 0
    # The SDP should find a lower bound close to 0
    @test res.objective isa Real
    @test res.objective ≥ -1e-4  # Should be ≈ 0 (vacuum energy)
end

@testset "Bose-Hubbard Chain Ground State" begin
    # Bose-Hubbard Model on a Chain Lattice (OBC)
    #
    # H = -J Σ_{<ij>} (a†_i a_j + a†_j a_i) + (U/2) Σ_i n_i(n_i - 1)
    #
    # where:
    #   <i,j> denotes neighboring sites on the chain (OBC)
    #   a†_i, a_i are bosonic creation and annihilation operators
    #   n_i = a†_i a_i is the number operator on site i
    #   J = hopping amplitude
    #   U = on-site interaction strength
    #
    # BROKEN: The SDP relaxation is unbounded without particle number constraints.
    # Bosonic systems (unlike fermionic) have no natural occupation limit, so the
    # optimization diverges to -∞. Need to add either:
    #   1. Particle number constraint (Σ n_i = N), or
    #   2. Occupation cutoff (n_i ≤ n_max), or
    #   3. Chemical potential to bound the problem

    M = 4  # Number of sites
    J = 1.0
    U = 1.0

    # Create bosonic variables
    registry, (a, a_dag) = create_bosonic_variables(1:M)

    # Number operators
    n = [a_dag[k] * a[k] for k in 1:M]

    # Hopping term: -J Σ_{<ij>} (a†_i a_j + a†_j a_i)
    ham_hop = sum(-J * (a_dag[k] * a[k+1] + a_dag[k+1] * a[k]) for k in 1:M-1)

    # Interaction term: (U/2) Σ_i n_i(n_i - 1)
    ham_int = sum((U / 2) * (n[k] * n[k] - n[k]) for k in 1:M)

    # Total Hamiltonian
    ham = ham_hop + ham_int

    # Create optimization problem
    pop = polyopt(ham, registry)

    # Solve with order-2 relaxation
    solver_config = SolverConfig(optimizer=SOLVER, order=2)
    res = cs_nctssos(pop, solver_config)

    # Refine with higher-order iterations
    res = cs_nctssos_higher(pop, res, solver_config)

    # Expected ground state energy
    E0_exact = -6.911424152031056

    println("Ground state energy lower bound: ", res.objective)
    println("Exact ground state energy: ", E0_exact)

    @test_broken res.objective ≈ E0_exact atol = 1e-4
end
