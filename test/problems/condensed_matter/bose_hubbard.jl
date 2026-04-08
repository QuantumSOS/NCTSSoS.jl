# Bose-Hubbard Model Tests
# Tests for the Bose-Hubbard model on chain lattices.
#
# The Bose-Hubbard Hamiltonian is:
#   H = -t Σ_{<i,j>} (b†_i b_j + b†_j b_i) + (U/2) Σ_i n̂_i(n̂_i - 1) + μ Σ_i n̂_i
#
# where:
#   <i,j> denotes neighboring sites on the lattice
#   b†_i, b_i are bosonic creation and annihilation operators
#   n̂_i = b†_i b_i is the number operator on site i
#   t = hopping amplitude
#   U = on-site interaction strength
#   μ = occupation penalty (equivalently, minus the usual chemical potential)
#
# The bosonic operators satisfy the canonical commutation relations (CCR):
#   [b_i, b†_j] = δ_{ij}
#   [b_i, b_j] = 0
#   [b†_i, b†_j] = 0

using NCTSSoS, Test
using JuMP, LinearAlgebra
using NCTSSoS: simplify, degree  # Disambiguate from JuMP/Graphs

const BOSE_HUBBARD_EXPECTATIONS_PATH = "expectations/bose_hubbard.toml"

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

function _two_site_bose_hubbard_sector_matrix(t_hop::Real, U::Real)
    # Fixed-N=2 basis ordered as |2,0⟩, |1,1⟩, |0,2⟩.
    # The hopping term contributes -√2 t between adjacent basis states and the
    # on-site interaction contributes U on the doubly occupied states.
    hop = -sqrt(2.0) * float(t_hop)
    return Hermitian(
        Float64[
            U hop 0.0
            hop 0.0 hop
            0.0 hop U
        ]
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
        @test first(n1_poly.terms[1]) == 1.0
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
        degree4_terms = filter(t -> degree(last(t)) == 4, simplified.terms)
        @test length(degree4_terms) == 1
        @test first(degree4_terms[1]) == 1.0

        # Check for the b†₁ b₁ term (degree 2)
        degree2_terms = filter(t -> degree(last(t)) == 2, simplified.terms)
        @test length(degree2_terms) == 1
        @test first(degree2_terms[1]) == 1.0
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
        identity_term = filter(t -> degree(last(t)) == 0, lhs.terms)
        @test length(identity_term) == 1
        @test first(identity_term[1]) == 1.0

        # Find the b† b term (degree 2) with coefficient 1
        number_term = filter(t -> degree(last(t)) == 2, lhs.terms)
        @test length(number_term) == 1
        @test first(number_term[1]) == 1.0

        # rhs should have 1 term: b† b
        @test length(rhs.terms) == 1
        @test first(rhs.terms[1]) == 1.0
        @test degree(last(rhs.terms[1])) == 2
    end

    @testset "Different-site commutator: [b_i, b†_j] = 0 for i ≠ j" begin
        # b₁ b₂† = b₂† b₁ (they commute, no delta term)
        lhs = simplify(b[1] * b_dag[2])
        rhs = simplify(b_dag[2] * b[1])

        # Both should be single terms with coefficient 1
        @test length(lhs.terms) == 1
        @test length(rhs.terms) == 1
        @test first(lhs.terms[1]) == 1.0
        @test first(rhs.terms[1]) == 1.0

        # Both should have the same normal-ordered form: b₂† b₁
        @test last(lhs.terms[1]) == last(rhs.terms[1])
    end

    @testset "Annihilators commute: [b_i, b_j] = 0" begin
        # b₁ b₂ = b₂ b₁
        m12 = simplify(b[1] * b[2])
        m21 = simplify(b[2] * b[1])

        @test length(m12.terms) == 1
        @test length(m21.terms) == 1
        @test last(m12.terms[1]) == last(m21.terms[1])
    end

    @testset "Creators commute: [b†_i, b†_j] = 0" begin
        # b₁† b₂† = b₂† b₁†
        m12 = simplify(b_dag[1] * b_dag[2])
        m21 = simplify(b_dag[2] * b_dag[1])

        @test length(m12.terms) == 1
        @test length(m21.terms) == 1
        @test last(m12.terms[1]) == last(m21.terms[1])
    end
end

@testset "Bose-Hubbard Ground State (vacuum reference)" begin
    # With a large positive occupation penalty μ > 2t, the vacuum is the ground state
    # because adding any particle costs more energy than the hopping can provide.
    # In this limit, E₀ = 0 (vacuum energy).

    N = 2
    t_hop = 1.0
    μ = 10.0  # Large occupation penalty → vacuum is ground state

    registry, (b, b_dag) = create_bosonic_variables(1:N)

    # Number operators
    n = [b_dag[k] * b[k] for k in 1:N]

    # Hamiltonian: hopping + occupation penalty (no interaction)
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

@testset "Two-site Bose-Hubbard model (canonical N = 2 sector)" begin
    # This is the smallest nontrivial bosonic condensed-matter regression with an
    # exact fixed-particle-number reference. We impose total particle number
    # N = 2 with a one-sided moment equality constraint so the SDP solves the
    # same canonical-sector problem as the 3×3 diagonalization below.

    N = 2
    t_hop = 1.0
    U = 2.0
    # Oracle provenance:
    # - objective: exact 3×3 fixed-N sector diagonalization below
    # - sides/nuniq: reviewed local implementation oracles recorded in
    #   `test/data/expectations/bose_hubbard.toml`
    oracle = expectations_oracle(BOSE_HUBBARD_EXPECTATIONS_PATH, "two_site_n2_u2_order2")

    registry, (b, b_dag) = create_bosonic_variables(1:N)
    n = [b_dag[k] * b[k] for k in 1:N]

    ham_hop = -t_hop * (b_dag[1] * b[2] + b_dag[2] * b[1])
    ham_int = (U / 2) * sum(n[k] * n[k] - n[k] for k in 1:N)
    ham = ham_hop + ham_int

    total_number = 1.0 * sum(n)
    canonical_sector = total_number - 2.0 * one(ham)

    exact_hamiltonian = _two_site_bose_hubbard_sector_matrix(t_hop, U)
    exact_e0 = eigmin(exact_hamiltonian)

    @test exact_e0 ≈ 1.0 - sqrt(5.0) atol = 1e-12

    pop = polyopt(ham, registry; moment_eq_constraints=[canonical_sector])
    result = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=2))

    @test result.objective ≈ oracle.opt atol = 1e-6
    @test result.objective ≈ exact_e0 atol = 1e-6
    @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
    @test result.n_unique_moment_matrix_elements == oracle.nuniq
end
