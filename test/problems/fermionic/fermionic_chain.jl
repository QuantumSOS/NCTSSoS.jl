using NCTSSoS, Test

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "standalone_setup.jl"))

# Fermionic Chain Hamiltonian Test
#
# H = j_c Σ_{k=1}^{N-1} (a_k† a_{k+1} + a_{k+1}† a_k)      # hopping (bulk)
#   - j_c (a_N† a_1 + a_1† a_N)                            # hopping (anti-periodic BC)
#   + j_c Σ_{k=1}^{N-1} (a_k† a_{k+1}† + a_{k+1} a_k)      # pairing (bulk)
#   - j_c (a_N† a_1† + a_1 a_N)                            # pairing (anti-periodic BC)
#   - hN + 2h Σ_{k=1}^N a_k† a_k                           # chemical potential
#
# With j_c = 1, h = 0, this is the fermionic representation of the XY model
# with anti-periodic boundary conditions.
#
# For N=4, the expected ground state energy is E₀ = -4.

@testset "Fermionic Chain Ground State (N=4, anti-periodic BC)" begin
    T = Float64  # FermionicAlgebra uses real coefficients
    N = 4
    j_c = 1.0
    h = 0.0

    # Create fermionic variables (a = annihilation, a_dag = creation)
    registry, (a, a_dag) = create_fermionic_variables(1:N)

    # Build Hamiltonian with anti-periodic boundary conditions
    # Hopping terms: bulk + boundary
    ham_hop = sum(j_c * (a_dag[k] * a[k+1] + a_dag[k+1] * a[k]) for k in 1:N-1)
    ham_hop -= j_c * (a_dag[N] * a[1] + a_dag[1] * a[N])  # anti-periodic BC

    # Pairing terms: bulk + boundary
    ham_pair = sum(j_c * (a_dag[k] * a_dag[k+1] + a[k+1] * a[k]) for k in 1:N-1)
    ham_pair -= j_c * (a_dag[N] * a_dag[1] + a[1] * a[N])  # anti-periodic BC

    # Chemical potential terms (with h = 0, these vanish)
    ham_chem = -h * N + 2h * sum(a_dag[k] * a[k] for k in 1:N)

    # Total Hamiltonian
    ham = ham_hop + ham_pair + ham_chem

    # Create optimization problem
    pop = polyopt(ham, registry)

    # Solve with order-2 relaxation
    solver_config = SolverConfig(optimizer=SOLVER, order=2)
    res = cs_nctssos(pop, solver_config)

    # Refine with higher-order iterations
    res = cs_nctssos_higher(pop, res, solver_config)

    # Expected ground state energy
    E0_exact = -4.0

    println("Ground state energy lower bound: ", res.objective)
    println("Exact ground state energy: ", E0_exact)

    @test_broken res.objective ≈ E0_exact atol = 1e-5
end
