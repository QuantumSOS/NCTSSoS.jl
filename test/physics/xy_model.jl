# XY Model Tests
# ===============
# XY model ground state energy tests using Pauli algebra.

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "setup.jl"))

@testset "XY Model Ground State (N=4, periodic BC)" begin
    T = ComplexF64
    N = 4
    j_c = 1.0

    # Create Pauli variables (σx, σy, σz for each site)
    registry, (x, y, z) = create_pauli_variables(1:N)

    # Build XY Hamiltonian with periodic boundary conditions
    # H = (j_c/4) Σᵢ (σxᵢ σxᵢ₊₁ + σyᵢ σyᵢ₊₁)
    ham = sum(T(j_c / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y] for i in 1:N)

    # Create optimization problem (Pauli simplification is automatic)
    pop = polyopt(ham, registry)

    # Solve with order-2 relaxation
    solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())
    res = cs_nctssos(pop, solver_config)

    # Refine with higher-order iterations
    res = cs_nctssos_higher(pop, res, solver_config)

    # Exact ground state energy for XY model with periodic BC
    # E₀ = -√2 for N=4, j_c=1
    E0_exact = -sqrt(2.0)

    @test res.objective ≈ E0_exact atol = 1e-5
end
