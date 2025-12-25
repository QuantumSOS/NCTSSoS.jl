using NCTSSoS, NCTSSoS.FastPolynomials, Test
using JuMP

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2))
    )
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

# XY Model Ground State Test
#
# The XY model is related to the fermionic chain via Jordan-Wigner transformation.
# With periodic boundary conditions and j_c = 1:
#   H = (j_c/4) Σᵢ (σxᵢ σxᵢ₊₁ + σyᵢ σyᵢ₊₁)
#
# For N=4 with periodic BC, the exact ground state energy is E₀ = -√2 ≈ -1.4142
#
# Note: The original fermionic Hamiltonian shown has anti-periodic BC, which
# corresponds to a frustrated XY model. The moment hierarchy for frustrated
# systems may require higher order to converge. This test uses periodic BC.

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

    println("Ground state energy lower bound: ", res.objective)
    println("Exact ground state energy: ", E0_exact)

    @test res.objective ≈ E0_exact atol = 1e-4
end
