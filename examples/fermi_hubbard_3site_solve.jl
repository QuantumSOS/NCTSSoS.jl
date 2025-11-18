#!/usr/bin/env julia

"""
Solve the Fermi-Hubbard Interaction Term for 3-Site Chain

This script demonstrates solving the interaction part of the Fermi-Hubbard model
using the NCTSSoS polynomial optimization framework with fermionic algebra.

Note: The full Hubbard Hamiltonian (hopping + interaction) requires special
symmetrization techniques. This example focuses on the interaction term which
demonstrates the fermionic algebra rules and SDP solving.

For the full Hamiltonian, see:
- examples/fermi_hubbard_3site.jl (construction)
- scripts/compute_exact_hubbard.jl (exact diagonalization)
- docs/src/examples/literate/hubbard_mpskit_groundstate.jl (MPS methods)
"""

using NCTSSoS, MosekTools

function solve_hubbard_interaction()
    println("="^70)
    println("Solving Hubbard Interaction Term: 3-Site Chain with U=1")
    println("="^70)

    # Parameters
    t = 1.0
    U = 1.0
    L = 3

    println("\nSystem Parameters:")
    println("  Lattice: 1D chain with $L sites")
    println("  Hopping parameter: t = $t")
    println("  Interaction strength: U = $U")

    # Create fermionic algebra system (6 modes: 3 sites × 2 spins)
    N_modes = 2 * L
    sys = fermionic_algebra(N_modes)
    c, c_dag = sys.variables

    println("  Total fermionic modes: $N_modes")
    println("  Fermionic constraints: $(length(sys.equality_constraints))")

    # Construct the Interaction Hamiltonian
    # Mode indexing: 1:L for spin-up, (L+1):2L for spin-down
    println("\nConstructing interaction Hamiltonian...")
    println("  H_int = U Σ_i n_i↑ n_i↓")

    # Interaction term: H_int = U Σ_i n_i↑ n_i↓
    # Use symmetrized number operators as shown in test/algebra_constructors.jl:492-494
    H = ComplexF64(0.0)
    for site in 1:L
        i_up = site
        i_dn = site + L
        # Symmetrized number operators
        n_up = ComplexF64(0.5) * (c_dag[i_up] * c[i_up] + c[i_up] * c_dag[i_up])
        n_dn = ComplexF64(0.5) * (c_dag[i_dn] * c[i_dn] + c[i_dn] * c_dag[i_dn])
        H += U * n_up * n_dn
    end

    println("  Hamiltonian has $(length(H.monos)) monomials")

    # Create polynomial optimization problem
    println("\nSetting up optimization problem...")
    pop = cpolyopt(H, sys)

    # Configure solver
    println("Configuring SDP solver...")
    solver_config = SolverConfig(
        optimizer = Mosek.Optimizer,
        order = 1  # Relaxation order (1 is sufficient for interaction term)
    )

    println("  Solver: Mosek")
    println("  Relaxation order: 1")

    # Solve for minimum interaction energy
    println("\nSolving for minimum interaction energy...")
    println("  (This may take a few moments)")

    result = cs_nctssos(pop, solver_config)

    println("\n" * "="^70)
    println("RESULTS")
    println("="^70)
    println("\nMinimum Interaction Energy:")
    println("  E_int = $(result.objective)")
    println("  E_int/site = $(result.objective / L)")

    println("\n" * "="^70)
    println("Physical Interpretation:")
    println("="^70)
    println("\nFor the Hubbard interaction term with U=$U:")
    println("  - Minimum energy occurs when double occupancy is minimized")
    println("  - Ground state avoids having spin-up and spin-down on same site")
    println("  - Result should be non-negative (E ≥ 0)")

    println("\n" * "="^70)
    println("Algebraic Rules (Enforced by fermionic_algebra):")
    println("="^70)
    println("\n1. Anti-commutation: {cᵢ, cⱼ} = 0, {c†ᵢ, c†ⱼ} = 0")
    println("2. Canonical: {cᵢ, c†ⱼ} = δᵢⱼ")
    println("3. Nilpotency: cᵢ² = 0, (c†ᵢ)² = 0")
    println("4. Pauli Exclusion: No two fermions in same state")

    println("\n" * "="^70)
    println("Done!")
    println("="^70)

    return result
end

# Run the solver
solve_hubbard_interaction()
