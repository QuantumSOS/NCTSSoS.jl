#!/usr/bin/env julia

"""
Fermi-Hubbard Model for 3-Site Chain with U=1, t=1

This script demonstrates solving the Fermi-Hubbard model on a 3-site chain lattice
using the fermionic algebra interface in NCTSSoS.jl.

## Mathematical Background

The Fermi-Hubbard Hamiltonian describes interacting fermions on a lattice:

H = -t Σ_⟨i,j⟩,σ (c†ᵢσ cⱼσ + c†ⱼσ cᵢσ) + U Σᵢ nᵢ↑ nᵢ↓

where:
- t = 1 is the hopping parameter (kinetic energy)
- U = 1 is the on-site interaction strength (Coulomb repulsion)
- c†ᵢσ, cᵢσ are creation/annihilation operators for site i, spin σ
- nᵢσ = c†ᵢσ cᵢσ is the number operator
- ⟨i,j⟩ denotes nearest-neighbor pairs

For a 3-site chain, we have nearest-neighbor pairs: (1,2) and (2,3)

## Fermionic Algebra Rules (Automatically Enforced)

The fermionic_algebra constructor automatically enforces:

1. Anti-commutation: {cᵢ, cⱼ} = 0, {c†ᵢ, c†ⱼ} = 0
2. Canonical anti-commutation: {cᵢ, c†ⱼ} = δᵢⱼ
3. Nilpotency: cᵢ² = 0, (c†ᵢ)² = 0

Total constraints: 2N² + 3N = 2(6)² + 3(6) = 90 constraints for 6 modes
"""

using NCTSSoS

println("="^70)
println("Fermi-Hubbard Model: 3-Site Chain with U=1, t=1")
println("="^70)

# Parameters
const t = 1.0  # Hopping parameter
const U = 1.0  # On-site interaction
const L = 2    # Number of lattice sites

println("\nSystem Parameters:")
println("  Lattice: 1D chain with $L sites")
println("  Hopping parameter: t = $t")
println("  Interaction strength: U = $U")
println("  Total fermionic modes: $(2*L) (spin-up + spin-down)")

# Create fermionic algebra system
# We need 2L modes: L sites × 2 spins
# Mode indexing:
#   - Modes 1,2: site 1 (spin-up, spin-down)
#   - Modes 3,4: site 2 (spin-up, spin-down)
#   - Modes 5,6: site 3 (spin-up, spin-down)
N_modes = 2 * L
sys = fermionic_algebra(N_modes)
c, c_dag = sys.variables

println("\nFermionic Algebra System Created:")
println("  Number of modes: $N_modes")
println("  Number of constraints: $(length(sys.equality_constraints))")
println("  Expected constraints (2N² + 3N): $(2*N_modes^2 + 3*N_modes)")

# Helper function: mode index for site i and spin σ
# σ = 1 for spin-up, σ = 2 for spin-down
mode_index(site::Int, spin::Int) = 2*(site-1) + spin

println("\nMode Indexing:")
for i in 1:L
    idx_up = mode_index(i, 1)
    idx_dn = mode_index(i, 2)
    println("  Site $i: spin-up → mode $idx_up, spin-down → mode $idx_dn")
end

# Construct the Hamiltonian
println("\nConstructing Hamiltonian...")

# Part 1: Hopping term
# H_hop = -t Σ_⟨i,j⟩,σ (c†ᵢσ cⱼσ + c†ⱼσ cᵢσ)
# For a chain: nearest neighbors are (1,2) and (2,3)

println("\n1. Hopping Terms:")
hopping_terms = []
for site in 1:(L-1)
    next_site = site + 1

    # Spin-up hopping
    i_up = mode_index(site, 1)
    j_up = mode_index(next_site, 1)

    # Spin-down hopping
    i_dn = mode_index(site, 2)
    j_dn = mode_index(next_site, 2)

    println("  Sites $site ↔ $next_site:")
    println("    Spin-up:   -t(c†[$i_up]c[$j_up] + c†[$j_up]c[$i_up])")
    println("    Spin-down: -t(c†[$i_dn]c[$j_dn] + c†[$j_dn]c[$i_dn])")

    # Collect hopping terms for both spins
    push!(hopping_terms, ComplexF64(-t) * (c_dag[i_up] * c[j_up] + c_dag[j_up] * c[i_up]))
    push!(hopping_terms, ComplexF64(-t) * (c_dag[i_dn] * c[j_dn] + c_dag[j_dn] * c[i_dn]))
end

# Part 2: Interaction term
# H_int = U Σᵢ nᵢ↑ nᵢ↓ = U Σᵢ (c†ᵢ↑ cᵢ↑)(c†ᵢ↓ cᵢ↓)

println("\n2. On-Site Interaction Terms:")
interaction_terms = []
for site in 1:L
    i_up = mode_index(site, 1)
    i_dn = mode_index(site, 2)

    # Number operators (already Hermitian)
    n_up = c_dag[i_up] * c[i_up]
    n_dn = c_dag[i_dn] * c[i_dn]

    println("  Site $site: U·n[$i_up]·n[$i_dn] = U·(c†[$i_up]c[$i_up])·(c†[$i_dn]c[$i_dn])")

    # Collect interaction term
    push!(interaction_terms, ComplexF64(U) * n_up * n_dn)
end

# Combine all terms
H = sum(hopping_terms) + sum(interaction_terms)

println("\n" * "="^70)
println("Hamiltonian Construction Complete")
println("="^70)

println("\nHamiltonian Summary:")
println("  Total terms: Hopping ($(2*(L-1)) pairs × 2 spins) + Interaction ($L sites)")
println("  H = H_hopping + H_interaction")
println("\nFull Hamiltonian Expression:")
println("  H = -t Σ_⟨i,j⟩,σ (c†ᵢσ cⱼσ + c†ⱼσ cᵢσ) + U Σᵢ nᵢ↑ nᵢ↓")

# Display the polynomial structure
println("\nPolynomial Details:")
println("  Number of monomials: $(length(H.monos))")
# Note: The polynomial contains hopping and interaction terms

println("\n" * "="^70)
println("Next Steps:")
println("="^70)
println("\nTo solve for the ground state energy:")

using MosekTools

# Create polynomial optimization problem
pop = cpolyopt(H, sys)

using NCTSSoS.FastPolynomials: simplify, star, terms
sa = sys.simplify_algo
sum(c * simplify(m, sa) for (c, m) in terms(H)) - sum(conj(c) * simplify(star(m), sa) for (c, m) in terms(H))

# currently, it does not know conj(c) = c^†

# Configure solver
solver_config = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 2  # Relaxation order (increase for higher accuracy)
)

# Solve
result = cs_nctssos(pop, solver_config)
println("Ground state energy: ", result.objective)

