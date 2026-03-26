# Demo: 4-site Fermi-Hubbard model — ground state energy and correlations
#
# Persona: condensed matter physicist, no optimization background.
# Tests: full pipeline (problem ID → algebra → formulation → solve),
#        FermionicAlgebra, correlation extraction from moments.
#
# Problem:
#   H = -t Σᵢ(c†ᵢσ cᵢ₊₁σ + h.c.) + U Σᵢ nᵢ↑ nᵢ↓
#   4 sites, periodic boundary conditions, t=1, U=4, half-filling.
#
# Reference: DMRG gives E₀ ≈ -2.83 for these parameters.

using NCTSSoS
using Clarabel

# --- Parameters ---
N = 4        # number of sites
t = 1.0      # hopping
U = 4.0      # on-site interaction

# --- Variables ---
# Spin-up and spin-down fermions on each site share one registry but use
# disjoint fermionic modes, so nᵢ↑ nᵢ↓ stays quartic instead of collapsing to nᵢ.
registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:N),
    ("c_dn", 1:N),
])

# --- Hamiltonian ---
# Hopping: -t Σᵢσ (c†ᵢσ cᵢ₊₁σ + h.c.), periodic boundary
ham_hop_up = -t * sum(
    c_up_dag[i] * c_up[mod1(i+1, N)] + c_up_dag[mod1(i+1, N)] * c_up[i]
    for i in 1:N
)
ham_hop_dn = -t * sum(
    c_dn_dag[i] * c_dn[mod1(i+1, N)] + c_dn_dag[mod1(i+1, N)] * c_dn[i]
    for i in 1:N
)

# Interaction: U Σᵢ nᵢ↑ nᵢ↓ where nᵢσ = c†ᵢσ cᵢσ
ham_int = U * sum(
    (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
    for i in 1:N
)

# Total Hamiltonian
ham = ham_hop_up + ham_hop_dn + ham_int

# --- Solve ---
# 4 sites × 2 spins = 8 fermionic modes → 16 variables (c + c†)
# order=2 is the standard starting point
pop = polyopt(ham, registry)
config = SolverConfig(optimizer=Clarabel.Optimizer, order=2)
result = cs_nctssos(pop, config)

println("Ground state energy bound: ", result.objective)
println("Expected (DMRG):           ≈ -2.83")

# --- Correlation: ⟨Sᶻ₁ Sᶻ₂⟩ ---
# Sᶻᵢ = (nᵢ↑ - nᵢ↓)/2 = (c†ᵢ↑cᵢ↑ - c†ᵢ↓cᵢ↓)/2
# ⟨Sᶻ₁ Sᶻ₂⟩ can be read from the solved moment map.
# After solving, the moment map contains ⟨m⟩ for each monomial m in the basis.
#
# To extract:
#   Sz1_Sz2 = (1/4) * (⟨n₁↑n₂↑⟩ - ⟨n₁↑n₂↓⟩ - ⟨n₁↓n₂↑⟩ + ⟨n₁↓n₂↓⟩)
#
# This requires accessing result's moment data — see the postprocessing guide
# in the polyopt-guide skill for the low-level moment map API.
