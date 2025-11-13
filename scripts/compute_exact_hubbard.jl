#!/usr/bin/env julia

"""
Compute exact ground state energies for 2-site Fermi-Hubbard model
using direct matrix construction and LinearAlgebra

This provides baseline reference values for testing the NCTSSoS fermionic algebra implementation.
"""

using LinearAlgebra
using Printf

println("="^70)
println("Exact Diagonalization: 2-Site Fermi-Hubbard Model")
println("="^70)

"""
Build the Hamiltonian matrix for 2-site Fermi-Hubbard model.

Basis states (16 total for 2 sites with spin):
|n₁↑, n₁↓, n₂↑, n₂↓⟩ where nᵢσ ∈ {0,1}

Indexing: state |n₁↑, n₁↓, n₂↑, n₂↓⟩ → index = 1 + n₁↑ + 2*n₁↓ + 4*n₂↑ + 8*n₂↓
"""
function hubbard_2site_hamiltonian(t::Float64, U::Float64; periodic::Bool=true)
    # 16 basis states for 2 sites × 2 spins
    dim = 16
    H = zeros(Float64, dim, dim)

    # Helper: convert (n1up, n1dn, n2up, n2dn) to index
    to_index(n1up, n1dn, n2up, n2dn) = 1 + n1up + 2*n1dn + 4*n2up + 8*n2dn

    # Iterate over all basis states
    for idx in 1:dim
        # Decode state
        i = idx - 1
        n1up = i & 1
        n1dn = (i >> 1) & 1
        n2up = (i >> 2) & 1
        n2dn = (i >> 3) & 1

        # Interaction term: U Σᵢ nᵢ↑ nᵢ↓
        H[idx, idx] += U * n1up * n1dn
        H[idx, idx] += U * n2up * n2dn

        # Hopping term: -t Σσ (c†₁σ c₂σ + c†₂σ c₁σ)
        # For spin-up: c†₁↑ c₂↑
        if n1up == 0 && n2up == 1
            # Apply c†₁↑ c₂↑ to state
            sign = (-1)^(n1dn)  # Fermionic sign from anticommutation
            new_idx = to_index(1, n1dn, 0, n2dn)
            H[new_idx, idx] += -t * sign
        end

        # For spin-up: c†₂↑ c₁↑
        if n2up == 0 && n1up == 1
            sign = (-1)^(n1dn + n2dn)
            new_idx = to_index(0, n1dn, 1, n2dn)
            H[new_idx, idx] += -t * sign
        end

        # For spin-down: c†₁↓ c₂↓
        if n1dn == 0 && n2dn == 1
            sign = (-1)^(n1up + n2up)
            new_idx = to_index(n1up, 1, n2up, 0)
            H[new_idx, idx] += -t * sign
        end

        # For spin-down: c†₂↓ c₁↓
        if n2dn == 0 && n1dn == 1
            sign = (-1)^(n1up + n2up + n2dn)
            new_idx = to_index(n1up, 0, n2up, 1)
            H[new_idx, idx] += -t * sign
        end
    end

    return Hermitian(H)
end

println("\n1. System: 2-site Fermi-Hubbard chain (periodic boundary conditions)")
println("   Hilbert space dimension: 16 (4 states per site)")
println("   Basis: |n₁↑, n₁↓, n₂↑, n₂↓⟩ with nᵢσ ∈ {0,1}")

# Parameters to test
t_values = [1.0]  # Hopping parameter
U_values = [0.0, 1.0, 4.0, 8.0]  # Interaction strength

println("\n2. Computing ground state energies...")
println("-"^70)
@printf("%-8s %-8s %-15s %-15s\n", "t", "U", "E₀", "E₀/site")
println("-"^70)

# Store results
results = Dict()

for t in t_values
    for U in U_values
        H = hubbard_2site_hamiltonian(t, U)

        # Diagonalize
        eigenvalues = eigvals(H)

        # Ground state energy
        E0 = minimum(eigenvalues)
        E0_per_site = E0 / 2

        results[(t, U)] = (E0=E0, E0_per_site=E0_per_site, eigenvalues=eigenvalues)

        @printf("%-8.1f %-8.1f %-15.6f %-15.6f\n", t, U, E0, E0_per_site)
    end
end

println("-"^70)

# Summary section
println("\n3. Summary of Exact Ground State Energies")
println("="^70)
println("\nFor NCTSSoS Testing:")
println("-"^70)

for U in U_values
    E0 = results[(1.0, U)].E0
    println("  U = $U*t: E₀ = $E0")
end

println("\n" * "="^70)
println("Physical Interpretation:")
println("="^70)

E0_free = results[(1.0, 0.0)].E0
E0_U4 = results[(1.0, 4.0)].E0
E0_U8 = results[(1.0, 8.0)].E0

println("\n• U=0 (Free fermions):")
println("  E₀ = $E0_free")
println("  Particles delocalize, occupy bonding orbital")

println("\n• U=4t (Intermediate interaction):")
println("  E₀ = $E0_U4")
println("  Competition between hopping and Coulomb repulsion")

println("\n• U=8t (Strong interaction):")
println("  E₀ = $E0_U8")
println("  Approaching Mott insulator, particles localize")

println("\n" * "="^70)
println("Done! These values can be used as exact benchmarks.")
println("="^70)
