#!/usr/bin/env julia

"""
Exact Diagonalization: 3-Site Fermi-Hubbard Model with Matrix Representation

This script demonstrates the matrix representation of fermionic creation and annihilation
operators and uses exact diagonalization to solve the 3-site Fermi-Hubbard model.

## Matrix Representation of Fermionic Operators

For a system with N fermionic modes, the Fock space has dimension 2^N.
Basis states are |n₁, n₂, ..., nₙ⟩ where nᵢ ∈ {0, 1}.

### Creation Operator c†ᵢ:
- Acts on |...nᵢ...⟩:
  * If nᵢ = 0: c†ᵢ|...0...⟩ = (-1)^Sᵢ |...1...⟩  where Sᵢ = Σⱼ<ᵢ nⱼ
  * If nᵢ = 1: c†ᵢ|...1...⟩ = 0 (Pauli exclusion)

### Annihilation Operator cᵢ:
- Acts on |...nᵢ...⟩:
  * If nᵢ = 1: cᵢ|...1...⟩ = (-1)^Sᵢ |...0...⟩
  * If nᵢ = 0: cᵢ|...0...⟩ = 0

The sign factor (-1)^Sᵢ ensures fermionic anticommutation relations:
{cᵢ, c†ⱼ} = δᵢⱼ and {cᵢ, cⱼ} = 0

## For 3-Site Hubbard Model:
- 6 fermionic modes: 3 sites × 2 spins
- Hilbert space dimension: 2^6 = 64
- Basis: |n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩
"""

using LinearAlgebra
using Printf

println("="^70)
println("Exact Diagonalization: 3-Site Fermi-Hubbard Model")
println("="^70)

"""
    hubbard_3site_hamiltonian(t::Float64, U::Float64)

Build the Hamiltonian matrix for 3-site Fermi-Hubbard model using exact
matrix representation of fermionic operators.

# Arguments
- `t::Float64`: Hopping parameter
- `U::Float64`: On-site interaction strength

# Basis States
64 total states for 3 sites × 2 spins:
|n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩ where nᵢσ ∈ {0,1}

# Indexing
Binary indexing: idx = 1 + n₁↑ + 2*n₁↓ + 4*n₂↑ + 8*n₂↓ + 16*n₃↑ + 32*n₃↓

# Returns
Hermitian matrix representing H = H_hop + H_int where:
- H_hop = -t Σ_⟨i,j⟩,σ (c†ᵢσ cⱼσ + h.c.)
- H_int = U Σᵢ nᵢ↑ nᵢ↓
"""
function hubbard_3site_hamiltonian(t::Float64, U::Float64)
    # 64 basis states for 3 sites × 2 spins
    dim = 64
    H = zeros(Float64, dim, dim)

    # Helper: convert occupation numbers to basis index
    # State |n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩ → index
    to_index(n1up, n1dn, n2up, n2dn, n3up, n3dn) =
        1 + n1up + 2*n1dn + 4*n2up + 8*n2dn + 16*n3up + 32*n3dn

    # Iterate over all basis states
    for idx in 1:dim
        # Decode state from index (binary representation)
        i = idx - 1
        n1up = i & 1
        n1dn = (i >> 1) & 1
        n2up = (i >> 2) & 1
        n2dn = (i >> 3) & 1
        n3up = (i >> 4) & 1
        n3dn = (i >> 5) & 1

        # ====================================================================
        # INTERACTION TERM: U Σᵢ nᵢ↑ nᵢ↓
        # ====================================================================
        # Diagonal term - measures double occupancy at each site
        H[idx, idx] += U * n1up * n1dn  # Site 1
        H[idx, idx] += U * n2up * n2dn  # Site 2
        H[idx, idx] += U * n3up * n3dn  # Site 3

        # ====================================================================
        # HOPPING TERM: -t Σ_⟨i,j⟩,σ (c†ᵢσ cⱼσ + c†ⱼσ cᵢσ)
        # ====================================================================
        # For a chain: nearest neighbors are (1,2) and (2,3)

        # --- Bond 1-2: Spin-up hopping ---
        # c†₁↑ c₂↑: removes electron from site 2↑, creates at site 1↑
        if n1up == 0 && n2up == 1
            # Fermionic sign: count occupied states to the left of mode 0 (site 1↑)
            # Modes to the left: none
            # Then count between mode 0 and mode 2 (site 2↑): n₁↓
            sign = (-1)^(n1dn)
            new_idx = to_index(1, n1dn, 0, n2dn, n3up, n3dn)
            H[new_idx, idx] += -t * sign
        end

        # c†₂↑ c₁↑: removes electron from site 1↑, creates at site 2↑
        if n2up == 0 && n1up == 1
            # Count occupied states between mode 0 (site 1↑) and mode 2 (site 2↑)
            sign = (-1)^(n1dn)
            new_idx = to_index(0, n1dn, 1, n2dn, n3up, n3dn)
            H[new_idx, idx] += -t * sign
        end

        # --- Bond 1-2: Spin-down hopping ---
        # c†₁↓ c₂↓: removes electron from site 2↓, creates at site 1↓
        if n1dn == 0 && n2dn == 1
            # Count occupied states to the left of mode 1 (site 1↓): n₁↑
            # Then count between mode 1 and mode 3 (site 2↓): n₂↑
            sign = (-1)^(n1up + n2up)
            new_idx = to_index(n1up, 1, n2up, 0, n3up, n3dn)
            H[new_idx, idx] += -t * sign
        end

        # c†₂↓ c₁↓: removes electron from site 1↓, creates at site 2↓
        if n2dn == 0 && n1dn == 1
            sign = (-1)^(n1up + n2up)
            new_idx = to_index(n1up, 0, n2up, 1, n3up, n3dn)
            H[new_idx, idx] += -t * sign
        end

        # --- Bond 2-3: Spin-up hopping ---
        # c†₂↑ c₃↑: removes electron from site 3↑, creates at site 2↑
        if n2up == 0 && n3up == 1
            # Count occupied states between mode 2 (site 2↑) and mode 4 (site 3↑)
            # These are: n₂↓
            sign = (-1)^(n2dn)
            new_idx = to_index(n1up, n1dn, 1, n2dn, 0, n3dn)
            H[new_idx, idx] += -t * sign
        end

        # c†₃↑ c₂↑: removes electron from site 2↑, creates at site 3↑
        if n3up == 0 && n2up == 1
            sign = (-1)^(n2dn)
            new_idx = to_index(n1up, n1dn, 0, n2dn, 1, n3dn)
            H[new_idx, idx] += -t * sign
        end

        # --- Bond 2-3: Spin-down hopping ---
        # c†₂↓ c₃↓: removes electron from site 3↓, creates at site 2↓
        if n2dn == 0 && n3dn == 1
            # Count occupied states between mode 3 (site 2↓) and mode 5 (site 3↓)
            # These are: n₃↑
            sign = (-1)^(n3up)
            new_idx = to_index(n1up, n1dn, n2up, 1, n3up, 0)
            H[new_idx, idx] += -t * sign
        end

        # c†₃↓ c₂↓: removes electron from site 2↓, creates at site 3↓
        if n3dn == 0 && n2dn == 1
            sign = (-1)^(n3up)
            new_idx = to_index(n1up, n1dn, n2up, 0, n3up, 1)
            H[new_idx, idx] += -t * sign
        end
    end

    return Hermitian(H)
end

println("\n1. System: 3-site Fermi-Hubbard chain")
println("   Hilbert space dimension: 64 (2^6 basis states)")
println("   Basis: |n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩ with nᵢσ ∈ {0,1}")
println("   Bonds: (1,2) and (2,3)")

# Parameters
t = 1.0
U = 1.0

println("\n2. Parameters:")
println("   Hopping parameter: t = $t")
println("   Interaction strength: U = $U")

# Build Hamiltonian
println("\n3. Building Hamiltonian matrix...")
H = hubbard_3site_hamiltonian(t, U)
println("   Matrix size: $(size(H, 1)) × $(size(H, 2))")
println("   Matrix type: $(typeof(H))")
println("   Is Hermitian: $(ishermitian(H))")

# Diagonalize
println("\n4. Diagonalizing Hamiltonian...")
eigenvalues = eigvals(H)

# Find ground state
E0 = minimum(eigenvalues)
E0_per_site = E0 / 3

println("\n" * "="^70)
println("RESULTS")
println("="^70)
println("\nGround State Energy:")
@printf("  E₀       = %.10f\n", E0)
@printf("  E₀/site  = %.10f\n", E0_per_site)

# Additional analysis
println("\n" * "="^70)
println("Energy Spectrum Analysis")
println("="^70)

# Sort eigenvalues
sorted_eigs = sort(eigenvalues)

println("\nLowest 10 eigenvalues:")
for i in 1:min(10, length(sorted_eigs))
    @printf("  E[%2d] = %12.8f\n", i-1, sorted_eigs[i])
end

println("\nHighest 5 eigenvalues:")
for i in (length(sorted_eigs)-4):length(sorted_eigs)
    @printf("  E[%2d] = %12.8f\n", i-1, sorted_eigs[i])
end

# Energy gap
gap = sorted_eigs[2] - sorted_eigs[1]
println("\nEnergy gap (E₁ - E₀): $gap")

println("\n" * "="^70)
println("Physical Interpretation")
println("="^70)
println("\nFor U=$U, t=$t (intermediate coupling regime):")
println("  • Ground state balances kinetic and interaction energy")
println("  • Energy gap indicates nature of excitations")
println("  • Negative energy indicates bound state formation")

println("\n" * "="^70)
println("Matrix Representation Details")
println("="^70)
println("\nFermionic Operator Matrix Construction:")
println("  • Basis states: |n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩")
println("  • Binary indexing: idx = Σᵢ nᵢ * 2^i")
println("  • Sign factors: (-1)^(# fermions to the left)")
println("  • Anticommutation: {cᵢ, c†ⱼ} = δᵢⱼ")
println("  • Pauli exclusion: c†ᵢ|1⟩ = 0, cᵢ|0⟩ = 0")

println("\n" * "="^70)
println("Done!")
println("="^70)
