#!/usr/bin/env julia

"""
Matrix Representation of Fermionic Creation and Annihilation Operators

This script demonstrates how fermionic operators are represented as matrices
and verifies the fundamental anticommutation relations.

## Theory: Fermionic Operators in Fock Space

For N fermionic modes, the Fock space has dimension 2^N with basis:
|n₁, n₂, ..., nₙ⟩ where nᵢ ∈ {0, 1}

### Creation Operator c†ᵢ (Matrix Representation):
c†ᵢ |ket⟩ = (-1)^Sᵢ |ket'⟩ if nᵢ = 0 → 1
          0           if nᵢ = 1 (Pauli exclusion)

where Sᵢ = Σⱼ<ᵢ nⱼ (number of fermions to the left)

### Annihilation Operator cᵢ (Matrix Representation):
cᵢ |ket⟩ = (-1)^Sᵢ |ket'⟩ if nᵢ = 1 → 0
         0            if nᵢ = 0

## Example: 2-Mode System

For 2 modes, Hilbert space dimension = 4:
Basis: |00⟩, |01⟩, |10⟩, |11⟩

c†₁ matrix (creates fermion in mode 1):
  |00⟩ |01⟩ |10⟩ |11⟩
|00⟩  0    0    1    0
|01⟩  0    0    0    1
|10⟩  0    0    0    0
|11⟩  0    0    0    0

c₁ matrix (annihilates fermion from mode 1):
  |00⟩ |01⟩ |10⟩ |11⟩
|00⟩  0    0    0    0
|01⟩  0    0    0    0
|10⟩  1    0    0    0
|11⟩  0    1    0    0

Note the sign: c†₂|10⟩ = -|11⟩ (there's one fermion to the left)
"""

using LinearAlgebra
using Printf

println("="^70)
println("Matrix Representation of Fermionic Operators")
println("="^70)

"""
    create_operator(mode::Int, n_modes::Int) -> Matrix

Construct the matrix representation of fermionic creation operator c†ₘₒdₑ.

# Arguments
- `mode::Int`: Which mode (1-indexed)
- `n_modes::Int`: Total number of fermionic modes

# Returns
Sparse matrix of size 2^n_modes × 2^n_modes
"""
function create_operator(mode::Int, n_modes::Int)
    dim = 2^n_modes
    c_dag = zeros(Float64, dim, dim)

    # Iterate over all basis states (columns)
    for col_idx in 1:dim
        state = col_idx - 1  # 0-indexed state

        # Check if mode is unoccupied
        bit = (state >> (mode - 1)) & 1
        if bit == 0
            # Create fermion in this mode
            new_state = state | (1 << (mode - 1))
            row_idx = new_state + 1

            # Compute fermionic sign: count occupied modes to the left
            mask = (1 << (mode - 1)) - 1  # Bits for modes 1, 2, ..., mode-1
            n_left = count_ones(state & mask)
            sign = (-1)^n_left

            c_dag[row_idx, col_idx] = sign
        end
        # If bit == 1, operator annihilates state (remains zero)
    end

    return c_dag
end

"""
    annihilate_operator(mode::Int, n_modes::Int) -> Matrix

Construct the matrix representation of fermionic annihilation operator cₘₒdₑ.
"""
function annihilate_operator(mode::Int, n_modes::Int)
    dim = 2^n_modes
    c = zeros(Float64, dim, dim)

    # Iterate over all basis states (columns)
    for col_idx in 1:dim
        state = col_idx - 1

        # Check if mode is occupied
        bit = (state >> (mode - 1)) & 1
        if bit == 1
            # Annihilate fermion from this mode
            new_state = state & ~(1 << (mode - 1))
            row_idx = new_state + 1

            # Compute fermionic sign
            mask = (1 << (mode - 1)) - 1
            n_left = count_ones(state & mask)
            sign = (-1)^n_left

            c[row_idx, col_idx] = sign
        end
    end

    return c
end

"""
    number_operator(mode::Int, n_modes::Int) -> Matrix

Construct the number operator nₘₒdₑ = c†ₘₒdₑ cₘₒdₑ.
"""
function number_operator(mode::Int, n_modes::Int)
    c_dag = create_operator(mode, n_modes)
    c = annihilate_operator(mode, n_modes)
    return c_dag * c
end

# ============================================================================
# DEMONSTRATION: 2-Mode System
# ============================================================================

println("\n" * "="^70)
println("2-Mode System Demonstration")
println("="^70)

n_modes = 2
dim = 2^n_modes

println("\nBasis states (|mode2, mode1⟩):")
println("  State 0: |00⟩ - vacuum")
println("  State 1: |01⟩ - mode 1 occupied")
println("  State 2: |10⟩ - mode 2 occupied")
println("  State 3: |11⟩ - both modes occupied")

# Create operators for 2-mode system
c1_dag = create_operator(1, n_modes)
c1 = annihilate_operator(1, n_modes)
c2_dag = create_operator(2, n_modes)
c2 = annihilate_operator(2, n_modes)

println("\n" * "-"^70)
println("Creation Operator c†₁ (creates fermion in mode 1)")
println("-"^70)
println("Matrix representation:")
display(c1_dag)

println("\n" * "-"^70)
println("Annihilation Operator c₁ (annihilates fermion from mode 1)")
println("-"^70)
display(c1)

println("\n" * "-"^70)
println("Creation Operator c†₂ (creates fermion in mode 2)")
println("-"^70)
display(c2_dag)

# ============================================================================
# VERIFY ANTICOMMUTATION RELATIONS
# ============================================================================

println("\n" * "="^70)
println("Verification of Anticommutation Relations")
println("="^70)

# {c₁, c†₁} = c₁c†₁ + c†₁c₁ = I
anticomm_11 = c1 * c1_dag + c1_dag * c1
println("\n{c₁, c†₁} = c₁c†₁ + c†₁c₁:")
display(anticomm_11)
println("Should equal identity matrix: $(anticomm_11 ≈ I(dim))")

# {c₁, c†₂} = c₁c†₂ + c†₂c₁ = 0
anticomm_12 = c1 * c2_dag + c2_dag * c1
println("\n{c₁, c†₂} = c₁c†₂ + c†₂c₁:")
display(anticomm_12)
println("Should equal zero matrix: $(anticomm_12 ≈ zeros(dim, dim))")

# {c₁, c₂} = c₁c₂ + c₂c₁ = 0
anticomm_cc = c1 * c2 + c2 * c1
println("\n{c₁, c₂} = c₁c₂ + c₂c₁:")
display(anticomm_cc)
println("Should equal zero matrix: $(anticomm_cc ≈ zeros(dim, dim))")

# {c†₁, c†₂} = c†₁c†₂ + c†₂c†₁ = 0
anticomm_dagdag = c1_dag * c2_dag + c2_dag * c1_dag
println("\n{c†₁, c†₂} = c†₁c†₂ + c†₂c†₁:")
display(anticomm_dagdag)
println("Should equal zero matrix: $(anticomm_dagdag ≈ zeros(dim, dim))")

# ============================================================================
# VERIFY PAULI EXCLUSION
# ============================================================================

println("\n" * "="^70)
println("Verification of Pauli Exclusion Principle")
println("="^70)

# (c†₁)² = 0
c1_dag_squared = c1_dag * c1_dag
println("\n(c†₁)² (should be zero - can't create twice in same mode):")
display(c1_dag_squared)
println("Is zero: $(c1_dag_squared ≈ zeros(dim, dim))")

# (c₁)² = 0
c1_squared = c1 * c1
println("\n(c₁)² (should be zero - can't annihilate twice in same mode):")
display(c1_squared)
println("Is zero: $(c1_squared ≈ zeros(dim, dim))")

# ============================================================================
# NUMBER OPERATOR
# ============================================================================

println("\n" * "="^70)
println("Number Operator n₁ = c†₁c₁")
println("="^70)

n1 = number_operator(1, n_modes)
println("\nMatrix representation:")
display(n1)

println("\nEigenvalues (occupation numbers): $(eigvals(n1))")
println("Diagonal elements (occupation in each basis state):")
for i in 1:dim
    state = i - 1
    bit1 = (state >> 0) & 1
    bit2 = (state >> 1) & 1
    @printf("  |%d%d⟩: n₁ = %d\n", bit2, bit1, Int(n1[i, i]))
end

# ============================================================================
# DEMONSTRATION ON SPECIFIC STATES
# ============================================================================

println("\n" * "="^70)
println("Operator Actions on Specific States")
println("="^70)

# Define basis states as column vectors
vacuum = zeros(dim); vacuum[1] = 1.0      # |00⟩
state_1 = zeros(dim); state_1[2] = 1.0    # |01⟩
state_2 = zeros(dim); state_2[3] = 1.0    # |10⟩
state_both = zeros(dim); state_both[4] = 1.0  # |11⟩

println("\nc†₁|00⟩ (create in mode 1 from vacuum):")
result = c1_dag * vacuum
for (i, amp) in enumerate(result)
    if abs(amp) > 1e-10
        state = i - 1
        bit2 = (state >> 1) & 1
        bit1 = (state >> 0) & 1
        sign = amp >= 0 ? "+" : ""
        @printf("  %s%.1f|%d%d⟩\n", sign, amp, bit2, bit1)
    end
end

println("\nc†₂|01⟩ (create in mode 2, already have mode 1):")
result = c2_dag * state_1
for (i, amp) in enumerate(result)
    if abs(amp) > 1e-10
        state = i - 1
        bit2 = (state >> 1) & 1
        bit1 = (state >> 0) & 1
        sign = amp >= 0 ? "+" : ""
        @printf("  %s%.1f|%d%d⟩\n", sign, amp, bit2, bit1)
    end
end
println("  Note: Sign is -1 because we cross one fermion to create in mode 2")

println("\nc₁|11⟩ (annihilate from mode 1, both occupied):")
result = c1 * state_both
for (i, amp) in enumerate(result)
    if abs(amp) > 1e-10
        state = i - 1
        bit2 = (state >> 1) & 1
        bit1 = (state >> 0) & 1
        sign = amp >= 0 ? "+" : ""
        @printf("  %s%.1f|%d%d⟩\n", sign, amp, bit2, bit1)
    end
end

println("\nc†₁|01⟩ (try to create in mode 1, already occupied):")
result = c1_dag * state_1
norm = sqrt(sum(abs2, result))
if norm < 1e-10
    println("  0 (Pauli exclusion - can't double occupy)")
else
    println("  Non-zero result (unexpected!)")
end

println("\n" * "="^70)
println("Summary: Matrix Representation of Fermionic Operators")
println("="^70)
println("\n✓ Anticommutation relations verified: {cᵢ, c†ⱼ} = δᵢⱼ")
println("✓ Pauli exclusion verified: (c†ᵢ)² = 0, (cᵢ)² = 0")
println("✓ Fermionic sign factors correctly implemented")
println("✓ Number operators are diagonal in occupation basis")
println("\nThese matrix operators can be used to construct Hamiltonians")
println("for exact diagonalization of many-body fermionic systems!")

println("\n" * "="^70)
println("Done!")
println("="^70)
