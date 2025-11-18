# Exact Diagonalization of the 3-Site Fermi-Hubbard Model

This directory contains a complete implementation of the Fermi-Hubbard model using exact diagonalization with explicit matrix representations of fermionic creation and annihilation operators.

## Overview

The Fermi-Hubbard model describes interacting fermions on a lattice:

```
H = -t Σ_⟨i,j⟩,σ (c†ᵢσ cⱼσ + c†ⱼσ cᵢσ) + U Σᵢ nᵢ↑ nᵢ↓
```

where:
- `t = 1.0` is the hopping parameter (kinetic energy)
- `U` is the on-site interaction strength (Coulomb repulsion)
- `c†ᵢσ, cᵢσ` are fermionic creation/annihilation operators
- `nᵢσ = c†ᵢσ cᵢσ` is the number operator

## Files

### 1. `fermionic_operator_matrices.jl` - Fundamental Operator Matrices

Demonstrates the matrix representation of fermionic operators in Fock space.

**Key Features:**
- Constructs explicit matrices for creation operator `c†ᵢ`
- Constructs explicit matrices for annihilation operator `cᵢ`
- Verifies anticommutation relations: `{cᵢ, c†ⱼ} = δᵢⱼ`
- Verifies Pauli exclusion: `(c†ᵢ)² = 0`
- Shows fermionic sign factors in action

**Example Output:**
```julia
Creation Operator c†₁ (2-mode system):
  |00⟩ |01⟩ |10⟩ |11⟩
|00⟩  0    0    1    0
|01⟩  0    0    0    1
|10⟩  0    0    0    0
|11⟩  0    0    0    0
```

**Run:**
```bash
julia --project examples/fermionic_operator_matrices.jl
```

### 2. `fermi_hubbard_3site_exact.jl` - Single Solution (U=1, t=1)

Solves the 3-site Hubbard model for the specific case U=1, t=1.

**System:**
- 3 sites, 2 spins → 6 fermionic modes
- Hilbert space dimension: 2^6 = 64
- Basis states: `|n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩`

**Results (U=1, t=1):**
```
Ground State Energy: E₀ = -2.3700171956
Energy per site:     E₀/site = -0.7900057319
Energy gap:          Δ = 0.2649229228
```

**Run:**
```bash
julia --project examples/fermi_hubbard_3site_exact.jl
```

### 3. `fermi_hubbard_3site_comparison.jl` - Comprehensive Analysis

Solves for multiple values of U to show the transition from free fermions to strong coupling.

**Parameter Scan:**
- U/t = 0.0 (free fermions)
- U/t = 0.5 (weak coupling)
- U/t = 1.0 (intermediate)
- U/t = 2.0
- U/t = 4.0
- U/t = 8.0 (strong coupling)

**Results Summary:**

| U/t  | E₀          | E₀/site     | Gap        | Physical Regime          |
|------|-------------|-------------|------------|--------------------------|
| 0.0  | -2.6762432  | -0.8920811  | 0.000000   | Free fermions            |
| 0.5  | -2.5092869  | -0.8364290  | 0.135114   | Weak coupling            |
| 1.0  | -2.3700172  | -0.7900057  | 0.264923   | Intermediate             |
| 2.0  | -2.1588286  | -0.7196095  | 0.493377   | Strong correlations      |
| 4.0  | -1.9097602  | -0.6365867  | 0.389001   | Approaching Mott         |
| 8.0  | -1.6991948  | -0.5663983  | 0.226975   | Mott insulator regime    |

**Run:**
```bash
julia --project examples/fermi_hubbard_3site_comparison.jl
```

## Matrix Representation Theory

### Fermionic Operators in Fock Space

For N fermionic modes, the Hilbert space has dimension 2^N with basis states `|n₁, n₂, ..., nₙ⟩` where `nᵢ ∈ {0, 1}`.

#### Creation Operator c†ᵢ

Acts on basis state `|...nᵢ...⟩`:
- If `nᵢ = 0`: `c†ᵢ|...0...⟩ = (-1)^Sᵢ |...1...⟩` where `Sᵢ = Σⱼ<ᵢ nⱼ`
- If `nᵢ = 1`: `c†ᵢ|...1...⟩ = 0` (Pauli exclusion)

#### Annihilation Operator cᵢ

Acts on basis state `|...nᵢ...⟩`:
- If `nᵢ = 1`: `cᵢ|...1...⟩ = (-1)^Sᵢ |...0...⟩`
- If `nᵢ = 0`: `cᵢ|...0...⟩ = 0`

#### Fermionic Sign Factor

The sign `(-1)^Sᵢ` ensures anticommutation relations:
- `{cᵢ, c†ⱼ} = cᵢc†ⱼ + c†ⱼcᵢ = δᵢⱼ`
- `{cᵢ, cⱼ} = 0`
- `{c†ᵢ, c†ⱼ} = 0`

### 3-Site Hubbard Model Details

#### Basis Indexing

Binary encoding: `|n₁↑, n₁↓, n₂↑, n₂↓, n₃↑, n₃↓⟩`

```
index = 1 + n₁↑ + 2·n₁↓ + 4·n₂↑ + 8·n₂↓ + 16·n₃↑ + 32·n₃↓
```

Examples:
- `|000000⟩` (vacuum) → index 1
- `|100000⟩` (one spin-up at site 1) → index 2
- `|111111⟩` (fully occupied) → index 64

#### Hamiltonian Matrix Structure

The Hamiltonian is a 64×64 Hermitian matrix:

**Diagonal elements:** Interaction energy
```
H[idx, idx] = U·(n₁↑·n₁↓ + n₂↑·n₂↓ + n₃↑·n₃↓)
```

**Off-diagonal elements:** Hopping terms
```
H[new_idx, idx] = -t · sign
```
where `sign = (-1)^(# of fermions between hopping sites)`

## Algebraic Rules Implementation

### 1. Anticommutation Relations

```julia
{c₁, c†₁} = I  (identity matrix)
{c₁, c†₂} = 0  (zero matrix)
{c₁, c₂}  = 0
{c†₁, c†₂} = 0
```

**Implementation:** Sign factor `(-1)^Sᵢ` where `Sᵢ` counts fermions to the left.

### 2. Pauli Exclusion Principle

```julia
(c†ᵢ)² = 0  (can't create twice in same mode)
(cᵢ)²  = 0  (can't annihilate twice from same mode)
```

**Implementation:** Conditional matrix elements (check occupation before acting).

### 3. Number Operator

```julia
nᵢ = c†ᵢ cᵢ
```

Diagonal matrix with eigenvalues 0 or 1 (occupation number).

### 4. Hermiticity

```julia
H = H†  (ensures real eigenvalues)
```

## Physical Interpretation

### U = 0 (Free Fermions)
- No interaction between electrons
- Pure kinetic energy minimization
- Electrons delocalize into molecular orbitals
- Ground state is 4-fold degenerate

### U = 1 (Intermediate Coupling)
- Balanced kinetic vs interaction energy
- Partial delocalization
- Beginning to avoid double occupancy
- Ground state is non-degenerate

### U → ∞ (Mott Insulator)
- Interaction energy dominates
- Double occupancy heavily suppressed
- Electrons localize to avoid interaction
- Charge gap opens

## Comparison with Symbolic Approach

This exact diagonalization approach complements the symbolic fermionic algebra approach in `fermi_hubbard_3site.jl`:

| Feature | Exact Diagonalization | Symbolic Algebra |
|---------|----------------------|------------------|
| Method | Matrix operators | Polynomial constraints |
| Basis | Fock states | Operator polynomials |
| Size | 64×64 matrix | 90 constraints |
| Output | Eigenvalues/eigenvectors | Optimization bounds |
| Accuracy | Numerically exact | SDP relaxation bounds |

Both approaches enforce the same algebraic rules:
- Anticommutation: `{cᵢ, c†ⱼ} = δᵢⱼ`
- Nilpotency: `cᵢ² = 0`
- Pauli exclusion: `(c†ᵢ)² = 0`

## Key Findings

1. **Matrix representation successfully implements fermionic algebra**
   - Sign factors enforce anticommutation
   - Conditional elements enforce Pauli exclusion

2. **Ground state energy evolution with U**
   - Increases (becomes less negative) as U increases
   - Smooth transition from metallic to insulating

3. **Energy gap behavior**
   - Non-monotonic: peaks around U ≈ 2t
   - Decreases at strong coupling (Mott regime)

4. **Computational efficiency**
   - Exact solution for small systems (≤6 modes)
   - Scales exponentially: 2^N Hilbert space dimension
   - 3 sites (64 states) diagonalizes in milliseconds

## References

1. Hubbard, J. (1963). "Electron Correlations in Narrow Energy Bands"
2. Lieb, E. H., & Wu, F. Y. (1968). "Absence of Mott Transition" - Phys. Rev. Lett.
3. Essler et al. (2005). "The One-Dimensional Hubbard Model"

## See Also

- `fermi_hubbard_3site.jl` - Symbolic algebra approach
- `scripts/compute_exact_hubbard.jl` - 2-site exact solutions
- `docs/src/examples/literate/hubbard_mpskit_groundstate.jl` - MPS methods
