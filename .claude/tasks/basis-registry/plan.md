# Implementation Plan: Registry-Aware Basis Generation

## Objective
Rewrite `basis.jl` to use `VariableRegistry` for type consistency and correct variable indices.

## Current Problems
1. `basis.jl` defaults to `T=Int`, but `@ncpolyvar` uses `UInt64` via VariableRegistry
2. Generates indices 1:n instead of actual variable indices from registry
3. Type mismatch when comparing/sorting monomials from different sources

## New API Design

### Primary Functions
```julia
# Replace old API entirely
get_ncbasis(registry::VariableRegistry{A,T}, d::Int) -> Vector{Monomial{A,T}}
get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) -> Vector{Monomial{A,T}}
```

### Internal Dispatch (algebra-specific)
```julia
_generate_canonical_words_deg(::Type{PauliAlgebra}, indices::Vector{T}, d) where T
_generate_canonical_words_deg(::Type{FermionicAlgebra}, indices::Vector{T}, d) where T
_generate_canonical_words_deg(::Type{BosonicAlgebra}, indices::Vector{T}, d) where T
_generate_canonical_words_deg(::Type{UnipotentAlgebra}, indices::Vector{T}, d) where T
_generate_canonical_words_deg(::Type{ProjectorAlgebra}, indices::Vector{T}, d) where T
_generate_canonical_words_deg(::Type{NonCommutativeAlgebra}, indices::Vector{T}, d) where T
```

## Core Algorithm: Canonical-Direct Generation

**CRITICAL**: Generate ONLY canonical-form words directly. Do NOT generate all words then simplify.

### Why NOT "generate all → simplify → dedupe"?
1. **Degree Mixing**: Simplification produces sums of different degrees (e.g., a₁a₁† → 1 - a₁†a₁)
   - These belong in DIFFERENT basis blocks for SOS optimization
   - Including degree-0 terms in degree-2 basis breaks moment matrix structure
2. **Exponential Redundancy**: Most non-canonical words reduce to existing canonical forms
3. **Linear Independence**: Including both non-canonical and canonical forms violates basis requirements

### Correct Algorithm
```julia
function _generate_basis_deg(::Type{A}, indices::Vector{T}, d::Int) where {A,T}
    d == 0 && return [Monomial{A}(T[])]  # identity
    d < 0 && return Monomial{A,T}[]

    # Directly generate canonical words of degree d (algebra-specific)
    canonical_words = _generate_canonical_words_deg(A, indices, d)
    return sort!([Monomial{A}(w) for w in canonical_words])
end
```

## Algebra-Specific Canonical Generation

### NonCommutativeAlgebra (simplest case)
- **Canonical form**: All words (no simplification rules)
- **Algorithm**: Generate all n^d words of degree d
```julia
function _generate_canonical_words_deg(::Type{NonCommutativeAlgebra}, indices::Vector{T}, d::Int) where T
    # All words are canonical - no reduction
    return _generate_all_words(indices, d)
end
```

### PauliAlgebra (UInt16, self-adjoint)
- **Canonical form**: No σᵢ² patterns (since σᵢ² = I reduces degree)
- **Indices**: unsigned, represent σx, σy, σz at each site
- **Algorithm**: Generate words avoiding consecutive identical indices
```julia
function _generate_canonical_words_deg(::Type{PauliAlgebra}, indices::Vector{T}, d::Int) where T
    # Generate words without consecutive repeats (σᵢ² = I)
    return _generate_no_consecutive_repeats(indices, d)
end
```

### UnipotentAlgebra (UInt16, self-adjoint)
- **Canonical form**: No consecutive repeats (U² = I)
- **Algorithm**: Same as Pauli
```julia
function _generate_canonical_words_deg(::Type{UnipotentAlgebra}, indices::Vector{T}, d::Int) where T
    return _generate_no_consecutive_repeats(indices, d)
end
```

### ProjectorAlgebra (UInt16, self-adjoint)
- **Canonical form**: No consecutive repeats (P² = P reduces degree)
- **Algorithm**: Same as Pauli/Unipotent
```julia
function _generate_canonical_words_deg(::Type{ProjectorAlgebra}, indices::Vector{T}, d::Int) where T
    return _generate_no_consecutive_repeats(indices, d)
end
```

### FermionicAlgebra (Int32, signed)
- **Canonical form**: Normal order (creators before annihilators, sorted by mode index)
- **Indices**: +i = annihilation aᵢ, -i = creation aᵢ†
- **Constraints**:
  - No repeated indices (Pauli exclusion: aᵢaᵢ = 0, aᵢ†aᵢ† = 0)
  - Creators (negative indices) come first, then annihilators (positive)
  - Within each group, sorted by |index|
```julia
function _generate_canonical_words_deg(::Type{FermionicAlgebra}, indices::Vector{T}, d::Int) where T
    # Extract modes (unsigned values)
    modes = unique(abs.(indices))

    # Generate all subsets of creators and annihilators totaling degree d
    # Format: [−m₁, −m₂, ..., +n₁, +n₂, ...] where mᵢ < mⱼ and nᵢ < nⱼ
    return _generate_fermionic_normal_order(modes, d)
end
```

### BosonicAlgebra (Int32, signed)
- **Canonical form**: Normal order (creators before annihilators)
- **Indices**: +i = annihilation bᵢ, -i = creation bᵢ†
- **Constraints**:
  - Repeated indices allowed (no Pauli exclusion for bosons)
  - Creators (negative) come first, then annihilators (positive)
```julia
function _generate_canonical_words_deg(::Type{BosonicAlgebra}, indices::Vector{T}, d::Int) where T
    modes = unique(abs.(indices))
    # Generate multisets with repetition allowed
    return _generate_bosonic_normal_order(modes, d)
end
```

## Implementation Steps

### Step 1: Add algebra type to VariableRegistry ✓ (validated)
- [ ] Change `VariableRegistry{T}` to `VariableRegistry{A,T}` where A is AlgebraType
- [ ] Update all `create_*_variables` functions:
  - `create_pauli_variables` → returns `VariableRegistry{PauliAlgebra,T}`
  - `create_fermionic_variables` → returns `VariableRegistry{FermionicAlgebra,T}`
  - `create_bosonic_variables` → returns `VariableRegistry{BosonicAlgebra,T}`
  - `create_projector_variables` → returns `VariableRegistry{ProjectorAlgebra,T}`
  - `create_unipotent_variables` → returns `VariableRegistry{UnipotentAlgebra,T}`
  - `create_noncommutative_variables` → returns `VariableRegistry{NonCommutativeAlgebra,T}`
- [ ] Update internal helpers: `_create_physical_variables`, `_create_noncommutative_variables`
- [ ] Add helper functions:
  ```julia
  algebra_type(::VariableRegistry{A,T}) where {A,T} = A
  index_type(::VariableRegistry{A,T}) where {A,T} = T
  ```

### Step 2: Rewrite `get_ncbasis` and `get_ncbasis_deg`
- [ ] Change signature to take `VariableRegistry{A,T}`
- [ ] Use `indices(registry)` to get all variable indices
- [ ] Infer both A and T from registry type parameters
- [ ] Call `_generate_basis_deg(A, indices, d)`

### Step 3: Implement canonical word generators
- [ ] `_generate_all_words(indices, d)` - for NonCommutativeAlgebra
- [ ] `_generate_no_consecutive_repeats(indices, d)` - for Pauli/Unipotent/Projector
- [ ] `_generate_fermionic_normal_order(modes, d)` - for FermionicAlgebra
- [ ] `_generate_bosonic_normal_order(modes, d)` - for BosonicAlgebra

### Step 4: Remove old API
- [ ] Remove `filter_constraint` parameter (replaced by algebra dispatch)
- [ ] Remove `T` type parameter (inferred from registry)
- [ ] Remove `n` parameter (from registry)
- [ ] Remove or update `monomials` alias
- [ ] Keep `has_consecutive_repeats` as internal utility

### Step 5: Update tests
- [ ] Rewrite tests to use registry-based API
- [ ] Unskip the 25 testsets that were skipped due to type mismatch
- [ ] Add tests for fermionic/bosonic normal-order generation
- [ ] Verify degree separation (degree-d basis contains ONLY degree-d monomials)

### Step 6: Update dependent code
- [ ] Find all usages of old `get_ncbasis` API
- [ ] Update to use registry-based API

## Resolved Questions
1. ~~Should we generate ALL words then simplify, or generate canonical directly?~~
   **Answer**: Generate canonical directly (validated by lead-researcher)

2. ~~For FermionicAlgebra, should we include mixed-degree terms from simplification?~~
   **Answer**: No. Degree-d basis contains ONLY degree-d canonical monomials.

3. ~~Does adding algebra type to VariableRegistry break anything?~~
   **Answer**: No. All construction goes through `create_*_variables` functions.

## Open Questions
1. Should `has_consecutive_repeats` stay as a public utility?
2. How to handle registry with multiple variable types (e.g., mixed Pauli + Fermionic)?

## Testing Strategy
- Unit tests for each `_generate_canonical_words_deg` specialization
- Integration tests comparing with NCTSSOS `get_ncbasis` output
- Verify type consistency with `@ncpolyvar` created variables
- Verify degree purity (no mixed-degree bases)

## References
- Wang & Magron (2021): NCTSSOS.jl
- Wittek (2015): Ncpol2sdpa - canonical word generation
- Helton & McCullough (2004): Noncommutative SOS theory
