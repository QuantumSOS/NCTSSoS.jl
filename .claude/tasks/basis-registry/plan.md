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
# Returns Vector{Term} - each word simplified to its canonical form with coefficient
get_ncbasis(registry::VariableRegistry{A,T}, d::Int) -> Vector{Term{A,T}}
get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) -> Vector{Term{A,T}}
```

**Key insight**: Return `Vector{Term}` because `simplify(Monomial) → Term`.
The linear relationship between redundant words and canonical forms is preserved in the coefficients.

### Internal Dispatch (algebra-specific)
```julia
_generate_basis_deg(::Type{A}, indices::Vector{T}, d) where {A,T}
```

## Core Algorithm: Generate All → Simplify

```julia
function _generate_basis_deg(::Type{A}, indices::Vector{T}, d::Int) where {A,T}
    d == 0 && return [Term(one(Float64), Monomial{A}(T[]))]
    d < 0 && return Term{A,T}[]

    # Generate ALL words of degree d
    all_words = _generate_all_words(indices, d)

    # Simplify each word → Term (preserves linear relationship)
    terms = Term{A,T}[]
    for word in all_words
        mono = Monomial{A}(word)
        simplified = simplify(mono)  # Returns Term or Vector{Term}
        append!(terms, simplified)
    end

    return terms
end
```

**Why this approach:**
1. The mapping from redundant words to canonical forms encodes SOS constraints
2. `Term` preserves coefficients from simplification (e.g., anticommutation signs)
3. Mixed-degree terms from simplification are included (needed for constraint relations)

## Algebra-Specific Behavior

### NonCommutativeAlgebra
- `simplify` is identity (no relations)
- All n^d words returned as Terms with coeff=1

### PauliAlgebra
- σᵢ² = I (reduces degree)
- σᵢσⱼ = -σⱼσᵢ for i≠j (sign changes)
- Returns canonical forms with coefficients

### UnipotentAlgebra / ProjectorAlgebra
- U² = I / P² = P (reduces degree)
- Returns simplified terms

### FermionicAlgebra
- Anticommutation: aᵢaⱼ† = δᵢⱼ - aⱼ†aᵢ
- Pauli exclusion: aᵢaᵢ = 0
- Returns terms including mixed-degree (important for constraints!)

### BosonicAlgebra
- Commutation: bᵢbⱼ† = δᵢⱼ + bⱼ†bᵢ
- Returns normal-ordered terms with lower-degree contributions

## Implementation Steps

### Step 1: Add algebra type to VariableRegistry ✓ DONE
- [x] Change `VariableRegistry{T}` to `VariableRegistry{A,T}`
- [x] Update all `create_*_variables` functions
- [x] Add `algebra_type()` and `index_type()` helpers

### Step 2: Rewrite `get_ncbasis` and `get_ncbasis_deg`
- [ ] Change signature to take `VariableRegistry{A,T}`
- [ ] Return `Vector{Term{A,T}}` instead of `Vector{Monomial{A,T}}`
- [ ] Use `indices(registry)` to get variable indices
- [ ] Call `_generate_basis_deg(A, indices, d)`

### Step 3: Implement core basis generation
- [ ] `_generate_all_words(indices::Vector{T}, d::Int)` - all words of degree d
- [ ] `_generate_basis_deg(::Type{A}, indices, d)` - generate + simplify

### Step 4: Remove old API
- [ ] Remove `filter_constraint` parameter (automatic via simplify)
- [ ] Remove `T` type parameter (inferred from registry)
- [ ] Remove `n` parameter (from registry)
- [ ] Update `monomials` alias

### Step 5: Update tests
- [ ] Rewrite tests to use registry-based API
- [ ] Unskip the 25 testsets that were skipped due to type mismatch
- [ ] Add tests verifying simplification relationships

### Step 6: Update dependent code
- [ ] Find all usages of old `get_ncbasis` API
- [ ] Update to use registry-based API

## Testing Strategy
- Unit tests for `_generate_all_words`
- Tests verifying simplification preserves linear relationships
- Integration tests comparing with NCTSSOS `get_ncbasis` output
- Type consistency tests with `@ncpolyvar` created variables

## References
- Wang & Magron (2021): NCTSSOS.jl
- Helton & McCullough (2004): Noncommutative SOS theory
