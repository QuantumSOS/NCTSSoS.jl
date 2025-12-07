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

### Step 2: Rewrite `get_ncbasis` and `get_ncbasis_deg` ✓ DONE
- [x] Change signature to take `VariableRegistry{A,T}`
- [x] Return `Vector{Term}` with ComplexF64 coefficients
- [x] Use `indices(registry)` to get variable indices
- [x] Call simplify for each generated word

### Step 3: Implement core basis generation ✓ DONE
- [x] `_generate_all_words(indices::Vector{T}, d::Int)` - all words of degree d
- [x] Simplify each monomial and collect terms

### Step 4: Remove old API ✓ DONE (fastpoly_test scope)
- [x] Remove old get_ncbasis(::Type{A}, n, d) from basis.jl
- [x] Remove old get_ncbasis_deg(::Type{A}, n, d) from basis.jl
- [x] Remove _generate_words(n, d, T) helper
- [x] Remove monomials alias
- [x] Remove has_consecutive_repeats
- [x] Remove filter_constraint parameter
- [ ] Update utils.jl get_basis wrappers (DEFERRED)
- [ ] Update dependent code (sparse.jl, gns.jl, etc.) (DEFERRED)

### Step 5: Update tests ✓ DONE (fastpoly_test scope)
- [x] Rewrite test/fastpoly_test/basis.jl to use registry API
- [x] Update test/fastpoly_test/utils.jl
- [x] Update test/fastpoly_test/simplify.jl
- [x] Update test/fastpoly_test/allocations.jl
- [ ] Update other test files using old API (DEFERRED - test/sparse.jl etc.)

### Step 6: Update dependent code (DEFERRED)
- [ ] src/sparse.jl
- [ ] src/gns.jl
- [ ] test/sos_solver.jl, test/moment_solver.jl, etc.
- [ ] docs/examples

**Note:** Steps 4-6 deferred for external code per user request. FastPolynomials tests (1017) all passing.

## Testing Strategy
- Unit tests for `_generate_all_words`
- Tests verifying simplification preserves linear relationships
- Integration tests comparing with NCTSSOS `get_ncbasis` output
- Type consistency tests with `@ncpolyvar` created variables

## References
- Wang & Magron (2021): NCTSSOS.jl
- Helton & McCullough (2004): Noncommutative SOS theory
