# Implementation Plan: Registry-Aware Basis Generation

## Objective
Rewrite `basis.jl` to use `VariableRegistry` for type consistency and correct variable indices.

## STATUS: COMPLETE ✅

## Final API Design

### Primary Functions
```julia
# Returns Vector{Polynomial} - each polynomial is the simplified form of one input monomial
get_ncbasis(registry::VariableRegistry{A,T}, d::Int) -> Vector{Polynomial{A,T,ComplexF64}}
get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) -> Vector{Polynomial{A,T,ComplexF64}}
```

**Key insight**: Return `Vector{Polynomial}` to preserve the 1-to-1 mapping between input words and their simplified forms. For fermionic/bosonic algebras, a single input word can simplify to a multi-term polynomial.

## Core Algorithm: Generate All → Simplify → Wrap in Polynomial

```julia
function get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) where {A,T}
    # Generate ALL words of degree d
    all_words = _generate_all_words(indices(registry), d)

    # Each word becomes one polynomial (the simplified form)
    result = Polynomial{A,T,ComplexF64}[]
    for word in all_words
        mono = Monomial{A}(word)
        simplified = simplify(mono)  # Returns Term or Vector{Term}
        push!(result, Polynomial(simplified))  # Wrap in Polynomial
    end
    return result
end
```

## Algebra-Specific Behavior

| Algebra | simplify(Monomial) | Output Polynomial |
|---------|-------------------|-------------------|
| NonCommutative | Term (coeff=1) | 1-term polynomial |
| Pauli | Term (coeff may be ±1, ±i) | 1-term polynomial |
| Unipotent | Term (U²=I reduces) | 1-term polynomial |
| Projector | Term (P²=P reduces) | 1-term polynomial |
| Fermionic | Vector{Term} (anticommutation) | Multi-term polynomial |
| Bosonic | Vector{Term} (commutation) | Multi-term polynomial |

## Implementation Steps - ALL COMPLETE ✅

### Step 1: Add algebra type to VariableRegistry ✓
- [x] Change `VariableRegistry{T}` to `VariableRegistry{A,T}`
- [x] Update all `create_*_variables` functions
- [x] Add `algebra_type()` and `index_type()` helpers

### Step 2: Rewrite `get_ncbasis` and `get_ncbasis_deg` ✓
- [x] Change signature to take `VariableRegistry{A,T}`
- [x] Return `Vector{Polynomial}` with ComplexF64 coefficients
- [x] Use `indices(registry)` to get variable indices
- [x] Call simplify for each generated word
- [x] Wrap each simplified result in a Polynomial

### Step 3: Implement core basis generation ✓
- [x] `_generate_all_words(indices::Vector{T}, d::Int)` - all words of degree d
- [x] Simplify each monomial and wrap in Polynomial

### Step 4: Remove old API ✓ (fastpoly_test scope)
- [x] Remove old get_ncbasis(::Type{A}, n, d) from basis.jl
- [x] Remove old get_ncbasis_deg(::Type{A}, n, d) from basis.jl
- [x] Remove deprecated helpers

### Step 5: Update tests ✓ (fastpoly_test scope)
- [x] Rewrite test/fastpoly_test/basis.jl to use registry API
- [x] Update test/fastpoly_test/utils.jl
- [x] Update test/fastpoly_test/simplify.jl
- [x] Update test/fastpoly_test/allocations.jl

### Step 6: Update dependent code (DEFERRED - separate task)
- [ ] src/sparse.jl
- [ ] src/gns.jl
- [ ] External tests

## References
- Wang & Magron (2021): NCTSSOS.jl
- Helton & McCullough (2004): Noncommutative SOS theory
