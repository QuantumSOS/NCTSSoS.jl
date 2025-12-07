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
_generate_basis(::Type{PauliAlgebra}, indices::Vector{T}, d) where T
_generate_basis(::Type{FermionicAlgebra}, indices::Vector{T}, d) where T
_generate_basis(::Type{BosonicAlgebra}, indices::Vector{T}, d) where T
_generate_basis(::Type{UnipotentAlgebra}, indices::Vector{T}, d) where T
_generate_basis(::Type{ProjectorAlgebra}, indices::Vector{T}, d) where T
_generate_basis(::Type{NonCommutativeAlgebra}, indices::Vector{T}, d) where T
```

## Algebra-Specific Rules

### Core Algorithm (ALL algebras)
```julia
function _generate_basis(::Type{A}, indices, d)
    all_words = generate_all_words(indices, d)
    canonical_monomials = Set{Monomial}()
    for word in all_words
        terms = simplify(Monomial{A}(word))
        for term in terms
            push!(canonical_monomials, term.monomial)
        end
    end
    return sort(collect(canonical_monomials))
end
```

Key insight: Generate ALL words, then simplify each to get canonical forms.
Collect unique canonical monomials (deduplication via Set).

### PauliAlgebra (UInt16, self-adjoint)
- Indices: unsigned, represent σx, σy, σz at each site
- Simplification: σᵢ² = I (reduces degree), cyclic products
- Result: canonical monomials after Pauli simplification

### UnipotentAlgebra (UInt16, self-adjoint)
- Indices: unsigned with site encoding
- Simplification: U² = I (reduces degree)
- Result: canonical monomials without consecutive repeats

### ProjectorAlgebra (UInt16, self-adjoint)
- Indices: unsigned with site encoding
- Simplification: P² = P (idempotent)
- Result: canonical monomials without consecutive repeats

### FermionicAlgebra (Int32, signed)
- Indices: +i = annihilation aᵢ, -i = creation aᵢ†
- Registry already contains both signs
- Simplification: anticommutation → normal order + lower degree terms
- Pauli exclusion: aᵢaᵢ = 0, aᵢ†aᵢ† = 0
- Result: collect all canonical monomials from simplification

### BosonicAlgebra (Int32, signed)
- Indices: +i = annihilation bᵢ, -i = creation bᵢ†
- Registry already contains both signs
- Simplification: commutation → normal order + lower degree terms
- No exclusion (can have repeated indices)
- Result: collect all canonical monomials from simplification

### NonCommutativeAlgebra (configurable)
- Generic non-commutative variables
- Simplification: none (identity)
- Result: all words as-is

## Implementation Steps

### Step 1: Add algebra type to VariableRegistry
- [ ] Change `VariableRegistry{T}` to `VariableRegistry{A,T}` where A is AlgebraType
- [ ] Update all `create_*_variables` functions to use new signature:
  - `create_pauli_variables` → returns `VariableRegistry{PauliAlgebra,T}`
  - `create_fermionic_variables` → returns `VariableRegistry{FermionicAlgebra,T}`
  - `create_bosonic_variables` → returns `VariableRegistry{BosonicAlgebra,T}`
  - `create_projector_variables` → returns `VariableRegistry{ProjectorAlgebra,T}`
  - `create_unipotent_variables` → returns `VariableRegistry{UnipotentAlgebra,T}`
  - `create_noncommutative_variables` → returns `VariableRegistry{NonCommutativeAlgebra,T}`
- [ ] Update internal helpers: `_create_physical_variables`, `_create_noncommutative_variables`

### Step 2: Rewrite `get_ncbasis` and `get_ncbasis_deg`
- [ ] Change signature to take `VariableRegistry{A,T}`
- [ ] Use `indices(registry)` to get all variable indices
- [ ] Infer both A and T from registry type parameters
- [ ] Call unified `_generate_basis_deg` with algebra dispatch

### Step 3: Implement core basis generation
- [ ] `_generate_all_words(indices::Vector{T}, d::Int)` - generate all words of degree d
- [ ] `_generate_basis_deg(::Type{A}, indices, d)` - core algorithm:
  ```julia
  function _generate_basis_deg(::Type{A}, indices::Vector{T}, d::Int) where {A,T}
      d == 0 && return [Monomial{A}(T[])]  # identity
      d < 0 && return Monomial{A,T}[]

      all_words = _generate_all_words(indices, d)
      canonical_monomials = Set{Monomial{A,T}}()

      for word in all_words
          terms = simplify(Monomial{A}(word))
          for term in terms
              push!(canonical_monomials, term.monomial)
          end
      end

      return sort!(collect(canonical_monomials))
  end
  ```

### Step 4: Handle algebra-specific optimizations (optional)
- [ ] For NonCommutativeAlgebra: skip simplify (identity operation)
- [ ] For algebras with simple rules: consider direct generation if more efficient

### Step 5: Remove old API
- [ ] Remove `filter_constraint` parameter (automatic via simplify)
- [ ] Remove `T` type parameter (inferred from registry)
- [ ] Remove `n` parameter (from registry)
- [ ] Remove `monomials` alias (or update to use registry)
- [ ] Remove `has_consecutive_repeats` if no longer needed

### Step 6: Update tests
- [ ] Rewrite tests to use registry-based API
- [ ] Unskip the 25 testsets that were skipped due to type mismatch
- [ ] Add tests for signed algebra basis generation

### Step 7: Update dependent code
- [ ] Find all usages of old `get_ncbasis` API
- [ ] Update to use registry-based API

## Open Questions
1. Should `has_consecutive_repeats` stay as a public utility?
2. For FermionicAlgebra, should we generate ALL combinations then filter, or generate in normal order directly?
3. How to handle registry with multiple variable types (e.g., mixed Pauli + Fermionic)?

## Testing Strategy
- Unit tests for each `_generate_basis` specialization
- Integration tests comparing with NCTSSOS `get_ncbasis` output
- Verify type consistency with `@ncpolyvar` created variables
