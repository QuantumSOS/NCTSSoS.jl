# FastPolynomials Review Plan

Review in this order (follows include dependencies and logical flow).

**Verification Command:** After any modification, run:
```bash
make test-FastPoly
```

---

## Layer 1: Type System Foundation

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 1 | `src/algebra_types.jl` | Algebra type hierarchy + index encoding | `AlgebraType`, `NonCommutativeAlgebra`, `PauliAlgebra`, `FermionicAlgebra`, `BosonicAlgebra`, `ProjectorAlgebra`, `UnipotentAlgebra`, `encode_index`, `decode_site`, `decode_operator_id` |
| 2 | `src/variable_registry.jl` | Symbol ↔ index mapping | `VariableRegistry`, `create_*_variables()`, `symbols`, `indices`, `subregistry` |

**Review Focus:**
- [x] Algebra type hierarchy is complete and well-documented
- [x] Index encoding correctly packs site + operator_id into unsigned integers
- [x] Variable creation functions return correct types for each algebra
- [x] Registry bidirectional mapping is consistent

**Review Status:** ✅ APPROVED (2025-01-22)

**Changes Made:**
- Removed unused `algebra_type` and `index_type` functions
- Fixed `subregistry` type stability with parametric type `VT<:AbstractVector{T}`
- Removed unused `support` and `is_symmetric` from `polynomial.jl`

**Layer 1 Complete** ✅
- `algebra_types.jl`: APPROVED
- `variable_registry.jl`: APPROVED

---

## Layer 2: Core Algebraic Types

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 3 | `src/monomial.jl` | Monomial representation | `Monomial{A,T}`, `degree`, `adjoint`, `*`, `isless`, `isone` |
| 4 | `src/term.jl` | Coefficient + Monomial pair | `Term{M,C}`, arithmetic ops, `iszero`, `isone` |
| 5 | `src/polynomial.jl` | Sorted term collection | `Polynomial{A,T,C}`, `coefficients`, `monomials`, `terms`, `variable_indices`, arithmetic |

**Review Focus:**
- [x] Monomial is immutable (no hash field, hash computed on demand)
- [x] Monomial ordering (`isless`) is consistent and total
- [x] Term arithmetic handles coefficient types correctly
- [x] Polynomial invariants: sorted terms, no duplicates, no zero coefficients
- [x] Polynomial arithmetic preserves invariants

**Review Status:** ✅ APPROVED

**Changes Made (monomial.jl - 2025-01-XX):**
- Made `Monomial` immutable (`struct` instead of `mutable struct`)
- Removed `hash` field - hash now computed on demand via `Base.hash(m.word, h)`
- Removed `update_hash!` function (no longer needed)
- Removed `star`, `star!`, `adjoint!` - consolidated to just `Base.adjoint`
- Simplified `Base.:(==)` - direct word comparison without hash caching
- Updated simplification functions to return new monomials instead of mutating
- Added physics notation notes to `adjoint` docstring (dagger †, star *)

**Changes Made (term.jl - 2025-01-22):**
- Made `Term` immutable (`struct` instead of `mutable struct`) for performance
  - 0 heap allocations for arithmetic (was 32 bytes per operation)
  - ~5x faster term creation due to stack allocation/inlining
  - Monomial reference is shared (no deep copy ever occurs)
- Added `Base.hash` to fix equality/hash contract violation
- Generalized `one`/`zero` to work with any `AbstractMonomial` (including `ComposedMonomial`)
- Added `one`/`isone` for `ComposedMonomial` in composed_monomial.jl
- Fixed complex coefficient display with parentheses: `(1.0 + 2.0im) * [1, 2]`
- Kept `iterate` for destructuring support used in sos_solver.jl

**Changes Made (polynomial.jl - 2025-01-22):**
- Added `default_coeff_type` trait in algebra_types.jl (ComplexF64 for Pauli, Float64 for others)
- Updated `Polynomial(m::Monomial)` to use `default_coeff_type(A)` instead of hardcoded Float64
- Added `simplify(p::Polynomial)` function that normalizes each monomial using algebra-specific rules
- Refactored polynomial multiplication to use raw concatenation + simplify (decoupled arithmetic from simplification)
- Removed old `_add_simplified_terms!` functions (replaced by `_collect_simplified_terms!` in simplify)
- Changed `degree(zero_polynomial)` from -1 to -Inf to preserve algebraic identity `deg(p*q) = deg(p) + deg(q)`
- Removed `maxdegree` alias - all usages converted to `degree`

**Layer 2 Complete** ✅
- `monomial.jl`: APPROVED
- `term.jl`: APPROVED
- `polynomial.jl`: APPROVED

---

## Layer 3: Algebra-Specific Simplification

| # | File | Algebra | Rule | Returns |
|---|------|---------|------|---------|
| 6 | `src/simplification/projector.jl` | `ProjectorAlgebra` | P² = P | `Monomial` |
| 7 | `src/simplification/unipotent.jl` | `UnipotentAlgebra` | U² = I | `Monomial` |
| 8 | `src/simplification/noncommutative.jl` | `NonCommutativeAlgebra` | site commutation only | `Monomial` |
| 9 | `src/simplification/pauli.jl` | `PauliAlgebra` | σᵢ² = I, XY=iZ | `Term` |
| 10 | `src/simplification/fermionic.jl` | `FermionicAlgebra` | {aᵢ, aⱼ†} = δᵢⱼ | `Polynomial` |
| 11 | `src/simplification/bosonic.jl` | `BosonicAlgebra` | [cᵢ, cⱼ†] = δᵢⱼ | `Polynomial` |

**Review Focus:**
- [x] Each `simplify` returns the documented type
- [x] `simplify` returns new monomial (immutable design)
- [x] `simplify!` only available when return type matches input (NC, Projector, Unipotent)
- [x] No `simplify!` for Pauli (returns Term), Fermionic/Bosonic (returns Polynomial)
- [x] Algebraic rules are mathematically correct
- [x] Edge cases: empty monomial, single operator, identity

**Review Status:** ✅ APPROVED (2025-01-23)

**Changes Made (simplification - 2025-01-22/23):**

*noncommutative.jl, projector.jl, unipotent.jl:*
- Added docstrings for `simplify!` functions
- Reordered functions to standard pattern: internal → `simplify!` → `simplify`
- Optimized `unipotent.jl`: use single `deleteat!(word, i:i+1)` instead of two calls

*pauli.jl:*
- Removed misleading `simplify!` (returns `Term`, not `Monomial`)
- Kept only `simplify` with clear docstring explaining return type

*fermionic.jl:*
- Removed `simplify!` (returns `Polynomial`)
- Removed `_fermi_mode` wrapper, use `_operator_mode` directly
- Use shared `normal_order_key` for sorting
- Use shared `combine_like_terms` (was `_combine_like_terms_fermi`)
- **Fixed `iszero` bug**: Changed from sequential duplicate check to net flux algorithm
  - Old: checked consecutive same-type operators (missed `a₁ a₂ a₁ = 0`)
  - New: `|annihilations - creations| >= 2` per mode → nilpotent

*bosonic.jl:*
- Removed `simplify!` (returns `Polynomial`)
- Removed `bosonic_mode` and `bosonic_sort_key` wrappers
- Use `_operator_mode` and shared helpers from utils.jl

*utils.jl - New shared helpers:*
- `normal_order_key(op)` - sort key for normal ordering (creation < annihilation, then by mode)
- `find_first_out_of_order(word)` - find first position out of normal order
- `is_normal_ordered(word)` - check if word is in normal order (uses `find_first_out_of_order`)
- `combine_like_terms(terms)` - generic term combining for any algebra

**Simplify API Pattern:**

| Algebra | `simplify` | `simplify!` | Returns | Reason |
|---------|------------|-------------|---------|--------|
| NonCommutative | ✅ | ✅ | Monomial | No type change |
| Projector | ✅ | ✅ | Monomial | No type change |
| Unipotent | ✅ | ✅ | Monomial | No type change |
| Pauli | ✅ | ❌ | Term | Complex phase |
| Fermionic | ✅ | ❌ | Polynomial | Wick contractions |
| Bosonic | ✅ | ❌ | Polynomial | Rook numbers |

**Layer 3 Complete** ✅

---

## Layer 4: Advanced Monomial Types

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 12 | `src/composed_monomial.jl` | Tensor products across algebras | `ComposedMonomial{Ts}`, `simplify` → `Vector{Term}` |

**Review Focus:**
- [x] ComposedMonomial correctly handles mixed algebra types
- [x] Simplification dispatches to each component's algebra
- [x] Cartesian product expansion is correct for multi-term algebras

**Review Status:** ✅ APPROVED (2025-01-23)

**Changes Made (composed_monomial.jl - 2025-01-23):**
- Removed precomputed `hash` field - hash now computed on demand
- Removed `_compute_composed_hash` helper function
- Added unified iteration protocol for `Monomial`, `Term`, `Polynomial` yielding `(coefficient, monomial)` pairs
- Added `coeff_type` trait for compile-time coefficient type inference
- Removed dead legacy code:
  - `_infer_coef_type(component_terms::Tuple)` - replaced by `_infer_coef_type_from_types`
  - `_to_term_vector(::Monomial)`, `_to_term_vector(::Term)`, `_to_term_vector(::Vector{<:Term})`, `_to_term_vector(::Polynomial)` - iteration protocol eliminates need for conversion
- Simplified `_expand_simplified_components` to use iteration protocol directly via `_cartesian_product_iter!`
- Updated `simplify` to always return `Vector{Term}` for consistent API

**Layer 4 Complete** ✅

---

## Layer 5: Canonicalization & Basis

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 13 | `src/canonicalization.jl` | Ordering & normalization | `symmetric_canon`, `cyclic_canon`, `cyclic_symmetric_canon`, `canonicalize` |
| 14 | `src/basis.jl` | Basis generation for SDP | `get_ncbasis`, `get_ncbasis_deg` |

**Review Focus:**
- [ ] `symmetric_canon`: min(word, reverse(word))
- [ ] `cyclic_canon`: lexicographically smallest rotation
- [ ] `cyclic_symmetric_canon`: combines both for trace states
- [ ] Basis generation produces correct degree-bounded monomials
- [ ] No duplicate basis elements

---

## Layer 6: State Polynomial System (GNS)

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 15 | `src/state_types.jl` | State type hierarchy | `StateType`, `Arbitrary`, `MaxEntangled` |
| 16 | `src/state_word.jl` | State expectation products | `StateWord`, `NCStateWord`, `simplify(NCStateWord)`, `expval`, `ς`, `tr`, `get_state_basis` |
| 17 | `src/state_polynomial.jl` | Polynomials over state words | `StatePolynomial`, `NCStatePolynomial`, `expval`, arithmetic |

**Review Focus:**
- [ ] StateWord invariants: sorted, involution-canonicalized
- [ ] MaxEntangled uses cyclic-symmetric canonicalization
- [ ] Arbitrary uses symmetric (involution) canonicalization
- [ ] NCStateWord correctly separates state part from operator part
- [ ] `expval` collapses NCStateWord to StateWord correctly
- [ ] State polynomial arithmetic preserves invariants

---

## Layer 7: Utilities & Module Entry

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 18 | `src/utils.jl` | Helper functions | `sorted_union`, `sorted_unique`, `neat_dot`, `_neat_dot3` |
| 19 | `src/FastPolynomials.jl` | Module definition | `module FastPolynomials`, includes, exports |

**Review Focus:**
- [ ] All public functions are exported
- [ ] Include order matches dependency order
- [ ] No missing exports for documented API
- [ ] Utility functions are correctly implemented

---

## Key Invariants Checklist

| Invariant | Where Enforced |
|-----------|----------------|
| Polynomial terms are sorted | `Polynomial` constructor |
| Polynomial has no duplicate monomials | `Polynomial` constructor (combines like terms) |
| Polynomial has no zero coefficients | `Polynomial` constructor (filters zeros) |
| Monomial is immutable | `struct Monomial` (not `mutable struct`) |
| StateWord monomials are sorted | `StateWord` constructor |
| StateWord monomials are canonicalized | `_state_canon` in constructor |
| Index encoding: site in lower bits | `encode_index` / `decode_site` |

---

## Simplification Return Type Reference

| Algebra | `simplify` Returns | Reason |
|---------|-------------------|--------|
| `ProjectorAlgebra` | `Monomial` | No coefficient change (P²=P) |
| `UnipotentAlgebra` | `Monomial` | No coefficient change (U²=I removes pairs) |
| `NonCommutativeAlgebra` | `Monomial` | No coefficient change (reordering only) |
| `PauliAlgebra` | `Term` | Complex phase from cyclic products (XY=iZ) |
| `FermionicAlgebra` | `Polynomial` | Multi-term from Wick contractions |
| `BosonicAlgebra` | `Polynomial` | Multi-term from commutator corrections |

---

## Quick Reference

```bash
# Run FastPolynomials tests only
make test-FastPoly

# Run all tests (slower)
make test

# Run specific test file
julia --project -e 'include("test/fastpoly_test/monomials.jl")'
```
