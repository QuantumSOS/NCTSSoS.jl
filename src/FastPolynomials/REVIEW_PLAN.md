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
- [ ] Monomial hash is precomputed and updated correctly after mutation
- [ ] Monomial ordering (`isless`) is consistent and total
- [ ] Term arithmetic handles coefficient types correctly
- [ ] Polynomial invariants: sorted terms, no duplicates, no zero coefficients
- [ ] Polynomial arithmetic preserves invariants

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
- [ ] Each `simplify` returns the documented type
- [ ] `simplify!` mutates correctly and updates hash
- [ ] `simplify` (non-mutating) preserves original
- [ ] Algebraic rules are mathematically correct
- [ ] Edge cases: empty monomial, single operator, identity

---

## Layer 4: Advanced Monomial Types

| # | File | Purpose | Key Types/Functions |
|---|------|---------|---------------------|
| 12 | `src/composed_monomial.jl` | Tensor products across algebras | `ComposedMonomial{Ts}`, `simplify` → `Vector{Term}` |

**Review Focus:**
- [ ] ComposedMonomial correctly handles mixed algebra types
- [ ] Simplification dispatches to each component's algebra
- [ ] Cartesian product expansion is correct for multi-term algebras

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
| Monomial hash is current | `update_hash!` after mutation |
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
