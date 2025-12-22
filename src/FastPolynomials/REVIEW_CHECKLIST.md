# FastPolynomials Review Checklist

Review components in this order (foundational → abstract):

## Phase 1: Core Infrastructure
- [ ] `src/utils.jl` — Utility functions used throughout
- [ ] `src/algebra_types.jl` — Six algebra type definitions (NonCommutative, Pauli, Fermionic, Bosonic, Projector, Unipotent)
- [ ] `src/variable_registry.jl` — Symbol ↔ integer index mapping with algebra-specific encoding

## Phase 2: Monomial Layer
- [ ] `src/monomial.jl` — Core `Monomial{A,T}` type
- [ ] `src/canonicalization.jl` — Monomial ordering and normalization
- [ ] `src/composed_monomial.jl` — Composite monomial structures

## Phase 3: Term & Polynomial Layer
- [ ] `src/term.jl` — `Term{M,C}` = Monomial + Coefficient
- [ ] `src/polynomial.jl` — `Polynomial{A,T,C}` = sorted collection of Terms
- [ ] `src/basis.jl` — Basis construction for moment matrices

## Phase 4: Simplification Rules (one per algebra)
- [ ] `src/simplification/noncommutative.jl` — Base case (no special rules)
- [ ] `src/simplification/pauli.jl` — σᵢ² = I, anticommutation, cyclic phases
- [ ] `src/simplification/fermionic.jl` — {aᵢ, aⱼ†} = δᵢⱼ, normal ordering
- [ ] `src/simplification/bosonic.jl` — [cᵢ, cⱼ†] = δᵢⱼ, commutation corrections
- [ ] `src/simplification/projector.jl` — Pᵢ² = Pᵢ idempotency
- [ ] `src/simplification/unipotent.jl` — U² = I involutory

## Phase 5: State Polynomials (GNS construction support)
- [ ] `src/state_types.jl` — State algebra definitions
- [ ] `src/state_word.jl` — State word representation
- [ ] `src/state_polynomial.jl` — Polynomials over state algebras

## Phase 6: Module Entry Point
- [ ] `src/FastPolynomials.jl` — Module definition, exports, includes

---

**Data Flow Reference:**
```
Variables → Monomials → simplify() → Terms → Polynomial
```

**Notes:**
- 
