# Task: simplify-return-types

## Request
Change `simplify!` and `simplify` return types:
1. `Monomial` for NonCommutative, Projector, Unipotent algebra
2. `Term` for Pauli algebra (no change)
3. `Polynomial` for Bosonic and Fermionic algebra
4. `ComposedMonomial` follows most complex algebra's return type

## Project State
- Branch: `fix/simplify-sort-behavior`
- Current implementation returns `Term` for simple algebras, `Vector{Term}` for Bosonic/Fermionic
- Key files:
  - `src/FastPolynomials/src/simplification/noncommutative.jl`
  - `src/FastPolynomials/src/simplification/pauli.jl`
  - `src/FastPolynomials/src/simplification/projector.jl`
  - `src/FastPolynomials/src/simplification/unipotent.jl`
  - `src/FastPolynomials/src/simplification/bosonic.jl`
  - `src/FastPolynomials/src/simplification/fermionic.jl`
  - `src/FastPolynomials/src/simplification/composed_monomial.jl`

## Decisions
- Proceed despite type instability concerns in generic code
- Will add conversion infrastructure to handle varying return types
- Will update all call sites systematically

## Progress
- [x] Research complete (lead-researcher)
- [ ] Implementation
- [ ] QA Review
- [ ] Documentation

## Handoff Summary
- Research identified MEDIUM-HIGH risk but user decided to proceed
- Key call sites: `_neat_dot3`, `_add_simplified_terms!`, `_to_term_vector`, basis generation
- Need to add `Polynomial(::Monomial)` and `Polynomial(::Polynomial)` constructors
- Bosonic/Fermionic: wrap `Vector{Term}` in `Polynomial` at return
