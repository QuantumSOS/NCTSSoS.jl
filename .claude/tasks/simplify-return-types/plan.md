# Implementation Plan: simplify-return-types

## Approach
Implement changes in phases: first add infrastructure (F007-F010), then modify simplify functions (F001-F006), finally run integration tests (F011-F012).

## Steps

### Phase 1: Infrastructure
1. [ ] F007: Add Polynomial constructors for Monomial/Polynomial inputs
   - Files: `src/FastPolynomials/src/types/polynomial.jl`

2. [ ] F008: Update `_add_simplified_terms!` overloads
   - Files: `src/FastPolynomials/src/types/polynomial.jl`

3. [ ] F009: Update `_to_term_vector` overloads
   - Files: `src/FastPolynomials/src/simplification/composed_monomial.jl`

4. [ ] F010: Update `_neat_dot3` and utils
   - Files: `src/FastPolynomials/src/utils.jl`

### Phase 2: Simple Algebras (return Monomial)
5. [ ] F001: NonCommutative simplify → Monomial
   - Files: `src/FastPolynomials/src/simplification/noncommutative.jl`

6. [ ] F002: Projector simplify → Monomial
   - Files: `src/FastPolynomials/src/simplification/projector.jl`

7. [ ] F003: Unipotent simplify → Monomial
   - Files: `src/FastPolynomials/src/simplification/unipotent.jl`

### Phase 3: Complex Algebras (return Polynomial)
8. [ ] F004: Bosonic simplify → Polynomial
   - Files: `src/FastPolynomials/src/simplification/bosonic.jl`

9. [ ] F005: Fermionic simplify → Polynomial
   - Files: `src/FastPolynomials/src/simplification/fermionic.jl`

### Phase 4: Composed Monomials
10. [ ] F006: ComposedMonomial simplify → appropriate type
    - Files: `src/FastPolynomials/src/simplification/composed_monomial.jl`

### Phase 5: Integration
11. [ ] F011: FastPolynomials tests pass
12. [ ] F012: Full test suite passes

## Testing Strategy
- Run `make test-FastPoly` after each phase
- Run `make test` after completing all phases
- Fix regressions immediately before proceeding

## Citations
- Research report: `.claude/tasks/simplify-return-types/context.md`
