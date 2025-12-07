# Task: fastpoly-review

## Request
Review all source files in `src/FastPolynomials/src/` and their corresponding test files. For each file:
1. Understand the code structure and design
2. Run QA analysis to identify missing test coverage
3. Add missing tests (especially error paths, edge cases, invariant validation)
4. Remove redundant/dead code if found
5. Simplify test imports

## Project State
- Branch: `fastpolynomial-redesign`
- FastPolynomials is a custom polynomial implementation for NCTSSoS
- Uses algebra-aware type system with multiple dispatch

## Progress
- [x] variable_registry.jl - Reviewed, tests added (61→114), removed `monomial_ring`
- [ ] algebra_types.jl
- [x] monomial.jl - Reviewed, tests added (29→66), removed buggy cross-type isless
- [ ] term.jl
- [ ] polynomial.jl
- [ ] composed_monomial.jl
- [ ] canonicalization.jl
- [ ] basis.jl
- [ ] simplification/*.jl (6 files)
- [ ] state_types.jl
- [ ] state_word.jl
- [ ] state_polynomial.jl
- [ ] utils.jl
- [ ] FastPolynomials.jl (main module)

## Decisions
- Remove redundant APIs (like `monomial_ring`)
- Only import internal `_` prefixed functions explicitly in tests
- Extract tests into separate files when they belong to different source files

## Handoff Summary
- **Completed**: variable_registry.jl, monomial.jl reviews
- **Key Finding**: basis.jl uses Int by default but @ncpolyvar uses UInt64 - needs rewrite
- **Blocker**: 25 testsets skipped due to type mismatch (see TODO in basis.jl)
- **Next Step**: Continue with term.jl or algebra_types.jl
