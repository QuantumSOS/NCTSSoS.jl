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
- [ ] monomial.jl
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
- **Completed**: variable_registry.jl full review
- **Commit**: `a027fba` - refactor(fastpoly): improve variable_registry tests and remove redundant API
- **Next Step**: Start with monomial.jl
