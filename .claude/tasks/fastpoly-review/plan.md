# Implementation Plan: fastpoly-review

## Approach
Systematically review each FastPolynomials source file with its corresponding test file. Use QA gatekeeper to identify coverage gaps, then add missing tests.

## Source Files → Test Files Mapping

| Source File | Test File | Status |
|-------------|-----------|--------|
| `variable_registry.jl` | `variables.jl` | ✅ Done |
| `algebra_types.jl` | `algebra_types.jl` | ✅ Done |
| `monomial.jl` | `monomials.jl` | ✅ Done |
| `term.jl` | `term.jl` | ✅ Done |
| `polynomial.jl` | `polynomial.jl` | ✅ Done |
| `composed_monomial.jl` | `composed_monomial.jl` | ✅ Done |
| `canonicalization.jl` | `canonicalization.jl` | ✅ Done |
| `basis.jl` | `basis.jl` | ✅ Done |
| `simplification/pauli.jl` | `simplify.jl` | ⬜ Pending |
| `simplification/fermionic.jl` | `simplify.jl` | ⬜ Pending |
| `simplification/bosonic.jl` | `simplify.jl` | ⬜ Pending |
| `simplification/projector.jl` | `simplify.jl` | ⬜ Pending |
| `simplification/unipotent.jl` | `simplify.jl` | ⬜ Pending |
| `simplification/noncommutative.jl` | `simplify.jl` | ⬜ Pending |
| `state_types.jl` | `state_word.jl` | ⬜ Pending |
| `state_word.jl` | `state_word.jl` | ⬜ Pending |
| `state_polynomial.jl` | `statepolynomial.jl` | ⬜ Pending |
| `utils.jl` | `utils.jl` | ⬜ Pending |
| `FastPolynomials.jl` | (module file) | ⬜ Pending |

## Additional Test Files
| Test File | Purpose |
|-----------|---------|
| `arithmetic.jl` | Polynomial arithmetic operations |
| `compare.jl` | Comparison operations |
| `allocations.jl` | Memory allocation tests |

## Review Workflow Per File
1. Read source file, understand types and functions
2. Read corresponding test file
3. Run QA gatekeeper to identify coverage gaps
4. Discuss findings with user
5. Add missing tests (prioritize: invariants, error paths, edge cases)
6. Remove redundant code if found
7. Simplify imports
8. Commit changes

## Priority Order (Recommended)
1. **Core types**: monomial.jl, term.jl, polynomial.jl
2. **Algebra**: algebra_types.jl, composed_monomial.jl
3. **Operations**: canonicalization.jl, basis.jl
4. **Simplification**: All 6 simplification files
5. **State**: state_types.jl, state_word.jl, state_polynomial.jl
6. **Utilities**: utils.jl
7. **Module**: FastPolynomials.jl (exports, structure)
