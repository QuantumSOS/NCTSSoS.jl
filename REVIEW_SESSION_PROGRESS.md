# CAS Review Session Progress

**Last Updated:** 2025-12-31
**Branch:** `refactor/state-symbol` (Section 8 review complete)

## Completed Sections

### Section 1: AlgebraType (Abstract Singleton) - COMPLETE
- All checklist items verified with tests in `test/polynomials/algebra_types.jl`
- Tests cover: `default_coeff_type`, `show`, all 6 concrete subtypes

### Section 2: VariableRegistry{A, T} - COMPLETE
- All checklist items verified
- Added 20 new tests for `subregistry` function in `test/polynomials/variables.jl`

### Section 3: Monomial{A, T} - COMPLETE
- All checklist items verified
- Added 13 new tests for `variable_indices` and `expval` in `test/polynomials/monomials.jl`

### Section 4: Term{M, C} - COMPLETE
- All 15 checklist items verified
- Added 1 new testset: `coeff_type` (11 tests) in `test/polynomials/term.jl`
- Existing test coverage was comprehensive; only `coeff_type` was missing tests
- All 1489 polynomial tests pass

### Section 5: Polynomial{A, T, C} - COMPLETE
- All 52 checklist items verified
- Added 3 new testsets (37 tests total) in `test/polynomials/polynomial.jl`:
  1. `isless(Polynomial, Polynomial) - Graded Lex Ordering` (10 tests)
  2. `convert(Polynomial{A,T,C2}, p) - Coefficient Type Conversion` (16 tests)
  3. `coeff_type - Coefficient Type Accessor` (11 tests)
- All 1526 polynomial tests pass

### Section 6: Simplification (per Algebra) - COMPLETE
- All 12 checklist items verified
- Added comprehensive NCTSSOS oracle tests in `test/polynomials/simplify.jl` (163 new tests total):
  1. Word-level canonicalization (symmetric, cyclic, cyclic-symmetric)
  2. Word-level constraint_reduce (unipotent, projector)
  3. Basis generation counts (n=1-4, d=0-3)
  4. Site-aware simplification (two-site UnipotentAlgebra, ProjectorAlgebra)
  5. PauliAlgebra Algebraic Identities
  6. FermionicAlgebra CAR Identities
  7. BosonicAlgebra CCR Identities
  8. Multi-site Commutation
- All 1675 polynomial tests pass

### Section 7: Index Encoding Utilities - COMPLETE
- All checklist items verified with thorough tests in `test/polynomials/algebra_types.jl`
- Tests cover: `site_bits`, `max_sites`, `max_operators`, `encode_index`, `decode_site`, `decode_operator_id`, `select_uint_type`

### Section 8: Algorithms - COMPLETE
- **Basis Generation** (`src/algorithms/basis.jl`):
  - `_generate_all_words` - All words of exact length d
  - `get_ncbasis_deg` - Degree-specific basis with simplification
  - `get_ncbasis` - Concatenate degrees 0:d
  - NCTSSOS oracle verification in `test/polynomials/simplify.jl`

- **Canonicalization** (`src/algorithms/canonicalization.jl`):
  - `symmetric_canon(word)` - min(word, reverse(word))
  - `cyclic_canon(word)` - Lexicographically minimal rotation
  - `cyclic_symmetric_canon(word)` - min(cyclic(word), cyclic(reverse(word)))
  - `canonicalize(m; cyclic=false|true)` - Monomial interface
  - `canonicalize(p; cyclic=false|true)` - Polynomial interface with term merging
  - NCTSSOS oracle verification with hardcoded expected values
  - Site-aware specializations documented (Unsigned type dispatch)

## Major Refactor Completed

### StateSymbol Refactor
- **Branch:** `refactor/state-symbol`
- **Commit:** `7c8fa9b`
- **Changes:**
  - Added `StateSymbol{ST,A,T}` as atomic state expectation type with auto-canonicalization
  - Refactored `StateWord` to store `state_syms::Vector{StateSymbol}` instead of `state_monos::Vector{Monomial}`
  - Exported `StateSymbol` from NCTSSoS
  - Updated all references from `state_monos` to `state_syms`
  - Added regression tests for StateSymbol canonicalization

## NCTSSOS Oracle Reference

**Local repo:** `/Users/yushengzhao/projects/NCTSSOS`

Key oracle functions:
- `NCTSSOS.get_ncbasis(n, d; ind=..., binary=false)` - basis generation
- `NCTSSOS._sym_canon(word)` - symmetric canonicalization
- `NCTSSOS._cyclic_canon(word)` - cyclic canonicalization
- `NCTSSOS.star(m)` - adjoint operation
- `NCTSSOS.constraint_reduce!` - unipotent/projector simplification

## Test Commands

```bash
# Polynomial tests only
julia --project -e 'using Pkg; Pkg.test(test_args=["--polynomials"])'

# Solver tests
julia --project -e 'using Pkg; Pkg.test(test_args=["--solvers"])'

# Full CI suite
julia --project -e 'using Pkg; Pkg.test()'

# Full local suite with Mosek
julia --project -e 'using Pkg; Pkg.test(test_args=["--local"])'
```

## Files Changed in This Session

1. `test/polynomials/variables.jl` - Added `subregistry` tests
2. `test/polynomials/monomials.jl` - Added `variable_indices`, `expval` tests
3. `src/states/word.jl` - Added `StateSymbol`, refactored `StateWord`
4. `src/NCTSSoS.jl` - Added `StateSymbol` export
5. `src/optimization/sparsity.jl` - Updated `state_monos` -> `state_syms`
6. `src/states/polynomial.jl` - Updated comment
7. `test/polynomials/state_word.jl` - Updated field refs, added StateSymbol tests
8. `test/polynomials/statepolynomial.jl` - Updated field ref
9. `CAS_REVIEW_CHECKLIST.md` - Updated checked items (Sections 1-8 complete)
10. `test/polynomials/term.jl` - Added `coeff_type` testset (11 tests)
11. `test/polynomials/polynomial.jl` - Added 3 testsets (37 tests)
12. `test/polynomials/simplify.jl` - Added NCTSSOS oracle tests (163 tests)
13. `test/polynomials/runtests.jl` - Removed stale include for nctssos_oracle.jl

## Summary

**All 8 sections of the CAS Review Checklist are now COMPLETE.**

The polynomial core (Sections 1-8) has been fully reviewed with:
- 1675 polynomial tests passing
- NCTSSOS oracle comparisons for word-level operations
- Comprehensive coverage of all 6 algebra types
- Site-aware specializations documented

The review verified that NCTSSoS produces structurally correct results matching NCTSSOS for:
- Variable creation and monomial operations
- Term and polynomial arithmetic
- Algebra-specific simplification rules
- Basis generation
- Canonicalization algorithms
