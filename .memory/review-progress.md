# Code Review Progress

## src/types/

- [x] `algebra.jl` — AlgebraType hierarchy, site encoding
- [ ] `monomial.jl`
  - [x] Lines 1-540: AbstractMonomial, NormalMonomial, helpers
  - [ ] Lines 541-1521: Monomial struct, constructors, multiplication
- [ ] `polynomial.jl`
- [x] `registry.jl`
- [x] `arithmetic.jl` — implemented _coeff_to_number, addressed all FIXMEs

## src/simplification/

- [x] `noncommutative.jl` — site-based simplification
- [x] `site_helpers.jl` — shared validation/sorting helpers
- [x] `projector.jl` — idempotency P²=P
- [x] `unipotent.jl` — involution U²=I
- [x] `pauli.jl` — cyclic products, two-pointer in-place simplify
- [x] `fermionic.jl` — Wick's theorem normal ordering
- [x] `bosonic.jl` — rook number algorithm for normal ordering

## src/states/

- [x] `types.jl`
- [ ] `polynomial.jl`
- [ ] `word.jl`

## src/optimization/

<!--- [ ] `gns.jl` skipped-->
- [ ] `sos.jl`
- [ ] `moment.jl`
- [ ] `sparsity.jl`
- [ ] `interface.jl`
- [x] `elimination.jl`
- [x] `problem.jl`

## src/algorithms/

- [x] `canonicalization.jl`
- [x] `basis.jl`

## src/util/

- [x] `helpers.jl` — reviewed, bugs fixed, docs updated
  - [ ] `_neat_dot3(NCStateWord, NCStateWord, NCStateWord)` — needs review after states/ review

## test/polynomials/

- [x] `algebra_types.jl` — DRY refactor with data-driven tests, verified expectations
