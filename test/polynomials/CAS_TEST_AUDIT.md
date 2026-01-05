# CAS test audit: `test/polynomials`

Owner: yushengzhao. Scope: feature coverage + math correctness, self-contained.

## Coverage map (tests → implementation)

- `test/polynomials/algebra_types.jl` → `src/types/algebra.jl` (`AlgebraType`, `coeff_type`, show).
- `test/polynomials/variables.jl` → `src/types/registry.jl` (registry invariants, indexing, symbol formatting, type selection, `create_*_variables`, `subregistry`).
- `test/polynomials/monomials.jl` → `src/types/monomial.jl` (construction, `*`, `^`, `one/zero/isone/iszero`, `hash`, `isless`, `adjoint`).
- `test/polynomials/term.jl` → `src/types/term.jl` (identity/zero semantics, scalar ops, iteration protocol).
- `test/polynomials/polynomial.jl` + `test/polynomials/arithmetic.jl` + `test/polynomials/compare.jl`
  - → `src/types/polynomial.jl` (invariants, ordering, arithmetic, iteration, conversions/promotions).
- `test/polynomials/simplify.jl`
  - → `src/simplification/*.jl` (NonCommutative, Pauli, Unipotent, Projector, Fermionic, Bosonic identities + helpers in `src/util/helpers.jl`).
- `test/polynomials/basis.jl` + `test/polynomials/state_basis.jl`
  - → `src/algorithms/basis.jl` (basis generation, counts, type stability).
- `test/polynomials/state_word.jl` + `test/polynomials/statepolynomial.jl`
  - → `src/states/word.jl`, `src/states/polynomial.jl` (state word canonicalization, ordering, multiplication, `expval`).
- `test/polynomials/canonicalization.jl`
  - → `src/algorithms/canonicalization.jl` (`symmetric_canon`, `cyclic_canon`, `cyclic_symmetric_canon`, `canonicalize`).
- `test/polynomials/utils.jl` → `src/util/helpers.jl` (shared helpers).
- `test/polynomials/allocations.jl` → regression on allocations for hot-path accessors.

## Correctness: key findings

### Bug: Pauli canonicalization accidentally used site-bit decoding

Root cause:
- `src/algorithms/canonicalization.jl` had site-aware `symmetric_canon(::Vector{<:Unsigned})` / `cyclic_symmetric_canon(::Vector{<:Unsigned})` implemented via `decode_site` (bit-packed site encoding).
- `symmetric_canon(::Monomial)` delegated to `symmetric_canon(m.word)`, so `Monomial{PauliAlgebra,UInt*}` would dispatch to the site-bit path.
- Pauli indices are *not* site-bit-encoded (they use `site = (idx-1)÷3+1`), so `decode_site` is not a valid site key for Pauli.

Fix:
- `src/algorithms/canonicalization.jl` now canonicalizes `Monomial` using the algebra type:
  - Pauli: site key `(idx-1)÷3+1` (stable; never reorders within-site).
  - Site-bit algebras: site key `decode_site(idx)`.

Regression checks:
- `test/polynomials/canonicalization.jl` now includes Pauli cases that would fail under the old behavior (intra-site indices must *not* be sorted by numeric value).

### Independent oracle checks (math truth, no external CAS)

- `test/polynomials/matrix_oracles.jl`:
  - Pauli: exact `2^n × 2^n` matrix model; exhaustively checks all 3-letter words on 2 sites against `simplify`.
  - Fermions: exact Jordan–Wigner matrix model; checks nilpotency + CAR cross-mode sign against `simplify`.

## Remaining gaps (feature-level, not line coverage)

- Bosonic CCR currently validated via algebraic identities in tests; no exact finite-dimensional representation exists without truncation artifacts (possible future: truncated Fock oracle + low-occupation subspace checks).
- Site-aware canonicalization is now monomial-only (algebra-aware), avoiding unsigned-word dispatch guessing.
