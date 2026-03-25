# nctssos-expert-gns-edge-cases

## Target Type
feature

## Target
GNS extraction, degenerate optima, and numerical edge cases

## Use Case
Can a specialist push edge cases in GNS extraction — degenerate optima, flatness detection, numerical tolerances — and get correct, precise results? Covers the workflows in `gns_construction_guide.jl`, `gns_optimizer_extraction.jl`, and `pauli_gns_construction.jl`, plus adversarial inputs not shown in examples.

## Expected Outcome
- Calls `gns_reconstruct(monomap, registry, H_deg; hankel_deg, method=:svd, atol, rtol)` and `gns_reconstruct(monomap, sparsity, H_deg; ...)` (sparse variant) and gets mathematically correct operator matrices.
- Verifies flatness: `rank(full_Hankel) == rank(principal_block)` via `gns.singular_values` and `gns.rank` / `gns.full_rank`. When flatness fails, the warning is emitted and `verify_gns` correctly reports the violation.
- Compares `:svd` vs `:cholesky` methods on the same problem; confirms they produce equivalent GNS representations (up to unitary equivalence).
- Constructs a problem with degenerate ground state (rank > 1 moment matrix), verifies `gns_reconstruct` extracts the correct rank and that `robustness_report(gns, hankel)` reports appropriate condition numbers.
- Feeds a moment map missing required entries and gets a structured `ArgumentError` listing the first 5 missing keys (not a cryptic `KeyError`).
- Tests sparse GNS with overlapping cliques and hits the documented `ArgumentError` about pairwise-disjoint requirement; understands this is an unimplemented amalgamation (Theorem 4.2).
- Pushes `atol` / `rtol` to extreme values (1e-15, 1e-3) on a near-rank-deficient Hankel matrix and verifies rank detection degrades gracefully (no silent wrong answers).
- Uses `verify_gns(gns, monomap, registry; poly=H, f_star=opt_value, atol=1e-6)` to check moment reproduction and operator feasibility.

## Agent

### Background
Senior researcher in operator algebras and noncommutative optimization (Dr. Lena). Deep background in NPA hierarchies, SOS/moment duality, and GNS construction theory (Theorem 4.5 in Burgdorf-Klep-Povh). Has published papers on flat extensions and moment problems. Reads `gns.jl`, `gns_cholesky.jl`, and `gns_diagnostics.jl` source code directly. Knows the mathematical guarantees: flat extension ⟹ finite-rank representation, and tests whether the code delivers them.

### Experience Level
Expert

### Decision Tendencies
- Skips examples; goes straight to `gns_reconstruct` method signatures and reads the implementation.
- First test: reproduce Example 4.23 from Burgdorf-Klep-Povh (as in `gns_optimizer_extraction.jl`) and verify operator matrices match the textbook.
- Then constructs adversarial inputs: rank-1 moment matrix with `atol=1e-15` (should detect rank=1), rank-2 matrix with `atol=1e-3` (should still detect rank=2 or warn).
- Checks that `(X + X') / 2` Hermiticity symmetrization (applied for `NonCommutativeAlgebra`, `PauliAlgebra`, `ProjectorAlgebra`, `UnipotentAlgebra`) doesn't mask asymmetry bugs.
- Tests the 3 overloads of `gns_reconstruct` (dense Hankel matrix, dense monomap, sparse monomap) to verify consistent results.
- Will attempt `gns_reconstruct` with `H_deg` too small for flatness and verify the warning fires.

### Quirks
- Distrusts default `atol=1e-8`; always passes explicit tolerances and checks `gns.singular_values` to verify the gap between kept and dropped values.
- Manually constructs Hankel matrices via `hankel_matrix(monomap, get_ncbasis(registry, d))` and compares eigenvalues against `gns.singular_values`.
- When `robustness_report` returns poor condition numbers, digs into whether it's the problem or the solver (runs same problem with Mosek and COSMO, compares).
- Files bug reports with: (a) minimal polynomial, (b) exact analytical moment values, (c) expected GNS rank and operators, (d) proof of why the code's answer is wrong.
- Checks whether the identity monomial is always in the selected basis after rank truncation (documented invariant in `gns.jl`).
