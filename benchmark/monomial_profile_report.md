# Monomial Implementation Profile — Findings

Date: 2026-06-12 · Julia 1.12.3 · branch `perf/monomial-allocations`
Harness: `benchmark/monomial_profile.jl` (re-runnable, no extra deps)

## Verdict

The **per-monomial primitives are well-tuned** (zero-alloc hash/compare, in-place
`simplify!`, trusted constructors). The **composition of them in hot loops is not**:
every product allocates 1–3 short-lived `Vector`s, producing 15 GiB of churn and
33–37% GC time for a single moment-relaxation build (Heisenberg n=10, order=3).
And for Pauli there is a structural ceiling: a heap `Vector{UInt8}` word is the
wrong representation when a bit-packed symplectic encoding gives O(1) zero-alloc
products.

## Micro-benchmarks (per op, typed-closure harness)

| operation | ns/op | B/op | notes |
|---|---|---|---|
| pauli `hash(m)` | ~0 | 0 | `_hash_word_vector` claim verified |
| pauli `==` / `isless` | ~2 | 0 | |
| pauli Dict lookup | 3.9 | 0 | |
| pauli `simplify!` (reused buffer, deg-9 raw) | 25.8 | 0 | in-place path is clean |
| pauli `_mul_term` (deg3 × deg3) | 26.8 | 64 | 1 alloc: combined word |
| pauli `simplify!(_neat_dot3(...))` moment kernel | 56.1 | 64 | 1 alloc per matrix element term |
| pauli full `*` → Polynomial | 35.7 | 144 | word + terms vec + wrapper |
| pauli `adjoint(m)` | 11.5 | 64 | |
| `_process_terms!` (90 terms) | 2510 | 2624 | sort scratch; acceptable |
| pauli `H*H` (27×27 terms) | 43757 | 82624 | ≈113 B/product, consistent |
| fermionic `*` (deg2 × deg2, PBW) | 1094 | 4768 | `Vector{Tuple{Int,Vector}}` expansion churn |

## End-to-end attribution (compute_sparsity + moment_relax, Heisenberg n=10, order=3)

Baseline: **19.6 s wall · 15.11 GiB alloc · 33% GC**.

Sampled allocation sites (share of NCTSSoS-src allocs):

1. `helpers.jl:157 _neat_dot3` — **~38%**. Fresh result `Vector` per product. Also
   dominates the *time* profile (the `similar` + reversed-view copy lines).
2. `moment.jl:147 _moment_key` — **~17%**. `copy(mono.word)` per moment-dict key.
3. `pauli.jl:214 simplify` — **~12%**. `copy(word)` on vectors that were already
   freshly allocated by `_neat_dot3` (pure waste). *(fixed, see below)*
4. `arithmetic.jl:82 _simplified_to_terms` — **~12%**. Allocates a 1-element
   `Vector` for every Monoid/TwistedGroup product result (always ≤1 term).
5. `sparsity.jl gsupp` — **~8%**. `last.(g.terms)` broadcast allocated inside a
   triple loop; `last.(terms)` again. *(fixed, see below)*

Also found: `get_term_sparsity_graph` sorted `activated_supp` but then used
**linear `in`** on it — O(n) membership per product instead of O(log n) `insorted`.

## Applied fixes

Round 1 — `src/optimization/sparsity.jl`:
- `gsupp`: hoist support monomials out of the triple loop; `simplify!` on owned
  vectors; push monomials directly instead of `last.()` broadcasts.
- `get_term_sparsity_graph`: `simplify!` on owned vectors; `insorted` binary
  search instead of linear `in` on the sorted support (O(n) → O(log n)).

Round 2:
- `src/types/arithmetic.jl`: new `_simplified_monomials` — tuple-returning
  (zero-alloc) sibling of `_simplified_to_terms` for Monoid/TwistedGroup, used
  by the three sparsity hot loops (`init_activated_supp`,
  `get_term_sparsity_graph`, `gsupp`).
- `src/optimization/moment.jl` `_moment_key`: generic path now calls
  `symmetric_canon(mono)` directly (drops a `StateSymbol` construction + an
  extra copy, value-identical); TwistedGroup/PBW path returns `mono.word`
  without copying (words are immutable by contract; keys must not be mutated).

Round 3 — scratch-buffer assembly kernel (kills the `_neat_dot3` share):
- `src/util/helpers.jl`: new `_neat_dot3!(buf, a, m, b)` / `neat_dot!(buf, a, b)`
  — write the raw adjoint(a)·m·b word into a reusable scratch buffer (resized
  to fit) instead of allocating a fresh `Vector` per product. `simplify!` then
  runs in place on the buffer (Pauli's in-place resize-on-shrink included).
- Aliasing contract: for Monoid/TwistedGroup algebras the simplified result IS
  the buffer, so words are **copied only at insertion time** (and after dedup
  probes, only on *first* insertion); PBW `simplify!` always returns fresh
  words, so no copy. Helpers: `_copy_if_buffer_aliased` (arithmetic.jl),
  `_push_scaled_buffered_terms!` / `_remember_buffered_monomials!` (moment.jl,
  the latter probes the Dict with a transient buffer-aliased wrapper and copies
  only on miss).
- Converted call sites: `_build_constraint_matrix`, `_moment_matrix_basis`,
  `_polynomial_total_basis` (both loops), `_append_moment_eq_constraints!`
  (moment.jl); `init_activated_supp`, `get_term_sparsity_graph` (membership
  checks now allocate nothing — transient probe + `insorted`), and
  `term_sparsity_graph_supp`/`gsupp` (sparsity.jl). `gsupp` now accumulates
  first occurrences directly (same iteration order and dedup semantics as the
  old `union(supp_sets...)`, bit-identical output, no per-pair vectors).

### Measured cumulative effect (Heisenberg n=10, order=3 build)

| | wall | alloc | GC |
|---|---|---|---|
| baseline | 19.6 s | 15.11 GiB | 33% |
| + round 1 | 18.8 s | 12.46 GiB | 37% |
| + round 2 | 16.7 s | 9.34 GiB | 33% |
| + round 3 | **11.2 s** | **5.77 GiB** | 25% |

Cumulative: **−43% wall, −62% allocations** vs baseline. After round 3 the
allocation profile is dominated by the unavoidable insertion-time word copies
(`_push_scaled_buffered_terms!`) and `_linearize_moment_polynomial` /
`_build_constraint_matrix` term-vector growth — i.e. data that actually
escapes. The `_neat_dot3` per-product churn is gone.

### Scale check: 2D Heisenberg 5×5 (25 sites, 40 bonds), order=2 build

Harness: `benchmark/heisenberg2d_build_bench.jl` (timing + structural SHA-256
fingerprint). Single 75-variable clique, total basis 1,089,526 monomials.
A/B vs `main` (merge-base `0136783`), 5 timed runs each after warmup:

| | wall (median) | alloc | speedup |
|---|---|---|---|
| `main` (0136783) | 38.98 s | 28.26 GiB | — |
| this branch | **10.53 s** | **4.53 GiB** | **3.7× wall, 6.2× alloc** |

Equivalence: SHA-256 fingerprint over cliques, term-sparsity supports, block
bases, total basis, `n_unique`, and every constraint-matrix polynomial is
**bit-identical** on both sides:
`42d2e647b6ed12d65b68ddf853e2fbca36fcf9d2e35867d4fe09f653de4975ea`.

### Verification

- SHA-256 digests of sparsity structures **and** full `MomentLinearData`
  (moments, adjoint keys, pivots, key_to_monomial, objective/constraint forms,
  total basis) are bit-identical before/after on Pauli, NC, and fermionic
  problems (re-verified after each round, including round 3).
- Full `make test`: 4458 passed, 2 broken (pre-existing), 1 failed —
  `E9 n=4 dense vs sparse` COSMO tolerance check (gap 2.6e-3 vs atol 1e-3).
  Confirmed pre-existing flake, not a regression: in a clean A/B env the failing
  comparison produces **bit-identical objectives** with and without the patch;
  the test file itself documents ~6e-3 drift across dependency resolutions.
  Same single flake after round 3.

## Remaining headroom, in order of return-on-effort

1. **Structural (the real ceiling, Pauli only):** symplectic bit-packed Pauli
   monomials — x-mask/z-mask in `UInt64`/`UInt128`. Product = XOR + popcount
   phase: O(1), zero-alloc, isbits (terms become isbits tuples → no GC tracking,
   integer hash/compare, much faster sorts and Dicts). This is the Stim /
   PauliStrings.jl representation; expect order-of-magnitude gains in symbolic
   assembly for spin problems. Fits behind `AbstractMonomial{PauliAlgebra}` as a
   parallel monomial type without disturbing the generic word path.
2. **Minor:** the GNS extraction path (`gns.jl`) still uses the
   `symmetric_canon(expval(...))` double-copy and `simplify(A, fresh_vector)`
   patterns; cold post-solve code, fix opportunistically.
