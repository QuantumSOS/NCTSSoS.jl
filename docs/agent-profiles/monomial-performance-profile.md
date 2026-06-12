# Monomial Performance Profile & Improved Design

Date: 2026-06-12 · Julia 1.12.3 · Apple Silicon · Pauli system, N = 10 sites,
Heisenberg-style Hamiltonian (27 terms), degree-2 moment basis (436 monomials).

## 1. Measured baseline

| Operation | time/op | alloc/op |
|---|---|---|
| `NormalMonomial * NormalMonomial` | 0.11 µs | 432 B (~10 allocations) |
| moment entry (`_neat_dot3` + `simplify` + terms) | 0.17 µs | 240 B |
| `NormalMonomial * Polynomial` (27 terms) | 12.6 µs | 76.9 KB |
| `Polynomial * Polynomial` (27×27) | 62.6 µs | 289.6 KB |
| moment-like loop, 80×80 basis × 27 terms | 16.9 ms | **35.9 MB** |
| `Dict` lookup sweep over 436-monomial basis | 27 µs | 27.4 KB (≈63 B/lookup) |
| `hash(::NormalMonomial)` | 24 ns | 16 B |
| `simplify(A, w)` vs `simplify!(A, w)` on fresh buffer | 181 vs 143 ns | 192 vs 112 B |

Allocation profile of the moment-like loop is dominated by four types:
`Vector{UInt8}` (word buffers), `Memory{UInt8}`, and the
`Vector{Tuple{ComplexF64,NormalMonomial}}` term vectors — i.e. **all transient
plumbing, none of it the actual result**.

## 2. Root causes

1. **Every monomial product funnels through full `Polynomial` construction.**
   `m1 * m2` does `vcat` → `simplify` → `_simplified_to_terms` (1-element
   `Vector` of tuples) → `Polynomial` → `_process_terms` (out-of-place `sort` +
   dedup). For Monoid/TwistedGroup algebras the product is *provably ≤ 1 term*,
   so the sort, the dedup, and two of the vectors are pure overhead.

2. **`simplify` copies buffers the caller already owns.**
   `simplify(A, w) = simplify!(A, copy(w))`. Hot paths call it on freshly
   allocated buffers (`vcat` in `arithmetic.jl`, `_neat_dot3` results in
   `moment.jl` lines ~471/535/555/576) — the copy is 100% waste. Measured:
   ~40 ns + 80 B per call.

3. **`NormalMonomial * Polynomial` accumulates quadratically.**
   `result = result + coef * prod_poly` inside the loop re-allocates, converts,
   concatenates and re-sorts the whole accumulated polynomial *per term*:
   ~2.8 KB per term vs ~190 B needed.

4. **Constructor validation runs on trusted words.**
   The only `NormalMonomial` constructor does `any(iszero, word)` +
   `_validate_word(A, word)` — O(n) scans — even when the word comes straight
   out of `simplify`, which already guarantees canonical form.

5. **Generic array hashing allocates.** `hash(::Vector{UInt8})` costs ~22 ns +
   16 B in Julia 1.12. Moment assembly and basis dedup do millions of dict
   probes, so this shows up as steady GC pressure.

6. **`Vector{T}` word storage.** Each monomial = 2 heap objects
   (`Vector` + `Memory`); `==`/`hash`/`isless` all pointer-chase. Typical
   canonical words are 1–6 bytes — smaller than the pointer that points at them.

## 3. Improved design

### Tier 1 — hot-path hygiene (no representation change) — *prototyped & measured*

| Operation | current | prototype | speedup |
|---|---|---|---|
| `m * m` | 0.11 µs / 432 B | 0.03 µs / 96 B | **3.7× / 4.5×** |
| `m * p` (27 terms) | 13.0 µs / 76.9 KB | 1.26 µs / 5.1 KB | **10× / 15×** |
| `p * p` (27×27) | 67 µs / 290 KB | 46 µs / 119 KB | 1.5× / 2.4× |

Prototype outputs verified `==` to current results. Changes:

1. **Trusted constructor.** Add a private inner-constructor path
   (e.g. `NormalMonomial{A,T}(word, ::Val{:trusted})` or an internal
   `_unchecked_monomial`) that skips zero-check and `_validate_word`. Use it
   everywhere the word provably comes from `simplify`/`simplify!`
   (`_simplified_to_terms`, adjoint, moment assembly). Public constructors keep
   validating — "validate at the boundary, trust inside".

2. **Single-term multiply kernel.** `_mul_term(m1, m2) -> (coef, mono)` for
   `MonoidAlgebra`/`TwistedGroupAlgebra`: allocate the concatenated buffer once
   (`copyto!` ×2, no `vcat` of two copies), `simplify!` in place, wrap with the
   trusted constructor. No intermediate term-`Vector`, no `Polynomial`, no
   sort. `m1 * m2` (public API) wraps this kernel in a `Polynomial` at the end.
   PBW algebras keep the multi-term path.

3. **Rewrite `m * p`, `p * m`, `p * p` as: one pre-sized terms buffer →
   push all products → one `Polynomial(terms)` at the end.** Kills the
   quadratic accumulation (cause 3).

4. **`_process_terms`: `sort!` in place and skip entirely for ≤1 term.**
   Every caller passes a freshly built vector it owns; the defensive `sort`
   copy and the no-op dedup pass for singletons are waste.

5. **Adopt `simplify!` on owned buffers** in `arithmetic.jl` and the four
   `moment.jl` call sites; keep `simplify` (copying) as the safe public entry.

6. **Custom `Base.hash(::NormalMonomial)`**: hash the word bytes directly
   (e.g. chunk-wise `hash` over `UInt64` reads or FNV-1a) instead of generic
   `AbstractArray` hashing → 0 allocations per probe.

Scope: `src/types/monomial.jl`, `src/types/arithmetic.jl`,
`src/types/polynomial.jl`, `src/optimization/moment.jl` (+ simplifier files
gain nothing — they're already in-place). Existing test suite is the
regression net; results are bit-identical.

### Tier 2 — moment-assembly buffer reuse + interning

The 80×80 moment-like loop allocates 36 MB because every entry builds fresh
buffers and fresh monomials, although the set of *distinct canonical words*
saturates quickly.

- `_neat_dot3!(scratch, a, m, b)`: write the triple product into a reusable
  scratch buffer (one per assembly loop; thread-local if parallelized), then
  `simplify!` in place.
- **Intern canonical words**: `get!` into a `Dict{Vector{T},NormalMonomial}`
  (copying the scratch only on first sight), or better, key the moment-index
  dict directly by the scratch word using a non-allocating probe and only
  materialize a `NormalMonomial` for novel words.
- Expected effect: steady-state assembly allocates only for *new* canonical
  words → tens of KB instead of tens of MB; removes the 2.4% GC overhead and
  most of the 63 B/probe dict cost.

### Tier 3 — representation change (benchmark-gated, larger blast radius)

Only after Tiers 1–2, and only if profiles still show monomial plumbing:

- **Option A — inline small words.** Store words as an isbits
  small-buffer type (`SmallCollections.SmallVector{N,T}` or hand-rolled
  `NTuple{N,T}` + length, N ≈ 16 covering degree-8 relaxations of 2-byte
  indices). Makes `==`/`hash`/`isless` register operations and removes both
  heap objects per monomial. Cost: a capacity parameter (or a spill-to-heap
  hybrid to stay type-stable), and churn across `Polynomial`, states, and
  composed monomials. General win across all algebras.
- **Option B — symplectic Pauli words.** For `PauliAlgebra` specifically,
  represent a word as (x-mask, z-mask, phase) bit-vectors
  (`UInt64`/`BitVector` per mask). Products become XOR + popcount — O(1)
  instead of the current sort-based O(n log n) `simplify!`. This is the
  standard stabilizer-formalism trick (cf. QuantumClifford.jl) and would make
  condensed-matter workloads (the dominant Pauli use case) dramatically
  faster. Can be introduced as an internal assembly representation without
  changing the public `NormalMonomial{PauliAlgebra}` type.

### Non-goals / rejected

- **Hash-consing all monomials globally**: pointer-equality wins don't justify
  global mutable state + thread-safety burden; Tier-2 local interning gets the
  same benefit where it matters.
- **Making `Polynomial` lazy/unsorted**: the sorted-unique invariant is load-
  bearing for binary search and canonical forms; fix the *paths into* the
  invariant instead.



## 4. Triage (2026-06-12, branch `perf/monomial-allocations`)

Follow-up measurements to decide scope, plus an adversarial review pass.

### New evidence

| Measurement | Result |
|---|---|
| End-to-end `cs_nctssos`, Heisenberg N=6, order 2, COSMO | 30 s, 2.3 GB, GC 0.4%. Stack-attributed: COSMO ≈45% + 51% Base (solver linalg), JuMP/MOI/Dualization 4%, **NCTSSoS ≈0.1%** |
| Build-only (`compute_sparsity` + `moment_relax`) scaling | N=6: 0.05 s / 90 MB · N=10: 0.46 s / 738 MB · N=14: 2.6 s / 2.9 GB, **26% GC** · N=18: 8.3 s / 8.0 GB, **29% GC** — GC fraction climbs with N |
| Build-only allocation ownership (N=14, `Profile.Allocs`) | NCTSSoS assembly logic 59.7%, monomial/poly plumbing 34.7%, other NCTSSoS 5.6%, JuMP/MOI ≈0% — **94% NCTSSoS-owned** |
| Phase split (N=14) | `compute_sparsity` 0.21 s / 0.4 GB vs **`moment_relax` 2.7 s / 2.5 GB** — the villain is moment assembly |
| Interning keystone (distinct canonical words / triple products) | N=6: 0.94% · N=10: 3.4% — **≥96% interning hit rate, assumption holds** |

### Decisions

- **Tier 1 — GO now**, with guardrails from review:
  - Public `Polynomial(::Vector)` constructors stay **non-mutating**; the
    in-place `sort!` lives in an internal owned-buffer path
    (`_process_terms!`-style), used only at audited call sites.
  - Trusted constructor stays private; document per-algebra `simplify!`
    postconditions; never wrap a reusable scratch buffer without copying.
  - Custom hash must keep the contract: incorporate algebra type **and** the
    seed `h`, preserve `isequal`/`hash` consistency; add cross-algebra
    collision tests.
  - Differential property tests across **all** algebras (random words,
    duplicate-heavy term lists, zero coefficients, PBW expansions) — not just
    the Pauli/Heisenberg path the prototype used.
- **Tier 2 — GO, as a separate profile-directed sprint after Tier 1 lands.**
  Promotion is now evidence-backed (8 GB churn / 29% GC at N=18, 94%
  NCTSSoS-owned, ≥96% interning hit rate, `moment_relax` identified as the
  dominant phase). Guardrails:
  - Line-level allocation table of `moment_relax` *before* implementing;
    target the top sites, not the micro-paths by assumption.
  - Never store the scratch buffer as a `Dict` key — probe by scratch,
    store a copy on first sight (mutable-key footgun).
  - Track **peak RSS**, not just allocation churn; interning retains memory.
  - Acceptance gates: N=18 build ≤ 2 GB churn, GC < 10%, wall ≥ 2× faster,
    peak RSS not worse.
- **Tier 3 — defer implementation; define triggers instead of a flat no.**
  Re-profile after Tier 2. Trigger a representation *spike* (not a rewrite) if
  word storage/hash/`==` still exceeds ~25% of build allocations, or if N≥30
  Pauli workloads remain build-bound with the production solver (Mosek).
  Option B (symplectic Pauli) remains the only order-of-magnitude candidate
  and can be scoped as an internal assembly representation only.

### Verdict on “should we push for even less allocation?”

**Yes for the build path, no for the representation rewrite — yet.**
End-to-end wall clock is solver-dominated (NCTSSoS ≈0.1% of a COSMO solve),
so zero-allocation as a goal in itself buys nothing user-visible there. But the
build phase is 94% NCTSSoS-owned allocation with GC% climbing in N — that’s
real headroom and a scaling risk for N=30–50 targets. Tier 1 + Tier 2 attack
exactly that with measured justification; Tier 3 stays gated on post-Tier-2
profiles with the production solver.

### Evidence limits (kept honest)

- The end-to-end ratio is COSMO-specific and N=6-sized; Mosek shrinks solve
  time and shifts the ratio toward build. Re-measure with Mosek before any
  Tier 3 call.
- Allocation churn ≠ peak RSS; the 8 GB figure is GC pressure, not footprint.
- Profile sample attribution misses BLAS threads.
- The interning saturation was measured on Heisenberg-type supports; check one
  non-chain Hamiltonian during the Tier 2 sprint.

## 5. Reproduction

Profiling scripts used (not committed): `/tmp/mono_profile.jl`,
`/tmp/mono_profile2.jl`, `/tmp/mono_proto.jl` — Pauli N=10 setup as above;
`@timed` loops + `Profile.Allocs` with `sample_rate=0.1`. Triage scripts:
`/tmp/e2e_profile2.jl` (stack-attributed end-to-end), `/tmp/build_scaling.jl`,
`/tmp/build_allocs.jl`, `/tmp/triage_extra.jl` (phase split + saturation).

## 6. Incidental findings

- The inner `Polynomial{A,T,C}` constructor's docstring says it "enforces
  invariants: sorted, deduplicated, non-zero coefficients" — it does not; it
  stores terms as-given (the outer constructor processes). The comment should
  state the actual contract: *trusted path, caller guarantees invariants*.
- Two outer `Polynomial(::Vector{Tuple{C,NormalMonomial}})` constructors exist
  (one `C<:Number` constrained, one unconstrained `CIn` that skips
  `_process_terms`). The unconstrained one is unreachable for valid
  coefficient types (the struct bounds `C<:Number`, and the constrained
  method is more specific) — dead code; delete it while touching this file.
