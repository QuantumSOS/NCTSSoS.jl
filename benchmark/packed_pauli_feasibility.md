# Bit-Packed Symplectic Pauli Representation — Feasibility & Blast-Radius Analysis

Date: 2026-06-12 · branch `perf/monomial-allocations` (post round 3, commit `41632d6`)
Status: **analysis only — no implementation**. Companion to `monomial_profile_report.md`.

## Verdict

**Go — but not as originally proposed.** The "parallel monomial type behind
`AbstractMonomial{PauliAlgebra}`" framing is refuted by the audit (§3): it is a
type-system rewrite, not a slot-in. The viable design is a **packed scratch
domain**: an internal isbits `PackedPauli` used *only inside the six hot
kernels*, converting at the boundaries, with words interned so each unique
monomial allocates once. This captures the dominant costs (O(L) hashing, O(L)
compares, per-product word copies, `simplify!`'s per-product site sort) with
**zero public type changes, zero fixture changes, and bit-identical digest A/B
preserved**.

- Phase 1 (packed kernels): go. Estimated 2–3× on the build benchmark.
- Phase 2 (packed moment keys, `K = PackedPauli`): optional follow-up, separate
  audit; needs a key-canonicalizing digest variant.
- Full parallel monomial type: **no-go** (cost/risk wildly out of proportion).

## 1. What the audit found

### 1.1 The killer constraint: `Polynomial` hard-codes `NormalMonomial`

```julia
struct Polynomial{A,T,C} <: AbstractPolynomial{C}
    terms::Vector{Tuple{C,NormalMonomial{A,T}}}   # src/types/polynomial.jl:79
end
```

A `PackedPauliMonomial` cannot be stored in any `Polynomial` without changing
this field — i.e. adding a monomial type parameter to `Polynomial` (public
type, used everywhere as `Polynomial{A,T,C}`) or splitting the struct. On top
of that, **60+ method signatures** are constrained `M<:NormalMonomial{A,T}`
(sparsity, moment, states, GNS, symmetry). `AbstractMonomial{A,T}` exists as a
supertype, but the codebase is *not* written against it; it is written against
`NormalMonomial`. The "designed for this dispatch" claim in the proposal is
true for *algebra* dispatch (parameter `A`), false for *representation*
dispatch (the word layout).

Making `NormalMonomial`'s field parametric (`NormalMonomial{A,T,W}`) is worse:
`NormalMonomial{A,T}` stops being concrete, every container type churns.

### 1.2 Consumers of Pauli words (`.word` access / word semantics)

| consumer | hot? | affected by packed scratch domain? |
|---|---|---|
| `monomial.jl` hash/`==`/`isless`/`show`/`iterate`/`adjoint` | hot via kernels | no — `NormalMonomial` unchanged |
| `arithmetic.jl` `_mul_term`, `_simplified_*`, `_copy_if_buffer_aliased` | warm (user algebra) | no |
| `moment.jl` `_build_constraint_matrix`, `_moment_matrix_basis`, `_polynomial_total_basis` | **hot** | **specialized** |
| `sparsity.jl` `init_activated_supp`, `get_term_sparsity_graph`, `term_sparsity_graph_supp` | **hot** | **specialized** |
| `moment.jl` `_moment_key` (Pauli: aliases `mono.word`, K=`Vector{T}`) | warm | unchanged (phase 1); phase-2 target |
| `moment.jl:1117` `_append_moment_eq_constraints!` | hot only with moment-eq constraints | not specialized (stays correct, just unaccelerated) |
| `canonicalization.jl` `symmetric_canon(Pauli) = copy(word)` (identity rep) | cold | no |
| `gns.jl` (allocating `neat_dot`/`_neat_dot3`) | cold (post-solve) | no |
| `symmetry.jl`, `sympleq/` (Clifford recognizer; already has a symplectic tableau!) | cold (setup) | no |
| `states/` (`NCStateWord`, `StateSymbol`) | — | **Pauli explicitly unsupported** (`states/word.jl:212`) |
| `v2rdm_structured.jl`, `particle_number.jl` | — | fermionic-only |
| `composed.jl` (tensor products with Pauli factor) | warm | no — generic path |
| `registry.jl` `create_pauli_variables` (`T` = smallest UInt for `3n` indices) | cold | no |
| test fixtures (`test/data/expectations/*.toml`) | — | objective values only; no word bytes, no basis orderings |
| tests hard-coding `.word` bytes (`test/polynomials/*.jl`) | — | construct `NormalMonomial` directly; unchanged |

Side note: `src/sympleq/tableau.jl` already converts Pauli words ↔ binary
symplectic rows (`Matrix{UInt8}`, one byte per bit) for the Clifford-symmetry
recognizer. It is prior art *inside the repo* for the encoding and the `Y=iXZ`
η-convention, but it is matrix-shaped and cold-path; not reusable as the kernel
type, though the packed type should reuse its conventions.

### 1.3 Where the remaining 11.2 s / 5.77 GiB actually lives

Top CPU (`hash`→`_hash_word_vector`) is Dict/Set probes of `NormalMonomial`s in
exactly the six kernels above (`_remember_buffered_monomials!` probes, `seen`
sets in `gsupp!`/`init_activated_supp`, `insorted` binary searches with O(L)
`isless`, `_process_terms!`/`sort!` with O(L) compares). Top allocation
(`_push_scaled_buffered_terms!` insertion `copy(word)`) is in
`_build_constraint_matrix`. Also unaddressed by round 3: `simplify!`'s
per-product `InsertionSort` by site. All of these are *inside* the proposed
specialization boundary. What is *outside*: `_linearize_moment_polynomial`'s
`Dict{Vector{T},…}` churn (generic array hashing allocates on 1.12 → phase 2),
clique decomposition, JuMP build, solver.

## 2. Ordering — resolved exactly, no re-baselining

Current Pauli `isless`: **graded** (degree = word length first), then **lex on
the index vector**; indices are `3(site−1)+type+1` with X=0<Y=1<Z=2 and words
site-sorted, ≤1 letter/site.

Packed reproduction (machine-verified, §5): degree = `popcount(x|z)`. For equal
degree, the words share a prefix up to the first site where the masks differ,
`s = trailing_zeros((x₁⊻x₂)|(z₁⊻z₂))`. At `s`: if exactly one word has support
→ it is smaller; if both → smaller type key wins, with the monotone key
`2z − (x∧z)` (X→0, Y→1, Z→2). O(1) for ≤64 sites, short limb loop beyond.

Consequence: every sorted structure (`_sorted_basis_keys`, `sorted_union`,
`_process_terms!` order inside polynomials, `insorted` searches) is reproduced
**bit-identically**. The digest A/B methodology survives phase 1 unchanged.

Phase-2 note: `key_lt` on `Vector{T}` keys is a *different* order (ungraded
zip-lex + length tie-break). Also O(1)-reproducible from masks (first-diff-site
rule with an "is the absent word exhausted at sites ≥ s" prefix check via
`support >> s == 0`), so even packed keys can keep today's key order; only the
key *serialization* changes, so the digest script needs a
canonicalize-keys-to-words variant for phase 2.

## 3. Integration seam — three designs compared

| | A: parallel monomial type | B: parametric word field | C: packed scratch domain |
|---|---|---|---|
| `Polynomial` change | new type param (public break) | none, but `NormalMonomial{A,T}` non-concrete | **none** |
| signatures touched | 60+ | 100+ | **0 generic; 6 new specialized methods** |
| states/GNS/symmetry/printing | all need methods | all churn | **untouched** |
| fixtures/tests | re-baseline risk | re-baseline risk | **bit-identical** |
| captures hash/copy/sort wins | yes | yes | **yes (kernels are where they live)** |
| captures `Polynomial` storage win | yes | yes | partially (interned, shared words) |
| est. diff | ≥2000 lines modified | worse | ~1000 lines, almost all *added* |
| verdict | no-go | no-go | **go** |

### Design C concretely

New internal file `src/types/packed_pauli.jl`:

```julia
struct PackedPauli{K}            # isbits; K limbs of 64 sites
    x::NTuple{K,UInt64}
    z::NTuple{K,UInt64}
end
# product: x₃=x₁⊻x₂, z₃=z₁⊻z₂
# phase:  k = pc(x₁∧z₁) + pc(x₂∧z₂) + 2·pc(z₁∧x₂) − pc(x₃∧z₃)  (mod 4)
# adjoint = identity (canonical Pauli words are Hermitian; matches
#           `adjoint(::NormalMonomial{PauliAlgebra}) = copy(m)`)
# hash = O(1) integer mix; isless = graded comparator of §2
```

Six Pauli-specialized methods, dispatching on the *more specific* argument
types (`Polynomial{PauliAlgebra,T,C}` / `Vector{NormalMonomial{PauliAlgebra,T}}`),
generic methods untouched:

1. `_build_constraint_matrix` — pack basis + poly terms once; accumulate
   `(coef, PackedPauli)` per element; sort/merge with the §2 comparator
   (replaces `_polynomial_from_owned_terms!`'s O(L) compares); materialize
   `Polynomial` via an intern cache `Dict{PackedPauli,NormalMonomial}` so each
   unique word across the whole build allocates **once** and is shared by
   reference (words are immutable by documented contract — sharing is safe;
   today's per-product `copy` exists for buffer aliasing, not uniqueness).
2. `_moment_matrix_basis` — `Dict{PackedPauli,Nothing}` dedup, O(1) probes,
   sort packed, unpack/intern at the end.
3. `_polynomial_total_basis` — same.
4. `init_activated_supp` — packed diagonal products + packed `sorted_union`.
5. `get_term_sparsity_graph` — packed products, `insorted` on packed reps.
6. `term_sparsity_graph_supp` — packed `seen` set; **`out` keeps
   first-occurrence order**, which depends only on loop order, not hash order
   (loop order identical → output identical).

Triple product `a†·m·b`: `a† = a` (the generic path's reverse+stable-site-sort
is a no-op modulo phase, and distinct-site letters commute), so it's two packed
multiplications. The `phase_k ∈ 0:3` slots into the existing TwistedGroup
contract; the `0x04` zero sentinel is never produced — correct, since the Pauli
group has no zero divisors, and today's `simplify!(PauliAlgebra,…)` never
returns it either. `_conj_coef(PauliAlgebra, k) = (0x04 − k) & 0x03` acts on
the phase alone — representation-independent.

Pack width: choose `K = cld(n_sites, 64)` at kernel entry (max site scan of the
basis is O(B·L), trivial); a function barrier keeps the kernels concretely
typed. `K=1` covers n≤64 — every current test and benchmark.

## 4. Things that break / need care, and the handling

| item | handling |
|---|---|
| public API / types | nothing breaks (phase 1 adds internal type + 6 methods) |
| digest A/B | bit-identical in phase 1 (order reproduced exactly); phase 2 needs key-canonicalizing variant |
| Dict iteration order | hash values change, but no kernel output depends on hash order (outputs are sorted or first-occurrence-by-loop-order — audited) |
| n > 64 sites | `NTuple{K,UInt64}` limbs; explicit boundary tests at 64/65 |
| code duplication (6 kernels × 2 paths) | `invoke`-based equivalence tests in CI: run the generic method via `invoke` on small Pauli problems, compare structures element-wise against the specialized path |
| moment-eq-constraint path, GNS, symmetry, states, composed | intentionally not specialized; correct via generic path |
| word-vector sharing via interning | safe under the existing immutability contract (`_unchecked_monomial` docstring); `copy(m)` still deep-copies |
| solver-backed fixtures | untouched — same SDP, same orderings |

## 5. Verification already performed (this analysis)

Machine-checked against the live oracles in this repo (Julia, 200k randomized
trials each + exhaustive single-site table, **0 failures**):

- **Phase formula** `k = pc(x₁∧z₁)+pc(x₂∧z₂)+2·pc(z₁∧x₂)−pc(x₃∧z₃) mod 4`
  vs `simplify!(PauliAlgebra, vcat(wa,wb))` — words *and* phases match,
  including identity words, up to 60 sites.
- **Graded comparator** of §2 vs `isless(::NormalMonomial{PauliAlgebra},…)`.

## 6. Implementation plan (when approved)

1. **Phase 0** — land `PackedPauli{K}` + ops + unit/property tests (the §5
   scripts, promoted into `test/polynomials/packed_pauli.jl`, plus limb
   boundary and comparator-totality tests). No call sites. ~350 src / ~250 test
   lines.
2. **Phase 1a** — specialize the three sparsity kernels. Digest A/B
   (sparsity structures) bit-identical; narrow suites
   (`test/correlated_sparsity/`, `test/problems/condensed_matter/`).
3. **Phase 1b** — specialize the three moment kernels + intern cache. Digest
   A/B (full `MomentLinearData`) bit-identical; `make test`; re-run
   `benchmark/monomial_profile.jl`, update the report table. ~450–500 src lines
   across `moment.jl`/`sparsity.jl`.
4. **Phase 2 (optional, separate audit)** — `K = PackedPauli` moment keys:
   `_moment_key` specialization, `key_lt`/`key_isequal` methods reproducing
   today's vector-key order, `_available_moment_set`/`_missing_polynomial_monomials`
   (moment.jl:838–857) updates, audit of `lowering.jl` resolver and
   dualization/SOS key consumers, key-canonicalizing digest variant.

Expected: phase 1 ≈ 2–3× on the Heisenberg n=10 order=3 build (11.2 s → ~4–6 s),
phase 2 a further ~1.3–1.5×. Pauli-only; NC/projector/unipotent/PBW unaffected.
Do **not** promise order-of-magnitude end-to-end: clique decomposition, JuMP
build, and the solver are untouched (Amdahl).
