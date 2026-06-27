# Matching QMBCertify Performance Without Losing Generality

## Executive Summary

QMBCertify (Jie Wang) is a hand-tuned Julia script that certifies ground-state
properties of the Heisenberg model on periodic lattices.  It handles 1D chains
up to L=30 and 2D square lattices up to 4×4.  The companion paper
(arXiv:2604.01555) scales to N=100 (1D) and 16×16 (2D).

NCTSSoS.jl currently solves N=16 (sparse d=4, ~12 min) and N=24 (order-2
singlet, ~6 min) on the same 1D Heisenberg benchmark.  The bottleneck is
symbolic construction time, not the SDP solver.

This document analyzes *why* QMBCertify is fast, extracts the transferable
architectural lessons, and lays out a phased plan to achieve the same speed
in NCTSSoS.jl without sacrificing its generality across algebras, lattices,
and problem types.

---

## Part I: Why QMBCertify Is Fast

QMBCertify achieves its speed by skipping three layers of abstraction that
NCTSSoS pays for on every problem.

### 1. No symbolic moment matrix

NCTSSoS builds a full `Matrix{Polynomial}` — every entry is a polynomial in
moment variables.  Then it scans this matrix to discover unique moments, builds
an index table, and converts to JuMP `AffExpr`.  For an n×n block, that is
O(n²) polynomial objects, each requiring allocation, simplification, and
hashing.

QMBCertify never builds this matrix.  It maintains a flat vector
`cons = [AffExpr(0) for i in 1:length(tsupp)]` — one entry per unique moment
— and accumulates JuMP coefficients directly:

```julia
add_to_expression!(cons[Locb], coef, gram[l][j,k])
```

Each (row, col) pair contributes a scalar coefficient to a known moment index.
No polynomial object is ever created.

### 2. Analytic DFT instead of generic representation theory

NCTSSoS uses SymbolicWedderburn to discover irreducible representations of the
symmetry group on the basis.  For a cyclic group of order N on a basis of size
12,001, this is massively redundant — the answer is the discrete Fourier
transform, known analytically since the 18th century.

QMBCertify writes `cos(2πrk/L)` and `sin(2πrk/L)` inline.  No character
table, no projection matrices, no Wedderburn decomposition.  The
`eigen_circmat` function is a few lines of trigonometry.

### 3. In-place integer-array monomials

Every monomial in QMBCertify is a `Vector{UInt16}` — a sorted list of encoded
site×axis indices.  The `reduce!` function chain does all simplification by
mutating this array in place:

- `reduce1!` — bubble-sort by site
- `reduce3!` — cancel σ²=I pairs
- `reduce2!` — apply σ_x σ_y = iσ_z rules
- `reduce4` — minimize over translation/reflection/permutation orbits
- `isz` — kill sign-symmetry-zero monomials

No heap allocation per operation.  No hash table of canonical forms.  No
polynomial arithmetic.

### Quantitative consequence

| Operation | NCTSSoS (N=64, d=4) | QMBCertify equivalent |
|:--|:--|:--|
| Build basis | Type-parameterized `NormalMonomial`, sorted sets | `get_basis` returns `Vector{Vector{UInt16}}` |
| Simplify products | `simplify(PauliAlgebra, word)` with dispatch | `reduce!` mutating `UInt16[]` |
| Find orbits | SymbolicWedderburn group enumeration | `reduce4`: try all translations, take min |
| Block-diagonalize | Generic Wedderburn irrep decomposition (88 GiB) | Inline DFT coefficients (kilobytes) |
| Emit JuMP model | `MomentProblem` → `MomentLinearData` → `build_jump_model` | Direct `@variable` + `add_to_expression!` |
| **Total construction** | **1121 seconds, 88 GiB** | **Seconds, megabytes** |

### The tradeoff QMBCertify makes

QMBCertify only solves *one* problem: the Heisenberg model on periodic
lattices with Pauli operators.  Change the algebra (fermions, bosons), the
lattice (non-periodic, irregular graph), the constraints (projector identities,
Bell inequalities), or the sparsity pattern — and QMBCertify has nothing to
offer.  It is a 1500-line script that knows every answer in advance.

The speed difference is not an algorithm gap.  It is the difference between a
library and a script.

---

## Part II: Transferable Architectural Lessons

Each of QMBCertify's three speed advantages can be made available in a general
framework without sacrificing generality.  The pattern is: **keep the generic
path as fallback; add fast paths that activate when structure is recognized.**

### Lesson 1: Streaming accumulation (bypass symbolic matrices)

**Pattern:** Define a protocol where the bilinear expansion loop *visits* each
(row, col, monomial, coefficient) tuple and directly accumulates into the
output, instead of materializing a `Matrix{Polynomial}`.

```
enumerate_basis_pairs(basis, constraints) do row, col, mono, coef
    idx = moment_index(mono)       # lookup, not discovery
    accumulate!(output, idx, coef, block, row, col)
end
```

The `output` can be:
- A `MomentProblem` (generic path — current behavior)
- A raw JuMP model with `AffExpr` vector (fast path — QMBCertify style)

The enumeration logic is shared.  What changes is whether intermediate
polynomial objects are stored or scalar coefficients are streamed.

**This is the single biggest win.** It eliminates the dominant allocation:
millions of intermediate `Polynomial` objects that exist only to be
immediately decomposed.

### Lesson 2: Analytic fast paths for known group families

**Pattern:** Most physically relevant symmetry groups fall into a small
taxonomy with known representation theory:

| Group family | Analytic decomposition | Cost |
|:--|:--|:--|
| Cyclic Z_N | DFT | O(N) |
| Dihedral D_N | Real DFT + parity split | O(N) |
| Symmetric S_n | Young tableaux | Known |
| Direct products | Product of factors | Multiplicative |
| Arbitrary finite | Wedderburn (generic fallback) | O(\|G\| · \|basis\|²) |

When the user passes `pauli_site_permutation([2:N; 1])`, the system can
*recognize* this as a cyclic generator and dispatch to the DFT path
automatically.  The generic Wedderburn path stays as a fallback for groups
that don't match any known family.

This is multiple dispatch applied to the symmetry layer:

```julia
decompose(::CyclicGroup, basis)      # → DFT
decompose(::DihedralGroup, basis)     # → real DFT + parity
decompose(::AbstractGroup, basis)     # → SymbolicWedderburn fallback
```

### Lesson 3: Buffer reuse in hot loops

**Pattern:** The cost of `NormalMonomial{PauliAlgebra,T}` is not the type tag
— it is that NCTSSoS allocates a new `Vector` for every intermediate product.
The fix is buffer reuse, not abandoning types:

```julia
# Current: allocates a new Vector for every product
function simplify(::Type{PauliAlgebra}, word::Vector{T}) ... end

# Fast: reuse a thread-local buffer
function simplify!(buf::Vector{T}, ::Type{PauliAlgebra}, word::Vector{T}) ... end
```

QMBCertify gets this for free because `reduce!` mutates in-place.  NCTSSoS
can do the same inside its hot loops without changing the public API.

### Lesson 4: Precomputed product tables

For a fixed basis, the result of `reduce(b_i† · b_j)` → `(moment_index,
coefficient)` is deterministic.  Cache it once; reuse across all momentum
blocks that share the same orbit structure.

QMBCertify does this implicitly by precomputing `coe1`, `coe2`, `bi1`, `bi2`
arrays before entering the JuMP construction loop.

---

## Part III: What NOT To Do

- **Don't hardcode block indices.**  QMBCertify's `posepsd8!` has hardcoded
  arrays `blocks = [[1], [2,3,5,9,...], ...]`.  This is fragile and
  Heisenberg-specific.  Instead, compute SU(2) block decompositions from
  Clebsch–Gordan coupling at startup — a one-time O(2^k) cost that
  generalizes to any SU(2)-invariant problem.

- **Don't abandon the type system.**  `NormalMonomial{A,T}` with compile-time
  algebra dispatch is the right design.  The problem is allocation patterns in
  hot loops, not the type hierarchy.

- **Don't write a separate "Heisenberg solver."**  QMBCertify is a dead end
  architecturally — it can never grow beyond its one problem.  Instead, build
  *general mechanisms* (streaming accumulation, analytic group decomposition,
  in-place simplification) that happen to make the Heisenberg case fast.

---

## Part IV: Current State of NCTSSoS.jl

### What already exists

| Capability | Status | Location |
|:--|:--|:--|
| Pauli algebra + simplification | ✅ Production | `src/simplification/pauli.jl` |
| Translation-invariant relaxation | ✅ Works at small N | `src/optimization/pauli_chains.jl` |
| Contiguous chain basis | ✅ Correct sizes at N=100 | `pauli_contiguous_chain_basis` |
| Charge sector splitting | ✅ Order-2 | `PauliChargeSectorSpec` |
| Singlet constraints | ✅ Order-2 | `PauliSingletConstraintSpec` |
| Pauli sign symmetry | ✅ | `pauli_sign_symmetry` |
| Generic Wedderburn symmetry | ✅ But too slow for large N | `src/optimization/symmetry.jl` |
| Fermionic v2RDM | ✅ | `src/optimization/v2rdm_structured.jl` |

### What is missing

| Capability | Gap |
|:--|:--|
| Streaming JuMP emission | Builds full `Matrix{Polynomial}` intermediate |
| Analytic DFT for cyclic groups | Uses generic Wedderburn |
| Pauli k-RDM positivity | Only fermionic v2RDM exists |
| State optimality constraints | None |
| SU(2) constraints beyond order 2 | Only singlet equalities |
| Rigorous certification | None |
| 2D lattice symmetry | None |

### Measured bottlenecks

From `perf/2604_01555_heisenberg_gap.md`:

| N | d | Generic Wedderburn time | Allocation | Max block | Target max block |
|--:|--:|--:|--:|--:|--:|
| 16 | 4 | 8.8 s | 1.08 GiB | 30 | 31 |
| 32 | 4 | 75.7 s | 7.17 GiB | 30 | 31 |
| 64 | 4 | 1121 s | 88.1 GiB | 30 | 31 |

The block sizes are correct.  The construction cost is 100–1000× too high.

From the specialized translation path (already in the repo):

| N | d | Path | Max block (realified) |
|--:|--:|:--|--:|
| 20 | 4 | `pauli_translation_invariant_moment_relaxation` | 62 |

This path produces the right structure but currently rejects additional
constraints (RDM, state optimality).

---

## Part V: Phased Implementation Plan

### Benchmark targets

Two tiers, never mixed:

**Tier A — QMBCertify code parity:**
- 1D Heisenberg chain, L=30, same options as QMBCertify defaults
- Pass: `|E_NCTSSoS − E_QMBCertify| ≤ 1e-6`, runtime ≤ 2× QMBCertify, RSS ≤ 2× QMBCertify

**Tier B — Paper structural parity:**
- 1D Heisenberg chain, N=100, d=4, sparse contiguous basis
- Pass: basis = 12,001, orbit reps = 121, logical max block = 31, structural build ≤ 10s

Every milestone uses exact numerical targets.  "Competitive" is not allowed.

### Degree closure table

Before implementing any constraint, verify moment coverage:

| Feature | Required moments | Closure rule (d=4, h=2) |
|:--|:--|:--|
| Base moment PSD | u†v for u,v ∈ B_d | Generated by construction |
| Objective H | Pauli words in H | d ≥ 1 sufficient |
| k-site contiguous RDM | All Pauli strings on k sites | k ≤ 2d (i.e. k ≤ 8) |
| Linear state-opt ⟨[H,A]⟩=0 | Commutator words | test width a: a+h-1 ≤ 2d → a ≤ 7 |
| PSD state-opt A†[H,[H,A]]≥0 | Triple products | 2a+h ≤ 2d → a ≤ 3 |
| SU(2) equalities | No new moments | Existing moments suffice |

Implementation rule: add `closure_check(feature, basis)` that enumerates
required monomials after translation reduction and fails before model
construction if any are missing.

---

### Phase 0 — Lock references and success criteria

**Duration:** 2–4 days

**Actions:**
- Pin QMBCertify commit, NCTSSoS commit, Julia version, Mosek version, hardware
- Run QMBCertify on reference hardware for L=10,20,30 and record:
  objective, PSD block histogram, moment count, build time, solve time, peak RSS
- Write TOML reference files under `test/data/expectations/`:
  `heisenberg_qmbcertify_base.toml`, `heisenberg_qmbcertify_rdm.toml`

**Milestone:** Reference TOML files committed with exact numerical targets.

**Stop gate:** None.  This is prerequisite bookkeeping.

---

### Phase 1 — Streaming accumulation + analytic DFT

**Duration:** 1–2 weeks

**Goal:** N=100 d=4 structural decomposition in <10s; L=30 end-to-end solve
matching QMBCertify base bound.

**Actions:**

1. **Define `AbstractModelEmitter` protocol**
   - `StreamingJuMPEmitter`: accumulates directly into `AffExpr` vector
   - `SymbolicEmitter`: builds `MomentProblem` (current behavior, for
     validation)
   - Both implement `accumulate!(emitter, moment_idx, coef, block, row, col)`

2. **Implement analytic cyclic DFT decomposition**
   - Detect cyclic generators → bypass SymbolicWedderburn
   - Produce momentum-sector blocks with DFT coefficients
   - Handle k=0, k=N/2, and conjugate pairs for realification
   - Validate: match SymbolicWedderburn output at N=4,6,8

3. **Add in-place simplification buffers**
   - Thread-local reusable `Vector{T}` for `simplify!`
   - Avoid allocation in the bilinear expansion inner loop

4. **Precompute product tables**
   - Cache `reduce(b_i† · b_j)` → `(moment_index, coefficient)` per orbit
   - Reuse across all momentum blocks

**Milestones:**
- N=100 d=4: basis 12,001, orbit reps 121, max block 31, build ≤ 10s, ≤ 2 GiB
- N=8,10,12: streaming path optimum matches generic path within 1e-8
- L=30 base solve matches QMBCertify within 1e-6

**Stop gate:** If N=100 structural build exceeds 30s after optimization,
profile before proceeding to Phase 2.

---

### Phase 2 — RDM and state optimality constraints

**Duration:** 1–3 weeks

**Goal:** Match QMBCertify bound quality, not just structural size.

**Actions:**

1. **Extend streaming path to accept additional PSD blocks and linear constraints**
   - The current `pauli_translation_invariant_moment_relaxation` rejects all
     constraints; lift this restriction for the specific constraint types below

2. **Implement Pauli k-site contiguous RDM positivity**
   - RDM from global moments: ρ_S = 2^{-k} Σ_P y_P P
   - Decompose under U(1) magnetization: ρ = ⊕_m ρ^(m), each block PSD
   - Start with k=2,3,4; extend to k=8 only after closure verification
   - Generate block structure from Clebsch–Gordan coupling, not hardcoded tables

3. **Add first-order linear state optimality**
   - ⟨i[H,A]⟩ = 0 for A in a selected monomial family
   - Family: all monomials of test width a ≤ 7 (for d=4) that are
     translation-orbit-reduced and sign-invariant
   - Implement `closure_check` that rejects a > 2d - h + 1

4. **Add PSD state optimality (double commutator)**
   - A†[H,[H,A]] ≥ 0 for test width a ≤ 3 (at d=4)
   - Produces additional small PSD blocks per momentum sector
   - Implement `closure_check` that rejects 2a + h > 2d

**Milestones:**
- `closure_check` rejects invalid d/k/a combinations before model construction
- L=30 + RDM + state-opt matches QMBCertify bound within 1e-6
- Solver burden (time, RSS) measured separately from construction burden

**Gate:** Try N=100 with RDM only if L=30 extrapolates below 64 GiB peak RSS.
Otherwise Phase 3 is mandatory before scaling.

---

### Phase 3 — SU(2) from representation theory

**Duration:** 3–6 weeks

**Goal:** Further block reduction via continuous SU(2) symmetry.

**Actions:**

1. **Implement Clebsch–Gordan decomposition for Pauli operator basis**
   - Each non-identity Pauli transforms as spin-1 under SU(2)
   - Couple tensor products into total (J, m, α) basis
   - Moment matrix: M = ⊕_J (I_{2J+1} ⊗ X_J), with X_J ⪰ 0

2. **Decompose RDM blocks under SU(2)**
   - Physical Hilbert space: (ℂ²)^{⊗k} = ⊕_j V_j ⊗ M_j
   - PSD condition on multiplicity matrices M_j
   - Replaces hardcoded block indices with derived decomposition

3. **Extend singlet constraints beyond order 2**
   - Higher-degree SU(2) equalities from irreducible tensor decomposition
   - Not just axis equality (xx=yy=zz) but full singlet channel relations

**Milestones:**
- For order ≤ 2: generated SU(2) blocks reproduce existing singlet results
- CG transform unitarity residual ≤ 1e-12
- Off-block residual ≤ 1e-10
- L=30 + SU(2) runtime ≤ 2× QMBCertify or ≥ 2× improvement over Phase 2

**Stop gate:** Only claim N=100 numerical parity after this phase succeeds
at L=30.

---

### Phase 4 — Rigorous certification

**Duration:** 2–4 weeks

**Goal:** Computer-verified ground state energy bounds.

**Design decision:** Certify in the **dual**, not by projecting primal moments.

**Actions:**

1. **Extract reduced dual PSD blocks Q_i from solver**
2. **Verify the SOS identity:**
   H − λ = Σ_i v_i† Q_i v_i + (equality multipliers) + (residual)
3. **Symmetrize dual blocks by group averaging** (preserves PSD)
4. **Round to rationals; re-verify**
5. **Compute rigorous eigenvalue lower bounds via Arblib.jl**
6. **Report certified bound:** λ_cert = λ_num − ε, where ε is the certified
   residual norm

**Critical constraint:** Never project onto affine constraints if that can
leave the PSD cone.  The dual approach avoids this: group averaging of a PSD
matrix stays PSD.

**Milestones:**
- Small N (4,6,8): all dual blocks PSD with certified min eigenvalue margin
- L=30: certified bound gap ≤ 1e-6 from numerical bound
- Rational arithmetic certificate passes Arblib verification

---

### Phase 5 — 2D lattice (separate project)

**Duration:** 4–8+ weeks, after 1D succeeds

**Goal:** Match QMBCertify 4×4, then attempt paper 16×16.

**Actions:**

1. **Start with QMBCertify 4×4, not paper 16×16**
2. **Implement rectangular lattice symmetry spec:**
   - Translation Z_Lx × Z_Ly → 2D DFT
   - Point group D4 → momentum-star stabilizers → small finite groups per star
3. **Reuse RDM / SU(2) / certification from 1D phases**
4. **Validate against QMBCertify 4×4 within 1e-6**

**Milestones:**
- 4×4 bound within 1e-6 of QMBCertify
- Runtime/RSS ≤ 2× QMBCertify
- Only then attempt larger lattices

---

## Part VI: Architectural Summary

```
┌─────────────────────────────────────────┐
│  User API: polyopt, cs_nctssos, etc.    │  ← unchanged
├─────────────────────────────────────────┤
│  Problem recognition                    │  ← NEW
│  "Is this a periodic Pauli chain?"      │
│  "Is the group cyclic? dihedral?"       │
├──────────────┬──────────────────────────┤
│  Fast path   │  Generic path            │
│  - DFT       │  - SymbolicWedderburn    │
│  - streaming │  - MomentProblem         │
│  - in-place  │  - polynomial matrices   │
├──────────────┴──────────────────────────┤
│  Shared validation layer                │  ← NEW
│  "Fast path == generic path for N≤12"   │
├─────────────────────────────────────────┤
│  JuMP / Mosek                           │  ← unchanged
└─────────────────────────────────────────┘
```

The fast path is not a different algorithm.  It is a more efficient
implementation of the same mathematical relaxation, valid only when structure
is recognized, verified by automated small-N equivalence tests.

The generic path stays as the fallback for every problem QMBCertify cannot
touch: Bell inequalities, fermionic systems, arbitrary graphs, projector
algebras, state polynomials, and anything else the framework already handles.

---

## Part VII: Risk Register

| Risk | Likelihood | Impact | Mitigation |
|:--|:--|:--|:--|
| RDM constraints move bottleneck from construction to solver memory | High | High | Measure per-stage; don't attempt N=100 until L=30 extrapolation is safe |
| SU(2) block decomposition errors give wrong bounds | Medium | Critical | Generate from CG coupling with tests, never hardcode |
| Certification projection breaks PSD | Medium | High | Work in dual; group-average to symmetrize |
| Specialized path silently differs from generic | Medium | Critical | Mandatory small-N equivalence tests before opt-in |
| 2D D4 symmetry requires momentum-star theory, not simple DFT | High | Medium | Defer to Phase 5; separate project with dedicated design |
| N=100 solver memory exceeds available RAM even with all reductions | Medium | High | Phase 3 (SU(2)) must succeed before N=100 attempt |

---

## Part VIII: Priority Summary

| Priority | Phase | What it buys |
|:--|:--|:--|
| 1 | Phase 0 | Falsifiable targets; no ambiguity |
| 2 | Phase 1 | 100× construction speedup; unlocks N=40–100 structurally |
| 3 | Phase 2 | Bound quality parity with QMBCertify |
| 4 | Phase 3 | Solver memory reduction via SU(2); needed for N=100 |
| 5 | Phase 4 | Rigorous certification (independent of scaling) |
| 6 | Phase 5 | 2D lattice (separate project) |

Phases 1–2 together close the gap with QMBCertify's public code.  Phase 3 is
required to match the paper at N=100.  Phase 4 adds scientific rigor.  Phase 5
is a separate research effort.
