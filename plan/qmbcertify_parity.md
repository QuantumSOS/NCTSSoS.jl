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

**Cost attribution note:** The 1121s / 88 GiB at N=64 is NOT solely polynomial
materialization.  The scaling analysis (8.6× per N-doubling, from 8.8s at N=16
to 1121s at N=64) matches O(|G|×|basis|²) Wedderburn decomposition cost, not
O(|basis|²) polynomial construction.  Estimated breakdown: ~80% generic
Wedderburn decomposition, ~20% polynomial matrix construction/linearization.
Small-N no-solver stage attribution confirms the direction: at N=8,10,12 and
d=4, generic charge/spatial/sign construction spends most of its measured time
inside charge/spatial Wedderburn block construction, with finite group closure
secondary and basis/generator setup negligible.  Large-N attribution is still
an extrapolation until the same stages are measured at N=16+ on an isolated
machine.

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

### Lesson 1: Analytic fast paths for known group families

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

**Note:** `pauli_translation_invariant_moment_relaxation` in `pauli_chains.jl`
already implements momentum-sector DFT — it builds blocks with analytic DFT
coefficients for cyclic translation symmetry.  The work is extending this
existing path to accept additional constraints and to handle the full discrete
symmetry group (reflection, sign, axis permutation), not building DFT from
scratch.

### Lesson 2: Streaming accumulation (bypass symbolic matrices)

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

This optimization becomes significant after the DFT fast path eliminates the
Wedderburn bottleneck.  Standalone, it saves ~20% of allocation; combined with
DFT, it completes the transition from O(|G|×|basis|²) generic decomposition to
O(reduced_entries) direct emission.

**Note:** `_ConstraintMatrixEntryCache` in `symmetry.jl` now covers the Pauli
charge lazy path and the nontrivial generic symmetry path for all
`offblock_check` modes.  It caches and streams constraint matrix entries
without building full polynomial matrices, while the `:full` and `:randomized`
certificates still run.  The no-symmetry direct-linear path now covers
ordinary Monoid/Pauli/PBW polynomial problems automatically, and explicit
trivial signed finite actions now route through the no-symmetry branch by
default instead of generic symmetry reduction; the remaining gap is broader
analytic integration with the DFT fast path.  After the perf-harness guard
hardening, the focused
`test/relaxations/symmetry.jl` suite was rerun directly through easy-ssh and
passed 1088/1088 checks in about 2m8s wall time.

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
  generalizes to any SU(2)-invariant Pauli spin system.

- **Don't abandon the type system.**  `NormalMonomial{A,T}` with compile-time
  algebra dispatch is the right design.  The problem is allocation patterns in
  hot loops, not the type hierarchy.

- **Don't fork an external one-problem script.**  QMBCertify is a dead end
  architecturally — it can never grow beyond its one problem.  Instead, build
  *general mechanisms* (streaming accumulation, analytic group decomposition,
  in-place simplification) that happen to make the Heisenberg case fast.
  An internal, tested Pauli-chain fast backend (e.g., extending
  `PauliChainSymmetrySpec` or `pauli_translation_invariant_moment_relaxation`)
  is exactly the right architecture.  The mistake would be copying QMBCertify
  as a standalone module rather than building general mechanisms that the
  Pauli-chain backend uses.

---

## Part IV: Current State of NCTSSoS.jl

### What already exists

| Capability | Status | Location |
|:--|:--|:--|
| Pauli algebra + simplification | ✅ Production | `src/simplification/pauli.jl` |
| Translation-invariant relaxation | ✅ DFT path with RDM/state-opt hooks | `src/optimization/pauli_chains.jl` |
| Contiguous chain basis | ✅ Correct sizes at N=100 | `pauli_contiguous_chain_basis` |
| Charge sector splitting | ✅ Order-2 | `PauliChargeSectorSpec` |
| Singlet constraints | ✅ Order-2 | `PauliSingletConstraintSpec` |
| Pauli sign symmetry | ✅ | `pauli_sign_symmetry` |
| Reflection splitting in translation path | ✅ | `reflection_symmetry=true` |
| Global Pauli-axis rotations | ✅ Scalar equalities plus standalone base PSD axis-isotypic split | `axis_rotation_equalities=true`, `axis_rotation_symmetry=true` |
| Pauli k-RDM positivity | ✅ Full, U(1), and local SU(2) modes | `contiguous_rdm_k` |
| Linear and PSD state optimality | ✅ Small-N tested | `linear_state_opt_width`, `psd_state_opt_width` |
| Structural target accounting | ✅ Analytic N=100 translation/SU(2) targets plus RDM/state-opt block shapes, no large build | `pauli_translation_structural_targets`, `pauli_su2_contiguous_chain_structural_targets`, `pauli_su2_translation_orbit_structural_targets` |
| Pauli-chain fast-path recognition | ✅ Conservative profile + SolverConfig/`cs_nctssos` bridge | `pauli_chain_fast_path_profile`, `pauli_translation_invariant_nctssos(pop, config)`, `cs_nctssos(pop, config; ...)` |
| Certification provenance hooks | ✅ Dual block/value diagnostics for current fast path | `sos_dual_blocks`, `sos_dual_block_values`, `sos_dual_block_diagnostics`, `sos_dual_certificate_residual` |
| Generic Wedderburn symmetry | ✅ But too slow for large N | `src/optimization/symmetry.jl` |
| Fermionic v2RDM | ✅ | `src/optimization/v2rdm_structured.jl` |
| Lazy constraint entry cache | ✅ Pauli charge and nontrivial generic symmetry paths, all offblock_check modes | `_ConstraintMatrixEntryCache` in `symmetry.jl` |

After the focused harness hardening, `test/relaxations/sparsity.jl` was also
made directly runnable and passed 6/6 checks through easy-ssh in about 28.5s,
covering term-sparsity initialization and `PolyOptResult` structural fields.
The combined `test/relaxations/runtests.jl` suite then exposed two dense-GNS
standalone tolerance failures when the suite inherited a non-COSMO solver
fallback.  After giving `gns_pipeline.jl` its own Clarabel standalone fallback
and aligning the affected flatness/entrywise checks with the existing
`5e-5` GNS verification tolerance, the targeted GNS pipeline passed 100/100
checks with both the Mosek wrapper (37.6s) and standalone Clarabel fallback
(41.4s).  The full relaxation component suite then passed 7833/7833 checks
through easy-ssh in 8m50.4s.  The package-level `Pkg.test()` gate then passed
through easy-ssh with 11,672 passing checks and 2 intentional broken checks
(11,674 total) in 21m47.1s.

### What is missing

| Capability | Gap |
|:--|:--|
| Streaming JuMP emission | Supported ordinary no-symmetry Monoid/Pauli/PBW polynomial paths now select direct-linear data automatically; explicit trivial signed finite actions default to the same no-symmetry branch, while broader analytic integration remains pending |
| Automatic analytic DFT dispatch beyond supported Pauli-chain specs | `cs_nctssos` routes supported 1D translation/reflection/sign specs and finite global axis rotations; 2D and broader analytic group families still need reducers |
| Axis-rotation block reducer | Symbolic and direct-linear base moment PSD blocks can be emitted in the axis-isotypic basis, including reflected sectors; sign is subsumed by the axis group, while broader lattice/group families remain pending |
| Reviewed QMBCertify numerical runs | TOML fixtures pin required runs and structural targets; A0_L10 and A1_L20/A1_L30 are reviewed; A0/A2 larger rows stay blocked until an exploratory profile variant earns an accepted reviewed status |
| Full SU(2) moment-matrix reduction | Base translation moment PSD blocks are SU(2)-reduced, including all-sector real reflection splitting, moment-equality rows, closed-support full/U(1)/SU(2) RDM blocks, and PSD state-opt blocks; remaining gap is equality-heavy solver form and larger numerical validation |
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

This path produces the right structure and now accepts RDM/state-opt hooks;
A0_L10 and A1_L20/A1_L30 QMBCertify reference rows are reviewed so far, while
A0/A2 larger rows are blocked by the reviewed-status policy until a deliberate
profile variant earns an accepted status.

Small no-solver comparison from `perf/pauli_translation_compare.jl`
(2026-06-28, easy-ssh, one thread, N < 14 guard active):

| N | d | Generic charge/spatial/sign time | Specialized translation DFT time | Ratio |
|--:|--:|--:|--:|--:|
| 8 | 4 | 4.368377 s | 0.234115 s | 18.659× |
| 10 | 4 | 4.981501 s | 0.499372 s | 9.976× |
| 12 | 4 | 3.462360 s | 0.668133 s | 5.182× |

These are structural construction comparisons only: no SDP solver calls, and
no large-N run.  The monitored wall time for this script run was approximately
55 seconds.

The comparison script now also emits stage attribution.  A small no-solver
rerun at N=8,10,12 (2026-06-28, easy-ssh, one thread, N < 14 guard active)
measured:

| N | Generic total | Generic group closure | Generic charge/Wedderburn blocks | Translation total | Translation block assembly | Translation linearization |
|--:|--:|--:|--:|--:|--:|--:|
| 8 | 6.158898 s | 0.539741 s | 3.793574 s | 0.246363 s | 0.123529 s | 0.085401 s |
| 10 | 6.834882 s | 0.553580 s | 4.427924 s | 0.552811 s | 0.323981 s | 0.171238 s |
| 12 | 5.371578 s | 0.455162 s | 3.021858 s | 0.736382 s | 0.420429 s | 0.239546 s |

The run was no-solver and completed in roughly 61 seconds monitored wall time.

A guarded N=8 rerun after the direct-builder finalization cleanup
(2026-06-29, easy-ssh, one thread, N < 14 guard active) measured
6.476534 s for the generic charge/spatial/sign path and 0.248171 s for the
specialized translation DFT path, a 26.097× construction speedup.  The run was
no-solver and completed in 56.80 seconds monitored wall time.

A guarded N=8 rerun after the standalone axis/add-on composition and
SolverConfig routing updates (2026-06-30, easy-ssh, one thread, N < 14 guard
active) measured 6.587105 s for the generic charge/spatial/sign path and
0.248927 s for the specialized translation DFT path, a 26.462× construction
speedup.  The run was no-solver and completed in about 64 seconds monitored
wall time.

A guarded N=8 target-only axis profile after adding structural accounting for
`axis_rotation_symmetry=true` (2026-06-30, easy-ssh, one thread, N < 14 guard
active) reported 20 axis-isotypic PSD blocks, solver-facing max block 12, and
462 required axis equality rows without constructing a relaxation or calling a
solver.

A guarded N=8,10,12 rerun after the base SU(2) solve-boundary updates
(2026-06-29, easy-ssh, one thread, N < 14 guard active) remained no-solver and
completed in about 83 seconds monitored wall time.  It measured construction
ratios of 22.642×, 11.420×, and 6.513× respectively; the specialized path
reported product-cache hit rates of 0.797342, 0.831025, and 0.855107.

Small guarded no-solver smokes after the constructed RDM metric/provenance
updates (2026-06-29, easy-ssh, one thread, N < 14 guard active) still pass the
diagnostic scripts.  `perf/pauli_translation_compare.jl` at N=8, order=4,
reflection enabled measured 6.741488 s for the generic charge/spatial/sign
path and 0.266486 s for the specialized translation DFT path, a 25.298×
construction ratio.  `perf/pauli_sparse_chain_d4_blocks.jl` at N=8, d=4
reported the generic charge/spatial/sign structural block decomposition in
3.878721 s with max block 30 and no solver call.  The
`perf/pauli_charge_singlet_prep.jl` N=8 no-solver smoke reported
`moment_relax_symmetric` in 3.460524 s, largest singlet-prep block 12, and
690 reduced PSD scalar variables.

---

## Part V: Phased Implementation Plan

### Benchmark targets

Two tiers, never mixed:

**Tier A — QMBCertify code parity:**
- 1D Heisenberg chain, L=30
- Pin: QMBCertify commit, Julia version, Mosek version, hardware (CPU model,
  RAM, threads), BLAS
- QMBCertify options: document exact command, RDM k-values, state-opt options,
  solver tolerances
- Separate stage timings: basis construction, symmetry decomposition, model
  build, solver, total
- Sub-profiles: A0 (base moment only), A1 (+RDM), A2 (+state-opt),
  A3 (+SU(2))
- Pass: `|E_NCTSSoS − E_QMBCertify| ≤ 1e-6`, runtime ≤ 2× QMBCertify,
  RSS ≤ 2× QMBCertify
- Report: block histogram, moment count, PSD dimensions, objective

**Tier B — Paper structural parity:**
- 1D Heisenberg chain, N=100, d=4, sparse contiguous basis
- Define "structural build": basis generation + orbit reps + DFT block
  assembly (excludes JuMP model and solver)
- Distinguish logical complex block size from solver-facing realified block
  size (existing path shows max 62 realified vs 31 logical)
- Pass: basis = 12,001, orbit reps = 121, logical max block = 31,
  structural build ≤ 10s
- Report: full block histogram (not just max), total moment count, memory

Every milestone uses exact numerical targets.  "Competitive" is not allowed.

### Verified QMBCertify source snapshot

Local source checkout:

```text
/Users/exaclior/projects/QMBCertify
commit b18830a9460de4daa03013e389808d522c7823cf
date   2025-07-19 14:14:33 +0800
remote git@github.com:wangjie212/QMBCertify.git
```

Initial constraint inventory from that snapshot:

| QMBCertify option | Source location | Meaning |
|:--|:--|:--|
| `rdm=8`, `rdm=9`, `rdm=10` | `src/rdm_positivity.jl`, `src/sdp.jl` | contiguous RDM positivity on 8, 9, or 10 sites with hardcoded block decompositions |
| `lso=true` | `src/sdp.jl` | linear state optimality rows `ℓ([H,u]) = 0` for filtered generated monomials |
| `pso>0` | `src/sdp.jl`, `PSDstate_entry` | PSD state optimality blocks using the Araujo/Fawzi one-H condition |

Source-level details from commit `b18830a9460de4daa03013e389808d522c7823cf`:

- `GSB(...; lso=true, lol=Int(L/2), pso=3, rdm=false, extra=0, ...)`
  is the default API.  The examples override these defaults case-by-case.
- Chain `rdm` support is exactly `8`, `9`, or `10`; larger values only print
  "Adding rdm > 10 is not supported!" and do not add a PSD block.
- RDM support is not "all local Pauli strings".  Before adding the PSD block,
  `sdp.jl` adds only local words whose counts of X/Y/Z labels are all even,
  then reduces each word with `reduce4`.
- `posepsd8!`, `posepsd9!`, and `posepsd10!` build hardcoded real PSD block
  decompositions from `real(kron(Pauli...))` restricted to the listed index
  blocks.  The resulting block sizes are source data, not discovered by a
  generic representation routine at runtime.
- Chain `lso=true` calls `generate_mons(L, lol, rdm-1)`, filters candidates
  through `filter_mons`, and adds one free scalar variable per surviving
  commutator row.  `filter_mons` rejects candidates whose reduced commutator
  moments are not present in `tsupp` and uses random integer weights only as a
  row-deduplication hash, so reviewed A1/A2 runs must pin actual emitted row
  counts rather than infer them from `lol` alone.
- Chain `pso>0` sets `gb[i] = basis[i][length.(basis[i]) .<= pso]` and emits
  PSD state-opt blocks using `PSDstate_entry(ba1, ba2, h, coe, L)`, where a
  term is added only when exactly one endpoint of a Hamiltonian edge changes
  the candidate word under `count_reduce`.

The RDM block diagnostics `pauli_qmbcertify_rdm_blocks(8)` and
`pauli_qmbcertify_rdm_block_sizes(8)` reproduce the source-pinned `posepsd8!`
storage rows and sizes `[72, 64, 56]` from QMBCertify's even-axis local RDM
graph.  Repeated equivalent components are represented once, matching the
reusable PSD block variables in QMBCertify.  The direct-linear finite-axis
quotient path can now emit these QMBCertify-specific row groups as shared RDM
PSD blocks, where the surrounding moment basis uses the same
translation/reflection/global-axis identifications.  The symbolic constructor
and non-quotiented direct-linear paths still reject this decomposition.
Structural target payloads now aggregate the source-pinned QMBCertify RDM
reference metrics across k=8,9,10 as first-class fields
(`qmbcertify_rdm_reference_block_sizes`, max block, dense/symmetric entries,
bytes, and construction flag), while preserving the per-k provenance records.
The source-derived per-k RDM metric fixture now also pins dense/symmetric byte
counts for k=8,9,10; a red loader check failed on the missing byte fields in
14.9s, and the green loader check passed all 727 checks in 12.9s.
The N=100 full-RDM and QMBCertify-style RDM structural rows now pin aggregate
PSD dense/symmetric byte counts as well; the red loader check failed on the
missing fields in 13.4s, and the green loader check passed all 731 checks in
13.5s.
The target-only profile prints those aggregate metrics next to the ordinary
per-k reference block table, so no-solver structural runs expose the same
budget fields that the expectation fixtures pin.  The explicit translation
constructors and the `cs_nctssos(pop, config; ...)` SolverConfig bridge still
regress the unsafe structural-target-only routes, while the direct-linear
finite-axis quotient bridge now constructs a small `N=4`, order-2, `k=4`
QMBCertify RDM smoke with representative row-block sizes matching
`pauli_qmbcertify_rdm_block_sizes(4; ambient_sites=4)`.  The focused
Pauli-chain test file passed again with 2,210/2,210 translation-invariant
relaxation tests in 2m30.6s.
A high-k no-solver smoke on 2026-06-30 then exercised the source-pinned local
shape at `N=8`, order 4, `contiguous_rdm_k=8`,
`contiguous_rdm_decomposition=:qmbcertify`, reflected finite-axis quotient, and
no solver calls.  The default `contiguous_rdm_support=:closed` guard rejected a
missing moment key, so the intentional high-k profile was rerun with
`contiguous_rdm_support=:extend`.  That run constructed the expected
QMBCertify RDM blocks `[72, 64, 56]`, reported 346 post-quotient linear
moments, 27 PSD blocks, no zero rows, solve support `true`, and 43.424737593s
construction after warmup with peak process HWM about 3.08 GiB.  This is
correctness/provenance progress, not yet the final high-k performance target.
The same profile was rerun after a sampled CPU profile showed that the
QMBCertify branch still built the full 256x256 contiguous RDM block before
discarding it.  Moving the QMBCertify branch ahead of full-block construction
and batching the three source-pinned row groups preserved the same 346 linear
moments, 27 PSD blocks, `[72, 64, 56]` RDM blocks, and solve support, while the
measured construction dropped from a fresh 43.096912671s baseline to
7.021947403s after warmup (about 6.1x faster) and allocation dropped from
3.98 GiB to 2.01 GiB.  A follow-up registration-token cleanup removed the hot
`Tuple(key)` conversion for vector moment keys; the same profile then measured
6.677197971s construction with 1.88 GiB allocated and unchanged model counts.
The next fresh same-case baseline measured 6.394811728s construction with
1.88 GiB allocated.  Seeding the finite-axis quotient rewrite from the
quotient-owned representative key/monomial map, then caching the reusable
QMBCertify RDM row-group templates, reduced the same measured construction to
4.520138600s with 1.42 GiB allocated.  The model counts stayed unchanged:
346 post-quotient linear moments, 27 PSD blocks, `[72, 64, 56]` RDM blocks,
and solve support `true`.  Precomputing the finite-axis quotient map from the
union of base moment support and exact QMBCertify RDM add-on support then
allowed base blocks to emit directly in quotient coordinates, while only the
RDM add-on blocks are quotiented as they are appended.  The same N=8 no-solver
profile measured 2.235024531s construction with 904.40 MiB allocated, 1.8165s
block assembly, 0.3632s precompute, and the same model counts and solve
support.  A small value cache for repeated QMBCertify RDM add-on forms then
reduced the same warmed profile to 725.85 MiB allocated at 2.227067479s, with
unchanged counts and solve support; this is a memory-pressure improvement, not
a measurable warmed-time speedup.  Reusing the restricted axis/reflection
projectors already computed during group selection reduced the warmed profile
again to 2.057139341s with 717.30 MiB allocated; block assembly fell to
1.6429s.  The post-fix runs peaked at about 1.77-1.93 GiB process HWM.
A small public `cs_nctssos` smoke at `N=4`, order 2, `k=4`, direct-linear
finite-axis quotient, and Clarabel also passed, checking the solver-facing
QMBCertify block sizes and solve-support flag.

`perf/qmbcertify_reference_runs.jl` now provides a guarded harness for reviewed
QMBCertify A0/A1/A2 reference runs.  It reads the TOML profile definitions,
defaults to L=10, refuses L>13 unless explicitly overridden, and emits
copyable fixture environment/case rows with commit, Julia/Mosek environment,
objective, block histogram, moment count, parsed build/solve timings, solver
status, and RSS.  The harness rejects non-optimal solves instead of emitting
reviewed fixture rows; rejected candidates now print copyable
`failed_attempts` TOML rows with parsed solver status and stage timings.
Successful QMBCertify/Mosek output can omit explicit `termination status`
lines, so the harness infers `OPTIMAL` from an emitted `optimum = ...` row.
Reviewed rows now require parsed or inferred optimal status and reject both
explicit bad statuses and completely missing status evidence.
Reviewed solve rows also require NCTSSoS commit provenance; on the easy-ssh
workspace this means passing `NCTS_NCTSSOS_COMMIT` explicitly because `.git` is
not synced.  Missing commit provenance is allowed only for bootstrap-only or
explicitly unreviewed probe runs.
On the remote easy-ssh workspace it can bootstrap the pinned QMBCertify checkout
with `NCTS_QMB_BOOTSTRAP=true`; `NCTS_QMB_BOOTSTRAP_ONLY=true` verifies the
checkout/environment without launching a solve.  The harness is include-safe
so small profile probes can reuse its parsers and QMBCertify bootstrap helpers
without launching the default reviewed-run `main()`.  The companion
`perf/qmbcertify_profile_probe.jl` accepts TOML argument overrides through
`NCTS_QMB_ARG_OVERRIDES` or `NCTS_QMB_ARG_OVERRIDES_FILE` and emits
`profile_probes`/`failed_attempts` rows rather than reviewed `cases`.  Failed
probe rows keep the parent `case_id` separate from the probe label, so an
experimental option variant cannot be mistaken for a reviewed profile, and
they include compact effective-argument/override summaries for provenance.

A reviewed A0_L10 reference run is recorded in
`test/data/expectations/heisenberg_qmbcertify_base.toml` from the
2026-06-28 autodl easy-ssh environment.  The run used the pinned QMBCertify
commit, Julia 1.12.6, Mosek 11.2.0, one Julia thread, and recorded objective
`-0.4515446349992884`, moment count 527, max PSD block 84, build time
3.891732289s, solve time 4.184498559s, and peak RSS 1196720128 bytes.
A 2026-06-30 guarded rerun of the same A0_L10 profile on easy-ssh reproduced
the objective, moment count, and block histogram; the harness reported
3.941845042s build time, 4.429378152s solve time, 25.659613652s total wall
time, and peak RSS 1147973632 bytes.
A guarded A0_L20 run on 2026-07-06 used one Julia thread and
`NCTS_QMB_MOSEK_THREADS=1` on autodl.  It did not crowd the machine, finishing
inside the 10-minute guard in about 70s, but MOSEK returned
`termination_status=SLOW_PROGRESS` and `solution_status=FEASIBLE_POINT`.
That row is recorded as a failed attempt, not a reviewed reference.  A0_L30
should wait until an exploratory A0 variant first demonstrates an accepted
status under the explicit reviewed-status policy.

The source-pinned A1 example in `examples/example.jl` is an N=20 run.  The
A1_L10 candidate was an extrapolated small-probe target from that source line,
not a literal QMBCertify example row.  A1_L10 and A2_L10 were reproduced on
2026-06-29 and rejected by the guarded harness because QMBCertify/Mosek returned
`termination_status=DUAL_INFEASIBLE` and
`solution_status=INFEASIBILITY_CERTIFICATE`.  They are not reviewed reference
runs.  The failed solver statuses and stage timings are recorded under
`failed_attempts` in `test/data/expectations/heisenberg_qmbcertify_rdm.toml`;
the intended RDM/state-opt parity profiles need confirmation before recording
A1/A2 data.  A small A1_L10 probe with `lso=false` was also rejected with the
same solver status, so the RDM failure is not explained by QMBCertify's default
linear state-opt rows.

### Degree closure table

The PSD state-opt formula has been checked against QMBCertify source and
the paper text.  It is the Araujo/Fawzi one-H condition
`ℓ(vHw* - 1/2(Hvw* + vw*H))`, not a double-commutator condition.

Before implementing any constraint, verify moment coverage:

| Feature | Required moments | Closure rule (d=4, h=2) |
|:--|:--|:--|
| Base moment PSD | u†v for u,v ∈ B_d | Generated by construction |
| Objective H | Pauli words in H | d ≥ 1 sufficient |
| k-site contiguous RDM | All Pauli strings on k sites | k ≤ 2d (i.e. k ≤ 8) |
| Linear state-opt ⟨[H,A]⟩=0 | Commutator words | test width a: a+h-1 ≤ 2d → a ≤ 7 |
| PSD state-opt | `vHw* - 1/2(Hvw* + vw*H)` words | test width a: 2a+h ≤ 2d → a ≤ 3 |
| SU(2) equalities | No new moments | Existing moments suffice |

Implementation status: `closure_check(feature, basis)` enumerates required
monomials after translation reduction for explicit bases.  The translation
constructors now fail before block assembly using cheap degree bounds for the
standard contiguous basis and exact `closure_check` enumeration for
user-supplied bases.

---

### Phase 0 — Lock references and success criteria

**Duration:** 2–4 days

**Actions:**
- Pin QMBCertify commit, NCTSSoS commit, Julia version, Mosek version, hardware
- Run QMBCertify on reference hardware for L=10,20,30 and record:
  objective, PSD block histogram, moment count, build time, solve time, peak RSS
- Write TOML reference files under `test/data/expectations/`:
  `heisenberg_qmbcertify_base.toml`, `heisenberg_qmbcertify_rdm.toml`
- Extract exact QMBCertify constraint inventory: which RDM k-values, contiguous
  vs all subsets, state-opt formula, A-family, normalization, spin-sector
  assumptions, discrete symmetries, and solver default options

**Milestone:** Reference TOML files committed with exact numerical targets AND
constraint inventory documented.

**Current status:** The pinned checkout bootstrap is implemented in
`perf/qmbcertify_reference_runs.jl`, the small A0_L10 reference case is
recorded and loader-tested, and the source-level RDM/linear-state-opt/PSD
state-opt constraint inventory is documented and pinned in the expectation
fixtures.
The larger A0 plus A2 reference runs remain pending.  The source-derived
A1_L10/A2_L10 RDM candidates are recorded as failed attempts, including an A1
`lso=false` probe, not reviewed references.  A1_L10 is explicitly marked as
extrapolated from the source-pinned N=20 example, pending profile confirmation.
A guarded A1_L12 profile probe on 2026-06-30 then solved successfully under
the same source-pinned A1 command (`rdm=8`, `pso=0`, `lso=true`,
`lol=Int(L/2)`, `extra=9`) with objective/site `-0.4489492424398082`, 1,044
moments, max PSD block 114, 19.301671432s parsed build time, 5.550101014s
solve time, 41.909802836s harness wall time, and 1.30 GiB peak RSS.  It is
recorded as a `profile_probes` row, not a reviewed required case, and indicates
that the A1_L10 failure is likely a short-chain/profile artifact rather than a
universal A1 command failure.  A guarded A2_L12 profile probe also solved under
the source-derived state-opt command (`rdm=8`, `pso=3`, `lso=true`,
`lol=Int(L/2)`, `extra=9`) with objective/site `-0.4489491908788246`, 1,184
moments, max PSD block 114, 20.402644327s parsed build time, 7.782392123s
solve time, 45.123412898s harness wall time, and 1.51 GiB peak RSS.  Both L=12
rows remain probes, not reviewed required cases.  The expectation-loader
fixture checks now cover these rows with exact PSD block histograms,
inferred `OPTIMAL` termination status, objective totals, and environment
provenance.  The harness parser now treats a successful QMBCertify log with an
`optimum = ...` row and no printed non-optimal status as `OPTIMAL`, matching
QMBCertify's source behavior of printing termination/solution status only when
`termination_status(model) != MOI.OPTIMAL`.  The
standalone `test/expectations_loader.jl` file also includes the expectations
helper when run outside `test/runtests.jl`; it now tests the QMBCertify
size guard, probe-length parsing, inline/file override conflict, effective
probe argument override handling, and optimal-status inference.  It passed
459/459 checks through easy-ssh in 6.44s remote wall.
After the large-run guard was lifted with the constraint that any run likely to
stall autodl must be skipped or removed from the plan, the source-pinned A1
profile was run at L=20 and L=30 with one Julia thread and
`NCTS_QMB_MOSEK_THREADS=1`.  Both completed safely and are now reviewed
reference rows: A1_L20 solved in 89.364792735s total wall with 1.57 GiB peak
RSS, objective/site `-0.4452477259458568`, 3,534 moments, and max PSD block
114; A1_L30 solved in 249.631930213s total wall with 2.48 GiB peak RSS,
objective/site `-0.4441690264512694`, 7,405 moments, and max PSD block 114.
The matching A2_L20 run also completed safely inside the 15-minute guard, but
MOSEK returned `SLOW_PROGRESS`/`FEASIBLE_POINT` after 176.571691721s solve
time, so it is recorded as a failed attempt rather than a reviewed A2 row.
A2_L30 should not be launched until an exploratory A2 variant first
demonstrates an accepted status under the explicit reviewed-status policy; the
evidence does not show a machine-stall risk, but it does show the current
profile cannot produce a reviewed row at L=20.
The matching NCTSSoS finite-axis source-like construction probes were then run
with one Julia thread, no JuMP lowering, and no solver call.  A1_L20 completed
in 235.17s remote wall with 4.04 GiB peak RSS, A2_L20 in 239.06s with
3.90 GiB peak RSS, A1_L30 in 503.67s with 6.86 GiB peak RSS, and A2_L30 in
496.51s with 6.58 GiB peak RSS.  These runs show that L=30 construction is
safe on autodl under a hard guard and that block assembly remains the dominant
wall-time cost.  The next no-solve SOS-dual model-sizing pass then ran the
heavier A2 row: A2_L20 built in 242.63s remote wall with 3.58 GiB peak RSS,
153,454 dual variables, 50 PSD cones, and max block 120; A2_L30 built in
517.50s remote wall with 6.93 GiB peak RSS, 226,909 dual variables, 65 PSD
cones, and max block 120.  Those rows justified a smaller solver rung, not an
unbounded large solve.  A1_L20 then solved with Mosek in 377.39s remote wall,
118.795514017s solve time, 4.23 GiB peak RSS, `OPTIMAL` status, objective
`-9.092707010191038`, and SOS residual `3.52e-9`; compared with QMBCertify
A1_L20's reviewed total estimate `-8.904954518917137`, NCTSSoS is lower by
about `0.187752491273901`.  A2_L20 solved similarly in 384.37s remote wall,
127.170103655s solve time, 3.84 GiB peak RSS, `OPTIMAL` status, objective
`-9.092707017660949`, and SOS residual `1.41e-9`.  Finally, the reviewed
QMBCertify-comparable A1_L30 row solved in 1092.76s remote wall with
542.494610598s solve time, 7.09 GiB peak RSS, `OPTIMAL` status, objective
`-13.664267269723961`, and SOS residual `4.49e-10`.  Runtime and memory are
safe, but bound-quality parity is not closed: the reviewed QMBCertify A1_L30
total estimate is `-13.325070793538082`, so NCTSSoS is lower by about
`0.339196476185879`.  The same mismatch direction at L=20 and L=30 points to a
formulation or normalization gap, not an L=30 runtime artifact.  Do not launch
A2_L30 as a reviewed parity row until an exploratory QMBCertify A2 variant
first demonstrates an accepted status under the explicit policy; otherwise it
would only add runtime evidence, not a clean parity comparison.

**Stop gate:** None.  This is prerequisite bookkeeping.

---

### Phase 1A — Analytic DFT / orbit specialization (the big mover)

**Duration:** 1–2 weeks

**Goal:** N=100 d=4 structural decomposition in <10s.

This is the dominant optimization.  The N=64 scaling data shows
O(|G|×|basis|²) Wedderburn cost; replacing it with analytic DFT eliminates
~80% of construction time and allocation.

**Actions:**

1. **Extend existing `pauli_translation_invariant_moment_relaxation` in
   `pauli_chains.jl`**
   - It already builds momentum-sector blocks with DFT coefficients
   - It now accepts contiguous RDM, U(1)/SU(2) RDM, linear state-opt, and PSD
     state-opt hooks; reviewed QMBCertify parity runs are still pending

2. **Build translation orbit representatives from relative support patterns**

3. **Handle k=0, k=N/2, and conjugate-pair realification explicitly**

4. **Add full discrete symmetry profile beyond cyclic translation:**
   sign symmetry, conjugation, reflection, axis permutation — cyclic DFT
   alone likely produces max block ~62 (realified), not 31
   - Scalar axis-rotation moment equalities are now available through
     `axis_rotation_equalities=true`.  The finite global H/S axis-rotation
     generators are now exposed as sparse signed-permutation matrices on the
     translation-orbit basis via
     `pauli_axis_rotation_action_matrices`, together with their full
     24-element group closure and element orders.  Central irrep projectors
     and an orthonormal isotypic transform are exposed through
     `pauli_axis_rotation_irrep_projectors` and
     `pauli_axis_rotation_isotypic_transform`.  The symbolic and direct-linear
     base constructors now apply that transform to standalone base moment PSD
     blocks when `axis_rotation_symmetry=true`, including reflection-adapted
     sectors; first-order linear state-opt scalar rows, contiguous RDM PSD
     blocks, and PSD state-opt blocks compose with this standalone split as
     unsplit add-ons.
     `sign_symmetry=true` is treated as subsumed by the finite axis group
     rather than as a second split.  The constructors also add the required
     axis equality rows.  Equality rows are deduplicated up to scalar sign and
     covered through the symbolic, direct-linear, `SolverConfig`, and
     `cs_nctssos` bridges at small N, including complex HPSD construction mode.
     When `su2_symmetry=true` is explicit, finite axis rotations are accepted
     as redundant because they are subsumed by the SU(2) reducer.

5. **Add conservative fast-path recognition / SolverConfig bridge**
   - `pauli_chain_fast_path_profile` recognizes translation, reflection, sign,
     and global axis-rotation generators; complete global H/S axis-rotation
     specs route through the specialized backend, or as SU(2)-subsumed
     generators when `su2_symmetry=true`
   - `pauli_translation_invariant_nctssos(pop, SolverConfig(...))` routes
     supported Pauli-chain specs without invoking generic Wedderburn
   - `cs_nctssos` routes supported Pauli-chain specs to the translation backend
     and leaves unsupported specs on the generic path

6. **Design the group/decomposition API so 2D direct-product translations
   are not a rewrite later**

7. **Validate:** constructed translation reports match analytic structural
   targets at small N with `N > 2d`; compare objective values and structural
   bounds against the generic path, not bit-for-bit block histograms, because
   the generic real Wedderburn adapter and the direct DFT/HPSD report use
   different block conventions
   - A 2026-06-29 no-solver smoke of
     `perf/pauli_translation_compare.jl` at N=8, order=4, reflection enabled
     measured 6.74s for the generic charge/spatial/sign preparation path
     versus 0.266s for the specialized translation DFT path (~25.3x), with
     the shared-machine N < 14 guard active.

**Milestones:**
- N=100 d=4: basis 12,001, orbit reps 121, logical max block 31,
  structural build ≤ 10s, ≤ 2 GiB
- N=10,12: constructed DFT reports match analytic structural targets for
  translation and reflection variants.  N=16 remains outside the current
  shared-machine test envelope unless explicitly permitted.

**Stop gate:** If max block exceeds 31 after all discrete symmetries,
investigate which symmetry is missing before proceeding.

---

### Phase 1B — Orbit-level product caching + lazy emission

**Duration:** 3–5 days

**Goal:** Eliminate redundant product computations across momentum blocks.

**Current status:** Translation momentum entry construction now uses
`_TranslationProductCache` to cache the k-independent orbit-reduced product
terms for each `(row, col)` pair and apply momentum phases at emission time.
The cache owns a reusable scratch buffer for `simplify!`, and tests verify that
cached monomial keys survive buffer mutation.  Reports and structural targets
surface cache hits, misses, entries, lookups, and hit rate; constructed N=12
d=4 checks show substantial reuse, while the analytic N=100 d=4 target exceeds
the 90% hit-rate milestone.  PSD-state-opt add-on blocks now use the same
principle: unphased translated residual terms are cached once and reused across
momentum sectors.
The direct base-SU2+extending-SU2-RDM path now also reports its real
`_TranslationProductCache` counters instead of zeros, and its target-only
structural payload computes the same cache accounting.  A red N=4/order-1/k=3
smoke failed on missing cache hits in 51.25s; the green smoke matched
constructed and target stats (`hits=12`, `misses=10`, `lookups=22`) in
60.75s; and the focused Pauli-chain regression passed 2,821/2,821 tests in
3m33.7s test time and 386.84s remote wall time.  The corresponding reflected
N=4/order-1/k=3 diagnostic matched the same cache stats in 54.22s.

**Actions:**

1. **Cache orbit-reduced products:** for each orbit pair (i,j), store the
   k-independent part of reduce(b_i† · b_j) → (canonical_key,
   base_coefficient).  Momentum-dependent DFT phases are applied separately.
   - Cache key: (orbit_rep_i, orbit_rep_j)
   - Cache value: (canonical_moment_key, base_coefficient)
   - K-dependent: DFT phase factor e^{i2πk·r/N}, applied at emission time

2. **Buffer reuse:** thread-local reusable `Vector{T}` for `simplify!` to
   avoid allocation in inner loops

3. **Key ownership:** explicitly copy `Vector{T}` keys when storing in linear
   forms (mutable vector keys from scratch buffers will corrupt data otherwise)

**Milestones:**
- Product cache hit rate > 90% for N=100 d=4
- Per-block emission time dominated by DFT coefficient computation, not
  simplification

---

### Phase 1C — Streaming JuMP emission (if profiling justifies)

**Duration:** 3–5 days

**Goal:** Eliminate remaining `Matrix{Polynomial}` intermediate for the fast
path.

**Note:** Only pursue if post-1A/1B profiling shows polynomial materialization
is still a significant fraction of remaining cost.  The lazy path
(`_ConstraintMatrixEntryCache`) may already be sufficient.

**Current status:** `_ConstraintMatrixEntryCache` now has a source-backed
off-block reducer for the checked modes.  The Pauli charge branch in
`moment_relax_symmetric` uses this cache for `:full`, `:randomized`, and `:off`
instead of materializing a full `Matrix{Polynomial}` whenever off-block checks
are enabled.  A targeted red check failed on the missing cache reducer, then the
green no-solver two-site Pauli fixture verified that cache-backed `:full`,
`:randomized`, and `:off` reductions match the materialized path and preserve
the corrupted-basis off-block failures for the checked modes.  A separate
no-solver N=4/order-2 Pauli charge integration check verified that
`offblock_check=:full` matches the default charge-sector report.  These
easy-ssh runs stayed below N=14; the green helper check completed in about
34 seconds, and the charge integration check completed in about 30 seconds.
The same cache reducer was then lifted into the nontrivial generic symmetry
branch.  A second targeted red check failed on the missing generic cache helper;
the green no-solver two-site Pauli symmetry smoke verified cache-vs-materialized
generic reductions and production `moment_relax_symmetric` reports for
`:full`, `:randomized`, and `:off` in about 50 seconds.
As a small checked-mode sanity run, `perf/pauli_charge_singlet_prep.jl` with
`NCTS_PERF_NS=8` and `NCTS_PERF_OFFBLOCK_CHECK=randomized` completed through
easy-ssh in about 63 seconds remote wall, with no solver call.  The reported
N=8 `moment_relax_symmetric` phase was 3.493262s, allocated 621.45 MiB, and
produced 33 PSD blocks with largest block 12.

A `MomentLinearBuilder` shell exists for staged key/form
ownership and finalization into `MomentLinearData`.  Fermionic parity rows and
moment-equality rows now append directly into
`MomentLinearData.zero_constraints` when the PSD cache is already synchronized.
PSD-block JuMP lowering now consumes cached linear zero rows instead of
relinearizing symbolic `:Zero` constraints, and finalized `MomentLinearData`
can be lowered directly to PSD-block or moment-variable JuMP models without a
symbolic `MomentProblem` wrapper.  Finalized `MomentLinearData` can also be
dualized directly into an SOS model without a symbolic wrapper; regression
coverage includes a finalized fermionic HPSD block lowered and dualized without
first building a symbolic `MomentProblem`.
`pauli_translation_invariant_nctssos(...; direct_linear=true, ...)` and
supported `cs_nctssos(...; direct_linear=true, contiguous_rdm_k=...,
linear_state_opt_width=..., psd_state_opt_width=...)` fast-path calls now solve
the directly finalized linear cache and return that `MomentLinearData` in the
result.  The `cs_nctssos` bridge also forwards the Pauli fast-path options for
moment sectors, real-vs-complex moment matrices, contiguous RDM
decomposition/support, U(1)/SU(2) validation, and state-opt hooks; the explicit
`pauli_translation_invariant_nctssos(pop, SolverConfig(...); ...)` overload is
covered by the same direct-linear/RDM/state-opt regression path.  The ordinary
no-symmetry `cs_nctssos` route now selects direct linear data automatically for
supported Monoid/Pauli/PBW polynomial problems, and an explicit trivial signed
finite action, `SymmetrySpec(SignedPermutation())`, is treated as no symmetry
instead of sending it through the generic symmetry pipeline or rejecting it.  A
RED public-interface regression failed on the old blanket
`direct_linear=true` with `solver_config.symmetry` guard in 0:38.39 wall with
1,500,924 KiB peak RSS and no swaps.  The GREEN interface run passed in
2:33.96 wall with 2,103,872 KiB peak RSS and no swaps, matching the
2.5-3.5 minute / <2.5 GiB estimate.  The default route was then tightened too:
a RED public-interface regression showed that
`cs_nctssos(pop, SolverConfig(..., symmetry=SymmetrySpec(SignedPermutation())))`
still returned a symmetry report even though objective, block sizes, and moment
counts already matched the plain no-symmetry route.  That RED run failed in
0:46.34 wall with 1,413,592 KiB peak RSS and no swaps.  The GREEN run routed
the same trivial finite action through the plain no-symmetry branch and passed
in 2:38.29 wall with 2,138,344 KiB peak RSS and no swaps.  This narrows the
trivial-action streaming/materialization gap without changing nontrivial
symmetry behavior.  A route-probe RED regression then showed the ordinary
no-symmetry default still built a symbolic `MomentProblem`; it failed in 0:29.50
wall with 1,441,264 KiB peak RSS and no swaps.  The GREEN interface suite passed
after automatic direct-linear routing in 2:32.35 wall with 2,077,988 KiB peak
RSS and no swaps, matching the 2.5-3.5 minute / <2.5 GiB estimate.  No large
construction, L>=14 solve, or scaling run was used for this evidence.  The same
policy was then extended to `cs_nctssos_higher`: a RED route-probe regression
showed the higher-step path still handed a symbolic `MomentProblem` to
`solve_sdp`, failing in 0:29.49 wall with 1,438,608 KiB peak RSS and no swaps.
The GREEN interface suite passed after routing supported ordinary no-symmetry
higher steps through `_moment_relax_linear`, in 2:31.12 wall with 2,099,940 KiB
peak RSS and no swaps.  The same higher-step route now also treats an explicit
trivial signed finite action as no symmetry: a RED route-probe regression failed
on the old blanket `solver_config.symmetry` rejection in 0:30.12 wall with
1,416,740 KiB peak RSS and no swaps, and the GREEN interface suite passed in
2:32.04 wall with 2,086,048 KiB peak RSS and no swaps.
An API-consistency RED then showed `cs_nctssos_higher(...; direct_linear=true)`
still rejected the explicit keyword, failing in 0:30.70 wall with 1,430,184 KiB
peak RSS and no swaps; after threading the keyword through the same no-symmetry
direct-linear gate, the GREEN interface suite passed in 2:33.36 wall with
2,292,108 KiB peak RSS and no swaps.
The Pauli-chain translation backend now has a direct linear-entry emitter for base
momentum blocks, fixed-momentum reflection
blocks, and realified conjugate-pair reflection blocks.  Full contiguous RDM
blocks also emit directly for both closed and extending support modes,
including complex `real_moment_matrix=false` HPSD mode; U(1)- and
SU(2)-decomposed RDM blocks emit directly with their zero-row provenance in
both real and complex moment-matrix modes.
Translation-compatible scalar equality rows, scalar inequality blocks, and
moment-equality rows now emit directly through the builder as well, including
the complex `real_moment_matrix=false` form.
Linear state-opt rows emit directly as scalar zero constraints in both real and
complex moment-matrix modes, and PSD state-opt blocks emit directly as linear
PSD/HPSD blocks.  Small public solves cover direct-linear linear/PSD state-opt
parity against the symbolic path, including SOS dual residual/block
diagnostics, and a focused N=4 moment-linear regression covers direct complex
HPSD PSD-state-opt parity with native undoubled state-opt block sizes.
Target-only structural accounting is also covered for complex PSD-state-opt
blocks: the N=8, order-2 target keeps eight native 3x3 HPSD state-opt blocks
and counts Hermitian storage with dense `n^2` entries; the same target is
checked against the constructed complex report.  The focused Pauli-chain test
file passed again with 2,059/2,059 translation-invariant relaxation tests in
2m46.4s.
Complex linear-state-opt target accounting is now checked the same way: the
N=8, order-2 target's row count, zero-feature histogram, block sizes, and
HPSD storage counts are compared with the constructed complex report.  The
focused Pauli-chain test file passed again with 2,066/2,066
translation-invariant relaxation tests in 2m44.5s.
Complex `H^2` moment-equality target accounting is likewise checked against a
constructed N=8, order-2 complex report, including row count, zero-feature
histogram, block sizes, and HPSD storage counts.  The focused Pauli-chain test
file passed again with 2,074/2,074 translation-invariant relaxation tests in
2m46.0s.
Complex axis-rotation equality target accounting is now checked against a
constructed N=8, order-2 complex report, including equality-row count,
quotient statistics, zero-feature histogram, block sizes, and HPSD storage
counts.  The focused Pauli-chain test file passed again with 2,087/2,087
translation-invariant relaxation tests in 2m47.1s.
Combined complex add-on target accounting is now checked against a constructed
N=8, order-2 HPSD report with U(1) RDM zero rows, axis-rotation equalities,
linear state-opt rows, and `H^2` moment-equality rows all enabled.  That
regression exposed and fixed a structural-target undercount for complex U(1)
RDM zero rows: the direct-linear report stores real and imaginary scalar zero
rows even when the PSD block itself remains native HPSD.  The focused
Pauli-chain test file passed again with 2,098/2,098 translation-invariant
relaxation tests in 2m45.4s.
A narrower single-feature regression now checks the same complex U(1) RDM
accounting directly: the N=8, order-2 HPSD target has native block sizes
`[1, 3, 3, 1]`, 88 magnetization-offblock zero rows, and HPSD storage counts
matching the constructed report.  The focused Pauli-chain test file passed
again with 2,110/2,110 translation-invariant relaxation tests in 2m45.5s.
Complex SU(2)-decomposed RDM target accounting is now checked at the add-on
block boundary as well: the N=8, order-2 HPSD target reports native local
SU(2) RDM block sizes `[1, 1]`, while the constructed base
translation/SU(2) report confirms those appended RDM block sizes and still
exposes its separate incomplete SU(2) zero-row accounting.  The focused
Pauli-chain test file passed again with 2,123/2,123 translation-invariant
relaxation tests in 2m44.8s.
Complex full-RDM target accounting is now checked directly too: the N=8,
order-2 HPSD target keeps the full contiguous two-site RDM as one native 4x4
HPSD block, adds no scalar zero rows, and counts Hermitian storage with dense
`n^2` entries matching the constructed report.  The focused Pauli-chain test
file passed again with 2,135/2,135 translation-invariant relaxation tests in
2m46.7s.
Complex direct-linear extending-support full-RDM emission is now covered as
well: the N=8, order-2, k=5 HPSD path appends one native 32x32 contiguous RDM
block, adds no scalar RDM rows, matches structural target storage accounting,
and satisfies finalized `MomentLinearData` invariants.  Because this support
mode intentionally introduces moment keys outside the closed base PSD pivot
set, complex PSD-block lowering now has an explicit regression that the default
orphan policy rejects those keys while `orphan_policy=:free_variables` builds
the HPSD JuMP model.  The same contract is now covered through the public
`pauli_translation_invariant_nctssos(...; direct_linear=true,
real_moment_matrix=false, formulation=:psd_blocks, representation=:complex)`
wrapper and the `cs_nctssos(pop, config; ...)` SolverConfig bridge on N=4,
order-1, k=3 extending-support solves: the default orphan policy rejects free
support moments with an error that points to the explicit orphan-policy
options, while `orphan_policy=:free_variables` builds and solves the
direct-linear HPSD model; the explicit wrapper also matches the symbolic path
with the same options.  The focused Pauli-chain test file passed again with
2,188/2,188 translation-invariant relaxation tests in 2m54.5s.
Direct-linear complex HPSD U(1)- and SU(2)-decomposed RDM modes are also
covered through both the explicit translation API and the `cs_nctssos`
SolverConfig bridge.  Complex HPSD mode now rejects reflected
Pauli-chain requests with non-reflection-fixed momentum sectors early, while
still allowing fixed sectors such as `momenta=[0]`; structural targets follow
the same rule and report undoubled HPSD block sizes for those fixed sectors.
HPSD structural metrics now count `n^2` real scalar storage entries for an
`n×n` Hermitian block instead of using the real-symmetric triangular count.
Target-only structural accounting now also includes contiguous RDM block
shapes, PSD state-opt block shapes, and linear/PSD state-opt candidate counts;
linear state-opt row counts are exact for the default nearest-neighbor
Heisenberg target.  Constructed translation report metrics now expose linear
state-opt row counts and
contiguous RDM plus PSD state-opt constructed block sizes/counts, including
constructed RDM zero-row counts, and the no-solver profile prints those
metrics for constructed runs.
The internal builder-backed base relaxation matches the symbolic path at small
N; no-solver
N=12, d=4 direct-builder smokes produced 9,547 moments, 28 PSD blocks, max
block 62 without reflection, 9,547 moments, 56 PSD blocks, max block 44 with
reflection, and 9,547 base moments extending to 9,556 linear moments with a
2-site RDM block.  No-solver N=8, d=3 smokes produced 2,458 moments for
U(1)-decomposed RDM (88 zero rows, 9 PSD blocks, max block 80), linear
state-opt (9 zero rows, 5 PSD blocks, max block 80), and PSD state-opt (10 PSD
blocks, max block 80).  A no-solver N=8, d=2 SU(2)-RDM smoke produced 379
moments, 24 zero rows, 7 PSD blocks, and max block 26.
`perf/pauli_translation_profile.jl` now has
`NCTS_TRANSLATION_DIRECT_LINEAR=true` to profile the direct linear constructor
without solving; a remote N=4 smoke confirmed the symbolic and direct-linear
profile modes report matching structural histograms.
A small reflected order-4 no-solver profile at N=8,10,12 initially exposed an
O(repeated terms × registered moments) registration bottleneck in the
direct-linear constructor: the builder re-registered duplicate vector-valued
moment keys through a linear scan.  A per-construction registered-key token
cache and sparse reflected-row transforms reduced the same N=8,10,12 profile
to ~0.44s, 1.90s, 2.66s outer construction time, from ~4.91s, 41.65s,
67.38s.  For N=12, `block_assembly` is now ~1.64s and builder finalization is
~0.64s after avoiding redundant normalization/copy passes for staged linear
PSD blocks.
A fresh reflected direct-linear no-solver profile on 2026-06-29 (easy-ssh,
one thread, N < 14 guard active) measured 0.513s, 1.923s, and 2.690s outer
construction time at N=8,10,12.  The N=12 stage split was 1.647s block
assembly and 0.652s linearization, so the Phase 1C conclusion is unchanged:
direct-linear construction is useful for add-on constraints and provenance,
but should remain opt-in for base moment blocks.
The same small reflected base profile still favors the symbolic cached path
(~0.25s, 1.02s, 1.05s at N=8,10,12), so direct-linear construction should stay
opt-in for base moment blocks unless extra RDM/state-opt constraints change the
tradeoff.
The ordinary no-symmetry path now has an internal direct-linear construction
helper, `_moment_relax_linear`, that emits `MomentLinearData` without first
materializing symbolic constraint matrices for the covered real ordinary
polynomial cases, guarded fermionic parity cases, complex equality rows, and
complex one-sided moment-equality rows.  A red `moment_linear.jl` regression
failed on the missing helper in 23.1s; the green focused file passed after
implementation in about 114.5s, including parity against symbolic
`moment_relax` for a small no-symmetry relaxation with equality and inequality
constraints.
The helper is now exposed through the public no-symmetry `cs_nctssos` path for
ordinary `MonoidAlgebra` polynomial problems via `direct_linear=true`, and
later became the automatic default for supported ordinary no-symmetry
Monoid/Pauli/PBW polynomial paths and their `cs_nctssos_higher` continuation
steps.  State/trace, unsupported algebra families, generic-symmetry calls, and
Pauli-translation-specific add-on keywords remain explicitly fenced.  The
public regression first failed in 26.9s on the stale Pauli-fast-path keyword
classification, then `test/relaxations/interface.jl` passed in about 130s after
routing the no-symmetry branch to `_moment_relax_linear`.
The direct no-symmetry helper now also emits one-sided moment-equality rows
directly as scalar zero constraints, reusing the symbolic row-basis truncation
without staging `1x1 Matrix{Polynomial}` intermediates.  The focused red test
failed in 24.9s on the old explicit moment-equality guard; the green
`test/relaxations/moment_linear.jl` run passed in about 114.5s, including
parity between direct `MomentLinearData` and symbolic `MomentProblem.linear`.
It also emits fermionic parity zero rows directly before moment-equality rows.
The focused red test failed in about 30s on the old parity guard; after matching
the symbolic Hermitian real/skew split and zero-constraint indexing convention,
`test/relaxations/moment_linear.jl` passed in about 114.8s.  The public
no-symmetry `direct_linear=true` gate now admits guarded `FermionicAlgebra`
problems as well: the public red test failed in about 45.7s at the old
`MonoidAlgebra`-only gate, and `test/relaxations/interface.jl` passed in about
126.9s after the support check was widened.
Complex equality rows now split directly into the same matrix-wide Hermitian
and skew-Hermitian zero components as the symbolic constructor.  A focused red
test failed in about 37.6s on the old complex zero-cone guard; after replacing
the scalar-only zero append with a component-aware direct emitter,
`test/relaxations/moment_linear.jl` passed in about 115.4s.
Complex one-sided moment-equality rows now split directly as well.  The
focused red test failed in about 40.1s on the old complex moment-equality
guard; the first green attempt exposed the scalar split sign convention in
about 50s, and the corrected `test/relaxations/moment_linear.jl` run passed in
about 113.9s.
The public no-symmetry `direct_linear=true` gate now also admits guarded
`PauliAlgebra` problems, including non-Hermitian equality rows.  The public red
test failed in about 61s at the old gate; after widening the support check,
`test/relaxations/interface.jl` passed in about 127.8s while Pauli
translation-specific options remained fenced.
The public no-symmetry gate now admits `PBWAlgebra` problems as well.  A tiny
Bosonic number-operator red test failed in about 55.8s at the old support gate;
after widening the gate to PBW algebras, `test/relaxations/interface.jl` passed
in about 131.4s.  The focused `test/relaxations/moment_linear.jl` validation
also passed, including Bosonic PBW direct `MomentLinearData` parity with the
symbolic path.
After the perf-harness guard hardening, the public interface suite was rerun
through easy-ssh and passed again in about 128s, covering ordinary, fermionic,
bosonic, and Pauli no-symmetry `direct_linear=true` routes plus native
Hermitian SOS dualization.
The focused `test/relaxations/moment_linear.jl` suite was rerun through
easy-ssh after the same hardening and passed 306/306 checks with about 119s of
reported test time.  That run covers the direct no-symmetry builder, cached
linear zero rows, PSD-block and moment-variable lowering, direct SOS
dualization, Pauli translation linear entries, full/U(1)/SU(2) RDM add-ons,
linear and PSD state-opt blocks, and post-mutation cache attachment.
Real reflected PSD blocks are symmetrized at construction time before
moment-variable lowering so floating DFT/reflection noise does not produce
spurious asymmetric JuMP matrices.  For Pauli-chain translation fast paths,
`direct_linear=true` remains opt-in until it is consistently faster than the
symbolic translation constructor across the supported option matrix.

**Actions:**

1. **Define a `MomentLinearBuilder`** that accumulates raw (key, coefficient)
   pairs during basis enumeration, then finalizes: close adjoints, sort/dedup
   keys, discover pivots, construct immutable `MomentLinearData`
   - `MomentLinearData` is an invariant-checked final object — you cannot
     stream into it directly.  The builder handles the two-phase process:
     accumulate → finalize.

2. **Validate:** streamed `MomentLinearData` matches symbolic path for all
   algebra types (real, complex, PBW, fermionic) at small N

3. **Port post-construction mutations** (`_add_parity_constraints!`,
   `_add_moment_eq_constraints!`) to emit `ScalarLinearConstraint` directly
   instead of going through polynomial intermediates
   - `_add_parity_constraints!` and `_add_moment_eq_constraints!` now emit
     direct linear zero rows when the existing PSD cache is synchronized.
   - `_moment_relax_linear` now emits no-symmetry moment-equality rows directly
     during construction for real ordinary polynomial problems.
   - `_moment_relax_linear` now emits fermionic parity rows directly, including
     the symbolic path's Hermitian real/skew zero-row split.
   - Public no-symmetry `cs_nctssos(...; direct_linear=true)` now accepts the
     same guarded fermionic path.
   - `_moment_relax_linear` now splits complex equality matrices directly,
     preserving the symbolic path's constraint-index and origin conventions.
   - `_moment_relax_linear` now splits complex one-sided moment-equality rows
     directly, preserving the scalar zero-row sign convention used by the
     symbolic path.
   - Public no-symmetry `cs_nctssos(...; direct_linear=true)` now accepts
     guarded Pauli complex problems in addition to Monoid and Fermionic cases.
   - Public no-symmetry `cs_nctssos(...; direct_linear=true)` now accepts
     PBW-algebra problems, with a Bosonic smoke matching the symbolic path.
   - `_moment_relax_linear` now matches symbolic `moment_relax` for a Bosonic
     PBW number-operator case.

4. **Maintain full n×n linear-data convention internally** — SOS dualization
   iterates all (i,j) with factor-of-2-sensitive Hermitian logic.  Real
   JuMP PSD-block lowering may bind only the upper triangle at the final model
   boundary.

**Milestones:**
- N=8,10,12: streaming path optimum matches generic path within 1e-8 for
  all supported algebra types
- L=30 base solve matches QMBCertify A0 within 1e-6

---

### Phase 2A — Contiguous k=2 RDM positivity

**Duration:** 3–5 days

**Goal:** Add the simplest RDM constraint and validate infrastructure.

**Actions:**

1. **Implement 2-site contiguous RDM:** ρ_S = 2^{-2} Σ_P y_P P for all
   Pauli strings on 2 adjacent sites

2. **Decompose under U(1) magnetization** (Hamming weight sectors), NOT
   CG/SU(2) — that belongs in Phase 3

3. **Full 4×4 PSD first,** then U(1) block-diagonal (only valid if model
   enforces global U(1) symmetry)

4. **Implement `closure_check(k=2, basis)`** that verifies all required
   moment keys exist

5. **Extend the Phase 1A DFT path** to accept the additional PSD blocks

6. **For QMBCertify-style high-k RDMs, require an explicit support policy:**
   keep `contiguous_rdm_support=:closed` as the default guard, and use
   `:extend` only when the RDM block intentionally introduces new moment keys
   outside the base moment PSD support.
   - Structural targets now accept `contiguous_rdm_decomposition=:qmbcertify`
     with either support policy as source-pinned shape metadata.  Actual
     construction remains guarded by the direct-linear finite-axis quotient
     path, which owns the added moment keys.

**Milestones:**
- `closure_check` rejects k > 2d
- L=20 + k=2 RDM gives tighter bound than base moment relaxation

**Small-N note:** N=6 and N=8, order-2 Clarabel probes under the shared-machine
guard were solver-tolerance dominated (`ALMOST_OPTIMAL`) and did not give a
clean monotone tightening signal.  Do not use those small probes as Phase 2
quality evidence; the L=20/L=30 gates remain the meaningful checks.
A guarded direct-linear no-solver regression now covers intentional
`contiguous_rdm_support=:extend` moment registration for N=8, order=2,
`contiguous_rdm_k=5`.  The probe reports 379 base moments expanding to 811
linear moments with a 64x64 RDM PSD block, and
`test/relaxations/pauli_chains.jl` verifies the same direct-linear invariants.
The focused Pauli-chain test file passed again through easy-ssh with
2,008/2,008 translation-invariant relaxation tests in 2m44.1s.  This is
structural/builder coverage only, not Phase 2 bound-quality evidence.
The direct-linear builder path is also now covered for the constructed
U(1)-decomposed k=3 and k=4 contiguous RDM add-ons at N=8, order=2.  The
guarded no-solver probe reports k=3 blocks `[1, 6, 6, 1]` with 88
magnetization-offblock zero rows and k=4 blocks `[1, 8, 12, 8, 1]` with 372
zero rows; both have 379 linear moments and pass linear-data invariants.  The
focused Pauli-chain test file passed again with 2,028/2,028
translation-invariant relaxation tests in 2m45.9s.
The complex-HPSD direct-linear builder path is now covered for the same N=8,
order-2 U(1)-decomposed k=3 and k=4 RDM add-ons.  The guarded no-solver probe
keeps the native Hermitian block sizes (`[1, 3, 3, 1]` and `[1, 4, 6, 4, 1]`)
instead of realifying them to doubled PSD blocks, while preserving the same 88
and 372 magnetization-offblock zero rows and 379 linear moments.  The focused
Pauli-chain test file passed again with 2,052/2,052 translation-invariant
relaxation tests in 2m44.7s.
Target-only structural reports now accept QMBCertify source-pinned RDM block
shapes with `contiguous_rdm_support=:extend`; the target remains shape-only
and solve-blocked, while construction is owned by the direct-linear
finite-axis quotient path.  The focused Pauli-chain test file passed again
through a temporary test environment with 2,551/2,551
translation-invariant relaxation tests in 2m58.9s.  The QMBCertify RDM
expectation fixture now pins the same `:extend` support policy, and the
expectation-loader tests passed with 392/392 checks.
A small public `N=4`, order-2 solve now also covers QMBCertify RDM blocks
combined with `psd_state_opt_width=1` on the direct-linear finite-axis
quotient path; the report stays solve-supported, appends the expected
QMBCertify RDM block sizes, and includes PSD-state-opt blocks.  In the same
tiny finite-axis quotient setting, width-1 through width-3 linear state-opt
rows reduce to zero both with and without the QMBCertify RDM add-on, so those
small probes are not useful evidence for Phase 2B bound-quality parity.
Small construction-only probes at N=6 and N=8 show the same quotient behavior:
finite-axis quotient reports zero linear state-opt rows while preserving
solver-supported QMBCertify RDM blocks; at N=8, order=4, the plain quotient,
full-RDM quotient, and QMBCertify-RDM quotient all produce 300 quotient moment
keys and 850 forced-zero classes.
The direct-linear finite-axis quotient path now also constructs the
source-like small composition `N=8`, order 4, `contiguous_rdm_k=8`,
`contiguous_rdm_decomposition=:qmbcertify`, `contiguous_rdm_support=:extend`,
and `psd_state_opt_width=3` without a solver call.  That guarded probe reports
349 quotient moment keys, QMBCertify RDM blocks `[72, 64, 56]`, five PSD
state-opt blocks of size 18, zero scalar rows, linear-data invariants passing,
and solve support `true`.  Extending-support RDM profiles therefore also admit
new PSD-state-opt moment keys when the combined add-on profile needs them.
The same support-policy audit showed that `psd_state_opt_width=3` is not
closed by the base contiguous moment support at N=8/order 4: both the no-RDM
case and the closed-support full-RDM case reject the disconnected moment
`UInt8[0x01, 0x04, 0x0a, 0x10]`, while the full-RDM extending-support case
constructs successfully with 367 moments and five PSD-state-opt blocks of
size 18.  Shape-only fixtures now mark the full-RDM-plus-PSD-state target as
extending support instead of closed support.
An N=10 target-vs-construction audit then pinned the finite-axis quotient
storage convention: with reflection, finite-axis quotient, QMBCertify `k=8`
RDM, extending support, and `psd_state_opt_width=3`, the constructed report has
35 PSD blocks, 1,481 quotient moments, max block 120, RDM blocks
`[72, 64, 56]`, and six PSD-state-opt blocks all of size 18.  Structural
targets now expose `axis_rotation_quotient=true` so target-only add-on block
sizes follow that fully-realified quotient convention.
Precomputing the finite-axis quotient map from the union of base support,
QMBCertify RDM support, and translated PSD-state-opt add-on support then
allowed the same source-like `N=8`, order-4, `k=8`, reflected quotient profile
to keep base emission in quotient coordinates while quotienting add-on blocks
as they are appended.  The guarded no-solver construction measured
35.275823678s, preserved 349 quotient moments, 6,994 raw moment keys, 850
forced-zero classes, QMBCertify RDM blocks `[72, 64, 56]`, five PSD-state-opt
blocks of size 18, and solve support `true`.  The PSD-state-opt support helper
now mirrors the actual translated block emitter, so the N=8/order-4
closed-support check rejects an order-1 PSD-state-opt basis and reports 1,252
required width-2 moment keys instead of undercounting the support.  A focused
easy-ssh regression rerun passed the expectation fixtures (393/393) and the
Pauli-chain relaxation file (2,581/2,581) in 335.68s remote wall time.
The PSD-state-opt translated-term cache now stores the unphased
row/column/translation residual terms once and phases them per momentum
sector.  A narrow N=8 width-3 microprofile measured old PSD-state-opt block
emission at 0.802930416s versus 0.175521885s cache build plus 0.050991021s
cached emission for the five real momentum sectors.  The source-like N=8,
order-4, reflected finite-axis quotient profile then measured 1.857817185s
for the warmed construction after a compile-heavy first call, preserving 349
quotient moments, 6,994 raw moment keys, 850 forced-zero classes, 32 PSD
blocks, max block 120, QMBCertify RDM blocks `[72, 64, 56]`, five
PSD-state-opt blocks of size 18, and solve support `true`.  The focused
easy-ssh regression rerun passed the expectation fixtures (393/393) and the
Pauli-chain relaxation file (2,583/2,583) in 350.58s remote wall time.
The same source-like construction profile was then run at N=12, still below
the shared-machine N < 14 performance guard and still without an SDP solve.
It constructed in 37.930539342s, with 46.74s remote wall time, 2,158 quotient
moments, 41,122 raw moment keys, 4,773 forced-zero classes, 38 PSD blocks, max
block 120, QMBCertify RDM blocks `[72, 64, 56]`, seven PSD-state-opt blocks of
size 18, product-cache hit rate 0.8551068883610451, and solve support `true`.
This is construction/shape evidence for the same small L=12 family as the new
QMBCertify A1/A2 probes, not a bound-quality comparison.
After real PSD-block lowering and triangular real binding, the same N=12
source-like profile was lowered build-only, still with no solver call.  The
probe reported 2,158 moments, 2,157 free keys, 38 PSD blocks, max block 120,
QMBCertify RDM blocks `[72, 64, 56]`, seven PSD-state-opt blocks of size 18,
96,037 JuMP variables, and 93,880 scalar binding equalities.  Construction
took 8.638280823s, lowering 1.452863180s, remote wall 61.63s, peak HWM about
1.96 GiB, and the product-cache hit rate remained 0.8551068883610451.  This is
still build-only shape evidence; the free-key count makes a shared-machine
solver run inappropriate without further formulation work.
A solver-backed N=8 smoke of the same source-like profile with Clarabel's
moment-variable path was deliberately stopped after about five minutes when
process RSS reached roughly 43 GiB inside Clarabel/QDLDL factorization.  The
remote status marker was cleared after the kill.  This is not a bound-quality
result; it is negative solver-memory evidence and reinforces the Phase 2/3
warning that construction support does not imply this profile should be solved
with the generic real moment-variable formulation on shared hardware.
A follow-up probe with `perf/pauli_translation_solver_probe.jl` checked whether
the complex PSD-block lowering could provide a safer small solve.  The exact
QMBCertify RDM decomposition cannot take that route today: it requires
`axis_rotation_quotient=true`, while the finite-axis quotient currently
requires `real_moment_matrix=true`, and complex PSD-block lowering requires
complex Hermitian moment data.  A non-QMBCertify full-RDM control at N=6,
order 4, `contiguous_rdm_k=6`, `contiguous_rdm_support=:extend`, and
`psd_state_opt_width=3` did lower successfully with `real_moment_matrix=false`
and `formulation=:psd_blocks`: construction took 1.486793528s, lowering took
2.207393789s, remote wall was 54.92s, peak HWM was about 1.93 GiB, and the
model had 700 moments, 699 free keys, 31 Hermitian PSD blocks
(`9=>6, 30=>23, 31=>1, 64=>1`), 27,641 JuMP variables, and 26,242 complex
binding equalities.  Mosek was not available in the NCTSSoS project
environment for the same probe.  A Clarabel solve of that lowered N=6 control
model was stopped after roughly two minutes when RSS reached about 16.9 GiB
without a result.  This shows that PSD-block lowering fixes model construction
for the complex control case, but does not yet provide a safe solver-backed
bound-quality smoke for the source-like profile on shared hardware.
Real PSD-block lowering for `MomentLinearData` was then added so the
finite-axis quotient path can at least avoid the real moment-variable lowering
boundary.  The focused `test/relaxations/moment_linear.jl` file passed through
easy-ssh in 96.02s, covering the new real PSD-block route and free-orphan
handling.  The exact N=8 source-like finite-axis quotient profile
(`k=8`, `:qmbcertify`, `:extend`, `psd_state_opt_width=3`) now constructs and
lowers build-only with `formulation=:psd_blocks`, `representation=:real`, and
`orphan_policy=:free_variables`: construction took 6.448200598s, lowering took
1.365188000s, remote wall was 58.39s, peak HWM was about 1.93 GiB, and the
lowered model had 349 moments, 348 free keys, 32 PSD blocks, 64,846 JuMP
variables, and 127,514 scalar binding equalities.  This is useful construction
progress, but the large free-key/binding count means a solver-backed run is
still not safe to launch blindly on the shared machine.
Real PSD-block lowering now binds only the upper triangle of symmetric PSD
variables.  The focused red/green regression in `test/relaxations/lowering.jl`
first caught the old duplicate mirrored binding (`3` scalar equalities instead
of `2`) and then passed; the full lowering file passed in 44.48s, and
`test/relaxations/moment_linear.jl` passed in 95.93s.  The same N=8
source-like build-only probe now lowers in 1.236605337s with 243.60 MiB
allocated and 64,498 scalar binding equalities, down from 1.365188000s,
371.76 MiB, and 127,514 equalities.  The model still has 64,846 JuMP
variables and 348 free keys, so this halves the binding rows but does not yet
make the source-like solve safe on shared hardware.
An exploratory scaled-pivot audit checked whether those free keys are mostly
singleton block entries with non-unit scalar coefficients.  They are not:
in the N=8 source-like profile only 9 of 348 free keys appear as singleton
scaled-pivot candidates, across 130 singleton free-key entries.  Generalized
scalar pivots would therefore add lowering complexity for negligible reduction
in the source-like formulation size, so that route is deferred.
A tiny end-to-end N=4, order-2 finite-axis quotient smoke with QMBCertify
`k=4` RDM, `:extend` support, and `psd_state_opt_width=1` still solves through
the new real PSD-block lowering after the triangular-binding reduction.  It
lowers to 148 JuMP variables, 142 scalar binding equalities, and 21 PSD blocks,
then Clarabel returns `OPTIMAL` with objective `-1.9999999856185144`;
construction took 4.146486958s, lowering 1.105972836s, solve 9.779200728s,
remote wall 66.26s, and peak HWM stayed about 1.94 GiB.  This is a
solver-compatibility smoke only; it is too small to count as Phase 2
bound-quality parity evidence.
The same tiny public-wrapper route is now covered in
`test/relaxations/pauli_chains.jl`: `cs_nctssos(...; direct_linear=true,
axis_rotation_quotient=true, contiguous_rdm_decomposition=:qmbcertify,
formulation=:psd_blocks, representation=:real,
orphan_policy=:free_variables)` solves and matches the existing tiny
moment-variable result to `1e-6`.  Rerunning the Pauli-chain file through
easy-ssh with a Clarabel `SOLVER` passed 2,588/2,588 checks in 304.99s.
The symbolic `MomentProblem` entry point for real PSD-block lowering is also
covered now: `test/relaxations/lowering.jl` passed with a Clarabel `SOLVER` in
35.36s, including a real noncommutative PSD-block model and real free-orphan
handling.  After the perf-harness guard hardening, the same focused lowering
file was rerun directly through easy-ssh and passed 31/31 checks in about 25s
wall time after adding the same standalone `SOLVER` fallback pattern used by
the other relaxation-focused test files.  A no-solver audit of the N=8/order-4
finite-axis quotient plus
QMBCertify `k=8` RDM profile with `linear_state_opt_width=1:3` found
consistent zero LSO rows (`moments=346`, `zero=0`, forced-zero classes `850`)
in 53.69s remote wall.  That is evidence that the current finite-axis quotient
makes these small commutator rows redundant; it is not yet Phase 2B
bound-quality evidence.
A follow-up A1/A2-shaped construction probe then exercised the QMBCertify
width-7 linear-state-opt setting under the same shared-machine guard.  The
first run correctly exposed a coverage gap: width-7 LSO rows referenced
moments outside the base finite-axis quotient support.  Direct-linear
extending-support profiles now allow LSO rows to register those additional
free moment keys, matching the existing `:extend` policy for RDM and PSD
state-opt add-ons while preserving the closed-support coverage check.  The
guarded N=8/order-4 A1-like profile (`k=8`, `:qmbcertify`, `:extend`,
`linear_state_opt_width=7`, no PSD state-opt) now constructs in
10.828337185s, with 346 quotient moments, 345 free keys, 810 LSO zero rows,
27 PSD blocks, max block 120, and solve support `true`.  The A2-like profile
with `psd_state_opt_width=3` constructs in 4.498608772s after warmup, with
349 quotient moments, 348 free keys, the same 810 LSO rows, 32 PSD blocks,
max block 120, QMBCertify RDM blocks `[72, 64, 56]`, and five PSD-state-opt
blocks of size 18.  The focused Pauli-chain regression now includes this
width-7 source-like construction path and passed 2,599 checks in 313.71s
remote wall.  This is still construction/solve-support evidence only; it does
not promote the small probe to reviewed Phase 2 bound-quality parity.
The same A1/A2-shaped construction-only comparison was run at N=12, matching
the length of the existing guarded QMBCertify A1/A2 profile probes while still
staying below the N < 14 shared-machine guard.  The A1-like NCTSSoS profile
constructed in 32.392926516s with 2,143 quotient moments, 2,142 free keys,
810 LSO rows, 31 PSD blocks, max block 120, QMBCertify RDM blocks
`[72, 64, 56]`, solve support `true`, product-cache hit rate
0.8551068883610451, and peak process RSS below 1.9 GiB.  The A2-like profile
constructed in 26.364808459s after warmup with 2,158 quotient moments, 2,157
free keys, the same 810 LSO rows, 38 PSD blocks, max block 120, seven
PSD-state-opt blocks of size 18, and the same product-cache hit rate.  The
two-profile remote wall time was 95.36s.  These rows are the correct
small-L construction comparison against the existing QMBCertify L=12 probes;
they still do not provide an NCTSSoS objective because the source-like solver
formulation remains the active memory bottleneck.
A follow-up N=8 A2-like build-only lowering audit compared the two JuMP
formulations with the same guarded profile (`k=8`, `:qmbcertify`, `:extend`,
`linear_state_opt_width=7`, `psd_state_opt_width=3`) and no solver call.
The PSD-block formulation constructed in 8.973778323s, lowered in
1.245082414s, used 64,846 JuMP variables, 65,308 scalar binding/zero
equalities, 32 PSD-triangle cone constraints, 61.31s remote wall time, and
2,065,528 KiB peak RSS.  The moment-variable formulation constructed in
8.844380433s, lowered in 0.653069896s, used 349 JuMP variables, 811 scalar
equalities, the same 32 affine PSD-triangle cone constraints, 60.72s remote
wall time, and 2,026,072 KiB peak RSS.  Thus the width-7 linear-state-opt
rows add the expected 810 scalar rows, but the source-like solver blocker is
not another mirrored-binding duplication in PSD-block lowering.  The compact
moment-variable build still leaves affine PSD cones for the solver bridge and
the earlier stopped Clarabel run shows that canonicalization/factorization,
not construction, remains the unsafe step for this profile.
A direct SOS-dual build of the same N=8 A2-like source profile constructs in
38.406562867s and dualizes in 2.694283931s, with 65,308 dual variables, only
348 scalar equality constraints, the same 32 PSD-triangle cones, 48.35s remote
wall time, and 1,619,344 KiB peak RSS.  This confirms that explicit scalar
binding rows can be removed in the dual model.  However, a guarded Clarabel
solve of that dual model with a 240s timeout was terminated inside
Clarabel/QDLDL factorization after 241.75s and 42,058,780 KiB peak RSS, with
no objective.  Dualization therefore does not by itself make the N=8
source-like A2 solve safe on the shared machine; the next solver-memory step
must change the cone/formulation structure rather than merely primal-vs-dual
orientation.
Running the same N=8 dual A2-like model with Mosek through the docs
environment changes that conclusion for practical small-source solves:
construction took 38.248220356s, dualization 3.219861635s, Mosek solve
5.186839150s, remote wall 59.27s, peak RSS 1,909,064 KiB, status `OPTIMAL`,
and objective `-3.6510934246065627`.  This is the first safe source-like
NCTSSoS objective for the width-7 LSO plus width-3 PSD-state-opt profile under
the N < 14 guard.  It does not prove L=30/L=100 parity, but it shows the
construction and dual model are solver-usable with the production solver and
that the Clarabel failures should not be used as NCTSSoS formulation failures.
The solver probe harness now reproduces this path with
`NCTS_SOLVER_PROBE_DUALIZE=true`: the N=8 A2-like run through
`perf/pauli_translation_solver_probe.jl` and the docs Mosek environment
constructed in 8.921513553s, dualized in 1.184948720s, solved in
4.668893662s, used 65,308 dual variables, 348 scalar equalities, 32 PSD cones,
810 zero-dual variables, reached `OPTIMAL`, and returned the same objective
`-3.6510934246065627` in 71.82s remote wall time with 2,235,832 KiB peak RSS.
A tiny N=4 dual Mosek harness smoke also solves to
`-2.0000000135072087`, so the harness covers both the optimizer loading path
and the source-like dual model route.
The same harness then produced an N=12 A2-like source objective, still below
the N < 14 guard and using the same `k=8`, width-7 LSO, and width-3
PSD-state-opt settings.  Construction took 31.815833665s, dualization
1.517045922s, Mosek solve 8.364793684s, remote wall 101.24s, and peak RSS
3,152,464 KiB.  The model had 2,158 quotient moments, 2,157 free keys,
810 zero-dual variables, 38 PSD cones, 94,690 dual variables, max block 120,
QMBCertify RDM blocks `[72, 64, 56]`, seven PSD-state-opt blocks of size 18,
and product-cache hit rate 0.8551068883610451.  Mosek returned `OPTIMAL` with
objective `-5.392562452506382`.  This is still not L=30 parity, but it is a
strictly stronger small-L NCTSSoS objective check than the earlier
construction-only N=12 comparison.
The paired N=12 A1-like profile, with the same RDM and width-7 LSO settings
but without PSD state-opt blocks, also solves through the same dual/Mosek
harness.  Construction took 31.511600503s, dualization 1.498736917s, Mosek
solve 7.871238416s, remote wall 100.21s, and peak RSS 3,090,328 KiB.  The
model had 2,143 quotient moments, 2,142 free keys, 810 zero-dual variables,
31 PSD cones, 93,493 dual variables, max block 120, QMBCertify RDM blocks
`[72, 64, 56]`, no PSD-state-opt blocks, and the same product-cache hit rate.
Mosek returned `OPTIMAL` with objective `-5.392564604362093`.  The A1/A2
difference at N=12 is about `2.15e-6`, so these small-source runs should be
treated as solver/backend evidence and rough bound-quality smoke, not as a
strict monotonicity or reviewed QMBCertify parity claim.
The matching L=12 QMBCertify profile-probe totals are now linked directly in
the NCTSSoS A1/A2 fixture rows: NCTSSoS A1 is `-0.005173695084395` below the
QMBCertify A1 total estimate, and NCTSSoS A2 is `-0.005172161960486` below the
QMBCertify A2 total estimate.  These deltas are loader-checked as probe
comparisons only; a red loader check failed on the missing comparison fields
in 14.9s, and the green loader check passed all 739 checks in 13.0s.
The matching N=12 A0-like base profile, with no RDM, no LSO, and no
PSD-state-opt blocks, also solves through the same harness.  Construction took
6.648130367s, dualization 1.190284484s, Mosek solve 7.110387065s, remote wall
73.58s, and peak RSS 2,946,584 KiB.  The model had 1,621 quotient moments,
1,620 free keys, 28 PSD cones, 86,379 dual variables, max block 120, no RDM
or PSD-state-opt blocks, and product-cache hit rate 0.8551068883610451.  Mosek
returned `OPTIMAL` with objective `-5.47924827890285`.  The A0 to A1/A2 jump is
therefore about `0.08668` at N=12, while A1 and A2 are separated only at the
few-`1e-6` level on this small profile.
These three NCTSSoS N=12 source-like dual/Mosek rows are now pinned in
`test/data/expectations/heisenberg_qmbcertify_rdm.toml` under
`nctssos_source_like_solver_probes`.  The expectation-loader regression checks
the A0/A1/A2 objectives, model sizes, RDM/PSO block shapes, LSO row counts,
environment metadata, explicit solver-probe harness provenance, non-parity
comparison scope, and the A0→A1/A2 objective jump.  The fixture test passed
529/529 checks through easy-ssh in 12.9s.  These rows remain NCTSSoS backend
evidence, not reviewed QMBCertify parity rows.
The solver probe harness now has an opt-in
`NCTS_SOLVER_PROBE_EMIT_FIXTURE=true` mode that prints a copyable
`[[nctssos_source_like_solver_probes]]` row with profile, environment,
harness path, model mode, timing, moment counts, RDM/PSO block sizes, dual
variable counts, zero-dual counts, peak RSS, and product-cache hit rate.  A
red N=4/order-1 construct-only grep failed before the emitter existed; the
green construct-only smoke passed in about 55s remote wall, and a matching
N=4/order-1 SOS-dual no-solver smoke passed in about 50s remote wall.  This
keeps future small-N NCTSSoS source-like evidence reproducible without
hand-copying unstructured probe logs.
The same emitter now also prints a copyable
`[nctssos_solver_probe_environments.<id>]` block with Julia version, optional
Mosek package version from the active project, CPU model, RAM, Julia thread
count, BLAS vendor, project name, and harness path.  A red N=4/order-1
construct-only grep failed before the environment block existed; the green
environment-block smoke passed in about 60s remote wall with the N < 14 guard
active.
The environment block also carries optional NCTSSoS commit provenance: when
`NCTS_NCTSSOS_COMMIT` is supplied, the harness emits `nctssos_commit`.
A red N=4/order-1 construct-only grep with `NCTS_NCTSSOS_COMMIT=testcommit123`
failed before the field existed; the green smoke emitted the commit field and
passed in about 60s remote wall.
The solver-probe fixture row also preserves requested symmetry flags, actual
report symmetry flags, and solve-support/blocker metadata.  This matters for
SU(2) extended-support profiles, where `probe_su2_symmetry=true` does not by
itself prove which reducer was active.  A red N=4/order-1 SU(2)-RDM
construct-only grep failed because the TOML row omitted `su2_symmetry`; the
green smoke emitted `su2_symmetry=true`, `report_su2_rdm_symmetry=true`, and
`solve_supported=true` in about 62s remote wall.
The pinned N=12 source-like NCTSSoS fixture rows now require the same
symmetry/support fields.  A red expectation-loader check failed on missing
`reflection_symmetry`, `su2_symmetry`, `report_*_symmetry`, and
`solve_supported` keys; after adding the finite-axis reflected quotient
metadata to A0/A1/A2 rows, `test/expectations_loader.jl` passed 586/586 checks
through easy-ssh in about 9s remote wall.

---

### Phase 2B — Extended RDM (k=3,4) + linear state optimality

**Duration:** 1–2 weeks

**Goal:** Match QMBCertify bound quality.

**Actions:**

1. **Extend RDM to k=3,4 contiguous sites** with U(1) decomposition

2. **Add first-order linear state optimality:** ⟨i[H,A]⟩ = 0 for
   translation-orbit-reduced, sign-invariant monomials A of test width
   a ≤ 7 (at d=4)
   - Specify: Hermitian A only (anti-Hermitian gives trivially zero for
     Hermitian H)
   - Apply per momentum sector

3. **Implement `closure_check`** for state-opt that rejects a > 2d - h + 1

4. **Memory gate:** measure solver RSS at L=30 per constraint tier

**Milestones:**
- L=30 + k=4 RDM + linear state-opt matches QMBCertify A2 within 1e-6
- Solver burden (time, RSS) measured separately from construction

---

### Phase 2C — PSD state optimality

**Duration:** 3–5 days

**Goal:** Add the strongest state-opt constraint.

**Actions:**

1. **Add PSD state-opt blocks** for test width a ≤ 3 (at d=4), with entries
   `ℓ(vHw* - 1/2(Hvw* + vw*H))`.

2. **Produces additional small PSD blocks** per momentum sector.

**Gate:** Attempt N=100 with RDM only if L=30 extrapolates below 64 GiB peak
RSS.  Otherwise Phase 3 is mandatory before scaling.

---

### Phase 3 — SU(2) symmetry reduction

**Duration:** 6–12+ weeks (this is the hardest phase)

**Goal:** Further block reduction via continuous SU(2) symmetry.

**Current status:** Local RDM SU(2) decomposition is implemented, including
real algebraic Schur coefficient-domain histograms in the RDM metric payloads.
A support-complete Pauli-basis moment helper now builds reduced HPSD blocks with
SU(2) transform provenance, rejects objectives not covered by the reduced PSD
moments, and has small no-solver regression coverage for reduced block sizes.
Analytic full-basis and translation-orbit contiguous-chain SU(2) structural
targets now give N=100 block counts without constructing the full basis; the
translation-orbit target also has a reflection-combined shape estimate for the
future reducer, and these shape-only targets are pinned in the QMBCertify base
expectation fixture.  Fixed-support word reflection diagnostics verify, at
small support, that reversal is block-diagonal in the SU(2)-coupled basis and
recover the expected symmetric/antisymmetric multiplicities.  The low-level
support-complete SU(2) moment helper can now verify objective invariance itself
when passed `ops=(σx, σy, σz)`, so callers do not have to rely solely on a
manual `assume_su2_invariant=true` assertion.  Translation-orbit
representatives are deliberately not accepted by that helper: they are
support-complete only modulo translation, while the current SU(2) transform
groups by exact active support.  A separate
`pauli_su2_translation_orbit_basis_summary` diagnostic now validates that
modulo-translation support completeness and matches the analytic
translation-orbit structural targets.  The companion
`pauli_su2_translation_orbit_basis_reduction_plan` now provides stable
support-orbit multiplicity labels for the future translation/SU(2) reducer,
and `pauli_su2_translation_orbit_basis_reduction_diagnostics` reports the
same off-sector/copy accounting used by the exact-support SU(2) diagnostics.
Those reduction-plan blocks also carry coefficient-domain provenance for the
complex square-root algebraic Clebsch-Gordan transforms, and the metric payload
exposes those domains as per-block vectors together with transformed row-count
metrics.  `pauli_su2_translation_orbit_basis_transform_blocks` now builds the
concrete Clebsch-Gordan row transforms on the translation-orbit representative
basis for small structural checks, and the matching reduced-moment helper can
conjugate an orbit-basis moment block into SU(2) multiplicity blocks with
unitarity/off-block/copy residual diagnostics.  A low-level
`pauli_su2_translation_orbit_basis_moment_problem` wrapper now accepts a
caller-supplied translation-orbit polynomial PSD block and emits SU(2)-reduced
HPSD blocks with transform provenance; it deliberately does not build
translation momentum products itself.  The companion
`pauli_su2_translation_orbit_zero_momentum_blocks` helper now builds the
zero-momentum translation polynomial block on an orbit basis and reduces it
through the same SU(2) pipeline; the underlying
`pauli_su2_translation_orbit_momentum_blocks` helper handles arbitrary momentum
sectors and drops the identity row outside k=0, matching the translation fast
path sector split; its returned blocks carry the same stable momentum/spin
labels as the analytic target.  The
`pauli_su2_translation_orbit_momentum_block_bundle` helper now builds all
selected canonical momentum sectors in fast-path order and matches the
analytic target labels and logical block sizes on the small N=6 structural
case.  The `pauli_su2_translation_orbit_moment_problem` wrapper now reduces a
translation-invariant objective to orbit representatives, constructs every
selected momentum block, emits SU(2)-reduced HPSD constraints with the analytic
momentum/spin labels, and linearizes/builds a JuMP model on the same small N=6
case.  Fixed-sector reflection splitting is now implemented in
`pauli_su2_translation_orbit_momentum_blocks`, and real nonfixed
conjugate-pair reflection splitting is implemented through
`pauli_su2_translation_orbit_real_reflection_momentum_blocks`.  The
bundle/moment-problem wrappers now support the all-sector real reflected base
SU(2) path, while the complex-HPSD debug form still requires reflection-fixed
momenta.  The fixed-sector reducer rephases the SU(2) orbit transform by the
momentum shift between canonical support-axis words and actual translation
representatives; the k=N/2 structural case is covered by the small N=10, d=3
target.  The same rephase required the nonfixed antiunitary reflection split to
use transpose, not adjoint, when reducing the reflection action.
The public symbolic `pauli_translation_invariant_moment_relaxation` and
direct-linear `_pauli_translation_base_linear_relaxation` constructors now have
a guarded base-moment translation/SU(2) branch (`su2_symmetry=true`),
supporting both the public default realified PSD form and the complex-HPSD
debugging form; `sign_symmetry=true` is accepted and treated as subsumed by
SU(2) for this branch.  Explicit finite axis-rotation requests are accepted
under `su2_symmetry=true`; standalone symbolic and direct-linear base moment
construction also support axis-isotypic PSD splitting for both unreflected and
reflection-adapted 1D sectors, with sign symmetry treated as subsumed by the
axis group.  The new
`pauli_axis_rotation_action_matrices` diagnostic exposes the finite H/S
generator action and full 24-element group closure on translation-orbit
representatives as sparse signed-permutation matrices; an N=4/order-1 check
verifies the expected H/S images, group closure, element-order distribution,
and central projector ranks.  The N=4/order-2 isotypic transform has block
sizes `[2, 0, 2, 3, 6]`, is orthonormal, and block-diagonalizes all 24 group
actions; the standalone symbolic/direct-linear axis split reduces N=4/order-2
base PSD blocks from `[26, 24, 24]` to
`[2, 4, 6, 12, 1, 2, 6, 12, 1, 4, 6, 12]`, adds 80 axis equality rows, and a
small symbolic and direct-linear public solve gives objective
`-1.9999999624331184` against the ordinary
direct-linear value `-1.9999999950284961`.  A small N=6 no-solver public report check
confirms that sign on/off gives identical SU(2) block labels, block sizes, and
zero-row counts.  Small N=4 reflected-axis construction and solve smokes split
the order-2 base PSD blocks to
`[2, 4, 6, 6, 6, 1, 2, 3, 6, 1, 2, 3, 6, 1, 4, 6, 6, 6]`, add 80 axis
equality rows, and give matching symbolic/direct objectives near `-2`.
The reflected-axis reducer now detects higher-order sectors where reflection
mixes the finite-axis projectors; those sectors fall back to a single
`axis_reflection_mixed_family` block instead of applying an invalid over-split.
A guarded N=10/order-4 target-only profile with reflection, standalone axis
splitting, LSO width 7, and PSO width 2 now reports 32 PSD blocks, solver-facing
max block 120, 23,608 raw moment keys, 1,017 nonzero axis-quotient moment
classes, 2,969 forced-zero moment classes, 29,444 axis equality rows, and
819 LSO rows without constructing or solving an SDP.  A constructed no-solver
N=10/order-4 reflection+axis
profile without state-opt add-ons builds the relaxation in 3.780177270s after
warmup, and a follow-up report-metric run builds it in 3.946993596s after
warmup.  The constructed report now prints 26 PSD blocks, 23,608 moments,
1,017 nonzero axis-quotient moment classes, 2,969 forced-zero moment classes,
29,444 axis equality rows, and peak process HWM about 1.89 GiB.
Small N=4 `SolverConfig` bridge checks accept the reflected standalone route,
the SU(2)-subsumed route, and the automatic `cs_nctssos` direct-linear path.
The Pauli-chain relaxation test file was rerun on 2026-06-30 through
`easy-ssh` with Clarabel supplied as the direct-file `SOLVER`; the main
translation-invariant suite passed 1,995/1,995 tests in 2m23.8s, after the
early Pauli basis/axis/profile testsets also passed.  The reflected-axis
report histogram now preserves `translation_reflection_axis_irrep`
provenance instead of collapsing those blocks into the unreflected
`translation_axis_irrep` bucket.
The direct-linear path now has an opt-in `axis_rotation_quotient=true` path
that replaces finite-axis scalar equality rows by signed quotient moment keys.
For base direct-linear cases it emits already-quotiented PSD forms; custom
bases, scalar side constraints, RDM add-ons, and state-opt add-ons retain the
post-assembly rewrite fallback.  A small N=4/order-1 solve probe matched the
axis-equality objective (`-2.9999999991535025` vs
`-2.9999999991535033`) while reducing the direct-linear model from 7 moments
and 4 axis zero rows to 3 moments and 0 zero rows.  A small no-solver smoke
reported 7 raw axis moment keys, 3 quotient moment keys, and 0 zero rows, with
0.011332798s measured construction after warmup.
On the guarded N=10/order-4 reflection+axis no-solver profile, the quotient
preserves the expected model size reduction from 23,608 raw moment keys to
1,017 quotient moment keys, with 2,969 forced-zero classes and 0 zero rows.
Pre-quotient emission reduced measured construction from the earlier
post-assembly quotient rewrite time of 16.293457953s to 7.309214324s after
warmup; warmup construction was 14.589617s, peak process HWM was about 1.88
GiB, and product-cache hit rate remained 0.831025.  A temp-environment rerun
of `test/relaxations/pauli_chains.jl` through `easy-ssh` passed 1,995/1,995
tests in 2m45.3s after adding the direct test imports (`JuMP`, `COSMO`,
`LinearAlgebra`, `Clarabel`).  A follow-up public API regression verifies that
the symbolic `SolverConfig` wrapper and automatic `cs_nctssos` bridge reject
`axis_rotation_quotient=true` unless `direct_linear=true`; the focused
Pauli-chain file passed 2,198/2,198 translation-invariant relaxation tests in
2m31.2s.  The same public test block now pins the direct-linear quotient
argument guards for missing finite-axis constraints, `su2_symmetry=true`, and
`real_moment_matrix=false`; after correcting the reflected complex test setup,
the focused file passed 2,204/2,204 translation-invariant relaxation tests in
2m32.0s.
A guarded rerun of the N=10/order-4 reflected-axis direct-linear quotient
no-solver profile on 2026-06-30 measured 2.196610180s construction after
warmup, with 1019.57 MiB allocated, 1.701997271s in block assembly,
0.421757607s in precompute, the same 1,017 quotient moments, 2,969 forced-zero
classes, 26 PSD blocks, and 0.831025 product-cache hit rate.
The dense block transform now also skips near-zero left transform coefficients
before walking a source row.  The Pauli-chain regression file still passed
2,204/2,204 translation-invariant relaxation tests in 2m31.6s, and the same
N=10/order-4 profile measured 2.123942361s construction with 1.632159282s in
block assembly, unchanged model counts, and the same product-cache hit rate.
A follow-up sampled profile showed that the remaining quotient overhead was
dominated by repeated vector-key-to-tuple token conversion inside
`_axis_quotient_info` while rewriting every linear form.  The quotient map now
also stores a direct vector-key lookup table, avoiding that conversion in the
hot path while keeping the token map for diagnostics/fallback.  The guarded
N=10/order-4 reflection+axis no-solver profile then measured 5.758779242s
construction after warmup, with the same 23,608 raw moment keys, 1,017
quotient moment keys, 2,969 forced-zero classes, 0 zero rows, and 0.831025
product-cache hit rate; warmup construction was 12.869761s and peak process
HWM was about 1.87 GiB.  The N=4 solve parity probe still matched the
axis-equality objective (`-2.9999999991535025` vs
`-2.9999999991535033`) while reducing 7 moments/4 zero rows to 3 moments/0
zero rows.  A second profile showed another avoidable pre-quotient-base cost:
the constructor still collected raw axis-equality moment support even though
the quotient path no longer emits those equality rows.  Guarding that
collection behind the non-pre-quotient fallback reduced the same N=10/order-4
profile to 4.885385332s after warmup, with block assembly down to 4.410592183s
and the same 1,017 quotient moments/0 zero rows; warmup construction was
11.322285s and peak process HWM was about 1.70 GiB.  A targeted N=4 standalone
axis solve parity probe matched objectives (`-2.9999999898006817` vs
`-2.999999985654257`) while reducing 19 moments/20 zero rows to 3 moments/0
zero rows.  The focused Pauli-chain test rerun passed 1,995/1,995 tests in
2m47.4s.  A tiny-form fast path in
`_linear_moment_form_from_owned_pairs!` now skips generic sorting for one- and
two-term forms; `test/relaxations/moment_linear.jl` covers reversed two-term
ordering and duplicate cancellation.  This was correct but not material for
the N=10 reflected-axis quotient profile, which measured 4.874313453s after
warmup.  The larger quotient-rewrite win came from aggregating forms with more
than eight terms by quotient representative key before final sorting.  That
reduced the guarded N=10/order-4 reflection+axis profile to 4.250164624s after
warmup, with block assembly 3.771922241s, the same 23,608 raw moment keys,
1,017 quotient moments, 2,969 forced-zero classes, 0 zero rows, and 0.831025
product-cache hit rate; warmup construction was 11.463560s and peak process
HWM was about 1.87 GiB.  The focused Pauli-chain test rerun passed 1,995/1,995
tests in 2m46.5s.  Moving the pre-quotient rewrite earlier, immediately after
base momentum-entry construction and before reflection/axis block transforms,
reduced the same guarded N=10/order-4 reflection+axis profile to 2.738170506s
after warmup.  Block assembly dropped to 2.258745390s, with the same 23,608
raw moment keys, 1,017 quotient moments, 2,969 forced-zero classes, 0 zero
rows, 26 PSD blocks, and 0.831025 product-cache hit rate; warmup construction
was 10.117964s and peak process HWM was about 1.88 GiB.  A targeted N=4 solve
parity probe matched objectives (`-2.9999999841286167` vs
`-2.999999985948631`) while reducing 19 moments/20 zero rows to 3 moments/0
zero rows.  The focused Pauli-chain test rerun passed 1,995/1,995 tests in
2m45.1s.  The next refinement pushes the quotient map into product-term to
linear-form construction, so base entries no longer normalize raw moment forms
only to quotient-normalize them immediately afterward.  Non-fixed reflected
real sectors now force the full `2n` realified coordinate space before
reflection adaptation, preserving the transform-domain contract even when the
quotient removes all imaginary forms.  The valid N=10/order-4 reflection+axis
profile measured 2.083705552s after warmup, with block assembly 1.608836549s,
the same 23,608 raw moment keys, 1,017 quotient moments, 2,969 forced-zero
classes, 0 zero rows, 26 PSD blocks, and 0.831025 product-cache hit rate;
warmup construction was 9.115407s and peak process HWM was about 1.82 GiB.  A
targeted N=4 reflection solve parity probe matched objectives
(`-2.999999984403329` vs `-2.999999999437976`) while reducing 19 moments/20
zero rows to 3 moments/0 zero rows.  The focused Pauli-chain test rerun passed
1,995/1,995 tests in 2m47.6s.  A follow-up attempt to aggregate transformed
linear-block entries through a `Dict` reduced some normalization churn but
slowed the same N=10/order-4 profile to 2.317124653s, with block assembly
1.835490310s, so it was not retained.  After removing that experiment, the
guarded N=4 constructor smoke still reported 3 quotient moments, 0 zero rows,
and thirteen 1x1 PSD blocks; the N=10/order-4 reflected-axis no-solver profile
measured 2.166027934s after warmup, with block assembly 1.674878120s, 1.02 GiB
allocated, the same 1,017 quotient moments/2,969 forced-zero classes/26 PSD
blocks, and 0.831025 product-cache hit rate.  The focused Pauli-chain file
passed again in a temporary test environment with the direct test imports
(`COSMO`, `JuMP`, `Clarabel`), including the 1,995/1,995
translation-invariant relaxation testset in 2m44.3s.  Streaming cached
momentum-product terms directly into `LinearMomentForm`s removed another
temporary `(coefficient, monomial)` vector layer for base entries; this was a
small improvement, not a headline win.  The guarded N=10/order-4 reflected-axis
profile moved from 2.166027934s to 2.123480136s after warmup, with block
assembly 1.639242351s, unchanged model counts, and 0.831025 product-cache hit
rate.  The focused Pauli-chain file passed again, including the 1,995/1,995
translation-invariant relaxation testset in 2m44.0s.  The remaining deeper
Phase 1A/1B work is to remove direct-linear product/form materialization
overhead only where profiling shows a material payoff, not to add more axis
equality rows.
The real
base SU(2) PSD blocks now symmetrize floating CG-transform roundoff before
linear lowering, so the generated moment problem passes the JuMP model-build
symmetry check instead of failing on `~1e-16` antisymmetric noise.
Translation-compatible scalar equality rows and scalar inequality PSD blocks
now append through this guarded branch in both symbolic and direct-linear
constructors, including the all-sector reflected real PSD path and the
complex-HPSD debugging path.  First-order linear state-opt zero rows append
through the same branch.  A 2026-06-30 N=4 construction and solve smoke also
validated standalone axis PSD splitting plus linear state-opt rows in both
symbolic and direct-linear constructors: the default sign-filtered path added
80 axis-rotation equality rows and 3 linear state-opt rows without increasing
the direct-linear moment count beyond 70, and both wrappers solved to
`-1.9999999902657695`.
A follow-up 2026-06-30 N=4 solve smoke validated standalone axis base splitting
with unsplit PSD add-ons: `psd_state_opt_width=1, sign_symmetry=false` solved
to `-1.9999999526231225`, and `contiguous_rdm_k=2` solved to
`-1.9999999942112616`, with symbolic and direct-linear wrappers matching in
both cases.
Solve-level diagnostics on N=4 and N=6 showed that the base SU(2) branch is
sound only when off-sector and magnetic-copy rows are reintroduced as
Wigner-Eckart zero/copy constraints.  A naive axis-closure row patch and a
naive transformed off-block/copy patch were both rejected by small solves
(`N=4` exposed the overconstraint immediately).  The high-level
`pauli_translation_invariant_nctssos` wrapper now accepts the validated base
translation/SU(2) reports carrying those Wigner rows, while still rejecting
unvalidated RDM and PSD-state-opt combinations.  The contiguous RDM SU(2)
reducer is separate and already carries its own `spin_magnetic_offblock` and
`magnetic_copy` rows.
`translation_symmetry_profile` now exposes whether a report
has base translation/SU(2) moment reduction and whether that reduction is also
reflection-combined; reports also carry zero-row feature histograms so
linear-state-opt provenance is visible without rebuilding the moment problem.
The same profile distinguishes axis-rotation scalar equality rows from the
still-missing axis-rotation PSD block reducer, and `translation_report_metrics`
forwards that histogram for downstream structural checks.
`translation_solve_support` now gives a machine-checkable supported/blocker
status for translation reports; bare base translation/SU(2) moment reducers
with Wigner zero/copy rows, invariant scalar constraints, scalar
moment-equality rows, and axis-rotation equality rows are now solve-supported.
Linear state-opt rows, PSD state-opt blocks, and closed-support full, U(1),
and SU(2)-decomposed contiguous RDM blocks are also solve-supported after the
small-N public-wrapper checks below.  The blocker payload now also reports unsupported block
features and unsupported zero-row features, so scripts can distinguish the
current RDM or PSD-state-opt blockers without parsing prose.
`perf/pauli_translation_profile.jl` exposes the same lists in constructed
no-solver reports; an N=6, d=2 base-SU(2)+`H^2` moment-equality smoke was
promoted to solve-supported after a small public-wrapper check.  The same
profile now prints the
constructed zero-row feature histogram, so Wigner-Eckart rows, scalar
equalities, linear state-opt rows, and current blockers can be audited without
opening the serialized report object.  `translation_report_metrics` now also
exposes explicit row counters for scalar equalities, axis-rotation equalities,
moment equalities, linear state-opt rows, and the SU(2) Wigner zero/copy row
classes; the no-solver profile prints those counters.  A 2026-06-29 N=6,
d=2 base-SU(2)+`H^2` no-solver profile printed 1 moment-equality row and
1,104 SU(2) base zero rows split as 690 spin off-block, 352 magnetic
off-diagonal, and 62 magnetic-copy rows.  A 2026-06-30 public solve smoke
through both symbolic and direct-linear wrappers matched objectives in
`ALMOST_OPTIMAL` status and leaves no unsupported zero-row features.
`translation_report_metrics` also forwards the translation symmetry-profile
booleans, so scripts can distinguish base SU(2) moment reduction, full/U(1)/
SU(2) RDM add-ons, and ordinary translation constraints without reparsing block
labels.
Before closed-support SU(2)-RDM blocks were promoted onto the base SU(2)
moment reducer, an N=6, d=2 no-solver SU(2)-RDM profile correctly rejected the
default sign-reduced support because the k=2 RDM block required the missing
`σz₁` moment; rerunning with `NCTS_TRANSLATION_SIGN=false` printed separate
RDM-only SU(2) labels and no solver call.
The public translation-relaxation docstring now states the current contract:
base translation moment blocks are SU(2)-reduced for the guarded base branch,
closed-support full/U(1)/SU(2) RDM add-ons, and `psd_state_opt_width` add-ons.
Earlier N=6, d=2 no-solver profiles that printed ordinary
translation moment-block labels for SU(2)-compatible RDM/PSD add-ons are
superseded by the small public-wrapper solve checks below.
`translation_report_metrics` forwards solve-support status for profile scripts,
and ordinary analytic translation structural targets now carry the same
`structural_target_only` solve blocker as the SU(2) target payloads.  A
target-only N=12, d=4 no-solver profile printed `solve supported | false` and
`solve blocker | structural_target_only` for both ordinary translation and
SU(2) structural target tables.  Structural-target payloads also carry empty
unsupported block/zero-feature lists, and the target-only sparse-chain perf
script prints the same blocker fields, so analytic shape reports cannot be
mistaken for constructed-but-unsupported solver paths.  The SU(2) full-basis,
translation-orbit, and translation/reflection structural targets use the same
contract.  The QMBCertify expectation fixtures now record the same
`structural_target_only` blocker and empty unsupported-feature lists for every
shape-only structural target, including RDM/state-opt targets.
Ordinary analytic translation structural targets now also count target-only
axis-rotation scalar equality rows by enumerating orbit-product moment keys
without building site-space bases or PSD blocks; the profile prints those rows
as known add-on zero rows when `axis_rotation_equalities=true`.
The same target payload now counts U(1)-decomposed contiguous-RDM scalar zero
rows exactly; SU(2)-RDM scalar row counts stay construction-dependent because
the Schur transform can split a logical zero entry into real/imaginary rows.
For the default nearest-neighbor Heisenberg width, target-only linear state-opt
rows are counted by enumerating nonzero translated commutator rows without
building PSD blocks.
Target-only `H^2` moment-equality smoke rows are counted the same way, using
the existing translated moment-equality row truncation but stopping before
constraint materialization.
Known target-only zero rows are also exposed as feature/decomposition/reason
histograms so scripts can consume them without parsing profile table text.
The no-solver profile usage block now includes the verified small N=6 PSD
state-opt smoke and N=12 target-only smoke commands, keeping the documented
workflow inside the shared-machine guard.
The target-only profile scripts now enforce the same guard before printing
large structural targets: on 2026-07-02, easy-ssh fail-fast probes confirmed
that `N=14` is rejected for both `perf/pauli_translation_profile.jl` and
`perf/pauli_sparse_chain_d4_blocks.jl` unless the explicit large-run override
is set.  The scripts now also use `abspath(PROGRAM_FILE) == @__FILE__` entry
guards so the expectation loader can include them as harness modules; its
`Pauli Perf Harness Guards` testset covers both `_check_size_guard([14])` and
target-only `main()` fail-fast behavior for the two scripts.  The same entry
guard pattern now covers the remaining perf harnesses, and a follow-up
easy-ssh direct-script smoke confirmed that `heisenberg_mosek_scaling`,
`pauli_charge_singlet_prep`, `pauli_translation_compare`, and
`pauli_translation_solver_probe` also reject `N=14` before any large run.
The QMBCertify profile-probe entrypoint now has the same coverage: with
`NCTS_QMB_NS=14`, `probe_main()` rejects before bootstrap or loading
QMBCertify.  The expectation-loader regression passed through easy-ssh with the
QMBCertify guard testset at 11/11 checks in 3.7s and total remote wall about
66s.
Analytic `pauli_su2_translation_orbit_structural_targets` payloads now carry
their own structural-target solve blocker, so target-only SU(2) shape reports
cannot be confused with constructed SDP relaxations.  They also report the
pre-reduction dense entries, active Wigner-Eckart entries, reduced multiplicity
entries, off-block entries, and magnetic-copy entries; for the N=100, d=4
translation/reflection target this pins 1,335,414 off-block entries and 66,412
magnetic-copy entries as the Phase 3 constraint accounting budget.
The target payload is also sector-resolved (`su2_accounting_records`), and
tests assert that those records sum back to the aggregate dense/active/reduced
entry counts.  Ordinary translation, full-basis SU(2), and translation/SU(2)
structural target payloads also carry block coefficient-domain and exact
coefficient-domain vectors aligned with their block order, plus compact
histograms for downstream scripts and fixtures.
`perf/pauli_translation_profile.jl` now prints `solve supported`, the blocker
tag, and the blocker reason in no-solver reports, so structural targets and
unsupported SU(2) combinations cannot be mistaken for solve-ready fast paths.
These diagnostics also forward from `TranslationInvariantResult`, so solved
fast-path results can be inspected without manually extracting `.report`.
No-solver profile reports also print block coefficient-domain and exact
coefficient-domain histograms, making cyclotomic DFT and algebraic
Clebsch-Gordan provenance visible in the same output used for structural
checks, including target-only full-basis and translation-orbit SU(2)
structural sections.  The QMBCertify expectation fixtures pin the same compact
histograms without duplicating per-block vectors.
Moment equality constraints are now admitted in base translation/SU(2)
construction only when every row-reduced moment is present in the reduced SU(2)
moment basis; uncovered rows still fail the coverage check instead of
reintroducing removed off-sector/copy moments silently.  These extra
constraints remain construction metadata progress only until their reduced SDPs
have solve validation.  The
low-level `pauli_su2_translation_orbit_wigner_eckart_moment_problem` helper now
emits those spin off-block, magnetic off-diagonal, and magnetic-copy zero rows
for one or more selected complex translation/SU(2) momentum sectors, while
preserving the same reduced multiplicity PSD block labels.  The public complex
and real non-reflection construction paths now route through that helper and
preserve Wigner zero/copy provenance when scalar, moment-equality, and
realification rows are appended.  Reflected fixed-momentum and real
conjugate-pair construction paths now reuse the same Wigner zero/copy rows
before reflection splitting, and preserve that provenance in symbolic and
direct-linear forms.  A small direct `N=4, d=2` non-reflection solve probe
after adding the realified Wigner rows initially gave the wrong bound (`-1`
versus the ordinary translation relaxation's `-2`); this was traced to the
nonzero-momentum SU(2) rephase convention.  Using the conjugate momentum
phase in the translation-orbit SU(2) transform restores the N=4 ordinary/SU(2)
objective match to `≈3e-9`, and a regression now protects that sign.  The
reflected `N=4, d=2` public wrapper also matches ordinary translation to
`≈1e-9`, so the public wrapper now accepts the bare base non-reflected and
reflected SU(2) reducer.  Invariant scalar equality/inequality constraints were
validated on the same small chain against ordinary translation (`≈1e-11` or
better), so those scalar combinations are solve-supported too.  Scalar
axis-rotation equality rows were also validated on the same `N=4, d=2` chain
against ordinary translation in both symbolic and direct-linear base SU(2)
paths (`≈8e-10`), so that redundant-but-explicit constraint combination is now
solve-supported.  Closed-support `contiguous_rdm_decomposition=:su2` blocks now
combine with the base SU(2) moment reducer and match the full RDM path on the
same small chain.  PSD state-opt blocks now also combine with the base SU(2)
moment reducer and match the ordinary PSD state-opt path on the same small
chain; the combined SU(2)-RDM plus PSD state-opt add-on path is covered by the
same N=4 symbolic/direct-linear objective check.  Closed-support full and
U(1)-decomposed RDM blocks now also combine with the base SU(2) moment reducer
and match the ordinary full/U(1) RDM paths on the same small chain.  A
guarded N=8, d=2 no-solver profile (2026-06-30, easy-ssh, one thread, N < 14
guard active) measured 0.166585 s for the constructed SU(2)+full-RDM path
(16 PSD blocks, max solver-facing block 8, 1,372 zero rows) and 0.147508 s
for the constructed SU(2)+U(1)-RDM path (18 PSD blocks, max solver-facing
block 4, 1,392 zero rows); both reports were solve-supported.
First-order linear state-opt rows
(`linear_state_opt_width=2`)
were validated on the same small public-wrapper solve against ordinary
translation in both non-reflected and reflected modes (`≈5e-9` and `≈2e-10`),
so that combination is now solve-supported.  Moment-equality rows were promoted
after the N=6 public-wrapper symbolic/direct-linear check above.  Structural
tests now pin the
sign-invariant state-opt candidate filter: width 2 has 12 raw contiguous
Pauli tests and 3 Heisenberg-sign-invariant tests.  A small `N=4, d=2`
moment-equality probe with `H^2` rows was promoted after a small N=6 public
wrapper solve matched the symbolic and direct-linear objectives in
`ALMOST_OPTIMAL` status.  The public solve wrapper now has regression coverage
that accepts base SU(2)+moment-equality reports in both symbolic and
direct-linear modes.
The constructors' report labels/sizes match the analytic small-N targets.  The
analytic translation/SU(2) structural target now also emits stable block labels
for translation and
reflection-combined target blocks, including spin/parity provenance for
reflection-adapted sectors, plus label-size histograms for logical and
solver-facing block comparisons.  Larger numerical fast-path runs remain
pending.
An exploratory source-like N=8, order-4 construction probe with requested
SU(2) RDM decomposition,
`contiguous_rdm_k=8`, `contiguous_rdm_decomposition=:su2`,
`:extend` support, and `psd_state_opt_width=3` was stopped before completion
on shared hardware.  The run was estimated at 1--3 minutes, but after about
5 minutes it was still inside `_linear_moment_form_from_owned_pairs!` /
`_transform_linear_block`, had produced no shape output, and had climbed to
about 26 GB RSS.  This was construction-only and stayed under the N < 14 guard,
but it is not a safe shared-machine probe.  Phase 3 therefore still needs a
pre-lowering SU(2) linear-form reduction strategy before source-like RDM plus
PSD-state-opt SU(2) profiles can be used as solver-memory evidence.  The code
path also showed that `contiguous_rdm_support=:extend` currently bypasses the
base translation/SU(2) moment reducer, so this probe is not yet evidence for
base moment-block SU(2) solver-memory reduction.
A follow-up patch converts the dense local SU(2) Schur transform to the
Pauli linear sparse-row representation before applying it to the RDM linear
block.  The k=8 Schur transform has 8,820 nonzeros out of 65,536 entries
(13.5% density, max 70 nonzeros in a row), so this avoids rediscovering those
zeros during the transform.  The focused `test/relaxations/moment_linear.jl`
rerun passed in 122.47s and covers dense-vs-sparse transform equality.  The
same source-like N=8/order-4 SU(2)+RDM+PSD-state-opt probe was still stopped
after 4m52s at about 24.5 GB RSS with no shape output.  The remaining blocker
is therefore not merely dense Schur zero scanning; the full transformed
SU(2)-RDM zero/copy materialization still needs a more structural emission
strategy.
The direct-linear extended-support SU(2) RDM path now skips those RDM-only
zero/copy rows and emits just the reduced multiplicity PSD blocks.  A small
N=4/order-1/k=3 smoke reports two SU(2) RDM blocks, zero contiguous-RDM zero
rows, and 32 moments in 46.13s remote wall.  The focused
`test/relaxations/pauli_chains.jl` file passed 2,594 checks in 308.44s after
this change.  Contiguous RDM construction also now exploits the signed
permutation structure of Pauli products and canonicalizes one matrix entry at
a time instead of holding the whole raw `entry_terms` matrix.  This is
correctness-preserving and reduces avoidable peak overlap, but the guarded
N=8/order-4/k=8 SU(2) extended-RDM no-PSO probe still climbed to about
14.1 GiB RSS by 2m15s with no shape output and was stopped.  The remaining
implementation gap is direct reduced SU(2)-RDM emission without first storing
the full `2^k × 2^k` computational RDM linear block.
That gap is now closed for direct-linear `:su2` RDMs with
`contiguous_rdm_support=:extend`: the builder aggregates each local RDM entry
directly into translation-orbit moment keys and then emits the SU(2)
multiplicity PSD blocks, without materializing the full computational RDM
block.  A k=3 regression compares the direct reduced blocks against the old
full-RDM-plus-Schur-transform path up to floating-point tolerance, and the
focused `test/relaxations/pauli_chains.jl` file passed 2,596 checks in
308.79s remote wall.  The guarded N=8/order-4/k=8 construction-only probe with
SU(2) RDM decomposition, `:extend` support, and `psd_state_opt_width=3`
completed in 54.11s remote wall (46.329912733s measured construction, peak
RSS 2,277,040 KiB).  It reported 7,691 moments, 7,135 free keys, zero
contiguous-RDM zero rows, 15 PSD blocks, max block 242, RDM blocks
`[28, 56, 40, 14, 1]`, and five PSD-state-opt blocks of size 78.  This is
real construction evidence for the source-like SU(2) path under the N < 14
guard; it is not yet N=100 solver-memory evidence.
The solver probe harness now exposes `NCTS_SOLVER_PROBE_U1` and
`NCTS_SOLVER_PROBE_SU2`, so the same dual/Mosek route can exercise SU(2)
profiles reproducibly.  A tiny N=4 reflected base-SU(2) dual smoke emits
`probe_su2_symmetry=true`, builds 13 PSD cones with max block 4, and dualizes
without solving in 58.87s remote wall.  More importantly, a guarded N=8,
order-4, reflected SU(2) profile with `contiguous_rdm_k=8`,
`contiguous_rdm_decomposition=:su2`, `contiguous_rdm_support=:extend`, and
`psd_state_opt_width=3` now solves through the dual/Mosek harness.  It
constructed in 14.136314575s, dualized in 5.603045518s, solved in
19.746948661s, used 4,098 moments, 4,097 free keys, 50 PSD cones, 22,854 dual
variables, max block 56, SU(2) RDM blocks `[28, 56, 40, 14, 1]`, five
PSD-state-opt blocks of size 18, and returned `OPTIMAL` with objective
`-3.651093409167063` in 97.64s remote wall time with 5,107,360 KiB peak RSS.
This is small-N solver-memory evidence only, but it shows the SU(2) route
does reduce solver-facing block size and remains solvable with the production
solver on a source-like RDM+PSO profile.
Adding width-7 linear state-opt rows to that same N=8 SU(2) RDM+PSO profile
also solves.  Construction took 14.434792432s, dualization 5.889646127s,
Mosek solve 20.484275070s, remote wall 99.34s, and peak RSS 5,282,384 KiB.
The model kept 4,098 moments, 4,097 free keys, 50 PSD cones, max block 56, the
same SU(2) RDM and PSD-state-opt block sizes, and added 819 zero-dual rows.
Mosek returned `OPTIMAL` with objective `-3.6510934100794663`, effectively the
same small-N value as the SU(2) RDM+PSO run without LSO.  This gives a guarded
small-N solve for a closer SU(2) A2-like profile, still not a large-N parity
claim.
The two N=8 SU(2) solve rows are pinned in
`test/data/expectations/heisenberg_qmbcertify_rdm.toml` and checked by
`test/expectations_loader.jl`: the loader verifies the requested/report
symmetry split, SU(2) RDM and PSD-state-opt block sizes, zero-dual accounting,
and objective agreement between the LSO and no-LSO variants.
A guarded N=12 run of the same SU(2) RDM+PSO+LSO profile did not solve within
the shared-machine guard.  Construction completed in 46.361893644s and
dualization in 20.721531983s, producing 28,244 moments, 28,243 free keys,
819 zero-dual rows, 68 PSD cones, 31,455 dual variables, max block 56, SU(2)
RDM blocks `[28, 56, 40, 14, 1]`, and seven PSD-state-opt blocks of size 18.
Mosek was still in the conic optimizer when the 420s timeout terminated the
run; remote wall time was 424.86s and peak RSS reached 53,692,756 KiB.  This
is useful negative Phase 3 evidence: SU(2) reduces the largest PSD block
relative to the finite-axis N=12 A2-like run, but the current equality-heavy
dual form is not yet a safe solver-memory path at N=12.
The solver probe now prints the actual report symmetries separately from the
requested flags (`report_su2_moment_symmetry`, `report_su2_rdm_symmetry`, and
related metrics).  A tiny N=4/order-1 dual no-solver check confirms the
important distinction for `contiguous_rdm_support=:extend`: the profile reports
`probe_su2_symmetry=true` and `report_su2_rdm_symmetry=true`, but
`report_su2_moment_symmetry=false`.  The N=12 timeout is therefore evidence
against the current RDM-extended dual form, not evidence that the full base
translation/SU(2) moment reducer has reached the source-like profile.
The solver-probe fixture emitter now routes its execution-state fields through
a pure `_probe_execution_state` helper, and the expectation-loader checks the
construct-only, primal-lowered, dual-lowered, and solved cases directly instead
of only grepping the source.  The seven pinned NCTSSoS solver-probe TOML rows
now also carry `construction_only=false`, `model_built=true`, and `solved=true`,
and `test/expectations_loader.jl` verifies those fields in the reviewed fixture
block.  A tiny N=4/order-1 construct-only fixture smoke finished in about 65s
remote wall and emitted `termination_status="not_solved"`,
`construction_only=true`, `model_built=false`, and `solved=false`, preserving
the TOML provenance contract after the helper refactor.  The same fixture path
now emits the SU(2) base Wigner zero-row accounting split
(`su2_base_zero_row_count`, `su2_base_spin_offblock_row_count`,
`su2_base_magnetic_offdiag_row_count`, and
`su2_base_magnetic_copy_row_count`).  The seven pinned solved NCTSSoS
solver-probe rows carry the all-zero split for rows without base SU(2) moment
reduction, and the expectation loader checks split consistency.  A tiny
N=4/order-1 construct-only base-SU(2)+extended-RDM smoke finished in about 67s
remote wall and emitted the nonzero split `34 = 4 + 23 + 7`, proving the
machine-readable fixture path carries the Wigner accounting used by the Phase 3
memory-reduction work.  Solver-probe fixture rows now also emit sorted
`construction_stage_*_seconds` fields from the report timing dictionary, so
stage attribution is copyable into reviewed TOML instead of console-only.  A
tiny N=4/order-1 construct-only smoke finished in about 60s remote wall and
emitted `construction_stage_axis_diagnostics_seconds`,
`construction_stage_block_assembly_seconds`, `construction_stage_linearization_seconds`,
and the other small-stage fields in the fixture row.  The translation-profile
repeat summary now also prints a bounded `Repeat Stage Summary` table, sorted
by the largest repeated construction stages.  The red loader guard first failed
on the missing table, the green loader passed in about 61s, and a tiny
N=4/order-1 SU(2)+extended-RDM no-solver repeat smoke finished in 78s remote
wall while printing stable shape plus stage ranges such as `axis_diagnostics`,
`su2_linear_base`, and `su2_base_transforms`.
The direct-linear constructor now has an opt-in guarded
base-SU(2)+extended-SU(2)-RDM merge path (`base_su2_extend_rdm=true`): it
clones the finalized base translation/SU(2) linear data, keeps the
Wigner-Eckart zero/copy rows, and appends the direct reduced SU(2)-RDM blocks
without materializing the full computational RDM block.  The same merge also
appends direct-linear LSO rows and PSD state-opt blocks under the existing
`:extend` free-moment policy.  Tiny N=4/order-1 and N=4/order-2 dual no-solver
smokes through `perf/pauli_translation_solver_probe.jl` now report both
`report_su2_moment_symmetry=true` and `report_su2_rdm_symmetry=true` for
`contiguous_rdm_support=:extend`; the order-2 sign-disabled smoke also reports
PSD state-opt blocks `[3, 6, 6]`.  The focused Pauli-chain regression file
passed in a temp solver environment, including the 2,603-check
translation-invariant relaxation testset in 3m23.8s and total remote wall time
353.25s.  Replacing the Wigner zero-row append with a scalar-only
Hermitian/skew-Hermitian split preserved that small-N behavior: the N=4/order-1
dual no-solver opt-in smoke constructed in 6.366529048s, dualized in
1.018680418s, reported 52 zero-dual rows, kept both SU(2) report flags true,
and finished in 77.26s remote wall with 2,386,368 KiB peak RSS.  A guarded
N=8/order-4 source-like no-solver probe with the opt-in path reached the 420s
timeout before shape output while constructing symbolic base SU(2) Wigner rows
(previously 425.58s remote wall, 8,585,648 KiB peak RSS; the scalar-split rerun
again produced no `construct_seconds` before timeout and needed manual cleanup
after Julia ignored SIGTERM).  The base-SU(2) opt-in path now emits those
Wigner zero/copy rows directly as `LinearMomentForm`s for real moment matrices.
The same N=4/order-1 smoke constructed in 4.054818142s, preserved the previous
PSD histogram `[1 => 3, 2 => 3, 4 => 1]`, kept 52 zero-dual rows, and finished
in 76.05s remote wall with 2,443,508 KiB peak RSS.  The guarded N=8/order-4
source-like no-solver probe now reaches shape output: construction took
242.530448133s, dualization took 108.669236007s, peak RSS reached 28,576,876
KiB, and the model had 7,475 moments, 142,291 zero-dual rows, 57 PSD cones, max
block 56, SU(2) RDM blocks `[28, 56, 40, 14, 1]`, and PSD state-opt blocks
`[18, 18, 18, 18, 18]`.  The command still exited 124 because Julia was killed
during atexit GC cleanup at the 420s guard, after model-shape output had been
printed.  The opt-in path therefore closes the small-N semantic gap and clears
the symbolic Wigner construction timeout at N=8, but source-like sizes still
need zero-row volume and allocation reduction before this path is a safe
default.  The focused Pauli-chain regression file was rerun after the
direct-linear Wigner change and passed all 2,603 checks in 3m26.7s test time
and 356.50s remote wall time.  A generic proportional zero-row dedup pass was
rejected: it reduced the N=4/order-1 opt-in smoke to 25 zero-dual rows, but the
N=8 source-like profile timed out before shape output because normalized-form
key construction was too expensive.  The retained exact-duplicate zero-row
dedup pass is weaker but cheap enough: N=4/order-1 keeps the same PSD histogram
and reduces zero-dual rows from 52 to 37 while constructing in 4.139767196s and
finishing in 76.36s remote wall.  The N=8/order-4 source-like no-solver probe
then constructed in 243.464339453s, dualized in 80.927225967s, finished
successfully in 395.11s remote wall, and reduced zero-dual rows from 142,291
to 115,449, dual variables from 147,547 to 120,705, dualization allocation
from 25.50 GiB to 20.38 GiB, and peak RSS from 28,576,876 KiB to 26,009,052
KiB.  Focused regression was rerun after exact dedup and again passed all
2,603 checks in 3m30.8s test time and 360.96s remote wall.  The first
generation-side reduction is to stop emitting adjoint spin-offblock Wigner
zero rows in the direct-linear real base-SU(2) path: a cross block equal to
zero already implies its adjoint block is zero.  The N=4/order-1 opt-in smoke
kept the same PSD histogram, reduced zero-dual rows to 36, constructed in
3.965136016s, dualized in 0.962532841s, and finished in 75.69s remote wall.
The N=8/order-4 source-like no-solver probe then constructed in
204.244947965s, dualized in 63.587208961s, finished in 336.47s remote wall,
and reduced zero-dual rows to 89,447, dual variables to 94,703, construction
allocation to 58.62 GiB, dualization allocation to 17.17 GiB, and peak RSS to
21,177,936 KiB.  Focused regression was rerun after this generation-side row
skip and again passed all 2,603 checks in 3m27.7s test time and 357.63s
remote wall.  The solver probe now prints the SU(2) base zero-row split by
Wigner reason so these reductions are directly accountable.  Skipping mirrored
magnetic-offdiagonal blocks in the same direct-linear helper reduced the
N=4/order-1 opt-in smoke from 36 to 27 zero-dual rows, with the magnetic
offdiagonal component dropping from 27 to 18 and the PSD histogram unchanged.
On the N=8/order-4 source-like no-solver probe it constructed in
184.536594444s, dualized in 53.053934864s, finished in 298.58s remote wall,
and reduced zero-dual rows to 74,584, dual variables to 79,840, construction
allocation to 53.22 GiB, dualization allocation to 14.73 GiB, and peak RSS to
17,715,412 KiB.  Emitting only the Hermitian half of magnetic-copy matrix
equalities was neutral at N=4 but useful at N=8: the final combined probe
constructed in 183.688220761s, dualized in 57.364053855s, finished in 307.73s
remote wall, and reduced zero-dual rows to 71,726, dual variables to 76,982,
construction allocation to 51.90 GiB, dualization allocation to 14.32 GiB, and
peak RSS to 17,046,568 KiB.  Its SU(2) base zero-row split was 52,113
spin-offblock rows, 15,446 magnetic-offdiagonal rows, and 3,348 magnetic-copy
rows.  Focused regression was rerun after the combined generation-side row
skips and again passed all 2,603 checks in 3m44.2s test time and 377.29s
remote wall.  Further reduction should keep moving duplicate Wigner work out
of row generation; post-hoc dedup helps model size but does not materially
reduce construction time.  Construct-only probe mode was added next so base
construction can be profiled without paying JuMP lowering or dualization.  The
N=8/order-4 source-like construct-only profile finished in 238.89s remote
wall, constructed in 183.021349034s, allocated 51.91 GiB, and kept the same
71,726 zero rows and 57 PSD blocks.  The fine stage split puts the next
bottleneck inside base Wigner row construction:
`su2_base_wigner_rows=113.101290713s`, `su2_base_reflection=29.544035649s`,
`su2_extend_addons=30.934965250s`, `su2_base_finalize=5.753755040s`,
`su2_base_transforms=1.655845452s`, and `su2_base_entries=0.087172341s`.
Adding coarse Wigner sub-stage timers showed the spin-offblock path as the
largest subcomponent: before streaming, the N=8 construct-only Wigner split was
`spin_transform=40.628393113s`, `spin_append=23.347249046s`,
`diag_transform=29.137386146s`, `magnetic_copy_append=13.416433905s`, and
`magnetic_offdiag_append=6.423562891s`.  The spin-offblock zero rows now stream
each transformed entry directly into the zero-row builder instead of first
materializing the full cross-block `Matrix{LinearMomentForm}`.  The N=8
construct-only smoke kept the same 71,726 zero rows and 57 PSD blocks, while
construction moved to 172.652316195s and remote wall to 229.10s.  The
corresponding dual no-solver probe kept 76,982 dual variables, constructed in
180.494032944s, dualized in 59.326808053s, and finished in 304.28s remote
wall.  Focused regression was rerun after streaming and again passed all 2,603
checks in 3m30.6s test time and 362.99s remote wall.
Trying sparse transform rows for all Wigner transforms was rejected: it kept
the shape but regressed the N=8 construct-only profile to 180.722383944s
because diagonal block transforms became slower.  Keeping sparse rows only for
the spin-offblock streaming path is neutral-to-small-positive: the N=8
construct-only profile kept 71,726 zero rows and 57 PSD blocks, constructed in
171.718134489s, and finished in 226.23s remote wall; the N=8 dual no-solver
probe kept 76,982 dual variables, constructed in 180.683454316s, dualized in
49.180073368s, and finished in 292.46s remote wall.  Focused regression passed
again after this mixed sparse-spin change, with all 2,603 checks passing in
3m32.3s test time and 367.21s remote wall.  A follow-up attempt to skip
empty accumulated spin-stream entries was rejected: N=8 construct-only
regressed to 175.537839356s with unchanged shape, so the branch was removed.
Likewise, combining the real/imaginary zero-row coefficient split into one pass
was rejected: the N=4 smoke preserved shape, but N=8 construct-only regressed
to 177.957469529s with the same 71,726 zero rows and 57 PSD blocks, so the
helper was removed.  After reverting that helper, the N=4 smoke again preserved
the 27-row shape, and focused Pauli-chain regression passed all 2,603 checks in
3m04.5s test time and 320.80s remote wall.  Passing the existing builder
registration-token cache into Wigner zero-row emission is retained: it keeps
the same N=4 and N=8 shapes, drops the N=8 construct-only profile to
169.786254456s with 40.30 GiB allocated, and reduces the Wigner spin-stream
stage to 51.005597152s.  The N=8 dual no-solver probe preserved 76,982 dual
variables and 71,726 zero-dual rows, with construction at 168.582903731s,
dualization at 58.413544123s, and 292.56s remote wall.  Focused Pauli-chain
regression passed again with all 2,603 checks in 3m11.0s test time and 329.52s
remote wall.  A narrower attempt to split real/imaginary Wigner spin-stream
rows directly while bypassing the complex intermediate was rejected: N=8
construct-only preserved shape but regressed to 178.834641206s, with the
Wigner spin-stream stage rising to 72.722218911s.  The same registration-token
cache was then threaded into linear-state-opt zero-row registration; this was
neutral-to-small-positive on the N=8 construct-only source-like profile
(169.721296525s, same shape) and focused Pauli-chain regression passed again
with all 2,603 checks in 3m02.5s test time and 315.40s remote wall.  Applying
sparse transform rows only to SU(2) reflection-adapted linear blocks is
retained: it preserves the N=8 source-like shape, reduces construct-only time
to 160.651421227s, and moves the reflection stage to 24.115869384s.  The
corresponding dual no-solver probe preserved 76,982 dual variables and 71,726
zero-dual rows, with construction at 162.073673170s, dualization at
56.141713527s, and 283.98s remote wall.  Focused Pauli-chain regression passed
again with all 2,603 checks in 3m03.6s test time and 316.46s remote wall.  A
follow-up cleanup that stored sparse reflection row bases in the row-block
metadata was rejected: it preserved shape but measured 165.383954188s
construct-only, so the simpler on-demand sparse-row conversion was restored.
The SU(2) extend-addons stage now reports substage timings.  On the N=8
source-like construct-only profile, it preserved shape, constructed in
165.725184012s, and split `su2_extend_addons=32.777054289s` into
`su2_extend_clone=9.684063259s`, `su2_extend_rdm=18.622408677s`,
`su2_extend_finalize=4.212198799s`, `su2_extend_linear_state_opt=0.047296484s`,
and `su2_extend_psd_state_opt=0.042591793s`.  The next add-on work should
therefore target builder cloning/finalization or direct SU(2)-RDM emission, not
state-opt row generation.  The first clone-side fix is retained: the
extend-addons path now uses an append-only shallow builder clone for finalized
base linear blocks instead of deep-copying every base PSD and zero form.  The
N=8 construct-only source-like profile preserved shape, dropped construction
to 147.743282110s, reduced allocation to 35.28 GiB, and reduced peak HWM to
8.95 GiB; the add-on clone substage fell to 0.073619934s.  The N=8 dual
no-solver probe preserved 76,982 dual variables and 71,726 zero-dual rows, with
construction at 151.737425023s, dualization at 56.017400004s, and 269.17s
remote wall.  Focused Pauli-chain regression passed again with all 2,603
checks in 3m02.1s test time and 315.05s remote wall.
The direct SU(2)-RDM append path now reports setup, block-build, realification,
registration/prepare, and append substages.  On the N=8 construct-only
source-like profile, this split showed `su2_extend_rdm=10.159557910s`, dominated
by `su2_extend_rdm_block_build=8.319424908s`; realification and registration
were not the main cost.  A retained Pauli-code action cache factors each local
RDM code into flip mask, sign mask, and phase once per k, instead of
recomputing the k-site action inside every sparse-row contribution.  The N=8
construct-only profile preserved the exact shape, measured
145.556229052s construction, and reduced the RDM block-build substage to
7.449856869s (`su2_extend_rdm=9.336208748s`).  The matching N=8 dual
no-solver probe preserved 76,982 dual variables and 71,726 zero-dual rows, with
construction at 150.331159309s, dualization at 57.597293851s, and 269.29s
remote wall.  The first focused regression exposed an internal old-signature
helper call; a compatibility wrapper now builds the cache for direct tests, and
focused Pauli-chain regression passes again with all 2,603 checks in 3m02.3s
test time and 326.16s remote wall.
The Wigner spin-stream stage now reports its own form-building and zero-append
substage timings.  The instrumentation-only N=8 construct-only profile
preserved shape, constructed in 149.745031587s, allocated 35.29 GiB, and split
`su2_wigner_spin_stream=54.326066098s` into
`su2_wigner_spin_form_build=39.128453389s` and
`su2_wigner_spin_zero_append=15.196159922s`.  The retained follow-up fix adds a
size hint for the transformed spin-offblock form pair buffer before pushing
scaled source-entry terms.  The same N=8 construct-only profile preserved
shape, constructed in 147.676061251s, reduced allocation to 31.40 GiB, and
reduced the Wigner spin stream to 45.037742364s
(`su2_wigner_spin_form_build=33.294894025s`).  The matching N=8 dual no-solver
probe preserved 76,982 dual variables and 71,726 zero-dual rows, with
construction at 147.984362508s, dualization at 55.191441117s, and 264.06s
remote wall.  Focused Pauli-chain regression passes again with all 2,603
checks in 3m01.1s test time and 314.65s remote wall.
The dense linear-block transform now applies the same buffer-sizing rule before
forming transformed entries.  The N=8 construct-only source-like profile again
preserved shape, reduced construction to 144.562581851s, and lowered allocation
to 28.24 GiB; the Wigner diagonal-transform substage measured
30.855271270s.  The matching N=8 dual no-solver probe preserved 76,982 dual
variables and 71,726 zero-dual rows, with construction at 143.918617669s,
dualization at 57.085799868s, and 262.79s remote wall.  Focused Pauli-chain
regression passes again with all 2,603 checks in 3m01.3s test time and
313.36s remote wall.
The sparse linear-block transform now applies the same buffer-sizing rule,
covering reflection-adapted blocks and sparse-row add-on transforms.  The N=8
construct-only source-like profile preserved shape, reduced construction to
139.002709155s, lowered allocation to 24.81 GiB, and moved the reflection
stage to 26.795498037s.  The matching N=8 dual no-solver probe preserved
76,982 dual variables and 71,726 zero-dual rows, with construction at
141.752911207s, dualization at 55.653298373s, and 258.57s remote wall.
Focused Pauli-chain regression passes again with all 2,603 checks in 3m02.0s
test time and 314.78s remote wall.
The Wigner spin form builder now also splits accumulation from reduction/sort
attribution.  The attribution-only N=8 construct-only profile preserved shape,
constructed in 139.798434973s, and split
`su2_wigner_spin_form_build=34.769088684s` into
`su2_wigner_spin_form_accumulate=9.153249927s` and
`su2_wigner_spin_form_reduce=25.615838757s`, confirming reduction/deduplication
as the next local bottleneck.  A retained follow-up adds Vector-specific
`key_lt` and `key_isequal` fast paths for linear moment form reduction, avoiding
generic iterator comparison in the hot sort/dedup loop.  The N=8 construct-only
profile preserved shape, measured 138.376832558s construction and 24.81 GiB
allocation, with `su2_wigner_spin_form_reduce=25.035851618s`.  The matching N=8
dual no-solver probe preserved 76,982 dual variables and 71,726 zero-dual rows,
with construction at 138.265902468s, dualization at 56.753026344s, and 256.72s
remote wall.  Focused Pauli-chain regression passes again with all 2,603 checks
in 3m02.0s test time and 314.38s remote wall.
A follow-up attempt to skip sorting when accumulated reducer input is already
nondecreasing was rejected.  The targeted `test/relaxations/moment_linear.jl`
file passed in 98.77s remote wall, but the N=8 SU(2) construct-only profile
regressed to 139.617294498s construction with
`su2_wigner_spin_form_reduce=30.756852213s`, so the sorted-input prepass was
removed.  After removal, `test/relaxations/moment_linear.jl` passed again in
108.24s remote wall.
The preserving-order coefficient transforms now avoid re-sorting filtered
forms: real-part, imaginary-part, scalar-multiple, and conjugate forms are
constructed as trusted `LinearMomentForm`s because they retain the already
sorted, deduplicated key order.  The targeted `test/relaxations/moment_linear.jl`
file passed in 109.22s remote wall.  The N=8 construct-only source-like profile
preserved shape, measured 137.006697241s construction, allocated 24.81 GiB, and
kept 71,726 zero rows and 57 PSD blocks.  The matching N=8 dual no-solver probe
preserved 76,982 dual variables and 71,726 zero-dual rows, with construction at
135.256776051s, dualization at 57.331337278s, and 253.89s remote wall.  Focused
Pauli-chain regression passes again with all 2,603 checks in 3m01.4s test time
and 313.63s remote wall.
A follow-up builder overload lets `add_zero_constraint!` accept an already
reduced `LinearMomentForm{K,C}` directly.  It still owns/copies the form keys,
but avoids rebuilding and re-sorting forms that were just split into real or
imaginary coefficient views; the builder ownership regression now covers this
path.  The targeted `test/relaxations/moment_linear.jl` file passed in 108.78s
remote wall.  The N=8 construct-only source-like profile preserved shape, cut
construction to 131.489690944s, reduced construction allocation to 23.42 GiB,
and reduced `su2_wigner_spin_zero_append` to 9.821286111s.  The matching N=8
dual no-solver probe preserved 76,982 dual variables and 71,726 zero-dual rows,
with construction at 134.022120528s, dualization at 55.854745187s, and 251.28s
remote wall.  Focused Pauli-chain regression passes again with all 2,603 checks
in 3m04.7s test time and 317.98s remote wall.
A follow-up attempt to replace `_average_pauli_linear_forms` with a sorted
merge was rejected.  The targeted `test/relaxations/moment_linear.jl` file
passed in 108.89s remote wall, but the N=8 construct-only profile was
effectively neutral at 131.448511175s construction and did not improve the
targeted reflection stage, so the simpler concatenate/sort reducer was
restored.
The real SU(2) reflection path now symmetrizes each realified source block once
and transforms it with a symmetric sparse-row transform that computes only one
triangle.  This uses the identity
`sym(U * A * U') = U * sym(A) * U'` for real reflection row bases and preserves
the old output symmetrization semantics while avoiding duplicate transformed
entries.  The targeted `test/relaxations/moment_linear.jl` file passed in
108.62s remote wall.  The N=8 construct-only source-like profile preserved
shape, cut construction to 124.024315138s, reduced construction allocation to
21.84 GiB, and reduced the reflection stage to 15.109653361s.  The matching
N=8 dual no-solver probe preserved 76,982 dual variables and 71,726 zero-dual
rows, with construction at 122.840809390s, dualization at 57.985723035s, and
242.65s remote wall.  Focused Pauli-chain regression passes again with all
2,603 checks in 3m00.5s test time and 313.23s remote wall.
A follow-up attempt to combine the real/imaginary coefficient split for
real-builder zero rows was rejected again.  The targeted
`test/relaxations/moment_linear.jl` file passed in 108.99s remote wall, but the
N=8 construct-only profile was neutral-to-worse at 124.222105211s and
`su2_wigner_spin_zero_append` rose to 13.866331908s, so the separate trusted
real and imaginary split helpers were restored.
A small-form insertion-sort fast path for 3-4 term `LinearMomentForm`
construction was also rejected.  The targeted `test/relaxations/moment_linear.jl`
file passed in 109.52s remote wall, but the N=8 construct-only profile was
neutral-to-worse at 124.462966826s, so the generic sort/reduce path remains.
A follow-up attempt to cache sparse spin-transform row handles and source-entry
lengths inside `_append_sparse_transformed_pauli_su2_zero_linear_block!` was
also rejected.  The N=8 construct-only profile preserved shape but regressed to
127.241774471s construction, with `su2_wigner_spin_zero_append` rising to
15.713764642s, so the simpler direct indexing loop was restored.
Passing row/column provenance separately to
`_add_translation_zero_linear_form!` so component labels could be assembled in
one merge was also rejected as noise.  The N=8 construct-only profile preserved
shape and measured 123.870184427s construction, but the target
`su2_wigner_spin_zero_append` stage was worse at 15.168898963s, so the existing
row/column label merge remains.
A `Vector{<:Unsigned}` key comparator specialization using direct `==` and `<`
for element comparison was rejected.  The targeted
`test/relaxations/moment_linear.jl` file passed in 109.49s remote wall, but the
N=8 construct-only profile preserved shape and regressed to 125.095278739s, so
the existing generic Vector fast path remains.
A bounded `Int8` state-sign table in `_pauli_rdm_code_actions` was rejected for
direct SU(2)-RDM block construction.  The targeted
`test/relaxations/moment_linear.jl` file passed in 109.43s remote wall, but the
N=8 construct-only profile preserved shape and regressed to 125.442903170s;
`su2_extend_rdm_block_build` worsened to 13.648459173s, so the cheaper
`count_ones(sign_mask & state)` path remains.
An internal trusted zero-constraint append path is retained for freshly
generated translation linear forms.  Public `add_zero_constraint!` still copies
caller-owned forms, preserving the existing mutation-safety contract, but
`_add_translation_zero_linear_form!` now uses the trusted path after registering
its generated real/imaginary forms.  The targeted
`test/relaxations/moment_linear.jl` file passed in 110.03s remote wall.  The
N=8 construct-only source-like profile preserved shape, cut construction to
100.032625567s, reduced allocation to 16.99 GiB, and lowered
`su2_wigner_spin_zero_append` to 5.863164530s.  The matching N=8 dual no-solver
probe preserved 76,982 dual variables, 71,726 zero-dual rows, and 57 PSD cones,
with construction at 99.433888686s, dualization at 51.500074183s, and 207.61s
remote wall.  Focused Pauli-chain regression passed all 2,603 checks in 3m01.6s
test time and 314.89s remote wall.
The Wigner spin sparse-transform path now uses a direct coefficient accumulator
for large generated forms and keeps the vector/sort reducer for small forms.
This targets the same N=8 bottleneck after the trusted zero-row append removed
the dominant copy cost; the accumulator size hint is capped so it does not
scale unboundedly with raw pair count.  The targeted
`test/relaxations/moment_linear.jl` file passed in 110.40s remote wall.  The
final N=8 construct-only source-like profile preserved shape, reduced
construction further to 91.068124893s, lowered allocation to 16.62 GiB, and
moved `su2_wigner_spin_form_reduce` from 23.440694175s to 8.573884865s while
`su2_wigner_spin_form_accumulate` rose to 13.143410653s; net
`su2_wigner_spin_stream` was 27.505764923s.  The matching N=8 dual
no-solver probe preserved 76,982 dual variables, 71,726 zero-dual rows, and 57
PSD cones, with construction at 91.387842977s, dualization at 52.279801098s,
and 201.44s remote wall.  Focused Pauli-chain regression passed all 2,603
checks in 3m00.7s test time and 313.02s remote wall.
The dense linear-block transform now uses the same large-form accumulator while
keeping the vector/sort path for small forms.  This targets the Wigner diagonal
transform stage after spin-stream reduction stopped dominating alone.  The
targeted `test/relaxations/moment_linear.jl` file passed in 110.11s remote
wall.  The N=8 construct-only source-like profile preserved shape, reduced
construction further to 83.485636931s, lowered allocation to 16.37 GiB, and
cut `su2_wigner_diag_transform` from 25.038815940s to 16.640285791s.  The
matching N=8 dual no-solver probe preserved 76,982 dual variables, 71,726
zero-dual rows, and 57 PSD cones, with construction at 83.792025404s,
dualization at 52.567356195s, and 193.77s remote wall.  Focused Pauli-chain
regression passed all 2,603 checks in 3m03.7s test time and 316.36s remote
wall.
The real symmetric sparse reflection transform now uses the same large-form
accumulator while preserving the triangle-only transform and mirror assignment.
The targeted `test/relaxations/moment_linear.jl` file passed in 109.81s remote
wall.  The N=8 construct-only source-like profile preserved shape, reduced
construction further to 75.415068466s, lowered allocation to 14.65 GiB, and cut
`su2_base_reflection` from 13.774774737s to 7.142088329s.  The matching N=8
dual no-solver probe preserved 76,982 dual variables, 71,726 zero-dual rows,
and 57 PSD cones, with construction at 76.865449786s, dualization at
51.486976939s, and 184.52s remote wall.  Focused Pauli-chain regression passed
all 2,603 checks in 3m00.9s test time and 313.24s remote wall.
Reusing one scratch `Dict` inside the Wigner spin sparse-transform loop was
rejected.  The targeted `test/relaxations/moment_linear.jl` file passed in
108.91s remote wall, and the N=8 construct-only profile preserved shape and
reduced allocation to 12.35 GiB, but total construction regressed to
76.518400791s against the 75.415068466s retained baseline, so per-form
accumulators remain.
Applying the large-form accumulator to the generic sparse linear-block
transform was also rejected as too broad for a noise-level speed change.  The
targeted `test/relaxations/moment_linear.jl` file passed in 110.33s remote
wall, and the N=8 construct-only profile preserved shape and reduced
allocation to 14.42 GiB, but construction was effectively neutral at
75.395349900s versus the 75.415068466s retained baseline, so the generic
sparse transform keeps its simpler vector/sort reducer.
A final current-tree verification after rejecting that broad helper change
kept the same shape.  The targeted `test/relaxations/moment_linear.jl` file
passed in 99.05s remote wall.  The N=8 construct-only profile reported
76.753755227s construction, 14.65 GiB allocated, 7,475 linear moments, 71,726
zero rows, and 57 PSD blocks in 140.21s remote wall.  The matching dual
no-solver probe preserved 76,982 dual variables, 71,726 zero-dual rows, and 57
PSD cones, with construction at 76.341669910s, dualization at 52.486021496s,
and 187.06s remote wall.  Focused Pauli-chain regression passed all 2,603
checks in 3m01.7s test time and 314.51s remote wall.
Precomputing the local RDM Pauli-code sign parity table was also rejected.
An isolated k=8 SU(2)-RDM block probe first measured per-sector block build
times of 2.819930078s, 5.581690690s, 1.352912341s, 0.022404288s, and
0.000491602s for spin2 `0,2,4,6,8`.  Replacing the innermost
`count_ones(sign_mask & state)` parity check with a precomputed table produced
2.645041679s, 5.647858169s, 1.344058004s, 0.023252523s, and 0.000493055s:
noise-level mixed movement with a slight regression in the dominant spin2=2
sector, so the simpler direct parity check remains.
Small SU(2)-RDM reduced blocks now use a guarded dense orbit accumulator
instead of one `Dict` per PSD entry when `orbit_count * block_dim^2` is at
most 8,000,000 slots; larger k=9/10-style blocks keep the existing sparse
fallback.  The same isolated k=8 SU(2)-RDM block probe reduced per-sector
block build times to 1.758103962s, 2.728100170s, 0.542031046s, 0.017904615s,
and 0.000407713s, with allocation drops in the three dominant sectors.  A new
small k=3 regression compares the dense and `Dict` reducers on the same SU(2)
RDM block.  The targeted `test/relaxations/moment_linear.jl` file passed in
110.94s remote wall.  The N=8 construct-only source-like profile preserved
shape and improved construction to 72.747815634s with 14.25 GiB allocated;
`su2_extend_rdm_block_build` fell to 3.552244089s and total
`su2_extend_rdm` to 5.598634750s.  The matching dual no-solver probe preserved
76,982 dual variables, 71,726 zero-dual rows, and 57 PSD cones, with
construction at 72.934934192s, dualization at 54.784999533s, and 184.14s
remote wall.  Focused Pauli-chain regression passed all 2,603 checks in
3m00.7s test time and 312.82s remote wall.
Replacing vector-key Wigner spin-form accumulators with a sector-local
integer key index was rejected.  The targeted
`test/relaxations/moment_linear.jl` file still passed in 111.61s remote wall,
but the N=8 construct-only profile regressed to 73.206488623s construction and
14.31 GiB allocated against the retained 72.747815634s / 14.25 GiB baseline.
`su2_wigner_spin_form_reduce` improved to 7.941908690s, but
`su2_wigner_spin_form_accumulate` rose to 14.526332913s and total
`su2_wigner_spin_stream` rose to 27.385756366s, so the vector-key sparse
accumulator remains.
The direct-linear Wigner diagonal path now computes only the selected
reference, magnetic-offdiagonal, and magnetic-copy transformed entries that
are actually emitted, instead of materializing each full transformed spin
block first.  This keeps the existing real-Hermitian half-row policy for the
direct SU(2) path.  The N=4/order-1 smoke preserved the retained 27 zero rows
and 6 PSD blocks, constructing in 2.054362574s and finishing in 55.34s remote
wall.  The N=8 construct-only source-like profile preserved the exact shape
with 7,475 linear moments, 71,726 zero rows, 57 PSD blocks, RDM blocks
`[28, 56, 40, 14, 1]`, and PSD state-opt blocks `[18, 18, 18, 18, 18]`;
construction dropped to 63.795308055s with 12.20 GiB allocated, and
`su2_wigner_diag_transform` fell to 0.901893215s while
`su2_base_wigner_rows` fell to 38.511953121s.  The matching N=8 dual
no-solver probe preserved 76,982 dual variables, 71,726 zero-dual rows, and 57
PSD cones, with construction at 63.690291206s, dualization at 52.210195542s,
and 172.79s remote wall.  A primitive selected-entry transform regression was
added to `test/relaxations/moment_linear.jl`, which passed in 100.52s remote
wall.  Focused Pauli-chain regression passed all 2,603 checks in 3m03.3s test
time and 315.61s remote wall.
Follow-up instrumentation split the selected magnetic-offdiagonal append stage
into form construction and zero-row registration without changing the model
shape.  The N=4/order-1 smoke kept 27 zero rows and 6 PSD blocks, constructed
in 2.029806576s, and finished in 65.22s remote wall.  The N=8 construct-only
source-like profile again preserved 7,475 linear moments, 71,726 zero rows,
and 57 PSD blocks; construction measured 64.019398478s with 12.20 GiB
allocated and 118.05s remote wall.  The new attribution shows
`su2_wigner_magnetic_offdiag_append=9.150819925s`, split into
`su2_wigner_magnetic_offdiag_form_build=7.899668198s`
(`form_accumulate=5.219825807s`, `form_reduce=2.679842391s`) and
`su2_wigner_magnetic_offdiag_zero_append=1.248437759s`.  The next offdiag
optimization should target selected-form construction, not builder
registration.
Raising the large-form threshold to force vector/sort reduction only for
magnetic-offdiagonal selected forms was rejected.  The N=4/order-1 smoke
preserved shape, but the N=8 construct-only source-like profile regressed to
66.563998048s construction and 12.28 GiB allocated against the retained
64.019398478s / 12.20 GiB instrumentation baseline.  The offdiag accumulate
substage fell to 1.706590501s, but reduce time rose to 8.962911057s and
`su2_wigner_magnetic_offdiag_form_build` worsened to 10.669501558s, so the
existing `Dict` accumulator threshold remains.
A current-tree verification on 2026-07-01 kept the retained shape and
confirmed the same bottleneck envelope: the guarded N=8 construct-only
source-like profile reported 7,475 linear moments, 71,726 zero rows, 57 PSD
blocks, max block 56, construction time 64.829821879s, 12.20 GiB allocated,
and 118.87s remote wall.  The Wigner split was
`su2_wigner_spin_stream=25.185556015s`,
`su2_wigner_magnetic_offdiag_append=9.779451668s`
(`form_accumulate=3.876902914s`, `form_reduce=4.662379546s`,
`zero_append=1.237338148s`), and
`su2_wigner_magnetic_copy_append=3.878267212s`.
After splitting magnetic-copy attribution into form construction, reference
subtraction, and zero-row registration, the same N=8 profile again preserved
the exact shape.  It constructed in 64.568291705s, allocated 12.20 GiB, and
finished in 118.63s remote wall.  The new magnetic-copy split was
`form_build=1.833220704s`, `form_subtract=0.559551151s`, and
`zero_append=0.328714777s`; the larger local target remains selected
magnetic-offdiagonal form construction at 9.317272850s.  The focused
Pauli-chain regression passed all 2,609 checks through easy-ssh in 321.02s
remote wall after this instrumentation change.
A follow-up attempt to cache the momentum-entry form lengths once per Wigner
pass was rejected.  The N=4 direct-linear SU(2/RDM smoke preserved shape, but
the N=8 source-like construct-only profile was neutral-to-worse overall:
construction measured 64.992571851s and 119.64s remote wall versus the
retained 64.568291705s / 118.63s baseline.  The offdiag form-build substage
moved down to 7.536831431s, but the total profile did not improve, so the
simpler direct `length(entries[a,b])` path remains.
The add-on finalization bucket was then split without changing model shape.
The guarded N=8 source-like construct-only profile measured 63.761812476s
construction and 118.20s remote wall; `su2_extend_finalize=3.707302157s` was
entirely `su2_extend_finalize_linear`, so the next add-on target is
`finalize!(builder)` ownership/finalization cost rather than zero-row dedup.
The focused Pauli-chain regression passed all 2,609 checks through easy-ssh in
325.22s remote wall after this timer split.
The generic `finalize!` routine now accepts optional internal timing hooks,
and the SU(2) extend path passes them through.  A guarded N=8 source-like
construct-only profile preserved the same model shape, constructed in
64.458800937s, and finished in 120.68s remote wall; the finalize split was
`su2_extend_finalize_construct_data=3.474024103s`,
`key_copy=0.382812471s`, `moment_index=0.001627026s`,
`adjoint_key=0.001801639s`, `pivots=0.000249656s`, and
`free_keys=0.000312103s`.  This points at `MomentLinearData` construction and
invariant checks, not pivot discovery.  The targeted
`test/relaxations/moment_linear.jl` file passed through easy-ssh in 101.21s
remote wall, and the focused Pauli-chain regression passed all 2,609 checks in
320.91s remote wall after the finalizer timing hooks.
The `MomentLinearData` invariant checks were then split under the same optional
timing hook.  The guarded N=8 source-like construct-only profile again
preserved the exact shape, constructed in 64.846282937s, and finished in
120.06s remote wall.  The finalization split was
`su2_extend_finalize_construct_data=3.533404839s`, dominated by
`construct_data_zero_constraints=3.025061185s` and
`construct_data_psd_blocks=0.503269642s`; key copy was 0.372524804s, while
moment indexing, adjoint-key setup, pivot discovery, free-key collection, and
key-to-monomial mapping were all negligible.  The next credible finalizer
target is therefore zero-row invariant/self-adjoint checking, not builder
ownership.  After this invariant split, the targeted
`test/relaxations/moment_linear.jl` file passed through easy-ssh in 103.82s
remote wall, and the focused Pauli-chain regression passed all 2,609 checks in
322.41s remote wall.
The zero-constraint invariant bucket was then split one level deeper.  The
targeted `test/relaxations/moment_linear.jl` file passed through easy-ssh in
111.64s remote wall.  The guarded N=8 source-like construct-only profile
preserved the exact shape, constructed in 64.298552927s, and finished in
118.37s remote wall.  The zero-constraint invariant split was
`kind=0.014628794s`, `moment_keys=2.952871206s`, and
`self_adjoint=0.017631821s`.  The remaining finalizer target is therefore
moment-key membership validation across the zero rows.
That target was kept as a conservative internal optimization: public
`add_zero_constraint!` rows still force full membership validation, while
translation zero rows appended through the registered-key internal path skip
only the repeated membership scan at finalization.  A regression now verifies
that unregistered public zero rows still fail during `finalize!`.  The targeted
`test/relaxations/moment_linear.jl` file passed through easy-ssh in 112.70s
remote wall.  The guarded N=8 source-like construct-only profile preserved the
exact shape, constructed in 59.922866072s, and finished in 114.20s remote wall;
`su2_extend_finalize` dropped to 1.689067023s and
`su2_extend_finalize_construct_data_zero_constraints` dropped to 0.833964539s.
The focused Pauli-chain regression passed all 2,609 checks in 321.84s remote
wall after the trusted-key optimization.
A final trusted-row cut skips the remaining repeated form-normalization scan
for internally registered zero rows; `ScalarLinearConstraint` already checks
those forms at insertion, while public rows still take the full validation
path.  The targeted `test/relaxations/moment_linear.jl` file passed through
easy-ssh in 111.88s remote wall.  The guarded N=8 source-like construct-only
profile preserved the exact shape, constructed in 59.038112594s, and finished
in 113.13s remote wall; `su2_extend_finalize` dropped to 0.886018329s and
`su2_extend_finalize_construct_data_zero_constraints` dropped to 0.029550123s.
The focused Pauli-chain regression passed all 2,609 checks in 317.78s remote
wall after this finalizer cut.
Trying to remove bounds checks and bind sparse transform rows in the selected
Wigner transform loops was rejected.  The targeted
`test/relaxations/moment_linear.jl` file passed through easy-ssh in 112.28s
remote wall, but the guarded N=8 source-like construct-only profile was
neutral-to-worse: it preserved the exact shape but constructed in
59.111406225s and finished in 113.75s remote wall versus the retained
59.038112594s / 113.13s baseline.  The code was restored to the simpler loop.
The base SU(2) finalization bucket was then split into zero-row deduplication
and `finalize!` internals.  The targeted `test/relaxations/moment_linear.jl`
file passed through easy-ssh in 112.95s remote wall.  The guarded N=8
source-like construct-only profile preserved the exact shape, constructed in
59.237385659s, and finished in 113.53s remote wall.  The split showed
`su2_base_finalize_deduplicate_zero_constraints=1.938521460s` versus
`su2_base_finalize_linear=0.449870604s`; the internal base `finalize!` split
was already small (`key_copy=0.313040146s`, `construct_data=0.133618559s`,
with zero-constraint checks at 0.030314956s).  The remaining base-finalize
target is therefore zero-row deduplication, not `MomentLinearData` finalizer
checks.
The zero-row dedup loop now avoids a separate `Set` membership probe before
insertion; it records the set length, calls `push!`, and keeps the row only
when the set grew.  The targeted `test/relaxations/moment_linear.jl` file
passed through easy-ssh in 111.73s remote wall.  The guarded N=8 source-like
construct-only profile preserved the exact shape, constructed in
59.138913104s, and finished in 113.64s remote wall; base zero-row dedup was
1.922271156s and base `finalize!` was 0.447387928s.  The focused Pauli-chain
regression passed all 2,609 checks in 330.40s remote wall after the dedup
change.
A temporary count probe of the same N=8 profile showed that the base zero-row
dedup pass saw 70,961 zero rows and kept 70,907 unique rows, so only 54 rows
were duplicates.  The count instrumentation was not retained because it would
misuse the timing dictionary, but the measurement is useful: the remaining
dedup cost is hashing nearly all unique zero-row forms, not removing a large
number of duplicate rows.
The Wigner row emitter now caches the magnetic-level row indices for each
SU(2) transform block instead of repeatedly scanning row labels with
`findall`.  This is a local allocation/lookup cleanup, not a replacement for
the sparse transformed-form accumulator.  The focused Pauli-chain regression
passed all 2,609 checks in 3m01.9s through easy-ssh.  The guarded N=8
source-like construct-only profile, now launched through the perf driver with
`NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true`, preserved the retained shape
(70,907 zero rows, 35 PSD blocks), constructed in 49.357733655s, and finished
the measured build in 49.385744290s; the full script still took roughly two
minutes because it warms up the same case first.  The Wigner row bucket remains
large at `su2_base_wigner_rows=38.159964066s`, with
`su2_wigner_spin_stream=27.894525454s` and
`su2_wigner_magnetic_offdiag_append=7.126991875s`, so the next real target is
still transformed-form construction rather than row-label lookup.
A follow-up attempt to cache selected sparse-transform row subsets by magnetic
level inside the magnetic-offdiag loop was rejected.  It preserved shape, but
the guarded N=8 source-like construct-only profile regressed slightly to
49.567424979s construction and 49.596365706s measured build wall time.  The
broad Wigner bucket fell in that noisy run, but magnetic-offdiag itself
worsened to `su2_wigner_magnetic_offdiag_append=7.411221196s`, so the extra
wrapper cache was removed.
A capped pair-count prepass for sparse transformed forms was also rejected.
The idea was to stop counting input terms once the capped `Dict` size hint was
already saturated.  It preserved shape and one run measured 49.179436113s
construction, but the confirmation run regressed to 50.495294950s construction
and 50.524314479s measured build wall time, with no reliable improvement in
the Wigner form-build buckets.  The exact pair-count loops were restored.
The Wigner zero-row append path now pre-registers each momentum sector's base
moment keys once and then skips per-row key registration for the generated
real/imaginary zero constraints.  This is a local bookkeeping optimization:
all Wigner rows are linear combinations of the already-built sector momentum
entries, so the zero-row appender does not need to rediscover those keys
70,907 times.  A small N=4/order-1/k=3 assertion smoke preserved the 27 SU(2)
base zero rows, 32 linear moments, six PSD blocks, and product-cache counters
in about 40s through easy-ssh.  On the guarded N=8/order-4 source-like
base-SU(2)+extending-SU(2)-RDM profile, measured construction improved from
48.828452855s to 44.307946176s with identical model counts; the Wigner bucket
fell from 38.664837418s to 33.904247516s, mostly by reducing zero-append time
(`spin_zero_append` 6.379749314s to 1.502470769s,
`magnetic_offdiag_zero_append` 1.278062326s to 0.491364483s,
`magnetic_copy_zero_append` 0.317884028s to 0.093392991s).  The remaining
dominant work is still Wigner form accumulation/reduction, not registration.
A follow-up attempt to use a special unique-pair finalizer for `Dict`
accumulator output was rejected: it preserved counts, but the same N=8 profile
regressed to 44.895987530s construction, so the general owned-pair reducer was
restored.
A two-stage sparse spin off-block transform was also rejected.  It reused
left-transformed partial forms across right transform rows and preserved the
N=4 smoke counts, but the N=8 source-like profile regressed to 45.924092980s
construction and 13.05 GiB allocated because the partial forms cost more than
the avoided recomputation.
A direct real/imaginary zero-row accumulator was rejected as unsafe.  It
looked attractive because the real builder currently constructs a complex
`LinearMomentForm` and then splits it, but splitting before complex-form
reduction changed tolerance behavior: the N=8 source-like profile emitted
70,909 zero rows instead of 70,907 and regressed to 47.176844358s
construction, so the complex-first reduction was restored.
After restoring only the registration optimization, the focused
`test/relaxations/pauli_chains.jl` file passed 2,861/2,861 checks in 3m10.3s
test time through easy-ssh with the direct-file Clarabel `SOLVER`.
After the perf-harness guard work, the same focused Pauli-chain file was rerun
through easy-ssh and passed 2,861/2,861 checks in 3m12.4s test time.  The
direct-file fallback now uses Clarabel when `TestUtils.jl` has not already
provided the CI COSMO `SOLVER`, so focused standalone runs do not require test
extras in the root project environment.  The perf/profile harnesses now label
the reported axis-orbit reduction ratio as an active diagnostic rather than a
missing block split; `test/expectations_loader.jl` pins that wording and passed
through easy-ssh with the new 4/4 axis-diagnostics checks and 764/764 fixture
checks in about 54s wall time.  The focused Pauli-chain testset was renamed to
match the same wording, and `test/relaxations/pauli_chains.jl` passed again
through easy-ssh with 2,861/2,861 checks in 3m10.8s test time.
`perf/pauli_translation_profile.jl` now also accepts
`NCTS_TRANSLATION_WARMUP=false`, so one-off small-N profile smokes can skip the
compile/warmup duplicate case.  The expectation-loader guard count rose to
21/21 and passed in about 54s wall time, and an N=4/order-1 no-warmup smoke
printed `warmup: false`, omitted the warmup section, and finished in about
65s wall time.
`perf/pauli_translation_solver_probe.jl` now accepts the explicit
`NCTS_SOLVER_PROBE_LOWER_MODEL=false` construct-only switch while preserving
`NCTS_SOLVER_PROBE_LOWER` as a compatibility alias.  Conflicting values are
rejected by the harness parser.  The expectation-loader guard count rose to
24/24 and passed in about 50s wall time, and an N=4/order-1 construct-only
smoke printed `probe_lower_model=false`, constructed in 1.150086127s, emitted
no `lower_seconds`, and finished in about 53s wall time.
The translation profile harness now fail-fast validates
`NCTS_TRANSLATION_BASE_SU2_EXTEND_RDM=true`: constructed profiles require
`NCTS_TRANSLATION_DIRECT_LINEAR=true`, and constructed plus target-only profiles
both require `NCTS_TRANSLATION_SU2=true`,
`NCTS_TRANSLATION_RDM_DECOMPOSITION=su2`, and
`NCTS_TRANSLATION_RDM_SUPPORT=extend`.  Target-only base-SU2+extend reports now
also reject `NCTS_TRANSLATION_REAL=false`, and constructed reflected profiles
reject the same unsupported complex/reflection combination before setup.  This
guard was added after a bad source-like profile command wasted 65s before the
constructor rejected `contiguous_rdm_support=:closed`.  The expectation-loader
guard count first rose to 36/36 for Pauli perf harnesses and passed 937 fixture
checks in 68s remote wall; after the real/reflection tightening, it rose again
to 44/44 and passed the same 937 fixture checks in 58s remote wall.
The solver-probe harness now has the matching fail-fast guard for
`NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true`: it rejects missing SU(2), non-SU(2)
RDM decomposition, non-extend support, and the unsupported complex reflected
constructor before building the Heisenberg case.  The intermediate
expectation-loader run after this solver-probe guard passed 41/41 Pauli perf
harness guards and the same 937 fixture checks in 58s remote wall.
The profile harness now also fail-fast rejects
`NCTS_TRANSLATION_SINGLET_EQUALITIES=true` unless `NCTS_TRANSLATION_SU2=true`,
matching the constructor's requirement before the harness prints setup output
or builds the basis.  The focused expectation-loader run rose to 48/48 Pauli
perf harness guards and still passed the same 937 fixture checks in 60s remote
wall.
The translation profile harness now also accepts `NCTS_TRANSLATION_REPEATS`,
so noisy small-N construction profiles can run multiple measured repeats in
one Julia process after a single optional warmup.  The focused
expectation-loader run rose to 50/50 Pauli perf harness guards and still
passed the same 937 fixture checks in about 70s monitor wall time; the Pauli
perf harness guard section itself took 50.6s.  An N=4/order-1 no-warmup repeat
smoke printed `repeats: 2`, `N = 4 (repeat 1/2)`, and
`N = 4 (repeat 2/2)`, and finished in 74s remote wall.
The repeat harness now also prints a compact repeat summary with construction
time range, outer wall range, allocation range, and a shape-stability flag.
The expectation-loader Pauli perf harness guards rose to 56/56 and passed in
62s remote wall.  A tiny N=4/order-1 no-warmup repeat smoke printed the
summary row with `shape stable = true`, shape
`linear=7, zero=0, psd=10, max=1, rdm=0, pso=0`, construction range
0.020346005s .. 2.557216615s, and finished in 75s remote wall.
The repeat-stage summary depth is now configurable with
`NCTS_TRANSLATION_REPEAT_STAGE_LIMIT`, defaulting to the previous top-six
stage rows.  A red expectation-loader check first failed on the missing
`stage_limit` keyword, then the green run passed the Pauli perf harness guards
with 67/67 checks in 57.801s remote wall.  A guarded N=8/order-4
base-SU2+extended-RDM repeat profile with the limit raised to 14 preserved
shape and finished in 2m31.543s, exposing the next stable buckets:
`su2_wigner_spin_form_reduce=3.784083072s..4.338025927s`,
`su2_extend_rdm_block_build=3.151126679s..4.323082508s`, and
`su2_base_reflection=1.327525842s..2.316803069s`.
The Wigner spin-stream timing now also records per spin-pair buckets such as
`su2_wigner_spin_pair_2_4` whenever construction stage timing is enabled.  A
red N=4 direct SU(2)+extended-RDM check first failed on the missing
spin-pair key in 54.996s remote wall, then the green check passed in 52.371s.
On the guarded N=8/order-4 source-like repeat profile, the new split preserved
shape and showed that `su2_wigner_spin_pair_2_4` is the largest spin-offblock
pair at 2.581409732s..3.426192454s, followed by
`su2_wigner_spin_pair_4_6=1.772064473s..1.874272433s` and
`su2_wigner_spin_pair_4_8=1.477088671s..1.537640608s`.  The next Wigner
optimization should therefore target spin-pair structure rather than another
global threshold tweak.  The broad `test/relaxations/pauli_chains.jl`
regression passed 2,926/2,926 checks in 6m12.398s remote wall after adding the
spin-pair timing regression.
The spin-pair timing now also splits each spin-pair bucket into
`*_form_build` and `*_zero_append` sub-buckets.  A broad full-file RED attempt
was deleted as a validation route after it hit a 360s timeout in an unrelated
earlier section before reaching the intended assertion.  The replacement
narrow N=4/order-1 SU(2)+extended-RDM smoke failed correctly on the missing
`su2_wigner_spin_pair_*_form_build` key in 45.73s wall with 1,580,156 KiB peak
RSS and no swaps.  After returning the local Wigner append helper's form-build
and zero-append counters to the spin-pair call site, the same smoke passed in
54.80s wall with 1,594,592 KiB peak RSS and no swaps.  This is diagnostic
instrumentation for the retained Wigner bottleneck; it does not restore any
large solved run.
The first N=8/order-4 no-warmup profile using that split ran with one Julia
thread, no solver, and no large-size override after telemetry showed 703 GiB
available memory and load averages `10.44, 10.53, 11.70` on 25 CPUs.  It
completed in 1:42.19 wall with 4,344,032 KiB peak RSS and no swaps.  The
profile preserved the retained shape (7,475 linear moments, 70,907 zero rows,
30 PSD blocks, solver-facing max block 56) and measured 27.220579080s report
construction time.  The new split showed `su2_wigner_spin_form_build =
6.874191143s`, with `su2_wigner_spin_form_reduce = 3.996971146s` and
`su2_wigner_spin_form_accumulate = 2.847984928s`; the hottest spin-pair
sub-buckets were `su2_wigner_spin_pair_2_4_form_build = 2.593889145s` and
`su2_wigner_spin_pair_4_6_zero_append = 1.379397784s`.
That profile identified one small reduction cleanup: accumulator updates kept
near-zero complex cancellations until finalization, even though the finalizer
later drops coefficients within the same tolerance.  A RED accumulator smoke
failed on a retained `[1] => 5e-13im` entry in 5.51s wall with 731,444 KiB peak
RSS.  The accumulator updater now deletes entries whose updated value is within
the configured absolute tolerance, matching the finalizer.  The same smoke
passed in 15.85s wall with 770,028 KiB peak RSS after recompilation, and the
N=4/order-1 SU(2)+extended-RDM Wigner timing smoke still passed in 44.69s wall
with 1,569,272 KiB peak RSS.  The focused `test/relaxations/moment_linear.jl`
file passed through easy-ssh in 2:14.86 wall with 1,592,872 KiB peak RSS and no
swaps.
The follow-up N=8/order-4 no-warmup profile ran after telemetry showed 702 GiB
available memory and load averages `11.99, 12.47, 12.40` on 25 CPUs.  It
completed in 1:42.88 wall with 4,209,780 KiB peak RSS and no swaps, again with
the same model shape.  The intended bucket improved:
`su2_wigner_spin_form_reduce` fell from 3.996971146s to 2.496154265s,
`su2_wigner_spin_form_build` from 6.874191143s to 5.382762404s, and
`su2_wigner_spin_stream` from 9.556272858s to 8.123054803s.  Overall
construction stayed effectively flat at 27.519791475s because unrelated
buckets moved in the opposite direction, so this is retained as a local
cleanup and diagnostic improvement, not an end-to-end speedup claim.
An attempted cold N=8/order-4 complex profile with
`base_su2_extend_rdm=true` first exposed a bad wrapper behavior: because the
direct SU(2) base builder was real-moment-matrix-only, the wrapper silently
fell back to symbolic Wigner polynomial construction.  That run was stopped
after about 9 minutes at roughly 11.5 GiB RSS with no construction row.  A
temporary fail-fast guard made the unsupported combination explicit: a red
focused Pauli-chain test first failed because no exception was thrown; the
green rerun passed 2,908/2,908 checks in 3m11.8s test time, and the exact N=8
profile smoke exited in under a minute with the guard message instead of
entering the symbolic path.
The next complex-direct attempt fixed the immediate finalizer blocker by
making complex zero-row self-adjoint coefficient lookup logarithmic in the
sorted form terms instead of linear.  The non-reflected complex direct
SU(2)+extend path now constructs native HPSD blocks at small N; the focused
Pauli-chain regression passed 2,915/2,915 checks in 3m13.9s test time.  The
first N=8/order-4 complex direct profile then completed, but was not yet a
performance win: after warmup, construction took 210.378778745s, allocated
62.29 GiB, and peaked at 10.99 GiB RSS.  It emitted 7,475 linear moments,
114,795 zero rows, 53 native HPSD blocks, solver-facing max block 28, and
eight momentum sectors.  Stage attribution showed the remaining finalizer
problem plainly: `su2_base_wigner_rows=126.274716808s`,
`su2_wigner_spin_zero_append=53.968146449s`, and the two complex finalizers
still spent about 35-36s each in zero-row self-adjoint checks.
A follow-up made generated zero rows carry `trusted_self_adjoint` provenance
and taught `MomentLinearData` finalization to skip the redundant self-adjoint
form scan for those rows.  The complex direct SU(2)+extend regression now
checks that generated zero rows are trusted; the focused Pauli-chain file
passed 2,916/2,916 checks in 3m13.9s test time, and the focused
`test/relaxations/moment_linear.jl` file passed through easy-ssh in about
2 minutes.  Re-running the N=8/order-4 complex direct profile after warmup
reduced construction to 133.768528996s with 62.30 GiB allocated and 10.79 GiB
peak RSS.  The finalizer bottleneck collapsed:
`su2_base_finalize=5.957180874s`,
`su2_extend_finalize=0.755851172s`, and the generated zero-row checks now
appear as millisecond-scale trusted-key/trusted-self-adjoint accounting.  The
remaining complex-direct cost moved to the Wigner row generator and zero-row
append:
`su2_base_wigner_rows=121.326396645s`,
`su2_wigner_spin_stream=88.294524939s`,
`su2_wigner_spin_zero_append=54.384828213s`, and
`su2_wigner_magnetic_offdiag_append=24.539301723s`.
A second follow-up cached stored-key-to-adjoint-key lookup while generating
Wigner zero rows and made the complex splitter respect the existing
`register_keys=false` path.  The focused Pauli-chain regression passed
2,916/2,916 checks in 3m16.0s test time.  The same N=8/order-4 complex direct
profile then measured 104.986241688s construction, 27.63 GiB allocated, and
10.64 GiB peak RSS.  A third follow-up replaced the temporary stored-form
materialization plus separate real/imag passes with a single-pass Hermitian
split, and strengthened the small complex direct SU(2)+extend regression so the
trusted generated rows are independently checked for self-adjointness.  The
focused Pauli-chain file passed 2,917/2,917 checks in 3m15.9s test time.  The
same N=8/order-4 complex direct profile then measured 86.867482055s
construction, 21.72 GiB allocated, and 10.50 GiB peak RSS.  This is a 2.42x
construction speedup and 2.87x allocation reduction relative to the first
completing complex direct profile, and a 1.21x construction speedup over the
adjoint-cache-only version.  A fourth follow-up reused scratch accumulators for
large sparse transformed forms; the focused Pauli-chain file passed 2,917/2,917
checks in 3m18.4s test time, and the N=8/order-4 profile measured
86.417257405s construction with 16.42 GiB allocated and 10.32 GiB peak RSS.
The scratch reuse is mainly an allocation win, not a time win.  A one-sort
real/imag tuple-reduction experiment was rejected because it raised allocation
back to 22.61 GiB without improving construction time.  The remaining stage
costs in the retained version are
`su2_base_wigner_rows=72.126907309s`,
`su2_wigner_spin_stream=51.167156497s`,
`su2_wigner_spin_form_build=36.047727340s`,
`su2_wigner_spin_zero_append=15.115140902s`, and
`su2_wigner_magnetic_offdiag_append=14.317868548s`.  Therefore complex direct
support is now correct formulation support and materially faster, but it is
still not the performance path.
A supported real-moment N=8/order-4 `base_su2_extend_rdm=true` profile then
refreshed the retained bottleneck: after warmup, construction took
44.569070366s and allocated 11.50 GiB, with 6,721 base moments, 7,475 linear
moments, 71,726 zero rows, 35 PSD blocks, solver-facing max block 56, and
product-cache hit rate 0.797342.  Stage attribution still points at Wigner
form construction: `su2_base_wigner_rows=34.028817843s`,
`su2_wigner_spin_form_build=22.902204269s`, `su2_wigner_spin_stream=24.188404878s`,
and `su2_wigner_magnetic_offdiag_form_build=5.861787273s`, while
`su2_extend_rdm=5.564329353s` and finalizers are smaller.  The next real
Phase 3 implementation target is therefore Wigner form accumulation/reduction
and generated zero-row append, not RDM append plumbing or finalizer
self-adjoint scans.
Accumulator-backed linear-form finalization now skips the full duplicate
combining reducer when the accumulator has already enforced unique keys, while
falling back to the general reducer if duplicate adjacent keys are ever seen.
Exact accumulator cancellations are also deleted immediately instead of being
left for final zero filtering.  A narrower follow-up routes `Dict`
accumulator output through a no-scan distinct-key finalizer, keeping the
fallback only for generic claimed-unique pair vectors.  The focused
`moment_linear.jl` file passed through easy-ssh in about 2m10s, including the
new accumulator/fallback, distinct-key, and exact-cancellation regressions.  A
guarded N=8/order-4 reflected
`base_su2_extend_rdm=true` construct-only profile, run twice in one Julia
process, preserved the source-like shape exactly: 7,475 linear moments, 71,726
zero rows, 57 PSD blocks, max block 56, SU(2) RDM blocks
`[28, 56, 40, 14, 1]`, and five PSD-state-opt blocks of size 18.  The warmed
second construction measured 46.383921277s and 8.40 GiB allocated, with
`su2_base_wigner_rows=31.424227183s`,
`su2_wigner_spin_form_build=19.199129498s`, and
`su2_wigner_spin_form_reduce=8.337007465s`.  This is a reducer/allocation win,
not yet a clear total construction-time win against the prior 44.569s warm
profile.  The focused Pauli-chain regression also passed 2,919/2,919 checks
through easy-ssh; the main translation-invariant testset reported 3m17.2s.
Zero-row deduplication now also pre-sizes its exact-form `Set`, avoiding
growth churn while keeping the same exact duplicate semantics.  The guarded
N=8/order-4 reflected `base_su2_extend_rdm=true` construct-only profile again
preserved the exact source-like shape; the warmed second construction measured
45.718648585s and 8.40 GiB allocated, with
`su2_base_finalize_deduplicate_zero_constraints=1.959705876s`,
`su2_base_wigner_rows=31.585523974s`, and the same 71,726 zero rows.  The
dedup substage itself stayed about flat, so this is a harmless allocation
hygiene step rather than a material Phase 3 speedup.  The focused Pauli-chain
regression passed again with 2,919/2,919 checks; the main translation testset
reported 3m20.9s.
Magnetic-copy form subtraction now merges the two already-sorted linear forms
instead of concatenating and re-sorting them, while preserving the old
per-input-term tolerance behavior.  The focused `moment_linear.jl` regression
passed in about 2m00s with a direct subtract-helper test.  The guarded
N=8/order-4 reflected `base_su2_extend_rdm=true` construct-only profile again
preserved the exact source-like shape; the warmed second construction measured
46.770832903s and 8.31 GiB allocated, with
`su2_wigner_magnetic_copy_form_subtract=0.117749770s` instead of the earlier
0.54--1.32s range, `su2_base_wigner_rows=31.579934917s`, and the same 71,726
zero rows.  The focused Pauli-chain regression passed again with 2,919/2,919
checks; the main translation testset reported 3m17.6s.
A Wigner form-builder attribution split then separated the pair-count sizing
pass from actual accumulation and reduction.  The guarded N=8/order-4
construct-only profile again preserved the exact source-like shape; the warmed
second construction measured 47.123981439s and 8.31 GiB allocated, with
`su2_wigner_spin_form_hint=0.066502831s` versus
`su2_wigner_spin_form_accumulate=10.860447010s` and
`su2_wigner_spin_form_reduce=7.773885925s`.  Therefore the pair-count pass is
not a meaningful target.  A one-pass real/imaginary zero-row splitter was
also rejected: the focused `moment_linear.jl` and Pauli-chain regressions
passed, but the warmed N=8/order-4 profile regressed to 48.131598626s.  Raising
the Wigner large-form accumulator threshold from 128 to 512 was rejected for
the same reason: it preserved shape and reduced allocation slightly to 8.26
GiB, but warmed construction was 47.266779428s and spin/offdiagonal reduction
got worse.  The original threshold and separate real/imaginary split were
restored; the final restored-code N=8/order-4 profile preserved the same shape
and measured 46.756872552s with 8.31 GiB allocated, with
`su2_base_wigner_rows=31.612013925s`,
`su2_wigner_spin_form_build=19.720440765s`, and
`su2_wigner_magnetic_copy_form_subtract=0.117364185s`.
Temporary Wigner form-size counters then showed why threshold tuning failed:
the spin path built 26,105 forms, 25,590 of them through the `Dict` path, from
126,079,984 input terms to 37,404,914 output terms; the magnetic-offdiag path
built 7,725 forms, 7,565 of them through the `Dict` path, from 39,619,270
input terms to 11,529,572 output terms.  The maximum hinted input form size
was 12,800 terms and the maximum reduced output form was 2,690 terms.  The
next accepted patch therefore avoided repeated `sizehint!` calls on reused
scratch accumulators instead of changing the accumulator threshold.  The final
no-counter guarded N=8/order-4 profile preserved the exact source-like shape
and measured 45.763406804s with 6.59 GiB allocated, with
`su2_base_wigner_rows=29.141005832s` and
`su2_wigner_spin_form_build=18.195563637s`.
An exact-length `Vector{Pair}` fill path for accumulator finalization was
rejected: the focused `moment_linear.jl` regression passed, but the warmed
N=8/order-4 profile regressed to 46.082027255s and
`su2_wigner_spin_form_reduce=9.270289261s`.  The previous push-with-sizehint
finalizer was restored.
Deferring exact-zero cancellation deletion in the hot `Dict` accumulator was
also rejected: it preserved shape and total time was noisy, but the actual
Wigner form buckets got worse, so immediate deletion remains.  The accepted
per-term accumulation patch instead replaces complex `abs(value) <= atol`
checks with the equivalent `abs2(value) <= atol^2`, avoiding a square root in
the 100M+-term Wigner accumulation loop while keeping the real path unchanged.
The focused `moment_linear.jl` regression passed in about 2m07s with explicit
boundary tests for the helper, and two guarded N=8/order-4 profiles preserved
the exact source-like shape.  Their warmed constructions measured
43.676842513s and 45.148464159s with 6.59 GiB allocated; the second run
reported `su2_base_wigner_rows=28.837198801s`,
`su2_wigner_spin_form_accumulate=9.649552747s`, and
`su2_wigner_magnetic_offdiag_form_accumulate=3.088151665s`.  The focused
Pauli-chain regression passed again with 2,919/2,919 checks; the main
translation testset reported 3m20.3s.
Precomputing `atol^2` once per scaled-form insertion call was tested and
rejected.  The focused `moment_linear.jl` regression still passed in about
2m08s, but the guarded two-run N=8/order-4 profile did not improve the hot
Wigner buckets: the warmed construction measured 44.968868180s with 6.59 GiB
allocated, `su2_wigner_spin_form_accumulate=9.674107647s`,
`su2_wigner_spin_form_build=18.613274586s`, and
`su2_wigner_magnetic_offdiag_form_accumulate=3.088721288s`.  Shape was
unchanged (`7475` linear moments, `71726` zero rows, `57` PSD blocks, max
block `56`), but the extra helper complexity had no measurable payoff, so the
patch was reverted.
Skipping empty transformed Wigner forms before calling the zero-row append
helper was accepted.  This is semantically neutral because the append helper
already emits no row for an empty form, but it avoids unnecessary real/imag
splitting and registration checks after cancellations.  The focused
`pauli_chains.jl` regression passed 2,919/2,919 checks; the main translation
testset reported 3m21.7s, while the full remote job elapsed 387s including
setup/precompile overhead.  The guarded two-run N=8/order-4 profile preserved
the exact source-like shape; the warmed construction measured 44.603565132s
with 6.59 GiB allocated, `su2_base_wigner_rows=28.338168254s`,
`su2_wigner_spin_form_build=17.565127474s`,
`su2_wigner_spin_form_accumulate=9.532774902s`,
`su2_wigner_spin_form_reduce=7.975324067s`, and
`su2_wigner_magnetic_offdiag_form_build=5.274989277s`.
A direct `Pair` comparator for the reducer sort was tested and rejected.  The
shared `moment_linear.jl` regression passed in 138s, but the guarded two-run
N=8/order-4 profile worsened the target Wigner buckets despite a noisy-good
total time: the warmed construction measured 44.067682913s, but
`su2_base_wigner_rows=29.865375268s`,
`su2_wigner_spin_form_accumulate=9.712992323s`,
`su2_wigner_spin_form_reduce=8.663904254s`, and
`su2_wigner_magnetic_offdiag_form_build=5.981465234s`.  Shape was unchanged,
but the original `by=x -> x.first` sort path was restored.
Using the same `abs2 <= atol^2` predicate in accumulator finalization was also
tested and rejected.  The focused `moment_linear.jl` regression passed in 140s,
but the exact guarded N=8/order-4 source-like profile
(`contiguous_rdm_k=8`, SU(2)-RDM support extension, reflected sectors, and
`psd_state_opt_width=3`) preserved shape while worsening the hot Wigner buckets:
the warmed construction measured 45.194249860s with 6.56 GiB allocated,
`su2_base_wigner_rows=29.281911337s`,
`su2_wigner_spin_form_build=18.689043556s`,
`su2_wigner_spin_form_reduce=8.883557491s`, and
`su2_wigner_magnetic_offdiag_form_build=6.254492396s`.  The previous
accumulator-finalizer predicate was restored.
The same predicate swap was then tested in the sorted magnetic-copy
`_subtract_linear_forms` helper.  The focused `moment_linear.jl` regression
passed in 128s after adding complex-boundary checks, but the guarded
N=8/order-4 source-like profile preserved shape while worsening warmed
construction to 47.545115402s with 6.56 GiB allocated.  The old subtraction
predicate was restored; this helper is not a meaningful Wigner target now that
magnetic-copy subtraction is already around 0.1s.
Skipping the final zero-pruning scan for tolerance-pruned `Dict` accumulator
output was also tested and rejected.  The helper regression passed in 139s,
but the same guarded N=8/order-4 profile again preserved shape while measuring
47.183849435s warmed construction with 6.56 GiB allocated.  Since the Wigner
form buckets did not beat the retained baseline, the general distinct-pair
finalizer remains in use.
The sparse Wigner transform loops now skip empty source-entry forms before
calling the scaled-form accumulator.  This removes no-op calls without changing
the transformed form algebra.  The guarded N=8/order-4 source-like profile
preserved the exact shape (7,475 linear moments, 70,907 zero rows, 57 PSD
blocks, max block 56) and measured 44.830253848s warmed construction with
6.56 GiB allocated; this is a small cleanup, not a headline speedup.  The
focused Pauli-chain regression passed 2,925/2,925 checks in 3m20.2s test time
and 375s remote wall.
A follow-up fast return for transformed row pairs with zero source terms was
tested and rejected.  It preserved the same source-like shape, but the guarded
N=8/order-4 profile worsened to 46.823975476s warmed construction with
`su2_base_wigner_rows=31.272867990s`,
`su2_wigner_spin_form_build=18.130518457s`, and
`su2_wigner_magnetic_offdiag_form_build=6.787247615s`.  The early-return patch
was reverted; the retained cleanup is only the per-entry empty-form skip.
After adding repeat support to the translation profile harness, the retained
code was re-baselined on the same guarded N=8/order-4 source-like profile in
one Julia process.  The warmup construction measured 55.656982s; measured
repeats took 48.814211s and 48.650460s, both with 6.56 GiB allocated.  Shape
was stable across repeats at 7,475 linear moments, 70,907 zero rows, 57 PSD
blocks, solver-facing max block 56, five SU(2)-RDM blocks, and five
PSD-state-opt blocks.  The Wigner bucket remained dominant:
`su2_base_wigner_rows` was 31.397780593s / 31.209051128s,
`su2_wigner_spin_form_build` was 18.181837714s / 19.625646538s, and
`su2_wigner_magnetic_offdiag_form_build` was 6.850589098s / 6.560630475s.
The generic sparse block-transform helpers now mirror the Wigner-specific
cleanup by skipping empty source-entry forms before calling the scaled-form
accumulator.  This is a semantic no-op and primarily protects non-Wigner
translation/RDM/state-opt transforms from the same avoidable empty-form work.
The focused `moment_linear.jl` regression passed through easy-ssh in 149s,
including an empty-entry dense-vs-sparse transform check, and the focused
Pauli-chain regression passed 2,925/2,925 checks in 3m34.0s test time
(400s remote wall).
The Wigner transformed-form path now indexes source moment-entry keys once per
Wigner call, accumulates and reduces transformed forms on dense integer keys,
and maps back to canonical moment keys only when appending real/imaginary zero
rows or storing the reduced PSD reference block.  This targets the measured
form accumulation/reduction bucket without changing public `MomentLinearData`
keys.  A red helper test first failed on the missing indexed-entry helper in
about 11s; the helper-level `moment_linear.jl` regression then passed in about
142s.  A tiny N=4/order-1 SU(2)+extended-RDM no-solver smoke preserved
19 linear moments, 27 zero rows, and six PSD blocks in 80.5s remote wall.  On
a guarded N=8/order-4 reflected SU(2)+extended-RDM+PSD-state-opt repeat profile
(`psd_state_opt_width=3`), shape was stable at 7,691 linear moments, 70,907
zero rows, 57 PSD blocks, solver-facing max block 78, five SU(2)-RDM blocks,
and five PSD-state-opt blocks.  Measured repeats constructed in
29.165932615s and 27.737109328s with 7.10 GiB allocated; the Wigner bucket
fell to `su2_base_wigner_rows=13.205942586s/12.431789633s`,
`su2_wigner_spin_form_build=7.204083064s/6.477031461s`, and
`su2_wigner_magnetic_offdiag_form_build=2.004770605s/2.391255187s`.  The
focused `moment_linear.jl` regression passed again in about 138s after the
Wigner integration, and `test/relaxations/pauli_chains.jl` passed 2,925/2,925
checks with test-reported time 3m22.7s.
The direct SU(2)-RDM block path now exposes a block-build timing split before
further optimization.  A red `moment_linear.jl` check first failed on the
missing `stage_times_ns` keyword; the focused file then passed in about 126s,
including dense-vs-dict timing parity.  A tiny N=4/order-1
SU(2)+extended-RDM no-solver smoke preserved 19 linear moments, 27 zero rows,
and six PSD blocks in 68.8s remote wall while printing
`su2_extend_rdm_block_accumulate`,
`su2_extend_rdm_block_finalize`, and
`su2_extend_rdm_block_build`.  On the guarded N=8/order-4 reflected
SU(2)+extended-RDM+PSD-state-opt profile, the run exited 0 in 1m45.292s remote
wall and preserved 7,691 linear moments, 70,907 zero rows, 57 PSD blocks,
solver-facing max block 78, five SU(2)-RDM blocks, and five PSD-state-opt
blocks.  The RDM bucket measured `su2_extend_rdm=6.232251855s`, with
`su2_extend_rdm_block_build=4.340646595s`,
`su2_extend_rdm_block_accumulate=3.256644622s`, and
`su2_extend_rdm_block_finalize=1.038504237s`; the next RDM target is therefore
the dense accumulation loop, not final entry materialization.
The dense SU(2)-RDM accumulator now stores dense scratch arrays as
`entry_count x orbit_count`, matching the hot inner loop's entry-varying access
pattern.  The focused `moment_linear.jl` regression, including the
dense-vs-dict RDM check, passed in 2m22.314s remote wall.  On the same guarded
N=8/order-4 reflected SU(2)+extended-RDM+PSD-state-opt no-solver profile, the
run exited 0 in 1m41.563s remote wall, preserved the same 7,691 linear moments,
70,907 zero rows, 57 PSD blocks, solver-facing max block 78, five SU(2)-RDM
blocks, and five PSD-state-opt blocks, and reduced reported construction from
35.620843398s to 34.493630680s.  The RDM bucket improved from
`su2_extend_rdm=6.232251855s` to `5.193881751s`, with
`su2_extend_rdm_block_build` down from `4.340646595s` to `3.440338998s` and
`su2_extend_rdm_block_accumulate` down from `3.256644622s` to `2.293019486s`;
`su2_extend_rdm_block_finalize` rose slightly to `1.104490463s`.
Replacing the accumulator reducer's complex `abs(coef) <= atol` check with the
shared squared-norm `_within_abs_atol` predicate was tested and rejected.  The
focused `moment_linear.jl` regression still passed in 2m22.933s remote wall,
but the same guarded N=8 profile regressed to 35.981040219s construction:
`su2_wigner_spin_form_reduce` improved to 4.056195350s, but total Wigner and
RDM staging worsened (`su2_base_wigner_rows=14.884942972s`,
`su2_extend_rdm=5.538883312s`, `su2_extend_rdm_block_build=4.082760916s`).
The reducer predicate change was reverted; the retained RDM improvement is
only the dense scratch storage layout.
A combined indexed real/imag projection helper for Wigner zero-row append was
also tested and rejected.  A red primitive check failed on the missing helper
in 23.931s, and the implemented helper passed the focused `moment_linear.jl`
file in 2m34.173s.  The guarded N=8 profile then regressed to
36.057665536s construction despite local reductions in some form-reduce
sub-buckets: `su2_wigner_spin_form_build=6.896804289s`,
`su2_wigner_spin_form_reduce=3.848063238s`, and
`su2_wigner_magnetic_offdiag_form_reduce=0.890205135s`, while zero append and
other staging worsened (`su2_wigner_spin_zero_append=1.902531783s`,
`su2_wigner_magnetic_copy_zero_append=0.940777909s`,
`su2_extend_rdm=5.595083496s`).  The helper and its test were reverted.
Raising the Wigner sparse-form accumulator threshold from 128 to 512 indexed
terms is retained.  The focused `moment_linear.jl` regression passed in
2m33.888s remote wall.  A single guarded N=8 profile was neutral on total
construction (34.726238323s versus the prior 34.493630680s single-run
layout-only measurement) but improved the target Wigner buckets:
`su2_base_wigner_rows=13.287616402s`,
`su2_wigner_spin_form_build=6.755963445s`,
`su2_wigner_spin_form_reduce=3.809173949s`, and
`su2_wigner_magnetic_offdiag_form_reduce=0.842658986s`.  The warmed repeat
profile settled the tie: repeats constructed in 27.614597928s and
28.014444064s with stable shape, while the main stage ranges were
`su2_base_wigner_rows=12.361250308s..12.678834718s`,
`su2_wigner_spin_stream=8.548904135s..8.702133920s`,
`su2_wigner_spin_form_build=5.528483278s..6.421527513s`, and
`su2_extend_rdm=5.234199464s..5.864203154s`.  A follow-up 1024-term threshold
also passed the focused `moment_linear.jl` file in 2m21.035s remote wall and
improved the same warmed profile slightly further: repeats constructed in
27.428677922s and 27.441530104s with stable shape, while stage ranges were
`su2_base_wigner_rows=12.166824051s..12.590421396s`,
`su2_wigner_spin_stream=8.426647256s..8.628046412s`,
`su2_wigner_spin_form_build=5.871597445s..6.264049218s`, and
`su2_extend_rdm=5.111810590s..5.686848042s`.  The retained threshold is
therefore 1024: it keeps the simpler vector/sort reducer for more medium-sized
indexed Wigner forms and reserves the `Dict` accumulator for larger forms.  The
broader Pauli-chain regression then passed 2,925/2,925 checks in 6m14.534s
remote wall at the 512 threshold, then passed again at the retained 1024
threshold with 2,925/2,925 checks in 6m15.010s remote wall.
An online Wigner zero-row dedup experiment was rejected.  The red test failed
as expected on the missing cache keyword and the focused `moment_linear.jl`
green run passed in 2m20.440s, but the guarded N=8 profile regressed:
construction rose to 32.508567382s, `su2_base_wigner_rows` rose to
16.960150379s, and the final zero-row shape changed from 70,907 to 70,919
when final dedup was skipped.  The experiment was reverted.
A follow-up 2048-term Wigner sparse-form accumulator threshold was also
rejected.  It passed the focused `moment_linear.jl` file in 2m33.414s and kept
the 70,907 zero-row shape, but the guarded N=8 profile regressed from the
1024-threshold single-run baseline: construction rose from 27.751051523s to
28.084812006s and `su2_base_wigner_rows` rose from 12.799602772s to
13.030600339s.  The retained threshold remains 1024.
An SU(2)-RDM parity-sign lookup table was also rejected.  The focused
`moment_linear.jl` file passed in 2m22.125s, and one single-run profile looked
promising on the RDM sub-bucket (`su2_extend_rdm=4.456959360s`,
`su2_extend_rdm_block_build=3.017365193s`).  The warmed two-repeat profile
settled it as noise: construction ranged 27.115666496s..27.955448092s versus
the retained 1024-threshold baseline 27.428677922s..27.441530104s, while
`su2_extend_rdm=4.851732128s..5.941539919s` was not a durable improvement.
The lookup was reverted.
Temporary reflection sub-stage instrumentation showed that reflected SU(2)
base construction is dominated by the real symmetric linear block transform:
`su2_base_reflection_transform=5.090019013s` out of
`su2_base_reflection=5.495615474s`, with reflection matrix construction,
row-basis eigensolves, and metadata all negligible.  Raising that transform's
Dict-accumulator threshold to 1024 was rejected: a single guarded N=8 run
constructed in 27.474193006s, but the warmed two-repeat profile ranged
27.050036025s..28.002468597s, slightly worse on average than the retained
baseline.  The threshold change and temporary instrumentation were reverted.
The retained reflection improvement instead indexes linear-form keys before
the reflection-adaptation transform, applies the sparse block transform on
dense integer keys, and rekeys once after adaptation.  The focused
`moment_linear.jl` regression passed in 2m20.276s.  A single guarded N=8
profile preserved the model shape (`linear=7691`, `zero=70907`, `psd=57`,
`max=78`, `rdm=5`, `pso=5`) while reducing construction to 25.207561442s and
`su2_base_reflection` to 1.457387816s, down from the 1024-threshold single-run
baseline of 27.751051523s construction and 4.998605001s reflection.  The
warmed two-repeat profile confirmed the result with construction
24.868588662s..25.490104162s and the same shape, versus the retained
pre-indexing baseline 27.428677922s..27.441530104s.  The broad
`test/relaxations/pauli_chains.jl` regression then passed 2,925/2,925 checks
in 6m9.617s remote wall.
After the reflection change, a small implementation audit found two Wigner
hot-path checks still using the old hardcoded 128-term accumulator threshold
even though the retained Wigner constant was 1024.  Wiring those Wigner helpers
to `_PAULI_SU2_WIGNER_SPARSE_FORM_ACCUMULATOR_THRESHOLD` kept the generic and
real-reflection transform thresholds unchanged.  The focused
`moment_linear.jl` regression passed in 2m19.490s, and the guarded N=8
warmup+two-repeat profile preserved shape while improving construction to
24.289958669s..24.625747897s.  The broad
`test/relaxations/pauli_chains.jl` regression passed 2,925/2,925 checks in
6m13.321s remote wall after this threshold wiring fix.
A follow-up attempt to flatten SU(2)-RDM sparse-transform right rows before
the dense/dict accumulation loops was rejected.  The focused
`moment_linear.jl` regression passed in 2m19.485s, but the guarded N=8
warmup+two-repeat profile regressed to
24.506089488s..25.124581025s versus the retained
24.289958669s..24.625747897s range, so the simpler direct transform-row loops
were restored.
A second SU(2)-RDM finalization cleanup was also rejected.  Replacing the
duplicate-aware RDM entry reducer with the distinct-pair reducer passed the
focused `moment_linear.jl` regression in 2m23.862s, including the dense-vs-dict
RDM equivalence check, but the guarded N=8 warmup+two-repeat profile regressed
to 24.965777448s..25.305375236s.  The duplicate-aware finalizer was restored.
An `Int`-key specialization for the accumulator-backed distinct-pair finalizer
was also rejected.  The focused `moment_linear.jl` regression passed in
2m19.713s, and the first guarded N=8 repeat profile looked locally promising
for `su2_wigner_spin_form_reduce` but was inconclusive on total construction
(24.392169311s..25.337338326s).  A confirmation repeat profile regressed
clearly to 25.937127445s..26.038073476s, so the generic distinct-pair finalizer
was restored.
Per-magnetic-sector sparse transform row selections are now cached once per
SU(2) block in the Wigner zero-row loop instead of being rebuilt for every
off-diagonal pair and magnetic-copy row.  This preserves row labels and
algebraic forms; it only removes repeated row-selection allocation.  The first
focused Pauli-chain regression run exposed stale source-base LSO shape checks:
the live sign-canonical row key halves the old raw constructor counts
(`340 -> 170`, `672 -> 336`).  After updating the checks to assert that
dedup relation, the same focused file passed 3,041/3,041 checks through
easy-ssh in 7:58.29 wall with 4,430,076 KiB peak RSS and no swaps.  Because
that focused verification already consumed the run budget and no open question
required profiling, the optional N=8 construct-only profile was deleted from
this pass.
The cache container was then temporarily tightened from a `Dict{Int,...}` to a
dense per-`m2` vector because the spin ladder is ordered and contiguous.  A
corrected small N=4/order-1 no-solver Wigner smoke preserved 27 zero rows and
six PSD blocks in 43.78s wall with 1,558,280 KiB peak RSS; the first attempted
smoke had failed fast in 18.80s because it used the symbolic public relaxation
with the intentionally direct-linear-only `base_su2_extend_rdm=true` keyword.
The full focused Pauli-chain regression then passed 3,041/3,041 checks in
7:54.83 wall with 4,686,068 KiB peak RSS and no swaps.  A subsequent guarded
N=8/order-4 reflected SU(2)+extended-RDM+PSD-state repeat profile showed that
the dense-vector cache did not improve construction: repeats measured
24.624093013s..25.894671860s wall with 7.23 GiB allocated and 4,474,920 KiB
peak RSS, versus the retained pre-vector-cache baseline of about
24.289958669s..24.625747897s.  The dense-vector cache was therefore rejected
and the measured `Dict` cache restored.
After the later Wigner/reflection changes, the earlier rejected SU(2)-RDM
parity lookup was retested against the current exact profile and is now
retained.  The direct block builder uses a per-`k` Boolean state/sign parity
table instead of recomputing `count_ones(sign_mask & state)` inside the
accumulation loop.  For the current `k=8` source-like RDM this table is
`256×256` Boolean entries, so the memory cost is negligible relative to the
existing multi-GiB construction profiles.  The focused `moment_linear.jl`
regression passed in 2:25.10 wall with 1,525,812 KiB peak RSS and no swaps.
A guarded N=8/order-4 reflected SU(2)+extended-RDM+PSD-state profile without
LSO then held the same 70,907 zero-row shape, constructed in
23.596493173s..23.636027156s, and finished in 2:38.43 wall with 4,322,468 KiB
peak RSS.  The exact LSO+PSO source-like profile also stayed below the run
estimate: construction was 23.880491967s..24.163750365s, remote wall was
2:40.11, peak RSS was 4,385,012 KiB, and there were no swaps.  This is a small
retained construction win, not evidence to restore any L=20/L=30/N=100 run.
The Wigner sparse-transform zero-row emitter now skips reducer construction
when the selected source block contributes zero terms.  This is a semantic
no-op because an empty source term set cannot become nonzero after scaling.
The focused `moment_linear.jl` regression passed again in 2:24.52 wall with
1,558,484 KiB peak RSS and no swaps.  The exact N=8/order-4 reflected
SU(2)+extended-RDM+LSO+PSD-state profile then preserved the shape
(`linear=7475`, `zero=71726`, `psd=57`, `max=56`, `rdm=5`, `pso=5`) and
constructed in 23.833881947s..24.048666308s.  The measured spin-stream stage
fell to 7.800025884s..8.353884644s, from the preceding exact-profile
8.521549380s..9.280233222s range.  Remote wall was 2:40.09 and peak RSS was
4,382,852 KiB with no swaps.  This keeps the optimization local to the guarded
N=8 construction path; no L=20/L=30/N=100 run is reopened by this result.
A broader focused Pauli-chain regression was considered as the next
verification step for the two retained hot-path edits, but the live remote
gate rejected it: autodl reported load averages `25.68, 15.11, 13.26` on a
25-core host, with memory headroom still available.  The run was therefore
deleted from this pass's live queue instead of launched under CPU saturation.
A unit-level regression now also pins the Wigner empty-source fast path in the
linear-form primitive testset: empty selected sparse-transform inputs return an
empty form, while the neighboring non-empty case still scales normally.  The
focused `moment_linear.jl` file passed with that guard in 2:14.17 wall with
1,537,620 KiB peak RSS and no swaps.
After load dropped to `11.30, 13.90, 13.21`, the broader focused Pauli-chain
regression was restored to the live queue with an 8-9 minute / <5 GiB estimate
and a 900s timeout.  It passed 3,041/3,041 checks in 7:55.46 wall with
4,524,504 KiB peak RSS and no swaps, matching the estimate.
A single-term transformed-form fast path was tested and rejected.  The focused
`moment_linear.jl` file passed in 2:24.14 wall with 1,523,092 KiB peak RSS,
but the exact N=8 profile regressed against the retained empty-source baseline:
construction widened to 23.885966428s..24.627815924s and spin-stream to
8.428082616s..8.601362258s.  The profile finished in 2:41.42 wall with
4,411,776 KiB peak RSS and no swaps.  The single-term specialization was
removed; the empty-source skip remains retained.
The retained base-SU(2)+extending-SU(2)-RDM direct path is now reachable
through the public `SolverConfig`/`cs_nctssos` Pauli translation fast-path
bridge via `base_su2_extend_rdm=true`.  The flag is direct-linear-only:
symbolic dispatch rejects `base_su2_extend_rdm=true` unless
`direct_linear=true`, avoiding accidental keyword leakage into
`pauli_translation_invariant_moment_relaxation`.  The focused Pauli-chain
regression passed all 2,619 checks in 3m02.4s through easy-ssh; the new checks
cover both successful `cs_nctssos(...; direct_linear=true, su2_symmetry=true,
base_su2_extend_rdm=true, contiguous_rdm_decomposition=:su2,
contiguous_rdm_support=:extend)` routing and the symbolic rejection.
Target-only structural accounting now accepts the same
`base_su2_extend_rdm=true` flag for the non-reflected direct-linear shape:
base moment blocks come from the translation/SU(2) orbit target, extending
SU(2)-RDM PSD sizes match the constructed report's realification policy, and
the target carries SU(2) base dense/active/reduced entry accounting.  The perf
target-only script forwards and prints the flag instead of silently dropping it.
A small N=4/order-1/k=3 diagnostic matched target/report PSD sizes exactly
(`[2, 2, 2, 2, 4, 1]`, RDM `[4, 1]`) in about 47s through easy-ssh, and the
focused Pauli-chain regression then passed all 2,634 checks in 3m03.5s.
A small N=8/order-2 target-only profile with the same flag completed in about
26s and printed the expected base SU(2) accounting plus the extending
SU(2)-RDM solver-facing histogram `[1 => 1, 4 => 1]`.
The same structural target mode now covers reflected base translation/SU(2)
sectors by reusing the reflected orbit target: a small N=4/order-1/k=3
diagnostic matched reflected target/report PSD sizes
`[2, 2, 1, 1, 2, 4, 1]` and the `translation_su2_reflection` feature
histogram in about 51s, and the focused Pauli-chain regression passed all
2,642 checks in 3m03.4s.
A source-like N=10/order-4 reflected target-only profile with
`contiguous_rdm_k=8` completed in about 27s, reported 62 total PSD blocks, max
solver-facing block 56, extending SU(2)-RDM blocks `[28, 56, 40, 14, 1]`, and
12 reflected base-SU(2) accounting records.  This remains structural
shape/accounting evidence only, not a constructed or solved SDP.
The real solver-facing RDM block vector is now a public metric:
`pauli_su2_rdm_metrics(8).real_psd_block_sizes == [28, 56, 40, 14, 1]`, and a
small smoke matched that vector against the reflected source-like structural
target in about 13s.  The focused Pauli-chain regression then passed all 2,651
checks in 3m07.2s.
Target-only base-SU(2)+extending-SU(2)-RDM accounting now also reports
`su2_base_zero_row_budget`, a conservative Wigner zero-row upper bound computed
from structural SU(2) accounting as `2*offblock_entry_count + copy_entry_count`.
The same target now exposes the off-block and magnetic-copy budget components
separately, so profile output can distinguish the two structural sources
without constructing transformed rows.  This is intentionally a budget, not an
exact constructed-row count: exact row counts require building the transformed
linear forms.  A red N=4/order-1/k=3 target-only smoke failed on the missing
field in 21.54s; the green target-only smoke reported budget 54 in 30.74s; the
profile target-only smoke printed the new Wigner budget in 37.98s; and a
constructed N=4 smoke verified budget 54 >= 27 constructed SU(2) base zero rows
in 50.15s.  The focused Pauli-chain regression then passed all 2,814 checks in
3m35.0s test time and 409.45s remote wall time.  A follow-up red smoke failed
on the missing split fields in about 11s; the green target-only split smoke
passed in about 31s; the profile-output smoke printed both split budget rows in
about 33s; and a constructed N=4 comparison confirmed the split budgets dominate
the actual spin/offdiagonal and magnetic-copy row counts in about 46s.  The
reviewed N=100 structural target fixture now pins those split budgets for the
translation and translation-plus-reflection targets; a red loader check failed
on the missing keys in 14.9s, and the green loader check passed all 693 checks
in 12.9s after adding the fixture values.  The same fixture now pins
singlet-channel counts and equality-row budgets for the full-basis,
translation-orbit, and translation-plus-reflection SU(2) N=100 targets; the
red loader check failed on the missing keys in 14.9s, and the green loader
check passed all 702 checks in 12.9s after adding the values.  A follow-up
fixture pass now pins SU(2) storage byte counts and full-basis accounting
closure as well; the red loader check failed on missing storage/accounting
keys in 14.9s, and the green loader check passed all 711 checks in 12.9s.
The full-basis fixture also now pins transformed SU(2) block sizes and
solver-facing real PSD storage fields; the red loader check failed on the
missing fields in 15.0s, and the green loader check passed all 721 checks in
13.2s.
SOS dualization now has an opt-in native Hermitian cone mode,
`sos_dualize(...; hermitian_representation=:native)`, while preserving the
existing real-lift default.  The focused `test/relaxations/moment_linear.jl`
file passed in about 120s after adding the native branch, including a tiny
HPSD solve and residual check.  The solver probe now forwards
`NCTS_SOLVER_PROBE_SOS_HERMITIAN=native`; a tiny N=4/order-1 dual no-solver
smoke built 13 native Hermitian PSD cones, 13 JuMP variables, and no zero-dual
variables in about 60s remote wall, and a construct-only fixture smoke emitted
the new `sos_hermitian_representation` field in about 62s.  This is
formulation evidence only, not yet evidence that the N=8/N=12 source-like
SU(2) memory bottleneck is solved.
A larger non-reflected N=8/order-4 SU(2)-RDM+PSD-state+LSO no-solver probe
first exposed that complex-affine native assembly was too slow: construction
completed in 6.888807061s, but dual assembly produced no shape output after
roughly 210s and was interrupted.  Rewriting native assembly to use the same
real/imaginary coefficient arrays as the lifted path fixed that bottleneck.
The focused `moment_linear.jl` file passed again in about 128s.  On the same
N=8 no-solver profile, native Hermitian dualization then completed in
10.503018589s with 31,758 JuMP variables and 45 Hermitian PSD cones, while
the default real lift took 12.967106417s with 63,800 JuMP variables and 45
real PSD cones.  SOS dual variables now also skip JuMP string-name creation;
the focused file passed again in about 126s, and a refreshed N=8 native probe
preserved the same model size with 10.483842612s dualization.  Native
Hermitian cones are therefore a valid opt-in model-size reduction for complex
HPSD probes, while the real-valued base-SU(2)+extend path remains unaffected.
The public `solve_sdp(...; dualize=true, sos_hermitian_representation=:native)`
route now forwards the same option for ordinary `MomentProblem` and
`MomentLinearData` inputs; a red interface regression failed on the missing
keyword in about 115s, and the green `test/relaxations/interface.jl` run
passed in about 144s while preserving default state/trace dualization.
The high-level generic `cs_nctssos` and `cs_nctssos_higher` paths now forward
the same option to `solve_sdp`.  A red high-level regression failed on the
missing `cs_nctssos` keyword in about 115s; the green interface file passed in
about 131s after switching the smoke to a solver-stable Pauli case and keeping
state/trace default routing intact.
The Pauli translation fast-path bridge now forwards
`sos_hermitian_representation=:native` through the public `SolverConfig`
overload down to `solve_sdp`, instead of leaking the keyword into moment
construction or silently ignoring it.  A red focused check first exposed the
leak as an unsupported construction keyword.  The green fast-path check uses a
small N=4 non-reflected complex-moment chain and verifies that the resulting
dual model contains a `HermitianPositiveSemidefiniteConeTriangle`.  The
focused Pauli-chain file then printed an all-pass summary, 2,861/2,861 checks
in 3m10.9s through easy-ssh; after the summary the Julia process kept spinning
and was manually cleaned up, so this is green test-summary evidence rather than
a clean process-exit signal.  After the perf-harness diagnostics cleanup, the
same focused file was rerun through easy-ssh and exited cleanly with
2,861/2,861 checks passing in 3m11.6s test time and roughly 6 minutes remote
wall time.
The native-Hermitian path now also has a solved small source-like probe.  A
guarded N=8/order-4 complex HPSD SU(2)-RDM+PSD-state+LSO run with
`NCTS_SOLVER_PROBE_SOS_HERMITIAN=native` constructed in 6.636262698s,
dualized in 9.290325384s, solved with Mosek in 41.911939314s, and finished in
about 119s monitored wall time with 8.15 GiB peak RSS.  The model had 4,098
moments, 819 zero-dual rows, 45 native Hermitian PSD cones, 31,758 JuMP
variables, max native block 31, objective `-3.6510934109179667`, and max SOS
residual `2.09e-9`.  This row is pinned as a separate native-cone solver probe
fixture; the expectation loader checks its native cone representation, block
sizes `[14, 28, 20, 7, 1]`, PSD-state-opt blocks `fill(9, 8)`, max block 31,
and objective agreement with the corresponding real-lift small solve.  The
loader passed 872/872 fixture checks through easy-ssh in about 47s, and the
focused `test/relaxations/interface.jl` regression passed through easy-ssh in
about 1m51s, including the native `solve_sdp` and high-level `cs_nctssos`
dualization checks.
The matching non-reflected complex real-lift row was then replayed with the
same N=8/order-4 SU(2)-RDM+PSD-state+LSO profile: construction took
6.496921562s, dualization 12.148392384s, Mosek solve 21.065537317s, and the
run finished in about 98s monitored wall time with 5.60 GiB peak RSS.  The
real-lift model kept the same 4,098 moments, 819 zero-dual rows, 45 PSD cones,
max block 31, block sizes `[14, 28, 20, 7, 1]`, and PSD-state-opt blocks
`fill(9, 8)`, but required 63,800 JuMP variables.  The expectation fixture now
pins this paired row beside the native row, checks the objective agreement to
within `1e-8`, and verifies that the native formulation cuts solver variables
at least in half (`31,758 <= 63,800 / 2`).  The refreshed expectation loader
passed 937/937 checks through easy-ssh in about 47s.

**Prerequisites:**
- SU(2) irreducible tensor operators for Pauli operator space
- Schur-Weyl / total-spin decomposition for physical RDM Hilbert space
- Clebsch-Gordan coefficients with exact or interval arithmetic
- Wigner-Eckart theorem for selection rules
- Compatibility: SU(2) commutes with translations/reflections, but implemented
  bases must be aligned

#### Phase 3A — SU(2) decomposition of local RDMs

**Actions:**

1. **Physical Hilbert space:** (C²)^{⊗k} = ⊕_j V_j ⊗ M_j

2. **PSD condition on multiplicity matrices M_j**

3. **Scope:** Pauli spin systems only (not arbitrary algebras)

4. **Basis/multiplicity sizes for k ≤ 8:** document explicitly

   Implemented local-RDM Schur-Weyl multiplicity blocks:

   | k | `(spin2, multiplicity)` blocks | real PSD target sizes |
   |:--|:--|:--|
   | 0 | `[(0, 1)]` | `[1]` |
   | 1 | `[(1, 1)]` | `[1]` |
   | 2 | `[(0, 1), (2, 1)]` | `[1, 1]` |
   | 3 | `[(1, 2), (3, 1)]` | `[4, 1]` |
   | 4 | `[(0, 2), (2, 3), (4, 1)]` | `[4, 6, 1]` |
   | 5 | `[(1, 5), (3, 4), (5, 1)]` | `[10, 8, 1]` |
   | 6 | `[(0, 5), (2, 9), (4, 5), (6, 1)]` | `[10, 18, 10, 1]` |
   | 7 | `[(1, 14), (3, 14), (5, 6), (7, 1)]` | `[28, 28, 12, 1]` |
   | 8 | `[(0, 14), (2, 28), (4, 20), (6, 7), (8, 1)]` | `[28, 56, 40, 14, 1]` |

   The table is pinned by `pauli_su2_rdm_blocks(k)` tests for `k=0:8` and
   exposed through `pauli_su2_rdm_metrics(k).real_psd_block_sizes`.  Real PSD
   target sizes double only nontrivial multiplicity spaces; singleton
   multiplicity blocks stay scalar-sized when their reduced block is real.

5. **Test:** reproduce existing order-2 singlet results for the singlet
   (J=0) sector

   Current coverage pins the two-qubit `J=0` Schur column as
   `[0, -1/sqrt(2), 1/sqrt(2), 0]` and checks that the triplet columns form
   the complementary projector.  The 4-site Heisenberg
   charge/spatial/singlet problem suite continues to cover the public
   `PauliSingletConstraintSpec` path.

#### Phase 3B — SU(2) irreducible tensor decomposition of moment matrices

**Actions:**

1. **Each non-identity Pauli transforms as spin-1 under SU(2)**

2. **Couple tensor products into total (J, m, α) basis** (natural growth is
   3^k support patterns with identities and multiplicities)

3. **Moment matrix:** M = ⊕_J (I_{2J+1} ⊗ X_J), with X_J ⪰ 0

4. **Exact or interval CG coefficients** (not just Float64 — certification
   needs this)

   The recursive spin-1/2 coupling step now exposes the exact squared
   coefficient data as integer numerators over a shared integer denominator,
   with a separate sign.  The numerical Schur transform still uses Float64
   square roots, but the triplet and singlet two-qubit coefficients are pinned
   as exact data before conversion.  A 2026-07-01 easy-ssh smoke printed the
   exact triplet/singlet tuples and the expected singlet Schur column, and the
   focused Pauli-chain regression passed 2,655/2,655 tests in 3m01.2s.  Full
   interval/exact transform emission remains future certification work.
   The public `pauli_su2_schur_diagnostics(k)` now reports `J_z`, `J^2`, and
   aggregate residuals directly alongside transform unitarity, and the k=4/k=8
   regression checks pin the Phase 3B residual tolerances.  The k=8 diagnostic
   smoke printed max residual `3.11e-15`, and the focused Pauli-chain
   regression passed 2,661/2,661 tests in 3m06.1s through easy-ssh.

#### Phase 3C — Higher-degree singlet/channel equalities

**Actions:**

1. **Beyond axis equality (xx=yy=zz):** full singlet channel relations from
   irreducible tensor decomposition

   `pauli_su2_word_singlet_channels(support_size)` now exposes the spin-0
   rows of the fixed-support Pauli-word SU(2) transform.  This pins the
   irreducible-tensor singlet channel data needed by future higher-degree
   scalar equality generation; support-size 2 is tested as the invariant
   `XX + YY + ZZ` channel up to global phase, and support-size 4 exposes the
   expected three singlet channels.  The same selector is now available for
   ordinary support-complete Pauli bases through
   `pauli_su2_basis_singlet_channels(basis)`, preserving support labels and
   coefficient-domain provenance.  The first equality-emission contract is now
   exposed through `pauli_su2_word_singlet_channel_equalities`,
   `pauli_su2_basis_singlet_channel_equalities`, and
   `pauli_su2_translation_orbit_singlet_channel_equalities`: these return the
   non-singlet SU(2) rows, sparse column forms, and orthogonality residuals
   against the singlet channel.  The support-complete low-level
   `pauli_su2_basis_moment_problem` path can now append those rows as opt-in
   `:Zero` constraints with `PauliSU2SingletChannelEqualityOrigin`
   provenance; the N=4/order-2 smoke maps 44 complex candidate rows to 76
   Hermitian zero components.  The low-level
   `pauli_su2_translation_orbit_basis_moment_problem` path can now append the
   translation-orbit candidates with matching provenance; the N=6/order-2
   smoke maps 11 complex candidate rows to 19 Hermitian zero components.  The
   high-level symbolic translation/SU(2) Wigner-Eckart path can now consume
   the same candidates through `singlet_channel_equalities=true`; the
   N=6/order-2 public smoke adds 119 singlet-channel zero rows while preserving
   solve-support validation.  The reflected symbolic and direct-linear paths
   now accept the same option for both reflection-fixed and all-sector real
   reflection reductions.  A 2026-07-01 easy-ssh smoke passed in about
   7s with support-2 off-channel leakage `1.67e-15`; the focused Pauli-chain
   regression passed 2,688/2,688 tests in 3m03.3s.  The support-complete
   basis selector smoke passed in about 9s on the 49-word order-2 basis with
   residual `3.33e-16`; the focused Pauli-chain regression passed
   2,714/2,714 tests in 3m04.6s.  The equality-candidate smoke passed in 13s
   across fixed-support, support-complete, and translation-orbit bases, and
   the focused regression passed 2,765/2,765 tests in 3m06.1s.  Structural
   targets and the profile script now also report the singlet-channel equality
   candidate row count; the N=4/N=6 smoke checks passed in 25s and the
   target-only profile output smoke found the new rows, followed by a focused
   regression passing 2,768/2,768 tests in 3m20.2s.  The support-complete and
   translation-orbit equality candidates now also carry `basis_forms`, i.e.
   sparse `(basis word, coefficient)` rows ready for a future zero-row emitter;
   the support-complete and translation-orbit smokes passed in 17s and 5s,
   and the focused regression passed 2,772/2,772 tests in 3m40.0s.  The
   support-complete zero-row emitter smoke passed in 19.2s, and the focused
   Pauli-chain regression passed 2,782/2,782 tests in 3m59.3s.

2. **Must be compatible with translation/momentum sector structure**

   `pauli_su2_translation_orbit_singlet_channels(basis, n_sites; momentum)`
   now selects the spin-0 rows from the translation-orbit SU(2) transform
   blocks, including optional momentum rephasing.  The small N=6 structural
   regression pins the two singlet channels from the identity and two-site
   support orbit, verifies orthonormality after momentum rephasing, and keeps
   support-orbit labels available for future equality provenance.  The
   translation/SU(2) structural target now reports the same singlet-channel
   count and support-size breakdown (`[0 => 1, 2 => 1]` for the N=6, order-2
   target), so profiles can account for future equality rows before they are
   emitted.  The translation-orbit smoke passed in about 17s with residuals
   `3.33e-16`, and the focused Pauli-chain regression passed 2,704/2,704
   tests in 3m04.2s.  After adding the structural target hook, the same
   focused regression passed 2,716/2,716 tests in 3m05.9s.  The high-level
   `pauli_translation_structural_targets(...; base_su2_extend_rdm=true)`
   target now forwards the base SU(2) singlet-channel count and support-size
   breakdown for profiling scripts; the N=4 smoke passed in 27s, and the
   focused regression passed 2,718/2,718 tests in 3m03.8s.  The full
   contiguous-chain SU(2) structural target now reports the same
   singlet-channel accounting; the N=6 smoke passed in 11s, and the focused
   regression passed 2,720/2,720 tests in 3m04.5s.  The target-only profile
   script now prints these singlet-channel rows; an N=4 buffered output smoke
   completed and found the base/support rows in about 30s.  The
   translation-orbit zero-row emitter smoke passed in 20.9s, and the focused
   Pauli-chain regression passed 2,792/2,792 tests in 3m39.2s.  The public
   symbolic translation/SU(2) singlet-equality smoke passed in 30.3s, and the
   focused Pauli-chain regression passed 2,798/2,798 tests in 3m36.2s.  The
   reflected fixed-sector public smoke passed in 47.24s with 19 singlet
   equality zero rows, the all-sector reflected public smoke passed in 37.99s
   with 78 singlet equality zero rows, and the focused Pauli-chain regression
   passed 2,803/2,803 tests in 3m35.3s test time (392.90s remote wall).  The
   direct-linear reflected fixed-sector smoke passed in 54.13s with 19
   singlet equality zero rows, the direct-linear all-sector reflected smoke
   passed in 47.24s with 78 singlet equality zero rows, and the focused
   Pauli-chain regression passed 2,808/2,808 tests in 3m36.2s test time
   (392.75s remote wall).  `translation_report_metrics` now exposes
   `su2_singlet_channel_equality_row_count`; the N=6 metric smoke matched
   78/78 rows in 47.40s, and the focused Pauli-chain regression passed
   2,810/2,810 tests in 3m35.4s test time (390.88s remote wall).  The
   profile script now accepts `NCTS_TRANSLATION_SINGLET_EQUALITIES=true`; an
   N=4 no-solver smoke completed in 68.27s and printed 58 constructed
   singlet-channel equality zero rows.  The SOS zero-dual provenance path now
   preserves those labels: an N=6 no-solver dualization smoke found 78
   singlet equality zero duals among 1,182 zero duals in 40.41s, and the
   focused Pauli-chain regression passed 2,812/2,812 tests in 3m35.0s test
   time (378.87s remote wall).
   The base-SU2 plus extending SU(2)-RDM structural target now also accepts
   requested finite-axis equality rows instead of rejecting the combination:
   the red N=4/order-1 target smoke failed in 14.4s on the old guard, the
   green smoke passed with positive axis-row accounting in 17.7s test time
   (28.4s remote wall), and the focused Pauli-chain regression passed
   2,826/2,826 tests in 3m33.5s.
   The same target now accepts the Heisenberg `H^2` moment-equality smoke for
   base-SU2 plus extending SU(2)-RDM accounting: the red N=6/order-2 target
   smoke failed in 12.4s on the stale guard, the green target smoke passed in
   16.9s test time (24.5s remote wall), the constructed N=6/order-2 comparison
   matched block sizes and moment-equality row counts in 49.0s, and the
   focused Pauli-chain regression passed 2,835/2,835 tests in 3m34.9s.
   The public `cs_nctssos`/`SolverConfig` Pauli fast path now forwards
   `singlet_channel_equalities` and `singlet_channel_atol` into the
   translation/SU(2) direct-linear builder.  A high-level N=4 direct-linear
   solve preserved `MomentLinearData`, emitted singlet-channel equality rows
   with provenance, and matched the low-level direct SU(2) objective within
   `1e-6`; `test/relaxations/interface.jl` passed 107/107 checks in 152s
   remote wall, and the focused Pauli-chain regression passed 2,925/2,925
   checks in 3m21.2s test time (378s remote wall).

**Milestones:**
- 3A: reproduce existing order-2 singlet block sizes and bounds
- 3B: CG transform unitarity residual ≤ 1e-12, off-block residual ≤ 1e-10
- Overall: L=30 + SU(2) runtime ≤ 2× QMBCertify or ≥ 2× improvement over
  Phase 2

**Stop gate:** Only claim N=100 numerical parity after this phase succeeds
at L=30.

---

### Phase 4 — Certification provenance hooks

**Duration:** 1 week (hooks only; full certification is a separate project)

**Goal:** Ensure the fast path preserves enough metadata for future rigorous
certification.

**Design decision:** Certify in the **dual**, not by projecting primal moments.

**Current status:** The translation fast path records block labels, logical row
labels, DFT/reflection/RDM/SU(2) transform descriptors, zero-row origins, and
coefficient-domain metadata in `MomentLinearData`.  SOS dualization preserves
that metadata through `sos_dual_blocks`, `sos_dual_block_values`,
`sos_zero_duals`, `sos_zero_dual_values`, `sos_dual_block_diagnostics`, and
`sos_dual_certificate_residual`.  Small-N tests cover translation, reflection,
RDM, U(1)/SU(2) RDM, base translation/SU(2), base translation/SU(2) plus
axis-rotation equality rows, closed-support full/U(1)/SU(2) RDM add-ons,
linear state-opt, and PSD state-opt residuals to solver tolerance.  A small
N=4 provenance check for the newly promoted SU(2)+full/U(1)-RDM paths gave SOS
certificate residuals below `3e-14`.  The direct `MomentLinearData` path is
now also covered for the combined SU(2) moment + SU(2)-RDM + PSD-state-opt
profile: the test dualizes that fast-path model, checks the SOS certificate
residual, verifies diagnostic block labels and native cone sizes against the
translation report, and confirms both SU(2 RDM and PSD-state-opt block labels
survive extraction.  The focused Pauli-chain regression with this provenance
check passed all 2,609 checks through easy-ssh in 316.58s remote wall.
Zero-origin histograms work for both symbolic `MomentProblem` wrappers and
direct `MomentLinearData`.  Translation reports now also expose solve-support
status and blocker provenance, so unsupported reducer combinations are not
accidentally treated as solved certificate candidates.
SOS helpers now also expose `sos_zero_dual_diagnostics`, so zero/equality
dual multipliers carry origin labels, row kinds, coefficient-domain tags,
form sizes, coefficient magnitude, and numerical values alongside the PSD
block diagnostics.  A direct `MomentLinearData` smoke pins this public helper
on a solved zero-row model, including the `nothing` domain default for rows
without tagged labels.
PSD block diagnostics now surface `coefficient_domain` and
`exact_coefficient_domain` directly, instead of requiring certificate tooling
to unpack backend-specific transform records.  The Pauli-chain SOS diagnostics
regression checks these fields against the translation report for DFT and
reflection DFT blocks.
The base-SU2 plus extending SU(2)-RDM direct path also checks that report
coefficient-domain arrays are derived from the actual block transform metadata
stored on `MomentLinearData`, including reflected sectors.  The N=4 metadata
smoke completed in about 36s through easy-ssh, and the focused Pauli-chain
regression passed 2,671/2,671 tests in 3m03.5s.
Singlet-channel equality zero rows now preserve their specific
`PauliSU2SingletChannelEqualityOrigin` through the symbolic
translation/SU(2) realification bridge and through the direct-linear zero-row
helper.  A red N=6/order-2 reflected fixed-sector no-solver smoke failed on
the missing origin type in 44.32s; the green smoke preserved 19 singlet
origins and 19 report-counted singlet rows in 55.18s.  The focused
Pauli-chain regression then passed 2,815/2,815 tests in 3m33.1s test time and
387.29s remote wall time.
Those singlet-channel equality origins now also carry the transform
coefficient-domain metadata (`:complex_algebraic_float64` and
`:complex_sqrt_rational`) rather than only row/momentum labels.  A red
N=6/order-2 reflected fixed-sector smoke failed on the missing fields in
43.99s; the green smoke preserved 19 singlet origins with the expected domain
tags in 57.25s; and the focused Pauli-chain regression passed 2,816/2,816
tests in 3m38.0s test time and 395.98s remote wall time.
Zero-dual diagnostics now surface those same singlet-channel domain tags
directly from the origin label.  A red direct-helper regression failed on the
missing fields in about 46s; the green `moment_linear.jl` check passed in
about 102s, and an N=4/order-2 no-solver Pauli smoke confirmed the singlet
diagnostic domain tags in about 17s.
The lower-level `sos_zero_duals` handle records and `sos_zero_dual_values`
records now expose the same direct `label`, `coefficient_domain`, and
`exact_coefficient_domain` fields, so certificate extraction code does not need
to unpack origin internals for equality multipliers.  A red helper regression
failed on the missing fields in about 68s; the green `moment_linear.jl` check
passed in about 104s; and an N=6/order-2 no-solver Pauli smoke confirmed
singlet-channel zero-dual domain tags in about 34s.
SOS residual reports now also identify the checked coefficient count,
identity moment, worst residual moment, and worst residual value in addition to
the max norm and full residual dictionary.  A red direct-linear regression
failed on the missing fields in about 80s, and the green `moment_linear.jl`
check passed in about 102s across both real and Hermitian direct SOS residual
paths.  After the perf-harness guard hardening, the focused
`test/relaxations/sos.jl` suite was rerun through easy-ssh and passed 9/9
checks in 30.8s.  The focused `test/relaxations/dualization.jl` suite was
also rerun through easy-ssh and passed 3/3 moment-vs-SOS equivalence checks in
33.6s.
The source-like solver probe harness now computes and emits certificate
residual fields for solved SOS-dual rows as well.  A tiny N=4 dual/Mosek
emitter smoke reported max residual `7.66e-9` in about 89s end-to-end.  The
five pinned NCTSSoS source-like solver rows now carry loader-checked
`sos_certificate_*` fields and the `sos_dual_certificate_residual` source
helper: the N=12 A0/A1/A2 rows replayed under the N<14 guard with max
residuals `1.79e-9`, `9.64e-10`, and `9.50e-9`; the N=8 SU(2)-RDM plus
PSD-state-opt rows replayed with max residuals `5.28e-11` and `3.58e-10`.
The red fixture-loader check failed on the missing residual fields in 14.9s,
and the green loader check passed 759/759 checks in 13.0s.  A follow-up red
loader check failed on the missing residual-source provenance in 15.5s, and
the green loader check passed 764/764 checks in 12.9s.  After adding
CI-level perf guard coverage, the same loader passed 815/815 checks through
easy-ssh in about 54s, including 20/20 `Pauli Perf Harness Guards` and 8/8
`Perf Harness Entry Guards`.
After tightening the solver-probe construct-only env switch, the loader passed
823/823 checks through easy-ssh in about 50s, including 24/24
`Pauli Perf Harness Guards` and the same 8/8 `Perf Harness Entry Guards`.
`translation_linear_provenance(...)` now exposes the Phase-4 hook data already
carried by the translation fast path: PSD block labels, solver block sizes,
logical row labels, raw row labels, transform descriptors, coefficient-domain
tags, and zero/equality-row origins with term counts.  This is deliberately a
provenance hook, not interval certification.  The focused
`test/relaxations/pauli_chains.jl` file passed 2,877/2,877 checks through
easy-ssh after adding the helper tests; the main translation testset reported
3m11.9s, with remote wall time roughly 5-6 minutes including precompile and
monitor overhead.
SOS dual block value and diagnostic records now also expose
`transform_family` directly, matching the provenance helper and avoiding
another transform-internal unpacking step for certificate tooling.  The generic
direct `MomentLinearData` helper path pins `transform_family === nothing`, and
translation/reflection dual diagnostics pin the reported families against the
stored block transforms.  The focused `test/relaxations/moment_linear.jl` file
passed through easy-ssh in about 2m30s, and the focused
`test/relaxations/pauli_chains.jl` file passed 2,881/2,881 checks through
easy-ssh; the main translation testset reported 3m12.6s, with remote wall time
about 5m25s.
SOS zero-dual handle, value, and diagnostic records now expose direct
`feature`, `decomposition`, `reason`, and `term_count` fields, so equality
multiplier extraction no longer has to unpack origin-label internals for common
classification.  The direct helper tests pin the `nothing` defaults for
untagged rows, and the Pauli-chain provenance tests pin scalar equality plus
domain-tagged SU(2) singlet-channel equality rows.  The focused
`test/relaxations/moment_linear.jl` file passed through easy-ssh in about
2m30s, and the focused `test/relaxations/pauli_chains.jl` file passed
2,896/2,896 checks through easy-ssh; the main translation testset reported
3m12.8s, with remote wall time about 5m20s.
`sos_dual_certificate_diagnostics(...)` now bundles the numerical residual,
PSD block diagnostics, and zero/equality dual diagnostics into one public
record for future certificate extraction plumbing.  The direct helper test
pins the bundled residual/block/zero counts and worst-residual fields, while
the Pauli-chain scalar-equality smoke checks the same helper on a solved
translation SOS dual.  The focused `test/relaxations/moment_linear.jl` file
passed through easy-ssh in about 2m30s with the direct SOS helper testset at
55/55 checks, and the focused `test/relaxations/pauli_chains.jl` file passed
2,907/2,907 checks through easy-ssh; the main translation testset reported
3m11.4s, with remote wall time about 5m20s.
The solver-probe fixture emitter now also writes the combined certificate
diagnostics source and summary fields for solved SOS-dual rows.  An empty
zero-dual base case first failed with a `Vector{Any}` dispatch hole in
`sos_zero_dual_diagnostics`; the green focused `moment_linear.jl` rerun passed
in about 2m07s after broadening the public helper to accept empty diagnostic
vectors.  A tiny N=4/order-1 dual/Mosek fixture smoke then emitted
`sos_certificate_diagnostics_source`, PSD/zero-dual counts, native Hermitian
residual, minimum native eigenvalue, and maximum zero-dual value with
`OPTIMAL` status and max SOS residual `6.91e-13` in about 70s end-to-end.
The five pinned source-like solver-probe rows were then replayed on
2026-07-02 under the N<14 guard and refreshed with the same diagnostics fields
under Julia 1.12.6 and the emitted MOSEK 11.2.2 runtime.  The A0/A1/A2 L=12
rows completed in 85.35s, 111.42s, and 109.45s remote wall time with
residuals `1.79e-9`, `9.64e-10`, and `9.50e-9`; the two N=8
SU(2)-RDM+PSD-state-opt rows completed in 105.59s and 108.30s with residuals
`5.28e-11` and `3.58e-10`.  The refreshed loader now checks the diagnostics
source, PSD and zero-dual counts, native Hermitian residuals, native minimum
eigenvalues, and zero-dual magnitudes; it passed 806/806 fixture checks
through easy-ssh in about 52s.
The solver-probe environment emitter now derives the MOSEK runtime version
from `Mosek.getversion()` instead of package dependency metadata, so future
replay rows no longer need a manually patched `mosek_version`.  A red harness
guard first failed on the missing helper, the green loader passed 806/806
fixture checks in about 52s, and a tiny N=4/order-1 construct-only
Mosek-emitter smoke verified that the real environment fixture now emits
`mosek_version = "11.2.2"` in about 61s without a solver call.  After updating
the refreshed NCTSSoS solver-probe environment to that emitted value, the
loader passed 806/806 checks again in about 47s.
The solver-probe fixture emitter now also marks row execution state explicitly
with `construction_only`, `model_built`, and `solved`, so construct-only rows
do not have to be inferred from missing lowering/dualization/solve timing
fields.  The expectation-loader source guard rose to 14/14 solver-probe
fixture-field checks after adding the SU(2) base Wigner zero-row split fields,
and the loader passed 993/993 fixture checks through easy-ssh in about 60s.
A tiny N=4/order-1 construct-only fixture smoke printed
`termination_status = "not_solved"`, `construction_only = true`,
`model_built = false`, `solved = false`, `linear_moments = 7`, and
`psd_cones = 13` in 73s remote wall.  A second tiny base-SU(2)+extended-RDM
construct-only smoke printed the nonzero Wigner split
`su2_base_zero_row_count = 34`,
`su2_base_spin_offblock_row_count = 4`,
`su2_base_magnetic_offdiag_row_count = 23`, and
`su2_base_magnetic_copy_row_count = 7` in about 67s remote wall.  The same
fixture emitter now serializes report construction-stage timings as sorted
`construction_stage_*_seconds` TOML fields; the helper-level guard raised
Pauli perf harness checks to 61/61, and a tiny N=4/order-1 construct-only smoke
emitted stage fields including `construction_stage_axis_diagnostics_seconds`,
`construction_stage_block_assembly_seconds`, and
`construction_stage_linearization_seconds` in about 60s remote wall.
The translation-profile harness now carries report stage timings into repeated
run summaries as a bounded `Repeat Stage Summary` table.  The red guard failed
on the absent table, the green loader passed 993/993 fixture checks in about
61s, and a tiny N=4/order-1 SU(2)+extended-RDM no-solver repeat smoke printed
the new stage-range table in 78s remote wall.
The axis-quotient rewrite now preserves `trusted_self_adjoint` on rewritten
zero constraints instead of silently dropping that provenance bit.  A synthetic
builder regression pins the rewrite behavior, and the focused Pauli-chain file
passed 2,919/2,919 checks through easy-ssh; the main translation testset
reported 3m17.2s.
A follow-up touched-suite verification sweep on 2026-07-02 kept all runs below
the shared-machine `N < 14` guard and passed the current shared regression
surface through easy-ssh: `test/expectations_loader.jl` passed 993/993 checks
in 58.387s, `test/relaxations/pauli_chains.jl` passed 2,926/2,926 checks in
6m11.695s, `test/relaxations/moment_linear.jl` passed in 2m8.883s,
`test/relaxations/interface.jl` passed in 2m20.744s,
`test/relaxations/lowering.jl` passed in 29.647s, `test/relaxations/sos.jl`
passed in 36.875s, `test/relaxations/dualization.jl` passed in 42.308s,
`test/relaxations/sparsity.jl` passed in 30.750s,
`test/relaxations/symmetry.jl` passed in 2m23.550s, and
`test/relaxations/gns_pipeline.jl` passed in 46.204s.  This is strong
small-N implementation evidence; the full `test/relaxations/runtests.jl`
entrypoint then passed 7,941/7,941 checks in 9m12.395s.  These regressions are
not a substitute for the deferred L=30/N=100 parity gates.  A package-level
`Pkg.test()` run through easy-ssh also passed with 12,075 checks and 2 broken
tests in 23m4.475s remote wall after Julia resolved the test-target
environment.  A direct `test/runtests.jl` invocation is not the right package
gate here because it omits test-target extras such as Aqua, Documenter,
ExplicitImports, and COSMO from the active environment.
Completion is no longer blocked by a blanket `N >= 14` ban, but every larger
run now needs an explicit stall-risk check before launch.  The first cautious
large-run pass on 2026-07-06 found that QMBCertify A1_L20 and A1_L30 are safe
and reviewed on autodl with one thread, while A0_L20 and A2_L20 are safe but
not reviewed because MOSEK returned `SLOW_PROGRESS`/`FEASIBLE_POINT`.  The
matching NCTSSoS finite-axis A1/A2 construction-only probes at L=20 and L=30
also completed safely under one-thread hard guards, with L=30 peak RSS below
7 GiB and remote wall below 8.5 minutes.  A2_L20 and A2_L30 no-solve
SOS-dual model-sizing also completed safely, and A1_L20/A2_L20 solved with
Mosek inside 15-minute guards.  The reviewed-reference A1_L30 NCTSSoS solve
also completed safely inside a 25-minute guard, but A1 failed the bound-quality
parity target at both L=20 and L=30: L=20 objective `-9.092707010191038`
versus QMBCertify total `-8.904954518917137`, and L=30 objective
`-13.664267269723961` versus QMBCertify total `-13.325070793538082`.  A
follow-up formulation audit at L=20 showed that the old source-like NCTSSoS
linear state-opt rows were not QMBCertify rows: QMBCertify generates 1817
candidate monomials, keeps 198 filtered row images, while the old contiguous
width-7 path emitted 810 rows with only 46 row images overlapping after
QMBCertify-style reduction.  The new explicit
`linear_state_opt_mode=:qmbcertify` path uses the source candidate family and
support-filter semantics; it reduces the L=20 NCTSSoS row count to 408 and
solves safely with one Julia thread and one Mosek thread in 145.71s solve time
after 160.95s construction, peak RSS 3.90 GiB.  The objective remains
`-9.092707064688685`, essentially unchanged, so the A1 bound-quality gap is
not caused primarily by the contiguous LSO row family.  The next diagnostic
mismatch was the base basis/block layer: QMBCertify A1_L20 had base block max
114, while the then-active finite-axis NCTSSoS base reported 18 solver-facing
blocks of size 120.  This is historical diagnostic evidence; the later
source-base constructor closes this no-solver base histogram/support-quotient
gap, so it is not a reason to spend more L=20/L=30 solver time by itself.
A second L=20 no-solver basis audit made the diagnostic concrete.  QMBCertify's
A1 basis has 105 orbit families, with word-length histogram
`[1 => 1, 2 => 50, 3 => 13, 4 => 41]`; only 60 of those families lie in the
NCTSSoS contiguous order-4 representative set, 45 QMBCertify families are
outside it, and 61 NCTSSoS contiguous representatives are absent from
QMBCertify's base basis.  Three QMBCertify families have translation orbit
length 10 at L=20 (`σx₁σx₁₁`, `σy₁σy₁₁`, `σz₁σz₁₁`), while the current
analytic translation backend explicitly supports only the identity orbit and
full-length translation orbits.  Passing the exact QMBCertify union basis to
the NCTSSoS backend therefore fails before model construction with that
short-orbit guard.  The audit completed safely in 4m36.74s wall with 3.67 GiB
peak RSS, one Julia thread, no model lowering, and no solver call.
The QMBCertify chain basis is now a first-class structural target through
`pauli_qmbcertify_chain_basis(ops, order; extra=...)`.  The focused regression
pins the L=20/order-4/`extra=9` source basis lengths `[1140, 960]`, family
counts `[57, 48]`, unique basis count 2051, base block sizes
`[58, 114×9, 57, 48, 96×9, 48]`, word-length histogram
`[1 => 1, 2 => 50, 3 => 13, 4 => 41]`, and translation-orbit histogram
`[10 => 3, 20 => 102]`.  The focused Pauli-chain regression passed 2,941/2,941
checks through easy-ssh in 3m57.2s test time.
The same helper now reports overcomplete row accounting: 2100 source rows,
2050 unique non-identity rows, 50 duplicate source rows, unique row counts
`[1110, 940]` by parity, duplicate row counts `[30, 20]` by parity, three
short-orbit families, and 30 rows coming specifically from short-orbit
overcompletion.  The focused Pauli-chain regression passed 2,948/2,948 checks
through easy-ssh in 3m50.1s test time after this accounting was added.
The helper also pins the QMBCertify base-block metrics directly: 22 base
blocks, max block 114, total block side 2101, dense entries 211,129,
symmetric entries 106,615, and histogram
`[48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]`.  The focused Pauli-chain
regression passed 2,954/2,954 checks through easy-ssh in 3m56.4s main-testset
time after these metrics were added.  The larger-run gate estimated 4-7
minutes and <3.5 GiB RSS; actual monitored wall was about seven minutes,
with peak sampled RSS 2.80 GiB, so the run stayed inside the safe envelope.
The base helper now also exposes explicit QMBCertify source block descriptors:
parity label, momentum, family count, solver block size, identity-row flag,
and realification flag for all 22 blocks.  A red/green structural check through
easy-ssh confirmed the new field on the L=20/order-4/`extra=9` target; the
green run took about 20s including precompilation and sampled at 587 MiB RSS.
The broader focused Pauli-chain file was not rerun immediately because the
fresh load gate reported `nproc=25` with load average `45.09, 41.87, 44.02`.
That deferral was later cleared by the current focused Pauli-chain run: after
the load dropped, the suite passed 3,051/3,051 checks in 8:09.27 wall with
4,996,112 KiB peak RSS and no swaps.
`pauli_translation_structural_targets` now accepts `qmbcertify_base_extra` as a
target-only reference knob.  For L=20/order-4/`extra=9`, it reports the same
QMBCertify base reference block sizes and histogram while leaving
`qmbcertify_base_reference_active=false`, so the report does not pretend the
active contiguous-base construction has been replaced.  The red/green
structural-target path now also accepts `qmbcertify_base_construct=true` as a
no-solver active-source-base target.  In that mode, L=20/order-4/`extra=9`
reports source-row basis size 2101, orbit-family size 106, active
`qmbcertify_base_reference_active=true`, and active PSD block histogram
`[48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]`.  A RED guard failed in
10.77s wall because the keyword was missing; the GREEN guard passed under a
60s timeout after adding the target switch, and a follow-up no-solver
reference-vs-active guard passed under a 45s timeout.  Both checks used one
Julia thread and no solver.  The fresh load gate was still overcommitted at
28.36 on 25 CPUs, so the full focused Pauli-chain file was deferred for that
pass; the later 3,051-check focused Pauli-chain run above covers the current
target path.
`perf/pauli_translation_profile.jl` now allows
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true` in target-only mode and
passes it through to the active source-base structural target.  The stale RED
source guard failed because the profile still required
`NCTS_TRANSLATION_TARGET_ONLY=false`; the GREEN source guard passed after
removing that restriction and threading the keyword into `_print_target_only`.
The actual L=20/order-4/`extra=9` target-only profile then passed through
easy-ssh with `NCTS_TRANSLATION_ALLOW_LARGE=true`, one Julia thread, no solver,
and a 90s timeout.  It printed basis size 2101, orbit basis size 106,
22 active QMBCertify base PSD blocks, max block 114, symmetric entries
106,615, and active construction `true`.  The first rerun without
`NCTS_TRANSLATION_ALLOW_LARGE=true` stopped on the intended size guard, not on
the removed target-only restriction.
The active source-base structural target now also supports exact target-only
RDM and PSD-state-opt block-shape accounting.  RDM add-ons are accepted only
with `contiguous_rdm_decomposition=:qmbcertify`,
`contiguous_rdm_support=:extend`, and `real_moment_matrix=true`; linear
state-opt remains construction-only because its source row count is filtered
against the actual moment basis.  A RED no-solver guard failed in 24.89s wall
because the active source-base target still rejected RDM add-ons.  The GREEN
guard passed under a 90s timeout after adding the RDM/PSO target formulas:
for L=20/order-4/`extra=9`, RDM8 reports PSD blocks `[72, 64, 56]`, PSO3
reports PSD blocks `[36, 72×9, 36, 28, 56×9, 28]`, and the source-family PSO
candidate count is 64.  A target-only profile with active source base, RDM8,
and PSO3 then passed through easy-ssh with one Julia thread, no solver, and a
120s timeout.  It printed 47 total PSD blocks, max block 114, symmetric
entries 153,079, 3 RDM blocks, 22 PSO blocks, and zero known add-on equality
rows.  The fresh pressure gate before these checks was `44.75` load on
25 CPUs, so the full focused Pauli-chain file was deferred for that pass; the
later 3,051-check focused Pauli-chain run above covers the current target path.
The target-only profile now mirrors the solver-probe defaults for the
source-base RDM path: with
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true`, an omitted
`NCTS_TRANSLATION_RDM_DECOMPOSITION` defaults to `qmbcertify` and an omitted
`NCTS_TRANSLATION_RDM_SUPPORT` defaults to `extend`.  A RED source guard failed
before those defaults existed and passed after the profile parser was changed.
The behavioral L=20/order-4/`extra=9` target-only profile then ran without
setting either RDM env var, under a preflight load of `32.29` on 25 CPUs.  The
run was estimated at 45-75s and <2.5 GiB RSS; actual `/usr/bin/time -v` wall
was 45.93s with 2.34 GiB max RSS, one Julia thread, no model construction, and
no solver.  It printed `contiguous_rdm_decomposition=qmbcertify`,
`contiguous_rdm_support=extend`, and RDM blocks `[72, 64, 56]`.  This remains
a small target-only ergonomics check; broad suites stay deleted from the live
queue while the remote load exceeds its CPU count.
The easy-ssh check for this target-only field passed with exit code 0, sampled
RSS 0.85 GiB, and no solver call under the same high-load gate.
`perf/pauli_translation_profile.jl` now exposes this reference through
`NCTS_TRANSLATION_QMBCERTIFY_BASE_EXTRA`; the target-only L=20/order-4 profile
prints a "QMBCertify base reference" section with family counts `[57, 48]`,
22 blocks, max block 114, dense entries 211,129, symmetric entries 106,615,
histogram `[48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]`, and
`active construction=false`.  The red profile-output check failed before this
section existed; the green check passed through easy-ssh with exit code 0,
sampled RSS about 2.0 GiB, one Julia thread, no model construction, and no
solver call.
The source-base support inventory is now pinned as well.  Mirroring the
QMBCertify chain `tsupp` collection for the L=20/order-4/`extra=9` base yields
37,298 nonzero source product rows, 18,232 sign-zero rows, 958 diagonal
nonzero rows, 36,340 off-diagonal nonzero rows, and 3,420 unique reduced
support words with length histogram `[0 => 1, 2 => 10, 4 => 575, 6 => 1816,
8 => 1018]`.  These fields are available both from
`pauli_qmbcertify_chain_basis` and from the
`pauli_translation_structural_targets(...; qmbcertify_base_extra=9)` reference
payload.  Red/green easy-ssh checks passed for both entry points with one
Julia thread, no model construction, and no solver call; the structural-target
green check sampled 1.2 GiB RSS.  The target-only profile now prints the same
support inventory in its "QMBCertify base reference" section, including
37,298 nonzero rows, 18,232 zero rows, 3,420 unique support words, and the
word-length histogram.  The profile-output RED check failed before those lines
existed; the GREEN check passed through easy-ssh in 46.05s wall time with
2.09 GiB max RSS under `JULIA_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, a
180s timeout, and no solver/model construction.  The active broad-test and
constructed-scaling runs remain deleted from the live run plan unless a fresh
load gate shows they will not crowd the remote cgroup quota.
The first no-solver QMBCertify source-base constructor is now active as an
internal direct-linear builder.  It uses the source family rows, the
short-orbit/overcomplete block accounting, and the same `reduce4` support
quotient used by the source audit rather than ordinary translation orbit
representatives.  A TDD RED check failed on the missing constructor in 12.56s
with 800 MiB RSS.  The first GREEN attempt correctly exposed a bug: fixed
momentum `k=0` was generically realified to 116 rows instead of QMBCertify's
58-row real block.  After projecting fixed momenta to their real part and
realifying only conjugate-pair sectors, the L=20/order-4/`extra=9` no-solver
construction passed through easy-ssh in 34.15s wall time with 919 MiB max RSS:
22 PSD blocks, max block 114, 3,420 linear moments, zero scalar constraints,
and product-cache hit rate 0.9074.  This run built the source-base linear PSD
blocks but still did not build a JuMP model or call a solver.
`perf/pauli_translation_profile.jl` now exposes the same path with
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true`.  A small N=4 RED profile
check showed that the flag was previously ignored and fell through the old
path, taking 78.42s and 1.96 GiB RSS before the expected grep failure.  The
GREEN N=4 profile route passed after adding the flag and also caught a small-N
edge case: conjugate-pair sectors must be fully realified even when the entries
happen to be real, because QMBCertify's source block rule still doubles those
rows.  The repeatable L=20/order-4/`extra=9` profile command then passed with
`NCTS_TRANSLATION_ALLOW_LARGE=true`,
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true`, and
`NCTS_TRANSLATION_QMBCERTIFY_BASE_EXTRA=9`: the profile reported basis size
2101, orbit/family basis size 106, 22 PSD blocks, max block 114, 3,420 moments,
zero scalar constraints, and 0.907400 product-cache hit rate.  The constructor
phase itself took 3.874s and allocated 1.75 GiB; the full profiled process took
85.92s wall with 2.26 GiB max RSS.  This is repeatable no-solver construction
evidence, not bound-quality or solver evidence.

Conclusion: the source-base no-solver construction now matches QMBCertify's
base block histogram and support quotient, and the public solve wrapper can
route a QMBCertify-base solve through the source-base linear data with
`qmbcertify_base_construct=true`.  The TDD RED wrapper check failed in
8.80s with 795 MiB RSS because the new keyword was forwarded to the old
constructor.  After adding the explicit route, an N=4/order-1 Clarabel smoke
passed through easy-ssh in 83.29s with 1.98 GiB max RSS: objective
`-1.9999999996`, six linear moments, and PSD block sizes `[4, 6, 3, 3, 6, 3]`.
A guarded dense-vs-QMBCertify-base comparison at the same size then took
94.29s and 2.24 GiB RSS and deliberately failed the parity assertion:
QMBCertify base gave `-1.9999999996` while the dense order-1 relaxation gave
`-2.9999999946`.  This is correct negative evidence for this slice: the
source-base route is a construction/solve route, not a bound-quality parity
route.  PSD state-optimality add-ons are still source-like rather than fully
integrated with the new QMBCertify base builder.
The same source-base solve route is now exposed through the generic
`cs_nctssos(pop, SolverConfig(...); qmbcertify_base_construct=true)` Pauli
translation bridge.  A RED full-file Pauli-chain regression first exposed that
the `SolverConfig` wrapper leaked profile-derived `sign_symmetry`,
`reflection_symmetry`, `axis_rotation_symmetry`, and `check_invariance`
keywords into the source-base constructor; it failed in 7:47.75 wall with
4,613,776 KiB peak RSS and no swaps.  After threading the QMBCertify source-base
keywords through the shared fast-path option bundle and suppressing those
profile-derived keywords on the source-base route, the focused Pauli-chain suite
passed in 7:56.86 wall with 4,648,664 KiB peak RSS and no swaps.  The broader
interface suite also passed in 2:27.00 wall with 2,057,500 KiB peak RSS and no
swaps.  No L>=14 solve, L=20 replay, or L=30 replay was used for this bridge
evidence.
The bridge now also rejects unsupported source-base fast-path options before
construction; a small regression pins that `momenta` cannot be combined with
`qmbcertify_base_construct=true`, and the lower constructor's unsupported-key
message lists the supported source-base knobs accurately.  The focused
Pauli-chain suite passed after this guard in 8:09.27 wall with 4,996,112 KiB
peak RSS and no swaps.
The source-base constructor now also accepts scalar equality rows that are
compatible with QMBCertify's support quotient.  The RED focused Pauli-chain run
failed on the old equality rejection in 6:32.30 wall with 2,965,416 KiB peak
RSS and no swaps.  The first GREEN attempt deliberately exposed the reducer
boundary: generic translation reduction kept `σyσy` and `σzσz` equality
moments outside the source-base support and failed in 6:46.93 wall with
3,212,800 KiB peak RSS and no swaps.  After switching source-base scalar
equality rows to the same QMBCertify support reduction used for the objective,
the targeted L=20/order-4/`extra=9` no-solver constructor smoke passed in
32.17s with 991,588 KiB peak RSS and no swaps, and the focused Pauli-chain
suite passed 3,055/3,055 checks in 7:56.65 wall with 4,612,860 KiB peak RSS and
no swaps.
The source-base constructor now also accepts covered scalar inequality rows as
1x1 localizing PSD blocks using the same support reduction.  The RED targeted
constructor run failed on the old inequality rejection in 18.99s with
1,040,160 KiB peak RSS and no swaps.  The GREEN targeted constructor smoke
passed in 31.16s with 984,800 KiB peak RSS and no swaps, adding one PSD block
and no new moments.  The focused Pauli-chain suite then passed 3,061/3,061
checks in 8:03.91 wall with 4,784,408 KiB peak RSS and no swaps.
Moment-equality constraints remain intentionally rejected by this source-base
constructor until a separate row-basis formulation slice earns them.
A no-solver closure probe confirmed that the naive source-family row basis is
not sufficient for `H^2` moment-equality rows: for L=20/order-4/`extra=9`, the
106 candidate source rows reduced through the QMBCertify support quotient gave
1 covered row, 54 zero rows, and 51 rows with moments outside the 3,420-moment
source-base support; the first missing row was `σx₁σx₂`.  The probe itself ran
in 25.39s with 1,070,320 KiB peak RSS and no swaps.  This deletes the naive
"reuse source-family rows" moment-equality implementation from the live plan.

Do not spend more L=20 or L=30 solver time on the old base formulation or on
the base-only QMBCertify route to chase the A1 bound gap.  Current solver
evidence tops out at N=4, where even the tiny checks use roughly 2 GiB RSS and
80-95s wall time under high remote load.  The active L=20/L=30 solver parity
runs are therefore deleted from the live run plan until no-solver formulation
checks show that the missing add-on rows are represented and a fresh load/RSS
gate predicts the run will not crowd the remote cgroup.
The previous absolute "do not try N >= 14" cap is replaced by this load-gated
policy only: every larger candidate run must state the scientific question,
expected wall/RSS/load risk, and a stop/delete condition before launch.  If the
evidence predicts a stall or cgroup crowding, the run is removed or downsized
instead of queued.
The first state-opt slice is now in place for linear state-opt rows.  A TDD RED
check for `_pauli_qmbcertify_chain_base_linear_relaxation(...;
linear_state_opt_width=7)` failed in 8.29s with 792 MiB RSS because the
constructor did not accept the keyword.  The first implementation attempt reused
the ordinary translation-reduced LSO helper and correctly failed: at
L=20/order-4/`extra=9`, ordinary translation reduction had 1,817 unique
nonzero qmb LSO rows but covered zero source-base moments.  The root cause was
the reduction mismatch.  Rechecking the same rows with QMBCertify's support
quotient gave 1,380 unique nonzero reduced rows and 340 rows covered by the
source-base moment set.  The GREEN no-solver constructor check passed through
easy-ssh in 35.44s with 1.04 GiB max RSS: 340 qmbcertify-mode LSO rows,
3,420 linear moments, and max PSD block 114.  The profile route was then
verified with
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true`,
`NCTS_TRANSLATION_QMBCERTIFY_BASE_EXTRA=9`, and
`NCTS_TRANSLATION_LSO_WIDTH=7`; it reported 340 linear zero constraints,
1.116s internal `linear_state_opt` time, 4.974s report construction time,
81.61s process wall time, and 2.40 GiB max RSS.  This remains no-solver
formulation evidence, not bound-quality parity evidence.
The source-base structural target now counts qmbcertify-mode LSO rows without
building the linear constructor.  A RED target probe failed on the old
"do not estimate linear state-opt rows" guard in 19.14s with 1,356,440 KiB peak
RSS and no swaps.  The GREEN target-only probe reported 1,817 qmbcertify LSO
candidates, 170 covered sign-canonical rows, 22 base PSD blocks, and 3,420 base
support moments in 36.25s with 1,382,572 KiB peak RSS and no swaps.  The
focused Pauli-chain suite then passed 3,066/3,066 checks in 8:08.90 wall with
4,779,496 KiB peak RSS and no swaps.  This replaces the previous need to build
the full L=20 source-base constructor merely to know the LSO target count.
The first source-base RDM slice is now in place as well.  A TDD RED check for
`_pauli_qmbcertify_chain_base_linear_relaxation(...; contiguous_rdm_k=8)`
failed in 21.01s with 808 MiB RSS because the constructor did not accept the
keyword.  Before the GREEN run, an easy-ssh load smoke compiled the changed
package in 17.63s with 776 MiB max RSS.  The L=20/order-4/`extra=9` no-solver
constructor check then passed in 25.11s with 1.07 GiB max RSS: 3,534 linear
moments, 25 PSD blocks, QMBCertify RDM block sizes `[72, 64, 56]`, and max PSD
block 114.  The no-solver profile route was verified with
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true`,
`NCTS_TRANSLATION_QMBCERTIFY_BASE_EXTRA=9`,
`NCTS_TRANSLATION_RDM_K=8`,
`NCTS_TRANSLATION_RDM_DECOMPOSITION=qmbcertify`, and
`NCTS_TRANSLATION_RDM_SUPPORT=extend`; it reported 5.65s internal construction,
1.80s internal `contiguous_rdm` time, 89.28s process wall time, and 2.21 GiB
max RSS.  This is no-solver RDM formulation evidence only; it does not reopen
the L=20/L=30 solver queue.
The source-base PSD state-opt slice is now in place too.  A TDD RED check for
`_pauli_qmbcertify_chain_base_linear_relaxation(...; psd_state_opt_width=3)`
failed in 6.34s with 681 MiB RSS because the constructor did not accept the
keyword.  The GREEN no-solver constructor check passed in 128.23s with
1.46 GiB max RSS: 4,980 linear moments, 44 PSD blocks, source-family
QMBCertify PSO block sizes
`[36, 72, 72, 72, 72, 72, 72, 72, 72, 72, 36, 28, 56, 56, 56, 56, 56, 56, 56, 56, 56, 28]`,
and max PSD block 114.  The first no-solver profile route was verified with
`NCTS_TRANSLATION_QMBCERTIFY_BASE_CONSTRUCT=true`,
`NCTS_TRANSLATION_QMBCERTIFY_BASE_EXTRA=9`, and
`NCTS_TRANSLATION_PSO_WIDTH=3`; it reported 103.85s internal construction,
99.25s internal `psd_state_opt` time, 188.55s process wall time, and 2.50 GiB
max RSS.  Total allocation was high (68.19 GiB), so this is acceptable
no-solver formulation evidence but not a reason to start combined
RDM+LSO+PSO or solver-backed L=20/L=30 runs on the shared machine.
The PSO hot loop then gained a source-reduced term cache shared across
momentum sectors.  A RED symbol check failed in 4.85s with 700 MiB RSS before
the cache existed; the GREEN cached-vs-uncached block equivalence check passed
in 23.15s with 920 MiB RSS.  The standalone PSO constructor then improved to
36.59s wall, 1.37 GiB RSS, and 12.45s internal `psd_state_opt` time, preserving
4,980 linear moments, 44 PSD blocks, and max PSD block 114.
The full source-base combined shape is now validated without a solver:
`contiguous_rdm_k=8`, `linear_state_opt_width=7`, and
`psd_state_opt_width=3` passed in 41.51s with 1.65 GiB max RSS.  It produced
5,034 linear moments, 47 PSD blocks, QMBCertify RDM blocks `[72, 64, 56]`,
the source-family PSO blocks above, 672 qmbcertify-mode LSO rows, and max PSD
block 114.  The matching no-solver profile route reported 18.72s internal
construction, 10.72s `psd_state_opt`, 1.57s `linear_state_opt`, 1.86s
`contiguous_rdm`, 102.85s process wall time, and 2.28 GiB max RSS.  This
reopens only further no-solver formulation audits; L=20/L=30 solver parity
runs remain deleted until a separate solver-risk gate is justified.
The target-only source-base structural counter then exposed a narrower
accounting bug: the combined RDM+LSO+PSO target initially counted only 170
sign-canonical LSO rows because it checked LSO coverage against the base
support, not the support extended by the RDM and PSD state-opt blocks.  The
counter now builds the same extended support before counting qmbcertify-mode
LSO rows.  The GREEN target-only rerun reported 336 LSO rows, 1,817 LSO
candidates, 64 PSO candidates, 47 PSD blocks, max block 114, and RDM blocks
`[72, 64, 56]` in 52.31s wall with 1,313,800 KiB peak RSS and no swaps.  This
matched the 35-70s / <2.5 GiB estimate and used no solver or model lowering.
The focused Pauli-chain regression suite then passed 3,074/3,074 checks in
8:28.45 wall with 4,867,580 KiB peak RSS and no swaps.  The narrow
expectation-loader verifier passed in 1:10.40 wall with 2,403,484 KiB peak RSS
and no swaps, including the parity-plan run gates and 2,708 TOML fixture
checks.
The qmb source-base route now also has a small Phase-4 SOS dual provenance
regression.  The N=4/order-1 qmb base smoke dualized the source-base
`MomentLinearData`, solved the SOS-dual helper model, and confirmed six
`qmbcertify_base` dual blocks with `:translation_dft` transform family,
`:cyclotomic_float64` coefficient domains, `:cyclotomic` exact domains, and
certificate residual `4.891920202254596e-16`.  The targeted smoke ran in
1:08.50 wall with 2,073,636 KiB peak RSS and no swaps.  The focused
Pauli-chain regression then passed 3,085/3,085 checks in 8:28.87 wall with
4,856,120 KiB peak RSS and no swaps.
The same source-base SOS hook now covers scalar constraint provenance.  A tiny
N=4/order-1 qmb source-base smoke with one scalar equality and one scalar
inequality solved the SOS-dual helper model with max certificate residual
`1.2495481604722566e-15`, preserved the scalar-equality zero-dual label, and
preserved the scalar-inequality PSD-block label.  It ran in 34.98s wall with
1,049,836 KiB peak RSS and no swaps.  The focused Pauli-chain regression then
passed 3,096/3,096 checks in 8:30.17 wall with 4,761,924 KiB peak RSS and no
swaps.
The no-solver L=20 combined qmb add-on regression now also checks
`translation_linear_provenance` directly: PSD block records must preserve the
qmbcertify base, contiguous-RDM, and PSD state-opt labels/sizes/logical row
labels, while zero-constraint records must preserve qmbcertify LSO labels.  The
focused Pauli-chain regression passed 3,104/3,104 checks in 8:31.09 wall with
4,816,436 KiB peak RSS and no swaps.  This matched the 8-10 min / ~5 GiB
estimate on `autodl` after the pre-run load gate reported load averages near
10-12 and 653 GiB available RAM.
`perf/pauli_translation_solver_probe.jl` now has the same source-base opt-in
for construction-only and future model-sizing rows:
`NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true`,
`NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_EXTRA`, and
`NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_THREE_TYPE`.  A RED construction-only
probe with `N=6` failed as intended because the flag was ignored and the old
direct-linear path ran; the run took 70.36s wall with 2.33 GiB peak RSS.  After
wiring the source-base route, `N=6` and `N=8` source-base probes were removed
as validation sizes because the QMBCertify source family generates short-orbit
words such as `[2, 3]` there, which simplify with a non-real Pauli phase and
fail before useful shape output.  The valid L=20 construction-only solver
probe then passed with one thread, no lowering, and no solver: 79.46s process
wall, 2.31 GiB peak RSS, 3.98s internal construction, 3,420 linear moments,
22 PSD blocks, max block 114, and product-cache hit rate 0.9074.  This keeps
source-base solver-probe construction rows in scope, but it does not justify
combined A2 model lowering or any L=20/L=30 solver run.  Those runs remain
deleted until a smaller no-solve model-sizing rung is explicitly estimated and
passes a fresh load/RSS gate.
The solver probe now prints static model-size estimates immediately after
construction and before any JuMP lowering:
`estimated_model_variables`, `estimated_psd_cone_scalar_variables`,
`estimated_scalar_equalities_upper_bound`, `estimated_zero_dual_variables`,
`estimated_dense_schur_bytes`, `estimated_free_orphan_variables`,
`estimated_aux_orphan_blocks`, `estimated_lowering_would_error`, and
`estimated_risk_gate_status`.  The probe also prints
`estimated_model_size_gate_status`, computed from the dense-Schur proxy and the
same memory-headroom fraction used by the large-run guard.  It now also emits
the model-size gate inputs as fixture fields:
`estimated_model_size_gate_estimated_rss_bytes`,
`estimated_model_size_gate_mem_available_bytes`, and
`estimated_model_size_gate_max_rss_fraction`, so a copied no-solve row carries
the memory evidence behind the status.  The existing source-base model-sizing
rows were promoted with review-time `MemAvailable` telemetry
`746746384384` bytes and the default `0.8` headroom fraction.  A RED TOML guard
failed on the missing fields in 3.24s with 428,600 KiB RSS; the GREEN
source/TOML guard passed in 1.50s with 339,364 KiB RSS.  A single harness-include
helper smoke, with no model construction or solver call, checked the gate tuple
in 5.98s with 715,376 KiB RSS.  The gate tuple also carries
`estimated_model_size_gate_reason`; for `ok` fixture rows this is empty, while
blocked rows explain the memory-headroom failure.  A RED TOML guard failed on
the missing reason field in 3.25s with 428,200 KiB RSS; the GREEN source/TOML
guard passed in 1.50s with 337,692 KiB RSS, and a helper smoke checked the
blocked-headroom reason in 5.93s with 714,656 KiB RSS.  The
explicit risk-gate status is `"blocked_lowering_orphan_policy"` when the
selected primal lowering would error on free orphan keys, and `"ok"` otherwise.
The estimator is intentionally a
pre-launch risk gate rather than a replacement for actual JuMP counts.  A RED
loader check failed on the missing helper in 67.99s with 2.35 GiB peak RSS;
after implementation, the loader passed 1,778 checks in 69.64s with 2.34 GiB
peak RSS.  The L=20/order-4/`extra=9` source-base construction-only probe then
printed the new estimates without lowering or solving: 106,615 estimated
primal PSD-block variables, 106,616 scalar binding/equality rows as an upper
bound, and `estimated_lowering_would_error=true` under the default
`orphan_policy=:error` because the source-base construction has 3,391 free
keys.  That is explicit negative evidence against a blind primal model-lowering
run; any future model-sizing run must first choose and document an orphan
policy or a dualization route whose static estimate clears the load/RSS gate.
The probe itself completed in 77.36s wall with 2.32 GiB peak RSS.  A matching
construction-only SOS-dual estimate, still with `NCTS_SOLVER_PROBE_LOWER_MODEL=false`
and no solver, completed in 80.13s wall with 2.32 GiB peak RSS and reported
`estimated_model_mode=sos_dual`, 106,615 estimated dual variables, 3,419
coefficient equality rows, and `estimated_lowering_would_error=false`.  That
identifies the dual route as the next admissible no-solve model-sizing rung,
but it is not a solve or bound-quality result.  The actual no-solve SOS-dual
model-sizing rung then passed under the same load discipline with
`NCTS_SOLVER_PROBE_LOWER_MODEL=true` and `NCTS_SOLVER_PROBE_SOLVE=false`:
dualization took 1.955824559s, allocated 331.92 MiB, and built a JuMP model
with 106,615 variables, 22 PSD cones, zero scalar zero-dual variables, and
3,419 scalar equality rows.  The full process took 85.24s wall with 2.31 GiB
peak RSS.  This confirms the static estimate for the base source route and
promotes only no-solve A1/A2 source-base model sizing, not solver execution.
The A1 source-base no-solve SOS-dual model-sizing rung then passed with
`contiguous_rdm_k=8` and `linear_state_opt_width=7`, still with no PSO and no
solver.  It reported 3,534 moments, 396 zero rows, 25 PSD cones, RDM blocks
`[72, 64, 56]`, 113,315 estimated SOS-dual variables, 112,919 PSD-cone scalar
variables, 396 zero-dual variables, 3,533 coefficient equality rows, and
`estimated_lowering_would_error=false`.  Actual dualization took
1.942718103s, allocated 387.88 MiB, and built a JuMP model with 113,315
variables, 25 PSD cones, 396 zero-dual variables, and 3,533 equality rows.
The full process took 87.19s wall with 2.30 GiB peak RSS.  The preflight load
was already 50.71 on 25 CPUs, so this is bounded model-size evidence under a
hard timeout, not permission to launch another large run under the same load.
The A2 source-base static estimate then passed with
`contiguous_rdm_k=8`, `linear_state_opt_width=7`, and
`psd_state_opt_width=3`, still with no lowering and no solver.  It reported
5,034 moments, 672 zero rows, 47 PSD cones, RDM blocks `[72, 64, 56]`, the
source-family PSO block vector, 153,751 estimated SOS-dual variables,
153,079 PSD-cone scalar variables, 672 zero-dual variables, 5,033 coefficient
equality rows, and `estimated_lowering_would_error=false`; the process took
92.39s wall with 2.32 GiB peak RSS.  The actual no-solve A2 SOS-dual
model-sizing rung then matched that estimate: dualization took 3.304339117s,
allocated 861.06 MiB, and built a JuMP model with 153,751 variables, 47 PSD
cones, 672 zero-dual variables, and 5,033 equality rows.  The full process
took 101.13s wall with 2.33 GiB peak RSS.  This is model-size evidence only;
it does not authorize an A2_L20 or A2_L30 solve.
Future source-base solve probes must set `NCTS_SOLVER_PROBE_QMBCERTIFY_PROBE_ID`
and `NCTS_SOLVER_PROBE_QMBCERTIFY_OBJECTIVE_TOTAL_ESTIMATE`, so the emitted TOML
row carries the reviewed reference row and
`objective_minus_qmbcertify_total_estimate`.  A RED source-only guard failed on
the missing harness fields in 1.85s with 418,548 KiB peak RSS, and a GREEN
source-only guard passed after adding the optional reference metadata path in
0.22s with 275,388 KiB peak RSS.  A separate plan-phrase RED check failed in
1.81s with 415,400 KiB peak RSS before this requirement was documented.  This
metadata is now also enforced before emitting solved source-base fixture rows:
construction-only, no-solve, and non-source-base rows may omit it, but solved
`qmbcertify_base_construct=true` rows reject without the pair.  It is a
precondition for restoring an A1 source-base solve; without it the run would
answer only "did the solver finish?", not the bound-quality question.
The narrow loader later verified that enforcement: 2,713/2,713 checks passed in
1:15.22 wall with 2,870,928 KiB peak RSS and no swaps after a preflight gate of
load `13.71, 14.21, 14.26` and 675 GiB available memory.
The persisted-fixture side now enforces the same contract: a solved
`qmbcertify_base_construct=true` row must carry a nonempty
`qmbcertify_probe_id`, a `qmbcertify_objective_total_estimate`, and an
`objective_minus_qmbcertify_total_estimate` consistent with the row's objective.
A RED loader contract test failed in 10.80s with 886,504 KiB peak RSS because
the missing-reference synthetic row was accepted.  After adding the fixture
loader check, the same narrow loader passed in 1:10.70 wall with 2,313,064 KiB
peak RSS and no swaps after a preflight gate of load `13.21, 15.35, 14.89` and
669 GiB available memory.
Fresh telemetry on 2026-07-08 01:03 CST reported load averages
`8.76, 11.89, 13.04` on 25 CPUs with 703 GiB available memory, so the narrow
expectation loader was safe to rerun.  It passed in 1:08.61 wall with
2,497,028 KiB peak RSS, including 93/93 Pauli perf harness checks, 30/30
solver-probe fixture-field checks, 14/14 plan-run-gate checks, and 1,855
fixture checks.
The guarded A1 L=20 source-base SOS-dual solve then ran with one Julia thread,
one Mosek thread, `NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=900`, and
`NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=8`.  The launch guard saw load averages
`12.41, 12.89, 13.24` on 25 CPUs with 704 GiB available memory and allowed the
run.  Actual remote wall was 2:22.82 with 2,918,064 KiB peak RSS, well below the
estimate.  Construction took 6.746505033s, SOS-dualization took 1.745493998s,
and Mosek solved in 35.947141334s with `OPTIMAL` status.  The objective was
`-8.905006401419794`; compared with reviewed QMBCertify `A1_L20` total
`-8.904954518917137`, the delta is `-5.188250265675265e-5`.  This nearly closes
the old finite-axis A1 L=20 gap (`-0.187752491273901`) and shows the source-base
formulation fixes the dominant bound-quality mismatch at L=20, but it still
misses the 1e-6 parity target and does not prove L=30 parity.
The new source-base solve row is fixture-pinned in
`heisenberg_qmbcertify_rdm.toml`.  A RED stdlib TOML guard failed on the missing
13th solver-probe row in 3.01s with 428,688 KiB peak RSS; after adding the row
and environment, the same TOML guard passed in 1.48s with 339,092 KiB peak RSS.
A fresh telemetry check on 2026-07-08 01:12 CST reported load averages
`9.55, 10.57, 12.08` on 25 CPUs with 703 GiB available memory, so the narrow
expectation loader was safe to rerun.  It passed in 1:08.00 wall with
2,464,168 KiB peak RSS, including 93/93 Pauli perf harness checks, 30/30
solver-probe fixture-field checks, 14/14 plan-run-gate checks, 101/101 pressure
guard checks, and 1,924 fixture checks.
The L=30 A1 source-base no-solve SOS-dual model-sizing rung also cleared the
shared-machine gate.  It launched with load averages `10.81, 11.21, 12.13` on
25 CPUs, `NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=900`, and
`NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=8`; actual remote wall was 1:37.16 with
2,562,612 KiB peak RSS.  Construction took 16.723815959s, SOS-dualization took
2.731044266s, and no optimizer was called.  The model-size gate was `ok` with
169,600 dual variables, 35 PSD cones, 626 zero-dual variables, 7,404 scalar
equality rows, max block 114, and dense-Schur proxy 438,553,728 bytes.  This
cleared solve-risk/model-size evidence for an L=30 source-base A1 solve, but it
was not bound-quality evidence by itself.  A RED stdlib TOML guard failed
on the missing 14th solver-probe row in 3.14s with 423,856 KiB peak RSS; the
GREEN TOML guard passed after fixture promotion in 1.51s with 336,212 KiB peak
RSS.  With fresh telemetry at load averages `12.69, 11.57, 12.04` on 25 CPUs
and 689 GiB available memory, the narrow expectation loader was safe to rerun.
It passed in 1:09.26 wall with 2,800,004 KiB peak RSS, including 93/93 Pauli
perf harness checks, 30/30 solver-probe fixture-field checks, 14/14
plan-run-gate checks, 101/101 pressure guard checks, and 1,982 fixture checks.
The guarded A1 L=30 source-base SOS-dual solve then ran with one Julia thread,
one Mosek thread, `NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=900`, and
`NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=8`.  The launch guard saw load averages
`14.01, 12.39, 12.25` on 25 CPUs with 681 GiB available memory and allowed the
run.  Actual remote wall was 4:17.87 with 4,161,712 KiB peak RSS, still below
the estimate.  Construction took 18.10667775s, SOS-dualization took
2.759845413s, and Mosek solved in 148.466465284s with `OPTIMAL` status.  The
objective was `-13.326719349737361`; compared with reviewed QMBCertify
`A1_L30` total `-13.325070793538082`, the delta is
`-0.001648556199279838`.  This reduces the old finite-axis A1 L=30 gap
(`-0.339196476185879`) by about two orders of magnitude, but it still misses
the 1e-6 parity target.  Do not repeat this solve unless a formulation change
explains the residual gap first.  The new fixture row passed a TOML-only guard
in 1.49s with 339,076 KiB peak RSS.  After the plan wording was updated, a
second narrow source/TOML guard passed in 1.49s with 336,680 KiB peak RSS.  A
fresh telemetry check on 2026-07-08 01:29 CST reported load averages
`13.07, 12.92, 12.60` on 25 CPUs with 674 GiB available memory, so the narrow
expectation loader was safe to rerun.  It passed in 1:08.38 wall with
2,480,664 KiB peak RSS, including 93/93 Pauli perf harness checks, 30/30
solver-probe fixture-field checks, 14/14 plan-run-gate checks, 101/101 pressure
guard checks, and 2,051 fixture checks.
The formulation audit now follows the same distinction: it prints explicit
source-base metrics from `_pauli_qmbcertify_chain_base_linear_relaxation` instead
of relying only on the older finite-axis comparison.  A RED source guard caught
the missing `nctssos_source_base_*` output path in 1.85s with 418,296 KiB peak
RSS.  The first L=10 no-solver audit smoke then exposed a real audit bug: raw
QMBCertify words were converted directly into `NormalMonomial`, which failed on
repeated site operators after 1:08.00 wall with 2,224,348 KiB peak RSS.  The
next smoke exposed the phase variant of the same bug after 1:09.42 wall with
2,215,480 KiB peak RSS; the audit representative helper now simplifies raw
QMBCertify words and drops the scalar phase for set accounting only.  The L=10
audit smoke then completed in 2:04.87 wall with 2,573,552 KiB peak RSS, no
solver call, and `nctssos_source_base_ok=false` because this short case is not a
valid source-base parity row.  It still proved that the audit exits cleanly and
keeps small invalid source-base cases as diagnostics, not evidence.  The audit
also now prints a separate qmbcertify-mode LSO row comparison instead of
comparing QMBCertify rows only against the contiguous-test generator.  A source
guard for that output passed in 0.20s with 257,624 KiB peak RSS; a no-main audit
include check passed in 5.73s with 702,772 KiB peak RSS; and a targeted runtime
check for qmbcertify-mode LSO generation plus raw-word simplification passed in
6.88s with 732,116 KiB peak RSS.  With fresh telemetry at load averages
`11.92, 12.89, 13.11` on 25 CPUs and 664 GiB available memory, the narrow
expectation loader was safe to rerun; it passed in 1:08.05 wall with
2,786,896 KiB peak RSS, including 37/37 perf harness entry guards and 2,051
fixture checks.
The audit then gained `NCTS_FORMULATION_AUDIT_DIRECT_COMPARE=false`, so the
source-base check can skip the stale finite-axis construction.  A RED helper
check for sign-canonical qmbcertify LSO row keys failed in 6.73s with
747,480 KiB peak RSS; the GREEN helper check passed after adding
`_qmbcertify_linear_state_opt_row_key` in 15.85s with 771,936 KiB peak RSS.  A
guarded L=20 source-base-only audit launched at load averages
`10.76, 12.20, 12.80` on 25 CPUs and completed in 1:28.55 wall with
2,208,904 KiB peak RSS.  It confirmed that sign-canonical LSO dedup reduces the
source-base A1 L=20 zero rows from 396 to 198, exactly matching QMBCertify's
198 filtered LSO rows, while leaving moments, RDM blocks, base block histogram,
and max block unchanged.  The guarded L=20 source-base Mosek solve then ran
with one Julia thread, one Mosek thread, `NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=900`,
and `NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=8`; it completed in 2:21.68 wall with
2,899,104 KiB peak RSS.  Dual variables dropped from 113,315 to 113,117 and
zero-dual variables from 396 to 198, but the objective stayed
`-8.905006401419794`, so the `A1_L20` delta remains
`-5.188250265675265e-5`.  At that point, the residual A1 L=20 source-base
bound-quality gap was therefore not caused by duplicated ± LSO rows.  The
refreshed fixture row passed a TOML-only guard in 1.48s with 337,000 KiB peak
RSS after a whitespace-sensitive plan phrase assertion was corrected.  With
fresh telemetry at load averages
`9.11, 10.26, 11.67` on 25 CPUs and 705 GiB available memory, the narrow
expectation loader was safe to rerun; it passed in 1:08.89 wall with
2,476,268 KiB peak RSS, including 41/41 perf harness entry guards and 2,051
fixture checks.
The same sign-dedup refresh was replayed at L=30 only after a fresh pressure
gate reported load averages `12.47, 10.74, 11.64` on 25 CPUs and 703 GiB
available memory.  The run used one Julia thread, one Mosek thread,
`NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=900`, and
`NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=8`; it completed in 4:19.74 wall with
4,243,784 KiB peak RSS.  LSO rows dropped from 626 to 313, dual variables
dropped from 169,600 to 169,287, and zero-dual variables dropped from 626 to
313, but the objective stayed `-13.326719349737361`; the `A1_L30` delta
therefore remains `-0.001648556199279838`.  At that point, the remaining A1
source-base bound-quality gaps were not caused by duplicated ± LSO rows at L=20
or L=30.  Do not repeat these source-base A1 solves unless a new formulation
change moves the bound or a cheaper no-solver audit predicts a different model.
The refreshed fixture row passed a TOML-only guard in 1.56s with 337,452 KiB
peak RSS.  With fresh telemetry at load averages `9.67, 11.71, 12.08` on
25 CPUs and 704 GiB available memory, the narrow expectation loader was safe to
rerun; it passed in
1:09.31 wall with 2,498,764 KiB peak RSS, including 41/41 perf harness entry
guards and 2,051 fixture checks.
A source-only normalization guard now pins the source-base identity cross-term
scale `1/sqrt(L)` and the fixed-momentum QMBCertify real-block path.  The guard
passed in 0.22s with 258,488 KiB peak RSS and did not import NCTSSoS, build a
model, or call a solver.  This rules out accidental deletion of the known
QMBCertify block-normalization convention, but it does not explain the residual
bound gap.  After adding that guard, the narrow expectation loader passed in
1:07.50 wall with 2,323,212 KiB peak RSS, including 43/43 perf harness entry
guards and 2,051 fixture checks.
The formulation audit now also prints an objective-coefficient comparison:
NCTSSoS's reduced total Hamiltonian coefficients are compared with
QMBCertify's per-site objective coefficients scaled by `L`.  A targeted L=20
A1 helper check imported NCTSSoS but did not load QMBCertify, construct
source-base blocks, build a model, or call a solver; it passed in 8.72s with
780,996 KiB peak RSS and `max_abs_delta=0.0`.  The residual A1 gap is therefore
not an objective scaling/sign mismatch.  After adding the audit output guards,
the narrow expectation loader passed in 1:09.02 wall with 2,867,872 KiB peak
RSS, including 46/46 perf harness entry guards and 2,051 fixture checks.
The formulation audit now also has a source-entry check for QMBCertify RDM
blocks.  For L=20 and k=8, the parsed `posepsd8!` row blocks match
`pauli_qmbcertify_rdm_blocks(8; ambient_sites=20)` exactly, and an independent
bitwise Pauli-entry formula matches NCTSSoS's source-base RDM linear entries
after the harmless global `2^-8` scale.  The persistent helper check compared
12,416 RDM entries with max delta `0.0` in 12.02s wall and 946,140 KiB peak
RSS, with no model construction or solver call.  The residual A1 gap is
therefore not an RDM row-orientation or coefficient-scaling mismatch.
After adding the persistent audit guard, the narrow expectation loader passed
in 1:07.99 wall with 2,385,216 KiB peak RSS, including 93/93 Pauli perf harness
guards, 53/53 perf harness entry guards, 101/101 pressure guards, 14/14
plan-run gates, and 2,051 fixture checks.
An integrated L=20 A1 formulation audit with
`NCTS_FORMULATION_AUDIT_DIRECT_COMPARE=false` then ran under an explicit 180s /
3 GiB guard after telemetry reported load averages `17.83, 12.91, 12.81` on 25
CPUs.  It completed in 1:30.63 wall with 2,206,144 KiB peak RSS, compared
12,416 source RDM entries with max delta `0.0`, and preserved 3,534 source
moments, 198 source-base LSO rows, RDM blocks `[72, 64, 56]`, and source-base
block histogram `[48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]`.  This is
still no-solver formulation evidence; it does not justify repeating L=20/L=30
A1 solves without a new formulation change.
The source-base audit then gained a base PSD entry comparison against
QMBCertify's chain `reduce!(...; realify=true)` entries.  A broad
`test/expectations_loader.jl` include verifier was deleted from the live queue
after it hit a 45s timeout and 2,515,028 KiB peak RSS because it executes the
full loader rather than a narrow source guard.  The valid L=12 base-entry smoke
first compared 61,771 entries and exposed 37,666 mismatches in about 54s /
1.12 GiB peak RSS, proving a small counterexample was enough and deleting any
L=20 entry-comparison run from the active ladder.  The production fix carries
Pauli simplification phases through `_qmbcertify_chain_support_term` and uses
QMBCertify-style realification for source-base products; the audit fix drops
near-zero source coefficients both before and after source-form aggregation.
Intermediate L=12 smokes reduced mismatches to 29,404 and then 8,076; two
block-local diagnostics showed the remaining samples were cancellation
residuals (`1.22e-16` and `-6.66e-16`) that NCTSSoS already cleaned at
`1e-12`.  After telemetry reported load averages `9.53, 11.87, 13.24` on 25
CPUs with ample memory, the final L=12 source-base entry smoke completed in
54.19s wall with 1,130,588 KiB peak RSS and reported
`compared=61771`, `match=true`, max delta `1.2434497875801753e-14`, and
zero mismatches.  This is the decisive small verification for the phase-loss
bug; it does not reopen L=20/L=30 no-solver or solver runs.
After that small counterexample was green, the same phase-carry fix justified
one guarded L=20 no-solver audit, still without lowering or solving.  Telemetry
reported load averages `11.26, 12.17, 13.05` on 25 CPUs with 739 GiB available
memory, so the run was launched with one Julia thread, a 240s hard timeout, and
an explicit 180s / 3 GiB estimate.  The large-run pressure gate accepted the
run at load `11.61` on 25 CPUs.  It completed in 1:33.30 wall with
2,264,004 KiB peak RSS, matching the estimate.  The audit preserved the A1 L=20
source-base structure (3,534 linear moments, 198 LSO rows, RDM blocks
`[72, 64, 56]`, base block histogram `[48 => 2, 57 => 1, 58 => 1, 96 => 9, 114 => 9]`)
and now compared all 106,615 source-base PSD entries with
`nctssos_source_base_entry_match=true`, zero mismatches, and max delta
`2.1760371282653068e-14`.  The RDM source-entry check still compared 12,416
entries with max delta `0.0`.  This is the first L=20 no-solver evidence that
the source-base PSD entries themselves match QMBCertify after the phase fix.
It justifies considering a separately gated L=20 A1 solve replay, but it does
not by itself authorize any L=30 or A2 run.
That separately gated A1 L=20 source-base solve replay was then launched only
after telemetry reported load averages `12.82, 14.10, 13.74` on 25 CPUs with
about 738 GiB available memory.  The run used one Julia thread, one Mosek
thread, a 600s timeout, and explicit 300s / 4 GiB estimates; the pressure gate
accepted it at load `10.64, 13.21, 13.45` on 25 CPUs.  It completed in 2:22.57
remote wall with 2,883,240 KiB peak RSS, below estimate.  Construction took
6.750735315s, SOS-dualization took 1.767317094s, and Mosek solved in
35.763360136s with `OPTIMAL` status.  The objective was
`-8.904954795042105`, so the delta against reviewed QMBCertify `A1_L20`
(`-8.904954518917137`) is `-2.7612496822371213e-7`, inside the 1e-6 parity
target.  This closes the A1 L=20 bound-quality parity gap for the source-base
route after phase realification.  It does not authorize L=30 or A2 replay
without a separate no-solver/formulation reason and fresh load gate.  The first
full fixture-loader rerun after promoting this row completed in 1:11.74 wall
with 2,547,248 KiB peak RSS but failed one stale policy assertion that assumed
all source-base solver rows were either unreviewed or not parity.  After
updating that guard to allow the reviewed A1 L=20 parity row only, the same
loader passed in 1:10.47 wall with 2,849,284 KiB peak RSS, including 36/36
solver-probe fixture-field checks and 2,051 fixture checks.
After the L=20 solve replay met parity, the same phase-realification change
justified one L=30 no-solver source-base entry audit before considering any
L=30 solve replay.  Fresh telemetry reported load averages
`15.02, 15.77, 14.43` with about 704 GiB available memory; the audit used one
Julia thread, no optimizer, a 300s hard timeout, and explicit 240s / 4 GiB
estimates.  The pressure gate accepted it at load `20.91, 17.12, 14.90` on
25 CPUs.  It completed in 1:54.25 wall with 2,421,048 KiB peak RSS, below
estimate.  The audit preserved the A1 L=30 source-base structure (7,405 linear
moments, 313 LSO rows, RDM blocks `[72, 64, 56]`, base block histogram
`[48 => 2, 57 => 1, 58 => 1, 96 => 14, 114 => 14]`) and compared all 162,670
source-base PSD entries with `nctssos_source_base_entry_match=true`, zero
mismatches, and max delta `6.772360450213455e-14`.  The RDM source-entry check
still compared 12,416 entries with max delta `0.0`.  This is sufficient
no-solver formulation evidence to consider one separately gated A1 L=30 solve
replay.  It does not authorize A2 or N=100 numerical runs.
That separately gated A1 L=30 source-base solve replay was then launched only
after fresh telemetry reported load averages `19.50, 18.80, 15.94` with about
704 GiB available memory.  The run used one Julia thread, one Mosek thread, a
900s timeout, and explicit 360s / 5 GiB estimates; the pressure gate accepted
it at load `17.10, 18.28, 15.84` on 25 CPUs.  It completed in 4:09.77 remote
wall with 4,241,080 KiB peak RSS, below estimate.  Construction took
17.461991117s, SOS-dualization took 2.684905077s, and Mosek solved in
139.546716090s with `OPTIMAL` status.  The objective was
`-13.325070927620173`, so the delta against reviewed QMBCertify `A1_L30`
(`-13.325070793538082`) is `-1.3408209120768788e-7`, inside the 1e-6 parity
target.  This closes the A1 L=30 bound-quality parity gap for the source-base
route after phase realification.  It does not authorize A2 or N=100 numerical
runs.  The first full fixture-loader rerun after promoting this row completed
in 1:10.32 wall with 2,541,672 KiB peak RSS but failed on missing provenance
for the renamed L=30 environment key.  After adding the matching environment
entry, the same loader passed in 1:08.45 wall with 2,455,480 KiB peak RSS,
including 36/36 solver-probe fixture-field checks and 2,051 fixture checks.
The solver probe can now set QMBCertify-matching Mosek tolerances explicitly
through `NCTS_SOLVER_PROBE_MOSEK_TOL_PFEAS`,
`NCTS_SOLVER_PROBE_MOSEK_TOL_DFEAS`, and
`NCTS_SOLVER_PROBE_MOSEK_TOL_REL_GAP`, and emits those values in future
environment fixtures when present.  A source-only guard for the new knobs
passed in 0.21s with 275,332 KiB peak RSS.  A first include-level check without
developing the local package into the docs environment failed before reaching
the new code; the retry using the same docs
`Pkg.develop(PackageSpec(path=pwd()))` setup as solver probes constructed the
Mosek optimizer wrapper in 10.63s with 1,026,136 KiB peak RSS.  That
tolerance-knob-only replay remains optional evidence, not queued work; the
separate phase-fix A1 L=20 and L=30 solve replays above are the live parity
evidence.  The
narrow expectation loader passed after this change in 1:08.91 wall with
2,499,308 KiB peak RSS, including
36/36 solver-probe fixture-field checks and 2,051 fixture checks.
The qmbcertify-mode LSO append order now matches QMBCertify's `filter_mons`
more closely: support coverage is checked before a row enters the sign-canonical
dedup set.  This prevents an uncovered duplicate from suppressing a later
covered row.  A source-only guard for the ordering passed in 0.20s with
258,352 KiB peak RSS.  A guarded L=20 source-base no-solver shape check then
ran with one Julia thread after telemetry reported load averages
`12.89, 13.78, 13.02` on 25 CPUs and 705 GiB available memory; it passed in
36.72s with 1,138,596 KiB peak RSS, preserving 3,534 moments, 3,505 free keys,
198 zero constraints, 198 LSO rows, RDM blocks `[72, 64, 56]`, and max block
114.
After this patch, the narrow expectation loader passed in 1:09.83 wall with
2,864,404 KiB peak RSS, including 48/48 perf harness entry guards, 36/36
solver-probe fixture-field checks, and 2,051 fixture checks.
A subsequent full expectation-loader verifier after the lifted-run-policy plan
update passed in 1:08.13 wall with 2,448,140 KiB peak RSS, including 93/93
Pauli perf harness guards, 48/48 perf entry guards, 101/101 pressure guards,
14/14 plan-run gates, and 2,051 fixture checks.
The source-base rows now also pin an 8-byte dense-Schur storage proxy from the
equality-row upper bound: 93,516,488 bytes for A0, 99,856,712 bytes for A1, and
202,648,712 bytes for A2.  This is only a memory-risk proxy; it is not a solver
allocation measurement and does not replace the load/RSS gate.
The solver probe now turns that proxy into a post-construction
`estimated_model_size_gate_status`.  For intentional large runs, lowering or
solving is refused when this status is not `"ok"`, so a large construction-only
probe can still delete the next model-build rung before it crowds the remote.
Structural target rows now also carry `estimated_model_size_gate_status`; because
they do not construct a scalar-equality upper bound or dense-Schur proxy, the
status is `"blocked_missing_scalar_equality_estimate"`.  Treat this as a hard
defer/delete signal for lowering or solve attempts; block-storage targets alone
are not solver-memory evidence.
The RED stdlib-only TOML/source guard failed on the missing structural target
field in 3.00s with 429,884 KiB RSS.  After adding the status to the structural
target APIs and fixtures, the same guard passed in 1.52s with 338,036 KiB RSS;
the expanded guard including plan text passed in 1.50s with 339,788 KiB RSS.
A target-only package smoke at N=10, d=2 then checked the ordinary translation,
SU(2) translation-orbit, and SU(2) full-basis target APIs in 32.97s with
1,230,520 KiB RSS.  No model was constructed and no solver was called.
The target-only profile printer now emits the same model-size gate status and
reason beside block-storage metrics, so profile output cannot silently promote
structural storage to solve evidence.  A RED source guard failed on the missing
profile string in 1.85s with 417,156 KiB RSS; the GREEN source guard passed in
0.22s with 257,216 KiB RSS.  A target-only N=10, d=2 profile smoke then grepped
the rendered output for the status in 43.50s with 2,246,808 KiB RSS, with no
model construction and no solver call.
Constructed no-solver reports now expose the analogous SOS-dual
`estimated_sos_dual_scalar_equalities_upper_bound` and
`estimated_sos_dual_dense_schur_bytes` fields from the report's linear moment
count.  The formula mirrors the solver-probe estimator: real PSD reports use
`linear_moment_count - 1`, while complex/HPSD reports use
`2 * linear_moment_count - 1`; zero-dual rows remain variables, not equality
rows.  A RED source guard failed on the missing fields/profile labels in 1.82s
with 417,008 KiB RSS; the GREEN source guard passed in 0.20s with 257,964 KiB
RSS.  A no-solver N=4/order-1 package smoke checked both real and complex report
metrics in 47.25s with 1,449,868 KiB RSS.  A constructed-profile output smoke
then grepped the rendered summary for the new SOS-dual equality and dense-Schur
rows in 78.86s with 2,294,592 KiB RSS.  Neither smoke built a JuMP model or
called a solver.  The constructed profile now also prints a `SOS-dual model-size
gate status` computed from that dense-Schur proxy and current remote memory
headroom.  A RED source guard failed on the missing status label in 1.84s with
419,536 KiB RSS; the GREEN source guard passed in 0.20s with 257,224 KiB RSS.
The no-solver N=4/order-1 profile smoke then grepped the rendered status row in
78.38s with 2,132,880 KiB RSS.  The constructed profile now prints the same
model-size gate reason and telemetry inputs beside the status:
`SOS-dual model-size gate reason`, estimated RSS bytes, `MemAvailable` bytes,
and max RSS fraction.  A RED source guard failed on the missing labels in 1.82s
with 417,896 KiB RSS; the GREEN source guard passed in 0.23s with 277,408 KiB
RSS.  A helper smoke checked the blocked-headroom reason without constructing a
relaxation in 4.22s with 629,356 KiB RSS.  After a fresh safe-run telemetry
check reported load averages `10.92, 13.56, 15.75` on 25 CPUs with 695 GiB
available memory, a no-solver N=4/order-1 rendered-profile smoke grepped all
four new rows in 76.75s with 2,121,360 KiB RSS.
The solver-probe risk gate is now fixture-visible as a plain status string.
A RED source guard failed because `estimated_risk_gate_status` and the
`blocked_lowering_orphan_policy` label were absent; the GREEN source guard
passed after printing the status in probe output and emitted TOML rows.  This
was a source-only harness check under a 30s timeout, one Julia thread, and no
solver-probe include, because the remote load gate was `35.69` on 25 CPUs.
The existing source-base A0/A1/A2 SOS-dual fixture rows now carry the promoted
status too: `estimated_risk_gate_status = "ok"` beside
`estimated_lowering_would_error = false`.  A RED fixture-level TOML parse failed
on the missing key; the GREEN narrow parse passed after updating the three rows.
The check used only stdlib TOML under a 30s timeout, one Julia thread, and no
package load, because a fresh remote gate reported `44.79` load on 25 CPUs.
The base/A1/A2 source-base no-solve model-size rows are now fixture-pinned in
`test/data/expectations/heisenberg_qmbcertify_rdm.toml` under
`nctssos_qmbcertify_source_base_model_probes`, separate from the older
finite-axis/source-like rows.  The narrow loader check passed through easy-ssh
with 1,812 fixture checks in 71.49s wall and 2.36 GiB peak RSS, using one Julia
thread, no optimizer call, and a 180s timeout.  The preflight load was 65.02 on
25 CPUs, so this verification does not reopen the larger run queue.  A fresh
pressure check immediately afterward still reported load average 49.80 on
25 CPUs, while the existing evidence puts `test/relaxations/runtests.jl` near
9m12s and package-level `Pkg.test()` near 23m4s.  Those broad verification
runs are therefore deleted from the current live queue until a fresh load gate
shows enough headroom; running them now would crowd autodl.
The QMBCertify reviewed-reference status policy is now explicit in
`perf/qmbcertify_reference_runs.jl`: only `OPTIMAL` and `ALMOST_OPTIMAL`
are accepted for reviewed rows.  `SLOW_PROGRESS`/`FEASIBLE_POINT` remains a
failed attempt, not a reviewed reference.  A TDD RED check failed on the
missing policy helper in 3.12s with 486 MiB peak RSS; the GREEN check passed in
1.17s with 386 MiB peak RSS.  The loader integration rerun then passed through
easy-ssh with 1,812 fixture checks and 12/12 QMBCertify provenance checks in
70.45s wall with 2.73 GiB peak RSS.  The preflight load was 61.14 on 25 CPUs,
so no broader verification or solver run follows from this.
The same accepted-status policy is now fixture-pinned in both
`heisenberg_qmbcertify_base.toml` and `heisenberg_qmbcertify_rdm.toml` under
`reviewed_reference_status_policy`, so row review can be audited from data
without reading the harness code.  A TDD RED fixture parse failed on the missing
table in 0.05s with 12 MiB peak RSS; the GREEN parse passed in 0.06s with
13 MiB peak RSS.  The loader integration rerun then passed through easy-ssh
with 1,818 fixture checks and 12/12 QMBCertify provenance checks in 69.94s wall
with 2.73 GiB peak RSS.  The preflight load was 36.18 on 25 CPUs, so broader
verification remains behind the fresh-load gate.
The reference harness now consumes that fixture-backed policy by profile before
accepting reviewed rows.  A TDD RED check failed on the missing
`_reviewed_solve_status_policy(::String)` method in 2.76s with 491 MiB peak
RSS; the GREEN check passed in 1.26s with 387 MiB peak RSS and verified that
`SLOW_PROGRESS`/`FEASIBLE_POINT` is still rejected under the A2 fixture policy.
The loader integration rerun then passed with 1,818 fixture checks and 14/14
QMBCertify provenance checks in 71.14s wall with 2.74 GiB peak RSS.  The
preflight load was 62.43 on 25 CPUs, so this still authorizes no broader run.
`perf/qmbcertify_profile_probe.jl` now keeps probe labels such as `A1_probe`
for emitted rows while passing the parent fixture profile into the reviewed
status-policy lookup.  A RED source guard failed in 2.17s wall because the
probe did not pass `status_policy_profile`; the GREEN guard passed in 2.16s
wall after adding `status_policy_profile=profile` and the shared `_run_one`
keyword.  A no-solver Julia include check then passed in 6.78s wall under a
45s timeout with one Julia thread.  The latest remote pressure gate reported
load average 70.97 on 25 CPUs and 647 GiB available memory, so this evidence
only justifies short no-solver checks; it deletes broad Julia suites and
solver-backed probes from the current live queue until the CPU gate improves.
The same profile-probe path now prints and emits the parent profile's accepted
termination-status policy.  A RED harness check failed on the missing
`_probe_accepted_statuses("A2")` helper; the GREEN helper check passed after
threading the policy through the probe header and both `profile_probes` and
`failed_attempts` rows.  No-solver rendered-row checks then confirmed both
emitted probe TOML tables contain `accepted_termination_statuses = ["OPTIMAL",
"ALMOST_OPTIMAL"]`.  These checks used one Julia thread, 45s hard timeouts, no
QMBCertify load, no model construction, and no solver call.
The profile-probe bootstrap path now also returns before loading or bootstrapping
QMBCertify, matching the reviewed-reference harness behavior.  A RED source
order guard failed because `_load_qmbcertify` still preceded the
`bootstrap_only` branch; the GREEN source guard passed after moving the load
below the bootstrap return.  A bootstrap behavior check then ran with an
intentionally missing `NCTS_QMBCERTIFY_PATH`, `NCTS_QMB_BOOTSTRAP_ONLY=true`,
one Julia thread, and a 45s timeout; it emitted the environment fixture and
exited 0 without creating the missing path.  This keeps exploratory preflight
checks cheap under high remote load.
The A2_L30 solve stays removed from the active run queue until a deliberate
QMBCertify A2 profile variant produces an accepted reviewed status under this
policy and the matching no-solver/model-size NCTSSoS source-base row clears the
explicit solve-risk gate.  The no-solver base histogram/support-quotient
mismatch is already closed, so it is not the blocker for re-adding that solve.
A target-only
N=100 structural profile
with SU(2) RDM, width-7 LSO, and width-3 PSD state-opt completed in about 54s
without constructing a model or calling a solver; it keeps the shape/accounting
target in scope but is not numerical parity evidence.  The remaining unproven
items are open scientific questions, not queued runs: reviewed QMBCertify
A0_L20/A0_L30 and A2 reference rows, N=100 constructed/numerical parity claims,
and the L=30 SU(2) runtime gate after solver-formulation changes.  The first
native L=30 SU(2) solve attempt timed out at 30 minutes, so the solved row is
deleted from the live plan until formulation evidence changes the estimate.
A2 still needs reviewed QMBCertify profile-status evidence before any solve.  A
structural target whose `estimated_model_size_gate_status` is
`"blocked_missing_scalar_equality_estimate"` is not enough solve-risk evidence.
After A1 L=20/L=30 source-base parity closed, the first L=30 SU(2) gate was a
target-only profile, not lowering or solving.  Fresh telemetry reported load
averages `10.28, 12.65, 13.89` with about 704 GiB available memory, so the run
used one Julia thread, no model construction, no solver, a 180s timeout, and
explicit 90s / 2 GiB estimates.  The pressure gate accepted it at load
`10.01, 12.52, 13.83` on 25 CPUs.  It completed in 46.59s wall with
2,211,856 KiB peak RSS.  The structural target reported 85 solver-facing PSD
blocks, overall solver-facing max block 62, 819 known add-on zero rows, SU(2)
RDM blocks `[56, 40, 28, 14, 1]`, PSD state-opt solver-facing blocks
`[18 => 15, 9 => 1]`, and SU(2) translation-orbit solver-facing max block 22.
It still reported `estimated_model_size_gate_status =
"blocked_missing_scalar_equality_estimate"`, so it does not authorize lowering
or solving by itself; the next valid rung is a separately gated no-solve
solver-probe/model-sizing run.
The separately gated L=30 SU(2) no-solve solver-probe/model-sizing pair then
ran with one Julia thread and one Mosek-thread model backend, still with
`SOLVE=false`.  The native-Hermitian row launched after telemetry reported load
averages `9.99, 11.82, 13.47` and about 704 GiB available memory; the run used
a 600s timeout with explicit 300s / 12 GiB estimates and was accepted at load
`10.66, 11.88, 13.46` on 25 CPUs.  It completed in 3:56.47 wall with
6,179,184 KiB peak RSS.  Construction took 100.034857014s, native dualization
took 49.760605312s, and the no-solve dual model had 155 native Hermitian PSD
cones, 112,740 dual variables, 122,385 scalar equalities, max block 31, 819
zero-dual rows, and dense-Schur proxy 119,824,705,800 bytes with
`estimated_model_size_gate_status = "ok"`.
The matching real-lift no-solve row launched after telemetry reported load
averages `17.16, 15.73, 14.76` and about 705 GiB available memory; the run used
a 900s timeout with explicit 480s / 24 GiB estimates and was accepted at load
`15.22, 15.37, 14.66` on 25 CPUs.  It completed in 4:10.78 wall with
8,121,032 KiB peak RSS.  Construction took 97.611421667s, real-lift
dualization took 65.356786939s, and the no-solve dual model had the same 155
PSD cones, 228,602 dual variables, 122,385 scalar equalities, max block 31, 819
zero-dual rows, and the same dense-Schur proxy gate status.  Native Hermitian
therefore cuts the L=30 SU(2) no-solve dual variable count by about 50.7% and
reduces peak RSS by about 1.85 GiB.  This is model-size/runtime evidence for
choosing native before any L=30 SU(2) solve, not numerical bound-quality
evidence.
One separately gated native L=30 SU(2) solve candidate was then attempted only
after telemetry reported load averages `9.89, 14.14, 14.50` with about 703 GiB
available memory.  The run used one Julia thread, one Mosek thread, a 1800s
timeout, and explicit 1200s / 160 GiB estimates; the pressure gate accepted it
at load `9.39, 13.73, 14.36` on 25 CPUs.  Construction and native dualization
matched the no-solve row (99.946889994s construction, 47.821555527s
dualization, 112,740 dual variables, 122,385 scalar equalities, dense-Schur
proxy 119,824,705,800 bytes), but Mosek did not return a solver status before
the timeout.  The job exited 124 after 30:03.89 wall with 81,381,828 KiB peak
RSS and no swaps.  This did not stall the remote, but it fails the L=30 SU(2)
runtime gate.  Delete L=30 SU(2) solved rows from the active queue until the
solver formulation changes, e.g. by reducing the 122,385 scalar equalities or
otherwise avoiding the dense-Schur bottleneck.  Keep only the target-only and
no-solve model-size rows as current SU(2) evidence.
Do not launch any of those larger runs unless the existing
timing/RSS evidence predicts that the run will not crowd the machine; if
evidence says a run would stall autodl, remove or downgrade that run from this
plan instead of trying it.
The previous blanket N>=14 performance-test ban is lifted, but only as a
load-gated policy: each run must carry a wall-time/RSS/load estimate before
launch, and any run whose estimate predicts crowding is deleted or downgraded
from the live plan rather than attempted.
The large-run model-sizing template now uses the emitted
`estimated_risk_gate_status`, not just `estimated_lowering_would_error`, to
decide whether lowering is allowed.  A RED source guard failed on the missing
phrase; the first GREEN attempt exposed a line-wrap issue; the final source-only
guard passed through easy-ssh in 0.35s wall with 298,636 KiB peak RSS, one
Julia thread, and a 30s hard timeout.  The preflight gate was load
`27.47, 33.49, 37.15` on 25 CPUs with 680 GiB available memory, so the check
stayed source-only and no package-load, solver, model-construction, or broad
test run was launched.
The perf harnesses now enforce the same policy before intentional large runs:
`perf/shared_load_guard.jl` provides a shared load-average gate, and
`pauli_translation_profile.jl`, `pauli_translation_solver_probe.jl`, and
`qmbcertify_reference_runs.jl` call it when their `ALLOW_LARGE` flag is set.
A RED source check failed because the shared guard file did not exist.  The
GREEN source/wiring check passed through easy-ssh in 0.40s wall with
302,436 KiB peak RSS.  A forced-overload behavior check then set the guard CPU
count to one, saw `blocked_overloaded_remote` at load `45.03`, and confirmed
that `NCTS_LOAD_GUARD_ALLOW_OVERCOMMITTED=true` bypasses the guard for isolated
hardware only; it passed in 0.65s wall with 319,408 KiB peak RSS.  These checks
loaded only the guard helper and source files, not NCTSSoS, QMBCertify, JuMP,
or a solver.
The guard was then moved earlier in the harness lifecycle so it runs before
heavy imports.  A RED source-order check failed because
`pauli_translation_solver_probe.jl` included `shared_load_guard.jl` after
`using NCTSSoS`.  After moving the pre-import large-run checks above the
NCTSSoS/JuMP/Pkg imports, the GREEN source-order check passed through easy-ssh
in 0.40s wall with 297,540 KiB peak RSS.  A forced-block smoke with
`NCTS_TRANSLATION_ALLOW_LARGE=true` and `NCTS_TRANSLATION_LOAD_GUARD_CPUS=1`
reported `blocked_overloaded_remote` at load `31.76` and stopped from
`_translation_profile_preimport_large_run_pressure_guard`, before the profile
harness imports NCTSSoS.
The pressure gate now checks the 1-, 5-, and 15-minute load averages rather
than only the one-minute value.  A RED guard-only check failed because
`_ncts_large_run_pressure_status` did not accept `load5` or `load15`.  The
GREEN guard-only check passed through easy-ssh in 0.55s wall with 308,596 KiB
peak RSS, covering both 5-minute and 15-minute overload cases.  A follow-up
forced-block smoke with `NCTS_TRANSLATION_ALLOW_LARGE=true` and
`NCTS_TRANSLATION_LOAD_GUARD_CPUS=25` stopped before NCTSSoS imports and
reported `blocked_overloaded_remote` at load averages `42.52, 38.95, 38.74`;
the smoke took 2.47s wall with 420,260 KiB peak RSS.
The guard also now derives its CPU denominator from cgroup v2 `cpu.max` before
falling back to `nproc` or Julia's `Sys.CPU_THREADS`, because autodl reports
`nproc=25` while `/proc/self/status` exposes `Cpus_allowed_list=0-207`.  A RED
guard-only check failed because `_ncts_cpu_max_quota_count` did not exist; the
GREEN guard-only check parsed `2500000 100000` as 25 CPUs and passed in 0.62s
wall with 320,416 KiB peak RSS.  A real pre-import smoke with
`NCTS_TRANSLATION_ALLOW_LARGE=true` and no manual CPU override then inferred
`large_run_pressure_gate_cpus=25.0`, reported load averages
`37.56, 44.24, 41.29`, and blocked before importing NCTSSoS.
Explicit load-guard overrides are now auditable too.  A RED source check
failed because override mode returned before printing any status.  After moving
the override decision below load/cgroup collection, the guard prints
`large_run_pressure_gate_status=override_overcommitted_remote` and
`large_run_pressure_gate_override=true` before returning.  The shared guard now
also prints `large_run_pressure_gate_reason`, so a blocked or overridden run log
contains the human-readable reason as well as the status code.  A RED source
guard failed on the missing reason field in 1.83s with 419,908 KiB RSS; the
GREEN source guard passed in 0.21s with 273,416 KiB RSS.  A first behavior check
used an invalid `IOBuffer` stdout capture and failed in 2.58s with 447,896 KiB
RSS; the corrected temp-file capture verified the printed reason in 0.43s with
296,840 KiB RSS.  Override logs now also print
`large_run_pressure_gate_blocked_status`, preserving the original blocked status
as a machine-readable field even when the emitted status is
`override_overcommitted_remote`.  A RED source guard failed on the missing field
in 1.80s with 419,396 KiB RSS; the GREEN source/parser guard passed in 0.21s
with 273,212 KiB RSS, and a synthetic print check passed in 0.46s with
308,624 KiB RSS.  The earlier override guard-only verification passed through
easy-ssh in 0.79s wall with 340,660 KiB peak RSS; no package or solver imports
were loaded.
The large-run guard now prints independent load, estimate, and memory
sub-statuses: `large_run_pressure_gate_load_status`,
`large_run_pressure_gate_estimate_status`, and
`large_run_pressure_gate_memory_status`.  This prevents an overloaded remote
from hiding a missing wall-time/RSS estimate in the audit log.  A RED
source-only guard failed on the missing output field in 1.84s with
416,120 KiB RSS; the GREEN source-only guard passed in 0.22s with 275,116 KiB
RSS.  A synthetic override capture then forced an overload while omitting both
estimates and verified that the log still reported
`blocked_missing_wall_estimate` as the estimate sub-status; it passed in 1.02s
with 332,932 KiB RSS and no NCTSSoS import.  After a fresh telemetry check
reported load averages `10.25, 11.62, 13.45` on 25 CPUs with 695 GiB available
memory, the full `test/expectations_loader.jl` gate was safe to rerun; it passed
in 1:07.48 wall with 2,471,020 KiB peak RSS, including 89/89 pressure-guard
checks and 1,855 fixture checks.
The estimate sub-status is now split further into wall-time and RSS estimate
statuses: `large_run_pressure_gate_wall_estimate_status` and
`large_run_pressure_gate_rss_estimate_status`.  The aggregate
`large_run_pressure_gate_estimate_status` remains for compatibility, but the log
now shows both missing estimates when both are absent.  A RED source guard failed
on the missing fields in 1.87s with 416,920 KiB RSS; the GREEN source guard
passed in 0.23s with 275,556 KiB RSS.  A synthetic override capture then omitted
both estimates and verified that the same log contained
`blocked_missing_wall_estimate` and `blocked_missing_rss_estimate`; it passed in
1.04s with 332,220 KiB RSS and no NCTSSoS import.  A fresh telemetry check
reported load averages `16.16, 14.45, 14.21` on 25 CPUs with 703 GiB available
memory, so the full expectation loader was safe to rerun; it passed in
1:06.14 wall with 2,323,488 KiB peak RSS, including 91/91 pressure-guard checks
and 1,855 fixture checks.
The guard now prints per-substatus reason fields:
`large_run_pressure_gate_load_reason`,
`large_run_pressure_gate_estimate_reason`,
`large_run_pressure_gate_wall_estimate_reason`,
`large_run_pressure_gate_rss_estimate_reason`, and
`large_run_pressure_gate_memory_reason`.  This makes a blocked log auditable
when load, wall/RSS estimates, and memory headroom have independent outcomes.
A RED source-only check failed on the missing plan phrase in 1.83s with
420,580 KiB peak RSS, and a RED print-signature check failed on the missing
reason keywords in 2.36s with 438,328 KiB peak RSS.  The direct print-path
GREEN check passed in 0.39s with 292,076 KiB peak RSS.  A forced-overload,
missing-estimate guard capture then verified the computed load, estimate,
wall-estimate, RSS-estimate, and memory reason fields in 1.04s with 333,816 KiB
peak RSS.  These checks included only `perf/shared_load_guard.jl`; no NCTSSoS
import, solver call, or model construction was launched.  A fresh telemetry
check on 2026-07-08 00:57 CST reported load averages `9.77, 11.49, 13.03` on
25 CPUs with 704 GiB available memory, so the narrow expectation loader was
safe to rerun.  It passed in 1:07.95 wall with 2,840,892 KiB peak RSS, including
101/101 pressure-guard checks, 13/13 plan-run-gate checks, and 1,855 fixture
checks.
Missing load-average telemetry is now blocking by default.  A RED guard-only
check failed because all-missing load values returned `unknown_no_loadavg`;
the GREEN check now returns `blocked_missing_loadavg`, while still allowing a
non-missing under-threshold load window.  That guard-only verification passed
through easy-ssh in 0.37s wall with 303,296 KiB peak RSS.  Explicit override
remains the only path that can proceed without load telemetry, and it is
audited as above.
The same pre-import guard now requires explicit wall-time and RSS estimates
for any intentional large run.  A RED guard-only check failed because
`_ncts_large_run_estimate_status` was missing; the GREEN guard-only check
passed through easy-ssh in 0.49s wall with 310,072 KiB peak RSS, requiring
`blocked_missing_wall_estimate` and `blocked_missing_rss_estimate` before a
large harness may proceed.  A forced memory-block smoke with
`NCTS_TRANSLATION_ESTIMATED_WALL_SECONDS=60`,
`NCTS_TRANSLATION_ESTIMATED_RSS_GIB=1000000`, and a deliberately permissive
load threshold reported `blocked_insufficient_memory` before importing
NCTSSoS, with load averages `17.39, 20.38, 21.03`, about 646.6 GiB
available memory, and 427,808 KiB peak RSS in 2.47s wall.  A separate
forced-block smoke with no estimates reported `blocked_missing_wall_estimate`
before importing NCTSSoS, with load averages `20.17, 21.11, 21.28`,
about 646.7 GiB available memory, and 427,984 KiB peak RSS in 2.56s wall.
This proves the lifted `N>=14` policy is still gated by explicit run estimates
and current machine telemetry; runs without those estimates stay deleted from
the live queue.
The remaining perf harnesses with large-size overrides now use the same shared
gate before heavy imports: `pauli_translation_compare.jl`,
`heisenberg_mosek_scaling.jl`, `pauli_charge_singlet_prep.jl`,
`pauli_sparse_chain_d4_blocks.jl`, and `qmbcertify_formulation_audit.jl`.
A RED source check failed on the missing `shared_load_guard.jl` include in
`pauli_translation_compare.jl`; after wiring the guard, the GREEN source-order
check passed through easy-ssh in 0.22s wall with 273,436 KiB peak RSS.  Forced
missing-estimate smokes then blocked before heavy imports for the comparison,
solver-backed scaling, and formulation-audit harnesses: 2.42s / 428,852 KiB,
2.50s / 434,544 KiB, and 2.44s / 428,792 KiB respectively.  Top-level harness
includes now use `Base.include(@__MODULE__, ...)`, so synthetic-module include
tests still work; the focused include-side-effect check passed in 1.41s wall
with 402,144 KiB peak RSS.  The narrow expectation-loader verifier then passed
2,047/2,047 checks through easy-ssh in 66.20s wall with 2,479,352 KiB peak RSS,
while running only guard/fixture checks and no solver-backed or large
construction run.
The harness usage comments were then tightened so `ALLOW_LARGE=true` is never
documented by itself as sufficient: each large-run recipe now also names
wall-time/RSS estimates and safe load/memory telemetry.  A remote source-only
header guard checked the eight guarded harnesses for estimate/telemetry wording
and trailing whitespace; it passed in 2.42s wall after an over-strict first
checker tripped only on comment line wrapping.  No package import, model
construction, or solver call was launched.
A follow-up behavior regression now proves the same policy at runtime:
`ALLOW_LARGE=true` without wall/RSS estimates is rejected for the QMBCertify
reference harness, translation profile, sparse/charge perf harnesses,
translation comparison, and solver probe.  The narrow
`test/expectations_loader.jl` verifier passed through easy-ssh in 1:21.03 wall
with 2,695,480 KiB peak RSS and no swaps, including 103/103 Pauli perf-harness
guard checks and 2,245 TOML fixture checks.  The run used one Julia thread, no
solver-backed or large construction job, and the actual wall/RSS matched the
70-85s / <3.2 GiB estimate.
The source check was then broadened to require estimate and telemetry wording
for all eight guarded harnesses, including the solver-heavy
`heisenberg_mosek_scaling.jl` and formulation-audit scripts without importing
their heavy dependencies.  The first broadened check failed only because two
comments wrapped `safe load/memory telemetry` across lines; after normalizing
comment wrapping, the same narrow verifier passed in 1:08.87 wall with
2,468,136 KiB peak RSS and no swaps, including 117/117 large-run pressure-guard
checks and the same 2,245 TOML fixture checks.
The exploratory `qmbcertify_profile_probe.jl` path now documents the same
pressure gate and has behavior coverage too: a large profile probe with
`NCTS_QMB_ALLOW_LARGE=true` but no wall/RSS estimates rejects before loading
QMBCertify.  The narrow verifier passed in 1:10.37 wall with 2,852,508 KiB
peak RSS and no swaps, including 18/18 QMBCertify harness guard checks.
Wall-time estimates are now actionable too, not merely printed.  A RED
guard-only check failed because `_ncts_large_run_estimate_status` did not
accept `max_wall_seconds`; the GREEN helper check passed in 0.46s wall with
309,496 KiB peak RSS.  A forced excessive-wall smoke then set
`NCTS_TRANSLATION_ESTIMATED_WALL_SECONDS=3600` with the default
`NCTS_LOAD_GUARD_MAX_WALL_SECONDS=1800`; it reported
`blocked_excessive_wall_estimate` before importing NCTSSoS, in 2.57s wall with
428,192 KiB peak RSS.  The expectation-loader verifier was rerun after this
policy change and passed 2,050/2,050 checks through easy-ssh in 70.79s wall
with 2,489,708 KiB peak RSS.
The memory gate now enforces shared-machine headroom as well as absolute
availability.  A RED helper check failed because
`_ncts_memory_pressure_status` did not accept `max_rss_fraction`; the GREEN
helper check passed in 0.48s wall with 310,608 KiB peak RSS.  A forced
headroom smoke then set `NCTS_TRANSLATION_ESTIMATED_RSS_GIB=600` while autodl
reported about 644 GiB available memory.  The guard reported
`blocked_insufficient_memory_headroom` before importing NCTSSoS, in 2.42s wall
with 434,112 KiB peak RSS, because the estimate would consume more than the
default `NCTS_LOAD_GUARD_MAX_RSS_FRACTION=0.8`.  This closes the gap where a
run could be below `MemAvailable` but still crowd the shared machine.
After adding the dense-Schur proxy field, a broad helper attempt that included
the full expectation loader hit its 45s timeout with 2,633,512 KiB peak RSS.
That check is deleted from the current live verification queue under the
lifted policy.  The replacement ladder is source-only, TOML-only, parser-only,
and real-source estimator-slice checks until a fresh load/RSS gate justifies
loading the full package or rerunning the full loader.
The follow-up model-size gate used that ladder: a RED source guard failed on
the missing `estimated_model_size_gate_status` field in 0.94s with 417,724 KiB
peak RSS; the GREEN source, parser, real-source estimator-slice, and TOML checks
then passed in 0.34s/297,336 KiB, 0.37s/302,404 KiB, 0.59s/329,032 KiB, and
1.69s/349,640 KiB respectively.  No package load, model construction, or solver
call was part of that verification.
Full interval certification remains deferred.

The old blanket size ceiling is lifted.  Larger probes are allowed only as
run-specific candidates, not as default queue entries: estimate wall time,
peak RSS, current load, and the open question answered before launch.  If that
estimate predicts crowding, delete or downsize the run instead of carrying it
as pending work.

Deleted larger probes, to restore only after the load/RSS gate improves, with
one thread and hard timeouts:

These entries are deleted from the active run plan until fresh evidence proves
both scientific necessity and remote safety.  Fresh evidence means current
load/RSS telemetry, a run-specific wall-time/RSS estimate, and a clear statement
of which open plan question the run answers.  If the estimate predicts crowding,
leave the entry deleted or downsize it; do not launch it and do not keep it as
queued work.
The matching fixture rows now mark these cases with
`run_queue_status = "deleted_until_evidence_gate"` plus an explicit
`restore_gate`, so "missing reviewed run" cannot silently become queued work.

1. **Resolve QMBCertify reference rows:** A1_L20 and A1_L30 are reviewed.
   A0_L20 and A2_L20 completed without machine pressure but returned
   `SLOW_PROGRESS`/`FEASIBLE_POINT`, which the explicit reviewed-status policy
   rejects.  Do not retry A0_L30 or A2_L30 as reviewed rows until an
   exploratory profile variant first demonstrates an accepted status.  Use
   `perf/qmbcertify_profile_probe.jl` with `NCTS_QMB_ARG_OVERRIDES` or
   `NCTS_QMB_ARG_OVERRIDES_FILE` for those exploratory variants before
   promoting any row to a reviewed fixture.

2. **Resolve NCTSSoS bound-quality gap:** L=20/L=30 construction-only rows,
   A2_L20/A2_L30 no-solve SOS-dual rows, A1_L20/A2_L20 solved rows, and the
   A1_L30 solved row are now recorded.  Runtime/RSS are safe for these rungs,
   but the old finite-axis A1 path was not bound-quality parity at either L=20
   or L=30.  The source-base constructor now matches QMBCertify's base block
   histogram and support quotient, including the short-orbit/overcomplete source
   basis.  After source-base phase realification, the guarded A1 L=20 and L=30
   solve replays reached objective deltas `-2.7612496822371213e-7` and
   `-1.3408209120768788e-7` against reviewed QMBCertify `A1_L20` and `A1_L30`,
   respectively, both inside the 1e-6 parity target.  The old finite-axis L=30
   gap was about `-0.339196476185879`, and the pre-phase-realification
   source-base L=30 gap was `-0.001648556199279838`.  Do not repeat A1 L=20 or
   A1 L=30 solves unless a new formulation change or regression needs them, and
   do not queue any A2 solve unless a separate no-solver/formulation reason plus
   fresh load gate justify it.  A2 remains blocked by the reviewed-status
   policy.  Do not run
   A2_L30 as a reviewed parity row until an exploratory QMBCertify A2 variant
   first demonstrates an accepted status under the explicit policy.  A
   stale-plan RED guard caught the obsolete remaining
   base-mismatch instruction in 1.83s with 417,508 KiB RSS; the rewritten
   source-only guard passed in 0.21s with 260,420 KiB RSS.  A follow-up
   stale-plan RED guard caught the older "remaining dominant formulation
   mismatch" wording in 1.85s with 417,816 KiB RSS; the corrected historical
   wording passed the same source-only guard in 0.21s with 260,532 KiB RSS.
   A second stale-plan RED guard caught the obsolete A2_L30 base-mismatch
   dependency in 1.87s with 418,228 KiB RSS; the corrected reviewed-status plus
   solve-risk dependency passed the source-only guard in 0.22s with
   260,008 KiB RSS.

3. **Re-run scaling gates only under a fresh load gate:** the immediate L=30
   constructed scaling profile is removed from the active run queue on
   2026-07-07 because `autodl` reported `nproc=25` and load average
   `30.79, 40.25, 44.20`.  That evidence says even a one-thread Julia run
   would add pressure to an already overcommitted machine.  Re-add a large
   constructed profile only after a fresh `uptime`/`free -h` check shows enough
   CPU and memory headroom, and still run with one thread plus a hard timeout.
   The SU(2) RDM/LSO/PSO N=100 target-only gate has already passed; any
   further N=100 probe stays `NCTS_TRANSLATION_TARGET_ONLY=true` until a
   constructed profile has fresh safe-run evidence.
   The L=30 SU(2) native no-solve model-sizing row is safe and keeps the
   dense-Schur proxy at about 120 GiB, but the first native L=30 SU(2) solved
row timed out after 30:03.89 with 81,381,828 KiB peak RSS.  Do not repeat
that solved row until a formulation change reduces the 122,385 scalar
equalities or otherwise changes the solver memory/time model.  The timeout
did not stall the remote, but it failed the Phase 3 runtime gate.
After this SU(2) timeout was recorded and the solved row was deleted from the
active queue, the narrow expectation loader passed in 1:10.20 wall with
2,501,168 KiB peak RSS, including 14/14 parity-plan run gates and 2,051
fixture checks.
The SOS-dual coefficient-equation builder now skips structurally zero affine
coefficient rows, and the solver-probe model-size estimator now counts two real
moment variables for complex moment data lowered as `formulation =
:moment_variables, representation = :real`.  A guarded N=8 native-Hermitian
SU(2) no-solve check verified the zero-row skip without solving: it completed
in 1:44.08 wall with 2,991,520 KiB peak RSS and built 8,139 actual scalar
equalities, down from the conservative 8,195 coefficient rows.  This correctness
fix is real but too small to explain the L=30 timeout.
The first primal moment-variable SU(2) no-solve checks were more promising for
the equality bottleneck.  At N=8, the primal no-solve row completed in 1:34.64
wall with 4,497,772 KiB peak RSS, 8,196 moment variables, 821 scalar
equalities, and dense-Schur proxy 5,392,328 bytes.  At L=30, the matching
primal no-solve row completed in 3:55.11 wall with 10,108,084 KiB peak RSS,
122,386 moment variables, 821 scalar equalities, and the same 5,392,328-byte
dense-Schur proxy, compared with the native SOS-dual row's 122,385 scalar
equalities and 119,824,705,800-byte proxy.  However, the guarded N=8 primal
solve then took 4:25.47 wall with 11,420,812 KiB peak RSS, returning
`OPTIMAL` objective `-3.6510929335842697`, within about `4.8e-7` of the
existing N=8 SU(2) native-dual objective.  The primal formulation is therefore
numerically consistent and removes the scalar-equality bottleneck, but its
small-N solve time is worse than the native SOS-dual solve.  Do not jump
directly to an L=30 primal solve; first run an intermediate-size primal solve
under a fresh load/RSS gate to establish scaling.
The intermediate N=12 primal moment-variable no-solve row passed that gate: it
completed in 2:36.27 wall with 9,409,188 KiB peak RSS, 56,488 moment
variables, 821 scalar equalities, 65 PSD blocks, max PSD block 31, and the same
5,392,328-byte dense-Schur proxy.  The corresponding guarded N=12 primal solve
then timed out at 20:02.64 wall with 26,226,188 KiB peak RSS and no swaps,
inside the remote memory-safety envelope but outside the useful scaling
envelope.  Delete N >= 14, L=20, and L=30 primal solved probes from the active
run queue under the current formulation.  Keep only no-solve sizing probes or
smaller numerical smoke checks until another formulation change reduces Mosek's
factorization cost.
The solver probe harness now enforces that deletion: solved primal
`formulation = :moment_variables` SU(2) RDM/LSO/PSO probes at N >= 14 are
refused by default before construction starts.  No-solve sizing, SOS-dual
probes, N <= 12 smoke solves, and an explicit post-formulation-change override
remain available, but the override must be paired with fresh wall/RSS/load
evidence.
The persisted fixture loader now enforces the same stop signal: a solved
primal `moment_variables` SU(2) RDM/LSO/PSO row at `length >= 14` is invalid,
while no-solve sizing rows and smaller smoke solves remain valid.  A RED
synthetic fixture contract failed in 10.94s with 855,800 KiB peak RSS because
the deleted solved row was accepted; after adding the loader guard, the narrow
loader passed in 1:11.02 wall with 2,342,820 KiB peak RSS and no swaps after a
preflight gate of load `10.48, 13.67, 14.39` and 652 GiB available memory.
The probe fixture emitter now also records `formulation` and `representation`
for future rows.  A tiny N=4/order-1 construct-only fixture smoke emitted
`formulation = "psd_blocks"` and `representation = "complex"` in 1:15.27 wall
with 2,267,460 KiB peak RSS and no swaps; no model lowering or solver call was
performed.
The N=12 primal timeout row now also carries `run_queue_status =
"deleted_until_evidence_gate"` and
`restore_gate = "post_formulation_change_and_fresh_wall_rss_load_estimate"`,
so the larger deleted primal solved probes are fixture-visible rather than
only described in prose.  The loader contract requires timeout rows to carry
that deletion metadata.  The narrow verifier passed 2,713/2,713 checks in
1:12.72 wall with 2,497,536 KiB peak RSS and no swaps after a preflight gate of
load `10.90, 14.26, 14.21` on 25 CPUs and 652 GiB available memory.
The QMBCertify reviewed-reference harness now enforces the same deleted-queue
policy for A0/A2 rows: if a requested pending case has
`run_queue_status = "deleted_until_evidence_gate"`, the harness refuses the
run unless `NCTS_QMB_ACCEPTED_PROFILE_VARIANT_ID` names the accepted exploratory
profile variant that restored it.  The regression keeps reviewed A1 rows
unchanged and verifies that A0/A2 deleted rows reject by default.  The narrow
loader passed 2,713/2,713 checks in 1:10.41 wall with 2,767,460 KiB peak RSS
and no swaps after a preflight gate of load `12.25, 13.79, 14.06` on 25 CPUs
and 654 GiB available memory.
The restore override is now fixture-checked rather than string-checked: the
named evidence row must exist in `profile_probes`, have `status = "probe_run"`,
match the requested `parent_profile`, carry an accepted termination status, and
name an actual exploratory variant with non-`none` `overrides_summary`.
Because the restore gate also depends on fresh run-risk evidence, the same
probe must carry positive `estimated_wall_seconds`, `estimated_rss_gib`,
`total_wall_seconds`, and `peak_rss_bytes`.  New QMBCertify reference/profile
rows emitted by the harness print the estimate fields from the same env vars
used by the large-run pressure guard.  The regression has been extended to
reject missing evidence, wrong-profile evidence, the existing default A2 L=12
probe, a synthetic override probe with missing RSS evidence, and a synthetic
override probe with missing estimate evidence; a synthetic accepted A2 override
probe exercises the positive path without adding fake reviewed data to the
fixture.  After autodl cooled from load `21.48, 17.39, 15.13` to
`15.70, 16.42, 14.90`, the estimate-field tightening was verified by the
narrow loader: 2,713/2,713 checks passed in 1:11.95 wall with 2,806,380 KiB
peak RSS and no swaps.  The earlier
intermediate loader gate passed 2,713/2,713 checks in 1:10.72 wall with
2,497,844 KiB peak RSS and no swaps after a preflight gate of load
`12.92, 11.96, 12.95` and 673 GiB available memory, before the override-summary
tightening.  The stricter guard is verified by the later loader run below.
The exploratory QMBCertify profile harness now blocks deleted large profile
probes without `NCTS_QMB_ARG_OVERRIDES` or `NCTS_QMB_ARG_OVERRIDES_FILE`; a
large probe must be an actual variant, not a disguised rerun of a rejected
reviewed row.  The regression accepts reviewed A1 large rows and small A2
probes without overrides, rejects A2_L20 with no overrides, and accepts A2_L20
when a variant override is declared.  The narrow loader passed 2,713/2,713
checks in 1:12.10 wall with 2,812,076 KiB peak RSS and no swaps after a
preflight gate of load `12.24, 12.79, 13.06` and 653 GiB available memory.
No QMBCertify profile solve or N>=14 NCTSSoS run was launched.
After the override-summary tightening, the same narrow loader passed
2,713/2,713 checks in 1:10.98 wall with 2,845,712 KiB peak RSS and no swaps
after a preflight gate of load `10.77, 12.49, 12.93` and 652 GiB available
memory.  During the loader's synthetic pressure-gate probes the one-minute load
briefly reached `26.2`; no substantive run was launched under that higher load.
After requiring positive wall/RSS fields on restore evidence, the same narrow
loader passed 2,713/2,713 checks in 1:11.95 wall with 2,786,764 KiB peak RSS
and no swaps after a preflight gate of load `13.03, 14.09, 13.55` and 651 GiB
available memory.  During that loader's synthetic pressure-gate probes the
one-minute load reached `28.08`, so no heavier verifier or performance run was
queued afterward.

4. **Re-verify after fixture promotion:** the earlier dense-Schur fixture update
   used only narrow source/TOML/parser/slice checks because a broad loader
   include had timed out at 45.47s and 2.63 GiB peak RSS under worse load.  After
   fresh telemetry on 2026-07-08 00:35 CST reported load averages
   `16.34, 15.12, 15.56` on 25 CPUs with 696 GiB available memory, the full
   `test/expectations_loader.jl` gate was re-added with one Julia thread and a
   150s hard timeout.  It passed in 1:06.58 wall with 2,487,556 KiB peak RSS,
   including the harness guard testsets and 1,855 fixture checks.  Keep
   `test/relaxations/runtests.jl` and package-level `Pkg.test()` deferred until
   a fresh load/RSS gate justifies those longer one-thread runs.  These checks
   prove regression coverage; they do not replace the reference/runtime evidence
   from the larger probes.  After the deleted-run-plan wording was added, a
   2026-07-08 00:39 CST telemetry check reported load averages
   `9.48, 12.75, 14.44` on 25 CPUs with 695 GiB available memory.  That was
   safe for the same loader gate, which passed in 1:07.78 wall with
   2,316,892 KiB peak RSS.  It was not treated as a reason to restore L=20/L=30
   solves because no new scientific target or model-size evidence required
   those runs.  A follow-up loader guard now enforces execution-state metadata
   across all NCTSSoS solver/model/construction probe fixture arrays:
   `construction_only`, `model_built`, `solved`, and `termination_status` must
   be explicit, solved rows must carry an objective and solve time, and timeout
   rows must record `timeout_exit_status = 124` plus the deleted larger solved
   probe lengths.  After fresh telemetry reported load averages
   `19.87, 12.94, 13.19` on 25 CPUs with 704 GiB available memory, the same
   one-thread loader gate passed in 1:11.29 wall with 2,485,508 KiB peak RSS
   and no swaps.  This strengthens the deleted-run policy only; it does not
   restore any L=20/L=30/N>=14 solved run.
   A full current-tree verification sweep then ran under the same lifted
   load-gated policy with one Julia thread and hard timeouts.  The narrow
   expectation loader passed in 1:10.39 wall with 2,495,876 KiB peak RSS and no
   swaps; `test/relaxations/pauli_chains.jl` passed in 7:56.74 wall with
   4,618,772 KiB peak RSS and no swaps; `test/relaxations/moment_linear.jl`
   passed in 2:15.88 wall with 1,414,224 KiB peak RSS and no swaps; the touched
   relaxation batch (`interface.jl`, `lowering.jl`, `sos.jl`, `sparsity.jl`,
   `symmetry.jl`, and `gns_pipeline.jl`) passed in 7:20.59 wall with
   2,441,944 KiB peak RSS and no swaps.  Package-level `Pkg.test()` was then
   launched only after telemetry reported load averages `12.60, 14.22, 14.18`
   on 25 CPUs with 678 GiB available memory; it passed 14,239 checks with 2
   broken checks in 26:03.65 wall, 5,543,700 KiB peak RSS, and no swaps.  The
   package run matched the 24-28 minute / <6 GiB estimate.  These are
   regression and policy-gate checks only; they do not restore any deleted
   L=20/L=30/N>=14 solved or profiling run without a new scientific question,
   fresh telemetry, and a run-specific wall/RSS estimate.

An earlier telemetry check on 2026-07-08 00:16 CST took 2.15s and reported load
averages `24.32, 19.02, 17.65` on 25 CPUs with 631 GiB available memory.  That
was acceptable for the short source/TOML guards above, but not enough margin to
restore multi-minute L=20/L=30 solves, constructed scaling profiles, broad
loader checks, or package-wide tests to the live queue at that point.  Later
loader checks were restored only after fresher telemetry showed better
headroom; the L=20/L=30/N>=14 solved runs remained deleted.

Large-run command templates.  These templates are historical command recipes,
not a run queue.  Use them only after the gates change and the corresponding
entry is explicitly restored:

```bash
# QMBCertify reference/status probes.  Set NCTS_NCTSSOS_COMMIT to the local
# NCTSSoS commit being compared; easy-ssh does not guarantee a remote .git.
# A1_L20/A1_L30 are already reviewed; A0/A2 larger rows first need an
# exploratory variant that demonstrates an accepted reviewed status after
# SLOW_PROGRESS/FEASIBLE_POINT at L=20.
NCTS_QMB_ALLOW_LARGE=true \
NCTS_QMB_ENVIRONMENT_ID=<reviewed-environment-id> \
NCTS_NCTSSOS_COMMIT=<local-nctssos-commit> \
NCTS_QMB_PROFILE=<A0-or-A2> \
NCTS_QMB_NS=<20-or-30> \
NCTS_QMB_MOSEK_THREADS=1 \
NCTS_QMB_ESTIMATED_WALL_SECONDS=<estimate-before-launch> \
NCTS_QMB_ESTIMATED_RSS_GIB=<estimate-before-launch> \
julia --project=. --startup-file=no perf/qmbcertify_reference_runs.jl

# Use perf/qmbcertify_profile_probe.jl plus NCTS_QMB_ARG_OVERRIDES for
# exploratory variants that should not yet become reviewed cases.

# NCTSSoS source-base or finite-axis source-like rows.  Use
# Clarabel-labeled no-solve model sizing in the root project; use the docs
# project only when actually solving with Mosek.  For the QMBCertify
# source-base route, keep N=20 construction-only as the first valid rung; do
# not use N=6/N=8 source-base probes as validation evidence.
NCTS_SOLVER_PROBE_ALLOW_LARGE=true \
NCTS_SOLVER_PROBE_N=20 \
NCTS_SOLVER_PROBE_ORDER=4 \
NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_CONSTRUCT=true \
NCTS_SOLVER_PROBE_QMBCERTIFY_BASE_EXTRA=9 \
NCTS_SOLVER_PROBE_LOWER_MODEL=false \
NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true \
NCTS_SOLVER_PROBE_EMIT_FIXTURE=true \
NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=<estimate-before-launch> \
NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=<estimate-before-launch> \
julia --project=. --startup-file=no perf/pauli_translation_solver_probe.jl

# Finite-axis source-like model-sizing template.  Do not promote to a solve
# until the static estimates are reviewed and the no-solve model-sizing row
# passes a fresh load/RSS gate.  If `estimated_lowering_would_error=true`,
# fix the selected formulation/orphan policy before launching lowering.
# If `estimated_risk_gate_status` is `blocked_lowering_orphan_policy`, fix the
# selected formulation/orphan policy or choose the SOS-dual route before
# launching lowering.
# If `estimated_model_size_gate_status` is not `ok`, delete or downsize the
# next lowering/solve rung before launching it.
NCTS_SOLVER_PROBE_ALLOW_LARGE=true \
NCTS_SOLVER_PROBE_N=30 \
NCTS_SOLVER_PROBE_ORDER=4 \
NCTS_SOLVER_PROBE_DUALIZE=true \
NCTS_SOLVER_PROBE_OPTIMIZER=clarabel \
NCTS_SOLVER_PROBE_SOLVE=false \
NCTS_SOLVER_PROBE_EMIT_FIXTURE=true \
NCTS_SOLVER_PROBE_REFLECTION=true \
NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=true \
NCTS_SOLVER_PROBE_AXIS_SYMMETRY=true \
NCTS_SOLVER_PROBE_AXIS_EQUALITIES=true \
NCTS_SOLVER_PROBE_AXIS_QUOTIENT=true \
NCTS_SOLVER_PROBE_ESTIMATED_WALL_SECONDS=<estimate-before-launch> \
NCTS_SOLVER_PROBE_ESTIMATED_RSS_GIB=<estimate-before-launch> \
julia --project=. --startup-file=no perf/pauli_translation_solver_probe.jl

# For selected Mosek solves after model sizing passes, switch to:
NCTS_SOLVER_PROBE_OPTIMIZER=mosek \
NCTS_SOLVER_PROBE_MOSEK_THREADS=1 \
NCTS_SOLVER_PROBE_SOLVE=true \
julia --project=docs --startup-file=no -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); include("perf/pauli_translation_solver_probe.jl"); main()'

# A1 deltas: insert these env vars before the `julia ... solver_probe` command
# above to add QMBCertify RDM blocks and width-7 linear state-opt.
NCTS_SOLVER_PROBE_RDM_K=8 \
NCTS_SOLVER_PROBE_RDM_DECOMPOSITION=qmbcertify \
NCTS_SOLVER_PROBE_RDM_SUPPORT=extend \
NCTS_SOLVER_PROBE_LSO_WIDTH=7 \
NCTS_SOLVER_PROBE_LSO_MODE=qmbcertify

# A2 delta: use the A1 deltas plus this env var before the same command.
NCTS_SOLVER_PROBE_PSO_WIDTH=3

# SU(2) memory gate.  Compare real-lift and native Hermitian variants before
# promoting a large solved row.  Insert these env vars before the same command.
NCTS_SOLVER_PROBE_SU2=true \
NCTS_SOLVER_PROBE_BASE_SU2_EXTEND_RDM=true \
NCTS_SOLVER_PROBE_SU2_MOMENT_QUOTIENT=true \
NCTS_SOLVER_PROBE_RDM_K=8 \
NCTS_SOLVER_PROBE_RDM_DECOMPOSITION=su2 \
NCTS_SOLVER_PROBE_RDM_SUPPORT=extend \
NCTS_SOLVER_PROBE_LSO_WIDTH=7 \
NCTS_SOLVER_PROBE_PSO_WIDTH=3 \
NCTS_SOLVER_PROBE_REAL_MOMENT_MATRIX=false \
NCTS_SOLVER_PROBE_REFLECTION=false \
NCTS_SOLVER_PROBE_SOS_HERMITIAN=native \
NCTS_SOLVER_PROBE_SOS_EQUATION_SCALING=max_abs \
NCTS_SOLVER_PROBE_SOS_EQUATION_SCALING_FLOOR=0.1

# N=100 structural gate: target-only first, no constructed model or solve.
NCTS_TRANSLATION_ALLOW_LARGE=true \
NCTS_TRANSLATION_TARGET_ONLY=true \
NCTS_TRANSLATION_NS=100 \
NCTS_TRANSLATION_ORDER=4 \
NCTS_TRANSLATION_SU2=true \
NCTS_TRANSLATION_RDM_K=8 \
NCTS_TRANSLATION_RDM_DECOMPOSITION=su2 \
NCTS_TRANSLATION_RDM_SUPPORT=extend \
NCTS_TRANSLATION_LSO_WIDTH=7 \
NCTS_TRANSLATION_PSO_WIDTH=3 \
NCTS_TRANSLATION_ESTIMATED_WALL_SECONDS=<estimate-before-launch> \
NCTS_TRANSLATION_ESTIMATED_RSS_GIB=<estimate-before-launch> \
julia --project=. --startup-file=no perf/pauli_translation_profile.jl
```

**Actions:**

1. **Streaming/DFT path records (hook in place):** block-to-basis mapping,
   DFT/CG transform descriptors, equality multiplier provenance, constraint
   origins

2. **Dual block extraction helper (hook in place):** given solved model,
   extract reduced dual PSD blocks Q_i with block labels

3. **Coefficient-domain documentation (hook in place):** original Pauli basis
   (rational), DFT basis (algebraic/irrational), CG basis (algebraic).  Naive
   rational rounding of DFT/CG coefficients is NOT valid for certification.

4. **Small-N SOS identity check (numerical hook in place):** extract dual
   blocks, verify SOS identity H − λ = Σ v_i† Q_i v_i + (residual) to solver
   tolerance

**Coefficient-domain note:** A future certificate extractor must keep three
domains separate.  Pauli word coefficients in the original operator basis are
rational up to the conventional `1/4` Hamiltonian normalization.  Momentum
blocks introduce DFT phases, so their row transforms live in cyclotomic
extensions (`cos(2πk/N)`, `sin(2πk/N)`), even when the realified solver block
contains only real floating-point coefficients.  SU(2) blocks will introduce
Clebsch-Gordan coefficients, typically square-root algebraic numbers.  A
certificate may round the final residual with interval bounds, but it must not
pretend DFT or CG transform coefficients are rational data.

Full certification (interval arithmetic, Arblib, certified residual bounds) is
deferred to a separate plan.

### Phase 3 completion evidence — 2026-07-11

The exact singlet quotient is now implemented on the direct translation/SU(2)
path.  It uses deterministic thin Clebsch-Gordan singlet rows through support
eight, exact sparse-integer rank computation for projected coordinate masks,
conditioned canonical pivot charts, sparse rewriting of the complete
`MomentLinearData`, and cone-only Wigner-Eckart reduction.  Direct and custom
inputs reject incomplete coordinate support.  The generated translation path
may use a projected mask only with explicit descriptor and per-orbit
provenance.  Quotient mode does not materialize the raw Wigner or singlet
equality families; it preserves the reduced PSD cones and records the quotient
chart and residual diagnostics.

The load-gated ladder used the complex/native SU(2)+k=8 RDM+LSO7+PSO3 profile,
one Julia thread, one Mosek thread, and no gate override.  The accepted L=30
run emitted fixture
`NCTSSOS_A2_SU2_MOMENT_QUOTIENT_L30_source_like_sos_dual_mosek_2026_07_11`
with these measurements:

- raw moments `94_129`, quotient moments `3_529` (3.7491% retained);
- `185` native Hermitian cones, maximum block `28`;
- `12_001` dual variables and `4_528` scalar equalities;
- Mosek `OPTIMAL`, primal and dual `FEASIBLE_POINT`, relative gap
  `1.4046251926234205e-8`;
- objective `-13.324596480957192`;
- SOS coefficient residual `3.3181295178152936e-7`, minimum native eigenvalue
  `-1.9242165196601706e-9`;
- construction `421.477801331s`, solve `281.124462594s`, harness wall
  `778.215219242s`, remote wall `821.14s`;
- peak RSS `12_801_622_016` bytes, no swap, inside the 15 GiB estimate.

The objective differs by `5.96901113e-7` from the earlier exact same-profile
L=30 feasible result `-13.324597077958305`, satisfying the `1e-6`
symmetry-equivalent objective gate.  The older result stalled; the accepted
row above is the reviewed baseline because capped row scaling plus the recorded
Mosek tolerances returned strict `OPTIMAL` while preserving a sub-`1e-6`
certificate residual.  The A1 QMBCertify value is not used as an objective
reference: A1 omits the A2 PSO profile and is not symmetry-equivalent.

The accepted wall is more than twice as fast as the recorded 30-minute Phase-2
NCTSSoS run.  N=8, N=10, N=12, and L=20 rungs established the safe progression;
the final loader independently checks the L=30 quotient counts, solver status,
gap, residual, eigenvalue, model size, wall, and RSS.  Phase 4 interval
certification and Phase 5 two-dimensional work remain separate projects and
are not claimed here.

After independent review, adversarial regressions were added for missing
patterns, projected-mask provenance, singular charts, row permutation,
mislabeled nonzero invariant rows, and concurrent Clebsch-Gordan cache access.
The corrected implementation rewrites every candidate zero row before dropping
it and protects shared rank/CG caches with locks.  A final exact-rank L=30
no-solve/dualization recheck reproduced the solved fixture's model exactly:
`94_129` raw moments, `3_529` quotient moments, 185 cones (maximum block 28),
`12_001` dual variables, and `4_528` scalar equalities.  Construction took
420.5785s; remote wall was 8:48.20 with 13,071,404 KiB peak RSS and no swap.
Because the solver model is identical, the strict-`OPTIMAL` fixture above
remains the authoritative numerical result.  The focused quotient suite passed
2,804/2,804 checks, the independent re-review reported no actionable findings,
and the final package run passed 17,378 checks with two intentional broken
checks and zero failures in 26:08.0 test time (26:31.28 wrapper wall), using
6,287,336 KiB peak RSS and no swap.

**Milestones:**
- Dual extraction works for all Phase 1-3 model variants
- Small-N SOS identity verification passes to solver tolerance

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
| Cyclic DFT alone fails to reach max block 31 (need reflection/sign/axis-permutation too) | High | High | Phase 1A includes full discrete symmetry profile; validate block histogram, not just max block |
| RDM constraints move bottleneck from construction to solver memory | High | High | Measure per-stage; don't attempt N=100 until L=30 extrapolation is safe |
| SU(2) block decomposition errors give wrong bounds | Medium | Critical | Generate from CG coupling with tests, never hardcode; exact/interval CG coefficients |
| Specialized path silently differs from generic | Medium | Critical | Mandatory small-N equivalence tests before opt-in; test at N > 2d to avoid short-orbit collisions |
| Naive product-table cache blows memory at 12001² pairs | Medium | High | Cache at orbit level, not basis-pair level; monitor cache size |
| State-opt closure rules are wrong or drift from QMBCertify | Medium | High | Keep the source-pinned Araujo/Fawzi formula in tests; verify any future variant against QMBCertify source first |
| DFT realification bugs at k=0, k=N/2, or short-orbit stabilizers | Medium | High | Explicit test cases for each special momentum sector |
| Fast-path recognition false positives (wrong dispatch) | Medium | High | Precise recognizer spec with rejection tests; no accidental fast-path activation |
| Logical vs realified block-size mismatch | Medium | Medium | Track and report both; existing path shows 62 realified vs 31 logical |
| Thread-local simplification buffers introduce aliasing/nondeterminism | Low | High | Explicit key-copy ownership rules; mutation-corruption tests |
| Streaming bypasses MomentProblem, breaks dual/certification provenance | Medium | High | Phase 4 provenance hooks required before streaming can replace symbolic path |
| DFT/CG irrational coefficients prevent rational certification | High | Medium | Document coefficient domain; plan interval arithmetic, not naive rounding |
| N=100 solver memory exceeds available RAM even with all reductions | Medium | High | Phase 3 (SU(2)) must succeed before N=100 attempt |
| U(1)/SU(2) assumptions exclude true ground state for some N/boundary/model combinations | Low | High | Document which Hamiltonians/boundaries are covered; test edge cases |
| Integration risk with existing pauli_chains.jl specialized path | Medium | Medium | Extend existing path rather than building parallel infrastructure |
| 2D D4 symmetry requires momentum-star theory, not simple DFT | High | Medium | Defer to Phase 5; Phase 1 designs API to accommodate 2D later |
| QMBCertify parity target mismatched due to normalization, tolerances, or option differences | Medium | Medium | Phase 0 extracts exact options; TOML reference files pin all parameters |

---

## Part VIII: Priority Summary

| Priority | Phase | What it buys | Expected speedup (standalone) |
|:--|:--|:--|:--|
| 1 | Phase 0 | Falsifiable targets; constraint inventory | — |
| 2 | Phase 1A | Eliminate Wedderburn bottleneck (~80% of cost) | ~5× |
| 3 | Phase 1B | Eliminate redundant product computation | ~1.5× incremental |
| 4 | Phase 1C | Eliminate remaining polynomial intermediate | ~1.3× incremental (only if profiling justifies) |
| 5 | Phase 2A-C | Bound quality parity with QMBCertify | N/A (accuracy, not speed) |
| 6 | Phase 3 | Solver memory reduction via SU(2); needed for N=100 | Depends on block reduction |
| 7 | Phase 4 | Certification provenance (hooks only) | N/A |
| 8 | Phase 5 | 2D lattice (separate project) | N/A |

Phases 1A–1C together target ~50-200× construction speedup (profiling needed
to confirm).  Phases 2A-C close the bound quality gap.  Phase 3 is required
for N=100 numerical runs.  Phase 4 preserves the option for future rigorous
certification.  Phase 5 is a separate research effort.
