# Performance gap against arXiv:2604.01555v1

Reference source: `/Users/exaclior/MyBrain/YushengBrain/references/sources/2604.01555v1/`.

## Verdict

The current generic symmetry path is **not on-par** with the paper for the 1D sparse degree-4 Heisenberg-chain target. It gets the right qualitative block-size reduction, but it spends far too much time and allocation constructing that reduction.

The paper/QMBCertify implementation uses a specialized translation-orbit + DFT construction. NCTSSoS currently pushes the sparse chain through generic Clifford/Wedderburn machinery, which scales badly with chain length.

## Paper target

For the 1D Heisenberg chain sparse basis with `r=1`, `d=4`, `N=100`, the paper reports:

| stage | maximal block / basis size |
|:--|--:|
| original SDP | 8,127,090,301 |
| after Pauli equalities | 322,029,976 |
| after sparse contiguous basis | 12,001 |
| after symmetry | 31 |

Relevant paper mechanisms:

- sparse contiguous monomial basis;
- sign/conjugate/permutation/mirror symmetries;
- translation block diagonalization by explicit DFT;
- realification/conjugate-pair handling for momentum sectors;
- additional RDM/state-optimality strengthening for final numerical bounds.

## Current NCTSSoS measurements

Remote host: `autodl` via `easy-ssh`.
Julia: 1.12.6.
CPU reported: Intel Xeon Platinum 8470Q.
Solver calls: none for structural benchmarks.

### Existing order-2 charge/spatial/singlet prep path

Command artifact: `perf/results/pauli_charge_singlet_prep_N8_12_16.md`.

| N | dense half-basis | group order | PSD blocks | largest block | reduced PSD scalar vars | `moment_relax_symmetric` |
|--:|--:|--:|--:|--:|--:|--:|
| 8 | 277 | 16 | 33 | 12 | 690 | 3.708 s |
| 12 | 631 | 24 | 43 | 18 | 2,204 | 7.297 s |
| 16 | 1,129 | 32 | 53 | 24 | 5,108 | 22.369 s |

This path is acceptable for small order-2 tests.

### Sparse degree-4 chain structural path

Command artifacts:

- `perf/results/pauli_sparse_chain_d4_blocks_N16_32.md`
- `perf/results/pauli_sparse_chain_d4_blocks_N64.md`

Those artifacts include an `N=4` warmup. Ignore its printed large-N sparse-basis formula comparison; periodic windows collide when `N ≤ d`, so the paper formula is only meaningful for the measured `N=16,32,64` rows below.

| N | sparse basis | group order | PSD blocks | largest block | target max block | decomposition time | allocations |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 16 | 1,921 | 64 | 95 | 30 | 31 | 8.781 s | 1.08 GiB |
| 32 | 3,841 | 128 | 167 | 30 | 31 | 75.713 s | 7.17 GiB |
| 64 | 7,681 | 256 | 311 | 30 | 31 | 1,121.018 s | 88.08 GiB |

The block sizes are in the right ballpark. The construction cost is not.

I did **not** submit `N=100`: based on `N=64`, the generic decomposition would be a long, allocation-heavy run before any SDP model or Mosek solve. That would measure stubbornness, not engineering.

## Root cause

The original charge-sector transform built dense charge rows. That was fixed by storing charge rows sparsely before block multiplication.

Remaining bottleneck: `_pauli_charge_transform_groups` still delegates each charge sector to generic `SymbolicWedderburn.symmetry_adapted_basis` for the full spatial/sign group. For a periodic chain this is the wrong abstraction. The reference code does not ask a generic representation package to rediscover Fourier modes; it constructs orbit representatives and DFT blocks directly.

## Changes made here

- Added sparse charge-transform construction in `src/optimization/symmetry.jl`.
- Added sparse-aware `_sparse_transform_rows(::SparseMatrixCSC)`.
- Added structural benchmark harness: `perf/pauli_sparse_chain_d4_blocks.jl`.
- Captured remote benchmark outputs under `perf/results/`.

## 2026-08-13: Mosek formulation memory follow-up

Phase-boundary profiling of the specialized TI DFT path found that the dual
SOS formulation is a Mosek memory hog for this problem family. With
`dualize=true`, `sos_dualize` creates one scalar equality per unique moment,
coupling entries across many small PSD-variable blocks. At `N=100`, this is
about 46,000 equalities over about 200 PSD blocks of size at most 30. Mosek's
solve-phase peak for this form grows approximately as `N³`: controlled runs at
`N=16,24,32` increased peak RSS by 0.7, 1.5, and 7.4 GiB and took 12, 33, and
103 seconds, respectively.

The primal moment formulation (`dualize=false`) instead uses free scalar moment
variables with PSD LMI constraints. Its solve-phase memory is approximately
linear in chain length at order 4, about 0.07 GiB per site over the measured
range. All moment-form rows below used fresh processes, Mosek with 16 threads,
and feasibility and relative-gap tolerances of `1e-8`.

| N | dual SOS bound | dual peak RSS | dual time | moment LMI bound | moment peak RSS | moment time |
|--:|--:|--:|--:|--:|--:|--:|
| 32 | -14.2306124253 | 8.11 GiB (+7.40) | 102.9 s | -14.2306029477 | 2.69 GiB (+1.98) | 48.1 s |
| 48 | not run | not run | not run | -21.3325283611 | 4.08 GiB (+3.37) | 69.8 s |
| 64 | not run | not run | not run | -28.4368184925 | 5.31 GiB (+4.60) | 110.0 s |
| 100 | -44.4239941277 | ~270 GiB | 1640.5 s | -44.4239799523 | 8.98 GiB (+8.27) | 237.5 s |

The `N=100` headline is therefore **~270 GiB and 1640 seconds → 9.0 GiB and
237.5 seconds**. The two bounds differ by `1.4e-5` absolute, or `3.2e-7`
relative, which is at solver-tolerance scale.

The reusable profiler in `perf/ti_memory_profile.jl` records the symbolic
build, SOS dualization, symbolic-data drop, and Mosek solve boundaries using
`VmRSS`/`VmHWM` from `/proc/self/status` and `Base.summarysize` attribution.

This is not Julia-side retention: at `N=32`, the complete `MomentProblem`
(including constraint and linear-form caches) occupied 0.053 GiB by
`Base.summarysize`, while the JuMP SOS model occupied 0.016 GiB, and dropping
all symbolic data before `optimize!` freed effectively no RSS. The leverage is
the solver formulation, not retaining fewer symbolic structures or reducing
Mosek from 16 to 4 threads. The TI entry point already defaults to
`dualize=false`; the Mosek performance harness now follows that default while
retaining an environment knob for explicit A/B runs.

## Verification run

Remote via `easy-ssh`, Mosek preferred where solver-backed tests were needed.

| check | result | actual runtime |
|:--|:--|--:|
| `test/relaxations/pauli_chains.jl` | pass | ~16 s after precompile |
| `test/problems/condensed_matter/heisenberg_symmetry.jl` with explicit Mosek `SOLVER` | pass | ~67 s |
| `perf/pauli_sparse_chain_d4_blocks.jl`, N=4 smoke | pass | <1 min |

## Recommended next implementation step

Do **not** keep trying to tune the generic Wedderburn path for this paper target.

Add an explicit opt-in 1D chain specialization, e.g. `PauliChainSymmetrySpec`, with strict validation:

1. Requires `pauli_contiguous_chain_basis(..., d; periodic=true)`-compatible basis.
2. Requires translation/reflection/sign-compatible chain symmetries.
3. Builds translation orbits by relative support pattern.
4. Constructs DFT momentum-sector transforms directly.
5. Handles `k=0`, `k=N/2`, and conjugate momentum pairs explicitly.
6. Adds small-`N` equivalence tests against the generic path before using it at scale.

Success gates for claiming parity with the paper:

- `N=100,d=4` sparse basis size `12,001`.
- largest PSD block exactly `31` (or a documented stronger reduction).
- structural decomposition time comfortably below a few seconds to low tens of seconds, not tens of minutes.
- end-to-end moment/SOS model construction measured separately from Mosek solve time.
- final bounds checked against the paper tables for at least one small case before scaling.

## 2026-08-13: TI RDM and state-optimality strengthening

The specialized TI path now has three opt-in strengthening mechanisms from
arXiv:2604.01555:

- `rdm_levels=[k]` adds the contiguous `k`-site reduced-density-matrix PSD
  constraint, split into fixed-magnetization blocks.  Its full complex Pauli
  representation is assembled directly in those blocks and realified exactly.
- `state_optimality=:linear` adds supported, symmetry-inequivalent instances of
  `ℓ([H,u]) = 0` for the full periodic Hamiltonian.
- `state_optimality=:linear_psd` also adds the PSD state-optimality matrix on
  the degree-3 contiguous basis plus two-site words through
  `state_optimality_range`.  This remains opt-in because the paper reports
  numerical trouble from this condition for some J1-J2 models.

The implementation also has an opt-in global Pauli-axis permutation quotient
for the isotropic XXX objective.  Together with sign and conjugation symmetry,
it identifies equivalent x/y/z moments and keeps the trivial sign sector plus
one representative nontrivial sector.  The objective and supplied basis are
checked for invariance before this quotient is used.

The Hamiltonian convention agrees with the paper: the default
`heisenberg_chain_hamiltonian` coupling is `1/4`, so it returns
`Σᵢ Sᵢ⋅Sᵢ₊₁ = (1/4)Σᵢ σᵢ⋅σᵢ₊₁`.  The solver objective is the full periodic
chain energy and the table below divides it by `N`.  For example, the exact
`N=10` energy is `-4.515446354492043`, or `-0.4515446354492043` per spin.

QMBCertify's published setup uses a 10-site RDM, full-chain linear
commutators, and a degree-3 PSD state-optimality basis with two-site words
through separation 5.  The exact k=10 model did not fit this 31 GiB orb: even
after the axis quotient it was OOM-killed during the Mosek phase after about
4 minutes.  The practical runs below therefore use k=9.  Separation 5 suffices
at `N=10,20`; increasing the PSD basis to separation 10 recovers the requested
`1e-5` agreement at `N=30`, but does not replace the missing 10-site RDM at
`N=100`.

All rows use the moment-LMI formulation (`dualize=false`), Mosek 11 with 8
threads, and `1e-8` primal/dual feasibility and relative-gap tolerances.  The
reported value is JuMP's dual objective bound; where Mosek stopped with
`SLOW_PROGRESS`, both primal and dual statuses were `FEASIBLE_POINT` and the
primal objective agreed closely with the bound.

| N | RDM k | PSD range | paper SDP New / spin | this implementation / spin | absolute difference | wall time | peak RSS |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 10 | 9 | 5 | -0.4515446 | -0.4515432554 | 1.34e-6 | 140.2 s | 7.02 GiB |
| 20 | 9 | 5 | -0.4452196 | -0.4452226662 | 3.07e-6 | 264.9 s | 8.82 GiB |
| 30 | 9 | 10 | -0.4440668 | -0.4440728648 | 6.06e-6 | 303.0 s | 11.60 GiB |
| 100 | 9 | 10 | -0.4432378 | -0.4434115129 | 1.74e-4 | 1582.3 s | 28.07 GiB |

Thus the practical model meets the `1e-5` comparison target through `N=30`,
but not at `N=100`.  At `N=100` it nevertheless reduces the DMRG-relative gap
from 0.228% for the unstrengthened TI moment relaxation to 0.0411%, improving
on the paper's 0.0820% `SDP Old` gap.  The `N=100` solve took 26.4 minutes,
longer than extrapolation from the smaller runs, and reached a peak of 28.07
GiB.  Its full objective/bound were `-44.3411527203`/`-44.3411512925`, with 213
PSD blocks, largest side 126, and 42,797 total moments.

As a separate validity check, the exact `N=10` ground state satisfies all 14
retained PSD state-optimality blocks: the smallest eigenvalue was
`-2.04e-15`, and the largest linear-state-optimality residual was `3.86e-16`.
This distinguishes the known Mosek numerical sensitivity from a symbolic
sign or adjoint error in the new constraints.

## 2026-08-14: exact paper configuration (k=10) on a large-memory host

The exact QMBCertify configuration — `N=100`, 10-site RDM, `linear_psd` state
optimality with separation-5 two-site words, and the axis-permutation
quotient — was run on a 128-core, ~1 TiB host (`6xa800`) where the model fits.
Settings: moment-LMI form, Mosek with 32 threads, `1e-8` feasibility and
relative-gap tolerances, logging off.

```text
RESULT N=100 order=4 rdm=10 state=linear_psd state_range=5 axis_symmetry=true
objective=-44.3291247035 per_site=-0.4432912470
objective_bound=-44.3291243612 bound_per_site=-0.4432912436
status=SLOW_PROGRESS primal=FEASIBLE_POINT dual=FEASIBLE_POINT
wall=4693.3s max_block=252 n_blocks=214 unique_moments=31485
```

Peak RSS was 245.5 GiB (`Maximum resident set size: 257408656 kB`); wall time
78.2 minutes.

| configuration | bound / spin | vs paper `-0.4432378` | DMRG-relative gap |
|:--|--:|--:|--:|
| paper SDP New (k=10) | -0.4432378 | — | 0.0019% |
| this k=10 run | -0.4432912436 | 5.34e-5 | 0.0139% |
| practical k=9, range 10 | -0.4434115129 | 1.74e-4 | 0.0411% |
| paper SDP Old | — | — | 0.0820% |

The k=10 model tightens the bound by 3.2x over the practical k=9 model but
still misses the `1e-5` parity target by about 5x.  Mosek terminated with
`SLOW_PROGRESS` (both primal and dual statuses `FEASIBLE_POINT`, primal and
dual objectives agreeing to `3.4e-7` absolute), so the remaining gap may be
numerical stalling rather than model weakness.  A diagnostic rerun with 64
threads, tolerances relaxed to `1e-7`, and Mosek logging enabled is used to
distinguish the two; the harness exposes the tolerance through
`NCTS_MOSEK_TOL` (default `1e-8`).
