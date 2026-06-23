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
