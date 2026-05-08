# H₂/Nk=2 no explicit particle-number constraint comparison

Generated from `output/phase2/pivot_coverage_audit.md`, `output/phase2/pivot_coverage_audit_no_N.md`, `output/phase2/h2_nk2/summary.md`, and `output/phase2/h2_nk2_no_N/summary.md`.

## Pivot coverage

| run | canonical moments | pivots | orphans | orphan % | zero constraints |
|---|---:|---:|---:|---:|---:|
| with explicit `N̂ - N` | 9,393 | 2,913 | 6,480 | 68.9875 | 983 |
| without explicit `N̂ - N` | 2,913 | 2,913 | 0 | 0.0000 | 6 |

The orphan source is gone. No attribution probe was needed.

## BPSDP bridge / solve metadata

| quantity | with `N̂ - N` | without `N̂ - N` |
|---|---:|---:|
| termination | `NUMERICAL_ERROR` | `ITERATION_LIMIT` |
| raw BPSDP reason | `cg_failure` | `max_iter` |
| BPSDP blocks | 221 | 18 |
| scalar 1×1 PSD blocks | 0 with `:aux_psd_free` (406,706 historically with `:free_variables`) | 0 |
| largest block dim | 65 | 65 |
| `A` shape (`dual_rows × primal_dim`) | 17,475 × 238,380 | 30,394 × 18,113 |
| primal dim | 238,380 | 18,113 |
| dual rows | 17,475 | 30,394 |
| A nnz | 170,704 | 123,415 |
| outer iterations | 0 | 5,000 |
| inner CG iterations | 210 | 88,278 |
| objective gap | `Inf` | 5.236154039920393e-4 |

Dropping the redundant constraint changes the failure mode: BPSDP no longer dies immediately in CG. It runs normally until the configured 5,000-iteration cap.

## A*A' diagnostic changes

| operator | facial-reduction blocks before | facial-reduction blocks after | excess-kernel sum before | excess-kernel sum after | max excess before | max excess after |
|---|---:|---:|---:|---:|---:|---:|
| `A_eq` | 18/18 | 18/18 | 816 | 64 | 56 | 8 |
| `A_full` | 6/18 | 2/18 | 2,499 | 927 | 599 | 421 |

`A_eq` is still singular in every block, but the huge particle-number localizing-row kernel is gone. The surviving rank deficiency is now genuine trace/spin/PQG algebra, not orphan-variable fallout.

## Per-block excess-kernel comparison

### `A_eq`

| block | excess before | excess after | decision before | decision after |
|---|---:|---:|---|---|
| `twoD_K0_twoSz0` | 56 | 2 | Facial reduction needed | Facial reduction needed |
| `twoD_K0_twoSz2` | 28 | 2 | Facial reduction needed | Facial reduction needed |
| `twoD_K0_twoSzm2` | 28 | 2 | Facial reduction needed | Facial reduction needed |
| `twoD_K1_twoSz0` | 56 | 2 | Facial reduction needed | Facial reduction needed |
| `twoD_K1_twoSz2` | 30 | 2 | Facial reduction needed | Facial reduction needed |
| `twoD_K1_twoSzm2` | 30 | 2 | Facial reduction needed | Facial reduction needed |
| `twoG_K0_twoSz0` | 54 | 2 | Facial reduction needed | Facial reduction needed |
| `twoG_K0_twoSz2` | 56 | 2 | Facial reduction needed | Facial reduction needed |
| `twoG_K0_twoSzm2` | 56 | 2 | Facial reduction needed | Facial reduction needed |
| `twoG_K1_twoSz0` | 54 | 2 | Facial reduction needed | Facial reduction needed |
| `twoG_K1_twoSz2` | 56 | 2 | Facial reduction needed | Facial reduction needed |
| `twoG_K1_twoSzm2` | 56 | 2 | Facial reduction needed | Facial reduction needed |
| `twoQ_K0_twoSz0` | 56 | 4 | Facial reduction needed | Facial reduction needed |
| `twoQ_K0_twoSz2` | 36 | 8 | Facial reduction needed | Facial reduction needed |
| `twoQ_K0_twoSzm2` | 36 | 8 | Facial reduction needed | Facial reduction needed |
| `twoQ_K1_twoSz0` | 56 | 4 | Facial reduction needed | Facial reduction needed |
| `twoQ_K1_twoSz2` | 36 | 8 | Facial reduction needed | Facial reduction needed |
| `twoQ_K1_twoSzm2` | 36 | 8 | Facial reduction needed | Facial reduction needed |

### `A_full`

| block | excess before | excess after | decision before | decision after |
|---|---:|---:|---|---|
| `twoD_K0_twoSz0` | 158 | 6 | Direct Cholesky | Direct Cholesky |
| `twoD_K0_twoSz2` | 80 | 4 | Facial reduction needed | Direct Cholesky |
| `twoD_K0_twoSzm2` | 80 | 4 | Facial reduction needed | Direct Cholesky |
| `twoD_K1_twoSz0` | 158 | 6 | Direct Cholesky | Direct Cholesky |
| `twoD_K1_twoSz2` | 80 | 4 | Facial reduction needed | Direct Cholesky |
| `twoD_K1_twoSzm2` | 80 | 4 | Facial reduction needed | Direct Cholesky |
| `twoG_K0_twoSz0` | 599 | 421 | Facial reduction needed | Facial reduction needed |
| `twoG_K0_twoSz2` | 110 | 6 | Direct Cholesky | Direct Cholesky |
| `twoG_K0_twoSzm2` | 110 | 6 | Direct Cholesky | Direct Cholesky |
| `twoG_K1_twoSz0` | 556 | 402 | Facial reduction needed | Facial reduction needed |
| `twoG_K1_twoSz2` | 108 | 4 | Direct Cholesky | Direct Cholesky |
| `twoG_K1_twoSzm2` | 108 | 4 | Direct Cholesky | Direct Cholesky |
| `twoQ_K0_twoSz0` | 60 | 8 | Direct Cholesky | Direct Cholesky |
| `twoQ_K0_twoSz2` | 38 | 10 | Direct Cholesky | Direct Cholesky |
| `twoQ_K0_twoSzm2` | 38 | 10 | Direct Cholesky | Direct Cholesky |
| `twoQ_K1_twoSz0` | 60 | 8 | Direct Cholesky | Direct Cholesky |
| `twoQ_K1_twoSz2` | 38 | 10 | Direct Cholesky | Direct Cholesky |
| `twoQ_K1_twoSzm2` | 38 | 10 | Direct Cholesky | Direct Cholesky |

## Bridge smoke

`NCTSSOS_RUN_H2_NK2_BRIDGE_SMOKE=1 julia --startup-file=no --project=/tmp/h2_noN_env.LQrrIs test/relaxations/h2_nk2_bridge_smoke.jl` passed on HAI: 10/10 tests in 40.5 s.
