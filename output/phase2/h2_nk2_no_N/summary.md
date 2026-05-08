# Phase 2 H₂/Nk=2 A A* diagnostic summary

Generated: `2026-05-07T22:34:02`
Input dir: `output/phase2/h2_nk2_no_N`
Working RAM threshold for direct Cholesky: 111.441 GiB

## BPSDP solve check

| quantity | value |
|---|---:|
| termination | `ITERATION_LIMIT` |
| BPSDP outer iterations | `5000` |
| BPSDP inner iterations | `88278` |
| BPSDP objective gap | `0.0005236154039920393` |

## Decision counts

| decision | count |
|---|---:|
| Direct Cholesky | 16 |
| Facial reduction needed | 20 |

## Instance-level decision

- `A_eq`: 18/18 blocks classified as `Facial reduction needed`.
- `A_full`: 16/18 blocks classified as `Direct Cholesky`; 2/18 still require facial reduction/presolve by the Phase-2 rule.
- BPSDP PSD-block solve status: `ITERATION_LIMIT`.
- Verdict: row/nullspace elimination or facial-reduction-style presolve is a blocking workstream before treating H₄/Nk=2 as a plain BPSDP/Mazziotti-CG run. Direct Cholesky remains viable only for the cleaned `A_full` subblocks that survive that presolve.
- Causality note: no immediate BPSDP `cg_failure` was observed; any remaining singularity is residual algebraic/scaling structure rather than orphan-variable fallout.

## Per-block decision table

| block | A | m | n | nnz(A) | nnz(AA*) | σ₁ | rank@1e-12 | ker | pred ker | excess | κ_range | gap | AMD nnz(L) | CG est | decision |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| twoD_K0_twoSz0 | A_full | 2088 | 1552 | 2324 | 3888 | 196.204 | 1552 | 536 | 530 | 6 | 196.204 | 5.446e+12 | 2998 | 126.066 | Direct Cholesky |
| twoD_K0_twoSz0 | A_eq | 8 | 1552 | 244 | 32 | 194.365 | 4 | 4 | 2 | 2 | 10.922 | 3.985e+30 | 20 | 29.744 | Facial reduction needed |
| twoD_K0_twoSz2 | A_full | 306 | 222 | 372 | 690 | 74.000 | 222 | 84 | 80 | 4 | 74.000 | 4.215e+13 | 505 | 77.421 | Direct Cholesky |
| twoD_K0_twoSz2 | A_eq | 6 | 222 | 72 | 18 | 72.000 | 2 | 4 | 2 | 2 | 1.000 | 1.013e+15 | 12 | 9.000 | Facial reduction needed |
| twoD_K0_twoSzm2 | A_full | 306 | 222 | 372 | 690 | 74.000 | 222 | 84 | 80 | 4 | 74.000 | 4.215e+13 | 505 | 77.421 | Direct Cholesky |
| twoD_K0_twoSzm2 | A_eq | 6 | 222 | 72 | 18 | 72.000 | 2 | 4 | 2 | 2 | 1.000 | 1.013e+15 | 12 | 9.000 | Facial reduction needed |
| twoD_K1_twoSz0 | A_full | 2088 | 1552 | 2320 | 3860 | 194.000 | 1552 | 536 | 530 | 6 | 194.000 | 3.431e+12 | 2981 | 125.355 | Direct Cholesky |
| twoD_K1_twoSz0 | A_eq | 8 | 1552 | 240 | 20 | 192.000 | 4 | 4 | 2 | 2 | 12.000 | 9.382e+13 | 14 | 31.177 | Facial reduction needed |
| twoD_K1_twoSz2 | A_full | 534 | 392 | 624 | 1106 | 98.000 | 392 | 142 | 138 | 4 | 98.000 | 1.273e+13 | 827 | 89.095 | Direct Cholesky |
| twoD_K1_twoSz2 | A_eq | 6 | 392 | 96 | 18 | 96.000 | 2 | 4 | 2 | 2 | 1.000 | 1.126e+15 | 12 | 9.000 | Facial reduction needed |
| twoD_K1_twoSzm2 | A_full | 534 | 392 | 624 | 1106 | 98.000 | 392 | 142 | 138 | 4 | 98.000 | 1.273e+13 | 827 | 89.095 | Direct Cholesky |
| twoD_K1_twoSzm2 | A_eq | 6 | 392 | 96 | 18 | 96.000 | 2 | 4 | 2 | 2 | 1.000 | 1.126e+15 | 12 | 9.000 | Facial reduction needed |
| twoG_K0_twoSz0 | A_full | 8529 | 3486 | 9041 | 28129 | 1.386e+04 | 3486 | 5043 | 4622 | 421 | 1.259e+04 | 7.123e+10 | 18951 | 1009.935 | Facial reduction needed |
| twoG_K0_twoSz0 | A_eq | 14 | 3486 | 862 | 98 | 1.386e+04 | 12 | 2 | 0 | 2 | 3.323e+04 | 1.039e+12 | 56 | 1640.508 | Facial reduction needed |
| twoG_K0_twoSz2 | A_full | 2088 | 1604 | 2692 | 5360 | 4006.010 | 1554 | 534 | 528 | 6 | 4006.010 | 9.964e+11 | 3739 | 569.637 | Direct Cholesky |
| twoG_K0_twoSz2 | A_eq | 8 | 1604 | 324 | 32 | 3996.712 | 6 | 2 | 0 | 2 | 134.208 | 2.055e+14 | 20 | 104.263 | Facial reduction needed |
| twoG_K0_twoSzm2 | A_full | 2088 | 1604 | 2692 | 5360 | 4006.010 | 1554 | 534 | 528 | 6 | 4006.010 | 9.964e+11 | 3739 | 569.637 | Direct Cholesky |
| twoG_K0_twoSzm2 | A_eq | 8 | 1604 | 324 | 32 | 3996.712 | 6 | 2 | 0 | 2 | 134.208 | 2.055e+14 | 20 | 104.263 | Facial reduction needed |
| twoG_K1_twoSz0 | A_full | 8266 | 3064 | 8256 | 26334 | 7852.974 | 3040 | 5226 | 4824 | 402 | 1.028e+04 | 3.445e+11 | 18608 | 912.499 | Facial reduction needed |
| twoG_K1_twoSz0 | A_eq | 10 | 3064 | 320 | 30 | 7843.507 | 8 | 2 | 0 | 2 | 490.219 | 3.325e+13 | 20 | 199.268 | Facial reduction needed |
| twoG_K1_twoSz2 | A_full | 2086 | 1604 | 2592 | 5074 | 4005.933 | 1554 | 532 | 528 | 4 | 4005.933 | 6.434e+11 | 3580 | 569.632 | Direct Cholesky |
| twoG_K1_twoSz2 | A_eq | 6 | 1604 | 224 | 18 | 3996.636 | 4 | 2 | 0 | 2 | 92.164 | 2.250e+14 | 12 | 86.402 | Facial reduction needed |
| twoG_K1_twoSzm2 | A_full | 2086 | 1604 | 2592 | 5074 | 4005.933 | 1554 | 532 | 528 | 4 | 4005.933 | 6.434e+11 | 3580 | 569.632 | Direct Cholesky |
| twoG_K1_twoSzm2 | A_eq | 6 | 1604 | 224 | 18 | 3996.636 | 4 | 2 | 0 | 2 | 92.164 | 2.250e+14 | 12 | 86.402 | Facial reduction needed |
| twoQ_K0_twoSz0 | A_full | 2094 | 1658 | 3074 | 11362 | 1.353e+04 | 1556 | 538 | 530 | 8 | 1.353e+04 | 1.527e+11 | 6728 | 1046.824 | Direct Cholesky |
| twoQ_K0_twoSz0 | A_eq | 14 | 1658 | 322 | 98 | 1.346e+04 | 8 | 6 | 2 | 4 | 8607.038 | 1.443e+12 | 56 | 834.967 | Facial reduction needed |
| twoQ_K0_twoSz2 | A_full | 314 | 276 | 622 | 1862 | 9566.275 | 226 | 88 | 78 | 10 | 9566.275 | 5.237e+11 | 1052 | 880.266 | Direct Cholesky |
| twoQ_K0_twoSz2 | A_eq | 14 | 276 | 118 | 98 | 9529.393 | 6 | 8 | 0 | 8 | 1238.254 | 7.247e+12 | 56 | 316.699 | Facial reduction needed |
| twoQ_K0_twoSzm2 | A_full | 314 | 276 | 622 | 1862 | 9566.275 | 226 | 88 | 78 | 10 | 9566.275 | 3.826e+11 | 1052 | 880.266 | Direct Cholesky |
| twoQ_K0_twoSzm2 | A_eq | 14 | 276 | 118 | 98 | 9529.393 | 6 | 8 | 0 | 8 | 1238.254 | 5.468e+12 | 56 | 316.699 | Facial reduction needed |
| twoQ_K1_twoSz0 | A_full | 2094 | 1658 | 3070 | 11394 | 1.353e+04 | 1556 | 538 | 530 | 8 | 1.353e+04 | 5.979e+10 | 6744 | 1046.838 | Direct Cholesky |
| twoQ_K1_twoSz0 | A_eq | 14 | 1658 | 318 | 98 | 1.346e+04 | 8 | 6 | 2 | 4 | 9358.889 | 4.076e+12 | 56 | 870.672 | Facial reduction needed |
| twoQ_K1_twoSz2 | A_full | 542 | 446 | 1006 | 3634 | 9596.419 | 396 | 146 | 136 | 10 | 9596.419 | 1.344e+12 | 2040 | 881.652 | Direct Cholesky |
| twoQ_K1_twoSz2 | A_eq | 14 | 446 | 142 | 98 | 9547.336 | 6 | 8 | 0 | 8 | 1120.417 | 8.063e+15 | 56 | 301.254 | Facial reduction needed |
| twoQ_K1_twoSzm2 | A_full | 542 | 446 | 1006 | 3634 | 9596.419 | 396 | 146 | 136 | 10 | 9596.419 | 7.752e+11 | 2040 | 881.652 | Direct Cholesky |
| twoQ_K1_twoSzm2 | A_eq | 14 | 446 | 142 | 98 | 9547.336 | 6 | 8 | 0 | 8 | 1120.417 | 5.217e+15 | 56 | 301.254 | Facial reduction needed |

## Moment-matrix rank at the BPSDP optimum

No optimum was obtained (`termination_status = ITERATION_LIMIT`), so primal moment-rank validation is unavailable for this run.

| block | rank@1e-10 | rank@1e-12 | rank@1e-14 | min positive eig@1e-12 |
|---|---:|---:|---:|---:|

## Notes

- `pred ker` is the construction-level predicted redundancy: exact zero rows plus exact duplicate rows after local block restriction.
- `excess = empirical kernel@1e-12 - pred ker`; this is the number to chase for missed algebraic/facial-reduction structure.
- METIS fill was not used; AMD fill is CHOLMOD's ordering on `AA* + δI`.
- Histogram PNGs live under `plots/sigma_histogram_<block>_<operator>.png`.
