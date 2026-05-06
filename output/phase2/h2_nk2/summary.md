# Phase 2 H₂/Nk=2 A A* diagnostic summary

Generated: `2026-05-05T21:24:14`
Input dir: `output/phase2/h2_nk2`
Working RAM threshold for direct Cholesky: 111.441 GiB

## BPSDP solve check

| quantity | value |
|---|---:|
| termination | `NUMERICAL_ERROR` |
| BPSDP outer iterations | `0` |
| BPSDP inner iterations | `210` |
| BPSDP objective gap | `Inf` |

## Decision counts

| decision | count |
|---|---:|
| Direct Cholesky | 12 |
| Facial reduction needed | 24 |

## Instance-level decision

- `A_eq`: 18/18 blocks classified as `Facial reduction needed`.
- `A_full`: 12/18 blocks classified as `Direct Cholesky`; 6/18 still require facial reduction/presolve by the Phase-2 rule.
- BPSDP PSD-block solve status: `NUMERICAL_ERROR`.
- Verdict: row/nullspace elimination or facial-reduction-style presolve is a blocking workstream before treating H₄/Nk=2 as a plain BPSDP/Mazziotti-CG run. Direct Cholesky remains viable only for the cleaned `A_full` subblocks that survive that presolve.
- Causality note: the BPSDP `cg_failure` is consistent with the measured singular/equality-kernel structure, but this table alone does not prove facial reduction is the only cause; scaling and implementation checks remain follow-up work.

## Per-block decision table

| block | A | m | n | nnz(A) | nnz(AA*) | σ₁ | rank@1e-12 | ker | pred ker | excess | κ_range | gap | AMD nnz(L) | CG est | decision |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| twoD_K0_twoSz0 | A_full | 2312 | 1552 | 3604 | 8528 | 204.195 | 1552 | 760 | 602 | 158 | 204.195 | 7.610e+12 | 5687 | 128.607 | Direct Cholesky |
| twoD_K0_twoSz0 | A_eq | 232 | 1552 | 1524 | 1216 | 202.348 | 102 | 130 | 74 | 56 | 50.587 | 3.465e+14 | 658 | 64.012 | Facial reduction needed |
| twoD_K0_twoSz2 | A_full | 418 | 222 | 708 | 2146 | 80.000 | 222 | 196 | 116 | 80 | 80.000 | 1.474e+14 | 1324 | 80.498 | Facial reduction needed |
| twoD_K0_twoSz2 | A_eq | 118 | 222 | 408 | 562 | 78.000 | 52 | 66 | 38 | 28 | 39.000 | 2.346e+13 | 292 | 56.205 | Facial reduction needed |
| twoD_K0_twoSzm2 | A_full | 418 | 222 | 708 | 2146 | 80.000 | 222 | 196 | 116 | 80 | 80.000 | 1.474e+14 | 1324 | 80.498 | Facial reduction needed |
| twoD_K0_twoSzm2 | A_eq | 118 | 222 | 408 | 562 | 78.000 | 52 | 66 | 38 | 28 | 39.000 | 2.346e+13 | 292 | 56.205 | Facial reduction needed |
| twoD_K1_twoSz0 | A_full | 2312 | 1552 | 3600 | 8436 | 202.000 | 1552 | 760 | 602 | 158 | 202.000 | 8.035e+12 | 5641 | 127.914 | Direct Cholesky |
| twoD_K1_twoSz0 | A_eq | 232 | 1552 | 1520 | 1140 | 200.000 | 102 | 130 | 74 | 56 | 50.000 | 2.011e+13 | 626 | 63.640 | Facial reduction needed |
| twoD_K1_twoSz2 | A_full | 646 | 392 | 1264 | 3394 | 106.000 | 392 | 254 | 174 | 80 | 106.000 | 4.321e+13 | 2155 | 92.661 | Facial reduction needed |
| twoD_K1_twoSz2 | A_eq | 118 | 392 | 736 | 578 | 104.000 | 50 | 68 | 38 | 30 | 26.000 | 4.691e+13 | 318 | 45.891 | Facial reduction needed |
| twoD_K1_twoSzm2 | A_full | 646 | 392 | 1264 | 3394 | 106.000 | 392 | 254 | 174 | 80 | 106.000 | 4.321e+13 | 2155 | 92.661 | Facial reduction needed |
| twoD_K1_twoSzm2 | A_eq | 118 | 392 | 736 | 578 | 104.000 | 50 | 68 | 38 | 30 | 26.000 | 4.691e+13 | 318 | 45.891 | Facial reduction needed |
| twoG_K0_twoSz0 | A_full | 8755 | 3486 | 15347 | 60709 | 1.391e+04 | 3486 | 5269 | 4670 | 599 | 1.249e+04 | 8.730e+10 | 35536 | 1005.834 | Facial reduction needed |
| twoG_K0_twoSz0 | A_eq | 240 | 3486 | 7168 | 1792 | 1.390e+04 | 138 | 102 | 48 | 54 | 3.332e+04 | 2.777e+11 | 840 | 1642.791 | Facial reduction needed |
| twoG_K0_twoSz2 | A_full | 2314 | 1604 | 4148 | 10234 | 4024.963 | 1604 | 710 | 600 | 110 | 7.223e+04 | 2.165e+10 | 6599 | 2418.724 | Direct Cholesky |
| twoG_K0_twoSz2 | A_eq | 234 | 1604 | 1780 | 1258 | 4015.661 | 106 | 128 | 72 | 56 | 1619.739 | 1.099e+14 | 702 | 362.214 | Facial reduction needed |
| twoG_K0_twoSzm2 | A_full | 2314 | 1604 | 4148 | 10234 | 4024.963 | 1604 | 710 | 600 | 110 | 7.223e+04 | 2.310e+10 | 6597 | 2418.724 | Direct Cholesky |
| twoG_K0_twoSzm2 | A_eq | 234 | 1604 | 1780 | 1258 | 4015.661 | 106 | 128 | 72 | 56 | 1619.739 | 1.982e+13 | 678 | 362.214 | Facial reduction needed |
| twoG_K1_twoSz0 | A_full | 8492 | 3064 | 10272 | 34824 | 7879.974 | 3064 | 5428 | 4872 | 556 | 9.028e+04 | 4.726e+10 | 23464 | 2704.134 | Facial reduction needed |
| twoG_K1_twoSz0 | A_eq | 236 | 3064 | 2336 | 1224 | 7870.500 | 134 | 102 | 48 | 54 | 1967.625 | 7.114e+12 | 674 | 399.221 | Facial reduction needed |
| twoG_K1_twoSz2 | A_full | 2312 | 1604 | 4048 | 9884 | 4024.886 | 1604 | 708 | 600 | 108 | 7.222e+04 | 2.522e+10 | 6431 | 2418.701 | Direct Cholesky |
| twoG_K1_twoSz2 | A_eq | 232 | 1604 | 1680 | 1180 | 4015.584 | 104 | 128 | 72 | 56 | 1619.708 | 1.510e+12 | 648 | 362.210 | Facial reduction needed |
| twoG_K1_twoSzm2 | A_full | 2312 | 1604 | 4048 | 9884 | 4024.886 | 1604 | 708 | 600 | 108 | 7.222e+04 | 2.658e+10 | 6431 | 2418.701 | Direct Cholesky |
| twoG_K1_twoSzm2 | A_eq | 232 | 1604 | 1680 | 1180 | 4015.584 | 104 | 128 | 72 | 56 | 1619.708 | 6.595e+12 | 648 | 362.210 | Facial reduction needed |
| twoQ_K0_twoSz0 | A_full | 2320 | 1658 | 4708 | 16288 | 1.357e+04 | 1658 | 662 | 602 | 60 | 2.434e+05 | 5.070e+09 | 9616 | 4440.644 | Direct Cholesky |
| twoQ_K0_twoSz0 | A_eq | 240 | 1658 | 1956 | 1376 | 1.350e+04 | 110 | 130 | 74 | 56 | 8138.666 | 5.548e+11 | 754 | 811.931 | Facial reduction needed |
| twoQ_K0_twoSz2 | A_full | 428 | 276 | 1136 | 3452 | 9594.372 | 276 | 152 | 114 | 38 | 1.336e+05 | 1.700e+10 | 1988 | 3290.018 | Direct Cholesky |
| twoQ_K0_twoSz2 | A_eq | 128 | 276 | 632 | 704 | 9557.488 | 56 | 72 | 36 | 36 | 1233.805 | 2.660e+12 | 368 | 316.130 | Facial reduction needed |
| twoQ_K0_twoSzm2 | A_full | 428 | 276 | 1136 | 3452 | 9594.372 | 276 | 152 | 114 | 38 | 1.336e+05 | 1.582e+10 | 1988 | 3290.018 | Direct Cholesky |
| twoQ_K0_twoSzm2 | A_eq | 128 | 276 | 632 | 704 | 9557.488 | 56 | 72 | 36 | 36 | 1233.805 | 3.785e+12 | 368 | 316.130 | Facial reduction needed |
| twoQ_K1_twoSz0 | A_full | 2320 | 1658 | 4704 | 16256 | 1.357e+04 | 1658 | 662 | 602 | 60 | 2.435e+05 | 1.113e+10 | 9600 | 4440.703 | Direct Cholesky |
| twoQ_K1_twoSz0 | A_eq | 240 | 1658 | 1952 | 1312 | 1.350e+04 | 110 | 130 | 74 | 56 | 8851.726 | 8.669e+11 | 722 | 846.753 | Facial reduction needed |
| twoQ_K1_twoSz2 | A_full | 656 | 446 | 1824 | 6080 | 9624.723 | 446 | 210 | 172 | 38 | 1.727e+05 | 4.111e+10 | 3584 | 3740.241 | Direct Cholesky |
| twoQ_K1_twoSz2 | A_eq | 128 | 446 | 960 | 720 | 9575.639 | 56 | 72 | 36 | 36 | 1103.703 | 1.175e+13 | 400 | 298.998 | Facial reduction needed |
| twoQ_K1_twoSzm2 | A_full | 656 | 446 | 1824 | 6080 | 9624.723 | 446 | 210 | 172 | 38 | 1.727e+05 | 1.420e+10 | 3584 | 3740.241 | Direct Cholesky |
| twoQ_K1_twoSzm2 | A_eq | 128 | 446 | 960 | 720 | 9575.639 | 56 | 72 | 36 | 36 | 1103.703 | 3.304e+12 | 400 | 298.998 | Facial reduction needed |

## Moment-matrix rank at the BPSDP optimum

No optimum was obtained (`termination_status = NUMERICAL_ERROR`), so primal moment-rank validation is unavailable for this run.

| block | rank@1e-10 | rank@1e-12 | rank@1e-14 | min positive eig@1e-12 |
|---|---:|---:|---:|---:|

## Notes

- `pred ker` is the construction-level predicted redundancy: exact zero rows plus exact duplicate rows after local block restriction.
- `excess = empirical kernel@1e-12 - pred ker`; this is the number to chase for missed algebraic/facial-reduction structure.
- METIS fill was not used; AMD fill is CHOLMOD's ordering on `AA* + δI`.
- Histogram PNGs live under `plots/sigma_histogram_<block>_<operator>.png`.
