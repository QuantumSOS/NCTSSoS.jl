# Bell Inequalities Test Results

**Test file:** `test/problems/bell_inequalities/bell_inequalities.jl`
**Solver:** MosekTools (`--local` flag)
**Date:** 2026-01-03
**Total runtime:** 458 minutes (~7.6 hours)

## Summary

| Status | Count |
|--------|-------|
| **Passed** | 86 |
| **Failed** | 2 |
| **Total** | 88 |

## Failing Tests

### 1. A77 (Borderline Failure)

| Field | Value |
|-------|-------|
| **Test ID** | A77 |
| **Actual** | -0.6655565218932911 |
| **Expected** | -0.6655582 |
| **Difference** | ~1.7e-6 |
| **Tolerance** | 1e-6 |
| **Runtime** | 3m37.6s |

**Analysis:** This is a borderline failure - the difference is only slightly outside the 1e-6 tolerance. Likely a numerical precision issue rather than a fundamental algorithmic problem.

---

### 2. A10 (Significant Failure)

| Field | Value |
|-------|-------|
| **Test ID** | A10 |
| **Actual** | -4.207106780806911 |
| **Expected** | -0.4158004 |
| **Difference** | ~3.79 (order of magnitude) |
| **Tolerance** | 1e-6 |
| **Runtime** | 41.3s |
| **Relaxation Order** | 3 (from test config: `d[57] = 3`) |

**Analysis:** This is a significant failure - the computed value is approximately 10x larger in magnitude than expected. This suggests a potential issue with:
- The relaxation order configuration
- The polynomial parsing for this specific instance
- A structural issue with this particular Bell inequality configuration

**Instance details (from test file):**
- `ms[9] = 4` (4 X variables)
- `ns[9] = 5` (5 Y variables)
- Equation: `-1*X[3]-2*X[4]-2*Y[1]-1*Y[2]-1*Y[3]+1*X[1]*Y[1]+1*X[2]*Y[1]+1*X[3]*Y[1]-1*X[1]*Y[2]+1*X[2]*Y[2]+1*X[3]*Y[2]+1*X[4]*Y[2]+1*X[1]*Y[3]+1*X[4]*Y[3]-1*X[2]*Y[4]+1*X[3]*Y[4]-1*X[2]*Y[5]+1*X[4]*Y[5]`

---

## All Test Results

| Test | Status | Runtime | Notes |
|------|--------|---------|-------|
| A2 | Pass | 11.7s | |
| A3 | Pass | 0.3s | |
| A4 | Pass | 0.9s | |
| A5 | Pass | 6.0s | |
| A6 | Pass | 6.9s | |
| A7 | Pass | 5.2s | |
| A8 | Pass | 0.0s | |
| A9 | Pass | 36.4s | |
| **A10** | **FAIL** | 41.3s | obj=-4.21, exp=-0.42 |
| A11 | Pass | 38.8s | |
| A12 | Pass | 0.1s | |
| A13 | Pass | 35.6s | |
| A14 | Pass | 2m02.3s | |
| A15 | Pass | 38.9s | |
| A16 | Pass | 39.3s | |
| A17 | Pass | 0.1s | |
| A18 | Pass | 0.1s | |
| A19 | Pass | 35.7s | |
| A20 | Pass | 41.5s | |
| A21 | Pass | 17m27.8s | |
| A22 | Pass | 38.5s | |
| A23 | Pass | 0.3s | |
| A24 | Pass | 0.2s | |
| A25 | Pass | 0.2s | |
| A26 | Pass | 0.2s | |
| A27 | Pass | 0.3s | |
| A28 | Pass | 0.3s | |
| A29 | Pass | 0.3s | |
| A30 | Pass | 0.2s | |
| A31 | Pass | 0.2s | |
| A32 | Pass | 3m59.4s | |
| A33 | Pass | 3m05.9s | |
| A34 | Pass | 3m29.6s | |
| A35 | Pass | 3m47.4s | |
| A36 | Pass | 3m20.0s | |
| A37 | Pass | 3m11.6s | |
| A38 | Pass | 3m15.2s | |
| A39 | Pass | 0.2s | |
| A40 | Pass | 0.2s | |
| A41 | Pass | 3m15.3s | |
| A42 | Pass | 0.2s | |
| A43 | Pass | 0.2s | |
| A44 | Pass | 3m50.4s | |
| A45 | Pass | 3m08.7s | |
| A46 | Pass | 20m00.1s | |
| A47 | Pass | 26m09.1s | |
| A48 | Pass | 17m51.8s | |
| A49 | Pass | 3m10.1s | |
| A50 | Pass | 3m11.5s | |
| A51 | Pass | 4m02.5s | |
| A52 | Pass | 3m30.3s | |
| A53 | Pass | 3m50.4s | |
| A54 | Pass | 3m50.4s | |
| A55 | Pass | 0.2s | |
| A56 | Pass | 6m25.2s | |
| A57 | Pass | 3m48.7s | |
| A58 | Pass | 3m47.8s | |
| A59 | Pass | 3m49.5s | |
| A60 | Pass | 7m34.7s | |
| A61 | Pass | 3m33.1s | |
| A62 | Pass | 42m57.5s | |
| A63 | Pass | 36m37.6s | |
| A64 | Pass | 20m02.9s | |
| A65 | Pass | 4m21.1s | |
| A66 | Pass | 3m09.9s | |
| A67 | Pass | 9m16.9s | |
| A68 | Pass | 40m19.5s | |
| A69 | Pass | 7m00.2s | |
| A70 | Pass | 0.2s | |
| A71 | Pass | 4m13.0s | |
| A72 | Pass | 3m28.4s | |
| A73 | Pass | 3m09.0s | |
| A74 | Pass | 3m47.3s | |
| A75 | Pass | 3m48.1s | |
| A76 | Pass | 2m49.9s | |
| **A77** | **FAIL** | 3m37.6s | obj=-0.6655565, exp=-0.6655582 |
| A78 | Pass | 0.2s | |
| A79 | Pass | 3m29.7s | |
| A80 | Pass | 9m42.9s | |
| A81 | Pass | 2m51.6s | |
| A82 | Pass | 35m31.7s | |
| A83 | Pass | 3m26.3s | |
| A84 | Pass | 20m59.8s | |
| A85 | Pass | 3m06.6s | |
| A86 | Pass | 3m08.2s | |
| A87 | Pass | 7m42.4s | |
| A88 | Pass | 0.5s | |
| A89 | Pass | 4m31.7s | |
