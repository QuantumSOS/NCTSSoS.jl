# Test Tolerance Improvement Plan

## Overview

This document outlines the plan to improve test tolerances (`atol`) to match the main branch. **All 4 cases should first be tested with Mosek** to determine if the original tolerances are achievable.

---

## Problem Summary

After tightening tolerances, we identified 4 problem cases:

| Case | Test | Main Branch | Current | COSMO Result | Location |
|------|------|-------------|---------|--------------|----------|
| 1 | CS TS Example | `1e-4` | `1e-3` | Not tested with Mosek yet | `test/solvers/moment.jl`, `test/solvers/sos.jl` |
| 2 | Example 1 Dense | `1e-6` | `1e-5` | ~2e-6 (near threshold) | `test/solvers/moment.jl` line ~127 |
| 3 | Example 1 Sparse | `1e-7` | `1e-4` | `-0.0035406` vs expected `-0.0035512` | `test/solvers/moment.jl` line ~139 |
| 4 | SOS Term Sparsity | `1e-5` | `1e-4` | `0.9976161` vs expected `0.9975306` | `test/solvers/sos.jl` line ~282 |

---

## Step 1: Test ALL 4 Cases with Mosek

**Goal**: Determine if the original main branch tolerances are achievable with Mosek.

### Commands to Run on a800

```bash
# First, sync the current branch to a800
cd /home/yushengzhao/NCTSSoS.jl
git fetch origin
git checkout verify/comprehensive-functionality-check  # or the current branch
git pull

# Run moment.jl tests (covers Cases 1, 2, 3)
LOCAL_TESTING=true julia --project -e '
include("test/setup.jl")
include("test/solvers/moment.jl")
'

# Run sos.jl tests (covers Cases 1, 4)
LOCAL_TESTING=true julia --project -e '
include("test/setup.jl")
include("test/solvers/sos.jl")
'
```

### Expected Outcomes

For each case, record the **actual objective value** with Mosek:

| Case | Test | Expected Value | Mosek Result | Error | Can Use Main Branch Tolerance? |
|------|------|----------------|--------------|-------|-------------------------------|
| 1 | CS TS Example | `3.011288` | `3.0112938825722892` | `5.88e-6` | ✅ YES (`1e-4`) |
| 2 | Example 1 Dense | `0.0` | `-3.23e-8` | `3.23e-8` | ✅ YES (`1e-6`) |
| 3 | Example 1 Sparse | `-0.0035512` | `-0.0035512666506813` | `6.67e-8` | ✅ YES (`1e-7`) |
| 4 | SOS Term Sparsity | `0.9975306427277915` | `0.997530500521655` | `1.42e-7` | ✅ YES (`1e-5`) |

**Tested on**: 2024-12-29, a800 server with Mosek 10.x

---

## Step 2: Decision Tree Based on Mosek Results

### ✅ RESULT: Mosek achieves ALL main branch tolerances

**All 4 Cases**: Pure COSMO solver precision issue.
- **Action**: Keep current looser tolerances for COSMO (CI)
- **Reason**: COSMO is a free, open-source solver with lower precision than commercial Mosek
- **No algorithm bugs detected** - the implementation is correct

The precision gap is expected:
- Mosek achieves `1e-6` to `1e-8` precision
- COSMO achieves `1e-3` to `1e-5` precision

---

## Step 3: ~~Investigate Term Sparsity~~ (NOT NEEDED)

**SKIPPED** - Mosek achieved all main branch tolerances, confirming no algorithm bugs.

---

## Step 4: Update Tolerances ✅ COMPLETED

Based on Mosek test results, using **Option A**: Keep looser tolerances for COSMO (CI).

### Final Tolerance Strategy

| Test | COSMO (CI) | Mosek (LOCAL_TESTING) | Reason |
|------|------------|----------------------|--------|
| CS TS Example | `1e-3` | `1e-4` achievable | Large sparse problem |
| Example 1 Dense | `1e-5` | `1e-6` achievable | Good precision |
| Example 1 Sparse | `1e-4` | `1e-7` achievable | Term sparsity relaxes bound |
| SOS Term Sparsity | `1e-4` | `1e-5` achievable | Term sparsity relaxes bound |

Tests are wrapped in `if LOCAL_TESTING` blocks where Mosek is required for precision.

---

## Checklist

- [x] Sync current branch to a800 (2024-12-29)
- [x] Run `test/solvers/moment.jl` with Mosek on a800 (ALL PASS)
- [x] Run `test/solvers/sos.jl` with Mosek on a800 (ALL PASS)
- [x] Record actual objective values for all 4 cases (see table above)
- [x] Determine if main branch tolerances are achievable (YES, all 4 cases)
- [x] ~~If NOT achievable with Mosek: investigate term sparsity implementation~~ (NOT NEEDED)
- [x] Update tolerances accordingly (keep current looser tolerances for COSMO)
- [x] Document findings (this document)

---

## Summary Table: Main Branch vs Current Tolerances

| File | Test | Main Branch | Current (COSMO) | Mosek Error | Status |
|------|------|-------------|-----------------|-------------|--------|
| `moment.jl` | CHSH Inequality | `1e-6` | `1e-6` | N/A | ✅ OK |
| `moment.jl` | CS TS Example | `1e-4` | `1e-3` | `5.88e-6` | ✅ COSMO limitation |
| `moment.jl` | Example 1 Dense | `1e-6` | `1e-5` | `3.23e-8` | ✅ COSMO limitation |
| `moment.jl` | Example 1 Sparse | `1e-7` | `1e-4` | `6.67e-8` | ✅ COSMO limitation |
| `sos.jl` | I_3322 | `1e-6` | `1e-6` | N/A | ✅ OK |
| `sos.jl` | CS TS Example | `1e-4` | `1e-3` | `5.88e-6` | ✅ COSMO limitation |
| `sos.jl` | SOS Term Sparsity | `1e-5` | `1e-4` | `1.42e-7` | ✅ COSMO limitation |

## Conclusion

**No algorithm bugs detected.** The tolerance differences between main branch and current implementation are entirely due to solver precision:

- **Mosek** (commercial): Achieves `1e-6` to `1e-8` precision
- **COSMO** (open-source): Achieves `1e-3` to `1e-5` precision

The current test tolerances are appropriate for CI (using COSMO/Clarabel) while still validating correctness. Users running with Mosek in LOCAL_TESTING mode will see tighter bounds.
