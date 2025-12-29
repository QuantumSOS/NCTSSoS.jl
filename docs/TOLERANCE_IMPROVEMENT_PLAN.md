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
| 1 | CS TS Example | `3.011288` | ___ | ___ | `1e-4`? |
| 2 | Example 1 Dense | `0.0` | ___ | ___ | `1e-6`? |
| 3 | Example 1 Sparse | `-0.0035512` | ___ | ___ | `1e-7`? |
| 4 | SOS Term Sparsity | `0.9975306427277915` | ___ | ___ | `1e-5`? |

---

## Step 2: Decision Tree Based on Mosek Results

### If Mosek achieves main branch tolerances:

**Cases 1 & 2**: These are pure COSMO solver precision issues.
- Keep current looser tolerances for COSMO (CI)
- Optionally add `if LOCAL_TESTING` blocks with tighter tolerances for Mosek

**Cases 3 & 4**: If Mosek achieves `1e-7` and `1e-5` respectively:
- COSMO is the bottleneck, not the algorithm
- Same solution as Cases 1 & 2

### If Mosek STILL fails to achieve main branch tolerances:

This indicates a **potential bug** in the term sparsity implementation. Proceed to Step 3.

---

## Step 3: Investigate Term Sparsity (Only if Mosek fails)

### Key Files
- `src/optimization/sparsity.jl` - term sparsity implementation
- `src/optimization/moment.jl` - moment matrix construction

### Likely Suspects

#### Suspect A: Canonicalization Mismatch in `init_activated_supp`

```julia
function init_activated_supp(partial_obj, cons, mom_mtx_bases)
    return sorted_union(
        symmetric_canon.(monomials(partial_obj)),  # ← only objective is canonicalized
        mapreduce(monomials, vcat, cons; init=M[]),  # ← NOT canonicalized
        [neat_dot(b, b) for b in mom_mtx_bases]      # ← NOT canonicalized
    )
end
```

**Fix**: Apply `symmetric_canon` consistently to all components.

#### Suspect B: Missing Iteration Loop

Current implementation does only **one round** of term sparsity. Reference may iterate to fixed point:

```julia
activated = init_activated_supp(...)
while true
    ts = iterate_term_sparse_supp(activated, ...)
    new_activated = union(activated, ts.term_sparse_graph_supp)
    if new_activated == activated
        break
    end
    activated = new_activated
end
```

### Validation: Compare with Reference NCTSSOS

```bash
ssh a800 "cd /home/yushengzhao/NCTSSOS && julia --project -e 'include(\"test/moment_solver.jl\")'"
```

---

## Step 4: Update Tolerances

Based on Mosek test results, update the test files:

### Option A: Mosek achieves all tolerances (solver precision issue)

Keep looser tolerances for COSMO, add comments:

```julia
# COSMO achieves ~1e-5, Mosek achieves 1e-7
@test isapprox(result.objective, expected, atol=1e-5)
```

### Option B: Mosek fails some tolerances (algorithm issue)

1. Fix the bug (canonicalization or iteration)
2. Re-run tests
3. Tighten tolerances to match main branch

---

## Checklist

- [ ] Sync current branch to a800
- [ ] Run `test/solvers/moment.jl` with Mosek on a800
- [ ] Run `test/solvers/sos.jl` with Mosek on a800
- [ ] Record actual objective values for all 4 cases
- [ ] Determine if main branch tolerances are achievable
- [ ] If NOT achievable with Mosek: investigate term sparsity implementation
- [ ] Update tolerances accordingly
- [ ] Document findings

---

## Summary Table: Main Branch vs Current Tolerances

| File | Test | Main Branch | Current | Action Needed |
|------|------|-------------|---------|---------------|
| `moment.jl` | CHSH Inequality | `1e-6` | `1e-6` | None |
| `moment.jl` | CS TS Example | `1e-4` | `1e-3` | Test with Mosek |
| `moment.jl` | Example 1 Dense | `1e-6` | `1e-5` | Test with Mosek |
| `moment.jl` | Example 1 Sparse | `1e-7` | `1e-4` | Test with Mosek |
| `sos.jl` | I_3322 | `1e-6` | `1e-6` | None |
| `sos.jl` | CS TS Example | `1e-4` | `1e-3` | Test with Mosek |
| `sos.jl` | SOS Term Sparsity | `1e-5` | `1e-4` | Test with Mosek |
