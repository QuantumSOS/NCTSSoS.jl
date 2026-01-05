# Task: Fix State Polynomial Term Sparsity SOS Dualization

## Request

Run tests locally with MosekTools and debug any inconsistency in state polynomial term sparsity optimization. DO NOT modify expected test results.

## Background

The NCTSSoS.jl package implements sparse noncommutative polynomial optimization. State polynomial optimization with term sparsity is producing incorrect results:

- **Expected:** -4.0 (bilocal network test)
- **Actual:** -1.76M (completely wrong)

## Root Cause (Identified)

The SOS dualization in `_sos_dualize_state` (sos.jl:252-330) loses StateWord contributions when term sparsity creates multiple smaller blocks:

1. Constraint matrices are built per block with local (row, col) indexing
2. Coefficient extraction searches the full `total_basis`
3. The basis_idx returned points to full basis, not block basis
4. Multiple blocks containing the same StateWord at different positions lose contributions

## Decision: Direct Block Iteration

**Rationale:** Match NCTSSOS's on-the-fly loop approach instead of pre-computing symbolic matrices and then extracting.

**Approach:**
1. Store block basis with each constraint matrix in `StateMomentProblem`
2. Iterate directly over (row, col) pairs per block in `_sos_dualize_state`
3. Convert NCStateWord → StateWord and accumulate contributions
4. Remove `_get_state_Cαj` (no longer needed)

## Handoff Summary

- **Status:** Plan completed, ready for implementation
- **Branch:** ys/rev/rev/polynomial
- **Next step:** Implement the 4 code changes in plan.md and run tests
