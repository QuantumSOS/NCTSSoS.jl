# Task: Debug Unphysical Moment Matrices

## Problem Statement
The user reports that they are currently generating unphysical moment matrices. This is a critical issue that affects the correctness of the GNS reconstruction and the entire optimization pipeline.

## What "Unphysical" Likely Means
A moment matrix is unphysical if:
1. **Not Hermitian**: H ≠ H† (for complex moment matrices)
2. **Not PSD**: H has negative eigenvalues (should be positive semidefinite)
3. **Incorrect structure**: Doesn't satisfy the Hankel structure (M[α,β] = μ[αβ])
4. **Violates constraints**: Doesn't satisfy the problem constraints (e.g., unipotency)

## Context
- We just validated that the complex-to-real SDP reformulation is mathematically correct
- The codebase uses the efficient dual formulation: H_R = (X1 + X2), H_I = (X3 - X3^T)
- The mermin_square_gns.jl example includes a division by 2: H_R = (X1 + X2)/2
- The issue might be in how the moment matrix is extracted from the solver

## Investigation Areas

### 1. Primal vs Dual Confusion
- Are we extracting from the primal or dual solution?
- The realified problem is the DUAL of the complex moment problem
- Should we use `value()` or `dual()`?

### 2. The Division by 2 Factor
From mermin_square_gns.jl:167-168:
```julia
H_R = (X1 + X2) / 2
H_I = (X3 - X3') / 2
```

Why divide by 2? Is this:
- Related to the dual formulation scaling?
- A correction factor from the paper?
- An error?

### 3. Sign Conventions
- Check if signs are correct in the reconstruction
- Verify X3 - X3^T vs X3^T - X3

### 4. Constraint Extraction
Current approach (mermin_square_gns.jl:148-152):
```julia
cons = all_constraints(model, include_variable_in_set_constraints=true)
X = value.(dual(cons[1]))
```

Questions:
- Is cons[1] always the right constraint?
- Should we use `value()` or `dual()`?
- Is the dualize=true flag causing issues?

### 5. Numerical Issues
- Is the solver converging properly?
- Are there numerical precision issues?
- Check eigenvalues of the reconstructed H

## Key Questions to Answer
1. What does "unphysical" mean specifically in this case? (Need examples)
2. Is H Hermitian?
3. Is H positive semidefinite? (Check eigenvalues)
4. Where does the /2 factor come from and is it correct?
5. Should we use `value(dual(cons[1]))` or something else?
6. Does the paper explain the dual-primal relationship clearly?

## Root Cause Identified

**BUG**: The examples use `H_R = (X1 + X2)/2` when the correct formula is `H_R = X1 + X2` (no division by 2).

**Consequence**: All moment matrix values are scaled down by factor of 2, causing:
- Constraint violations (unipotency fails: A² = I/2 instead of I)
- Game constraints fail (products off by factors of 1/8)
- Incorrect eigenvalues (all scaled by 1/2)

**Mathematical Justification**: From paper arxiv:2307.11599, the dual formulation (PSDP-R') directly defines H_R = X1 + X2. The relationship X = Y/2 in Theorem 2.2 is about equivalence between formulations, not about extracting H from X.

**Fix**: Remove `/2` from extraction formulas in both affected examples.

## Expected Outputs
1. ✓ Diagnosis of what makes the moment matrices unphysical
2. ✓ Identification of the root cause (extraction, reconstruction, scaling, etc.)
3. ✓ Proposed fix with mathematical justification
4. ✓ Test to verify the fix produces physical moment matrices

## Files to Investigate
- src/sos_solver.jl (dualization logic)
- docs/src/examples/literate/mermin_square_gns.jl (extraction example)
- src/gns.jl (reconstruction logic)
- Paper arxiv:2307.11599 (dual formulation theory)
