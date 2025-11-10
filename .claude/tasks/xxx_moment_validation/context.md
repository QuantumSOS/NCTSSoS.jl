# XXX Moment Matrix Validation Debugging

## Task Overview
Debug why the manually computed moment matrix (`H_manual`) from the exact quantum ground state doesn't match the moment matrix (`H_nctssos`) extracted from NCTSSoS.jl's SDP dual solution.

## Current Status
- **File**: `docs/src/examples/literate/xxx_moment_matrix_validation.jl`
- **Problem**: Moment matrices don't match
  - Diagonal elements: H_nctssos shows 2.0 vs H_manual shows 1.0 (without division by 2)
  - With division by 2: Diagonals match, but off-diagonals are completely wrong
  - Example: H_manual[X₁, X₂] = -0.666... (correct, verified with Yao) vs H_nctssos[X₁, X₂] ≈ 0.0

## Methodology
1. **Manual computation**: H[i,j] = ⟨ψ|basis[i]† basis[j]|ψ⟩ from exact Yao ground state
2. **NCTSSoS extraction**: H = (X1 + X2) + i(X3 - X3') from SDP dual matrix X
3. **Validation**: Compare element-wise with tolerance checks

## Key Questions
1. Is there a normalization factor issue in the dual extraction formula?
2. Does the basis ordering from `get_basis()` match NCTSSoS's internal ordering?
3. Is the extraction formula H = (X1 + X2) + i(X3 - X3') correct?
4. Should there be division by 2 somewhere in the formula?

## Reference Example
- `xxx_pauli_gns.jl` uses H_R = X1 + X2 (no division) but doesn't validate against exact values
- Same extraction approach but no cross-check with manual computation

## Expected Outcome
Identify root cause and provide correct extraction formula or basis transformation to make matrices match within numerical tolerance (~1e-6).

## Root Cause Analysis Complete

**BUG IDENTIFIED**: Missing factor of 1/2 in moment matrix extraction formula.

**Current formula (WRONG)**:
```julia
H_R = X1 + X2
H_I = X3 - X3'
```

**Correct formula**:
```julia
H_R = (X1 + X2) / 2
H_I = (X3 - X3') / 2
```

**Evidence**:
- H_nctssos trace = 182 (should be 91) - exactly 2× too large
- All diagonal elements = 2.0 (should be 1.0) - exactly 2× too large
- Root cause: Dual variable relationship X = Y/2 from arxiv:2307.11599 Theorem 2.2

**Fix locations**:
1. `docs/src/examples/literate/xxx_moment_matrix_validation.jl` line 250-254
2. `docs/src/examples/literate/xxx_pauli_gns.jl` line 133-137
3. Possibly `docs/src/examples/literate/mermin_square_gns.jl`

See `/Users/exaclior/projects/NCTSSoS.jl/.claude/tasks/xxx_moment_validation/plan.md` for complete analysis.
