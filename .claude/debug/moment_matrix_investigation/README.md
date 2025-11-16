# Moment Matrix Discrepancy Investigation

This directory contains diagnostic scripts created to investigate why the moment matrices H_manual and H_nctssos don't match element-wise in the xxx_moment_matrix_validation example.

## Investigation Scripts

### 1. debug_moment_indexing.jl
**Purpose**: Verify that moment matrix indexing M[i,j] = ⟨ψ| basis[i]† * basis[j] |ψ⟩ is correct

**Key Validations**:
- Confirms the indexing formula is implemented correctly
- Verifies H_manual[2,3] matches direct Yao calculation of ⟨X₁X₂⟩
- Shows internal consistency: H_manual[2,3] = H_manual[1,14] = -0.666...

### 2. debug_nctssos_values.jl
**Purpose**: Examine H_nctssos values and search for where expected values appear

**Key Findings**:
- H_nctssos[2,3] ≈ 0.0 (different from H_manual[2,3] = -0.666)
- The value -0.666 does appear in H_nctssos, but at different matrix positions
- Both matrices satisfy all physical constraints (Hermitian, PSD, normalized, correct energy)

### 3. debug_moment_matrix.jl
**Purpose**: Investigate basis construction and expval operations

**Key Validations**:
- Verifies basis construction is correct
- Confirms extraction formula H = (X₁ + X₂)/2 + i(X₃ - X₃')/2
- Shows no issues with basis ordering or simplification

## Conclusion

The investigation confirmed that the discrepancy is **expected behavior, not a bug**. The SDP finds a different quantum state than exact diagonalization, both with the same ground state energy E₀ = -8.0. This is due to:
- Ground state symmetries (SU(2) rotational symmetry in XXX model)
- SDP's freedom to choose among equivalent optimal solutions
- Multiple valid states in degenerate/symmetric subspace

See `MOMENT_MATRIX_DISCREPANCY_EXPLAINED.md` in the project root for comprehensive explanation.

## How to Run

```bash
cd /Users/exaclior/projects/NCTSSoS.jl
julia .claude/debug/moment_matrix_investigation/debug_moment_indexing.jl
julia .claude/debug/moment_matrix_investigation/debug_nctssos_values.jl
julia .claude/debug/moment_matrix_investigation/debug_moment_matrix.jl
```

## Date
2025-11-10
