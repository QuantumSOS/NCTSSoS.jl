# Review: src/optimization/sos.jl

**Lines**: 460
**Purpose**: SOS (Sum-of-Squares) dualization of moment problems

## Structure

```
SOSProblem{T}                    # Wrapper for JuMP model
get_Cαj()                        # Extract coefficient matrices
sos_dualize()                    # Entry point (dispatches on type)
├── _sos_dualize_real()          # Real algebras: PSD cone
├── _sos_dualize_hermitian()     # Complex algebras: 2n×2n embedding
└── _sos_dualize_state()         # State polynomials: PSD cone
```

## Dualization Algorithm

### Real Case (lines 152-218)

Standard SOS dual:
1. Create matrix variables `G_j` for each constraint
   - `:Zero` → SymmetricMatrixSpace
   - `:PSD` → PSDCone
2. Create scalar `b` (bound variable)
3. For each basis monomial α: `∑_j ∑_{k,l} C_α_jkl * G_j[k,l] = c_α - δ_{α,1} * b`
4. Objective: `max b`

### Hermitian Case (lines 343-459)

For complex algebras (Pauli, Fermionic, Bosonic):
1. Create 2n×2n PSD matrices (embedding)
2. Extract real/imag parts:
   ```julia
   X1 = dv[1:n, 1:n] + dv[n+1:2n, n+1:2n]      # Re(G)
   X2 = dv[n+1:2n, 1:n] - dv[1:n, n+1:2n]      # Im(G)
   ```
3. Separate real/imag constraint arrays
4. Complex arithmetic: `(c_re + i*c_im) * (X1 + i*X2)`

### State Case (lines 252-325)

Similar to real case but:
- Uses `StateWord` (via `expval`) for basis comparison
- Directly iterates over upper triangle with 2× factor for off-diagonal

## Quality Assessment

| Aspect | Rating | Notes |
|--------|--------|-------|
| Documentation | ✓✓✓ | Clear docstrings, math explanation |
| Correctness | ✓✓✓ | Hermitian embedding verified |
| Code structure | ✓✓ | Three separate implementations |
| Performance | ✓✓ | Sparse coefficient extraction good |

## Key Observations

### Strengths

1. **Sparse coefficient extraction** (`get_Cαj`, lines 61-83):
   Only stores non-zero coefficients in dictionary. Efficient for sparse problems.

2. **Hermitian embedding** (lines 363-382):
   ```julia
   X1 = top-left + bottom-right       # Real part
   X2 = bottom-left - top-right       # Imaginary part
   ```
   Correct extraction from 2n×2n embedding.

3. **Basis symmetrization** (line 174):
   ```julia
   symmetric_basis = _sorted_symmetric_basis(mp.total_basis)
   ```
   Properly handles fact that moment basis may not be symmetric-canonical.

4. **Upper triangle iteration** (lines 303-318):
   ```julia
   effective_coeff = (row == col) ? coeff : 2 * coeff
   ```
   Correctly handles symmetric matrix contributions.

### Concerns

1. **Hermitian math comments** (lines 430-438):
   The comment about complex multiplication is slightly misleading:
   ```julia
   # coef * G = (c_re + i*c_im) * (X1 + i*X2)
   # Real part: c_re*X1 - c_im*X2
   # Imag part: c_im*X1 + c_re*X2
   ```
   But the code uses negative signs for the dual. The code is correct, but the
   comment derivation is confusing without explaining the dual relationship.

2. **Three separate implementations**:
   - `_sos_dualize_real`
   - `_sos_dualize_hermitian`
   - `_sos_dualize_state`

   Some shared logic (basis creation, constraint loop) could be factored out.

3. **Missing StateWord in hermitian case**:
   `_sos_dualize_hermitian` doesn't handle StatePolyOpt with complex algebras.
   This might be intentional (state polys typically use real algebras).

## Dual Formulation Reference

**Primal (moment)**:
```
min  ⟨c, y⟩
s.t. y_1 = 1
     M(y) ≽ 0        (moment matrix PSD)
     L_g(y) ≽ 0      (localizing matrices PSD)
```

**Dual (SOS)**:
```
max  b
s.t. p(x) - b = ∑_j σ_j(x) * g_j(x)   (SOS decomposition)
     σ_j(x) ≽ 0                        (SOS multipliers)
```

The coefficient matching `∑_j C_αj * G_j = c_α - δ_{α,1} * b` enforces polynomial equality.

## Code Quality: Good

Correct implementation. The Hermitian embedding is the most complex part and is handled well. Code duplication could be reduced.
