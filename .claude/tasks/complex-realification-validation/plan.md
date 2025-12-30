# Complex SDP Realification Validation - Research Analysis

## Executive Summary

After thorough analysis of both the NCTSSoS.jl codebase and the paper arxiv:2307.11599, **I have found a significant discrepancy** in how complex Hermitian positive semidefinite matrices are converted to real positive semidefinite matrices. The codebase uses the **standard/popular approach** documented in equation (1.1) of the paper, while the paper proposes a **new, more efficient dual-based reformulation** (PSDP-R') that avoids extra affine constraints.

**Key Finding**: The codebase is using the older, less efficient reformulation (PSDP-R) rather than the paper's new efficient reformulation (PSDP-R'). This is not mathematically incorrect, but it is computationally inefficient.

---

## Section 1: Current Codebase Implementation

### 1.1 Location of Realification Logic

The complex to real SDP conversion is implemented in:
- **File**: `/Users/exaclior/projects/NCTSSoS.jl/src/sos_solver.jl`
- **Function**: `sos_dualize(cmp::ComplexMomentProblem)` (lines 92-146)

### 1.2 Mathematical Approach in Codebase

The codebase uses the **standard block matrix reformulation**:

For a complex Hermitian matrix `H = H_R + i·H_I ∈ ℂⁿˣⁿ`, the equivalence used is:

```
H ⪰ 0  ⟺  Y = [H_R  -H_I] ⪰ 0
              [H_I   H_R]
```

This produces a **2n×2n real symmetric matrix** from an **n×n complex Hermitian matrix**.

### 1.3 Implementation Details

**From lines 92-116** of `sos_solver.jl`:

```julia
function sos_dualize(cmp::ComplexMomentProblem{T,P}) where {T,P}
    dual_model = GenericModel{real(T)}()

    # Create 2n×2n PSD variables for each constraint
    dual_variables = map(cmp.constraints) do (type,cons)
        G_dim = size(cons,1)
        @variable(dual_model, [1:2*G_dim, 1:2*G_dim] in (type == :Zero ? SymmetricMatrixSpace() : PSDCone()))
    end
    dual_variable_dims = map(dual_variables) do dv
        size(dv, 1) ÷ 2
    end

    # Extract the block components
    # X1 + X2 and X3 - X3^T
    Xs = [[
            begin
                dim = dual_variable_dims[i]
                dv[1:dim, 1:dim] .+ dv[dim+1:2*dim, dim+1:2*dim]  # X1 + X2
            end for (i, dv) in enumerate(dual_variables)
        ], [
            begin
                dim = dual_variable_dims[i]
                dv[dim+1:2*dim, 1:dim] .- dv[1:dim, 1+dim:2*dim]  # X3 - X3^T
            end for (i, dv) in enumerate(dual_variables)
        ]
    ]
```

This creates the 2n×2n matrix:
```
Y = [X1  X3^T]
    [X3  X2 ]
```

Then extracts:
- Real part coefficient: `X1 + X2`
- Imaginary part coefficient: `X3 - X3^T`

**From lines 124-145**, the constraint assembly:

```julia
# real and imag parts of fα constraints
fα_constraints = [[zero(...) for _ in 1:length(symmetric_basis)],
                  [zero(...) for _ in 1:length(symmetric_basis)]]

for (coef,mono) in terms(cmp.objective)
    for (fα_constraints_part, part_func) in zip(fα_constraints, [real, imag])
        fα_constraints_part[searchsortedfirst(symmetric_basis, mono)] += part_func(coef)
    end
end

add_to_expression!(fα_constraints[1][1], -one(T), b)

for (i, (_,sdp_constraint)) in enumerate(cmp.constraints)
    Cαjs = get_Cαj(cmp.total_basis, sdp_constraint)
    for (ky, coef) in Cαjs
        # This is the key transformation logic
        for (X_part, coef_part, sign, part_func) in zip([1, 2, 2, 1], [1, 1, 2, 2], [-1, -1, -1, 1], [real, imag, real, imag])
            add_to_expression!(fα_constraints[coef_part][ky[1]],
                              sign*part_func(coef), Xs[X_part][i][ky[2], ky[3]])
        end
    end
end
@constraint(dual_model, fα_constraints[1] .== 0)
@constraint(dual_model, fα_constraints[2] .== 0)
```

### 1.4 Constraint Loop Analysis

The critical loop on **line 138** iterates through:
```julia
for (X_part, coef_part, sign, part_func) in
    zip([1, 2, 2, 1], [1, 1, 2, 2], [-1, -1, -1, 1], [real, imag, real, imag])
```

This produces 4 terms for each coefficient:
1. `X_part=1, coef_part=1, sign=-1, part_func=real`: `-real(coef) * Xs[1][i]` → `-real(coef) * (X1 + X2)`
2. `X_part=2, coef_part=1, sign=-1, part_func=imag`: `-imag(coef) * Xs[2][i]` → `-imag(coef) * (X3 - X3^T)`
3. `X_part=2, coef_part=2, sign=-1, part_func=real`: `-real(coef) * Xs[2][i]` → `-real(coef) * (X3 - X3^T)`
4. `X_part=1, coef_part=2, sign=+1, part_func=imag`: `+imag(coef) * Xs[1][i]` → `+imag(coef) * (X1 + X2)`

---

## Section 2: Paper's Approach (arxiv:2307.11599)

### 2.1 The Standard Reformulation (PSDP-R)

The paper describes the **standard popular reformulation** in Section 2, equation (1.1):

```
H ⪰ 0  ⟺  Y = [H_R  -H_I] ⪰ 0
              [H_I   H_R]
```

With **additional affine constraints** (equation 1.2):

```
Y_{i,j} = Y_{i+n,j+n},           i=1,...,n, j=i,...,n
Y_{i,j+n} + Y_{j,i+n} = 0,       i=1,...,n, j=i,...,n
```

This requires **n(n+1) extra affine constraints**, which is computationally expensive.

### 2.2 The New Efficient Reformulation (PSDP-R')

The paper's **main contribution** (Theorem 2.2, page 4) is a **dual-based reformulation** that avoids extra constraints.

For a complex Hermitian matrix `H = H_R + i·H_I`, define the dual variable:

```
X = [X1  X3^T]  ∈ S^{2n}
    [X3  X2 ]
```

**Key transformation**:
```
H_R = X1 + X2
H_I = X3 - X3^T
```

**Recovery** (Theorem 2.2):
If `X*` is optimal for (PSDP-R'), then:
```
H* = (X1* + X2*) + i(X3* - (X3*)^T)
```
is optimal for the original complex problem.

### 2.3 Mathematical Equivalence

**Theorem 2.2** proves that (PSDP-R') and (PSDP-R) share the same optimum by showing:

1. **Primal to Dual**: If `Y = [H_R -H_I; H_I H_R]` is feasible for (PSDP-R), then `X = Y/2` is feasible for (PSDP-R') with same objective value.

2. **Dual to Primal**: If `X = [X1 X3^T; X3 X2]` is feasible for (PSDP-R'), then construct:
   ```
   Y = [X1+X2    X3^T-X3]
       [X3-X3^T  X1+X2  ]
   ```
   This `Y` is PSD (proven via rotation matrix transformation) and feasible for (PSDP-R).

### 2.4 Complexity Comparison

From **Table 1** (page 8):

| Formulation | Matrix Size | Number of Affine Constraints |
|-------------|-------------|------------------------------|
| (HSOS-R) - Standard | 2ω_{s,d} | 2ω²_{s,d} + 2ω_{s,d} + Σᵢω_{s,d-dᵢ} |
| (HSOS-R') - New | 2ω_{s,d} | ω²_{s,d} |

The new reformulation **halves** the number of affine constraints!

### 2.5 Numerical Evidence

Table 2 (page 9) shows **2-7x speedup** for the new reformulation on various problem sizes.

---

## Section 3: Comparison and Discrepancy Analysis

### 3.1 Which Approach Does the Codebase Use?

Examining the code in `sos_solver.jl` lines 92-146:

**Evidence that codebase uses the STANDARD approach (not the paper's new approach)**:

1. **Matrix Structure**: The code creates `[1:2*G_dim, 1:2*G_dim]` PSD variables directly, suggesting the Y-form.

2. **Block Extraction**: Lines 105-116 extract `X1+X2` and `X3-X3^T`, which matches the paper's dual formulation extraction.

3. **No Extra Constraints**: The code does **NOT** add the n(n+1) extra affine constraints from equation (1.2).

**CONCLUSION**: The codebase is actually using **a hybrid approach** that appears closer to (PSDP-R') but may not be exploiting all its benefits.

### 3.2 Detailed Mathematical Comparison

Let me trace through exactly what the codebase does vs. the paper:

#### Codebase Constraint Construction

For a complex coefficient `c = c_R + i·c_I` in the moment matrix at position (α, row, col):

The loop creates:
```
fα_real[α] += -real(c) * (X1 + X2)[row,col] - imag(c) * (X3 - X3^T)[row,col]
fα_imag[α] += -real(c) * (X3 - X3^T)[row,col] + imag(c) * (X1 + X2)[row,col]
```

#### Paper's Formulation (HSOS-R')

From equation (3.3) on page 7, the paper constructs:
```
Σᵢ [AR_i(X1_i + X2_i)]_{β,γ} + [AI_i(X3_i - (X3_i)^T)]_{β,γ} + δ_{(β,γ),(0,0)}λ = [b_R]_{β,γ}
Σᵢ [AR_i(X3_i - (X3_i)^T)]_{β,γ} - [AI_i(X1_i + X2_i)]_{β,γ} = [b_I]_{β,γ}
```

Where `AR` and `AI` are linear operators extracting real and imaginary parts.

### 3.3 The Discrepancy

**THE APPROACHES MATCH!** The codebase is implementing the paper's efficient (PSDP-R') formulation, not the standard (PSDP-R) formulation!

**Key Evidence**:
1. The codebase uses `X1 + X2` and `X3 - X3^T` exactly as in (PSDP-R')
2. No extra n(n+1) affine constraints are added
3. The constraint structure matches equation (3.3) from the paper

### 3.4 Sign Convention Analysis

Let me verify the signs more carefully:

**Codebase** (line 138-140):
```julia
for (X_part, coef_part, sign, part_func) in zip([1, 2, 2, 1], [1, 1, 2, 2], [-1, -1, -1, 1], [real, imag, real, imag])
    add_to_expression!(fα_constraints[coef_part][ky[1]],
                      sign*part_func(coef), Xs[X_part][i][ky[2], ky[3]])
end
```

This produces:
- `fα_real`: `-real(coef)·(X1+X2) - imag(coef)·(X3-X3^T)`
- `fα_imag`: `-real(coef)·(X3-X3^T) + imag(coef)·(X1+X2)`

**Paper** (equation 3.3):
- Real constraint: `Σ[AR(X1+X2)] + [AI(X3-X3^T)] = bR`
- Imag constraint: `Σ[AR(X3-X3^T)] - [AI(X1+X2)] = bI`

The signs MATCH when accounting for the fact that:
- `AR(coef) = real(coef)`
- `AI(coef) = imag(coef)`
- The codebase moves terms to the left side (hence the negatives)

---

## Section 4: Final Assessment

### 4.1 Is the Codebase Correct?

**YES!** The codebase correctly implements the paper's efficient (PSDP-R') reformulation.

### 4.2 Why Might the User Suspect Incorrectness?

Possible reasons for suspicion:
1. **Lack of documentation**: The code doesn't cite the paper or explain which formulation it uses
2. **Complex sign manipulation**: The 4-tuple loop is not immediately obvious
3. **No validation tests**: No test comparing against a reference implementation

### 4.3 Potential Issues

While mathematically correct, there are some concerns:

1. **Code Clarity**: The 4-tuple `zip([1,2,2,1], [1,1,2,2], [-1,-1,-1,1], [real,imag,real,imag])` is opaque and hard to verify

2. **No Reference to Paper**: The code predates or doesn't reference the 2023 paper, suggesting independent derivation

3. **Implicit Block Structure**: The extraction of `X1+X2` and `X3-X3^T` is correct but uncommented

### 4.4 Validation Strategy

To validate correctness, one should:

1. **Test on Known Problems**: Compare results with a known-correct complex SDP solver (Hypatia, SeDuMi, etc.)

2. **Check Duality Gap**: Verify strong duality holds

3. **Test Recovery**: Verify that the recovered complex solution satisfies original complex constraints

4. **Numerical Stability**: Check that the real formulation doesn't introduce numerical errors

---

## Section 5: Recommendations

### 5.1 No Changes Needed to Math

The implementation is **mathematically correct** and uses the **state-of-the-art efficient formulation** from the 2023 paper.

### 5.2 Suggested Improvements for Clarity

If you want to improve the code (without changing functionality):

#### Improvement 1: Add Documentation

```julia
# This function implements the efficient complex-to-real SDP reformulation
# from "A more efficient reformulation of complex SDP as real SDP"
# (Wang, arxiv:2307.11599, 2023)
#
# Key idea: Instead of using H ⪰ 0 ⟺ [H_R -H_I; H_I H_R] ⪰ 0 (standard)
# we use the dual formulation X = [X1 X3^T; X3 X2] ⪰ 0 where:
#   H_R = X1 + X2
#   H_I = X3 - X3^T
# This avoids n(n+1) extra affine constraints.
function sos_dualize(cmp::ComplexMomentProblem{T,P}) where {T,P}
    ...
```

#### Improvement 2: Clarify the 4-Term Loop

Instead of the opaque 4-tuple, use named operations:

```julia
# Build constraints from complex coefficients
# For c = c_R + i·c_I, we need:
#   Real part: -c_R·(X1+X2) - c_I·(X3-X3^T)
#   Imag part: -c_R·(X3-X3^T) + c_I·(X1+X2)
for (ky, coef) in Cαjs
    c_R = real(coef)
    c_I = imag(coef)
    row, col = ky[2], ky[3]

    # Real part constraints
    add_to_expression!(fα_constraints[1][ky[1]], -c_R, Xs[1][i][row, col])  # -c_R·(X1+X2)
    add_to_expression!(fα_constraints[1][ky[1]], -c_I, Xs[2][i][row, col])  # -c_I·(X3-X3^T)

    # Imaginary part constraints
    add_to_expression!(fα_constraints[2][ky[1]], -c_R, Xs[2][i][row, col])  # -c_R·(X3-X3^T)
    add_to_expression!(fα_constraints[2][ky[1]], +c_I, Xs[1][i][row, col])  # +c_I·(X1+X2)
end
```

#### Improvement 3: Add Validation Test

Create a test that:
1. Solves a small complex SDP both ways (if Hypatia is available)
2. Compares the objectives
3. Verifies the recovered solution satisfies complex constraints

### 5.3 Performance Validation

To confirm the implementation achieves the paper's efficiency gains:

1. **Benchmark against standard formulation**: Implement (PSDP-R) with extra constraints and compare
2. **Count actual constraints generated**: Verify it matches ω²_{s,d} not 2ω²_{s,d} + 2ω_{s,d}
3. **Profile solver time**: Confirm speedup on larger problems

---

## Section 6: References and Code Locations

### Key Source Files

1. **Complex SDP Dualization**:
   - File: `src/sos_solver.jl`
   - Function: `sos_dualize(cmp::ComplexMomentProblem)` (lines 92-146)
   - Key logic: Lines 105-116 (block extraction), 138-140 (constraint assembly)

2. **Complex Moment Problem Construction**:
   - File: `src/complex_moment_solver.jl`
   - Function: `moment_relax(cpop::ComplexPolyOpt, ...)` (lines 9-43)
   - Function: `constrain_moment_matrix(...)` (lines 45-56)

3. **Complex Polynomial Optimization Setup**:
   - File: `src/pop.jl`
   - Struct: `ComplexPolyOpt` (lines 95-103)
   - Function: `cpolyopt(...)` (lines 106-126, 154-172)

4. **Interface**:
   - File: `src/interface.jl`
   - Function: `cs_nctssos(...)` (lines 73-97)
   - Note: Line 91 prevents non-dualized solving for complex problems

### Test Files Using Complex Optimization

1. `test/moment_solver.jl`: Lines 28, 62 use `cpolyopt`
2. `test/interface.jl`: Lines 24, 42, 59, 76 use `cpolyopt`
3. `test/heisenberg.jl`: Lines 26, 52, 73 use `cpolyopt`
4. `test/pop.jl`: Tests for `ComplexPolyOpt` constructor
5. `test/algebra_constructors.jl`: Lines 89, 113, 129, etc.

### Paper Reference

- **Title**: "A more efficient reformulation of complex SDP as real SDP"
- **Author**: Jie Wang
- **ArXiv**: [2307.11599v3](https://arxiv.org/abs/2307.11599)
- **Date**: July 2023 (v1), July 2024 (v3)
- **Key Equations**: (1.1), (2.1)-(2.3), Theorem 2.2, (3.3), Table 1

---

## Conclusion

**The NCTSSoS.jl codebase correctly implements the efficient complex-to-real SDP reformulation that matches the state-of-the-art approach described in the 2023 paper by Jie Wang.**

The implementation:
- ✅ Uses the dual-based formulation (PSDP-R')
- ✅ Avoids extra n(n+1) affine constraints
- ✅ Correctly extracts H_R = X1 + X2 and H_I = X3 - X3^T
- ✅ Properly assembles real and imaginary constraints
- ✅ Should achieve the efficiency gains reported in the paper

**No mathematical corrections are needed.** The only improvements recommended are for code clarity and validation testing.
