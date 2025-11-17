# Code Flow: Moment Matrix → SOS Polynomial Constraints

## Visual Code Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│  moment_solver.jl:39-102: moment_relax()                           │
│  Creates MomentProblem with moment variables y_α                    │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  moment_solver.jl:152-158: add_matrix_constraint!()                │
│                                                                      │
│  moment_mtx = [                                                     │
│      sum([coef * monomap[⟨b_i† · poly · b_j⟩]]                     │
│          for (coef, mono) in zip(coefficients(poly), ...))          │
│      for row_idx in basis, col_idx in basis                         │
│  ]                                                                   │
│                                                                      │
│  Each entry M[i,j] is a JuMP AffExpr with .terms field              │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  moment_solver.jl:158: Add constraint to model                     │
│                                                                      │
│  @constraint(model, moment_mtx in PSDCone())                        │
│                                                                      │
│  Stores as: VectorConstraint with .func field containing the       │
│  matrix entries as AffExpr objects                                  │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         │  Returns MomentProblem containing:
                         │  - model (with constraints)
                         │  - constraints (vector of ConstraintRefs)
                         │  - monomap (Dict: monomial → y variable)
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  sos_solver.jl:50: sos_dualize()                                   │
│  Input: MomentProblem                                               │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  sos_solver.jl:54-56: Create dual variables                        │
│                                                                      │
│  for each cons in moment_problem.constraints:                       │
│      G_j = @variable(dual_model, [1:dim, 1:dim] in                 │
│                      (cons.set == Zeros ? Symmetric : PSDCone()))   │
│                                                                      │
│  Primal: M ⪰ 0  →  Dual: G ⪰ 0                                     │
│  Primal: M = 0  →  Dual: G symmetric                               │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  sos_solver.jl:73-79: Initialize polynomial constraints            │
│                                                                      │
│  fα_constraints = [f_α for α in monomials]                         │
│  fα_constraints[1] -= b  # For identity monomial                   │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  sos_solver.jl:81-86: CRITICAL LOOP                                │
│                                                                      │
│  for (j, sdp_constraint) in enumerate(constraints):                │
│      ┌───────────────────────────────────────────────┐             │
│      │  Call get_Cαj()                              │             │
│      └──────────┬────────────────────────────────────┘             │
│                 │                                                   │
│                 ↓                                                   │
│      ┌──────────────────────────────────────────────────────────┐  │
│      │  sos_solver.jl:8-22: get_Cαj()                          │  │
│      │                                                           │  │
│      │  Input: constraint_object(sdp_constraint)                │  │
│      │     → VectorConstraint with .func = [M[1,1], M[1,2], ...] │  │
│      │                                                           │  │
│      │  for (position, cur_expr) in enumerate(constraint.func): │  │
│      │      # cur_expr is an AffExpr: a₁·y₁ + a₂·y₂ + ...       │  │
│      │      for (variable, coef) in cur_expr.terms:             │  │
│      │          α_idx = basis_dict[variable]                    │  │
│      │          row, col = position_to_matrix_indices(position) │  │
│      │          C_αj[(α_idx, row, col)] = coef                  │  │
│      │                                                           │  │
│      │  Returns: Dict{(α_idx, i, j) → coefficient}              │  │
│      └──────────┬────────────────────────────────────────────────  │
│                 │                                                   │
│                 ↓  C_αj returned                                    │
│      ┌──────────────────────────────────────────────────────────┐  │
│      │  Add contributions to polynomial constraints             │  │
│      │                                                           │  │
│      │  for ((α_idx, i, j), coef) in C_αj:                      │  │
│      │      fα_constraints[α_idx] -= coef * G_j[i,j]            │  │
│      │                                                           │  │
│      │  Accumulates: f_α -= Σ_{i,j} C_αj[α,i,j] · G_j[i,j]     │  │
│      └──────────────────────────────────────────────────────────┘  │
│                                                                      │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  sos_solver.jl:87: Enforce polynomial equalities                   │
│                                                                      │
│  @constraint(dual_model, fα_constraints .== 0)                     │
│                                                                      │
│  Final constraints for each α:                                      │
│    f_α - δ_{α,1}·b - Σ_j Σ_{i,k} C_αj[α,i,k]·G_j[i,k] = 0         │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ↓
┌─────────────────────────────────────────────────────────────────────┐
│  sos_solver.jl:90: Return SOSProblem                               │
│                                                                      │
│  Contains:                                                          │
│    - dual_model with objective: max b                               │
│    - Variables: G_j (PSD matrices) and b (scalar)                   │
│    - Constraints: polynomial equalities (one per monomial)          │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Detailed get_Cαj() Function Trace

### Input Structure

When `constraint_object(cons)` is called on a PSD constraint, it returns a `VectorConstraint`:

```julia
struct VectorConstraint{F, S, Shape}
    func::Vector{F}  # Vector of AffExpr, one per matrix entry
    set::S           # PSDCone or Zeros
    shape::Shape     # SquareMatrixShape
end
```

The `.func` field contains the matrix entries **flattened** in column-major order:
```
For 2×2 matrix [ M[1,1]  M[1,2] ]
               [ M[2,1]  M[2,2] ]

func = [M[1,1], M[2,1], M[1,2], M[2,2]]  # Column-major order
```

### Function Execution

```julia
function get_Cαj(basis_dict::Dict, localizing_mtx::VectorConstraint)
    dim = get_dim(localizing_mtx)  # Matrix dimension (e.g., 2 for 2×2)
    cis = CartesianIndices((dim, dim))  # [(1,1), (2,1), (1,2), (2,2)]

    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    # Iterate over flattened entries and their matrix positions
    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        # ci: CartesianIndex with .I[1]=row, .I[2]=col
        # cur_expr: AffExpr with .terms field

        # cur_expr.terms is an iterator of (variable, coefficient) pairs
        for (α, coeff) in cur_expr.terms
            # α: JuMP variable (e.g., y[3])
            # coeff: coefficient (e.g., 1.0 or -2.5)
            # basis_dict[α]: index of this variable in the basis (e.g., 3)

            α_idx = basis_dict[α]
            row = ci.I[1]
            col = ci.I[2]

            dictionary_of_keys[(α_idx, row, col)] = coeff
        end
    end

    return dictionary_of_keys
end
```

### Concrete Example

**Input:**
```julia
# Moment matrix with basis [1, x]
M = [y[1]  y[2];
     y[2]  y[3]]

# After @constraint(model, M in PSDCone())
# constraint_object(cons).func = [y[1], y[2], y[2], y[3]]  (column-major)

# Each entry is an AffExpr:
# y[1].terms = iterator yielding (y[1], 1.0)
# y[2].terms = iterator yielding (y[2], 1.0)
# etc.

basis_dict = Dict(y[1]=>1, y[2]=>2, y[3]=>3)
```

**Execution:**
```
Iteration 1: ci=(1,1), cur_expr=y[1]
  - cur_expr.terms yields (y[1], 1.0)
  - α_idx = basis_dict[y[1]] = 1
  - Store: (1, 1, 1) → 1.0

Iteration 2: ci=(2,1), cur_expr=y[2]
  - cur_expr.terms yields (y[2], 1.0)
  - α_idx = 2
  - Store: (2, 2, 1) → 1.0

Iteration 3: ci=(1,2), cur_expr=y[2]
  - α_idx = 2
  - Store: (2, 1, 2) → 1.0

Iteration 4: ci=(2,2), cur_expr=y[3]
  - α_idx = 3
  - Store: (3, 2, 2) → 1.0
```

**Output:**
```julia
C_αj = Dict(
    (1, 1, 1) => 1.0,
    (2, 2, 1) => 1.0,
    (2, 1, 2) => 1.0,
    (3, 2, 2) => 1.0
)
```

### With Non-trivial Entries

**Input:**
```julia
# Localizing matrix for g = 1 - x² with basis [1, x]
M_g = [y[1] - y[3]    y[2];
       y[2]           y[3]]

# After constraint, entries become:
# M_g[1,1].terms = iterator yielding (y[1], 1.0), (y[3], -1.0)
# M_g[1,2].terms = iterator yielding (y[2], 1.0)
# etc.
```

**Execution for M_g[1,1]:**
```
ci=(1,1), cur_expr = 1.0*y[1] - 1.0*y[3]

cur_expr.terms yields:
  - (y[1], 1.0)   → Store: (1, 1, 1) → 1.0
  - (y[3], -1.0)  → Store: (3, 1, 1) → -1.0
```

**Output:**
```julia
C_αj_g = Dict(
    (1, 1, 1) =>  1.0,   # y[1] in position [1,1]
    (3, 1, 1) => -1.0,   # y[3] in position [1,1] with negative coef
    (2, 2, 1) =>  1.0,   # y[2] in position [2,1]
    (2, 1, 2) =>  1.0,   # y[2] in position [1,2]
    (3, 2, 2) =>  1.0    # y[3] in position [2,2]
)
```

---

## How C_αj Builds Polynomial Constraints

### Accumulation Process

Starting with:
```julia
fα_constraints[1] = f_1 - b  # Objective coefficient for ⟨1⟩ minus b
fα_constraints[2] = f_2      # Objective coefficient for ⟨x⟩
fα_constraints[3] = f_3      # Objective coefficient for ⟨x²⟩
```

For constraint j with dual variable G_j and coefficients C_αj:
```julia
for ((α_idx, i, j), coef) in C_αj:
    fα_constraints[α_idx] -= coef * G_j[i,j]
end
```

**Example with two constraints:**

Constraint 0 (moment matrix): `G_0 ⪰ 0`
```
C_α0 = {(1,1,1)→1, (2,1,2)→1, (2,2,1)→1, (3,2,2)→1}

Adds to:
  fα_constraints[1] -= 1 * G_0[1,1]
  fα_constraints[2] -= 1 * G_0[1,2]
  fα_constraints[2] -= 1 * G_0[2,1]
  fα_constraints[3] -= 1 * G_0[2,2]
```

Constraint 1 (localizing matrix): `G_1 ⪰ 0`
```
C_α1 = {(1,1,1)→1, (3,1,1)→-1, (2,1,2)→1, (2,2,1)→1, (3,2,2)→1}

Adds to:
  fα_constraints[1] -= 1 * G_1[1,1]
  fα_constraints[3] -= (-1) * G_1[1,1]  = 1 * G_1[1,1]
  fα_constraints[2] -= 1 * G_1[1,2]
  fα_constraints[2] -= 1 * G_1[2,1]
  fα_constraints[3] -= 1 * G_1[2,2]
```

**Final constraints:**
```
fα_constraints[1] = f_1 - b - G_0[1,1] - G_1[1,1]
fα_constraints[2] = f_2 - G_0[1,2] - G_0[2,1] - G_1[1,2] - G_1[2,1]
fα_constraints[3] = f_3 - G_0[2,2] + G_1[1,1] - G_1[2,2]
```

Setting these equal to zero gives the dual polynomial constraints!

---

## Key Code Locations Reference

| Function | File | Lines | Purpose |
|----------|------|-------|---------|
| `moment_relax()` | `src/moment_solver.jl` | 39-102 | Creates primal moment problem |
| `add_matrix_constraint!()` | `src/moment_solver.jl` | 141-159 | Builds moment/localizing matrices |
| `sos_dualize()` | `src/sos_solver.jl` | 50-91 | Converts primal to dual |
| `get_Cαj()` | `src/sos_solver.jl` | 8-22 | Extracts coefficient structure |
| Dual variable creation | `src/sos_solver.jl` | 54-56 | Creates G matrices |
| Polynomial constraint init | `src/sos_solver.jl` | 73-79 | Sets up f_α - b |
| Polynomial constraint build | `src/sos_solver.jl` | 81-86 | Adds G contributions |
| Constraint enforcement | `src/sos_solver.jl` | 87 | Creates equalities |

---

## Summary

The conversion happens through **extracting the linear structure** of moment matrix entries:

1. **Moment matrices** have entries that are JuMP affine expressions
2. **`get_Cαj()`** iterates through these expressions and extracts which moment variables appear in which matrix positions
3. **Dual variables** `G_j` replace the primal PSD constraints
4. **Polynomial constraints** are built by accumulating the dual variable contributions weighted by the extracted coefficients
5. The result is a system of **linear equalities** in the dual variables that encode the original matrix constraints

The key insight: **The structure of the moment matrix becomes the coefficients in the polynomial constraints.**
