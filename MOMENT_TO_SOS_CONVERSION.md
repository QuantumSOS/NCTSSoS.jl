# How Moment Matrices Convert to SOS Polynomial Constraints

## The Core Question

**In the primal `MomentProblem`:** We have moment matrices in each clique/block made of JuMP variables.

**In the dual `SOSProblem`:** We need polynomial equality constraints.

**How does this conversion work?**

---

## Step-by-Step Conversion Process

### Step 1: Understand the Primal Moment Matrix Structure

In `moment_solver.jl:152-158`, moment matrices are created as:

```julia
moment_mtx = [
    sum([T_prom(coef) * monomap[simplify!(expval(_neat_dot3(row_idx, mono, col_idx)), sa)]
         for (coef, mono) in zip(coefficients(poly), monomials(poly))])
    for row_idx in local_basis, col_idx in local_basis
]
```

**What this creates:**
- A matrix where each entry `M[i,j]` is a **JuMP affine expression**
- Each expression is a linear combination of moment variables: `M[i,j] = Σ_α coef_α · y_α`

**Example:** For basis `[1, x]` and polynomial `poly = 1` (moment matrix):
```
M[1,1] = 1.0 * y_⟨1⟩         # ⟨1·1·1⟩ = ⟨1⟩
M[1,2] = 1.0 * y_⟨x⟩         # ⟨1·1·x⟩ = ⟨x⟩
M[2,1] = 1.0 * y_⟨x⟩         # ⟨x·1·1⟩ = ⟨x⟩  (symmetric)
M[2,2] = 1.0 * y_⟨x²⟩        # ⟨x·1·x⟩ = ⟨x²⟩
```

For a localizing matrix with `poly = 1 - x²`:
```
M_g[1,1] = 1.0*y_⟨1⟩ - 1.0*y_⟨x²⟩      # ⟨1·(1-x²)·1⟩
M_g[1,2] = 1.0*y_⟨x⟩ - 1.0*y_⟨x³⟩      # ⟨1·(1-x²)·x⟩
M_g[2,1] = 1.0*y_⟨x⟩ - 1.0*y_⟨x³⟩      # ⟨x·(1-x²)·1⟩
M_g[2,2] = 1.0*y_⟨x²⟩ - 1.0*y_⟨x⁴⟩     # ⟨x·(1-x²)·x⟩
```

**Key observation:** Each matrix entry is a **linear combination** of moment variables y_α.

**The primal constraint:**
```julia
@constraint(model, moment_mtx in PSDCone())
```

This adds: `M ⪰ 0` to the primal problem.

---

### Step 2: Extract the Coefficient Structure with `get_Cαj()`

**Location:** `sos_solver.jl:8-22`

```julia
function get_Cαj(basis_dict::Dict, localizing_mtx::VectorConstraint)
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        for (α, coeff) in cur_expr.terms  # <-- KEY LINE
            dictionary_of_keys[(basis_dict[α], ci.I[1], ci.I[2])] = coeff
        end
    end

    return dictionary_of_keys
end
```

**What this does:**

1. **Input:** A constraint whose matrix entries are JuMP affine expressions
   - Each expression has the form: `cur_expr = Σ coeff_k * variable_k`
   - Accessed via `cur_expr.terms` which gives `(variable, coefficient)` pairs

2. **Iteration:** For each position `(i,j)` in the matrix:
   - Look at the JuMP expression `M[i,j]`
   - Extract all terms: `(variable, coefficient)` pairs
   - For each variable that appears, record: `(monomial_index, row, col) → coefficient`

3. **Output:** Dictionary `C_αj` mapping `(α_idx, i, j) → coeff`
   - Meaning: "Moment variable y_α appears in position M[i,j] with coefficient coeff"

**Example:** Given the moment matrix from Step 1:
```
M[1,1] = 1.0 * y_1
M[1,2] = 1.0 * y_2
M[2,2] = 1.0 * y_4
```

`get_Cαj()` returns:
```
C_αj = {
    (1, 1, 1) → 1.0    # y_1 appears in M[1,1] with coef 1.0
    (2, 1, 2) → 1.0    # y_2 appears in M[1,2] with coef 1.0
    (2, 2, 1) → 1.0    # y_2 appears in M[2,1] with coef 1.0
    (4, 2, 2) → 1.0    # y_4 appears in M[2,2] with coef 1.0
}
```

For the localizing matrix:
```
M_g[1,1] = 1.0*y_1 - 1.0*y_4
```

`get_Cαj()` returns:
```
C_αj_g = {
    (1, 1, 1) →  1.0    # y_1 in M_g[1,1] with coef 1.0
    (4, 1, 1) → -1.0    # y_4 in M_g[1,1] with coef -1.0
    ...
}
```

**Mathematical interpretation:**

The coefficient dictionary `C_αj` encodes the decomposition:
```
M[i,j] = Σ_α C_αj[α, i, j] · y_α
```

This is the **key insight**: The moment matrix constraint can be written as a linear function of moment variables, parameterized by the coefficient structure C_αj.

---

### Step 3: Create Dual Variables

**Location:** `sos_solver.jl:54-56`

```julia
dual_variables = map(constraint_object.(moment_problem.constraints)) do cons
    G_dim = get_dim(cons)
    @variable(dual_model, [1:G_dim, 1:G_dim] in
              ((cons.set isa MOI.Zeros) ? SymmetricMatrixSpace() : PSDCone()))
end
```

**What happens:**

- For each primal constraint (moment or localizing matrix):
  - **Primal:** `M ⪰ 0` (PSD cone) → **Dual:** `G ⪰ 0` (PSD variable)
  - **Primal:** `M = 0` (zeros) → **Dual:** `G` symmetric (no PSD)

- The dimension of `G` matches the dimension of `M`

**Example:**
- Primal has: `M ⪰ 0` (3×3 moment matrix) and `M_g ⪰ 0` (3×3 localizing matrix)
- Dual gets: `G_0 ⪰ 0` (3×3) and `G_1 ⪰ 0` (3×3)

---

### Step 4: Initialize Polynomial Constraints

**Location:** `sos_solver.jl:73-79`

```julia
primal_objective_terms = objective_function(moment_problem.model).terms

fα_constraints = [GenericAffExpr{T,VariableRef}(
                      get(primal_objective_terms, α, zero(T)))
                  for α in symmetric_variables]

add_to_expression!(fα_constraints[1], -one(T), b)
```

**What this does:**

1. Extract the primal objective: `min Σ_α f_α · y_α`
2. Create one constraint expression per unique monomial α: `f_α`
3. For the identity monomial (α=1), subtract `b`: `f_1 - b`

**Example:** If objective is `f = 2·y_1 - y_4 - y_6` (representing 2 - x² - y²):
```
fα_constraints[1] = 2 - b     # For monomial ⟨1⟩
fα_constraints[2] = 0         # For monomial ⟨x⟩
fα_constraints[3] = 0         # For monomial ⟨y⟩
fα_constraints[4] = -1        # For monomial ⟨x²⟩
fα_constraints[5] = 0         # For monomial ⟨xy⟩
fα_constraints[6] = -1        # For monomial ⟨y²⟩
```

---

### Step 5: Add Dual Variable Contributions to Polynomial Constraints

**Location:** `sos_solver.jl:81-86`

```julia
for (i, sdp_constraint) in enumerate(moment_problem.constraints)
    Cαjs = get_Cαj(unsymmetrized_basis_vals_dict, constraint_object(sdp_constraint))
    for (ky, coef) in Cαjs
        add_to_expression!(fα_constraints[symmetrized_α2cons_dict[unsymmetrized_basis[ky[1]]]],
                           -coef, dual_variables[i][ky[2], ky[3]])
    end
end
```

**This is the CRITICAL step!**

For each primal constraint (each moment/localizing matrix):
1. Extract its coefficient structure: `C_αj = get_Cαj(...)`
2. For each entry `(α_idx, i, j) → coef` in `C_αj`:
   - Find the polynomial constraint for monomial α: `fα_constraints[α_idx]`
   - Add the term: `-coef * G[i,j]`

**Example:** Using C_αj from Step 2:

From moment matrix M (constraint index 0):
```
C_αj = {(1,1,1)→1, (2,1,2)→1, (2,2,1)→1, (4,2,2)→1}
```

This adds to polynomial constraints:
```
fα_constraints[1] -= 1.0 * G_0[1,1]  # From (1,1,1)→1
fα_constraints[2] -= 1.0 * G_0[1,2]  # From (2,1,2)→1
fα_constraints[2] -= 1.0 * G_0[2,1]  # From (2,2,1)→1
fα_constraints[4] -= 1.0 * G_0[2,2]  # From (4,2,2)→1
```

From localizing matrix M_g (constraint index 1):
```
C_αj_g = {(1,1,1)→1, (4,1,1)→-1, ...}
```

This adds:
```
fα_constraints[1] -= 1.0 * G_1[1,1]   # From (1,1,1)→1
fα_constraints[4] -= (-1.0) * G_1[1,1] = 1.0 * G_1[1,1]  # From (4,1,1)→-1
```

**After accumulating all contributions:**
```
fα_constraints[1] = (2 - b) - G_0[1,1] - G_1[1,1]
fα_constraints[2] = 0 - G_0[1,2] - G_0[2,1] - ...
fα_constraints[4] = -1 - G_0[2,2] + G_1[1,1] - ...
...
```

---

### Step 6: Enforce Polynomial Equalities

**Location:** `sos_solver.jl:87`

```julia
@constraint(dual_model, fα_constraints .== 0)
```

**Final dual constraints:**

For each monomial α:
```
f_α - δ_{α,1}·b = Σ_j Σ_{i,k} C_αj[α,i,k] · G_j[i,k]
```

Rearranged:
```
f_α - δ_{α,1}·b - Σ_j (G_j • C_αj) = 0
```

**These are linear equality constraints** in the dual variables `G_j[i,k]` and `b`.

---

## Complete Example Walkthrough

### Primal Problem

**Variables:** `y = [y_1, y_x, y_x²]` representing `[⟨1⟩, ⟨x⟩, ⟨x²⟩]`

**Objective:** `min 2·y_1 - y_x²` (i.e., `min 2 - x²`)

**Constraint:** Moment matrix with basis `[1, x]`
```
M = [ y_1   y_x  ]  ⪰ 0
    [ y_x   y_x² ]
```

**Normalization:** `y_1 = 1`

### Extracting C_αj

Matrix entries:
```
M[1,1] = y_1    → C_αj: (1, 1, 1) → 1
M[1,2] = y_x    → C_αj: (2, 1, 2) → 1
M[2,1] = y_x    → C_αj: (2, 2, 1) → 1
M[2,2] = y_x²   → C_αj: (3, 2, 2) → 1
```

### Dual Problem Setup

**Variables:**
- `G[1:2, 1:2] ⪰ 0` (dual matrix for M ⪰ 0)
- `b` (scalar bound)

**Initial polynomial constraints:**
```
fα_constraints[1] = 2 - b      # For ⟨1⟩
fα_constraints[2] = 0          # For ⟨x⟩
fα_constraints[3] = -1         # For ⟨x²⟩
```

**Add dual contributions** from C_αj:
```
fα_constraints[1] -= G[1,1]    # (1,1,1)→1
fα_constraints[2] -= G[1,2]    # (2,1,2)→1
fα_constraints[2] -= G[2,1]    # (2,2,1)→1
fα_constraints[3] -= G[2,2]    # (3,2,2)→1
```

**Final constraints:**
```
2 - b - G[1,1] = 0              ⟹ G[1,1] = 2 - b
0 - G[1,2] - G[2,1] = 0         ⟹ G[1,2] + G[2,1] = 0
-1 - G[2,2] = 0                 ⟹ G[2,2] = -1
```

**Problem:** `G[2,2] = -1` but `G ⪰ 0` requires `G[2,2] ≥ 0`!

This is **infeasible**, which means we need to lower `b` or add inequality constraints to make it feasible.

### With Inequality Constraint

Add constraint: `g = 1 - x² ≥ 0`

**Localizing matrix for g:**
```
M_g = [ y_1 - y_x²   y_x - y_x³  ]  ⪰ 0
      [ y_x - y_x³   y_x² - y_x⁴ ]
```

Extract C_αj_g:
```
M_g[1,1] = y_1 - y_x²  → (1,1,1)→1, (3,1,1)→-1
...
```

Create second dual variable: `G_1 ⪰ 0`

Updated constraints:
```
2 - b - G_0[1,1] - G_1[1,1] = 0
...
-1 - G_0[2,2] + G_1[1,1] = 0   ⟹ G_0[2,2] = G_1[1,1] - 1
```

Now if `G_1[1,1] ≥ 1`, we can have `G_0[2,2] ≥ 0`!

**Dual problem:**
```
maximize b
subject to:
  G_0[1,1] + G_1[1,1] = 2 - b
  G_0[2,2] - G_1[1,1] = -1
  (other constraints...)
  G_0 ⪰ 0
  G_1 ⪰ 0
```

The optimal `b* = 1` (achieved when `x = 0`, giving `f(0) = 2 - 0 = 2` constrained by `1 - 0 ≥ 0`).

---

## Summary: The Conversion Process

### Data Flow

```
Moment Matrix M (primal)
    ↓
Each entry M[i,j] is a JuMP AffExpr: Σ_α coef_α · y_α
    ↓
get_Cαj() extracts: (α_idx, i, j) → coef_α
    ↓
Create dual variable G (PSD matrix)
    ↓
For each monomial α, build polynomial constraint:
    f_α - δ_{α,1}·b = Σ_{i,j} C_αj[α,i,j] · G[i,j]
    ↓
Add as equality constraints to dual model
```

### Key Insights

1. **Primal constraint** `M ⪰ 0` becomes **dual variable** `G ⪰ 0`

2. **Structure of M** (which moments appear where) becomes **polynomial constraints**

3. **`get_Cαj()`** is the bridge: it extracts the linear structure of the moment matrix entries

4. **Each unique monomial** gets one polynomial equality constraint in the dual

5. **Objective coefficients** appear on the LHS of polynomial constraints

6. **Dual variable entries** `G[i,j]` appear on the RHS, weighted by coefficients from C_αj

### Why This Works (SDP Duality)

The primal problem is:
```
min  c'y
s.t. Σ_α A_α y_α ⪰ 0
     y_1 = 1
```

where `A_α` are the "coefficient matrices" - the matrix `A_α` has 1 in positions where monomial α appears.

The dual is:
```
max  ⟨1, Z⟩  = trace(Z)
s.t. ⟨A_α, Z⟩ = c_α  ∀α
     Z ⪰ 0
```

In NCTSSoS, this is reformulated as:
```
max  b
s.t. c_α - δ_{α,1}·b = ⟨A_α, Z⟩  ∀α
     Z ⪰ 0
```

The `C_αj` dictionary precisely encodes `⟨A_α, Z⟩ = Σ_{i,j} C_αj[α,i,j] · Z[i,j]`.
