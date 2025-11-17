# Complete Answer: Moment Matrix → SOS Polynomial Constraint Conversion

## Your Question

> "During converting `MomentProblem` to `SOSProblem` by `sos_dualize`, we have moment matrices in each clique and each block in terms of JuMP variables in `MomentProblem`'s `model`. How can we convert these matrices to constraints in the Sum of Square problem?"

## Short Answer

**The conversion happens through coefficient extraction:**

1. **Primal:** Moment matrices have entries that are **JuMP affine expressions** (linear combinations of moment variables y_α)

2. **Extract structure:** `get_Cαj()` extracts **which moment variable appears in which matrix position with which coefficient**

3. **Dual:** Each primal PSD constraint `M ⪰ 0` becomes a dual PSD variable `G ⪰ 0`

4. **Build constraints:** The extracted coefficients `C_αj` are used to build **polynomial equality constraints** relating the dual variables to the objective

5. **Result:** One polynomial equality constraint per unique monomial, encoding the sum-of-squares decomposition

---

## The Key Function: `get_Cαj()`

**Location:** `src/sos_solver.jl:8-22`

### What It Does

Takes a moment/localizing matrix constraint and extracts its linear structure:

**Input:** A matrix constraint where each entry `M[i,j]` is a JuMP affine expression:
```julia
M[i,j] = c₁·y_α₁ + c₂·y_α₂ + ... + cₖ·y_αₖ
```

**Output:** Dictionary mapping `(monomial_index, row, col) → coefficient`:
```julia
C_αj = {
    (α₁_idx, i, j) → c₁,
    (α₂_idx, i, j) → c₂,
    ...
}
```

**Meaning:** "Moment variable y_α appears in position M[i,j] with coefficient c"

### How It Works

```julia
function get_Cαj(basis_dict::Dict, localizing_mtx::VectorConstraint)
    dictionary_of_keys = Dict{Tuple{Int,Int,Int},T}()

    # Iterate over each matrix entry and its position
    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        # cur_expr is a JuMP AffExpr with .terms field
        # .terms yields (variable, coefficient) pairs

        for (α, coeff) in cur_expr.terms
            α_idx = basis_dict[α]  # Get monomial index
            row, col = ci.I[1], ci.I[2]  # Get matrix position

            # Record: "variable α appears at position (row,col) with coeff"
            dictionary_of_keys[(α_idx, row, col)] = coeff
        end
    end

    return dictionary_of_keys
end
```

### Concrete Example

**Moment matrix:**
```julia
M = [y[1]  y[2];     # [[⟨1⟩   ⟨x⟩  ]
     y[2]  y[3]]      #  [⟨x⟩   ⟨x²⟩]]
```

**After `get_Cαj()`:**
```julia
C_αj = {
    (1, 1, 1) → 1.0    # y[1] (⟨1⟩) appears at M[1,1] with coefficient 1
    (2, 1, 2) → 1.0    # y[2] (⟨x⟩) appears at M[1,2] with coefficient 1
    (2, 2, 1) → 1.0    # y[2] (⟨x⟩) appears at M[2,1] with coefficient 1
    (3, 2, 2) → 1.0    # y[3] (⟨x²⟩) appears at M[2,2] with coefficient 1
}
```

**Localizing matrix for g = 1 - x²:**
```julia
M_g = [y[1] - y[3]    y[2];
       y[2]           y[3]]
```

**After `get_Cαj()`:**
```julia
C_αj_g = {
    (1, 1, 1) →  1.0    # y[1] at M_g[1,1] with coef 1
    (3, 1, 1) → -1.0    # y[3] at M_g[1,1] with coef -1
    (2, 1, 2) →  1.0    # y[2] at M_g[1,2] with coef 1
    (2, 2, 1) →  1.0    # y[2] at M_g[2,1] with coef 1
    (3, 2, 2) →  1.0    # y[3] at M_g[2,2] with coef 1
}
```

---

## The Dualization Process (Step-by-Step)

### Step 1: Start with Primal Moment Problem

**Variables:** `y[α]` representing moments `⟨α⟩`

**Constraints:**
- Moment matrices: `M ⪰ 0` where `M[i,j] = Σ_α C_αj[α,i,j] · y[α]`
- Localizing matrices: `M_g ⪰ 0` for inequality constraints
- Equality matrices: `M_h = 0` for equality constraints

**Objective:** `min Σ_α f_α · y[α]`

### Step 2: Create Dual Variables

**Code:** `src/sos_solver.jl:54-56`

For each primal constraint:
```julia
dual_variables[j] = @variable(dual_model, [1:dim, 1:dim] in
    (constraint_is_equality ? SymmetricMatrixSpace() : PSDCone()))
```

- Primal: `M ⪰ 0` → Dual: `G ⪰ 0`
- Primal: `M = 0` → Dual: `G` symmetric (unconstrained sign)

### Step 3: Initialize Polynomial Constraints

**Code:** `src/sos_solver.jl:73-79`

One constraint per unique monomial:
```julia
fα_constraints[α] = f_α  # Objective coefficient for monomial α
fα_constraints[1] -= b   # Subtract bound variable from identity monomial
```

### Step 4: Extract Coefficients and Build Constraints

**Code:** `src/sos_solver.jl:81-86`

```julia
for (j, constraint) in enumerate(primal_constraints):
    # Extract coefficient structure
    C_αj = get_Cαj(basis_dict, constraint)

    # Add dual variable contributions to polynomial constraints
    for ((α_idx, i, k), coef) in C_αj:
        fα_constraints[α_idx] -= coef * dual_variables[j][i,k]
    end
end
```

**What this does:**

For each primal constraint with dual variable `G_j`, add to polynomial constraints:
```
fα_constraints[α] -= Σ_{i,k} C_αj[α, i, k] · G_j[i,k]
```

**Example:**

With `C_αj = {(1,1,1)→1, (2,1,2)→1, (3,2,2)→1}` and dual variable `G`:

```julia
fα_constraints[1] -= 1.0 * G[1,1]  # From (1,1,1)→1
fα_constraints[2] -= 1.0 * G[1,2]  # From (2,1,2)→1
fα_constraints[3] -= 1.0 * G[2,2]  # From (3,2,2)→1
```

### Step 5: Enforce Polynomial Equalities

**Code:** `src/sos_solver.jl:87`

```julia
@constraint(dual_model, fα_constraints .== 0)
```

**Final constraints** for each monomial α:
```
f_α - δ_{α,1}·b - Σ_j Σ_{i,k} C_αj[α,i,k]·G_j[i,k] = 0
```

Or rearranged:
```
f_α - δ_{α,1}·b = Σ_j (G_j • C_αj)
```

where `•` is the Frobenius inner product.

---

## Complete Example Walkthrough

### Problem Setup

**Primal:**
```
minimize  2 - x²

subject to:  1 - x² ≥ 0
             (moment relaxation order 1)
```

**Variables:** `y = [y_1, y_x, y_x²]` for `[⟨1⟩, ⟨x⟩, ⟨x²⟩]`

**Constraints:**

Moment matrix (basis `[1, x]`):
```
M = [y_1   y_x ]  ⪰ 0
    [y_x   y_x²]
```

Localizing matrix for `g = 1 - x²`:
```
M_g = [y_1 - y_x²   y_x - y_x³]  ⪰ 0
      [y_x - y_x³   y_x² - y_x⁴]
```

**Objective:** `2·y_1 - y_x² = 2 - x²`

### Coefficient Extraction

**From M:**
```julia
get_Cαj(M) = {
    (1, 1, 1) → 1.0,   # y_1 at M[1,1]
    (2, 1, 2) → 1.0,   # y_x at M[1,2]
    (2, 2, 1) → 1.0,   # y_x at M[2,1]
    (3, 2, 2) → 1.0    # y_x² at M[2,2]
}
```

**From M_g (first row only for brevity):**
```julia
get_Cαj(M_g) includes {
    (1, 1, 1) →  1.0,   # y_1 at M_g[1,1]
    (3, 1, 1) → -1.0,   # y_x² at M_g[1,1] (negative!)
    (2, 1, 2) →  1.0,   # y_x at M_g[1,2]
    ...
}
```

### Dual Problem Construction

**Variables:**
```julia
G_0[1:2, 1:2] ⪰ 0      # Dual of M ⪰ 0
G_1[1:2, 1:2] ⪰ 0      # Dual of M_g ⪰ 0
b                       # Bound variable
```

**Initialize polynomial constraints:**
```julia
fα_constraints[1] = 2 - b      # For ⟨1⟩
fα_constraints[2] = 0          # For ⟨x⟩
fα_constraints[3] = -1         # For ⟨x²⟩
```

**Add G_0 contributions:**
```julia
fα_constraints[1] -= 1.0 * G_0[1,1]  # (1,1,1)→1
fα_constraints[2] -= 1.0 * G_0[1,2]  # (2,1,2)→1
fα_constraints[2] -= 1.0 * G_0[2,1]  # (2,2,1)→1
fα_constraints[3] -= 1.0 * G_0[2,2]  # (3,2,2)→1
```

**Add G_1 contributions:**
```julia
fα_constraints[1] -= 1.0 * G_1[1,1]   # (1,1,1)→1
fα_constraints[3] -= (-1.0) * G_1[1,1]  # (3,1,1)→-1, so +=
fα_constraints[2] -= 1.0 * G_1[1,2]   # (2,1,2)→1
...
```

**Final constraints:**
```julia
Constraint for ⟨1⟩:
  2 - b - G_0[1,1] - G_1[1,1] = 0

Constraint for ⟨x⟩:
  0 - G_0[1,2] - G_0[2,1] - G_1[1,2] - G_1[2,1] = 0

Constraint for ⟨x²⟩:
  -1 - G_0[2,2] + G_1[1,1] - ... = 0
```

**Dual optimization problem:**
```julia
maximize  b

subject to:
  G_0[1,1] + G_1[1,1] = 2 - b
  G_0[1,2] + G_0[2,1] + G_1[1,2] + G_1[2,1] = 0
  G_0[2,2] - G_1[1,1] + ... = -1
  (additional constraints...)
  G_0 ⪰ 0
  G_1 ⪰ 0
```

---

## Why This Works: Mathematical Intuition

### Primal-Dual Relationship

The primal problem can be written as:
```
min  c'y
s.t. A(y) ⪰ 0
```

where `A(y) = Σ_α A_α · y_α` is a linear matrix function.

The `A_α` matrices are:
- `A_α[i,j] = 1` if moment α appears at position (i,j)
- `A_α[i,j] = 0` otherwise

The dual problem is:
```
max  b
s.t. ⟨A_α, Z⟩ = c_α - δ_{α,1}·b  ∀α
     Z ⪰ 0
```

where `⟨A_α, Z⟩ = Σ_{i,j} A_α[i,j] · Z[i,j]` is the Frobenius inner product.

**The C_αj dictionary precisely encodes the nonzero entries of A_α!**

```
C_αj[α, i, j] = coeff  ⟺  A_α[i,j] = coeff
```

### Sum-of-Squares Interpretation

The dual constraints enforce:
```
f(x) - b = Σ_j (σ_j(x)† · G_j · σ_j(x))
```

where `σ_j(x)` are vectors of monomials (the basis).

Since `G_j ⪰ 0`, each term `σ_j(x)† · G_j · σ_j(x) ≥ 0`, proving `f(x) ≥ b`.

The dual maximizes `b` to find the tightest lower bound provable by SOS.

---

## Summary Table

| Aspect | Primal (Moment) | Conversion | Dual (SOS) |
|--------|-----------------|------------|------------|
| **Variables** | `y_α` (moments) | → | `G_j[i,k]` (matrices), `b` (scalar) |
| **Constraint** | `M ⪰ 0` where `M[i,j] = Σ_α C_αj[α,i,j]·y_α` | `get_Cαj()` extracts `C_αj` | Polynomial equalities |
| **Structure** | PSD matrix constraint on moments | Coefficient extraction | Linear constraints on matrix variables |
| **Objective** | Minimize `Σ_α f_α · y_α` | → | Maximize `b` |
| **Final form** | Linear program over moments with PSD constraints | → | Linear program over matrices with polynomial constraints |

---

## Key Insight

**The moment matrix structure becomes the polynomial constraint coefficients.**

- **Primal:** Which moments appear where in the matrix determines the structure
- **Dual:** This structure becomes the coefficients in the polynomial equalities
- **Bridge:** `get_Cαj()` extracts this structure from JuMP expressions

The conversion is **exact** - both problems have the same optimal value by SDP duality!

---

## Files for Further Study

1. **MOMENT_TO_SOS_CONVERSION.md** - Detailed mathematical explanation
2. **CONVERSION_CODE_FLOW.md** - Visual code flow diagram
3. **src/sos_solver.jl:8-22** - `get_Cαj()` implementation
4. **src/sos_solver.jl:50-90** - `sos_dualize()` full implementation
5. **test/sos_solver.jl:66-92** - Test showing `get_Cαj()` in action
