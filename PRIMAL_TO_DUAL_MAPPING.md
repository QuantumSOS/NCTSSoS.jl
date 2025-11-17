# Mapping Primal JuMP Variables to Dual Constraints

## The Question

How can I map JuMP variables in `MomentProblem.model` (primal) to constraint indices in `SOSProblem.model` (dual)?

## Short Answer

The mapping is established through the **`symmetric_basis`** in `sos_solver.jl:68-89`:

```julia
# Primal variable at index i:
primal_var = symmetric_variables[i]  # Or: moment_problem.monomap[symmetric_basis[i]]

# Corresponding dual constraint at index i:
dual_constraint = dual_model[:coef_cons][i]

# The monomial:
monomial = symmetric_basis[i]
```

---

## Detailed Explanation

### Step 1: Build the Symmetric Basis

**Location:** `sos_solver.jl:66-71`

```julia
# Get all monomials from the primal problem's monomap
unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))

# Canonicalize and remove duplicates (handles α and α† mapping to same moment)
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(moment_problem.sa)))

# Get the JuMP variables corresponding to each unique monomial
symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)
```

**What this creates:**

- `symmetric_basis`: Vector of **unique monomials** (after canonicalization)
- `symmetric_variables`: Vector of **primal JuMP variables** corresponding to each monomial

**Example:**
```julia
# If monomap = {1 → y[1], x → y[2], x† → y[2], x² → y[3]}
# (Note: x and x† both map to y[2] for Hermitian symmetry)

unsymmetrized_basis = [1, x, x†, x²]

symmetric_basis = [1, x, x²]  # x† canonicalized to x, duplicate removed

symmetric_variables = [y[1], y[2], y[3]]  # The JuMP variables
```

### Step 2: Create Polynomial Constraints

**Location:** `sos_solver.jl:74`

```julia
fα_constraints = [GenericAffExpr{T,VariableRef}(
                      get(primal_objective_terms, α, zero(T)))
                  for α in symmetric_variables]
```

**What this creates:**

- One constraint expression per unique monomial
- Indexed by position in `symmetric_basis`
- `fα_constraints[i]` corresponds to `symmetric_basis[i]` and `symmetric_variables[i]`

### Step 3: Store in Dual Model

**Location:** `sos_solver.jl:89`

```julia
dual_model[:coef_cons] = @constraint(dual_model, fα_constraints .== 0)
```

**What this creates:**

- Stores the polynomial equality constraints with name `:coef_cons`
- Can be accessed as: `dual_model[:coef_cons]` (returns a vector of constraint references)
- `dual_model[:coef_cons][i]` is the constraint for `symmetric_basis[i]`

---

## The Complete Mapping

### Forward Mapping: Primal Variable → Dual Constraint

**Given:** A primal JuMP variable `y_var` from `moment_problem.model`

**Find:** The corresponding dual constraint index

**Method 1: Using the monomap (if you know the monomial)**

```julia
# If you know the monomial corresponding to y_var:
monomial = find_monomial_for_variable(y_var, moment_problem.monomap)  # Reverse lookup
canonical_monomial = canonicalize(monomial, moment_problem.sa)
constraint_idx = searchsortedfirst(symmetric_basis, canonical_monomial)
dual_constraint = dual_model[:coef_cons][constraint_idx]
```

**Method 2: Direct search in symmetric_variables**

```julia
# Find where y_var appears in symmetric_variables
constraint_idx = findfirst(==(y_var), symmetric_variables)
dual_constraint = dual_model[:coef_cons][constraint_idx]
monomial = symmetric_basis[constraint_idx]
```

### Reverse Mapping: Dual Constraint Index → Primal Variable

**Given:** A dual constraint index `i`

**Find:** The corresponding primal variable and monomial

```julia
# Direct indexing:
monomial = symmetric_basis[i]
primal_var = symmetric_variables[i]  # Or: moment_problem.monomap[monomial]
dual_constraint = dual_model[:coef_cons][i]
```

---

## Helper Dictionaries in the Code

### 1. `symmetrized_α2cons_dict`

**Location:** `sos_solver.jl:76`

```julia
symmetrized_α2cons_dict = Dict(
    zip(unsymmetrized_basis,
        map(x -> searchsortedfirst(symmetric_basis, canonicalize(x, moment_problem.sa)),
            unsymmetrized_basis))
)
```

**Purpose:** Maps **any monomial** (before canonicalization) to its constraint index

**Usage:**
```julia
# Given any monomial (possibly non-canonical):
monomial = x†  # Or any monomial from unsymmetrized_basis
constraint_idx = symmetrized_α2cons_dict[monomial]
dual_constraint = dual_model[:coef_cons][constraint_idx]
```

**Example:**
```julia
# If symmetric_basis = [1, x, x²]
# And x† canonicalizes to x (which is at index 2)

symmetrized_α2cons_dict = {
    1 → 1,    # Identity at index 1
    x → 2,    # x at index 2
    x† → 2,   # x† also maps to index 2 (same as x)
    x² → 3    # x² at index 3
}
```

### 2. `unsymmetrized_basis_vals_dict`

**Location:** `sos_solver.jl:78`

```julia
unsymmetrized_basis_vals_dict = Dict(
    zip(getindex.(Ref(moment_problem.monomap), unsymmetrized_basis),
        1:length(unsymmetrized_basis))
)
```

**Purpose:** Maps **primal JuMP variables** to their index in `unsymmetrized_basis`

**Usage:**
```julia
# Given a primal JuMP variable:
y_var = y[2]
unsymmetrized_idx = unsymmetrized_basis_vals_dict[y_var]
monomial = unsymmetrized_basis[unsymmetrized_idx]
```

**Note:** This is used internally in `get_Cαj()` to extract coefficient structure.

---

## Complete Example Code

### Setup

```julia
using NCTSSoS, JuMP

# Create a simple polynomial optimization problem
@ncpolyvar x
f = 2 - x^2
g = 1 - x^2

pop = polyopt(f; ineq_constraints=[g])

# Create moment problem (primal)
corr_sparsity = compute_correlative_sparsity(pop, ...)
term_sparsities = compute_term_sparsity(pop, corr_sparsity, ...)
moment_problem = moment_relax(pop, corr_sparsity, term_sparsities)

# Dualize to SOS problem
sos_problem = sos_dualize(moment_problem)
dual_model = sos_problem.model
```

### Access the Mapping Data Structures

Since `symmetric_basis` and `symmetric_variables` are local to `sos_dualize()`, you need to reconstruct them:

```julia
# Reconstruct the mapping (same logic as in sos_dualize)
sa = moment_problem.sa
unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

println("Symmetric basis (unique monomials):")
for (i, mono) in enumerate(symmetric_basis)
    println("  [$i] $mono")
end

println("\nPrimal variables (moment variables):")
for (i, var) in enumerate(symmetric_variables)
    println("  [$i] $var (represents ⟨$(symmetric_basis[i])⟩)")
end

println("\nDual constraints:")
for (i, cons) in enumerate(dual_model[:coef_cons])
    println("  [$i] $cons (polynomial equality for ⟨$(symmetric_basis[i])⟩)")
end
```

### Map Primal Variable to Dual Constraint

```julia
# Example: Find dual constraint for a specific primal variable
y_x² = moment_problem.monomap[x^2]  # Get the variable for ⟨x²⟩

# Method 1: Search in symmetric_variables
constraint_idx = findfirst(==(y_x²), symmetric_variables)
dual_constraint = dual_model[:coef_cons][constraint_idx]

println("Primal variable: $y_x²")
println("Represents moment: ⟨$(symmetric_basis[constraint_idx])⟩")
println("Dual constraint index: $constraint_idx")
println("Dual constraint: $dual_constraint")
```

### Map Monomial to Dual Constraint

```julia
# Given a monomial, find its dual constraint
monomial = x^2

# Canonicalize and find index
canonical_mono = canonicalize(monomial, sa)
constraint_idx = searchsortedfirst(symmetric_basis, canonical_mono)

# Get corresponding constraint and variable
dual_constraint = dual_model[:coef_cons][constraint_idx]
primal_var = symmetric_variables[constraint_idx]

println("Monomial: $monomial")
println("Canonical form: $canonical_mono")
println("Constraint index: $constraint_idx")
println("Primal variable: $primal_var")
println("Dual constraint: $dual_constraint")
```

### Iterate Over All Mappings

```julia
println("\n" * "="^70)
println("Complete Primal-Dual Mapping")
println("="^70)

for (i, (mono, pvar)) in enumerate(zip(symmetric_basis, symmetric_variables))
    dcons = dual_model[:coef_cons][i]
    println("\nIndex $i:")
    println("  Monomial:        $mono")
    println("  Primal variable: $pvar")
    println("  Dual constraint: $dcons")
end
```

---

## Programmatic Access: Building Mapping Functions

### Function 1: Get Dual Constraint Index from Primal Variable

```julia
function get_dual_constraint_index(
    primal_var::VariableRef,
    moment_problem::MomentProblem
)
    # Reconstruct symmetric_variables
    sa = moment_problem.sa
    unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
    symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

    # Find the variable
    idx = findfirst(==(primal_var), symmetric_variables)

    if idx === nothing
        error("Variable $primal_var not found in symmetric variables")
    end

    return idx
end

# Usage:
# idx = get_dual_constraint_index(y[2], moment_problem)
# dual_constraint = dual_model[:coef_cons][idx]
```

### Function 2: Get Primal Variable from Monomial

```julia
function get_primal_variable_for_monomial(
    monomial::M,
    moment_problem::MomentProblem
) where {M}
    canonical_mono = canonicalize(monomial, moment_problem.sa)
    return get(moment_problem.monomap, canonical_mono, nothing)
end

# Usage:
# primal_var = get_primal_variable_for_monomial(x^2, moment_problem)
```

### Function 3: Build Complete Mapping Dictionary

```julia
function build_primal_dual_mapping(
    moment_problem::MomentProblem,
    sos_problem::SOSProblem
)
    sa = moment_problem.sa
    unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
    symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

    dual_model = sos_problem.model

    mapping = Dict{VariableRef, NamedTuple{(:constraint_index, :monomial, :dual_constraint), Tuple{Int, Any, ConstraintRef}}}()

    for (i, (mono, pvar)) in enumerate(zip(symmetric_basis, symmetric_variables))
        mapping[pvar] = (
            constraint_index = i,
            monomial = mono,
            dual_constraint = dual_model[:coef_cons][i]
        )
    end

    return mapping
end

# Usage:
# mapping = build_primal_dual_mapping(moment_problem, sos_problem)
# info = mapping[y[2]]
# println("Constraint index: $(info.constraint_index)")
# println("Monomial: $(info.monomial)")
# println("Dual constraint: $(info.dual_constraint)")
```

---

## Key Insights

### 1. One-to-One Correspondence

There is a **one-to-one correspondence** between:
- Unique monomials (after canonicalization)
- Primal moment variables (for those monomials)
- Dual polynomial equality constraints

All indexed by the same position `i` in their respective arrays.

### 2. Hermitian Symmetry Handling

For Hermitian systems:
- Monomials `α` and `α†` may map to the **same** primal variable
- `canonicalize()` ensures they both map to the **same** dual constraint
- This is why `symmetric_basis` has fewer elements than `unsymmetrized_basis`

### 3. Constraint Naming

The dual constraints are stored with name `:coef_cons` in the dual model:
```julia
dual_model[:coef_cons]  # Vector of constraint references
```

This allows easy access after dualization.

### 4. The Mapping is Implicit

The code doesn't explicitly store a dictionary mapping primal variables to dual constraints. Instead:
- The **ordering** in `symmetric_basis` establishes the mapping
- Both primal variables and dual constraints use the **same indexing**

---

## Summary Table

| Index `i` | Monomial | Primal Variable | Dual Constraint | Constraint Meaning |
|-----------|----------|-----------------|-----------------|-------------------|
| 1 | `symmetric_basis[1]` | `symmetric_variables[1]` | `dual_model[:coef_cons][1]` | Polynomial equality for this monomial |
| 2 | `symmetric_basis[2]` | `symmetric_variables[2]` | `dual_model[:coef_cons][2]` | ... |
| ... | ... | ... | ... | ... |
| n | `symmetric_basis[n]` | `symmetric_variables[n]` | `dual_model[:coef_cons][n]` | ... |

---

## Code Locations Reference

| Component | File | Line | Description |
|-----------|------|------|-------------|
| Build `unsymmetrized_basis` | `src/sos_solver.jl` | 66 | Extract all monomials from monomap |
| Build `symmetric_basis` | `src/sos_solver.jl` | 68 | Canonicalize and remove duplicates |
| Build `symmetric_variables` | `src/sos_solver.jl` | 71 | Get primal JuMP variables |
| Create `fα_constraints` | `src/sos_solver.jl` | 74 | Initialize polynomial constraints |
| Store dual constraints | `src/sos_solver.jl` | 89 | Save as `:coef_cons` in dual model |
| Helper: `symmetrized_α2cons_dict` | `src/sos_solver.jl` | 76 | Map any monomial to constraint index |
| Helper: `unsymmetrized_basis_vals_dict` | `src/sos_solver.jl` | 78 | Map primal variable to basis index |

---

## Practical Usage Pattern

When working with primal and dual problems:

```julia
# 1. Create the problems
moment_problem = moment_relax(...)
sos_problem = sos_dualize(moment_problem)

# 2. Reconstruct the mapping basis
sa = moment_problem.sa
unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))
symmetric_basis = sorted_unique(canonicalize.(unsymmetrized_basis, Ref(sa)))
symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

# 3. Access mappings by index
for i in 1:length(symmetric_basis)
    monomial = symmetric_basis[i]
    primal_var = symmetric_variables[i]
    dual_constraint = sos_problem.model[:coef_cons][i]

    # Work with the mapping...
end
```

This gives you full control over the primal-dual correspondence!
