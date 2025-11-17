"""
Detailed Example: How Moment Matrix Constraints Become SOS Polynomial Constraints

This file demonstrates the conversion process step-by-step.
"""

using NCTSSoS, JuMP
using NCTSSoS: get_Cαj, constraint_object

# ============================================================================
# STEP 1: Create a Moment Matrix in the Primal Problem
# ============================================================================

println("="^70)
println("STEP 1: Create Moment Matrix (Primal)")
println("="^70)

# Create a simple primal model with moment variables
model = Model()
@variable(model, y[1:6])  # Moment variables: y_1, y_x, y_y, y_x², y_xy, y_y²
@constraint(model, y[1] == 1)  # Normalization

# Label the moments for clarity
moment_names = ["1", "x", "y", "x²", "xy", "y²"]
println("\nMoment variables (representing ⟨monomial⟩):")
for (i, name) in enumerate(moment_names)
    println("  y[$i] = ⟨$name⟩")
end

# ============================================================================
# STEP 2: Build Moment Matrix from JuMP Variables
# ============================================================================

println("\n" * "="^70)
println("STEP 2: Build Moment Matrix M with basis [1, x, y]")
println("="^70)

# Basis: [1, x, y]
# Moment matrix entries M[i,j] = ⟨b_i† · b_j⟩
# For real polynomials: b_i† = b_i

# M[1,1] = ⟨1·1⟩ = y[1]
# M[1,2] = ⟨1·x⟩ = y[2]
# M[1,3] = ⟨1·y⟩ = y[3]
# M[2,1] = ⟨x·1⟩ = y[2]
# M[2,2] = ⟨x·x⟩ = y[4]
# M[2,3] = ⟨x·y⟩ = y[5]
# M[3,1] = ⟨y·1⟩ = y[3]
# M[3,2] = ⟨y·x⟩ = y[5]
# M[3,3] = ⟨y·y⟩ = y[6]

M = [y[1] y[2] y[3];
     y[2] y[4] y[5];
     y[3] y[5] y[6]]

println("\nMoment matrix M (as JuMP expressions):")
println("M = ")
for i in 1:3
    print("    [ ")
    for j in 1:3
        print("$(M[i,j])  ")
    end
    println("]")
end

# Add PSD constraint
cons = @constraint(model, M in PSDCone())

println("\nConstraint added: M ⪰ 0")
println("Constraint type: ", typeof(constraint_object(cons)))

# ============================================================================
# STEP 3: Extract Coefficient Matrix C_αj using get_Cαj()
# ============================================================================

println("\n" * "="^70)
println("STEP 3: Extract Coefficient Matrix C_αj")
println("="^70)

# Create basis dictionary: JuMP variable → index
basis_dict = Dict(y[i] => i for i in 1:6)

println("\nBasis dictionary (JuMP variable → monomial index):")
for (i, name) in enumerate(moment_names)
    println("  y[$i] (⟨$name⟩) → index $i")
end

# Extract coefficients
C_αj = get_Cαj(basis_dict, constraint_object(cons))

println("\nCoefficient matrix C_αj:")
println("Maps: (monomial_index, row, col) → coefficient")
println("\nEntries in C_αj:")
for ((α_idx, i, j), coef) in sort(collect(C_αj), by=x->x[1])
    println("  C_αj[$(moment_names[α_idx]), $i, $j] = $coef")
end

# ============================================================================
# STEP 4: Understand the Meaning of C_αj
# ============================================================================

println("\n" * "="^70)
println("STEP 4: Interpret C_αj - Decomposing the Matrix")
println("="^70)

println("\nThe moment matrix can be written as:")
println("  M[i,j] = Σ_α C_αj[α, i, j] · y_α")
println("\nFor example:")
println("  M[1,1] = C_αj[1,1,1]·y[1] = 1·y[1]")
println("  M[1,2] = C_αj[2,1,2]·y[2] = 1·y[2]")
println("  M[2,2] = C_αj[4,2,2]·y[4] = 1·y[4]")

println("\nVerification - reconstruct M from C_αj:")
M_reconstructed = zeros(AffExpr, 3, 3)
for ((α_idx, i, j), coef) in C_αj
    M_reconstructed[i,j] += coef * y[α_idx]
end

println("M_reconstructed == M: ", M_reconstructed == M)

# ============================================================================
# STEP 5: Build Dual Problem - Create Dual Variables
# ============================================================================

println("\n" * "="^70)
println("STEP 5: Create Dual Variables G")
println("="^70)

dual_model = Model()

# For the primal constraint M ⪰ 0, create dual variable G ⪰ 0
@variable(dual_model, G[1:3, 1:3] in PSDCone())

println("\nDual variable G created:")
println("  G is a 3×3 PSD matrix variable")
println("  This corresponds to the primal constraint M ⪰ 0")

# Create bound variable
@variable(dual_model, b)
@objective(dual_model, Max, b)

println("\nObjective: maximize b")

# ============================================================================
# STEP 6: Build Polynomial Equality Constraints
# ============================================================================

println("\n" * "="^70)
println("STEP 6: Build Polynomial Equality Constraints")
println("="^70)

# Suppose objective is f = 2 - x² - y² = 2·y[1] - y[4] - y[6]
println("\nSuppose objective: f = 2 - x² - y²")
println("In terms of moments: f = 2·⟨1⟩ - ⟨x²⟩ - ⟨y²⟩ = 2·y[1] - y[4] - y[6]")

# Initialize polynomial constraints: f_α - δ_{α,1}·b
f_coefficients = [2.0, 0.0, 0.0, -1.0, 0.0, -1.0]  # [f_1, f_x, f_y, f_x², f_xy, f_y²]

fα_constraints = [AffExpr(f_coefficients[i]) for i in 1:6]
fα_constraints[1] -= b  # f_1 - b for identity monomial

println("\nInitial polynomial constraints (before adding dual contribution):")
for (i, name) in enumerate(moment_names)
    println("  f_α[$name] = $(fα_constraints[i])")
end

# Add contributions from dual variable: -Σ_{i,j} C_αj[α,i,j] · G[i,j]
println("\nAdding dual variable contributions: fα -= C_αj · G")
for ((α_idx, i, j), coef) in C_αj
    fα_constraints[α_idx] -= coef * G[i,j]
end

println("\nFinal polynomial equality constraints:")
for (i, name) in enumerate(moment_names)
    println("  Constraint for ⟨$name⟩: $(fα_constraints[i]) == 0")
end

# Add constraints to model
@constraint(dual_model, fα_constraints .== 0)

println("\nConstraints added to dual model!")

# ============================================================================
# STEP 7: Interpretation - What This Means
# ============================================================================

println("\n" * "="^70)
println("STEP 7: Interpretation")
println("="^70)

println("""
The polynomial equality constraints enforce:

  f(x) - b = G • C

where G • C means the Frobenius inner product:
  G • C = Σ_{i,j} G[i,j] · C[α,i,j]  for each monomial α

Breaking down the constraints:

1. For monomial ⟨1⟩:
   2 - b = 1·G[1,1] + 0·other_terms
   ⟹ We need G[1,1] = 2 - b

2. For monomial ⟨x⟩:
   0 = 1·G[1,2] + 1·G[2,1]
   ⟹ We need G[1,2] + G[2,1] = 0
   (But G is symmetric, so G[1,2] = G[2,1], thus G[1,2] = G[2,1] = 0)

3. For monomial ⟨x²⟩:
   -1 = 1·G[2,2]
   ⟹ We need G[2,2] = -1

   BUT WAIT! G ⪰ 0 means G[2,2] ≥ 0, so this is INFEASIBLE!
   This means the optimal b ≤ some value where feasibility is maintained.

The dual problem finds the maximum b such that we can find a PSD matrix G
satisfying these polynomial equalities, which proves f(x) - b ≥ 0.
""")

# ============================================================================
# STEP 8: Example with Localizing Matrix
# ============================================================================

println("="^70)
println("STEP 8: Localizing Matrix Example")
println("="^70)

println("\nFor inequality constraint g(x) ≥ 0, we create a localizing matrix:")
println("  M_g[i,j] = Σ_k (g)_k · y_{b_i† · m_k · b_j}")

println("\nExample: g(x,y) = 1 - x² - y² with basis [1, x]")
println("  g = 1·⟨1⟩ - 1·⟨x²⟩ - 1·⟨y²⟩")

# Localizing matrix for g with basis [1, x]
# M_g[1,1] = ⟨g·1·1⟩ = ⟨g⟩ = y[1] - y[4] - y[6]
# M_g[1,2] = ⟨g·1·x⟩ = ⟨g·x⟩ = y[2] - (would need ⟨x³⟩, ⟨xy²⟩...)
# Simplified for demonstration:

M_g = [y[1] - y[4]  y[2];
       y[2]         y[4]]

println("\nLocalizing matrix M_g (simplified):")
println("M_g = ")
for i in 1:2
    print("    [ ")
    for j in 1:2
        print("$(M_g[i,j])  ")
    end
    println("]")
end

cons_g = @constraint(model, M_g in PSDCone())

# Extract coefficients
basis_dict_simple = Dict(y[i] => i for i in 1:6)
C_αj_g = get_Cαj(basis_dict_simple, constraint_object(cons_g))

println("\nCoefficient matrix C_αj_g for localizing matrix:")
for ((α_idx, i, j), coef) in sort(collect(C_αj_g), by=x->x[1])
    println("  C_αj_g[$(moment_names[α_idx]), $i, $j] = $coef")
end

println("\nIn the dual, this creates another PSD variable G_g ⪰ 0")
println("and adds contributions to the polynomial constraints:")
println("  fα -= C_αj_g · G_g")

println("\n" * "="^70)
println("Summary: Moment Matrix → SOS Constraint Conversion")
println("="^70)
println("""
PRIMAL:
  - Variables: y_α (moments)
  - Constraints: M ⪰ 0 where M[i,j] = Σ_α C_αj[α,i,j] · y_α
  - Objective: min Σ_α f_α · y_α

DUAL:
  - Variables: G[i,j] (matrix) and b (scalar)
  - Constraints: f_α - δ_{α,1}·b = Σ_{i,j} C_αj[α,i,j] · G[i,j]  ∀α
                G ⪰ 0
  - Objective: max b

KEY FUNCTION: get_Cαj()
  - Input: Moment matrix M where entries are JuMP expressions
  - Output: Dictionary (α_idx, i, j) → coefficient
  - Meaning: Extracts which moment variable appears where with what coefficient

CONVERSION PROCESS:
  1. Primal PSD constraint M ⪰ 0 → Dual PSD variable G ⪰ 0
  2. Structure of M → Coefficient matrices C_αj via get_Cαj()
  3. Primal objective coefficients → Dual polynomial constraints
  4. Each monomial gets one polynomial equality constraint
  5. Dual constraints enforce: f(x) - b = sum of matrix products
""")
