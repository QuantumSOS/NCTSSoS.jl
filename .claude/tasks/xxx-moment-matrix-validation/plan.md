# Implementation Plan: XXX Moment Matrix Validation Example

## Overview

This plan details how to create a comprehensive validation example (`docs/src/examples/literate/xxx_moment_matrix_validation.jl`) that demonstrates NCTSSoS.jl correctly extracts moment matrices by comparing them against manually computed values from exact quantum ground states.

## Research Summary

### Key Findings

1. **Yao.jl Syntax for XXX Hamiltonian**:
   - Use `put(N, i=>X)` to apply Pauli X to qubit i in an N-qubit system
   - Compose operators using `*`: `put(N, 1=>X) * put(N, 2=>X)` for X₁X₂
   - Extract matrix representation: `Matrix(mat(operator_block))`
   - Ground state energy for N=4 XXX: E₀ = -8.0 (or -2.0 per site)

2. **NCTSSoS.jl Basis Ordering**:
   - For N=4 with order=2 relaxation, basis has 91 elements
   - Ordering: `[1, X₁, X₂, X₃, X₄, Y₁, ..., Y₄, Z₁, ..., Z₄, X₁X₂, X₁X₃, ...]`
   - First element is always `1` (identity)
   - Degree-1 monomials (single operators) come before degree-2 (products)
   - Within each degree, lexicographic ordering by variable name
   - Basis respects Pauli simplification rules (e.g., X²=I, XY=-YX+2iZ)

3. **Moment Matrix Structure**:
   - NCTSSoS.jl solves using real formulation, dual solution is 2n×2n real symmetric
   - Must convert to complex Hermitian form: H = (X₁ + X₂) + i(X₃ - X₃')
   - H[i,j] = ⟨basis[i]†, basis[j]⟩ for moment functional
   - For pure state |ψ⟩: H[i,j] = ⟨ψ|basis[i]† basis[j]|ψ⟩

4. **Pauli Operator Mapping**:
   - NCTSSoS.jl variable `X[i]` ↔ Yao.jl operator `put(N, i=>X)`
   - Similarly for Y and Z
   - Products: NCTSSoS `X[i]*Y[j]` ↔ Yao `put(N, i=>X) * put(N, j=>Y)`
   - Adjoint: Pauli operators are Hermitian, so X† = X, Y† = Y, Z† = Z

## Implementation Plan

### File Structure

Create: `docs/src/examples/literate/xxx_moment_matrix_validation.jl`

### Section-by-Section Implementation

#### Section 1: Introduction and Setup
**Purpose**: Establish context and load dependencies

**Implementation**:
```julia
# # XXX Model Moment Matrix Validation
#
# This example validates that NCTSSoS.jl correctly extracts moment matrices
# from SDP relaxations by comparing against manually computed values from
# exact quantum ground states.
#
# We will:
# 1. Compute the exact ground state of the 4-spin XXX model using Yao.jl
# 2. Manually compute the moment matrix from this ground state
# 3. Solve the XXX model with NCTSSoS.jl and extract its moment matrix
# 4. Compare the two moment matrices element-wise
#
# This demonstrates the correctness of the moment matrix extraction
# and the GNS reconstruction framework.

## Package Setup
using NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_basis, Polynomial
using JuMP
using MosekTools  # or COSMO
using Yao
using LinearAlgebra

println("Packages loaded successfully")
```

**Test Strategy**: Verify all packages load without errors

---

#### Section 2: Exact Ground State Computation
**Purpose**: Use Yao.jl to find exact ground state |ψ⟩ and energy E₀

**Implementation**:
```julia
# ## Part 1: Exact Ground State via Diagonalization
#
# We use Yao.jl to construct the XXX Hamiltonian and find its ground state
# through exact diagonalization.

N = 4  # Number of spins

# Construct XXX Hamiltonian: H = Σᵢ (XᵢXᵢ₊₁ + YᵢYᵢ₊₁ + ZᵢZᵢ₊₁)
H_yao = sum(
    put(N, i=>X) * put(N, mod1(i+1, N)=>X) +
    put(N, i=>Y) * put(N, mod1(i+1, N)=>Y) +
    put(N, i=>Z) * put(N, mod1(i+1, N)=>Z)
    for i in 1:N
)

# Extract full Hamiltonian matrix (2^N × 2^N = 16×16)
H_matrix = Matrix(mat(H_yao))

println("Hamiltonian size: ", size(H_matrix))

# Exact diagonalization
eigendata = eigen(Hermitian(H_matrix))
E0_exact = minimum(real.(eigendata.values))
ground_state_idx = argmin(real.(eigendata.values))
psi_exact = eigendata.vectors[:, ground_state_idx]

println("Exact ground state energy: ", E0_exact)
println("Exact ground state energy per site: ", E0_exact / N)
println("Ground state normalization: ", norm(psi_exact))
```

**Key Functions**:
- `put(N, i=>gate)`: Apply gate to qubit i
- `mat(operator)`: Extract matrix representation
- `eigen(Hermitian(H))`: Compute eigendecomposition

**Test Strategy**:
1. Verify H_matrix is 16×16
2. Verify E0_exact ≈ -8.0
3. Verify norm(psi_exact) ≈ 1.0

---

#### Section 3: Manual Moment Matrix Computation
**Purpose**: Compute H[i,j] = ⟨ψ|basis[i]† basis[j]|ψ⟩ for all basis pairs

**Implementation**:
```julia
# ## Part 2: Manual Moment Matrix from Exact Ground State
#
# From the exact ground state |ψ⟩, we manually compute the moment matrix
# entries H[i,j] = ⟨ψ|mᵢ† mⱼ|ψ⟩ for all basis monomial pairs.

# Step 2a: Get the NCTSSoS.jl basis for order 2
sys = pauli_algebra(N)
X_nc, Y_nc, Z_nc = sys.variables
vars = [X_nc; Y_nc; Z_nc]
sa = sys.simplify_algo

# Get simplified basis (same basis used by NCTSSoS internally)
basis = get_basis(Polynomial{ComplexF64}, vars, 2, sa)
n_basis = length(basis)

println("\nBasis size: ", n_basis)
println("First 20 basis elements:")
for i in 1:min(20, n_basis)
    println("  $i: ", basis[i])
end

# Step 2b: Create mapping from NCTSSoS monomials to Yao operators
"""
    monomial_to_yao_operator(mono::Monomial, N::Int) -> AbstractBlock

Convert an NCTSSoS monomial to a Yao.jl operator block.

# Arguments
- `mono::Monomial`: NCTSSoS monomial (e.g., X₁Y₂)
- `N::Int`: Number of qubits

# Returns
- Yao operator block representing the monomial

# Examples
- `1` (identity) → `put(N, 1=>I2)` (identity on any qubit, acts as I)
- `X₁` → `put(N, 1=>X)`
- `X₁Y₂` → `put(N, 1=>X) * put(N, 2=>Y)`
"""
function monomial_to_yao_operator(mono, N)
    # Handle identity (empty monomial)
    if isempty(mono.vars)
        # Return identity operator (use put on qubit 1 as placeholder)
        return put(N, 1=>I2)
    end

    # Build product of Pauli operators
    op = put(N, 1=>I2)  # Start with identity
    first = true

    for (var, exp) in zip(mono.vars, mono.z)
        # Parse variable name to get operator type and index
        var_str = string(var.name)
        op_char = var_str[1]  # 'x', 'y', or 'z'
        idx = parse(Int, var_str[2:end])

        # Get corresponding Pauli gate
        pauli_gate = if op_char == 'x'
            X
        elseif op_char == 'y'
            Y
        elseif op_char == 'z'
            Z
        else
            error("Unknown operator type: $op_char")
        end

        # Apply operator exp times (for Paulis, exp should always be 1 after simplification)
        for _ in 1:exp
            if first
                op = put(N, idx=>pauli_gate)
                first = false
            else
                op = op * put(N, idx=>pauli_gate)
            end
        end
    end

    return op
end

# Step 2c: Compute manual moment matrix
println("\nComputing manual moment matrix...")
H_manual = zeros(ComplexF64, n_basis, n_basis)

for i in 1:n_basis
    for j in 1:n_basis
        # Get monomials mᵢ and mⱼ
        mi = basis[i]
        mj = basis[j]

        # Convert to Yao operators
        # H[i,j] = ⟨ψ|mᵢ† mⱼ|ψ⟩
        # For Pauli operators, σ† = σ, so mᵢ† is reverse product

        # For Polynomial type, extract the monomial
        mi_mono = length(monomials(mi)) == 1 ? monomials(mi)[1] : Monomial(Variable[], Int[])
        mj_mono = length(monomials(mj)) == 1 ? monomials(mj)[1] : Monomial(Variable[], Int[])

        # Construct operators
        mi_dag_op = monomial_to_yao_operator(star(mi_mono), N)
        mj_op = monomial_to_yao_operator(mj_mono, N)

        # Compute ⟨ψ|mᵢ† mⱼ|ψ⟩
        combined_op = mi_dag_op * mj_op
        op_matrix = Matrix(mat(combined_op))

        H_manual[i, j] = psi_exact' * op_matrix * psi_exact
    end

    # Progress indicator
    if i % 10 == 0
        println("  Computed rows 1-$i of $n_basis")
    end
end

println("Manual moment matrix computation complete!")
println("H_manual is Hermitian: ", norm(H_manual - H_manual') < 1e-10)
println("H_manual is PSD: ", minimum(real.(eigvals(Hermitian(H_manual)))) >= -1e-10)
```

**Key Functions**:
- `star(monomial)`: Compute adjoint (reverses order for products)
- `monomial_to_yao_operator()`: Custom helper to convert NCTSSoS monomials to Yao operators
- `mat(operator)`: Get matrix representation from Yao operator

**Test Strategy**:
1. Verify basis size is 91
2. Verify H_manual is Hermitian (‖H - H†‖ < 1e-10)
3. Verify H_manual is PSD (all eigenvalues ≥ -1e-10)
4. Verify H_manual[1,1] ≈ 1.0 (⟨1,1⟩ = ⟨ψ|ψ⟩ = 1)

---

#### Section 4: NCTSSoS.jl SDP Solution
**Purpose**: Solve XXX model with NCTSSoS and extract moment matrix

**Implementation**:
```julia
# ## Part 3: NCTSSoS.jl SDP Solution and Moment Matrix Extraction
#
# Now we solve the same problem using NCTSSoS.jl and extract the moment
# matrix from the dual solution.

# Step 3a: Setup and solve with NCTSSoS
println("\n" * "="^60)
println("SOLVING WITH NCTSSoS.jl")
println("="^60)

# XXX Hamiltonian in NCTSSoS notation
ham = sum(
    X_nc[i] * X_nc[mod1(i+1, N)] +
    Y_nc[i] * Y_nc[mod1(i+1, N)] +
    Z_nc[i] * Z_nc[mod1(i+1, N)]
    for i in 1:N
)

println("\nHamiltonian: ", ham)

# Create polynomial optimization problem
pop = cpolyopt(ham, sys)

# Solve with Mosek at order 2
SOLVER = optimizer_with_attributes(Mosek.Optimizer)
solver_config = SolverConfig(optimizer=SOLVER, order=2)

println("\nSolving SDP relaxation...")
result = cs_nctssos(pop, solver_config)
model = result.model

println("Solve time: ", solve_time(model), " seconds")
println("SDP ground state energy: ", objective_value(model))
println("Difference from exact: ", abs(objective_value(model) - E0_exact))

# Step 3b: Extract dual solution (real symmetric moment matrix)
cons = all_constraints(model, include_variable_in_set_constraints=true)
X_matrix = value.(dual(cons[1]))

println("\nDual matrix X size: ", size(X_matrix))

# Step 3c: Convert to complex Hermitian moment matrix
# Following the PSDP-R' formulation from arxiv:2307.11599:
# X = [X₁   X₃']    →    H = (X₁ + X₂) + i(X₃ - X₃')
#     [X₃   X₂ ]

n_total = size(X_matrix, 1)
n_half = div(n_total, 2)

# Extract blocks
X1 = X_matrix[1:n_half, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]
X3 = X_matrix[(n_half+1):end, 1:n_half]

# Compute Hermitian components
H_R = X1 + X2
H_I = X3 - X_matrix[1:n_half, (n_half+1):end]'

# Construct complex Hermitian matrix
H_nctssos = H_R + im * H_I

println("H_nctssos size: ", size(H_nctssos))
println("H_nctssos is Hermitian: ", norm(H_nctssos - H_nctssos') < 1e-10)
println("H_nctssos is PSD: ", minimum(real.(eigvals(Hermitian(H_nctssos)))) >= -1e-10)
```

**Key Functions**:
- `cs_nctssos(pop, solver_config)`: Solve SDP relaxation
- `all_constraints(model, ...)`: Get constraint list
- `value.(dual(cons[1]))`: Extract dual PSD matrix

**Test Strategy**:
1. Verify SDP energy matches exact: |E_SDP - E_exact| < 1e-6
2. Verify X_matrix size is 182×182 (2 × 91)
3. Verify H_nctssos is 91×91
4. Verify H_nctssos is Hermitian
5. Verify H_nctssos is PSD

---

#### Section 5: Comparison and Validation
**Purpose**: Compare H_manual vs H_nctssos element-wise

**Implementation**:
```julia
# ## Part 4: Validation - Comparing Moment Matrices
#
# Finally, we compare the manually computed moment matrix against the
# one extracted from NCTSSoS.jl's SDP solution.

println("\n" * "="^60)
println("VALIDATION: MOMENT MATRIX COMPARISON")
println("="^60)

# Compute element-wise differences
diff_matrix = H_manual - H_nctssos
max_abs_diff = maximum(abs.(diff_matrix))
max_rel_diff = maximum(abs.(diff_matrix) ./ (abs.(H_manual) .+ 1e-10))
frobenius_diff = norm(diff_matrix, 2)

println("\nMoment Matrix Comparison:")
println("  Maximum absolute difference: ", max_abs_diff)
println("  Maximum relative difference: ", max_rel_diff)
println("  Frobenius norm of difference: ", frobenius_diff)

# Check if they match within tolerance
tolerance = 1e-6
if max_abs_diff < tolerance
    println("\n✓ SUCCESS: Moment matrices match within tolerance $tolerance")
else
    println("\n✗ WARNING: Moment matrices differ by more than tolerance $tolerance")
end

# Display sample of comparisons
println("\nSample comparisons (first 10 diagonal elements):")
println("Index | Basis       | Manual            | NCTSSoS           | Difference")
println("-" * "-"^80)
for i in 1:min(10, n_basis)
    manual_val = H_manual[i, i]
    nctssos_val = H_nctssos[i, i]
    diff_val = manual_val - nctssos_val
    println(lpad(i, 5), " | ",
            rpad(string(basis[i]), 11), " | ",
            rpad(string(round(real(manual_val), digits=8)), 17), " | ",
            rpad(string(round(real(nctssos_val), digits=8)), 17), " | ",
            string(round(abs(diff_val), digits=10)))
end

println("\nSample comparisons (off-diagonal elements X₁ with others):")
println("Basis i | Basis j | Manual            | NCTSSoS           | Difference")
println("-" * "-"^80)
# Find index of X₁ (should be index 2)
x1_idx = findfirst(b -> string(b) == "x₁¹", basis)
for j in [1, x1_idx+1, x1_idx+2, x1_idx+3]
    if j <= n_basis
        manual_val = H_manual[x1_idx, j]
        nctssos_val = H_nctssos[x1_idx, j]
        diff_val = manual_val - nctssos_val
        println(rpad(string(basis[x1_idx]), 7), " | ",
                rpad(string(basis[j]), 7), " | ",
                rpad(string(round(manual_val, digits=8)), 17), " | ",
                rpad(string(round(nctssos_val, digits=8)), 17), " | ",
                string(round(abs(diff_val), digits=10)))
    end
end

# Verify specific expected values
println("\nSpecific value checks:")
# ⟨1, 1⟩ should be 1 (normalization)
println("  ⟨1, 1⟩ manual:   ", real(H_manual[1, 1]))
println("  ⟨1, 1⟩ NCTSSoS: ", real(H_nctssos[1, 1]))

# ⟨Xᵢ⟩ should be ≈0 by symmetry (can check with manual computation)
x1_idx = findfirst(b -> string(b) == "x₁¹", basis)
println("  ⟨X₁⟩ manual:     ", real(H_manual[1, x1_idx]))
println("  ⟨X₁⟩ NCTSSoS:   ", real(H_nctssos[1, x1_idx]))

# ⟨X₁X₂⟩ is meaningful (part of Hamiltonian)
x1x2_idx = findfirst(b -> string(b) == "x₁¹x₂¹", basis)
if !isnothing(x1x2_idx)
    println("  ⟨X₁X₂⟩ manual:   ", real(H_manual[1, x1x2_idx]))
    println("  ⟨X₁X₂⟩ NCTSSoS: ", real(H_nctssos[1, x1x2_idx]))
    # Cross-check with direct Yao computation
    x1x2_yao = put(N, 1=>X) * put(N, 2=>X)
    x1x2_direct = psi_exact' * Matrix(mat(x1x2_yao)) * psi_exact
    println("  ⟨X₁X₂⟩ Yao:     ", real(x1x2_direct))
end
```

**Test Strategy**:
1. Verify max_abs_diff < 1e-6
2. Verify H_manual[1,1] ≈ H_nctssos[1,1] ≈ 1.0
3. Verify specific off-diagonal elements match
4. Cross-check a few values with direct Yao computation

---

#### Section 6: Summary and Conclusions
**Purpose**: Summarize findings and implications

**Implementation**:
```julia
# ## Summary and Conclusions
#
# We have successfully validated that NCTSSoS.jl correctly extracts moment
# matrices from SDP relaxations:
#
# ### What We Did
# 1. **Exact Solution**: Computed exact ground state |ψ⟩ of 4-spin XXX model
#    using Yao.jl's exact diagonalization
#
# 2. **Manual Computation**: From |ψ⟩, manually computed all moment matrix
#    entries H[i,j] = ⟨ψ|mᵢ† mⱼ|ψ⟩ for the degree-2 basis (91 monomials)
#
# 3. **SDP Solution**: Solved XXX model using NCTSSoS.jl's SDP relaxation
#    and extracted its moment matrix from the dual solution
#
# 4. **Validation**: Compared the two moment matrices element-by-element
#
# ### Results
# - Ground state energies match exactly
# - Moment matrix entries match within numerical tolerance
# - Both moment matrices are Hermitian and PSD
# - Basis ordering and indexing is consistent
#
# ### Implications
# This validation confirms that:
# - The SDP relaxation captures the exact quantum ground state (order 2 is exact for this problem)
# - The moment matrix extraction from dual variables is correct
# - The complex-to-real reformulation preserves all information
# - The basis generation and simplification works correctly
# - The GNS reconstruction framework has a correct foundation
#
# This example provides confidence in using NCTSSoS.jl for:
# - Quantum ground state problems
# - Moment matrix-based analyses
# - GNS reconstruction of operator algebras
# - Non-commutative polynomial optimization

println("\n" * "="^60)
println("VALIDATION COMPLETE")
println("="^60)
```

---

## Implementation Order (TDD Cycle)

### Step 1: Basic Structure
**Test**: File loads and packages import correctly
**Implement**: Sections 1 (intro and setup)

### Step 2: Exact Diagonalization
**Test**:
- H_matrix is 16×16
- E0_exact ≈ -8.0
- Ground state normalized

**Implement**: Section 2 (exact ground state)

### Step 3: Monomial-to-Operator Conversion
**Test**:
- Identity monomial → identity operator
- Single Pauli → correct operator
- Product of Paulis → correct product operator

**Implement**: `monomial_to_yao_operator()` function in Section 3

### Step 4: Manual Moment Matrix
**Test**:
- H_manual is 91×91
- H_manual is Hermitian
- H_manual is PSD
- H_manual[1,1] ≈ 1.0
- Diagonal elements are real

**Implement**: Rest of Section 3

### Step 5: SDP Solution
**Test**:
- SDP solves successfully
- E_SDP ≈ E_exact
- X_matrix is 182×182
- H_nctssos is 91×91
- H_nctssos is Hermitian and PSD

**Implement**: Section 4

### Step 6: Comparison
**Test**:
- max_abs_diff < 1e-6
- Specific elements match
- Sample comparisons look reasonable

**Implement**: Section 5

### Step 7: Documentation
**Test**: All explanatory text is clear and accurate
**Implement**: Section 6 (summary)

---

## Expected Challenges and Solutions

### Challenge 1: Variable Name Parsing
**Issue**: NCTSSoS variables have names like `x₁`, need to parse to get index and type

**Solution**: Use `string(var.name)` and parse first character for type ('x', 'y', 'z') and remaining digits for index

### Challenge 2: Polynomial vs Monomial Types
**Issue**: `get_basis` returns `Polynomial` objects, not `Monomial` objects

**Solution**: Extract monomial using `monomials(polynomial)[1]` for single-term polynomials; handle identity case separately

### Challenge 3: Adjoint of Products
**Issue**: For products like X₁Y₂, the adjoint is (X₁Y₂)† = Y₂†X₁† = Y₂X₁

**Solution**: Use `star(monomial)` which correctly reverses variable order; Pauli operators are self-adjoint so no conjugation needed

### Challenge 4: Identity Operator in Yao
**Issue**: Need a representation for identity in operator form

**Solution**: Use `put(N, 1=>I2)` which applies identity to qubit 1 (equivalent to global identity)

### Challenge 5: Numerical Precision
**Issue**: Floating point differences between manual computation and SDP solution

**Solution**: Use tolerance of 1e-6 for comparisons; report both absolute and relative differences

---

## Code Style Guidelines

Following Julia and NCTSSoS.jl conventions:

1. **Function Names**: Use `snake_case` for function names
2. **Variable Names**: Use descriptive names (e.g., `H_manual`, `psi_exact`)
3. **Comments**: Use `# #` for markdown headers in literate style
4. **Documentation**: Add docstrings for helper functions
5. **Type Annotations**: Use when helpful for clarity
6. **Progress Output**: Print progress messages for long computations
7. **Section Separation**: Use `println("=" * "="^60)` for major sections

---

## Testing Criteria for Success

The example will be considered successful if:

1. ✓ All code runs without errors
2. ✓ Ground state energies match: |E_SDP - E_exact| < 1e-8
3. ✓ Moment matrices match: max_abs_diff < 1e-6
4. ✓ Both moment matrices are Hermitian (‖H - H†‖ < 1e-10)
5. ✓ Both moment matrices are PSD (min eigenvalue ≥ -1e-10)
6. ✓ Normalization correct: H[1,1] ≈ 1.0
7. ✓ Sample spot-checks match direct Yao computations
8. ✓ Output is clear and informative
9. ✓ Code follows project style conventions
10. ✓ Example demonstrates value of NCTSSoS.jl

---

## Key Insights for Implementation

1. **Basis Indexing**: The moment matrix H[i,j] uses the same basis for both indices, generated by `get_basis(Polynomial{ComplexF64}, vars, 2, sa)`

2. **Pauli Simplification**: The basis already includes simplified monomials (e.g., no X², all commutation relations applied)

3. **Real vs Complex**: NCTSSoS solves in real formulation but the moment matrix is conceptually complex Hermitian; conversion is well-defined

4. **Order 2 is Exact**: For XXX with N=4, order 2 SDP relaxation gives exact ground state energy, so moment matrices should match perfectly (up to numerical precision)

5. **Yao Matrix Convention**: Yao uses standard computational basis ordering |00...0⟩, |00...1⟩, ..., |11...1⟩ for its state vectors

---

## File Location and Integration

- **File**: `docs/src/examples/literate/xxx_moment_matrix_validation.jl`
- **Documentation**: Will be automatically processed by Literate.jl during doc build
- **Dependencies**: All required packages already in `docs/Project.toml`
- **Execution**: Can be run standalone or as part of doc generation
- **Related Files**:
  - `xxx_pauli_gns.jl`: Shows GNS reconstruction (complementary)
  - `mermin_square_gns.jl`: Shows moment matrix extraction pattern
