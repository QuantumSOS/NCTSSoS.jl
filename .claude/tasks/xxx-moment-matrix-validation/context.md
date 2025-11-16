# Task: XXX Model Moment Matrix Validation Example

## Objective

Create a comprehensive example in `docs/src/examples/literate/` that validates NCTSSoS.jl's moment matrix extraction is correct for the XXX model with 4 spins.

## Requirements

The example must demonstrate:

1. **Exact Ground State Computation**
   - Use Yao.jl for exact diagonalization of the XXX Hamiltonian with 4 spins
   - Compute the exact ground state |ψ⟩ and ground state energy E₀

2. **Manual Moment Matrix Computation**
   - From the exact ground state |ψ⟩, manually compute expected moment matrix values
   - For basis monomials up to degree 2 (to match SDP relaxation order)
   - Compute ⟨m₁† m₂⟩ = ⟨ψ|m₁† m₂|ψ⟩ for all pairs of basis monomials

3. **NCTSSoS.jl Moment Matrix Extraction**
   - Use NCTSSoS.jl to solve the XXX model SDP relaxation (order 2)
   - Extract the moment matrix from the dual solution
   - Convert from real symmetric to complex Hermitian form

4. **Comparison and Validation**
   - Compare manual vs reconstructed moment matrix values
   - Report differences element-wise
   - Verify they match within numerical tolerance

## Project Context

### Existing Files
- `docs/src/examples/literate/xxx_pauli_gns.jl`: Existing XXX example that focuses on GNS reconstruction
- Other examples show various NCTSSoS.jl capabilities
- `src/gns.jl`: GNS reconstruction implementation

### XXX Model
- N=4 sites (4 spins)
- Hamiltonian: H = Σᵢ (Xᵢ Xᵢ₊₁ + Yᵢ Yᵢ₊₁ + Zᵢ Zᵢ₊₁) with periodic boundary conditions
- Uses Pauli operators satisfying: σᵢ² = I, {σᵢ, σⱼ} = 0, [σᵢ, σⱼ] = 2iεᵢⱼₖσₖ

### Available Dependencies
- Yao.jl: Already in docs/Project.toml, can be used for exact diagonalization
- NCTSSoS.jl: Main package
- JuMP, MosekTools/COSMO: SDP solvers
- LinearAlgebra: Standard library

## Key Technical Details

### Moment Matrix Extraction (from mermin_square_gns.jl)
```julia
# Extract dual solution X (real symmetric PSD matrix)
cons = all_constraints(model, include_variable_in_set_constraints=true)
X_matrix = value.(dual(cons[1]))

# Convert to complex Hermitian form
n_half = div(size(X_matrix, 1), 2)
X1 = X_matrix[1:n_half, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]
X3 = X_matrix[(n_half+1):end, 1:n_half]

H_R = X1 + X2
H_I = X3 - X_matrix[1:n_half, (n_half+1):end]'
H = H_R + im * H_I  # Complex Hermitian moment matrix
```

### Pauli Algebra Setup (from xxx_pauli_gns.jl)
```julia
N = 4
sys = pauli_algebra(N)
X, Y, Z = sys.variables
ham = sum(X[i] * X[mod1(i+1, N)] + Y[i] * Y[mod1(i+1, N)] + Z[i] * Z[mod1(i+1, N)] for i in 1:N)
pop = cpolyopt(ham, sys)
result = cs_nctssos(pop, solver_config)
```

## Deliverables

1. New file: `docs/src/examples/literate/xxx_moment_matrix_validation.jl`
2. The file should follow literate programming style with:
   - Markdown comments explaining each section (using `# #` for headers)
   - Clear separation of steps
   - Detailed verification output

## Success Criteria

- Manual and reconstructed moment matrices match within tolerance (e.g., 1e-6)
- Example runs successfully and produces clear validation output
- Code follows Julia best practices and existing example style
- Demonstrates NCTSSoS.jl correctly extracts moment matrices

---

## Research and Planning Summary (Sub-Agent Work)

### Key Research Findings

1. **Yao.jl Integration**:
   - Syntax: `put(N, i=>X)` for Pauli operators, compose with `*`
   - Ground state for N=4 XXX: E₀ = -8.0
   - Expectation values: ⟨ψ|O|ψ⟩ computed via `psi' * Matrix(mat(O)) * psi`

2. **NCTSSoS.jl Basis Structure**:
   - Order-2 relaxation produces 91 basis elements
   - Ordering: [1, X₁, ..., Z₄, X₁X₂, X₁X₃, ...] (lexicographic within each degree)
   - Basis respects Pauli algebra simplification

3. **Moment Matrix Extraction**:
   - Dual matrix X is 2n×2n real symmetric (n=91)
   - Convert to complex Hermitian: H = (X₁+X₂) + i(X₃-X₃')
   - H[i,j] = ⟨ψ|basis[i]† basis[j]|ψ⟩

### Implementation Architecture

The example is structured in 6 sections:
1. **Introduction & Setup**: Package loading and context
2. **Exact Diagonalization**: Yao.jl computes |ψ⟩ and E₀
3. **Manual Moment Matrix**: Compute H[i,j] from |ψ⟩ for all basis pairs
4. **NCTSSoS SDP Solution**: Solve and extract moment matrix from dual
5. **Validation**: Compare matrices element-wise
6. **Summary**: Document findings and implications

### Key Technical Components

**Helper Function**: `monomial_to_yao_operator(mono, N)`
- Converts NCTSSoS monomial (e.g., x₁¹y₂¹) to Yao operator (put(N,1=>X)*put(N,2=>Y))
- Handles identity, single operators, and products
- Uses variable name parsing to map x₁→X, y₁→Y, z₁→Z

**Validation Strategy**:
- Frobenius norm, max absolute/relative differences
- Sample diagonal and off-diagonal element comparisons
- Cross-check specific values with direct Yao computations

### Testing Approach (TDD)

1. Structure and imports
2. Exact diagonalization (verify E₀ ≈ -8.0, normalized |ψ⟩)
3. Monomial conversion function (unit tests for identity, single, product)
4. Manual moment matrix (Hermitian, PSD, H[1,1]=1)
5. SDP solution (E_SDP ≈ E_exact, proper matrix sizes)
6. Validation (max_abs_diff < 1e-6)
7. Documentation polish

### Expected Outcome

- Demonstrates NCTSSoS.jl correctly extracts moment matrices
- Validates SDP relaxation captures exact quantum ground state
- Provides confidence in GNS reconstruction framework
- Serves as comprehensive validation example for the package

**Status**: Plan complete. Ready for implementation by parent agent.
