# Moment Matrix Extraction Implementation Plan

## Executive Summary

This plan outlines the implementation of a feature to extract moment matrices from `PolyOptResult` models. The extraction method differs based on whether the model originated from a Moment problem (primal) or SOS problem (dual). The implementation will provide a clean API for users to access moment matrices needed for GNS reconstruction without manually constructing them.

## 1. Problem Analysis

### 1.1 Current State
- `PolyOptResult` only stores `model::GenericModel{T}`, not the full `MomentProblem` or `SOSProblem`
- Users need moment matrices for GNS reconstruction but must currently construct them manually
- The model type (Moment vs SOS) determines extraction method but isn't tracked

### 1.2 Key Challenge
Since `PolyOptResult` doesn't store whether the problem was dualized, we need to distinguish between:
- **Moment Problem**: Use `dual()` on PSD constraints to get moment matrix values
- **SOS Problem**: Use `dual()` on equality constraints to get coefficient values

### 1.3 Structural Information Available
The `PolyOptResult` already contains rich structural information:
- `corr_sparsity`: Defines cliques and constraint assignment
- `cliques_term_sparsities`: For each clique, contains `TermSparsity` objects with:
  - First element: moment matrix basis blocks
  - Remaining elements: localizing matrix basis blocks for each constraint

## 2. Distinguishing Moment vs SOS Models

### 2.1 Recommendation: Add Metadata Field

**Preferred Approach**: Add a `problem_type` field to `PolyOptResult`

**Rationale**:
1. **Reliability**: Model inspection is fragile and could break with JuMP changes
2. **Clarity**: Explicit metadata makes intent clear and code maintainable
3. **Performance**: No runtime inspection overhead
4. **Extensibility**: Easier to add new problem types in future

**Implementation**:
```julia
# src/interface.jl
@enum ProblemType MomentPrimal SOSDual

struct PolyOptResult{T,P,M}
    objective::T
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
    problem_type::ProblemType  # NEW FIELD
end
```

**Changes Required**:
- Modify `cs_nctssos()` to pass `dualize ? SOSDual : MomentPrimal`
- Modify `cs_nctssos_higher()` similarly
- Update `PolyOptResult` constructor calls

### 2.2 Alternative: Model Inspection (Not Recommended)

Could inspect model structure:
- Moment models: Have vector of PSD constraint refs stored separately
- SOS models: Have equality constraints with matrix variables

This is fragile and not recommended, but included for completeness:
```julia
function is_sos_model(model::GenericModel)
    # Check for equality constraints with matrix variables
    # This is fragile and version-dependent
    cons = all_constraints(model, include_variable_in_set_constraints=false)
    # Look for matrix variable patterns...
end
```

## 3. API Design

### 3.1 Primary Function: `moment_matrix()`

```julia
"""
    moment_matrix(result::PolyOptResult{T,P,M}; clique::Int=1) where {T,P,M}

Extract the moment matrix from a polynomial optimization result for a specific clique.

The moment matrix is a Hankel matrix indexed by a monomial basis, where entry (i,j)
corresponds to the moment ⟨basis[i] * basis[j]⟩. This matrix is used as input to
GNS reconstruction for extracting matrix representations of variables.

# Arguments
- `result::PolyOptResult`: Result from `cs_nctssos` containing the solved model
- `clique::Int=1`: Index of the clique (default: 1 for single-clique problems)

# Returns
- `MomentMatrix{T,M}`: Structure containing:
  - `matrix::Matrix{T}`: The moment matrix values
  - `basis::Vector{M}`: Monomial basis indexing the matrix
  - `clique_vars::Vector{Variable}`: Variables in this clique

# Example
```julia
@ncpolyvar x[1:2]
f = x[1]^2 + x[2]^2 - 1
pop = polyopt(f)
result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))

# Extract moment matrix
mm = moment_matrix(result)

# Use with GNS reconstruction
X_mats = reconstruct(mm.matrix, mm.clique_vars, 2)
```

# Notes
- For problems with term sparsity, the moment matrix may be block-structured
- The basis is derived from `result.cliques_term_sparsities[clique][1]`
- Extraction method depends on whether the problem was dualized (stored in `result.problem_type`)
"""
function moment_matrix(result::PolyOptResult{T,P,M}; clique::Int=1) where {T,P,M}
end
```

### 3.2 Return Type: `MomentMatrix` Struct

```julia
# src/gns.jl (or new src/moment_extraction.jl)
"""
Structure containing a moment matrix and its associated basis information.

# Fields
- `matrix::Matrix{T}`: The moment matrix values, where entry (i,j) = ⟨basis[i] * basis[j]⟩
- `basis::Vector{M}`: Monomial basis indexing rows and columns of the matrix
- `clique_vars::Vector{Variable}`: Variables appearing in this clique
- `clique_index::Int`: Index of the clique in the original problem
"""
struct MomentMatrix{T,M}
    matrix::Matrix{T}
    basis::Vector{M}
    clique_vars::Vector{Variable}
    clique_index::Int
end

function Base.show(io::IO, mm::MomentMatrix)
    println(io, "MomentMatrix for clique $(mm.clique_index)")
    println(io, "  Variables: $(mm.clique_vars)")
    println(io, "  Basis size: $(length(mm.basis))")
    println(io, "  Matrix size: $(size(mm.matrix))")
end
```

### 3.3 Multi-Clique Support

```julia
"""
    moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}

Extract moment matrices for all cliques in the optimization result.

# Returns
- `Vector{MomentMatrix{T,M}}`: One moment matrix per clique

# Example
```julia
# For a problem with multiple cliques
mms = moment_matrices(result)
for (i, mm) in enumerate(mms)
    println("Clique $i: variables ", mm.clique_vars)
end
```
"""
function moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}
    return [moment_matrix(result; clique=i) for i in 1:length(result.corr_sparsity.cliques)]
end
```

## 4. Handling Sparsity Structures

### 4.1 Correlative Sparsity (Cliques)

The moment matrix is extracted **per clique**:
- Each clique has its own variables and constraints
- Each clique has an independent moment matrix
- The API provides `clique::Int` parameter to select which one

### 4.2 Term Sparsity (Block Structure)

Within each clique, term sparsity creates **block-diagonal structure**:

```
cliques_term_sparsities[clique_idx][1].block_bases :: Vector{Vector{M}}
```

Each element is a basis for one diagonal block. We have two options:

**Option A: Reconstruct Full Block-Diagonal Matrix** (RECOMMENDED)
- Construct the full moment matrix from all blocks
- More intuitive for users familiar with standard formulation
- Easier integration with GNS reconstruction

**Option B: Return Individual Blocks**
- Return vector of smaller matrices
- More memory efficient but requires modified GNS reconstruction
- More complex API

**Recommendation**: Use Option A for simplicity and compatibility with existing `reconstruct()`.

### 4.3 Localizing Matrices (Future Extension)

Currently focus on moment matrices only. Localizing matrices can be added later:
```julia
# Future API
localizing_matrix(result, clique=1, constraint_idx=1)
```

## 5. Implementation Steps

### Step 1: Add `ProblemType` to `PolyOptResult`
**File**: `src/interface.jl`

```julia
# Add at top of file
@enum ProblemType MomentPrimal SOSDual

# Modify struct
struct PolyOptResult{T,P,M}
    objective::T
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
    problem_type::ProblemType
end

# Update cs_nctssos() around line 96
function cs_nctssos(...)
    # ... existing code ...
    return PolyOptResult(
        objective_value(problem_to_solve.model),
        corr_sparsity,
        cliques_term_sparsities,
        problem_to_solve.model,
        dualize ? SOSDual : MomentPrimal  # NEW
    )
end

# Update cs_nctssos_higher() around line 143
function cs_nctssos_higher(...)
    # ... existing code ...
    return PolyOptResult(
        objective_value(problem_to_solve.model),
        prev_res.corr_sparsity,
        cliques_term_sparsities,
        problem_to_solve.model,
        dualize ? SOSDual : MomentPrimal  # NEW
    )
end
```

### Step 2: Create `MomentMatrix` Struct
**File**: `src/gns.jl` (or new `src/moment_extraction.jl`)

```julia
struct MomentMatrix{T,M}
    matrix::Matrix{T}
    basis::Vector{M}
    clique_vars::Vector{Variable}
    clique_index::Int
end

function Base.show(io::IO, mm::MomentMatrix)
    println(io, "MomentMatrix for clique $(mm.clique_index)")
    println(io, "  Variables: $(mm.clique_vars)")
    println(io, "  Basis size: $(length(mm.basis))")
    println(io, "  Matrix size: $(size(mm.matrix))")
end
```

### Step 3: Implement Core Extraction Logic
**File**: `src/gns.jl` (or new `src/moment_extraction.jl`)

```julia
"""
Extract moment matrix values from a Moment (primal) problem model.

For Moment problems, the moment matrix constraint is:
    @constraint(model, moment_mtx in PSDCone())
where moment_mtx[i,j] is a JuMP expression with monomap variables.

The dual values of this constraint give the moment matrix entries.
"""
function extract_moment_matrix_primal(
    model::GenericModel{T},
    constraints::Vector{<:ConstraintRef},
    basis::Vector{M}
) where {T,M}
    # Constraints is a flattened vector of all PSD constraints
    # We need to find the one corresponding to this basis
    # The constraint matrix should be basis × basis

    n = length(basis)
    matrix = zeros(T, n, n)

    # Find the constraint with matching dimension
    for cons_ref in constraints
        cons_obj = constraint_object(cons_ref)
        if get_dim(cons_obj) == n
            # This is our moment matrix constraint
            dual_vals = dual(cons_ref)
            # dual_vals is a vectorized form, reshape to matrix
            matrix = reshape(dual_vals, n, n)
            break
        end
    end

    return matrix
end

"""
Extract moment matrix values from an SOS (dual) problem model.

For SOS problems, the equality constraints are:
    fα_constraints .== 0
where each fα contains coefficients from the dual matrix variables.

We need to reconstruct the moment matrix by:
1. Identifying which equality constraints correspond to our basis
2. Using dual values of these constraints to get moment values
"""
function extract_moment_matrix_dual(
    model::GenericModel{T},
    basis::Vector{M},
    monomap::Dict{M,<:AbstractJuMPScalar}  # Need this from somewhere
) where {T,M}
    # This is more complex because we need to map back from
    # the equality constraints to the moment matrix structure

    # Get all equality constraints
    eq_cons = all_constraints(model, GenericAffExpr{T,VariableRef}, MOI.EqualTo{T})

    n = length(basis)
    matrix = zeros(T, n, n)

    # For each entry (i,j) in the moment matrix:
    # The value is the dual of the constraint corresponding to basis[i] * basis[j]
    for i in 1:n, j in 1:n
        # Compute the monomial for this entry
        mono = simplify!(expval(_neat_dot3(basis[i], one(M), basis[j])), sa)

        # Find the constraint for this monomial
        # The dual value is the moment matrix entry
        # ... implementation depends on how we track constraint-to-monomial mapping
        matrix[i, j] = dual(...)  # Need to figure out constraint mapping
    end

    return Symmetric(matrix)  # Moment matrices are symmetric
end
```

### Step 4: Implement Main API Function
**File**: `src/gns.jl`

```julia
function moment_matrix(result::PolyOptResult{T,P,M}; clique::Int=1) where {T,P,M}
    # Validate inputs
    num_cliques = length(result.corr_sparsity.cliques)
    if clique < 1 || clique > num_cliques
        throw(ArgumentError("Invalid clique index $clique. Problem has $num_cliques cliques."))
    end

    # Get the term sparsity for this clique
    # First element is the moment matrix, rest are localizing matrices
    mom_term_sparsity = result.cliques_term_sparsities[clique][1]

    # Get the basis - need to reconstruct full basis from blocks
    full_basis = construct_full_basis(mom_term_sparsity.block_bases)

    # Extract matrix values based on problem type
    if result.problem_type == MomentPrimal
        matrix = extract_moment_matrix_primal(result.model, full_basis)
    else  # SOSDual
        matrix = extract_moment_matrix_dual(result.model, full_basis)
    end

    # Get clique variables
    clique_vars = result.corr_sparsity.cliques[clique]

    return MomentMatrix(matrix, full_basis, clique_vars, clique)
end

function moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}
    return [moment_matrix(result; clique=i)
            for i in 1:length(result.corr_sparsity.cliques)]
end

"""
Construct full basis from block-diagonal structure.

With term sparsity, the moment matrix is block-diagonal where each block
has its own basis. This function creates the full basis that indexes the
complete block-diagonal matrix.
"""
function construct_full_basis(block_bases::Vector{Vector{M}}) where {M}
    # If no term sparsity (single block), return it directly
    length(block_bases) == 1 && return block_bases[1]

    # Otherwise, concatenate all block bases
    # They should be disjoint or we need to handle overlaps
    return sorted_unique(reduce(vcat, block_bases))
end
```

### Step 5: Handle Term Sparsity Correctly
**File**: `src/gns.jl`

The challenge: `block_bases` contains potentially overlapping bases for diagonal blocks.

**Investigation needed**:
- Are blocks truly disjoint or do they overlap?
- Do we need to extract each block separately and assemble?

**Strategy A: Full Matrix Reconstruction** (if blocks can share basis elements)
```julia
function extract_with_term_sparsity(model, block_bases, problem_type)
    # Extract each block separately
    blocks = []
    for block_basis in block_bases
        if problem_type == MomentPrimal
            block_matrix = extract_moment_matrix_primal(model, block_basis)
        else
            block_matrix = extract_moment_matrix_dual(model, block_basis)
        end
        push!(blocks, (block_basis, block_matrix))
    end

    # Assemble into full basis and matrix
    full_basis = construct_full_basis(block_bases)
    n = length(full_basis)
    full_matrix = zeros(n, n)

    # Place each block in the appropriate location
    for (basis, block_mat) in blocks
        indices = [searchsortedfirst(full_basis, b) for b in basis]
        full_matrix[indices, indices] = block_mat
    end

    return full_matrix, full_basis
end
```

### Step 6: Update Module Exports
**File**: `src/NCTSSoS.jl`

```julia
export moment_matrix, moment_matrices, MomentMatrix
export ProblemType, MomentPrimal, SOSDual  # If users need these
```

### Step 7: Write Comprehensive Tests
**File**: `test/moment_extraction.jl` (NEW)

```julia
using Test, NCTSSoS
using Clarabel

@testset "Moment Matrix Extraction" begin
    @testset "Simple Unconstrained Problem - Primal" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        pop = polyopt(f)

        # Solve without dualization (Moment problem)
        result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer); dualize=false)

        # Extract moment matrix
        mm = moment_matrix(result)

        # Validate structure
        @test mm.clique_index == 1
        @test length(mm.clique_vars) == 2
        @test size(mm.matrix, 1) == size(mm.matrix, 2)
        @test issymmetric(mm.matrix)

        # Validate basis
        @test all(m -> variables(m) ⊆ mm.clique_vars, mm.basis)
    end

    @testset "Simple Unconstrained Problem - Dual" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        pop = polyopt(f)

        # Solve with dualization (SOS problem)
        result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer); dualize=true)

        # Extract moment matrix
        mm = moment_matrix(result)

        # Should produce similar results to primal
        @test mm.clique_index == 1
        @test issymmetric(mm.matrix)
    end

    @testset "Multi-Clique Problem" begin
        @ncpolyvar x[1:4]
        # Create problem that will have multiple cliques
        f = x[1] * x[2] + x[3] * x[4]
        g1 = 1 - x[1]^2 - x[2]^2
        g2 = 1 - x[3]^2 - x[4]^2

        pop = polyopt(f; ineq_constraints=[g1, g2])
        result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer, cs_algo=MF()))

        # Extract all moment matrices
        mms = moment_matrices(result)

        # Should have multiple cliques
        @test length(mms) >= 1

        # Each should be valid
        for mm in mms
            @test issymmetric(mm.matrix)
            @test size(mm.matrix, 1) == length(mm.basis)
        end
    end

    @testset "Integration with GNS Reconstruction" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2 - 1
        pop = polyopt(f)

        result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer, order=2))
        mm = moment_matrix(result)

        # Should be compatible with reconstruct()
        # Get the degree from the basis
        H_deg = maximum(degree.(mm.basis))

        # This should not throw
        X_mats = reconstruct(mm.matrix, mm.clique_vars, H_deg; atol=1e-6)

        @test length(X_mats) == length(mm.clique_vars)
    end

    @testset "Error Handling" begin
        @ncpolyvar x[1:2]
        f = x[1]^2 + x[2]^2
        pop = polyopt(f)
        result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer))

        # Invalid clique index
        @test_throws ArgumentError moment_matrix(result; clique=0)
        @test_throws ArgumentError moment_matrix(result; clique=100)
    end
end
```

### Step 8: Update Documentation
**File**: `docs/src/api.md` (or appropriate documentation file)

Add section:
```markdown
## Moment Matrix Extraction

After solving a polynomial optimization problem, you can extract the moment matrix
for use with GNS reconstruction or other analysis.

### Basic Usage

\`\`\`julia
@ncpolyvar x[1:2]
f = x[1]^2 + x[2]^2 - 1
pop = polyopt(f)
result = cs_nctssos(pop, SolverConfig(optimizer=Clarabel.Optimizer, order=2))

# Extract moment matrix
mm = moment_matrix(result)

# Access components
println("Matrix size: ", size(mm.matrix))
println("Basis: ", mm.basis)
println("Variables: ", mm.clique_vars)

# Use with GNS reconstruction
X_mats = reconstruct(mm.matrix, mm.clique_vars, 2)
\`\`\`

### Multi-Clique Problems

For problems with correlative sparsity, extract matrices for each clique:

\`\`\`julia
# Extract all moment matrices
mms = moment_matrices(result)

# Or extract specific clique
mm1 = moment_matrix(result; clique=1)
mm2 = moment_matrix(result; clique=2)
\`\`\`
```

## 6. Open Questions & Design Decisions

### 6.1 Constraint Tracking for SOS Problems

**Issue**: In SOS dual problems, we need to map equality constraint dual values back to moment matrix entries.

**Challenge**: The `sos_dualize()` function doesn't currently store the mapping from monomials to constraints.

**Solution Options**:

**Option A**: Store additional metadata in `PolyOptResult` (RECOMMENDED)
```julia
struct PolyOptResult{T,P,M}
    # ... existing fields ...
    problem_type::ProblemType
    monomap::Union{Nothing, Dict{M,<:AbstractJuMPScalar}}  # NEW
    # monomap is Nothing for SOSDual, filled for MomentPrimal
end
```

**Option B**: Reconstruct mapping from model structure (fragile)

**Option C**: Only support moment matrix extraction for primal problems (limiting)

**Recommendation**: Implement Option A. Store `monomap` from `MomentProblem` in `PolyOptResult` to enable extraction from both primal and dual.

### 6.2 Term Sparsity Block Assembly

**Question**: When term sparsity creates blocks, are they:
- Truly disjoint (no overlapping basis elements)?
- Partially overlapping?
- Hierarchical (blocks within blocks)?

**Investigation Needed**:
- Review `iterate_term_sparse_supp()` logic
- Check if `block_bases` can share elements
- Look at test cases with term sparsity

**Current Assumption**: Blocks are disjoint or can be placed in a block-diagonal matrix. Need to verify.

### 6.3 Complex-Valued Problems

**Question**: How to handle `ComplexPolyOpt` and `ComplexMomentProblem`?

The GNS reconstruction currently only handles real matrices. Complex moment problems need special handling:
- Real and imaginary parts stored separately
- Different constraint structure in `complex_moment_solver.jl`

**Recommendation**: Initially support only real problems. Add complex support in follow-up PR.

```julia
function moment_matrix(result::PolyOptResult{T,P,M}; clique::Int=1) where {T,P,M}
    if T <: Complex
        throw(ArgumentError("Moment matrix extraction for complex problems not yet implemented"))
    end
    # ... implementation for real case ...
end
```

## 7. Implementation Checklist

- [ ] **Step 1**: Add `ProblemType` enum and field to `PolyOptResult`
  - [ ] Modify struct definition
  - [ ] Update `cs_nctssos()` to pass problem type
  - [ ] Update `cs_nctssos_higher()` to pass problem type
  - [ ] Update `Base.show()` method if needed

- [ ] **Step 2**: Create `MomentMatrix` struct
  - [ ] Define struct with fields
  - [ ] Implement `Base.show()` method
  - [ ] Add to exports

- [ ] **Step 3**: Implement primal extraction
  - [ ] `extract_moment_matrix_primal()` function
  - [ ] Handle constraint identification
  - [ ] Handle dual value extraction

- [ ] **Step 4**: Implement dual extraction
  - [ ] Decide on monomap storage approach
  - [ ] `extract_moment_matrix_dual()` function
  - [ ] Map equality constraints to moment entries

- [ ] **Step 5**: Implement main API
  - [ ] `moment_matrix()` function with validation
  - [ ] `moment_matrices()` convenience function
  - [ ] `construct_full_basis()` helper
  - [ ] Handle term sparsity block structure

- [ ] **Step 6**: Write tests
  - [ ] Simple primal problem
  - [ ] Simple dual problem
  - [ ] Multi-clique problem
  - [ ] Term sparsity problem
  - [ ] Integration with GNS reconstruction
  - [ ] Error handling

- [ ] **Step 7**: Update documentation
  - [ ] API documentation
  - [ ] Usage examples
  - [ ] Integration with existing workflows

- [ ] **Step 8**: Handle edge cases
  - [ ] Empty cliques
  - [ ] Global constraints
  - [ ] Complex problems (error message)

## 8. File Structure Summary

**Modified Files**:
- `src/interface.jl`: Add `ProblemType` and modify `PolyOptResult`
- `src/gns.jl`: Add extraction functions and `MomentMatrix` struct
- `src/NCTSSoS.jl`: Update exports

**New Files**:
- `test/moment_extraction.jl`: Comprehensive tests

**Documentation Updates**:
- Update API docs with new functions
- Add examples to user guide

## 9. Testing Strategy

### 9.1 Unit Tests
- Test each extraction method independently
- Test basis construction with various sparsity patterns
- Test error handling for invalid inputs

### 9.2 Integration Tests
- Test full workflow: solve → extract → GNS reconstruction
- Compare primal and dual extraction results (should match)
- Test with real examples from physics (Heisenberg, Ising models)

### 9.3 Regression Tests
- Ensure existing functionality unchanged
- Check that `PolyOptResult` still serializes/displays correctly

## 10. Performance Considerations

### 10.1 Memory
- For large problems, moment matrices can be very large
- Consider lazy evaluation or views where possible
- Document memory requirements in API docs

### 10.2 Computation
- Dual value extraction is typically O(n²) for n×n matrix
- Block assembly is O(total_basis_size)
- Should be negligible compared to solving the SDP

### 10.3 Type Stability
- Ensure all functions are type-stable
- Use `@code_warntype` to verify
- Parametric types properly threaded through

## 11. Future Extensions

### 11.1 Localizing Matrix Extraction
```julia
localizing_matrix(result; clique=1, constraint_idx=1)
localizing_matrices(result; clique=1)
```

### 11.2 Complex Problem Support
Handle complex moment problems with proper splitting of real/imaginary parts.

### 11.3 Sparse Matrix Support
For very large problems, return sparse matrices when appropriate.

### 11.4 Validation Functions
```julia
validate_flat_extension(mm::MomentMatrix)
check_psd_property(mm::MomentMatrix)
```

### 11.5 Direct GNS Integration
```julia
# Convenience function
reconstruct_from_result(result::PolyOptResult; clique=1, atol=1e-3)
```

## 12. Known Limitations

1. **Initial Implementation**: Real problems only (no Complex support)
2. **Dual Extraction**: May require additional metadata storage
3. **Block Assembly**: Assumes disjoint or properly nested block structure
4. **Global Constraints**: Handling may need special attention
5. **Memory**: Large dense matrices may be prohibitive for huge problems

## 13. Success Criteria

The implementation will be considered successful when:

1. ✓ Users can extract moment matrices with a single function call
2. ✓ Extraction works for both primal and dual formulations
3. ✓ Multi-clique problems are properly handled
4. ✓ Term sparsity is correctly handled
5. ✓ Integration with `reconstruct()` works seamlessly
6. ✓ Comprehensive tests pass
7. ✓ Documentation is clear and includes examples
8. ✓ No breaking changes to existing API

## 14. Estimated Complexity

- **Code Changes**: ~200-300 lines of implementation
- **Tests**: ~150-200 lines
- **Documentation**: ~50-100 lines
- **Total Effort**: Medium complexity (2-3 days for experienced Julia developer)

Most complex parts:
1. Dual extraction with proper constraint mapping
2. Term sparsity block assembly
3. Ensuring type stability throughout

## 15. Alternative Approaches Considered

### Approach A: Store Full MomentProblem (Rejected)
- Store entire `MomentProblem` or `SOSProblem` in `PolyOptResult`
- **Pros**: Complete information available
- **Cons**: Memory overhead, unnecessary data retention

### Approach B: Reconstruct from Scratch (Rejected)
- Re-run moment relaxation to get constraints
- **Pros**: No structural changes needed
- **Cons**: Wasteful computation, requires original problem

### Approach C: Lazy Extraction (Future)
- Return a lazy object that extracts on-demand
- **Pros**: Memory efficient
- **Cons**: More complex API, harder to debug

## 16. Summary

This plan provides a comprehensive path to implementing moment matrix extraction from `PolyOptResult` objects. The key design decisions are:

1. **Add metadata**: Include `problem_type` field to distinguish primal/dual
2. **Clean API**: Single `moment_matrix()` function with optional clique selection
3. **Rich return type**: `MomentMatrix` struct with matrix, basis, and metadata
4. **Proper sparsity handling**: Support both correlative and term sparsity
5. **Phased implementation**: Start with real problems, extend to complex later

The implementation maintains backward compatibility while providing the essential functionality needed for GNS reconstruction and other moment matrix analyses.
