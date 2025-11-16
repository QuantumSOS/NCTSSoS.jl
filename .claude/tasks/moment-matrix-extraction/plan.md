# Implementation Plan: Moment Matrix Extraction

## Overview

This plan details how to add moment matrix extraction functionality to NCTSSoS.jl after solving relaxed sum of squares problems. The design follows a **linear, modular workflow** where:

1. **Solving phase**: `cs_nctssos()` solves the problem and stores the support list in `PolyOptResult`
2. **Extraction phase**: Users call `get_moment_matrices(result)` to extract moment matrices

**Key Design Principles**:
- **Separation of concerns**: Solving and extraction are independent operations
- **Zero API changes**: `cs_nctssos()` signature remains unchanged
- **Hierarchical structure**: `MomentSupport` mirrors cliques → blocks hierarchy
- **Semantic clarity**: Code structure explicitly represents problem sparsity
- **Efficient extraction**: Direct indexing via precomputed `dual_indices`, no lookups
- **Built once, used once**: Support structure created during constraint building, used during extraction
- **Type stable and performant**: Follows Julia best practices

**Core Data Structure**:
```julia
struct MomentSupport{M}
    cliques::Vector{CliqueSupport}   # Hierarchical: cliques → blocks → dual_indices
    global_support::Vector{M}         # Reference: sorted unique monomials
end
```

---

## Section 1: Design Decisions

### 1.1 Where to Store Moment Matrices

**Decision**: Create a separate standalone function that extracts and returns moment matrices. Do NOT modify `PolyOptResult`.

**Rationale**:
- **Better separation of concerns**: Solving and extraction are independent operations
- **Linear procedure**: Users explicitly call: `result = cs_nctssos(...); moments = get_moment_matrices(result)`
- **Zero overhead**: Users who don't need moments pay no cost (no extra fields, no memory)
- **More flexible**: Can extract moments at any time after solving, or not at all
- **Minimal code changes**: PolyOptResult remains unchanged
- **Functional programming**: Pure function that takes result and returns moments

**Proposed API**:
```julia
# Solve the problem
result = cs_nctssos(pop; solver_config=config)

# Extract moment matrices (separate step)
moment_matrices = get_moment_matrices(result)
```

**Alternative Considered**: Add `moment_matrices` field to PolyOptResult
- **Rejected**: Couples solving with extraction, adds complexity, less modular

### 1.2 Support List Storage

**Decision**: Store the support list in `PolyOptResult` during solving, not recompute during extraction.

**Rationale**:
- Computing basis products is actually time consuming for large problems
- Support list is already computed during `moment_relax()` / `sos_dualize()` when creating constraints
- Storing it reuses existing computation rather than duplicating work
- Negligible memory overhead compared to the JuMP model and solution data
- NCTSSOS stores it in `struct_data.tsupp` for this reason

**Implementation**:
- Add `support::Vector{M}` field to `PolyOptResult`
- Populate during constraint creation in `moment_relax()` / `sos_dualize()`
- Use stored support during moment extraction

### 1.3 Constraint Tracking

**Decision**: Use constraint names/tags in JuMP model to identify moment constraints.

**Rationale**:
- JuMP models already track constraints
- Can query dual values by constraint reference
- No need for additional bookkeeping structures
- Follows JuMP best practices
- Eg: `model[:cons] = @constraint `

**Approach**: In `moment_relax()` and `sos_dualize()`, name constraints consistently so they can be retrieved later.

### 1.4 Real vs Complex Support

**Decision**: Support both real and complex moment matrices using type parameters.

**Rationale**:
- NCTSSoS.jl already handles both via `ComplexPolyOpt` and `ComplexMomentProblem`
- Type parameter `T` in `PolyOptResult` already tracks numeric type
- Marginal additional complexity

---

## Section 2: Data Structures

### 2.1 Hierarchical Support Structure

**New Data Structures** for semantic clarity and efficient extraction:

```julia
"""
Support information for a single block within a clique.

Maps (row, col) indices in the block's moment matrix to the corresponding
index in the global dual variable vector.
"""
struct BlockSupport
    dual_indices::Matrix{Int}   # [i,j] = index into global dual_vars vector
end

"""
Support information for a single clique.

Contains support data for all blocks within this clique.
"""
struct CliqueSupport
    blocks::Vector{BlockSupport}  # One entry per block in the clique
end

"""
Complete hierarchical support structure for the optimization problem.

Mirrors the sparsity hierarchy: cliques → blocks → (row, col) basis pairs.
The global_support is kept for reference/debugging.
"""
struct MomentSupport{M}
    cliques::Vector{CliqueSupport}   # One entry per clique
    global_support::Vector{M}         # Flat, sorted, unique global support (optional, for reference)
end
```

**Key Design Points**:
- **`dual_indices`**: The essential information - maps matrix position to dual variable index
- **`monomials`**: NOT stored - only used during construction, not needed for extraction
- **Hierarchy**: Explicit structure mirrors cliques → blocks, making code self-documenting
- **`global_support`**: Optional, kept for debugging/inspection

### 2.2 Extended PolyOptResult

**Modification**: Add `moment_support` field to store hierarchical support structure.

```julia
struct PolyOptResult{T,P,M}
    objective::T
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
    moment_support::MomentSupport{M}  # NEW: Hierarchical support structure
end
```

**Field Details**:
- `moment_support`: Hierarchical mapping from (clique, block, row, col) to dual variable indices
- Built during `moment_relax()` / `sos_dualize()` when creating constraints
- Self-documenting structure that mirrors the problem's sparsity hierarchy
- Enables fast extraction via direct indexing: `dual_vars[dual_indices[i,j]]`

**Benefits of This Design**:
- **Semantic clarity**: Structure explicitly represents clique/block hierarchy
- **Efficient extraction**: Direct indexing, no recomputation or dictionary lookups
- **Self-documenting**: Code naturally navigates cliques → blocks → entries
- **Type-safe**: Compiler enforces correct hierarchy
- **Minimal storage**: Only store indices (Int), not monomials

### 2.3 Return Type for get_moment_matrices()

The extraction function returns:

```julia
Vector{Vector{Matrix{T}}}  # Hierarchical: cliques → blocks → matrices
```

**Structure**:
- Outer vector: One entry per clique
- Middle vector: One entry per block within the clique
- Inner matrix: Symmetric dense matrix of size `block_size × block_size`
- Type `T`: Matches the numeric type from PolyOptResult (Float64, ComplexF64, etc.)

**Preserves semantic structure**: Users can navigate results by clique and block, matching the problem structure.

---

## Section 3: Implementation Files and Functions

### 3.1 File: src/interface.jl

**Modifications**:
1. Update `PolyOptResult` struct (add `support::Vector{M}` field)
2. Update constructor to include `support` parameter
3. NO changes to `cs_nctssos()` function signature

**Struct Change**:
```julia
struct PolyOptResult{T,P,M}
    objective::T
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
    support::Vector{M}  # NEW
end
```

The `support` vector will be populated by `moment_relax()` / `sos_dualize()` and passed through when constructing the result.

### 3.2 File: src/moment_extraction.jl (NEW)

**Purpose**: Core extraction logic

**Functions**:

#### 3.2.1 Main Extraction Function
```julia
"""
    get_moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}

Extract moment matrices from a solved optimization problem.

# Arguments
- `result::PolyOptResult`: The solved optimization result containing the model and moment_support

# Returns
- `Vector{Vector{Matrix{T}}}`: Hierarchical moment matrices (cliques → blocks → matrices)

# Example
```julia
result = cs_nctssos(pop; solver_config=config)
moments = get_moment_matrices(result)

# Access moment matrix for clique i, block j
moment_matrix = moments[i][j]
```
"""
function get_moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}
    # 1. Extract dual variables from model constraints
    dual_vars = extract_dual_variables(result.model)

    # 2. Reconstruct moment matrices using hierarchical structure
    moment_mats = reconstruct_moment_matrices(
        dual_vars,
        result.moment_support
    )

    return moment_mats
end
```

#### 3.2.2 Dual Variable Extractor
```julia
"""
    extract_dual_variables(model::GenericModel{T}) where {T}

Extract dual variables from the solved JuMP model.

Returns a vector corresponding to the support list.
"""
function extract_dual_variables(model::GenericModel{T}) where {T}
    # Approach depends on how constraints are named/stored
    # Option 1: If constraints are stored with consistent naming
    con_refs = all_constraints(model, AffExpr, MOI.EqualTo{T})

    # Extract duals (note: negative sign as per NCTSSOS pattern)
    dual_vals = -dual.(con_refs)

    return dual_vals
end
```

#### 3.2.3 Moment Matrix Reconstructor
```julia
"""
    reconstruct_moment_matrices(dual_vars, moment_support)

Reconstruct block-wise moment matrices from dual variables using the hierarchical support structure.

# Arguments
- `dual_vars::Vector{T}`: Dual variable values from the solved optimization
- `moment_support::MomentSupport{M}`: Hierarchical mapping structure

# Returns
- `Vector{Vector{Matrix{T}}}`: Moment matrices organized by cliques → blocks
"""
function reconstruct_moment_matrices(
    dual_vars::Vector{T},
    moment_support::MomentSupport{M}
) where {T,M}

    moment_mats = Vector{Vector{Matrix{T}}}(undef, length(moment_support.cliques))

    # For each clique
    for (clq_idx, clique_support) in enumerate(moment_support.cliques)
        clique_mats = Vector{Matrix{T}}(undef, length(clique_support.blocks))

        # For each block in the clique
        for (blk_idx, block_support) in enumerate(clique_support.blocks)
            n = size(block_support.dual_indices, 1)
            mat = zeros(T, n, n)

            # Fill matrix using direct indexing (no recomputation!)
            for i in 1:n, j in 1:n
                dual_idx = block_support.dual_indices[i, j]
                mat[i, j] = dual_vars[dual_idx]
            end

            clique_mats[blk_idx] = mat
        end

        moment_mats[clq_idx] = clique_mats
    end

    return moment_mats
end
```

**Key simplification**: No monomial recomputation, no lookups - just direct indexing using the precomputed `dual_indices`!

#### 3.2.4 Build MomentSupport Structure (Critical Helper)

This function will be called during `moment_relax()` / `sos_dualize()` to build the hierarchical support structure.

```julia
"""
    build_moment_support(
        corr_sparsity::CorrelativeSparsity{P,M},
        cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
        sa::SimplifyAlgorithm
    ) where {P,M}

Build the hierarchical MomentSupport structure during constraint creation.

This function:
1. Collects all basis products from all cliques and blocks
2. Builds a global sorted unique support vector
3. Creates the hierarchical dual_indices mapping

# Returns
- `MomentSupport{M}`: Complete hierarchical structure with dual_indices
"""
function build_moment_support(
    corr_sparsity::CorrelativeSparsity{P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
    sa::SimplifyAlgorithm
) where {P,M}

    # Step 1: Collect all monomials from all blocks
    all_monomials = M[]

    for (clq_idx, clique_ts_vec) in enumerate(cliques_term_sparsities)
        for ts in clique_ts_vec
            for block_basis in ts.block_bases
                for i in 1:length(block_basis)
                    for j in 1:length(block_basis)
                        # Compute basis product
                        mono = adjoint_product(block_basis[i], block_basis[j])
                        mono = simplify!(mono, sa)
                        push!(all_monomials, mono)
                    end
                end
            end
        end
    end

    # Step 2: Build global support (sorted, unique)
    global_support = unique(sort(all_monomials))

    # Step 3: Create lookup dict for monomial → index
    support_map = Dict(mono => idx for (idx, mono) in enumerate(global_support))

    # Step 4: Build hierarchical dual_indices structure
    cliques = Vector{CliqueSupport}(undef, length(cliques_term_sparsities))

    for (clq_idx, clique_ts_vec) in enumerate(cliques_term_sparsities)
        blocks = Vector{BlockSupport}(undef, length(clique_ts_vec))

        for (ts_idx, ts) in enumerate(clique_ts_vec)
            # For simplicity, assuming one block per term sparsity
            # (need to verify this assumption with actual code structure)
            for (blk_idx, block_basis) in enumerate(ts.block_bases)
                n = length(block_basis)
                dual_indices = Matrix{Int}(undef, n, n)

                # Fill dual_indices matrix
                for i in 1:n, j in 1:n
                    mono = adjoint_product(block_basis[i], block_basis[j])
                    mono = simplify!(mono, sa)
                    dual_indices[i, j] = support_map[mono]
                end

                blocks[blk_idx] = BlockSupport(dual_indices)
            end
        end

        cliques[clq_idx] = CliqueSupport(blocks)
    end

    return MomentSupport(cliques, global_support)
end
```

#### 3.2.5 Other Helper Functions
```julia
"""
    adjoint_product(mono1::M, mono2::M) where {M<:Monomial}

Compute mono1^† * mono2 (adjoint of mono1 times mono2).
"""
function adjoint_product(mono1::M, mono2::M) where {M}
    # This depends on FastPolynomials internals
    # Likely: reverse mono1, then concatenate with mono2
    # Implementation TBD based on FastPolynomials API
end
```

### 3.3 File: src/moment_solver.jl

**Modification**: Update `moment_relax()` to build and return the MomentSupport structure.

**Current signature** (approximate):
```julia
function moment_relax(
    pop::PolyOpt,
    corr_sparsity::CorrelativeSparsity,
    cliques_term_sparsities::Vector{Vector{TermSparsity}}
) -> MomentProblem
```

**Modified to also return MomentSupport**:
```julia
function moment_relax(
    pop::PolyOpt,
    corr_sparsity::CorrelativeSparsity,
    cliques_term_sparsities::Vector{Vector{TermSparsity}},
    sa::SimplifyAlgorithm
) -> (MomentProblem, MomentSupport{M})
```

**Implementation**:
```julia
function moment_relax(...)
    # ... existing code to create MomentProblem ...

    # NEW: Build moment support structure
    moment_support = build_moment_support(
        corr_sparsity,
        cliques_term_sparsities,
        sa
    )

    return (moment_problem, moment_support)
end
```

### 3.4 File: src/sos_solver.jl

**Modification**: Update `sos_dualize()` to build and return the MomentSupport structure.

**This is the more important case** since `dualize=true` is the default in `cs_nctssos()`.

**Modified signature**:
```julia
function sos_dualize(
    moment_problem::MomentProblem,
    corr_sparsity::CorrelativeSparsity,
    cliques_term_sparsities::Vector{Vector{TermSparsity}},
    sa::SimplifyAlgorithm
) -> (SOSProblem, MomentSupport{M})
```

**Implementation**: Similar to `moment_relax()`, call `build_moment_support()` and return it alongside the SOSProblem.

**Note**: The support structure is the same regardless of primal/dual formulation - it always maps basis products to dual variable indices.

---

## Section 4: Integration Points

### 4.1 Linear Workflow

**User-facing workflow**:
```julia
# Step 1: Solve the optimization problem
result = cs_nctssos(pop; solver_config=config)

# Step 2: Extract moment matrices (separate function call)
moments = get_moment_matrices(result)

# Access hierarchically: moments[clique_idx][block_idx] is a matrix
moment_matrix_clique1_block1 = moments[1][1]
```

**Internal flow**:
```
1. Build correlative sparsity
2. Build term sparsities
3. Create moment problem (or dualize to SOS)
   └─> Build MomentSupport structure (hierarchical dual_indices + global_support)
4. Solve optimization
5. Return PolyOptResult (now includes moment_support field)

[User can then call get_moment_matrices on result]

6. get_moment_matrices() extracts dual variables from model
7. get_moment_matrices() reconstructs matrices via direct indexing using moment_support
8. Returns Vector{Vector{Matrix{T}}} (cliques → blocks → matrices)
```

**Key Points**:
- MomentSupport built **once** during constraint creation
- Contains dual_indices: precomputed mapping from (clique, block, i, j) → dual_var index
- Extraction is **fast and simple**: just index into dual_vars using dual_indices
- No recomputation, no lookups, no dictionary overhead

### 4.2 Constraint Naming Strategy

**Challenge**: Need to map dual variables back to support monomials.

**Solution**: During moment_relax() or sos_dualize(), ensure constraints are created in a consistent order that matches the support list.

**Implementation**:
- Modify `moment_relax()` in src/moment_solver.jl to track constraint order
- Store constraint references in order corresponding to support
- During extraction, query duals in same order

---

## Section 5: Handling Dual SOS Formulation

### 5.1 Challenge

The default formulation is `dualize=true`, which converts the moment problem to a dual SOS problem:
- Original: Minimize ⟨c, y⟩ subject to moment matrix PSD constraints
- Dual: Maximize b subject to polynomial equality constraints

In the dual, moment values appear as **dual variables** of the polynomial equality constraints, not as primal variables.

### 5.2 Solution

The extraction works the same way:
1. The dual SOS problem has affine constraints linking Gram matrices to polynomial coefficients
2. These constraints correspond to monomials in the support
3. After solving, extract dual variables of these constraints
4. Map back to moment matrix structure

**Key insight from NCTSSOS**: The pattern is identical whether using primal or dual formulation - always extract dual variables of the equality constraints.

---

## Section 6: Open Questions to Resolve

### 6.1 SimplifyAlgorithm Access

**Question**: Where is the `SimplifyAlgorithm` stored for later access during extraction?

**Options**:
1. Store in `PolyOptResult` as additional field
2. Store in `CorrelativeSparsity` structure
3. Pass as parameter to extraction function
4. Reconstruct from problem structure

**Recommended**: Option 1 (store in PolyOptResult) for simplicity.

### 6.2 Constraint Ordering

**Question**: How to ensure constraints are created in predictable order matching support list?

**Challenge**: Current code may not guarantee specific constraint ordering.

**Options**:
1. Modify constraint creation to explicitly track ordering
2. Store constraint references with their corresponding monomials
3. Build a mapping during moment_relax/sos_dualize

**Recommended**: Option 3 (build mapping during problem construction).

### 6.3 FastPolynomials API

**Question**: What's the correct way to compute adjoint products in FastPolynomials?

**Need to investigate**:
- How to get adjoint of a monomial
- How to multiply monomials
- Whether `neat_dot` or similar functions already exist

**Action**: Review FastPolynomials source code for monomial operations.

### 6.4 Multiple Blocks per Clique

**Question**: How are moment matrices organized when there are multiple blocks per clique?

**From exploration**:
- `TermSparsity` has `block_bases::Vector{Vector{M}}`
- Each clique can have multiple term sparsity structures
- Each term sparsity can have multiple blocks

**Decision needed**:
- One moment matrix per block? (NCTSSOS approach)
- One aggregated matrix per clique?
- Hierarchical structure?

**Recommended**: One matrix per block (matches NCTSSOS, finest granularity).

---

## Section 7: Implementation Steps

### Step 1: Understand FastPolynomials
- Read FastPolynomials source for monomial operations
- Identify adjoint and product functions
- Test basic monomial manipulation
- Understand how `neat_dot` works

### Step 2: Create MomentSupport Data Structures
- Define `BlockSupport`, `CliqueSupport`, and `MomentSupport{M}` structs in src/moment_extraction.jl
- Add docstrings explaining the hierarchy
- Ensure proper type parameters

### Step 3: Modify PolyOptResult
- Add `moment_support::MomentSupport{M}` field to PolyOptResult struct
- Update PolyOptResult constructor to accept moment_support parameter
- Ensure all call sites are updated

### Step 4: Implement build_moment_support()
- Create the function that builds the hierarchical MomentSupport structure
- Implements the algorithm: collect monomials → build global_support → create dual_indices
- Test on simple examples to verify correctness

### Step 5: Modify Problem Construction Functions
- Update `moment_relax()` to call `build_moment_support()` and return MomentSupport
- Update `sos_dualize()` to call `build_moment_support()` and return MomentSupport
- Update `cs_nctssos()` to pass moment_support through to PolyOptResult

### Step 6: Implement Dual Extraction
- Create `extract_dual_variables()` function in src/moment_extraction.jl
- Test on solved models
- Ensure correct ordering (matching global_support vector)

### Step 7: Implement Moment Reconstruction
- Create `reconstruct_moment_matrices()` function
- Use direct indexing: `dual_vars[dual_indices[i,j]]`
- Return hierarchical structure: Vector{Vector{Matrix{T}}}

### Step 8: Implement Main Extraction Function
- Create `get_moment_matrices()` as the public API
- Compose the extraction and reconstruction functions
- Handle edge cases

### Step 9: Export and Documentation
- Export `get_moment_matrices` from src/NCTSSoS.jl
- Add docstrings with examples
- Create example in docs/src/examples/ showing hierarchical access

### Step 10: Testing
- Create unit tests for each component:
  - `build_moment_support()` correctness
  - `dual_indices` mapping accuracy
  - `reconstruct_moment_matrices()` with known values
- Create integration test with known solution
- Validate against NCTSSOS on same problem
- Test both real and complex problems
- Test multi-clique sparse problems
- Test hierarchical access patterns

---

## Section 8: Testing Strategy

### 8.1 Unit Tests

**File**: test/test_moment_extraction.jl

Tests for:
1. Support list construction
2. Monomial adjoint products
3. Support lookup
4. Dual variable extraction
5. Matrix reconstruction (with known dual values)

### 8.2 Integration Tests

**Test Cases**:
1. Simple 2-variable problem with analytical solution
2. Pauli algebra problem (compare with GNS reconstruction)
3. Problem with sparsity (multiple cliques)
4. Complex-valued problem

### 8.3 Validation

**Approach**:
- Solve same problem in NCTSSOS
- Extract moment matrices from both
- Compare values (should match within solver tolerance)

---

## Section 9: Files to Create/Modify

### Files to Modify:

1. **src/interface.jl**
   - Modify `PolyOptResult` struct: add `support::Vector{M}` field
   - Update `PolyOptResult` constructor
   - Update `cs_nctssos()` to pass support from moment_relax/sos_dualize to PolyOptResult

2. **src/moment_solver.jl**
   - Modify `moment_relax()` to build and return support list
   - Update return type to return `(MomentProblem, Vector{M})`
   - Collect basis products during constraint creation

3. **src/sos_solver.jl**
   - Modify `sos_dualize()` to build and return support list
   - Update return type accordingly
   - Build support similar to NCTSSOS pattern

4. **src/NCTSSoS.jl**
   - Export `get_moment_matrices` function

### Files to Create:

1. **src/moment_extraction.jl** - Core extraction logic
   - `get_moment_matrices()` - Main public API
   - `extract_dual_variables()` - Extract duals from model
   - `reconstruct_moment_matrices()` - Build matrices from duals
   - `adjoint_product()` - Helper for monomial adjoint products
   - `get_simplify_algorithm()` - Helper to access SA

2. **test/test_moment_extraction.jl** - Unit and integration tests
   - Test support list correctness
   - Test dual extraction
   - Test matrix reconstruction
   - Integration tests with known solutions

### Files to Reference:

1. `src/FastPolynomials/src/FastPolynomials.jl` - Understand monomial API
2. `src/gns.jl` - Reference for working with moment matrices
3. `src/sparse.jl` - Understanding basis structures
4. `/Users/yushengzhao/projects/NCTSSOS/src/nccpop.jl` - Reference implementation patterns

---

## Section 10: Potential Challenges

### 10.1 Constraint Access

**Challenge**: JuMP may not preserve constraint ordering.

**Mitigation**: Store explicit mapping during problem construction.

### 10.2 Simplification Algorithm

**Challenge**: Need SA for monomial reduction, may not be readily accessible.

**Mitigation**: Store SA in result structure or pass as parameter.

### 10.3 Type Stability

**Challenge**: Supporting both real and complex while maintaining performance.

**Mitigation**: Use type parameters consistently, avoid type instabilities.

### 10.4 Memory Usage

**Challenge**: Large problems may have many blocks, large moment matrices. Storing support in PolyOptResult adds memory overhead.

**Mitigation**:
- Support vector is relatively small (one monomial per constraint)
- Users only call `get_moment_matrices()` when they need moments
- Moment matrices are not stored in PolyOptResult, only returned from function
- Users can discard the result object if they only need moments

---

## Section 11: Success Criteria

Implementation is successful when:

1. ✓ Can extract moment matrices from solved `PolyOptResult`
2. ✓ Matrices have correct dimensions (match block bases)
3. ✓ Matrices are symmetric
4. ✓ Values match dual variables from optimization
5. ✓ Works with both real and complex problems
6. ✓ Works with sparse (multi-clique) problems
7. ✓ All existing tests still pass
8. ✓ New tests achieve 100% coverage of new code
9. ✓ Results validated against NCTSSOS on test problems
10. ✓ Documentation is complete and clear

---

## Next Steps

1. **Review this plan** section by section with user
2. **Resolve open questions** (Section 6)
3. **Investigate FastPolynomials** API for monomial operations
4. **Begin implementation** following the steps in Section 7
