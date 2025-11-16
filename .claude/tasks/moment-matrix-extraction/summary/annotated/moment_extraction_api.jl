"""
SUMMARY OF CHANGES - Moment Matrix Extraction API

This NEW file implements the public API for extracting moment matrices from solved problems.

Major components:
1. get_moment_matrices() - Main public function
2. extract_dual_variables() - For dual SOS formulation
3. extract_primal_variables() - For primal moment formulation
4. reconstruct_moment_matrices() - Direct indexing reconstruction

Context: Part of moment matrix extraction implementation (.claude/tasks/moment-matrix-extraction/)
Plan reference: plan.md sections 3.2.1-3.2.3 (Extraction Functions)
"""

# Moment Matrix Extraction API
#
# This file contains the public API functions for extracting moment matrices
# from solved optimization problems. It is included after interface.jl to
# have access to the PolyOptResult type.

# ============================================================================
# Core Extraction Functions
# ============================================================================

"""
    get_moment_matrices(result::PolyOptResult{T,P,M}) where {T,P,M}

Extract moment matrices from a solved optimization problem.

Returns a hierarchical structure of moment matrices organized by cliques and blocks,
following the problem's sparsity pattern.

# Arguments
- `result::PolyOptResult`: The solved optimization result containing the JuMP model
  and moment support structure.

# Returns
- `Vector{Vector{Matrix{T}}}`: Hierarchical moment matrices where:
  - Outer vector: One entry per clique
  - Middle vector: One entry per block within the clique
  - Inner matrix: Dense `Matrix{T}` of size `block_size × block_size`

# Example
```julia
# Solve the problem
result = cs_nctssos(pop; solver_config=config)

# Extract moment matrices
moments = get_moment_matrices(result)

# Access moment matrix for clique 1, block 1
M11 = moments[1][1]

# Access specific entry
value = moments[i][j][row, col]
```

# Algorithm
1. Extract dual variables from the solved JuMP model's equality constraints
2. Use the precomputed `moment_support` structure to map dual variables to
   moment matrix entries via direct indexing
3. Construct symmetric moment matrices organized by cliques and blocks

# Notes
- The dual variables correspond to the affine constraints that link the Gram
  matrices (or moment matrices in primal form) to the polynomial basis.
- The function works for both primal moment formulations and dual SOS formulations.
- Matrices are returned as dense `Matrix{T}`, not wrapped in `Symmetric`, to
  allow flexible manipulation by the user.
"""
# CHANGED: NEW public API function for moment matrix extraction
# WHY: Linear workflow - users explicitly call this after solving
# DESIGN DECISION: Separate function vs embedded in PolyOptResult
#   - Chosen: Separate function for zero overhead when not needed
#   - Alternative: Automatic extraction during solving
#   - Trade-off: Extra function call, but better separation of concerns
# IMPLEMENTATION: Delegates to formulation-specific extractors
# TYPE SAFETY: Type parameter M2 added to PolyOptResult for canonicalized support type
function get_moment_matrices(result::PolyOptResult{T,P,M,M2}) where {T,P,M,M2}
    # CHANGED: Branch based on formulation type (primal vs dual)
    # WHY: Different extraction methods for different formulations
    #   - Dual SOS: Extract dual variables from equality constraints
    #   - Primal moment: Extract primal variable values
    # CRITICAL: is_dual flag set during solving (cs_nctssos parameter)
    # CONTEXT: Dual formulation is DEFAULT (dualize=true) in cs_nctssos
    if result.is_dual
        # Dual SOS formulation: extract dual variables from equality constraints
        # MATHEMATICAL RELATIONSHIP: Moment values = dual variables of polynomial constraints
        # WHY: In dual formulation, moments appear as duals, not primal variables
        moment_values = extract_dual_variables(result.model)
    else
        # Primal moment formulation: extract primal variable values
        # MATHEMATICAL RELATIONSHIP: Moment values = primal variables directly
        # WHY: In primal formulation, moments are the optimization variables
        moment_values = extract_primal_variables(result.model)
    end

    # CHANGED: Reconstruct moment matrices using hierarchical structure
    # WHY: Direct indexing - no recomputation, no lookups
    # PERFORMANCE: O(n_blocks × block_size²) pure array indexing
    # KEY INSIGHT: moment_support precomputed during solving, reused here
    moment_mats = reconstruct_moment_matrices(
        moment_values,
        result.moment_support
    )

    return moment_mats
end

"""
    extract_dual_variables(model::GenericModel{T}) where {T}

Extract dual variables from the solved JuMP model.

Queries all affine equality constraints in the model and retrieves their
dual values. The ordering matches the global support vector used during
constraint construction.

# Arguments
- `model::GenericModel{T}`: The solved JuMP model.

# Returns
- `Vector{T}`: Vector of dual variable values, ordered to match the global
  support vector.

# Notes
- Applies a negative sign to dual values following the NCTSSOS convention.
- Only extracts duals from affine equality constraints (`AffExpr == constant`).
- Used for dual SOS formulation.
"""
# CHANGED: NEW function to extract dual variables from JuMP model
# WHY: Dual SOS formulation stores moment values as dual variables
# CRITICAL INSIGHT: Dual variables of polynomial equality constraints = moment values
# MATHEMATICAL FOUNDATION: Lagrange multipliers of primal constraints = dual variables
# NCTSSOS PATTERN: Analyzed from nccpop.jl lines 586-597
function extract_dual_variables(model::GenericModel{T}) where {T}
    # CHANGED: Query all affine equality constraints
    # WHY: These constraints link Gram matrices to polynomial coefficients
    # IMPLEMENTATION: Uses JuMP's constraint reflection API
    # ORDERING: Constraint order matches global_support order (critical!)
    con_refs = all_constraints(model, AffExpr, MOI.EqualTo{T})

    # CHANGED: Extract duals with negative sign
    # WHY: NCTSSOS convention (see NCTSSOS/src/nccpop.jl line 593)
    # MATHEMATICAL REASON: Sign convention in dual formulation
    # CRITICAL: Must use negative sign or values will be wrong
    dual_vals = -dual.(con_refs)

    return dual_vals
end

"""
    extract_primal_variables(model::GenericModel{T}) where {T}

Extract primal variable values from the solved JuMP model.

For the primal moment formulation, moment values are stored as primal variables
in the JuMP model. This function extracts their values in the order they were created.

# Arguments
- `model::GenericModel{T}`: The solved JuMP model.

# Returns
- `Vector{T}`: Vector of primal variable values, ordered to match the global
  support vector.

# Notes
- Used for primal moment formulation.
- The ordering matches the total_basis order from moment_relax.
"""
# CHANGED: NEW function to extract primal variables from JuMP model
# WHY: Primal moment formulation stores moment values as primal variables
# IMPLEMENTATION: Simpler than dual extraction (no sign flip needed)
# ORDERING: Variable order matches total_basis order from moment_relax
# USAGE: Only called when dualize=false in cs_nctssos
function extract_primal_variables(model::GenericModel{T}) where {T}
    # CHANGED: Get all variables from the model
    # WHY: In primal formulation, moments ARE the optimization variables
    # ORDERING: all_variables() returns in creation order = total_basis order
    all_vars = all_variables(model)

    # CHANGED: Extract primal values
    # WHY: These are the solved moment values directly
    # NO SIGN FLIP: Unlike dual extraction, primal values used as-is
    primal_vals = value.(all_vars)

    return primal_vals
end

"""
    reconstruct_moment_matrices(
        moment_values::Vector{T},
        moment_support::MomentSupport{M}
    ) where {T,M}

Reconstruct block-wise moment matrices from moment values using the hierarchical
support structure.

# Arguments
- `moment_values::Vector{T}`: Moment values from the solved optimization. These can be
  either dual variables (for dual SOS formulation) or primal variables (for primal
  moment formulation).
- `moment_support::MomentSupport{M}`: Hierarchical mapping from matrix positions
  to moment value indices.

# Returns
- `Vector{Vector{Matrix{T}}}`: Moment matrices organized as cliques → blocks → matrices.

# Algorithm
For each clique and each block within that clique:
1. Allocate a matrix of the appropriate size
2. Fill entries using direct indexing: `mat[i,j] = moment_values[indices[i,j]]`
3. Return the hierarchical structure

# Performance
This function performs no monomial computations or dictionary lookups - only
direct array indexing, making it very efficient even for large problems.
"""
# CHANGED: NEW function to reconstruct matrices from moment values
# WHY: Core extraction logic - maps values to matrix structure
# ALGORITHM: Pure direct indexing - no computation, no lookups
# PERFORMANCE: O(n_blocks × block_size²) array indexing operations
# KEY INSIGHT: All hard work done during support construction; this is trivial
# DESIGN BENEFIT: Separation of concerns - support building vs matrix reconstruction
function reconstruct_moment_matrices(
    moment_values::Vector{T},
    moment_support::MomentSupport{M}
) where {T,M}

    # CHANGED: Pre-allocate hierarchical structure
    # WHY: Type stability and performance
    # STRUCTURE: Vector of vectors of matrices = cliques → blocks → matrices
    moment_mats = Vector{Vector{Matrix{T}}}(undef, length(moment_support.cliques))

    # CHANGED: Iterate over cliques
    # WHY: Process each clique's moment matrices independently
    # CONTEXT: Mirrors correlative sparsity structure
    for (clq_idx, clique_support) in enumerate(moment_support.cliques)
        clique_mats = Vector{Matrix{T}}(undef, length(clique_support.blocks))

        # CHANGED: Iterate over blocks within clique
        # WHY: Each block corresponds to one moment matrix
        # CONTEXT: Mirrors term sparsity structure (moment matrix blocks)
        for (blk_idx, block_support) in enumerate(clique_support.blocks)
            n = size(block_support.dual_indices, 1)
            mat = zeros(T, n, n)

            # CHANGED: Fill matrix using direct indexing
            # ALGORITHM: mat[i,j] = moment_values[dual_indices[i,j]]
            # PERFORMANCE: O(n²) array indexing - no recomputation!
            # KEY INSIGHT: dual_indices precomputed during support construction
            # CRITICAL SIMPLIFICATION vs NCTSSOS:
            #   - NCTSSOS: Recomputes basis products, looks up in dictionary
            #   - NCTSSoS: Direct indexing into precomputed mapping
            #   - Result: Much simpler and faster extraction
            for i in 1:n, j in 1:n
                val_idx = block_support.dual_indices[i, j]
                mat[i, j] = moment_values[val_idx]
            end

            # CHANGED: Store matrix (not wrapped in Symmetric)
            # WHY: Allow flexible user manipulation
            # NOTE: Matrices should be symmetric by construction, but not enforced
            clique_mats[blk_idx] = mat
        end

        moment_mats[clq_idx] = clique_mats
    end

    # RETURN: Hierarchical structure matching problem sparsity
    # BENEFIT: Users can navigate: moments[clique][block][row,col]
    return moment_mats
end
