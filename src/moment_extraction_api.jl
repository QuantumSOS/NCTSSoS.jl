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
function get_moment_matrices(result::PolyOptResult{T,P,M,M2}) where {T,P,M,M2}
    # 1. Extract values from model based on formulation type
    if result.is_dual
        # Dual SOS formulation: extract dual variables from equality constraints
        moment_values = extract_dual_variables(result.model)
    else
        # Primal moment formulation: extract primal variable values directly
        moment_values = extract_primal_variables(result.model)
    end

    # 2. Reconstruct moment matrices using hierarchical structure
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
- Extracts dual values directly without sign conversion.
- Only extracts duals from affine equality constraints (`AffExpr == constant`).
- Used for dual SOS formulation.
"""
function extract_dual_variables(model::GenericModel{T}) where {T}
    # Get all affine equality constraints
    con_refs = all_constraints(model, AffExpr, MOI.EqualTo{T})

    # Extract duals (NOTE: No negative sign needed for NCTSSoS.jl dual formulation)
    dual_vals = dual.(con_refs)

    return dual_vals
end

"""
    extract_primal_variables(model::GenericModel{T}) where {T}

Extract primal variable values from the solved JuMP model.

For the primal moment formulation, moment values are stored as primal variables
in the JuMP model. Since global_support = total_basis (both are canonical and identical),
we can extract values directly without reordering.

# Arguments
- `model::GenericModel{T}`: The solved JuMP model.

# Returns
- `Vector{T}`: Vector of primal variable values, ordered to match total_basis (= global_support).

# Notes
- Used for primal moment formulation.
- No reordering needed since total_basis is already canonical.
"""
function extract_primal_variables(model::GenericModel{T}) where {T}
    # Get all variables from the model (in total_basis = global_support order)
    all_vars = all_variables(model)

    # Extract primal values directly
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
function reconstruct_moment_matrices(
    moment_values::Vector{T},
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
                val_idx = block_support.dual_indices[i, j]
                mat[i, j] = moment_values[val_idx]
            end

            clique_mats[blk_idx] = mat
        end

        moment_mats[clq_idx] = clique_mats
    end

    return moment_mats
end

# Overload for PrimalMomentSupport - identical implementation
function reconstruct_moment_matrices(
    moment_values::Vector{T},
    moment_support::PrimalMomentSupport{M}
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
                val_idx = block_support.dual_indices[i, j]
                mat[i, j] = moment_values[val_idx]
            end

            clique_mats[blk_idx] = mat
        end

        moment_mats[clq_idx] = clique_mats
    end

    return moment_mats
end
