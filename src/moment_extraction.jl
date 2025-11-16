# Moment Matrix Extraction
#
# This file implements functionality to extract moment matrices from solved
# optimization problems in NCTSSoS.jl.
#
# The extraction follows a hierarchical structure that mirrors the problem's
# sparsity pattern: cliques → blocks → moment matrix entries.

# ============================================================================
# Data Structures
# ============================================================================

"""
    BlockSupport

Support information for a single block within a clique.

Maps `(row, col)` indices in the block's moment matrix to the corresponding
index in the global dual variable vector.

# Fields
- `dual_indices::Matrix{Int}`: A matrix where `dual_indices[i,j]` gives the
  index into the global dual variable vector for moment matrix entry `(i,j)`.

# Example
```julia
# For a 3×3 block
block_support = BlockSupport(reshape(1:9, 3, 3))
# block_support.dual_indices[1,1] = 1  # corresponds to dual_vars[1]
# block_support.dual_indices[2,3] = 6  # corresponds to dual_vars[6]
```
"""
struct BlockSupport
    dual_indices::Matrix{Int}
end

"""
    CliqueSupport

Support information for a single clique.

Contains support data for all blocks within this clique.

# Fields
- `blocks::Vector{BlockSupport}`: One `BlockSupport` entry per block in the clique.

# Example
```julia
# Clique with 2 blocks
clique_support = CliqueSupport([
    BlockSupport(block1_indices),
    BlockSupport(block2_indices)
])
```
"""
struct CliqueSupport
    blocks::Vector{BlockSupport}
end

"""
    MomentSupport{M}

Complete hierarchical support structure for the optimization problem.

This structure mirrors the problem's sparsity hierarchy:
`cliques → blocks → (row, col) basis pairs`, and provides efficient mapping
from moment matrix positions to dual variable indices.

# Fields
- `cliques::Vector{CliqueSupport}`: One entry per clique in the problem.
- `global_support::Vector{M}`: Flat, sorted, unique vector of all monomials
  appearing in the problem. This is kept for reference and debugging.

# Type Parameters
- `M`: Monomial type (e.g., `Monomial{Variable}`)

# Usage
The `MomentSupport` structure is built during problem construction and stored
in `PolyOptResult`. It enables efficient moment matrix extraction via direct
indexing into the dual variable vector.

# Example
```julia
# Access dual index for clique i, block j, entry (row, col)
dual_idx = moment_support.cliques[i].blocks[j].dual_indices[row, col]
moment_value = dual_vars[dual_idx]
```
"""
struct MomentSupport{M}
    cliques::Vector{CliqueSupport}
    global_support::Vector{M}
    primal_var_map::Vector{Int}  # Maps global_support index -> variable index (for primal)
end

# ============================================================================
# Note: Core extraction API functions (get_moment_matrices, extract_dual_variables,
# reconstruct_moment_matrices) are in moment_extraction_api.jl, included after
# interface.jl to access PolyOptResult type.
# ============================================================================

# ============================================================================
# Support Structure Construction
# ============================================================================

"""
    build_moment_support(
        corr_sparsity::CorrelativeSparsity{P,M},
        cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
        global_support::Vector{M},
        sa::SimplifyAlgorithm
    ) where {P,M}

Build the hierarchical MomentSupport structure during constraint creation.

This function is called during problem construction (in `moment_relax()` or
`sos_dualize()`) to create the mapping between moment matrix entries and
dual variable indices.

# Arguments
- `corr_sparsity::CorrelativeSparsity{P,M}`: Correlative sparsity information.
- `cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}`: Term sparsity
  information for each clique and block.
- `global_support::Vector{M}`: Pre-built global support vector matching the
  constraint ordering. This must be the canonicalized basis used in constraint
  construction.
- `sa::SimplifyAlgorithm`: Algorithm for simplifying monomials to canonical form.

# Returns
- `MomentSupport{M}`: Complete hierarchical structure with dual indices and
  global support.

# Algorithm
1. **Create index mapping**: Build a dictionary mapping each monomial in the
   global support to its index.
2. **Build hierarchical structure**: For each clique's MOMENT MATRIX block
   (the first TermSparsity in each clique), and for each basis pair `(i,j)`,
   compute `neat_dot(basis[i], basis[j])`, canonicalize it, and look up its
   index in the global support to populate the `dual_indices` matrix.

# Notes
- Only processes MOMENT MATRIX blocks (first TermSparsity), not localizing matrices.
- The global_support parameter ensures exact alignment with constraint construction.

# Performance
This function is called once during problem construction. The computational cost
is dominated by monomial simplification and canonicalization.
"""
function build_moment_support(
    corr_sparsity::CorrelativeSparsity{P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
    global_support::Vector{M2},
    sa::SimplifyAlgorithm;
    primal_var_map::Vector{Int}=Int[]  # Optional: for primal formulation
) where {P,M,M2}

    # Step 1: Create lookup dict for monomial → index
    support_map = Dict(mono => idx for (idx, mono) in enumerate(global_support))

    # Step 2: Build hierarchical dual_indices structure
    cliques = Vector{CliqueSupport}(undef, length(cliques_term_sparsities))

    for (clq_idx, clique_ts_vec) in enumerate(cliques_term_sparsities)
        blocks = BlockSupport[]

        # Only process the FIRST TermSparsity (moment matrix), not localizing matrices
        moment_matrix_ts = clique_ts_vec[1]

        for block_basis in moment_matrix_ts.block_bases
            n = length(block_basis)
            dual_indices = Matrix{Int}(undef, n, n)

            # Fill dual_indices matrix
            for i in 1:n, j in 1:n
                # Compute basis product matching moment_solver.jl:
                # expval(_neat_dot3(row, monomial, col))
                # For moment matrix (poly=1), this simplifies to expval(_neat_dot3(row, 1, col))
                mono = simplify(expval(neat_dot(block_basis[i], block_basis[j])), sa)
                # Canonicalize to match global_support
                mono = canonicalize(mono, sa)
                # Look up index in global support
                dual_indices[i, j] = support_map[mono]
            end

            push!(blocks, BlockSupport(dual_indices))
        end

        cliques[clq_idx] = CliqueSupport(blocks)
    end

    return MomentSupport(cliques, global_support, primal_var_map)
end
