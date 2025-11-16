"""
SUMMARY OF CHANGES - Moment Matrix Extraction Feature

This NEW file implements the core data structures and support building logic for moment matrix extraction.

Major components:
1. Hierarchical data structures: BlockSupport → CliqueSupport → MomentSupport
2. build_moment_support() function to construct the mapping during problem construction
3. Type-flexible design supporting both StateWord and NCStateWord monomials

Context: Part of moment matrix extraction implementation (.claude/tasks/moment-matrix-extraction/)
Plan reference: plan.md sections 2 (Data Structures) and 3.2.4 (Build MomentSupport)
"""

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

# CHANGED: New hierarchical data structure for moment matrix extraction
# WHY: Enables efficient extraction via direct indexing instead of recomputing basis products
# DESIGN DECISION: Hierarchical (cliques → blocks → indices) vs flat global support
#   - Chosen: Hierarchical for semantic clarity and self-documenting code
#   - Alternative: NCTSSOS uses flat support with algorithmic lookup
#   - Trade-off: Slightly more complex structure, but much clearer intent
# IMPACT: Core abstraction that mirrors problem sparsity structure

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
# IMPLEMENTATION NOTE: Only stores indices, not monomials
# WHY: Indices are what's actually needed for extraction; monomials only intermediate
# BENEFIT: Minimal storage overhead (~n² integers per block)
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
# DESIGN: Mirrors the clique-based sparsity structure
# WHY: Natural organization that matches how the problem is structured
# USAGE: Enables code like `moment_support.cliques[i].blocks[j].dual_indices[r,c]`
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
# TYPE PARAMETER: M can be StateWord, NCStateWord, or Monomial{Variable}
# WHY: Flexibility needed because canonicalization changes type (NCStateWord → StateWord)
# IMPACT: Separate M2 parameter in build_moment_support for type flexibility
struct MomentSupport{M}
    cliques::Vector{CliqueSupport}
    global_support::Vector{M}
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
# CHANGED: New function to build hierarchical support structure
# WHY: Reuse computation from constraint building instead of recomputing during extraction
# KEY INSIGHT: Support structure built ONCE during solving, used ONCE during extraction
# ALGORITHM: Matches basis product computation in moment_solver.jl and sos_solver.jl
#   - Uses _neat_dot3(basis[i], 1, basis[j]) to match moment matrix construction
#   - Applies simplify() and canonicalize() to match constraint ordering
#   - Looks up index in global_support dictionary for O(1) mapping
# CRITICAL: Only processes first TermSparsity (moment matrix), not localizing matrices
#   WHY: We only want to extract the moment matrix, not localizing matrices
# TYPE FLEXIBILITY: Separate M and M2 parameters
#   WHY: block_basis elements (M) may differ from canonicalized global_support (M2)
#   EXAMPLE: NCStateWord basis → StateWord canonicalized support
function build_moment_support(
    corr_sparsity::CorrelativeSparsity{P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}},
    global_support::Vector{M2},
    sa::SimplifyAlgorithm
) where {P,M,M2}

    # STEP 1: Create O(1) lookup for monomial → index mapping
    # WHY: Avoid O(n) binary search for each basis product lookup
    # PERFORMANCE: Dictionary creation O(n), lookups O(1) vs repeated binary search O(n log n)
    support_map = Dict(mono => idx for (idx, mono) in enumerate(global_support))

    # STEP 2: Build hierarchical structure mirroring cliques → blocks
    # ALLOCATION: Pre-allocate for type stability
    cliques = Vector{CliqueSupport}(undef, length(cliques_term_sparsities))

    for (clq_idx, clique_ts_vec) in enumerate(cliques_term_sparsities)
        # CHANGED: Use push! instead of pre-allocation
        # WHY: Don't know number of blocks ahead of time (only processing moment matrices)
        blocks = BlockSupport[]

        # CRITICAL: Only process the FIRST TermSparsity (moment matrix), not localizing matrices
        # WHY: Moment matrices correspond to dual variables; localizing matrices are auxiliary
        # CONTEXT: clique_ts_vec[1] = moment matrix, clique_ts_vec[2:end] = localizing matrices
        moment_matrix_ts = clique_ts_vec[1]

        for block_basis in moment_matrix_ts.block_bases
            n = length(block_basis)
            dual_indices = Matrix{Int}(undef, n, n)

            # CHANGED: Fill dual_indices matrix by computing basis products
            # ALGORITHM: Matches moment_solver.jl constraint construction
            #   Line 111: expval(_neat_dot3(row_idx, mono, col_idx))
            #   For moment matrix (poly=1), simplifies to: expval(_neat_dot3(row, 1, col))
            # WHY: Must use EXACT same computation to ensure indices align with constraints
            for i in 1:n, j in 1:n
                # Compute basis product matching moment_solver.jl:
                # expval(_neat_dot3(row, monomial, col))
                # For moment matrix (poly=1), this simplifies to expval(_neat_dot3(row, 1, col))
                # CRITICAL: simplify THEN canonicalize to match constraint construction order
                mono = simplify(expval(_neat_dot3(block_basis[i], one(block_basis[i]), block_basis[j])), sa)
                # Canonicalize to match global_support
                # WHY: global_support is canonicalized, so must canonicalize lookup key
                mono = canonicalize(mono, sa)
                # Look up index in global support
                # PERFORMANCE: O(1) dictionary lookup vs O(log n) binary search
                dual_indices[i, j] = support_map[mono]
            end

            push!(blocks, BlockSupport(dual_indices))
        end

        cliques[clq_idx] = CliqueSupport(blocks)
    end

    # RETURN: Complete hierarchical structure
    # BENEFIT: Self-documenting organization that mirrors problem structure
    return MomentSupport(cliques, global_support)
end
