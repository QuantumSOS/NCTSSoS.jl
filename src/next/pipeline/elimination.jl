"""
    NoElimination <: EliminationAlgorithm

Elimination strategy that treats the entire problem as a single dense block.

This strategy ignores sparsity structure and solves the full SDP without
clique decomposition. Use when the problem is dense or small enough that
sparsity exploitation provides no benefit.

See also: [`MF`](@ref), [`MMD`](@ref), [`AsIsElimination`](@ref)
"""
struct NoElimination <: EliminationAlgorithm end

"""
    AsIsElimination <: EliminationAlgorithm

Elimination strategy that uses maximal cliques without triangulation.

!!! warning
    Assumes the correlative sparsity graph is already chordal. If the graph
    is not chordal, the resulting cliques may not satisfy the Running
    Intersection Property (RIP), which could produce invalid bounds. Use
    `MF()` or `MMD()` for automatic chordal completion.

See also: [`MF`](@ref), [`MMD`](@ref), [`NoElimination`](@ref)
"""
struct AsIsElimination <: EliminationAlgorithm end

"""
    MaximalElimination <: EliminationAlgorithm

Elimination strategy that decomposes the graph into connected components.

Each connected component becomes a separate SDP block. This is the most
aggressive sparsity exploitation, useful when variables form disjoint groups.

See also: [`MF`](@ref), [`MMD`](@ref), [`NoElimination`](@ref)
"""
struct MaximalElimination <: EliminationAlgorithm end

"""
    NaturalElimination()

Chordal completion in ascending vertex order, with the clique-row selection of
ncpol2sdpa 1.14.0's NumPy `_generate_clique` backend. In correlative sparsity,
vertices include all registry variables (also unused ones), in registry-index
order. Objective monomials connect their supports separately; each constraint
connects its entire support. Use `cs_algo=NaturalElimination()` explicitly;
existing defaults are unchanged.

This computes the symbolic nonzero pattern of natural-order Cholesky, avoiding
upstream's random matrix weights and numerical cancellation. It does **not**
implement ncpol2sdpa's alternative `_generate_clique_alt`: when both CVXOPT and
CHOMPACK are installed, upstream selects approximate-minimum-degree (AMD)
ordering and CHOMPACK cliques.
Neither `MF()` nor `MMD()` is a promise of compatibility with that backend.
"""
struct NaturalElimination <: EliminationAlgorithm end

function cliquetree(graph::AbstractGraph{V}, ::NoElimination, snd::SupernodeType) where {V}
    return cliquetree(complete_graph(nv(graph)), BFS(), snd)
end

"""
    clique_decomp(G::SimpleGraph, clique_alg::EliminationAlgorithm)

Decomposes a graph into cliques using the specified elimination algorithm.

# Arguments
- `G::SimpleGraph`: Input graph to decompose
- `clique_alg::EliminationAlgorithm`: Algorithm for clique tree elimination

# Returns
- `Vector{Vector{Int}}`: Vector of cliques, each containing vertex indices
"""
function clique_decomp(G::SimpleGraph, clique_alg::EliminationAlgorithm)
    label, tree = cliquetree(G, alg=clique_alg)
    return map(x -> label[x], collect(Vector{Int}, tree))
end

function clique_decomp(G::SimpleGraph, ::MaximalElimination)
    return connected_components(G)
end

function clique_decomp(G::SimpleGraph, ::NaturalElimination)
    filled = copy(G)
    rows = Vector{Int}[]
    cliques = Vector{Int}[]
    for i in vertices(filled)
        later = filter(j -> j > i, neighbors(filled, i))
        row = vcat(i, later)
        # Upstream tests against *all* previous factor rows, not just the
        # retained cliques. Preserve both row order and vertex order.
        if !any(previous -> issubset(row, previous), rows)
            push!(cliques, row)
        end
        push!(rows, row)
        add_clique!(filled, later)
    end
    return cliques
end

"""
    clique_decomp(G::SimpleGraph, ::AsIsElimination)

Returns maximal cliques of the graph without triangulation.

!!! warning
    This method assumes the graph is already chordal. If the graph is not chordal,
    the resulting cliques may not satisfy the Running Intersection Property (RIP),
    which could lead to invalid SOS relaxation bounds. Use `MF()` or `MMD()` for
    automatic chordal completion on non-chordal graphs.
"""
function clique_decomp(G::SimpleGraph, ::AsIsElimination)
    is_chordal, _ = CheckChordal(G)
    if !is_chordal
        @warn "AsIsElimination: the sparsity graph is not chordal. " *
              "The resulting clique decomposition may violate the Running Intersection " *
              "Property and produce invalid bounds. Use `MF()` or `MMD()` instead."
    end
    return maximal_cliques(G)
end

# TODO: https://github.com/QuantumSOS/NCTSSoS.jl/issues/45, this only applied to TS
