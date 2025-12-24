struct NoElimination <: EliminationAlgorithm end
struct AsIsElimination <: EliminationAlgorithm end
struct MaximalElimination <: EliminationAlgorithm end

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
    return maximal_cliques(G)
end

# TODO: https://github.com/wangjie212/NCTSSoS.jl/issues/45, this only applied to TS
