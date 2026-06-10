# =============================================================================
# SympleQ Step 5: colour-preserving graph automorphisms
# =============================================================================

"""
    TermPermutation

Permutation of Pauli-term vertices. `images[i] == j` means term vertex `i` maps
to term vertex `j`.
"""
struct TermPermutation
    images::Vector{Int}

    function TermPermutation(images::AbstractVector{<:Integer})
        n = length(images)
        sorted = sort!(collect(Int, images))
        sorted == collect(1:n) || throw(ArgumentError("TermPermutation images must be a permutation of 1:$n."))
        return new(collect(Int, images))
    end
end

Base.length(p::TermPermutation) = length(p.images)
Base.getindex(p::TermPermutation, i::Integer) = p.images[i]
Base.:(==)(a::TermPermutation, b::TermPermutation) = a.images == b.images
Base.hash(p::TermPermutation, h::UInt) = hash(p.images, h)
Base.isone(p::TermPermutation) = all(i -> p.images[i] == i, eachindex(p.images))

function Base.show(io::IO, p::TermPermutation)
    print(io, "TermPermutation(", p.images, ")")
end

function _automorphism_order(g::SympleQGraph)
    n = nv(g.graph)
    color_counts = Dict{Any,Int}()
    for color in g.colors
        color_counts[color] = get(color_counts, color, 0) + 1
    end

    # Adjacency-first greedy order: after seeding with a most-constrained
    # vertex, always place next an unordered vertex with the largest number of
    # already-ordered neighbours (ties: smaller colour class, higher degree,
    # lower index). Each placement is then immediately constrained by edges to
    # previously placed vertices, which prunes the backtracking search far
    # earlier than a connectivity-blind order. The enumerated automorphism set
    # is unchanged — only the traversal order differs.
    key(v) = (color_counts[g.colors[v]], -Graphs.degree(g.graph, v), v)
    order = Int[]
    placed = falses(n)
    ordered_neighbors = zeros(Int, n)

    for _ in 1:n
        best = 0
        for v in 1:n
            placed[v] && continue
            if best == 0 ||
               ordered_neighbors[v] > ordered_neighbors[best] ||
               (ordered_neighbors[v] == ordered_neighbors[best] && key(v) < key(best))
                best = v
            end
        end
        push!(order, best)
        placed[best] = true
        for u in Graphs.neighbors(g.graph, best)
            ordered_neighbors[u] += 1
        end
    end

    return order
end

function _automorphism_candidates(g::SympleQGraph)
    n = nv(g.graph)
    return [
        [w for w in 1:n if g.colors[w] == g.colors[v] && Graphs.degree(g.graph, w) == Graphs.degree(g.graph, v)]
        for v in 1:n
    ]
end

function _automorphism_consistent(
    graph::SimpleGraph{Int},
    assigned::Vector{Int},
    v::Int,
    w::Int,
)
    for u in Graphs.vertices(graph)
        image_u = assigned[u]
        image_u == 0 && continue
        Graphs.has_edge(graph, v, u) == Graphs.has_edge(graph, w, image_u) || return false
    end
    return true
end

function _enumerate_color_preserving_automorphisms(
    g::SympleQGraph;
    max_automorphisms::Int=100_000,
)
    n = nv(g.graph)
    length(g.colors) == n || throw(ArgumentError("Graph colour vector length does not match vertex count."))

    order = _automorphism_order(g)
    candidates = _automorphism_candidates(g)
    assigned = zeros(Int, n)
    used = falses(n)
    automorphisms = Vector{Int}[]

    function search(depth::Int)
        length(automorphisms) >= max_automorphisms && return nothing
        if depth > n
            push!(automorphisms, copy(assigned))
            return nothing
        end

        v = order[depth]
        for w in candidates[v]
            used[w] && continue
            _automorphism_consistent(g.graph, assigned, v, w) || continue
            assigned[v] = w
            used[w] = true
            search(depth + 1)
            used[w] = false
            assigned[v] = 0
        end
        return nothing
    end

    search(1)
    length(automorphisms) >= max_automorphisms && @warn(
        "SympleQ automorphism enumeration hit max_automorphisms=$max_automorphisms; returned a truncated group."
    )
    return automorphisms
end

function _restrict_to_pauli_vertices(g::SympleQGraph, automorphism::Vector{Int})
    m = length(g.pauli_vertices)
    images = Vector{Int}(undef, m)
    pauli_lookup = Dict(vertex => i for (i, vertex) in enumerate(g.pauli_vertices))

    for (i, vertex) in enumerate(g.pauli_vertices)
        image_vertex = automorphism[vertex]
        images[i] = get(pauli_lookup, image_vertex, 0)
        images[i] == 0 && throw(ArgumentError(
            "Colour-preserving automorphism mapped Pauli vertex $vertex to non-Pauli vertex $image_vertex."
        ))
    end

    return TermPermutation(images)
end

"""
    automorphism_generators(graph::SympleQGraph; backend=:backtracking)

Return a deterministic generating set of colour-preserving automorphisms,
restricted to the Pauli-term vertices.

No extra dependency is pulled in for this MVP. `backend=:bliss` is accepted as a
request, but currently falls back to the deterministic in-tree backtracking
solver with a warning; this keeps non-SympleQ users free of graph-automorphism
binary dependencies.
"""
function automorphism_generators(
    graph::SympleQGraph;
    backend::Symbol=:backtracking,
    max_automorphisms::Int=100_000,
)
    if backend in (:auto, :backtracking)
        # native path
    elseif backend == :bliss
        @warn "Bliss backend is not wired yet; using the in-tree backtracking automorphism solver."
    else
        throw(ArgumentError("Unsupported SympleQ graph automorphism backend $backend."))
    end

    autos = _enumerate_color_preserving_automorphisms(graph; max_automorphisms)
    restricted = TermPermutation[]
    seen = Set{Tuple{Vararg{Int}}}()
    for auto in autos
        perm = _restrict_to_pauli_vertices(graph, auto)
        key = Tuple(perm.images)
        key in seen && continue
        push!(seen, key)
        isone(perm) || push!(restricted, perm)
    end

    return restricted
end
