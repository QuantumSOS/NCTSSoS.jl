# =============================================================================
# SympleQ Step 3: symplectic products and coloured anticommutation graph
# =============================================================================

"""
    SympleQGraph

Coloured graph used for SympleQ automorphism discovery.

The first `length(pauli_vertices)` vertices are Pauli-term vertices. Auxiliary
cycle vertices, when present, are stored after them and coloured separately so a
colour-preserving automorphism cannot confuse terms with cycle constraints.
"""
struct SympleQGraph
    graph::SimpleGraph{Int}
    colors::Vector{Any}
    pauli_vertices::Vector{Int}
    auxiliary_vertices::Vector{Int}
    cycles::Vector{Vector{Int}}
end

function Base.show(io::IO, g::SympleQGraph)
    print(
        io,
        "SympleQGraph($(nv(g.graph)) vertices, $(ne(g.graph)) edges, " *
        "$(length(g.pauli_vertices)) Pauli vertices, $(length(g.auxiliary_vertices)) cycle vertices)",
    )
end

@inline function _symplectic_product_rows(
    a::AbstractVector{<:Integer},
    b::AbstractVector{<:Integer},
)
    length(a) == length(b) || throw(DimensionMismatch("Symplectic rows have different lengths."))
    iseven(length(a)) || throw(ArgumentError("Symplectic rows must have even length."))
    n = length(a) ÷ 2
    acc = 0
    @inbounds for i in 1:n
        acc ⊻= (Int(a[i] & b[n + i]) & 1)
        acc ⊻= (Int(a[n + i] & b[i]) & 1)
    end
    return UInt8(acc & 1)
end

"""
    symplectic_product_matrix(tab::SymplecticTableau)

Return the binary matrix whose `(i,j)` entry is the Pauli symplectic product
`⟨pᵢ,pⱼ⟩ = xᵢ⋅zⱼ + zᵢ⋅xⱼ (mod 2)`.
"""
function symplectic_product_matrix(tab::SymplecticTableau)
    m = length(tab)
    products = zeros(UInt8, m, m)
    for j in 1:m, i in 1:(j - 1)
        value = _symplectic_product_rows(view(tab.paulis, i, :), view(tab.paulis, j, :))
        products[i, j] = value
        products[j, i] = value
    end
    return products
end

_term_color(coef, eta::Integer) = (:pauli, coef, Int(eta))

"""
    anticommutation_graph(tab::SymplecticTableau)

Build the colour-preserving graph for SympleQ Steps 1--3: Pauli-term vertices
are coloured by coefficient (and the canonical word phase `η`), with an edge for
each anticommuting pair.
"""
function anticommutation_graph(tab::SymplecticTableau)
    m = length(tab)
    graph = SimpleGraph(m)
    products = symplectic_product_matrix(tab)
    for j in 1:m, i in 1:(j - 1)
        products[i, j] == 1 && add_edge!(graph, i, j)
    end
    colors = Any[_term_color(tab.coeffs[i], tab.eta[i]) for i in 1:m]
    return SympleQGraph(graph, colors, collect(1:m), Int[], Vector{Int}[])
end
