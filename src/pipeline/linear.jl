import SparseArrays

struct DerivationSink{W,C}
    index::Dict{W,Int}
    row::Int
    rows::Vector{Int}
    cols::Vector{Int}
    coefs::Vector{C}
end
function push_term!(sink::DerivationSink{W,C}, c::Number, w::W) where {W,C}
    iszero(c) && return sink
    col = get(sink.index,w,0)
    col != 0 || throw(ArgumentError("The moment domain is not closed under the derivation; supply its linear closure."))
    push!(sink.rows,sink.row); push!(sink.cols,col); push!(sink.coefs,convert(C,c))
    sink
end

"""
    derivation_constraints(words, derivations; C=ComplexF64)

Sparse generator equations `A*y == 0` for invariant moments on an explicitly
closed domain. No group elements or moment matrix are formed. Exact generator
coefficients stay exact when `C` is exact; orthonormalization is not required.
`Symmetry(continuous=...)` joins these generator relations with finite orbit
moments and eliminates dependent coordinates. Continuous block splitting is
a separate step; moment elimination keeps each assembled block's shape.
"""
function derivation_constraints(words::AbstractVector{W}, derivations;
    C::Type=ComplexF64) where {W<:AbstractWord}
    index = Dict(w=>i for (i,w) in enumerate(words))
    length(index) == length(words) || throw(ArgumentError("The moment domain contains duplicate words."))
    rows,cols,coefs = Int[],Int[],C[]
    row = 0
    for d in derivations, w in words
        row += 1
        sink = DerivationSink(index,row,rows,cols,coefs)
        derive!(sink,d,w,one(C))
    end
    SparseArrays.sparse(rows,cols,coefs,row,length(words))
end
