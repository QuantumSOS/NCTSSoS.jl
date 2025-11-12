Base.cmp(a::Variable, b::Variable) = (a.name == b.name) ? 0 : (a.name < b.name ? -1 : 1)

Base.isless(a::Variable, b::Variable) = cmp(a, b) < 0

function Base.in(a::Variable, collection::Vector{Variable})
    return !isempty(searchsorted(collection, a))
end

function Base.cmp(a::Monomial, b::Monomial)
    degree(a) != degree(b) && return degree(a) < degree(b) ? -1 : 1
    a_idx, b_idx = 1, 1
    while a_idx <= length(a.vars) && b_idx <= length(b.vars)
        iszero(a.z[a_idx]) && (a_idx += 1; continue)
        iszero(b.z[b_idx]) && (b_idx += 1; continue)
        a.vars[a_idx] != b.vars[b_idx] && return cmp(a.vars[a_idx], b.vars[b_idx])
        # since degree(a) == degree(b) there must be something after a_idx
        a.z[a_idx] < b.z[b_idx] && return cmp(a.vars[a_idx + 1], b.vars[b_idx])
        # similarily here
        a.z[a_idx] > b.z[b_idx] && return cmp(a.vars[a_idx], b.vars[b_idx + 1])
        a_idx += 1
        b_idx += 1
    end
    # if I reach here, should trigger a warning
    return error("Comparison between monomials of equal degree, $a and $b failed")
end

function Base.in(a::Monomial, collection::Vector{Monomial})
    return !isempty(searchsorted(collection, a))
end

Base.isless(a::Monomial, b::Monomial) = cmp(a, b) == -1

Base.:(==)(a::Monomial, b::Monomial) = iszero(cmp(a, b))

function Base.:(==)(p::Polynomial{T}, q::Polynomial{T}) where {T}
    length(p.monos) != length(q.monos) && return false
    all(p.monos .== q.monos) || return false
    all(p.coeffs .== q.coeffs) || return false
    return true
end
