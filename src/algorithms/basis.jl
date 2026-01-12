@inline function _push_simplified_into_basis!(
    basis::Set{NormalMonomial{A,T}},
    ::Type{A},
    simplified::Vector{T},
) where {A<:MonoidAlgebra,T<:Integer}
    push!(basis, NormalMonomial{A,T}(simplified))
    return nothing
end

@inline function _push_simplified_into_basis!(
    basis::Set{NormalMonomial{A,T}},
    ::Type{A},
    simplified::Tuple{Vector{T},UInt8},
) where {A<:TwistedGroupAlgebra,T<:Integer}
    word, _phase_k = simplified
    push!(basis, NormalMonomial{A,T}(word))
    return nothing
end

@inline function _push_simplified_into_basis!(
    basis::Set{NormalMonomial{A,T}},
    ::Type{A},
    simplified::Vector{Tuple{Int,Vector{T}}},
) where {A<:PBWAlgebra,T<:Integer}
    for (c, word) in simplified
        iszero(c) && continue
        push!(basis, NormalMonomial{A,T}(word))
    end
    return nothing
end

function get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType,T<:Integer}
    d == 0 && return [one(NormalMonomial{A,T})]

    idxs = indices(registry)
    num_vars = length(idxs)
    num_vars == 0 && return NormalMonomial{A,T}[]

    basis_length = sum(num_vars^i for i in 0:d)
    words = [T[]]
    sizehint!(words, basis_length)

    for idx in idxs
        push!(words, [idx])
    end

    last_end_idx = 1

    @inbounds for i in 2:d
        last_end_idx += num_vars^(i - 1)
        for k in (num_vars^(i - 1)):-1:1
            for j in 1:num_vars
                push!(words, [words[last_end_idx - k + 1]; idxs[j]])
            end
        end
    end

    basis_set = Set{NormalMonomial{A,T}}()
    sizehint!(basis_set, basis_length)

    for word in words
        simplified = simplify!(A, word)
        _push_simplified_into_basis!(basis_set, A, simplified)
    end

    return sort!(collect(basis_set))
end
