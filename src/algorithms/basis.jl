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

"""
    get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A,T} -> Vector{NormalMonomial{A,T}}

Generate a basis of monomials up to degree `d` from the variables in `registry`.

This function enumerates all words of length 0 to `d` over the variable indices,
applies algebra-specific simplification, and returns the unique simplified monomials
in sorted order.

# Arguments
- `registry`: Variable registry containing the variable indices
- `d`: Maximum degree (inclusive)

# Returns
A sorted vector of unique `NormalMonomial{A,T}` elements up to degree `d`.

# Examples
```jldoctest
julia> using NCTSSoS

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> length(basis)  # 1 + 2 + 4 = 7 monomials
7

julia> degree.(basis)
7-element Vector{Int64}:
 0
 1
 1
 2
 2
 2
 2
```

For algebras with simplification (e.g., UnipotentAlgebra where x² = x):
```jldoctest
julia> using NCTSSoS

julia> reg, (u,) = create_unipotent_variables([("u", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> length(basis)  # Fewer due to x² = x simplification
5
```

See also: [`get_state_basis`](@ref), [`VariableRegistry`](@ref)
"""
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
