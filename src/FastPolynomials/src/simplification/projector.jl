"""
    Projector Algebra Simplification

Implements simplification for projector operators satisfying Pᵢ² = Pᵢ (idempotency).

# Algebraic Rules
- Pᵢ² = Pᵢ (idempotency: projectors square to themselves)
- Operators on different sites commute (sorted by site ascending)
- Operators on the same site are non-commutative (order preserved)
- Consecutive identical operators collapse: PᵢPᵢPᵢ → Pᵢ

# Index Encoding
Indices must be unsigned integers (`T<:Unsigned`) with bit-packed site information.
Use `encode_index(T, operator_id, site)` to create indices.

# Algorithm
1. Group operators by site (using `decode_site`)
2. Sort groups by site (ascending)
3. Within each site: remove consecutive duplicates (idempotency)
4. Concatenate sorted groups

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{ProjectorAlgebra}([idx1_s2, idx1_s1, idx1_s2]);

julia> result = simplify(m);

julia> result.word == [idx1_s1, idx1_s2]
true
```
"""

"""
    _simplify_projector_word!(word::Vector{T}) where {T<:Unsigned} -> Vector{T}

In-place site-aware simplification for projector algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, idempotency applies: consecutive identical operators collapse (P² = P).

Returns the simplified word vector (mutated in place).

# Algorithm
1. Stable sort by site (using `decode_site`)
2. Remove consecutive duplicates (idempotency)
"""
function _simplify_projector_word!(word::Vector{T}) where {T<:Unsigned}
    # Empty or single: nothing to simplify
    length(word) <= 1 && return word

    # Stable sort by site (operators on different sites commute, within-site order preserved)
    sort!(word, alg=InsertionSort, by=decode_site)

    # Apply P²=P: remove consecutive duplicates (keep first of each run)
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            # Consecutive identical: remove duplicate (P² = P)
            deleteat!(word, i + 1)
            # No backtrack needed - just keep checking current position
        else
            i += 1
        end
    end

    return word
end

"""
    simplify(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned} -> Monomial

Simplify a projector algebra monomial with site-aware commutation and idempotency.

Returns a new simplified Monomial (no coefficient changes, just reordering and collapsing).
The original monomial is unchanged.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> m = Monomial{ProjectorAlgebra}([idx1_s1, idx1_s1, idx1_s1]);

julia> result = simplify(m);

julia> result.word == [idx1_s1]
true

julia> length(m.word)  # Original unchanged
3
```
"""
function simplify(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}
    word_copy = copy(m.word)
    _simplify_projector_word!(word_copy)
    Monomial{ProjectorAlgebra,T}(word_copy)
end
