"""
    Projector Algebra Simplification

Implements simplification for projector operators satisfying Pᵢ² = Pᵢ (idempotency).

# Algebraic Rules
- Pᵢ² = Pᵢ (idempotency: projectors square to themselves)
- Operators on different sites commute (sorted by site ascending)
- Operators on the same site are non-commutative (order preserved)
- Consecutive identical operators collapse: Pᵢ^k → Pᵢ

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
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> word = simplify(ProjectorAlgebra, UInt16[idx1_s2, idx1_s1, idx1_s2]);

julia> word == UInt16[idx1_s1, idx1_s2]
true
```
"""

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_word(::Type{ProjectorAlgebra}, word::Vector{T}) where {T<:Unsigned}

Check that a projector word is in canonical form. Throws `ArgumentError` if invalid.

Canonical form requirements:
- Sites sorted in ascending order
- No consecutive identical operators (no P² terms)

This is used by `NormalMonomial{ProjectorAlgebra,T}` constructor to enforce invariants.
"""
function _validate_word(::Type{ProjectorAlgebra}, word::Vector{T}) where {T<:Unsigned}
    _is_site_sorted(word) || throw(ArgumentError("ProjectorAlgebra word must be site-sorted. Use `simplify` for raw words."))
    @inbounds for i in 1:length(word)-1
        word[i] == word[i+1] && throw(ArgumentError("ProjectorAlgebra word must not contain consecutive identical operators. Use `simplify` for raw words."))
    end
end

# =============================================================================
# Simplification
# =============================================================================

"""
    simplify!(::Type{ProjectorAlgebra}, word::Vector{T}) where {T<:Unsigned} -> Vector{T}

In-place site-aware simplification for projector algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, idempotency applies: consecutive identical operators collapse (P² = P).

Returns the simplified word vector (mutated in place).

# Algorithm
1. Stable sort by site (using `decode_site`)
2. Remove consecutive duplicates (idempotency)
"""
function simplify!(::Type{ProjectorAlgebra},word::Vector{T}) where {T<:Unsigned}
    filter!(!iszero, word)
    # Empty or single: nothing to simplify
    length(word) <= 1 && return word

    # Stable sort by site (operators on different sites commute, within-site order preserved)
    _stable_sort_by_site!(word)

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
    simplify(::Type{ProjectorAlgebra}, w::Vector{T}) where {T<:Unsigned} -> Vector{T}

Simplify a raw ProjectorAlgebra word into canonical form.

This is the primary entry point for ProjectorAlgebra simplification. Takes a raw word vector
and returns a simplified copy (sorted by site, idempotency applied).
"""
simplify(::Type{ProjectorAlgebra}, w::Vector{T}) where {T<:Unsigned} = simplify!(ProjectorAlgebra,copy(w))
