"""
    Unipotent Algebra Simplification (Site-Based)

Implements simplification for unipotent operators satisfying U² = I (identity).

# Algebraic Rules
- U² = I (squares to identity, operators are self-inverse)
- Consecutive identical operators cancel: UᵢUᵢ → I (removes pair)
- No cyclic products or cross-operator interactions (unlike Pauli)
- Operators on different sites commute (sorted by site ascending)
- Operators on the same site are non-commutative (order preserved)

# Index Encoding
Indices must be unsigned integers (`T<:Unsigned`) with bit-packed site information.
Use `encode_index(T, operator_id, site)` to create indices.

# Algorithm
1. Group operators by site (using `decode_site`)
2. Sort groups by site (ascending)
3. Within each site: apply U²=I (remove consecutive pairs via stack)
4. Concatenate sorted groups

# Examples
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> word = simplify(UnipotentAlgebra, UInt16[idx1_s2, idx1_s1, idx1_s2]);

julia> word == UInt16[idx1_s1]
true
```

U²=I pair removal:
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> word = simplify(UnipotentAlgebra, UInt16[idx1_s1, idx1_s1]);

julia> isempty(word)
true
```
"""

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_word(::Type{UnipotentAlgebra}, word::Vector{T}) where {T<:Unsigned}

Check that a unipotent word is in canonical form. Throws `ArgumentError` if invalid.

Canonical form requirements:
- Sites sorted in ascending order
- No consecutive identical operators (no U² terms - they should cancel to I)

This is used by `NormalMonomial{UnipotentAlgebra,T}` constructor to enforce invariants.
"""
function _validate_word(::Type{UnipotentAlgebra},word::Vector{T}) where {T<:Unsigned}
    _is_site_sorted(word) || throw(ArgumentError("UnipotentAlgebra word must be site-sorted. Use `simplify` for raw words."))
    @inbounds for i in 1:length(word)-1
        word[i] == word[i+1] && throw(ArgumentError("UnipotentAlgebra word must not contain consecutive identical operators. Use `simplify` for raw words."))
    end
end

# =============================================================================
# Simplification
# =============================================================================

"""
    simplify!(::Type{UnipotentAlgebra}, word::Vector{T}) where {T<:Unsigned} -> Vector{T}

In-place site-aware simplification for unipotent algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, U²=I applies: consecutive identical operators cancel (remove pairs).

Returns the simplified word vector (mutated in place).

# Algorithm
1. Stable sort by site (using `decode_site`)
2. Apply U²=I via stack (remove consecutive pairs)
"""
function simplify!(::Type{UnipotentAlgebra},word::Vector{T}) where {T<:Unsigned}
    # Empty or single: nothing to simplify
    length(word) <= 1 && return word

    # Stable sort by site (operators on different sites commute, within-site order preserved)
    _stable_sort_by_site!(word)

    # Apply U²=I: remove consecutive identical pairs with backtracking
    i = 1
    while i < length(word)
        if word[i] == word[i + 1]
            # Consecutive identical: remove both (U² = I)
            deleteat!(word, i:i+1)
            # Backtrack to catch cascading cancellations
            i > 1 && (i -= 1)
        else
            i += 1
        end
    end

    return word
end

"""
    simplify(::Type{UnipotentAlgebra}, word::Vector{T}) where {T<:Unsigned} -> Vector{T}

Simplify a raw UnipotentAlgebra word into canonical form.

This is the primary entry point for UnipotentAlgebra simplification. Takes a raw word vector
and returns a simplified copy (sorted by site, U²=I applied).
"""
simplify(::Type{UnipotentAlgebra}, word::Vector{T}) where {T<:Unsigned} =  simplify!(UnipotentAlgebra, copy(word))
