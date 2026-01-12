"""
    NonCommutative Algebra Simplification (Site-Based)

Implements simplification for non-commutative operators with site-based commutation.

# Algebraic Rules
- No simplification rules within a site (order preserved exactly)
- Operators on different sites commute (sorted by site ascending)
- Operators on the same site are non-commutative (order preserved)

# Index Encoding
Indices must be unsigned integers (`T<:Unsigned`) with bit-packed site information.
Use `encode_index(T, operator_id, site)` to create indices.

# Algorithm
1. Group operators by site (using `decode_site`)
2. Sort groups by site (ascending)
3. Within each site: preserve order exactly (no simplification)
4. Concatenate sorted groups

# Examples
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> word = simplify(NonCommutativeAlgebra, [idx1_s2, idx1_s1, idx1_s2]);

julia> word == [idx1_s1, idx1_s2, idx1_s2]
true
```
"""

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_word(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Unsigned}

Check that a noncommutative word is in canonical form. Throws `ArgumentError` if invalid.

Canonical form requirements:
- Sites sorted in ascending order (operators on different sites commute)

This is used by `NormalMonomial{NonCommutativeAlgebra,T}` constructor to enforce invariants.
"""
function _validate_word(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Unsigned}
    _is_site_sorted(word) || throw(ArgumentError("NonCommutativeAlgebra word must be site-sorted. Use `simplify` for raw words."))
end

_validate_word(::Type{NonCommutativeAlgebra}, ::Vector{T}) where {T<:Signed} = throw(ArgumentError(
    "NonCommutativeAlgebra indices must be unsigned (site-encoded). Use `encode_index` with `T<:Unsigned`."
))

# =============================================================================
# Simplification
# =============================================================================

"""
    simplify!(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Unsigned} -> Vector{T}

In-place site-aware simplification for non-commutative algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, order is preserved exactly (no simplification rules apply).

Returns the sorted word vector (mutated in place).

# Algorithm
1. Filter out zero indices
2. Stable sort by site (using `decode_site`)
3. Within each site: preserve order exactly
"""
function simplify!(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Unsigned}
    # Empty or single: nothing to simplify
    length(word) <= 1 && return word

    filter!(!iszero, word)
    # Stable sort by site: operators on different sites commute, within-site order preserved
    _stable_sort_by_site!(word)
    return word
end

"""
    simplify(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Unsigned} -> Vector{T}

Simplify a raw NonCommutative word into canonical form.

This is the primary entry point for NonCommutative simplification. Takes a raw word vector
and returns the simplified word vector (sorted by site).
"""
simplify(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Unsigned} = simplify!(NonCommutativeAlgebra, copy(word))

"""
    simplify(::Type{NonCommutativeAlgebra}, word::Vector{T}) where {T<:Signed}

Throws `ArgumentError`. Signed indices are not supported for `NonCommutativeAlgebra`.
Use `encode_index` with an unsigned integer type to construct site-encoded indices.
"""
simplify(::Type{NonCommutativeAlgebra}, ::Vector{T}) where {T<:Signed} = throw(ArgumentError(
    "NonCommutativeAlgebra indices must be unsigned (site-encoded). Use `encode_index` with `T<:Unsigned`."
))
