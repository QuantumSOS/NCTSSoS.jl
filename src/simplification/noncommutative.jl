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
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{NonCommutativeAlgebra}([idx1_s2, idx1_s1, idx1_s2]);

julia> result = simplify(m);

julia> result.word == [idx1_s1, idx1_s2, idx1_s2]
true
```
"""

"""
    _simplify_nc_word!(word::Vector{T}) where {T<:Unsigned} -> Vector{T}

In-place site-aware simplification for non-commutative algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, order is preserved exactly (no simplification rules apply).

Returns the sorted word vector (mutated in place).

# Algorithm
1. Stable sort by site (using `decode_site`)
2. Within each site: preserve order exactly
"""
function _simplify_nc_word!(word::Vector{T}) where {T<:Unsigned}
    # Empty or single: nothing to simplify
    length(word) <= 1 && return word

    # Stable sort by site: operators on different sites commute, within-site order preserved
    sort!(word, alg=Base.Sort.InsertionSort, by=decode_site)
    return word
end

"""
    simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned} -> Monomial

In-place simplification of a non-commutative algebra monomial.

Mutates the input monomial's word vector and returns the same monomial.

# Warning
This mutates the input monomial. Use `simplify` for a non-mutating version.
"""
function simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    _simplify_nc_word!(m.word)
    return m
end

"""
    simplify(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned} -> Monomial

Simplify a non-commutative algebra monomial with site-aware commutation.

Returns a new simplified Monomial (no coefficient changes, just reordering).
The original monomial is unchanged.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{NonCommutativeAlgebra}([idx1_s2, idx1_s1]);

julia> result = simplify(m);

julia> result.word == [idx1_s1, idx1_s2]
true

julia> m.word == [idx1_s2, idx1_s1]  # Original unchanged
true
```
"""
function simplify(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    word_copy = copy(m.word)
    _simplify_nc_word!(word_copy)
    Monomial{NonCommutativeAlgebra,T}(word_copy)
end
