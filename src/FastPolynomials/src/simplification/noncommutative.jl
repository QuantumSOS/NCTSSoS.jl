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

julia> t = simplify(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2, idx1_s2]
true
```
"""

"""
    simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned} -> Monomial

Site-aware in-place simplification for non-commutative algebra with encoded indices.

Operators on different sites commute and are sorted by site (ascending).
Within each site, order is preserved exactly (no simplification rules apply).

Returns the simplified Monomial (no coefficient changes, just reordering).

# Algorithm
1. Group operators by site (using `decode_site`)
2. Sort groups by site (ascending)
3. Within each site: preserve order exactly
4. Concatenate sorted groups

# Warning
This mutates the input monomial. Use `simplify` for a non-mutating version.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{NonCommutativeAlgebra}([idx1_s2, idx1_s1]);

julia> result = simplify!(m);

julia> result.word == [idx1_s1, idx1_s2]
true

julia> m.word == [idx1_s1, idx1_s2]  # Original was mutated
true
```
"""
function simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    word = m.word

    # Empty or single: nothing to simplify
    length(word) <= 1 && return m

    # Stable sort by site: operators on different sites commute, within-site order preserved
    sort!(word, alg=Base.Sort.InsertionSort, by=decode_site)

    # Create new monomial with correct hash (can't update hash in-place since Monomial is immutable)
    simplified_mono = Monomial{NonCommutativeAlgebra}(word)
    return simplified_mono
end

"""
    simplify(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned} -> Monomial

Simplify a non-commutative algebra monomial with site-aware commutation.

Non-mutating version - creates a copy and simplifies it.
Returns the simplified Monomial (no coefficient changes, just reordering).

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
    # Copy and delegate to simplify!
    m_copy = Monomial{NonCommutativeAlgebra,T}(copy(m.word), m.hash)
    simplify!(m_copy)
end
