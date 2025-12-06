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
    simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}

Site-aware in-place simplification for non-commutative algebra with encoded indices.

Operators on different sites commute and are sorted by site (ascending).
Within each site, order is preserved exactly (no simplification rules apply).

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

julia> t = simplify!(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true

julia> m.word == [idx1_s1, idx1_s2]  # Original was mutated
true
```
"""
function simplify!(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    word = m.word

    # Empty or single: nothing to simplify
    length(word) <= 1 && return Term(1.0, m)

    # Group by site using Dict to collect operators per site, preserving order within site
    site_groups = Dict{Int,Vector{T}}()

    for idx in word
        site = decode_site(idx)
        if !haskey(site_groups, site)
            site_groups[site] = T[]
        end
        push!(site_groups[site], idx)
    end

    # Sort sites (ascending order)
    sorted_sites = sort!(collect(keys(site_groups)))

    # Build result: for each site, just concatenate (no simplification within site)
    empty!(word)

    for site in sorted_sites
        ops = site_groups[site]
        append!(word, ops)
    end

    return Term(1.0, m)
end

"""
    simplify(m::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}

Simplify a non-commutative algebra monomial with site-aware commutation.

Non-mutating version - creates a copy and simplifies it.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{NonCommutativeAlgebra}([idx1_s2, idx1_s1]);

julia> t = simplify(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
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

"""
    Base.:*(m1::Monomial{NonCommutativeAlgebra,T}, m2::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}

Multiply two non-commutative monomials with site-aware simplification.

Site-encoded operators on different sites commute.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m1 = Monomial{NonCommutativeAlgebra}([idx1_s1]);

julia> m2 = Monomial{NonCommutativeAlgebra}([idx1_s2]);

julia> t = m1 * m2;

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true
```
"""
function Base.:*(m1::Monomial{NonCommutativeAlgebra,T}, m2::Monomial{NonCommutativeAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word

    # Handle empty cases
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Concatenate and simplify using site-aware simplify!
    result = Monomial{NonCommutativeAlgebra,T}(vcat(w1, w2), zero(UInt64))
    simplify!(result)
end
