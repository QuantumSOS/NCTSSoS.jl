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

julia> t = simplify(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true
```
"""

"""
    simplify!(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}

Site-aware in-place simplification for projector algebra with encoded indices.

Operators on different sites commute and are sorted by site (ascending).
Within each site, idempotency applies: consecutive identical operators collapse (P² = P).
Order of different operators within the same site is preserved (non-commutative within site).

# Algorithm
1. Group operators by site (using `decode_site`)
2. Sort groups by site (ascending)
3. Within each site: remove consecutive duplicates (idempotency)
4. Concatenate sorted groups

# Warning
This mutates the input monomial. Use `simplify` for a non-mutating version.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{ProjectorAlgebra}([idx1_s2, idx1_s1]);

julia> t = simplify!(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true

julia> m.word == [idx1_s1, idx1_s2]  # Original was mutated
true
```
"""
# TODO: need to look at NCTSSOS this algorithm seems bad
function simplify!(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}
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

    # Build result: for each site, apply idempotency (remove consecutive duplicates)
    empty!(word)

    for site in sorted_sites
        ops = site_groups[site]
        if !isempty(ops)
            prev = ops[1]
            push!(word, prev)
            for i in 2:length(ops)
                if ops[i] != prev
                    prev = ops[i]
                    push!(word, prev)
                end
            end
        end
    end

    return Term(1.0, m)
end

"""
    simplify(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}

Simplify a projector algebra monomial with site-aware commutation and idempotency.

Non-mutating version - creates a copy and simplifies it.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> m = Monomial{ProjectorAlgebra}([idx1_s1, idx1_s1, idx1_s1]);

julia> t = simplify(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1]
true

julia> length(m.word)  # Original unchanged
3
```
"""
function simplify(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}
    # Copy and delegate to simplify!
    m_copy = Monomial{ProjectorAlgebra,T}(copy(m.word), m.hash)
    simplify!(m_copy)
end

"""
    Base.:*(m1::Monomial{ProjectorAlgebra,T}, m2::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}

Multiply two projector monomials with site-aware simplification.

Site-encoded operators on different sites commute. Idempotency applies within sites.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m1 = Monomial{ProjectorAlgebra}([idx1_s1]);

julia> m2 = Monomial{ProjectorAlgebra}([idx1_s2]);

julia> t = m1 * m2;

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true
```
"""
function Base.:*(m1::Monomial{ProjectorAlgebra,T}, m2::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word

    # Handle empty cases
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Concatenate and simplify using site-aware simplify!
    result = Monomial{ProjectorAlgebra,T}(vcat(w1, w2), zero(UInt64))
    simplify!(result)
end
