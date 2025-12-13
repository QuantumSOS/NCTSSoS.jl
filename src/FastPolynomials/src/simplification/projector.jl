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
function simplify!(m::Monomial{ProjectorAlgebra,T}) where {T<:Unsigned}
    word = m.word

    # Empty or single: nothing to simplify
    length(word) <= 1 && return Term(1.0, m)

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

    # Create new monomial with correct hash (can't update hash in-place since Monomial is immutable)
    simplified_mono = Monomial{ProjectorAlgebra}(word)
    return Term(1.0, simplified_mono)
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
