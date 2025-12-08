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
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{UnipotentAlgebra}([idx1_s2, idx1_s1, idx1_s2]);

julia> t = simplify(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1]
true
```

U²=I pair removal:
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> m = Monomial{UnipotentAlgebra}([idx1_s1, idx1_s1]);

julia> t = simplify(m);

julia> isempty(t.monomial.word)
true
```
"""

"""
    simplify!(m::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}

Site-aware in-place simplification for unipotent algebra with encoded indices.

Operators on different sites commute and are sorted by site (ascending).
Within each site, U²=I applies: consecutive identical operators cancel (remove pairs).
Order of different operators within the same site is preserved (non-commutative within site).

# Algorithm
1. Group operators by site (using `decode_site`)
2. Sort groups by site (ascending)
3. Within each site: apply U²=I via stack (remove consecutive pairs)
4. Concatenate sorted groups

# Warning
This mutates the input monomial. Use `simplify` for a non-mutating version.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m = Monomial{UnipotentAlgebra}([idx1_s2, idx1_s1]);

julia> t = simplify!(m);

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true

julia> m.word == [idx1_s1, idx1_s2]  # Original was mutated
true
```
"""
function simplify!(m::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}
    word = m.word

    # Empty or single: nothing to simplify
    length(word) <= 1 && return Term(1.0, m)

    # Stable sort by site (operators on different sites commute, within-site order preserved)
    sort!(word; alg=InsertionSort, by=decode_site)

    # Apply U²=I: remove consecutive identical pairs with backtracking
    i = 1
    while i < length(word)
        if word[i] == word[i + 1]
            # Consecutive identical: remove both (U² = I)
            deleteat!(word, i)
            deleteat!(word, i)
            # Backtrack to catch cascading cancellations
            i > 1 && (i -= 1)
        else
            i += 1
        end
    end

    return Term(1.0, m)
end

"""
    simplify(m::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}

Simplify a unipotent algebra monomial with site-aware commutation and U²=I.

Non-mutating version - creates a copy and simplifies it.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> m = Monomial{UnipotentAlgebra}([idx1_s1, idx1_s1]);

julia> t = simplify(m);

julia> t.coefficient
1.0

julia> isempty(t.monomial.word)
true

julia> length(m.word)  # Original unchanged
2
```
"""
function simplify(m::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}
    # Copy and delegate to simplify!
    m_copy = Monomial{UnipotentAlgebra,T}(copy(m.word), m.hash)
    return simplify!(m_copy)
end
