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

    # Build result: for each site, apply U²=I (remove consecutive pairs via stack)
    empty!(word)

    for site in sorted_sites
        ops = site_groups[site]
        if !isempty(ops)
            # Stack-based approach for U²=I within site
            stack = T[]
            for op in ops
                if !isempty(stack) && stack[end] == op
                    # Same as stack top: pop (U² = I)
                    pop!(stack)
                else
                    # Different: push
                    push!(stack, op)
                end
            end
            append!(word, stack)
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
    simplify!(m_copy)
end

"""
    Base.:*(m1::Monomial{UnipotentAlgebra,T}, m2::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}

Multiply two unipotent monomials with site-aware simplification.

Site-encoded operators on different sites commute. U²=I applies within sites.

# Examples
```jldoctest
julia> using FastPolynomials

julia> using FastPolynomials: encode_index

julia> idx1_s1 = encode_index(UInt16, 1, 1);

julia> idx1_s2 = encode_index(UInt16, 1, 2);

julia> m1 = Monomial{UnipotentAlgebra}([idx1_s1]);

julia> m2 = Monomial{UnipotentAlgebra}([idx1_s2]);

julia> t = m1 * m2;

julia> t.coefficient
1.0

julia> t.monomial.word == [idx1_s1, idx1_s2]
true
```
"""
function Base.:*(m1::Monomial{UnipotentAlgebra,T}, m2::Monomial{UnipotentAlgebra,T}) where {T<:Unsigned}
    w1, w2 = m1.word, m2.word

    # Handle empty cases
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Concatenate and simplify using site-aware simplify!
    result = Monomial{UnipotentAlgebra,T}(vcat(w1, w2), zero(UInt64))
    simplify!(result)
end

"""
    Base.:*(m1::Monomial{UnipotentAlgebra,T}, m2::Monomial{UnipotentAlgebra,T}) where {T<:Signed}

Multiply two unipotent monomials with signed integer indices (legacy fallback).

For signed integers, we cannot use site-based encoding, so we fall back to
simple concatenation with consecutive duplicate removal (U²=I).

This is provided for backward compatibility with NCTSSoS legacy code that
uses `Int` indices instead of encoded unsigned indices.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{UnipotentAlgebra}(Int[1]);

julia> m2 = Monomial{UnipotentAlgebra}(Int[2]);

julia> t = m1 * m2;

julia> t.coefficient
1.0

julia> t.monomial.word == Int[1, 2]
true

julia> t2 = m1 * m1;  # U² = I

julia> isempty(t2.monomial.word)
true
```
"""
function Base.:*(m1::Monomial{UnipotentAlgebra,T}, m2::Monomial{UnipotentAlgebra,T}) where {T<:Signed}
    w1, w2 = m1.word, m2.word

    # Handle empty cases
    isempty(w1) && return Term(1.0, m2)
    isempty(w2) && return Term(1.0, m1)

    # Simple concatenation with U²=I (consecutive duplicate removal)
    result_word = T[]

    # Process w1
    for idx in w1
        if !isempty(result_word) && result_word[end] == idx
            pop!(result_word)  # U² = I
        else
            push!(result_word, idx)
        end
    end

    # Process w2
    for idx in w2
        if !isempty(result_word) && result_word[end] == idx
            pop!(result_word)  # U² = I
        else
            push!(result_word, idx)
        end
    end

    result = Monomial{UnipotentAlgebra,T}(result_word, zero(UInt64))
    return Term(1.0, result)
end

"""
    Base.adjoint(m::Monomial{UnipotentAlgebra,T}) where {T<:Signed}

Compute the adjoint of a unipotent monomial with signed indices.

Unipotent operators are self-adjoint (U = U†), so adjoint only reverses the word.
This overrides the default signed behavior which would also negate indices.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{UnipotentAlgebra}(Int[1, 2, 3]);

julia> adjoint(m).word
3-element Vector{Int64}:
 3
 2
 1
```
"""
function Base.adjoint(m::Monomial{UnipotentAlgebra,T}) where {T<:Signed}
    new_word = reverse(m.word)
    Monomial{UnipotentAlgebra}(new_word)
end
