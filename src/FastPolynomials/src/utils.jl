# =============================================================================
# Utility Functions for NCTSSoS Integration
# =============================================================================
#
# These utility functions are used by NCTSSoS and were part of the old
# FastPolynomials implementation. They are maintained here for backward
# compatibility with NCTSSoS source files.
#
# =============================================================================

"""
    sorted_union(args...) -> Vector

Compute the sorted union of multiple vectors. Returns a sorted vector
containing all unique elements from the input vectors.

# Examples
```jldoctest
julia> sorted_union([3, 1], [2, 1])
3-element Vector{Int64}:
 1
 2
 3

julia> sorted_union([1, 2], [2, 3], [3, 4])
4-element Vector{Int64}:
 1
 2
 3
 4
```
"""
sorted_union(xs...) = sort!(union(xs...))

"""
    sorted_unique(v::Vector) -> Vector

Return a sorted vector containing only unique elements from the input.

# Examples
```jldoctest
julia> sorted_unique([3, 1, 2, 1, 3])
3-element Vector{Int64}:
 1
 2
 3
```
"""
sorted_unique(v::Vector) = sort!(unique(v))

"""
    _neat_dot3(a::NCStateWord, m::NCStateWord, b::NCStateWord) -> NCStateWord

Compute adjoint(a) * m * b for NCStateWords.
This is a common pattern in moment matrix construction.

!!! note
    The returned NCStateWord is NOT simplified. Callers should call `simplify`
    explicitly if algebra-specific simplification is needed.

# Returns
An NCStateWord representing the triple product (unsimplified).
"""
function _neat_dot3(
    a::NCStateWord{ST,A,T}, m::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}
) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Concatenate state word parts directly (adjoint is identity due to Hermitian invariance)
    sw_prod = adjoint(a.sw) * m.sw * b.sw

    # Compute nc_word part via _neat_dot3 on Monomials (no simplification)
    nc_mono = _neat_dot3(a.nc_word, m.nc_word, b.nc_word)

    return NCStateWord(sw_prod, nc_mono)
end

"""
    simplify(ncsw::NCStateWord) -> NCStatePolynomial

Simplify an NCStateWord by applying algebra-specific simplification rules to its
nc_word (Monomial) part. Returns an NCStatePolynomial since simplification may
produce multiple terms (e.g., Pauli algebra phase factors).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{UnipotentAlgebra}(UInt[1, 1]);  # x₁²

julia> sw = StateWord{Arbitrary}([one(Monomial{UnipotentAlgebra,UInt})]);

julia> ncsw = NCStateWord(sw, m);

julia> result = simplify(ncsw);

julia> result isa NCStatePolynomial
true
```
"""
function simplify(ncsw::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Simplify the nc_word part (may produce multiple terms with phases)
    nc_poly = Polynomial(simplify(ncsw.nc_word))

    # Convert to NCStatePolynomial: each term gets the same StateWord
    coeffs = [t.coefficient for t in nc_poly.terms]
    ncsws = [NCStateWord(ncsw.sw, t.monomial) for t in nc_poly.terms]
    return NCStatePolynomial(coeffs, ncsws)
end

# Overload for regular Monomials (non-state context)
"""
    neat_dot(a::Monomial, b::Monomial) -> Monomial

Compute adjoint(a) * b for regular Monomials by concatenating words.

Returns a Monomial with the adjoint of a's word followed by b's word.

!!! note
    The returned Monomial is NOT simplified. Callers should call `simplify!`
    explicitly if algebra-specific simplification is needed.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3, 4]);

julia> result = neat_dot(m1, m2);

julia> result.word
4-element Vector{Int64}:
 -2
 -1
  3
  4
```
"""
function neat_dot(a::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # Preallocate result: reverse(a.word) ++ b.word (single allocation)
    n_a, n_b = length(a.word), length(b.word)
    result = Vector{T}(undef, n_a + n_b)

    # Copy reversed a.word (with negation for signed types)
    @inbounds for i in 1:n_a
        result[i] = T <: Signed ? -a.word[n_a - i + 1] : a.word[n_a - i + 1]
    end

    # Copy b.word
    @inbounds for i in 1:n_b
        result[n_a + i] = b.word[i]
    end

    return Monomial{A}(result)
end

"""
    _neat_dot3(a::Monomial, m::Monomial, b::Monomial) -> Monomial

Compute adjoint(a) * m * b for regular Monomials by concatenating words.

This is the three-argument form commonly used in moment matrix construction
where we need adjoint(row_index) * constraint_monomial * column_index.

!!! note
    The returned Monomial is NOT simplified. Callers should call `simplify`
    explicitly if algebra-specific simplification is needed.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> m3 = Monomial{PauliAlgebra}([4, 5]);

julia> result = _neat_dot3(m1, m2, m3);

julia> result.word
5-element Vector{Int64}:
 -2
 -1
  3
  4
  5
```
"""
function _neat_dot3(
    a::Monomial{A,T}, m::Monomial{A,T}, b::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    # Preallocate: reverse(a.word) ++ m.word ++ b.word (single allocation)
    n_a, n_m, n_b = length(a.word), length(m.word), length(b.word)
    result = Vector{T}(undef, n_a + n_m + n_b)

    # Copy reversed a.word (with negation for signed types)
    @inbounds for i in 1:n_a
        result[i] = T <: Signed ? -a.word[n_a - i + 1] : a.word[n_a - i + 1]
    end

    # Copy m.word
    @inbounds for i in 1:n_m
        result[n_a + i] = m.word[i]
    end

    # Copy b.word
    @inbounds for i in 1:n_b
        result[n_a + n_m + i] = b.word[i]
    end

    return Monomial{A}(result)
end
