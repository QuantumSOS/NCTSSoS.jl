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
    _neat_dot3(a::NCStateWord, m::NCStateWord, b::NCStateWord) -> NCStatePolynomial

Compute adjoint(a) * m * b for NCStateWords.
This is a common pattern in moment matrix construction.

# Returns
An NCStatePolynomial representing the simplified triple product.
The nc_word parts are simplified using algebra-specific rules (e.g., x²=I for unipotent).
"""
function _neat_dot3(a::NCStateWord{ST,A,T}, m::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Compute state word part: adjoint(a.sw) * m.sw * b.sw (commutative)
    sw_prod = adjoint(a.sw) * m.sw * b.sw

    # Compute nc_word part via _neat_dot3 on Monomials, then simplify
    nc_mono = _neat_dot3(a.nc_word, m.nc_word, b.nc_word)
    nc_poly = Polynomial(simplify(nc_mono))

    # Convert to NCStatePolynomial
    coeffs = [t.coefficient for t in nc_poly.terms]
    ncsws = [NCStateWord(sw_prod, t.monomial) for t in nc_poly.terms]
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

    Monomial{A}(result)
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

julia> m1 = Monomial{NonCommutativeAlgebra}(UInt16[1]);

julia> m2 = Monomial{NonCommutativeAlgebra}(UInt16[2]);

julia> m3 = Monomial{NonCommutativeAlgebra}(UInt16[3]);

julia> result = _neat_dot3(m1, m2, m3);

julia> result.word
3-element Vector{UInt16}:
 0x0001
 0x0002
 0x0003
```
"""
function _neat_dot3(a::Monomial{A,T}, m::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
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

    Monomial{A}(result)
end

# =============================================================================
# Polynomial to StatePolynomial Conversion (kept from legacy API)
# =============================================================================

"""
    ς(p::Polynomial{A,T,C}) -> StatePolynomial

Create a StatePolynomial from a Polynomial.
Converts each term's monomial to a StateWord{Arbitrary}.
"""
function ς(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return StatePolynomial(C[], StateWord{Arbitrary,A,T}[])

    state_words = [StateWord{Arbitrary}(t.monomial) for t in p.terms]
    coeffs = C[t.coefficient for t in p.terms]
    StatePolynomial(coeffs, state_words)
end

# =============================================================================
# Monomial Identity Fallback
# =============================================================================

"""
    Base.one(::Type{Monomial}) -> Monomial{NonCommutativeAlgebra,UInt}

Create the identity monomial (empty word) for the generic Monomial type.
This fallback is needed for code that uses `one(Monomial)` without type parameters.
"""
Base.one(::Type{Monomial}) = Monomial{NonCommutativeAlgebra}(UInt[])

# =============================================================================
# Internal AbstractPolynomial Type Alias
# =============================================================================

"""
    AbstractPolynomial{T}

Internal type alias used by NCTSSoS for type constraints.
Maps to Union of Polynomial types with any algebra type.

This is NOT exported - for internal use only.
"""
const AbstractPolynomial{T} = Union{
    Polynomial{<:AlgebraType,<:Integer,T},
    StatePolynomial{T,<:StateType,<:AlgebraType,<:Integer},
    NCStatePolynomial{T,<:StateType,<:AlgebraType,<:Integer}
}
