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
function sorted_union(args...)
    isempty(args) && return []
    result = unique(vcat(args...))
    sort!(result)
    return result
end

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
function sorted_unique(v::Vector)
    isempty(v) && return similar(v, 0)
    result = unique(v)
    sort!(result)
    return result
end

"""
    _neat_dot3(a::NCStateWord, m::M, b::NCStateWord) where M

Compute adjoint(a) * m * b for NCStateWords with a middle monomial/polynomial.
This is a common pattern in moment matrix construction.

The middle element `m` can be:
- A Monomial: converted to NCStateWord and multiplied
- An NCStateWord: direct multiplication
- A StateWord: converted to NCStateWord and multiplied

# Returns
An NCStateWord representing the triple product.
"""
function _neat_dot3(a::NCStateWord{ST,A,T}, m::Monomial{A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # adjoint(a) * m * b
    # Convert m to NCStateWord for multiplication
    m_ncsw = NCStateWord(one(StateWord{ST,A,T}), m)
    return adjoint(a) * m_ncsw * b
end

function _neat_dot3(a::NCStateWord{ST,A,T}, m::NCStateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    return adjoint(a) * m * b
end

function _neat_dot3(a::NCStateWord{ST,A,T}, sw::StateWord{ST,A,T}, b::NCStateWord{ST,A,T}) where {ST<:StateType,A<:AlgebraType,T<:Integer}
    # Convert StateWord to NCStateWord with identity nc_word
    m_ncsw = NCStateWord(sw)
    return adjoint(a) * m_ncsw * b
end

# Overload for regular Monomials (non-state context)
"""
    neat_dot(a::Monomial, b::Monomial) -> Monomial

Compute adjoint(a) * b for regular Monomials by concatenating words.

Returns a Monomial with the adjoint of a's word followed by b's word.
Does NOT apply simplification - callers should simplify! explicitly if needed.

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
    # adjoint(a) * b - concatenate adjoint(a).word with b.word
    adjoint(a) * b
end

"""
    _neat_dot3(a::Monomial, m::Monomial, b::Monomial) -> Polynomial

Compute adjoint(a) * m * b for regular Monomials with algebra-specific simplification.

Returns a Polynomial containing the simplified result. For algebras with non-trivial
simplification rules (e.g., Pauli algebra where σ² = I and σₓσᵧ = iσz), the
simplification is applied automatically.

This is the three-argument form commonly used in moment matrix construction
where we need adjoint(row_index) * constraint_monomial * column_index.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{NonCommutativeAlgebra}(UInt16[1]);

julia> m2 = Monomial{NonCommutativeAlgebra}(UInt16[2]);

julia> m3 = Monomial{NonCommutativeAlgebra}(UInt16[3]);

julia> result = _neat_dot3(m1, m2, m3);

julia> monomials(result)[1].word
3-element Vector{UInt16}:
 0x0001
 0x0002
 0x0003
```

For Pauli algebra, simplification is applied:
```julia
julia> a = Monomial{PauliAlgebra}(UInt8[1]);  # σx₁

julia> result = _neat_dot3(a, one(a), a);  # σx₁ * I * σx₁ = I

julia> isone(monomials(result)[1])  # Result is identity
true
```
"""
function _neat_dot3(a::Monomial{A,T}, m::Monomial{A,T}, b::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # adjoint(a) * m * b - concatenate words then simplify
    concatenated = adjoint(a) * m * b
    simplified_term = simplify(concatenated)
    return Polynomial(simplified_term)
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
