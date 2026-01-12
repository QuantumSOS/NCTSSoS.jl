# =============================================================================
# Utility Functions for NCTSSoS Integration
# =============================================================================
#
# These utility functions are used throughout NCTSSoS for polynomial algebra
# operations. They provide efficient helpers for sorted unions, dot products,
# and fermionic/bosonic operator handling.
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
    sorted_unique!(v::Vector) -> Vector

Return a sorted vector containing only unique elements from the input.
Modifies `v` in-place.

# Examples
```jldoctest
julia> sorted_unique!([3, 1, 2, 1, 3])
3-element Vector{Int64}:
 1
 2
 3
```
"""
sorted_unique!(v::Vector) = sort!(unique!(v))

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

    # Compute nc_word part via _neat_dot3 on NormalMonomials (no simplification)
    nc_mono = _neat_dot3(a.nc_word, m.nc_word, b.nc_word)

    return NCStateWord(sw_prod, nc_mono)
end

# Overload for regular Monomials (non-state context)
"""
    neat_dot(a::Vector{T}, b::Vector{T}) -> Vector{T}

Compute the raw word for adjoint(a) * b by concatenating words.

Returns a `Vector{T}` containing the adjoint of a's word followed by b's word.
The result is NOT simplified - callers must use `simplify(A, result)` if they need
a canonical `Monomial`.

For unsigned types (e.g., Pauli): adjoint just reverses the word.
For signed types (e.g., Fermionic): adjoint reverses and negates (creation ↔ annihilation).

# Examples
```jldoctest
julia> using NCTSSoS: neat_dot

julia> neat_dot(UInt16[1, 2], UInt16[3, 4])  # Unsigned: reversed
4-element Vector{UInt16}:
 0x0002
 0x0001
 0x0003
 0x0004

julia> neat_dot(Int32[1, 2], Int32[3, 4])  # Signed: reversed + negated
4-element Vector{Int32}:
 -2
 -1
  3
  4
```
"""
function neat_dot(a::Vector{T}, b::Vector{T}) where {T<:Unsigned}
    return append!(reverse(a), copy(b))
end

function neat_dot(a::Vector{T}, b::Vector{T}) where {T<:Signed}
    res = similar(a, length(a) + length(b))
    res[1:length(a)] .= .-@view(a[end:-1:1])
    res[length(a)+1:end] .= b
    return res
end

"""
    _neat_dot3(a::Vector{T}, m::Vector{T}, b::Vector{T}) -> Vector{T}

Compute the raw word for adjoint(a) * m * b by concatenating words.

This is the three-argument form commonly used in moment matrix construction
where we need adjoint(row_index) * constraint_monomial * column_index.

Returns a `Vector{T}` - callers must use `simplify(A, result)` if they need
a canonical `Monomial`.

Currently only implemented for unsigned integer types (e.g., Pauli algebra).

# Examples
```jldoctest
julia> using NCTSSoS: _neat_dot3

julia> _neat_dot3(UInt16[1, 2], UInt16[3], UInt16[4, 5])
5-element Vector{UInt16}:
 0x0002
 0x0001
 0x0003
 0x0004
 0x0005
```
"""
function _neat_dot3(
    a::Vector{T}, m::Vector{T}, b::Vector{T}
) where {T<:Unsigned}
    res = similar(a, length(a)+length(m)+length(b))
    res[1:length(a)] .= @view(a[end:-1:1])
    res[length(a)+1:length(a)+length(m)] .= m
    res[length(a)+length(m)+1:end] .= b
    return res
end

"""
    _neat_dot3(a::Vector{T}, m::Vector{T}, b::Vector{T}) where {T<:Signed} -> Vector{T}

Signed version for Fermionic/Bosonic algebras.
Adjoint reverses the word and negates operators (creation ↔ annihilation).
"""
function _neat_dot3(
    a::Vector{T}, m::Vector{T}, b::Vector{T}
) where {T<:Signed}
    res = similar(a, length(a)+length(m)+length(b))
    res[1:length(a)] .= .-@view(a[end:-1:1])  # adjoint: reverse + negate
    res[length(a)+1:length(a)+length(m)] .= m
    res[length(a)+length(m)+1:end] .= b
    return res
end

"""
    _neat_dot3(a::NormalMonomial{A,T}, m::NormalMonomial{A,T}, b::NormalMonomial{A,T}) where {A,T} -> Vector{T}

Convenience method that extracts words from NormalMonomials and delegates to the Vector version.
Returns the raw word Vector{T} for the caller to wrap in simplify().
"""
function _neat_dot3(
    a::NormalMonomial{A,T}, m::NormalMonomial{A,T}, b::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    return _neat_dot3(a.word, m.word, b.word)
end

# =============================================================================
# Shared Helper Functions for Fermionic/Bosonic Algebras
# =============================================================================

"""
    _is_creation(op::Integer) -> Bool

Check if an operator index represents a creation operator.
In the fermionic/bosonic encoding, creation operators have negative indices.

# Examples
```jldoctest
julia> using NCTSSoS: _is_creation

julia> _is_creation(-1)  # a₁† (creation)
true

julia> _is_creation(1)   # a₁ (annihilation)
false
```
"""
@inline _is_creation(op::T) where {T<:Signed} = op < 0

"""
    _operator_mode(op::T) where T<:Integer -> T

Extract the mode (site) index from a fermionic or bosonic operator.
The mode is the absolute value of the operator index.
Returns the same integer type as the input for type stability.

# Examples
```jldoctest
julia> using NCTSSoS: _operator_mode

julia> _operator_mode(-3)  # a₃† → mode 3
3

julia> _operator_mode(2)   # a₂ → mode 2
2
```
"""
@inline _operator_mode(op::T) where {T<:Signed} = abs(op)

"""
    normal_order_key(op::T) where T<:Integer -> Tuple{Int, T}

Compute the sort key for normal ordering of fermionic/bosonic operators.

Canonical order:
- Creation operators (negative) on the LEFT, sorted by mode in **descending** order
- Annihilation operators (positive) on the RIGHT, sorted by mode in ascending order

This choice makes `adjoint(m)` preserve the normal-form invariant without
introducing extra signs from re-sorting.

# Examples
```jldoctest
julia> using NCTSSoS: normal_order_key

julia> normal_order_key(-3)  # a₃† (creation)
(0, -3)

julia> normal_order_key(2)   # a₂ (annihilation)
(1, 2)
```
"""
@inline normal_order_key(op::T) where {T<:Signed} = (!_is_creation(op), op)

"""
    find_first_out_of_order(word::AbstractVector{T}) where T<:Integer -> Int

Find the first position where operators are out of normal order.
Returns the index `i` where `word[i]` and `word[i+1]` are out of order,
or 0 if the word is already in normal order.

Normal order means: all creators (negative) before annihilators (positive),
with creators sorted by mode descending and annihilators by mode ascending.

# Examples
```jldoctest
julia> using NCTSSoS: find_first_out_of_order

julia> find_first_out_of_order(Int32[-1, 1])  # c₁† c₁ (normal order)
0

julia> find_first_out_of_order(Int32[1, -1])  # c₁ c₁† (out of order at position 1)
1

julia> find_first_out_of_order(Int32[-1, -2]) # c₁† c₂† (out of order at position 1)
1
```
"""
function find_first_out_of_order(word::AbstractVector{T}) where {T<:Signed}
    @inbounds for i in 1:length(word)-1
        key_i = normal_order_key(word[i])
        key_i1 = normal_order_key(word[i+1])
        key_i > key_i1 && return i
    end
    return 0
end

"""
    is_normal_ordered(word::AbstractVector{T}) where T<:Integer -> Bool

Check if a word of fermionic/bosonic operators is in normal order.
Normal order means: all creators (negative) before annihilators (positive),
with creators sorted by mode descending and annihilators by mode ascending.

Implemented in terms of `find_first_out_of_order` for DRY.

# Examples
```jldoctest
julia> using NCTSSoS: is_normal_ordered

julia> is_normal_ordered(Int32[-2, -1, 1, 2])  # c₂† c₁† c₁ c₂ (normal)
true

julia> is_normal_ordered(Int32[1, -1])  # c₁ c₁† (not normal)
false
```
"""
is_normal_ordered(word::AbstractVector{T}) where {T<:Signed} = iszero(find_first_out_of_order(word))

@inline _sorted_symmetric_basis!(xs) = sorted_unique!(symmetric_canon.(xs))

@inline function _sorted_stateword_basis_from_ncsw!(xs)
    return sorted_unique!([symmetric_canon(expval(ncsw)) for ncsw in xs])
end

"""
    combine_like_terms(terms::Vector{Tuple{C,NormalMonomial{A,T}}}) where {A,T,C} -> Vector{Tuple{C,NormalMonomial{A,T}}}

Combine terms with identical monomials by summing their coefficients.
Filters out terms with zero coefficients.

    Internal utility for combining like terms by monomial word.
    Currently used by tests; kept as a general helper.

# Arguments
- `terms`: Vector of `(coefficient, monomial)` pairs to combine

# Returns
A new vector with like terms combined. If all terms cancel to zero,
returns a vector with a single zero term.
"""
function combine_like_terms(
    terms::Vector{Tuple{C,NormalMonomial{A,T}}}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Use inner constructor since we're reconstructing from existing valid words
    isempty(terms) && return [(zero(C), one(NormalMonomial{A,T}))]

    grouped = Dict{NormalMonomial{A,T},C}()
    for (coef, mono) in terms
        grouped[mono] = get(grouped, mono, zero(C)) + coef
    end

    result = Tuple{C,NormalMonomial{A,T}}[]
    for (m, coef) in grouped
        iszero(coef) && continue
        push!(result, (coef, m))
    end

    if isempty(result)
        push!(result, (zero(C), one(NormalMonomial{A,T})))
    end

    return result
end

# =============================================================================
# Pauli Algebra Helpers
# =============================================================================

"""
    _pauli_site(idx::Integer) -> Int

Extract the site number from a Pauli operator index.
Uses the encoding: site = (idx - 1) ÷ 3 + 1

# Examples
```jldoctest
julia> using NCTSSoS: _pauli_site

julia> _pauli_site(1)  # σₓ at site 1
1

julia> _pauli_site(4)  # σₓ at site 2
2

julia> _pauli_site(7)  # σₓ at site 3
3
```
"""
@inline _pauli_site(idx::Integer) = (idx - 1) ÷ 3 + 1

"""
    _pauli_type(idx::Integer) -> Int

Extract the Pauli type from an operator index.
Returns 0 for σₓ, 1 for σᵧ, 2 for σz.

# Examples
```jldoctest
julia> using NCTSSoS: _pauli_type

julia> _pauli_type(1)  # σₓ at site 1
0

julia> _pauli_type(2)  # σᵧ at site 1
1

julia> _pauli_type(3)  # σz at site 1
2
```
"""
@inline _pauli_type(idx::Integer) = (idx - 1) % 3

"""
    _pauli_index(site::Integer, type::Integer) -> Int

Compute the Pauli operator index from site and type.
Inverse of `_pauli_site` and `_pauli_type`.

# Examples
```jldoctest
julia> using NCTSSoS: _pauli_index

julia> _pauli_index(1, 0)  # σₓ at site 1
1

julia> _pauli_index(2, 1)  # σᵧ at site 2
5

julia> _pauli_index(3, 2)  # σz at site 3
9
```
"""
@inline _pauli_index(site::Integer, type::Integer) = (site - 1) * 3 + type + 1

# =============================================================================
# Phase Conversion
# =============================================================================

"""
    _phase_to_complex(phase::UInt8) -> ComplexF64

Convert a phase encoding to a complex number.
phase ∈ 0:3 maps to (im)^phase: 0→1, 1→i, 2→-1, 3→-i

Used in Pauli algebra where products accumulate phases as powers of i.

# Examples
```jldoctest
julia> using NCTSSoS: _phase_to_complex

julia> _phase_to_complex(UInt8(0))
1.0 + 0.0im

julia> _phase_to_complex(UInt8(1))
0.0 + 1.0im

julia> _phase_to_complex(UInt8(2))
-1.0 + 0.0im
```
"""
@inline function _phase_to_complex(phase::UInt8)
    phase == 0 && return ComplexF64(1.0, 0.0)
    phase == 1 && return ComplexF64(0.0, 1.0)
    phase == 2 && return ComplexF64(-1.0, 0.0)
    return ComplexF64(0.0, -1.0)  # phase == 3
end
