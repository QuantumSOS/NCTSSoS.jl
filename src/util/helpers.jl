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

    # Compute nc_word part via _neat_dot3 on NormalMonomials (no simplification)
    nc_mono = _neat_dot3(a.nc_word, m.nc_word, b.nc_word)

    return NCStateWord(sw_prod, nc_mono)
end

# Overload for regular Monomials (non-state context)
"""
    neat_dot(a::NormalMonomial, b::NormalMonomial) -> NormalMonomial

Compute adjoint(a) * b for regular Monomials by concatenating words.

Returns a NormalMonomial with the adjoint of a's word followed by b's word.

!!! note
    The returned NormalMonomial is NOT simplified. Callers should call `simplify!`
    explicitly if algebra-specific simplification is needed.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1, 2]);

julia> m2 = NormalMonomial{PauliAlgebra}([3, 4]);

julia> result = neat_dot(m1, m2);

julia> result.word
4-element Vector{Int64}:
 -2
 -1
  3
  4
```
"""
function neat_dot(a::NormalMonomial{A,T}, b::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
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

    return NormalMonomial{A,T}(result, _UNSAFE_NORMAL_MONOMIAL)
end

"""
    _neat_dot3(a::NormalMonomial, m::NormalMonomial, b::NormalMonomial) -> NormalMonomial

Compute adjoint(a) * m * b for regular Monomials by concatenating words.

This is the three-argument form commonly used in moment matrix construction
where we need adjoint(row_index) * constraint_monomial * column_index.

!!! note
    The returned NormalMonomial is NOT simplified. Callers should call `simplify`
    explicitly if algebra-specific simplification is needed.

# Examples
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: _neat_dot3

julia> m1 = NormalMonomial{PauliAlgebra}([1, 2]);

julia> m2 = NormalMonomial{PauliAlgebra}([3]);

julia> m3 = NormalMonomial{PauliAlgebra}([4, 5]);

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
    a::NormalMonomial{A,T}, m::NormalMonomial{A,T}, b::NormalMonomial{A,T}
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

    return NormalMonomial{A,T}(result, _UNSAFE_NORMAL_MONOMIAL)
end

# =============================================================================
# NormalMonomial Validation
# =============================================================================

"""
    validate(m::NormalMonomial{A,T}) where {A,T}

Validate that a NormalMonomial is in canonical form for its algebra type.
Throws `ArgumentError` if invalid.

# Canonical Forms by Algebra
- `PauliAlgebra`: ≤1 operator per site, sorted by site
- `FermionicAlgebra`: normal-ordered (creators left, creators by mode descending)
- `BosonicAlgebra`: normal-ordered (creators left, creators by mode descending)
- `ProjectorAlgebra`: sorted by site, no consecutive identical (no P²)
- `UnipotentAlgebra`: sorted by site, no consecutive identical (no U²)
- `NonCommutativeAlgebra`: sorted by site

For raw (non-canonical) input, use:
- `simplify(PauliAlgebra, word)` for Pauli algebra (returns a `Monomial`, iterable as `(c_internal, NormalMonomial)` pairs)
- `simplify(FermionicAlgebra, word)` / `simplify(BosonicAlgebra, word)` for PBW algebras

# Examples
```julia
julia> m = NormalMonomial{PauliAlgebra}([1, 2]);  # σx₁ and σy₁ on same site!

julia> validate(m)  # throws ArgumentError
ERROR: ArgumentError: Pauli word has multiple operators on site 1
```
"""
function validate(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    _validate_word!(A, m.word)
    return m
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
@inline _is_creation(op::T) where {T<:Integer} = op < 0

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
@inline _operator_mode(op::T) where {T<:Integer} = abs(op)

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
@inline function normal_order_key(op::T) where {T<:Integer}
    mode = _operator_mode(op)
    return _is_creation(op) ? (0, -mode) : (1, mode)
end

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
function find_first_out_of_order(word::AbstractVector{T}) where {T<:Integer}
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
is_normal_ordered(word::AbstractVector{T}) where {T<:Integer} = find_first_out_of_order(word) == 0

@inline _sorted_symmetric_basis(xs) = sorted_unique(symmetric_canon.(xs))

@inline function _sorted_stateword_basis_from_ncsw(xs)
    return sorted_unique([symmetric_canon(expval(ncsw)) for ncsw in xs])
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
