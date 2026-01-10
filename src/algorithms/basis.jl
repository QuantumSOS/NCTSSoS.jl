# =============================================================================
# Basis Generation for Non-Commutative Monomials
# =============================================================================
#
# Generates bases of monomials for non-commutative polynomial optimization.
# Uses VariableRegistry for type-consistent basis generation.
#
# Design:
# - Registry-aware: uses VariableRegistry{A,T} for proper index types
# - User-facing basis elements are `Monomial` (single-term, unit internal coefficient)
# - Algebra-dispatched: PBW expands; monoids stay one-monomial

"""
    _generate_all_words(indices::Vector{T}, d::Int) where {T<:Integer} -> Vector{Vector{T}}

Internal helper to generate all words of exactly length d using the provided indices.

# Arguments
- `indices`: Vector of variable indices to use
- `d`: Exact word length (degree)

# Returns
Vector of all possible words of length d using the given indices.

# Examples
```jldoctest
julia> using NCTSSoS: _generate_all_words

julia> words = _generate_all_words([1, 3, 5], 2);

julia> length(words)  # 3^2 = 9
9

julia> [1, 1] in words
true

julia> [3, 5] in words
true
```

Edge cases:
```jldoctest
julia> using NCTSSoS: _generate_all_words

julia> _generate_all_words([1, 2], 0)  # degree 0
1-element Vector{Vector{Int64}}:
 []

julia> _generate_all_words(Int[], 1)  # no indices
Vector{Int64}[]
```
"""
function _generate_all_words(indices::Vector{T}, d::Int) where {T<:Integer}
    d == 0 && return [T[]]
    isempty(indices) && return Vector{T}[]
    d == 1 && return [[idx] for idx in indices]

    result = Vector{T}[]
    for idx in indices
        for suffix in _generate_all_words(indices, d - 1)
            push!(result, vcat([idx], suffix))
        end
    end
    return result
end

"""
    get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}

Generate all simplified monomials of exactly degree d using variables from the registry.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Exact degree

# Returns
- `Vector{Monomial{A,T}}`: Canonical monomial *elements* produced by simplification.
  Each entry is a single-term `Monomial` with unit internal coefficient and a
  `NormalMonomial` word in normal form.

# Examples
```jldoctest
julia> using NCTSSoS

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis_deg(reg, 2);

julia> length(basis)  # 2^2 = 4 monomials
4

julia> all(m -> m isa Monomial, basis)
true
```

Degree 0 returns identity monomial:
```jldoctest
julia> using NCTSSoS

julia> reg, (x,) = create_noncommutative_variables([("x", 1:3)]);

julia> basis = get_ncbasis_deg(reg, 0);

julia> length(basis)
1

julia> isone(basis[1])
true
```
"""
function get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}
    idxs = indices(registry)

    monos = _get_ncbasis_deg(A, T, idxs, d)
    return [Monomial(m) for m in monos]
end

function _get_ncbasis_deg(
    ::Type{A}, ::Type{T}, idxs::Vector{T}, d::Int
) where {A<:PBWAlgebra, T<:Integer}
    d < 0 && return NormalMonomial{A,T}[]
    d == 0 && return [one(NormalMonomial{A,T})]

    all_words = _generate_all_words(idxs, d)
    result = NormalMonomial{A,T}[]

    for word in all_words
        # Use inner constructor to allow raw (non-canonical) words; simplify() canonicalizes.
        raw = NormalMonomial{A,T}(word, _UNSAFE_NORMAL_MONOMIAL)
        for (_, mono) in simplify(raw)
            push!(result, mono)
        end
    end

    unique!(sort!(result))
    return result
end

function _get_ncbasis_deg(
    ::Type{A}, ::Type{T}, idxs::Vector{T}, d::Int
) where {A<:Union{MonoidAlgebra,TwistedGroupAlgebra}, T<:Integer}
    d < 0 && return NormalMonomial{A,T}[]
    d == 0 && return [one(NormalMonomial{A,T})]

    all_words = _generate_all_words(idxs, d)
    result = NormalMonomial{A,T}[]

    for word in all_words
        # Monoid/TwistedGroup algebras simplify to a single monomial (up to an internal coefficient).
        raw = NormalMonomial{A,T}(word, _UNSAFE_NORMAL_MONOMIAL)
        for (_, mono) in simplify(raw)
            push!(result, mono)
        end
    end

    return result
end

"""
    get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}

Generate all simplified monomials up to and including degree d using variables from the registry.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Maximum degree (inclusive)

# Returns
- `Vector{Monomial{A,T}}`: Canonical monomial elements (single-term) appearing in simplification results

# Examples
```jldoctest
julia> using NCTSSoS

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> length(basis)  # 1 + 2 + 4 = 7
7

julia> all(m -> m isa Monomial, basis)
true
```

With unipotent algebra (U^2 = I simplification):
```jldoctest
julia> using NCTSSoS

julia> reg, (U,) = create_unipotent_variables([("U", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> all(m -> m isa Monomial, basis)
true
```
"""
function get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}
    idxs = indices(registry)

    result = NormalMonomial{A,T}[]
    for deg in 0:d
        append!(result, _get_ncbasis_deg(A, T, idxs, deg))
    end
    unique!(sort!(result))
    return [Monomial(m) for m in result]
end
