# =============================================================================
# Basis Generation for Non-Commutative Monomials
# =============================================================================
#
# Generates bases of monomials for non-commutative polynomial optimization.
# Uses VariableRegistry for type-consistent basis generation.
#
# Design:
# - Registry-aware: uses VariableRegistry{A,T} for proper index types
# - Returns Vector{Term}: handles simplification results with coefficients
# - Algebra-dispatched: simplification is algebra-specific

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
julia> using FastPolynomials: _generate_all_words

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
julia> using FastPolynomials: _generate_all_words

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

Generate all terms of exactly degree d using variables from the registry.

Returns `Vector{Term}` to properly handle simplification results, since
some algebras (like Pauli) may produce coefficients other than 1.0 during
simplification.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Exact degree

# Returns
- `Vector{Term{Monomial{A,T}, ComplexF64}}`: All simplified terms of degree d

# Examples
```jldoctest
julia> using FastPolynomials

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis_deg(reg, 2);

julia> length(basis)  # 2^2 = 4 monomials
4

julia> all(t -> t isa Term, basis)
true
```

Degree 0 returns identity term:
```jldoctest
julia> using FastPolynomials

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

    # Negative degree: return empty
    d < 0 && return Term{Monomial{A,T}, ComplexF64}[]

    # Degree 0: return identity term
    if d == 0
        identity_mono = Monomial{A}(T[])
        return [Term(one(ComplexF64), identity_mono)]
    end

    # Generate all words of length d using registry indices
    all_words = _generate_all_words(idxs, d)

    # Use ComplexF64 as coefficient type to handle all algebras
    # (PauliAlgebra returns ComplexF64, others return Float64 which promotes)
    terms = Term{Monomial{A,T}, ComplexF64}[]
    for word in all_words
        mono = Monomial{A}(word)
        simplified = simplify(mono)
        # Handle both Term and Vector{Term} return types from simplify
        if simplified isa Vector
            for term in simplified
                push!(terms, Term(ComplexF64(term.coefficient), term.monomial))
            end
        else
            push!(terms, Term(ComplexF64(simplified.coefficient), simplified.monomial))
        end
    end

    return terms
end

"""
    get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}

Generate all terms up to and including degree d using variables from the registry.

Returns `Vector{Term}` to properly handle simplification results.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Maximum degree (inclusive)

# Returns
- `Vector{Term{Monomial{A,T}, ComplexF64}}`: All simplified terms from degree 0 to d

# Examples
```jldoctest
julia> using FastPolynomials

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> length(basis)  # 1 + 2 + 4 = 7
7

julia> all(t -> t isa Term, basis)
true
```

With unipotent algebra (U^2 = I simplification):
```jldoctest
julia> using FastPolynomials

julia> reg, (U,) = create_unipotent_variables([("U", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> all(t -> t isa Term, basis)
true
```
"""
function get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}
    terms = Term{Monomial{A,T}, ComplexF64}[]
    for deg in 0:d
        append!(terms, get_ncbasis_deg(registry, deg))
    end
    return terms
end
