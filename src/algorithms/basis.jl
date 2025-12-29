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

Generate all simplified polynomials of exactly degree d using variables from the registry.

Returns `Vector{Polynomial}` where each element is the simplified form of one input monomial.
This preserves the 1-to-1 mapping between input words and output polynomials:
- NonCommutativeAlgebra: each polynomial is a single monomial (no simplification)
- PauliAlgebra: each polynomial is a weighted monomial (coefficient from simplification)
- FermionicAlgebra/BosonicAlgebra: each polynomial may have multiple terms

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Exact degree

# Returns
- `Vector{Polynomial{A,T,ComplexF64}}`: Simplified polynomials, one per input word

# Examples
```jldoctest
julia> using NCTSSoS

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis_deg(reg, 2);

julia> length(basis)  # 2^2 = 4 monomials
4

julia> all(p -> p isa Polynomial, basis)
true
```

Degree 0 returns identity polynomial:
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

    # Negative degree: return empty
    d < 0 && return Polynomial{A,T,ComplexF64}[]

    # Degree 0: return identity polynomial
    if d == 0
        identity_mono = Monomial{A}(T[])
        identity_term = Term(one(ComplexF64), identity_mono)
        return [Polynomial([identity_term])]
    end

    # Generate all words of length d using registry indices
    all_words = _generate_all_words(idxs, d)

    # Each word becomes one polynomial (the simplified form of that monomial)
    result = Polynomial{A,T,ComplexF64}[]
    for word in all_words
        mono = Monomial{A}(word)
        simplified = simplify(mono)
        # Handle different return types from simplify:
        # - Monomial: NC, Projector, Unipotent (coefficient is implicitly 1.0)
        # - Term: Pauli (has coefficient)
        # - Polynomial: Bosonic, Fermionic (multiple terms possible)
        if simplified isa Monomial
            term = Term(ComplexF64(1.0), simplified)
            push!(result, Polynomial([term]))
        elseif simplified isa Term
            term = Term(ComplexF64(simplified.coefficient), simplified.monomial)
            push!(result, Polynomial([term]))
        elseif simplified isa Polynomial
            # Convert to ComplexF64 coefficients
            terms = [Term(ComplexF64(t.coefficient), t.monomial) for t in simplified.terms]
            push!(result, Polynomial(terms))
        else
            # Legacy: handle Vector{Term} just in case
            terms = [Term(ComplexF64(t.coefficient), t.monomial) for t in simplified]
            push!(result, Polynomial(terms))
        end
    end

    return result
end

"""
    get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}

Generate all simplified polynomials up to and including degree d using variables from the registry.

Returns `Vector{Polynomial}` where each element is the simplified form of one input monomial.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Maximum degree (inclusive)

# Returns
- `Vector{Polynomial{A,T,ComplexF64}}`: Simplified polynomials from degree 0 to d

# Examples
```jldoctest
julia> using NCTSSoS

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> length(basis)  # 1 + 2 + 4 = 7
7

julia> all(p -> p isa Polynomial, basis)
true
```

With unipotent algebra (U^2 = I simplification):
```jldoctest
julia> using NCTSSoS

julia> reg, (U,) = create_unipotent_variables([("U", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> all(p -> p isa Polynomial, basis)
true
```
"""
function get_ncbasis(registry::VariableRegistry{A,T}, d::Int) where {A<:AlgebraType, T<:Integer}
    result = Polynomial{A,T,ComplexF64}[]
    for deg in 0:d
        append!(result, get_ncbasis_deg(registry, deg))
    end
    return result
end
