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
julia> using FastPolynomials

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis_deg(reg, 2);

julia> length(basis)  # 2^2 = 4 monomials
4

julia> all(p -> p isa Polynomial, basis)
true
```

Degree 0 returns identity polynomial:
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
julia> using FastPolynomials

julia> reg, (x,) = create_noncommutative_variables([("x", 1:2)]);

julia> basis = get_ncbasis(reg, 2);

julia> length(basis)  # 1 + 2 + 4 = 7
7

julia> all(p -> p isa Polynomial, basis)
true
```

With unipotent algebra (U^2 = I simplification):
```jldoctest
julia> using FastPolynomials

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

# =============================================================================
# Normal-Form Basis Generation for Fermionic/Bosonic Algebras
# =============================================================================
#
# For fermionic and bosonic algebras, we generate bases containing ONLY
# normal-form monomials. Normal form means:
#   - All creation operators (negative indices) on the LEFT
#   - All annihilation operators (positive indices) on the RIGHT
#   - Within each group: sorted by mode (absolute value) ascending
#
# This is more efficient than generating all words and simplifying:
#   - Fermionic: O(2^n) subsets vs O((2n)^d) all words
#   - Bosonic: O(multisets) vs O((2n)^d) all words
#
# Returns Vector{Monomial} (not Vector{Polynomial}) since normal-form
# basis elements ARE monomials by definition.

"""
    _get_modes(registry::VariableRegistry{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra}, T<:Signed}

Extract unique modes from a fermionic/bosonic registry.

In these algebras, indices are signed: positive for annihilation, negative for creation.
Both ±i refer to the same physical mode i.

# Returns
Sorted vector of unique mode indices (positive integers).
"""
function _get_modes(registry::VariableRegistry{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra}, T<:Signed}
    return sort!(unique!(abs.(indices(registry))))
end

"""
    _combinations(v::Vector{T}, k::Int) where T

Generate all k-combinations (subsets of size k) from vector v.

Returns vector of vectors, each containing k elements from v in order.
"""
function _combinations(v::Vector{T}, k::Int) where {T}
    k == 0 && return [T[]]
    k > length(v) && return Vector{T}[]
    k == length(v) && return [copy(v)]

    result = Vector{T}[]
    for (i, elem) in enumerate(v)
        for suffix in _combinations(v[i+1:end], k - 1)
            push!(result, vcat([elem], suffix))
        end
    end
    return result
end

"""
    _multisets(v::Vector{T}, k::Int) where T

Generate all k-multisets (combinations with replacement) from vector v.

Returns vector of vectors, each containing k elements from v (may repeat).
"""
function _multisets(v::Vector{T}, k::Int) where {T}
    k == 0 && return [T[]]
    isempty(v) && return Vector{T}[]

    result = Vector{T}[]
    for (i, elem) in enumerate(v)
        for suffix in _multisets(v[i:end], k - 1)
            push!(result, vcat([elem], suffix))
        end
    end
    return result
end

"""
    get_ncbasis_deg(registry::VariableRegistry{FermionicAlgebra,T}, d::Int) where {T<:Signed}

Generate normal-form fermionic monomials of exactly degree d.

Fermionic operators are nilpotent (a²=0, (a†)²=0), so each mode can appear
at most once as creation AND at most once as annihilation. Normal-form words
are enumerated as pairs (C, A) where:
  - C ⊆ modes: creation operators (appear as -c for c ∈ C)
  - A ⊆ modes: annihilation operators (appear as +a for a ∈ A)
  - |C| + |A| = d

# Returns
`Vector{Monomial{FermionicAlgebra,T}}` containing only normal-form monomials.
"""
function get_ncbasis_deg(registry::VariableRegistry{FermionicAlgebra,T}, d::Int) where {T<:Signed}
    d < 0 && return Monomial{FermionicAlgebra,T}[]

    modes = _get_modes(registry)
    n_modes = length(modes)

    # Degree 0: identity monomial
    if d == 0
        return [Monomial{FermionicAlgebra}(T[])]
    end

    result = Monomial{FermionicAlgebra,T}[]

    # Enumerate all (n_create, n_annihilate) splits where n_create + n_annihilate = d
    for n_create in 0:min(d, n_modes)
        n_annihilate = d - n_create
        n_annihilate > n_modes && continue

        # Enumerate all n_create-subsets of modes for creation
        for create_modes in _combinations(modes, n_create)
            # Enumerate all n_annihilate-subsets of modes for annihilation
            for ann_modes in _combinations(modes, n_annihilate)
                # Build normal-form word: sorted creators (negative) then sorted annihilators (positive)
                word = T[]
                # Creators: negative indices, sorted by mode
                append!(word, T[-m for m in create_modes])
                # Annihilators: positive indices, sorted by mode
                append!(word, T[m for m in ann_modes])
                push!(result, Monomial{FermionicAlgebra}(word))
            end
        end
    end

    return result
end

"""
    get_ncbasis_deg(registry::VariableRegistry{BosonicAlgebra,T}, d::Int) where {T<:Signed}

Generate normal-form bosonic monomials of exactly degree d.

Bosonic operators are NOT nilpotent (b² ≠ 0), so each mode can appear
multiple times. Normal-form words are enumerated as pairs (C, A) where:
  - C: multiset of creation operators (modes with repetition)
  - A: multiset of annihilation operators (modes with repetition)
  - |C| + |A| = d

# Returns
`Vector{Monomial{BosonicAlgebra,T}}` containing only normal-form monomials.
"""
function get_ncbasis_deg(registry::VariableRegistry{BosonicAlgebra,T}, d::Int) where {T<:Signed}
    d < 0 && return Monomial{BosonicAlgebra,T}[]

    modes = _get_modes(registry)

    # Degree 0: identity monomial
    if d == 0
        return [Monomial{BosonicAlgebra}(T[])]
    end

    result = Monomial{BosonicAlgebra,T}[]

    # Enumerate all (n_create, n_annihilate) splits where n_create + n_annihilate = d
    for n_create in 0:d
        n_annihilate = d - n_create

        # Enumerate all multisets of size n_create from modes (with replacement)
        for create_multiset in _multisets(modes, n_create)
            # Enumerate all multisets of size n_annihilate from modes
            for ann_multiset in _multisets(modes, n_annihilate)
                # Build normal-form word: sorted creators (negative) then sorted annihilators (positive)
                word = T[]
                # Creators: negative indices, already sorted
                append!(word, T[-m for m in create_multiset])
                # Annihilators: positive indices, already sorted
                append!(word, T[m for m in ann_multiset])
                push!(result, Monomial{BosonicAlgebra}(word))
            end
        end
    end

    return result
end

# Specialized get_ncbasis for fermionic/bosonic that returns Vector{Monomial}
function get_ncbasis(registry::VariableRegistry{FermionicAlgebra,T}, d::Int) where {T<:Signed}
    result = Monomial{FermionicAlgebra,T}[]
    for deg in 0:d
        append!(result, get_ncbasis_deg(registry, deg))
    end
    return result
end

function get_ncbasis(registry::VariableRegistry{BosonicAlgebra,T}, d::Int) where {T<:Signed}
    result = Monomial{BosonicAlgebra,T}[]
    for deg in 0:d
        append!(result, get_ncbasis_deg(registry, deg))
    end
    return result
end
