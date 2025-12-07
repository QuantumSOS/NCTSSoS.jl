# =============================================================================
# Basis Generation for Non-Commutative Monomials
# =============================================================================
#
# Generates bases of monomials for non-commutative polynomial optimization.
# Reference: NCTSSOS get_ncbasis in utils.jl
#
# Design:
# - Returns Vector{Monomial{A,T}} instead of Vector{Vector{UInt16}}
# - Algebra type A enables type-safe dispatch
# - Integer type T is configurable (default: Int)
# - Optional constraint filtering for algebras with simplification rules
#
# TODO: Rewrite basis generation to use VariableRegistry for type consistency.
#       Currently defaults to T=Int, but @ncpolyvar uses UInt64 via VariableRegistry.
#       This causes type mismatches when comparing/sorting monomials from different sources.
#       The fix is to make basis functions registry-aware so they produce compatible types.
#       Related: removed cross-type isless(::Monomial{A,T1}, ::Monomial{A,T2}) as it
#       used abs() which incorrectly treated creation/annihilation operators as equal.

"""
    has_consecutive_repeats(word::Vector) -> Bool

Check if a word vector has consecutive repeated elements.
Used for filtering monomials in algebras where consecutive repeats simplify
(e.g., UnipotentAlgebra where X_i * X_i = I).

# Examples
```jldoctest
julia> has_consecutive_repeats([1, 2, 3])
false

julia> has_consecutive_repeats([1, 1, 2])
true

julia> has_consecutive_repeats([1, 2, 2])
true
```
"""
function has_consecutive_repeats(word::Vector)
    for i in 1:(length(word) - 1)
        word[i] == word[i + 1] && return true
    end
    return false
end

"""
    _generate_words(n::Int, d::Int, ::Type{T}) where T<:Integer -> Vector{Vector{T}}

Internal helper to generate all words of exactly length d using indices 1:n.

# Arguments
- `n`: Number of variables (uses indices 1:n)
- `d`: Exact word length (degree)
- `T`: Integer type for the word elements

# Returns
Vector of all possible words of length d.
"""
function _generate_words(n::Int, d::Int, ::Type{T}) where {T<:Integer}
    d == 0 && return [T[]]
    n == 0 && return Vector{T}[]
    d == 1 && return [T[i] for i in 1:n]

    result = Vector{T}[]
    for i in 1:n
        for suffix in _generate_words(n, d - 1, T)
            push!(result, vcat(T[i], suffix))
        end
    end
    return result
end

# TODO: make this take variable registry this should fix the output vector type
"""
    get_ncbasis_deg(::Type{A}, n::Int, d::Int; T::Type{<:Integer}=Int,
                    filter_constraint::Bool=false) where A<:AlgebraType

Generate all monomials of algebra type A with n variables, of exactly degree d.

# Arguments
- `A`: Algebra type (PauliAlgebra, FermionicAlgebra, etc.)
- `n`: Number of variables (uses indices 1:n)
- `d`: Exact degree
- `T`: Integer type for the word (default: Int)
- `filter_constraint`: If true, filter out monomials that would simplify
  (e.g., consecutive repeats for UnipotentAlgebra). Default: false.

# Returns
- `Vector{Monomial{A,T}}`: All monomials of degree d, sorted lexicographically

# Examples
```jldoctest
julia> basis = get_ncbasis_deg(PauliAlgebra, 2, 2);

julia> length(basis)  # 2^2 = 4
4

julia> [m.word for m in basis]
4-element Vector{Vector{Int64}}:
 [1, 1]
 [1, 2]
 [2, 1]
 [2, 2]
```

Degree 0 returns only the identity:
```jldoctest
julia> basis = get_ncbasis_deg(PauliAlgebra, 3, 0);

julia> length(basis)
1

julia> isone(basis[1])
true
```
"""
function get_ncbasis_deg(
    ::Type{A}, n::Int, d::Int; T::Type{<:Integer}=Int, filter_constraint::Bool=false
) where {A<:AlgebraType}
    d == 0 && return [Monomial{A}(T[])]  # identity
    d < 0 && return Monomial{A,T}[]

    # Generate all words of length d using indices 1:n
    words = _generate_words(n, d, T)

    # Apply constraint filtering if requested
    if filter_constraint
        filter!(w -> !has_consecutive_repeats(w), words)
    end

    monos = [Monomial{A}(w) for w in words]
    sort!(monos)  # Uses isless for Monomial (degree-first, then lexicographic)
    return monos
end

"""
    get_ncbasis(::Type{A}, n::Int, d::Int; T::Type{<:Integer}=Int,
                filter_constraint::Bool=false) where A<:AlgebraType

Generate all monomials of algebra type A with n variables, up to and including degree d.
Returns a vector of monomials sorted by degree, then lexicographically.

# Arguments
- `A`: Algebra type (PauliAlgebra, FermionicAlgebra, etc.)
- `n`: Number of variables (uses indices 1:n)
- `d`: Maximum degree (inclusive)
- `T`: Integer type for the word (default: Int)
- `filter_constraint`: If true, filter out monomials that would simplify.
  Default: false.

# Returns
- `Vector{Monomial{A,T}}`: All monomials from degree 0 to d

# Count Formula
For n variables and maximum degree d, the total count is:
- Without filtering: 1 + n + n^2 + ... + n^d = (n^(d+1) - 1) / (n - 1)
- With filtering: depends on algebra constraints

# Examples
```jldoctest
julia> basis = get_ncbasis(PauliAlgebra, 2, 2);

julia> length(basis)  # 1 + 2 + 4 = 7
7
```

With constraint filtering for UnipotentAlgebra (removes consecutive repeats):
```jldoctest
julia> basis_unfiltered = get_ncbasis(UnipotentAlgebra, 2, 2);

julia> length(basis_unfiltered)  # 1 + 2 + 4 = 7
7

julia> basis_filtered = get_ncbasis(UnipotentAlgebra, 2, 2; filter_constraint=true);

julia> length(basis_filtered)  # 1 + 2 + 2 = 5 (removed [1,1] and [2,2])
5
```
"""
# TODO: need to generate basis more efficiently, see `NCTSSOS`
function get_ncbasis(
    ::Type{A}, n::Int, d::Int; T::Type{<:Integer}=Int, filter_constraint::Bool=false
) where {A<:AlgebraType}
    basis = Monomial{A,T}[]
    for deg in 0:d
        append!(basis, get_ncbasis_deg(A, n, deg; T=T, filter_constraint=filter_constraint))
    end
    return basis
end

"""
    monomials(::Type{A}, n::Int, d::Int; T::Type{<:Integer}=Int) where A<:AlgebraType

Generate all monomials of exactly degree d with n variables.
Alias for `get_ncbasis_deg` following MultivariatePolynomials.jl convention.

Note: This function is NOT exported to avoid conflict with the `monomials(::Polynomial)`
accessor which extracts monomials from a polynomial. Use `get_ncbasis_deg` explicitly
for basis generation.

# Arguments
- `A`: Algebra type
- `n`: Number of variables (uses indices 1:n)
- `d`: Exact degree
- `T`: Integer type for the word (default: Int)

# Returns
- `Vector{Monomial{A,T}}`: All monomials of degree d, sorted lexicographically

# Examples
```jldoctest
julia> using FastPolynomials

julia> basis = FastPolynomials.monomials(PauliAlgebra, 3, 2);

julia> length(basis)  # 3^2 = 9
9
```

See also: [`get_ncbasis_deg`](@ref), [`get_ncbasis`](@ref)
"""
monomials(::Type{A}, n::Int, d::Int; T::Type{<:Integer}=Int) where {A<:AlgebraType} =
    get_ncbasis_deg(A, n, d; T=T)

# =============================================================================
# Registry-Aware Basis Generation
# =============================================================================
#
# These functions generate bases using a VariableRegistry, producing properly
# typed monomials and returning Vector{Term} to handle simplification results.

"""
    _generate_all_words(indices::Vector{T}, d::Int) where {T<:Integer} -> Vector{Vector{T}}

Internal helper to generate all words of exactly length d using the provided indices.

Unlike `_generate_words(n, d, T)` which uses indices 1:n, this function accepts
an arbitrary vector of indices. This enables registry-aware basis generation
where indices may be non-contiguous (e.g., site-encoded unsigned indices).

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

This is the registry-aware version of `get_ncbasis_deg(::Type{A}, n, d)`.
It returns `Vector{Term}` to properly handle simplification results, since
some algebras (like Pauli) may produce coefficients other than 1.0 during
simplification.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Exact degree

# Returns
- `Vector{Term{Monomial{A,T}, C}}`: All simplified terms of degree d

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

This is the registry-aware version of `get_ncbasis(::Type{A}, n, d)`.
It returns `Vector{Term}` to properly handle simplification results.

# Arguments
- `registry`: Variable registry containing the indices to use
- `d`: Maximum degree (inclusive)

# Returns
- `Vector{Term{Monomial{A,T}, C}}`: All simplified terms from degree 0 to d

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
