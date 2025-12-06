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
    for i in 1:length(word)-1
        word[i] == word[i+1] && return true
    end
    false
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
    result
end

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
function get_ncbasis_deg(::Type{A}, n::Int, d::Int;
    T::Type{<:Integer}=Int,
    filter_constraint::Bool=false) where {A<:AlgebraType}
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
    monos
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
function get_ncbasis(::Type{A}, n::Int, d::Int;
    T::Type{<:Integer}=Int,
    filter_constraint::Bool=false) where {A<:AlgebraType}
    basis = Monomial{A,T}[]
    for deg in 0:d
        append!(basis, get_ncbasis_deg(A, n, deg; T=T, filter_constraint=filter_constraint))
    end
    basis
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
