# =============================================================================
# Canonicalization Algorithms
# =============================================================================
#
# This module provides functions to compute canonical forms of monomials and
# polynomials. Canonical forms are essential for:
# - Identifying equivalent monomials under symmetry transformations
# - Trace optimization (cyclic equivalence)
# - Polynomial simplification (combining like terms)
#
# Reference: NCTSSOS _sym_canon and _cyclic_canon in utils.jl
#
# IMPORTANT: For algebras with site-based commutation (ProjectorAlgebra,
# UnipotentAlgebra, NonCommutativeAlgebra with unsigned indices), the
# symmetric canonicalization must sort BOTH the word AND its reverse by site
# before comparing. This matches NCTSSOS's behavior where `_comm` is applied
# to both word and reverse(word) before taking the minimum.
# =============================================================================

# =============================================================================
# Symmetric Canonicalization
# =============================================================================

"""
    symmetric_canon(word::Vector{T}) where T

Return the lexicographically smaller of `word` and `reverse(word)`.
Used for identifying symmetric equivalence classes where a word and its
reversal represent the same physical quantity (e.g., in eigenvalue optimization).

The algorithm compares elements from both ends simultaneously, returning
as soon as a difference is found. This is O(n/2) in the best case.

!!! warning "Site-based algebras"
    For algebras with site-based commutation (using unsigned integer indices
    with encoded sites), use the `Monomial` dispatch instead. The raw word
    version does not handle site-based sorting of the reversed word.

# Arguments
- `word::Vector{T}`: A word (sequence of indices) representing a monomial

# Returns
- A copy of either `word` or `reverse(word)`, whichever is lexicographically smaller

# Examples
```jldoctest
julia> symmetric_canon([1, 2, 3])
3-element Vector{Int64}:
 1
 2
 3

julia> symmetric_canon([3, 2, 1])
3-element Vector{Int64}:
 1
 2
 3

julia> symmetric_canon([1, 2, 1])  # palindrome
3-element Vector{Int64}:
 1
 2
 1
```

See also: [`cyclic_canon`](@ref), [`cyclic_symmetric_canon`](@ref), [`canonicalize`](@ref)
"""
function symmetric_canon(word::Vector{T}) where {T}
    n = length(word)

    # Early return for empty or single-element words
    n <= 1 && return word

    # Compare from both ends, decide as soon as a difference is found
    for i in 1:div(n, 2)
        if word[i] < word[n+1-i]
            return word
        elseif word[i] > word[n+1-i]
            return reverse(word)
        end
        # If equal, continue to next pair
    end

    # All compared pairs are equal (palindrome or equal under reversal)
    return word
end

# =============================================================================
# Site-Awareness Trait System
# =============================================================================
#
# Uses Julian two-level dispatch pattern:
# 1. _site_trait(A, T) determines IF site-aware canonicalization applies
# 2. _site_extractor(A) determines HOW to compute site for sorting
# =============================================================================

"""Site-aware canonicalization trait - operators on different sites commute."""
struct SiteAware end

"""No site-aware canonicalization - use raw word algorithms."""
struct NoSiteAware end

# Trait dispatch: which algebra+index combinations support site-based commutation
@inline _site_trait(::Type{PauliAlgebra}, ::Type{<:Integer}) = SiteAware()
@inline _site_trait(::Type{NonCommutativeAlgebra}, ::Type{<:Unsigned}) = SiteAware()
@inline _site_trait(::Type{ProjectorAlgebra}, ::Type{<:Unsigned}) = SiteAware()
@inline _site_trait(::Type{UnipotentAlgebra}, ::Type{<:Unsigned}) = SiteAware()
@inline _site_trait(::Type{<:AlgebraType}, ::Type) = NoSiteAware()

# Site extractor: returns function to extract site from index
@inline _site_extractor(::Type{PauliAlgebra}) = idx -> (Int(idx) - 1) ÷ 3 + 1
@inline _site_extractor(::Type{<:AlgebraType}) = decode_site

# =============================================================================
# Site-Sorted Helpers
# =============================================================================

"""In-place sort word by site using stable InsertionSort."""
@inline function _sort_by_site!(word::Vector{T}, ::Type{A}) where {T,A<:AlgebraType}
    return sort!(word; by=_site_extractor(A), alg=InsertionSort)
end

"""Symmetric canon with site sorting: min(sorted_word, sorted_reverse)."""
function _symmetric_canon_site_sorted(word::Vector{T}, ::Type{A}) where {T,A<:AlgebraType}
    sorted_word = _sort_by_site!(copy(word), A)
    sorted_rev = _sort_by_site!(reverse(word), A)  # reverse already copies
    return min(sorted_word, sorted_rev)
end

"""
Cyclic canon with site sorting: for each rotation, sort by site, return smallest.

Algorithm: O(n) rotations × O(n log n) sort = O(n² log n)
"""
function _cyclic_canon_site_sorted(word::Vector{T}, ::Type{A}) where {T,A<:AlgebraType}
    n = length(word)
    n <= 1 && return _sort_by_site!(copy(word), A)

    best = _sort_by_site!(copy(word), A)
    for offset in 1:(n-1)
        rotation = [word[mod1(i + offset, n)] for i in 1:n]
        _sort_by_site!(rotation, A)
        if rotation < best
            best = rotation
        end
    end
    return best
end

"""Cyclic-symmetric canon with site sorting."""
function _cyclic_symmetric_canon_site_sorted(word::Vector{T}, ::Type{A}) where {T,A<:AlgebraType}
    return min(_cyclic_canon_site_sorted(word, A), _cyclic_canon_site_sorted(reverse(word), A))
end

"""
    symmetric_canon(m::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}

Return a new monomial with the symmetrically canonical word.
Preserves the algebra type.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([3, 2, 1]);

julia> symmetric_canon(m).word
3-element Vector{Int64}:
 1
 2
 3
```
"""
# Two-level dispatch: entry point dispatches on trait
symmetric_canon(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} =
    _symmetric_canon(_site_trait(A, T), m)

# NoSiteAware: use raw word algorithm
_symmetric_canon(::NoSiteAware, m::Monomial{A}) where {A} =
    Monomial{A}(symmetric_canon(m.word))

# SiteAware: use site-sorted algorithm
_symmetric_canon(::SiteAware, m::Monomial{A}) where {A} =
    Monomial{A}(_symmetric_canon_site_sorted(m.word, A))

# =============================================================================
# Cyclic Canonicalization
# =============================================================================

"""
    cyclic_canon(word::Vector{T}) where T

Return the lexicographically smallest cyclic rotation of `word`.
Used for trace optimization where cyclic permutations are equivalent
(since trace is invariant under cyclic permutation: tr(ABC) = tr(BCA) = tr(CAB)).

# Arguments
- `word::Vector{T}`: A word (sequence of indices) representing a monomial

# Returns
- A copy of the lexicographically smallest cyclic rotation

# Algorithm
Finds the optimal rotation offset first (O(n²) comparisons, O(1) allocations),
then constructs the result vector once. Time complexity: O(n²) comparisons.

# Examples
```jldoctest
julia> cyclic_canon([2, 3, 1])
3-element Vector{Int64}:
 1
 2
 3

julia> cyclic_canon([3, 1, 2])
3-element Vector{Int64}:
 1
 2
 3

julia> cyclic_canon([1, 2, 3])  # already canonical
3-element Vector{Int64}:
 1
 2
 3
```

See also: [`symmetric_canon`](@ref), [`cyclic_symmetric_canon`](@ref), [`canonicalize`](@ref)
"""
function cyclic_canon(word::Vector{T}) where {T}
    n = length(word)
    n <= 1 && return copy(word)

    best_offset = 0

    # Find the lexicographically smallest rotation by comparing offsets
    for offset in 1:(n-1)
        # Compare rotation at best_offset vs rotation at offset
        for i in 1:n
            idx_best = mod1(i + best_offset, n)
            idx_curr = mod1(i + offset, n)

            if word[idx_curr] < word[idx_best]
                best_offset = offset
                break
            elseif word[idx_curr] > word[idx_best]
                break
            end
            # If equal, continue comparing next position
        end
    end

    # Allocate result only once
    return [word[mod1(i + best_offset, n)] for i in 1:n]
end

"""
    cyclic_canon(m::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}

Return a new monomial with the cyclically canonical word.
Preserves the algebra type.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([2, 3, 1]);

julia> cyclic_canon(m).word
3-element Vector{Int64}:
 1
 2
 3
```
"""
# Two-level dispatch: entry point dispatches on trait
cyclic_canon(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} =
    _cyclic_canon(_site_trait(A, T), m)

# NoSiteAware: use raw word algorithm
_cyclic_canon(::NoSiteAware, m::Monomial{A}) where {A} =
    Monomial{A}(cyclic_canon(m.word))

# SiteAware: use site-sorted algorithm
_cyclic_canon(::SiteAware, m::Monomial{A}) where {A} =
    Monomial{A}(_cyclic_canon_site_sorted(m.word, A))

# =============================================================================
# Combined Cyclic-Symmetric Canonicalization
# =============================================================================

"""
    cyclic_symmetric_canon(word::Vector{T}) where T

Return the lexicographically smallest among all cyclic rotations of `word`
and all cyclic rotations of `reverse(word)`.

Used for trace optimization with both cyclic invariance (trace property)
and symmetric invariance (adjoint/transpose symmetry).

# Arguments
- `word::Vector{T}`: A word (sequence of indices) representing a monomial

# Returns
- The lexicographically smallest word among all cyclic rotations of both
  the original word and its reversal

# Examples
```jldoctest
julia> cyclic_symmetric_canon([3, 2, 1])
3-element Vector{Int64}:
 1
 2
 3

julia> cyclic_symmetric_canon([3, 1, 4, 2])
4-element Vector{Int64}:
 1
 3
 2
 4
```

See also: [`symmetric_canon`](@ref), [`cyclic_canon`](@ref), [`canonicalize`](@ref)
"""
function cyclic_symmetric_canon(word::Vector{T}) where {T}
    return min(cyclic_canon(word), cyclic_canon(reverse(word)))
end

"""
    cyclic_symmetric_canon(m::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}

Return a new monomial with the cyclic-symmetrically canonical word.
Preserves the algebra type.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([3, 2, 1]);

julia> cyclic_symmetric_canon(m).word
3-element Vector{Int64}:
 1
 2
 3
```
"""
# Two-level dispatch: entry point dispatches on trait
cyclic_symmetric_canon(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} =
    _cyclic_symmetric_canon(_site_trait(A, T), m)

# NoSiteAware: use raw word algorithm
_cyclic_symmetric_canon(::NoSiteAware, m::Monomial{A}) where {A} =
    Monomial{A}(cyclic_symmetric_canon(m.word))

# SiteAware: use site-sorted algorithm
_cyclic_symmetric_canon(::SiteAware, m::Monomial{A}) where {A} =
    Monomial{A}(_cyclic_symmetric_canon_site_sorted(m.word, A))

# =============================================================================
# Unified Canonicalize Interface
# =============================================================================

"""
    canonicalize(m::Monomial{A,T}; cyclic::Bool=false) where {A<:AlgebraType, T<:Integer}

Canonicalize a monomial using the specified method.

This function operates on the structural word level and assumes that algebraic
simplification (e.g., `simplify`) has already been applied if needed.

# Arguments
- `m::Monomial{A,T}`: The monomial to canonicalize
- `cyclic::Bool=false`: Canonicalization mode
  - `false` (default): Symmetric canonicalization only (min of word and reverse)
  - `true`: Cyclic-symmetric canonicalization (min of all rotations and their reversals)

# Returns
- A new monomial in canonical form

# Use Cases
- `cyclic=false`: Eigenvalue optimization (word and reverse are equivalent under adjoint)
- `cyclic=true`: Trace optimization (cyclic rotations are equivalent under trace)

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([3, 2, 1]);

julia> canonicalize(m).word  # symmetric (default)
3-element Vector{Int64}:
 1
 2
 3

julia> canonicalize(m; cyclic=true).word  # cyclic + symmetric
3-element Vector{Int64}:
 1
 2
 3
```

See also: [`symmetric_canon`](@ref), [`cyclic_canon`](@ref), [`cyclic_symmetric_canon`](@ref)
"""
function canonicalize(m::Monomial{A,T}; cyclic::Bool=false) where {A<:AlgebraType,T<:Integer}
    if cyclic
        return cyclic_symmetric_canon(m)
    else
        return symmetric_canon(m)
    end
end

"""
    canonicalize(p::Polynomial{A,T,C}; cyclic::Bool=false) where {A<:AlgebraType, T<:Integer, C<:Number}

Canonicalize each monomial in the polynomial and combine like terms.

This function:
1. Applies canonicalization to each monomial in the polynomial
2. Re-processes terms to combine those with equivalent canonical monomials
3. Removes any zero-coefficient terms resulting from cancellation

# Arguments
- `p::Polynomial{A,T,C}`: The polynomial to canonicalize
- `cyclic::Bool=false`: Canonicalization mode (see `canonicalize(::Monomial)`)

# Returns
- A new polynomial with canonical monomials and combined coefficients

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([3, 2, 1]);

julia> m2 = Monomial{PauliAlgebra}([1, 2, 3]);

julia> p = Polynomial([Term(1.0+0im, m1), Term(2.0+0im, m2)]);

julia> p_canon = canonicalize(p);

julia> length(terms(p_canon))  # Combined into one term
1

julia> coefficients(p_canon)[1]
3.0 + 0.0im
```

See also: [`canonicalize(::Monomial)`](@ref)
"""
function canonicalize(
    p::Polynomial{A,T,C}; cyclic::Bool=false
) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return zero(Polynomial{A,T,C})

    # Canonicalize each term's monomial
    new_terms = [Term(t.coefficient, canonicalize(t.monomial; cyclic=cyclic)) for t in p.terms]

    # Polynomial constructor handles sorting, deduplication, and zero removal
    Polynomial(new_terms)
end
