"""
    Pauli Algebra Simplification

Implements simplification for Pauli spin operators satisfying:
- σᵢ² = I (involution: each Pauli squares to identity)
- Different sites commute: [σᵢⱼ, σₖₗ] = 0 for j ≠ l
- Same site cyclic products: σₓσᵧ = iσz, σᵧσz = iσₓ, σzσₓ = iσᵧ

# Variable Encoding Convention
Variables are ordered by site first, then by type (x, y, z):
- Index 1: σx₁, Index 2: σy₁, Index 3: σz₁
- Index 4: σx₂, Index 5: σy₂, Index 6: σz₂

For index `idx`:
- site = (idx - 1) ÷ 3 + 1
- pauli_type = (idx - 1) % 3 (0=X, 1=Y, 2=Z)

# Algorithm
1. Sort by site using stable sort (sites commute, stable sort gives deterministic ordering)
2. Linear pass: reduce each site group to at most one Pauli operator
3. Track complex phase coefficient through reductions

# Examples
```jldoctest
julia> using NCTSSoS

julia> word = [1, 2];  # σx₁ σy₁ = iσz₁ (raw word; not in Pauli normal form)

julia> m = simplify(PauliAlgebra, word);  # canonical expansion as `Monomial`

julia> p = Polynomial(m);  # convert internal phase encoding to numeric coefficient

julia> coefficients(p)[1]
0.0 + 1.0im

julia> monomials(p)[1].word
1-element Vector{Int64}:
 3
```
"""

# Encoding helper functions (exported for testing)

"""
    _pauli_site(idx::Integer) -> Int

Extract site number from Pauli variable index.
Site = (idx - 1) ÷ 3 + 1
"""
@inline _pauli_site(idx::Integer) = (idx - 1) ÷ 3 + 1

"""
    _pauli_type(idx::Integer) -> Int

Extract Pauli type from variable index.
Type = (idx - 1) % 3, where 0=X, 1=Y, 2=Z
"""
@inline _pauli_type(idx::Integer) = (idx - 1) % 3

"""
    _pauli_index(site::Integer, type::Integer) -> Int

Create variable index from site and Pauli type.
Index = (site - 1) * 3 + type + 1
"""
@inline _pauli_index(site::Integer, type::Integer) = (site - 1) * 3 + type + 1

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_pauli_word!(word::Vector{T}) where {T<:Integer}

Check that a Pauli word is in canonical form. Throws `ArgumentError` if invalid.

Canonical form requirements:
- ≤1 operator per site (no σ² terms)
- Sites sorted in ascending order

This is used by `NormalMonomial{PauliAlgebra,T}` constructor to enforce invariants.
"""
function _validate_pauli_word!(word::Vector{T}) where {T<:Integer}
    length(word) <= 1 && return nothing

    prev_site = _pauli_site(word[1])
    for i in 2:length(word)
        curr_site = _pauli_site(word[i])
        if curr_site < prev_site
            throw(ArgumentError(
                "Pauli word not sorted by site: site $curr_site at index $i " *
                "comes after site $prev_site"
            ))
        elseif curr_site == prev_site
            throw(ArgumentError(
                "Pauli word has multiple operators on site $curr_site " *
                "(indices $(i-1) and $i). Use simplify(PauliAlgebra, word) for raw words."
            ))
        end
        prev_site = curr_site
    end
    return nothing
end

# Connect validation hook used by `NormalMonomial{A,T}` inner constructor.
_validate_word!(::Type{PauliAlgebra}, word::Vector{T}) where {T<:Integer} =
    _validate_pauli_word!(word)

# =============================================================================
# Cyclic product table
# =============================================================================
# Cyclic product table:
# σₐσᵦ = i * ε(a,b) * σ_c where c = (a+b) mod 3 if a≠b (sort of)
# More precisely:
#   XY = iZ, YZ = iX, ZX = iY (cyclic: phase = +i)
#   YX = -iZ, ZY = -iX, XZ = -iY (anti-cyclic: phase = -i)

"""
    _pauli_product(type1::Int, type2::Int) -> Tuple{UInt8, Int}

Compute the product of two Pauli operators on the SAME site.
Returns (phase_k_delta, result_type) where:
- phase_k_delta ∈ 0:3 encodes phase contribution as (im)^phase_k_delta
- result_type is 0, 1, or 2 (or -1 for identity)

If type1 == type2, returns (0, -1) to indicate identity (σᵢ² = I).

Phase encoding: 0 → 1, 1 → i, 2 → -1, 3 → -i
"""
@inline function _pauli_product(type1::Int, type2::Int)
    # Same type: σᵢ² = I (phase = 1 = (im)^0)
    type1 == type2 && return (UInt8(0), -1)

    # Different types: cyclic or anti-cyclic product
    # XY→Z, YZ→X, ZX→Y (cyclic, +i = (im)^1)
    # YX→Z, ZY→X, XZ→Y (anti-cyclic, -i = (im)^3)

    # Result type: the one that's neither type1 nor type2
    # For {0,1,2}, the third one is 3 - type1 - type2
    result_type = 3 - type1 - type2

    # Determine phase_k: 1 for cyclic (+i), 3 for anti-cyclic (-i)
    # Cyclic order: 0→1→2→0 (X→Y→Z→X)
    # (type2 - type1 + 3) % 3 == 1 means cyclic
    phase_k = if (type2 - type1 + 3) % 3 == 1
        UInt8(1)   # +i = (im)^1
    else
        UInt8(3)   # -i = (im)^3
    end

    return (phase_k, result_type)
end

"""
    _simplify_pauli_word!(word::Vector{T}) where {T<:Integer} -> Tuple{Vector{T}, UInt8}

In-place site-aware simplification for Pauli algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, Pauli product rules apply: σₓσᵧ = iσz, σᵢ² = I, etc.

Returns a tuple of (simplified_word, phase_k) where:
- simplified_word is a new vector with the reduced operators
- phase_k ∈ 0:3 encodes the accumulated phase as (im)^phase_k

Phase encoding: 0 → 1, 1 → i, 2 → -1, 3 → -i

# Algorithm
1. Stable sort by site (using `_pauli_site`)
2. Reduce each site group to at most one operator, tracking phase_k
"""
function _simplify_pauli_word!(word::Vector{T}) where {T<:Integer}
    phase_k = UInt8(0)  # Start with phase = 1 = (im)^0

    # Empty or single: nothing to simplify
    length(word) <= 1 && return (word, phase_k)

    # Stage 1: Sort by site (stable sort preserves within-site order for determinism)
    sort!(word, alg=InsertionSort, by=_pauli_site)

    # Stage 2: Reduce same-site groups in a single linear pass
    # After sorting, all operators on the same site are adjacent
    result = T[]
    i = 1
    while i <= length(word)
        current_site = _pauli_site(word[i])

        # Find extent of this site's group: word[i:j-1]
        j = i + 1
        while j <= length(word) && _pauli_site(word[j]) == current_site
            j += 1
        end

        # Reduce all operators in word[i:j-1] to at most one operator
        # current_type: -1 = identity accumulated, 0/1/2 = X/Y/Z accumulated
        current_type = _pauli_type(word[i])
        for k in (i + 1):(j - 1)
            next_type = _pauli_type(word[k])
            if current_type == -1
                # Identity * next = next (no phase change)
                current_type = next_type
            else
                delta_k, result_type = _pauli_product(current_type, next_type)
                phase_k = (phase_k + delta_k) % UInt8(4)  # Accumulate phase mod 4
                current_type = result_type  # Could be -1 (identity)
            end
        end

        # Append result for this site (skip if identity)
        if current_type != -1
            push!(result, T(_pauli_index(current_site, current_type)))
        end

        i = j  # Move to next site group
    end

    return (result, phase_k)
end

"""
    _phase_k_to_complex(phase_k::UInt8) -> ComplexF64

Convert phase_k encoding to complex number.
phase_k ∈ 0:3 maps to (im)^phase_k: 0→1, 1→i, 2→-1, 3→-i
"""
@inline function _phase_k_to_complex(phase_k::UInt8)
    phase_k == 0 && return ComplexF64(1.0, 0.0)
    phase_k == 1 && return ComplexF64(0.0, 1.0)
    phase_k == 2 && return ComplexF64(-1.0, 0.0)
    return ComplexF64(0.0, -1.0)  # phase_k == 3
end

"""
    simplify(m::NormalMonomial{PauliAlgebra,T}) where T -> Monomial

Simplify a Pauli algebra monomial.

Returns a Term containing the simplified monomial and accumulated phase coefficient.
The original monomial is unchanged.

# Algebraic Rules
- σᵢ² = I (involution)
- Different sites commute
- Same site cyclic products: σₓσᵧ = iσz, σᵧσz = iσₓ, σzσₓ = iσᵧ

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1, 2]);  # σx₁ σy₁

julia> t = simplify(m);

julia> t.coefficient
0.0 + 1.0im

julia> t.monomial.word
1-element Vector{Int64}:
 3

julia> m.word  # Original unchanged
2-element Vector{Int64}:
 1
 2
```
"""
function simplify(m::NormalMonomial{PauliAlgebra,T}) where {T}
    word_copy = copy(m.word)
    result, phase_k = _simplify_pauli_word!(word_copy)
    mono = NormalMonomial{PauliAlgebra,T}(result, _OWNED_NORMAL_MONOMIAL)
    return Monomial(phase_k, mono)
end

# =============================================================================
# Specialized Outer Constructor (validates, rejects non-canonical)
# =============================================================================

"""
    NormalMonomial{PauliAlgebra}(word::Vector{T}) where {T<:Integer}

Construct a Pauli monomial, validating that the input is in canonical form.

Throws `ArgumentError` if the word is not canonical. For non-canonical words,
use `simplify(PauliAlgebra, word)` which auto-canonicalizes and returns a `Monomial`.

Canonical form requirements:
- ≤1 operator per site (no σ² terms)
- Sites sorted in ascending order

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1]);  # σx₁ - canonical

julia> m.word
1-element Vector{Int64}:
 1

julia> NormalMonomial{PauliAlgebra}([1, 2])  # σx₁ σy₁ - NOT canonical (same site)
ERROR: ArgumentError: Pauli word has multiple operators on site 1 (indices 1 and 2). Use simplify(PauliAlgebra, word) for raw words.
```
"""
function NormalMonomial{PauliAlgebra}(word::Vector{T}) where {T<:Integer}
    word_filtered = filter(!iszero, word)
    return NormalMonomial{PauliAlgebra,T}(word_filtered, _OWNED_NORMAL_MONOMIAL)
end

function Base.:*(m1::NormalMonomial{PauliAlgebra,T}, m2::NormalMonomial{PauliAlgebra,T}) where {T<:Integer}
    return simplify(PauliAlgebra, vcat(m1.word, m2.word))
end
