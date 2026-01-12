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

julia> m = simplify(PauliAlgebra, word);  # canonical expansion

julia> p = Polynomial(m);  # convert internal phase encoding to numeric coefficient

julia> coefficients(p)[1]
0.0 + 1.0im

julia> monomials(p)[1].word
1-element Vector{Int64}:
 3
```
"""

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
function _validate_word(::Type{PauliAlgebra}, word::Vector{T}) where {T<:Integer}
    length(word) <= 1 && return nothing

    _is_site_sorted_pauli(word) || throw(ArgumentError("Pauli word not sorted by site"))

    prev_site = _pauli_site(word[1])
    @inbounds for i in 2:length(word)
        curr_site = _pauli_site(word[i])
        curr_site == prev_site &&
            throw(ArgumentError(
                "Pauli word has multiple operators on site $curr_site " *
                "(indices $(i-1) and $i). Use simplify(PauliAlgebra, word) for raw words."
            ))
    end
    return nothing
end

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
    simplify!(::Type{PauliAlgebra}, word::Vector{T}) where {T<:Integer} -> Tuple{Vector{T}, UInt8}

In-place site-aware simplification for Pauli algebra word vectors.

Operators on different sites commute and are sorted by site (ascending).
Within each site, Pauli product rules apply: σₓσᵧ = iσz, σᵢ² = I, etc.

Returns a tuple of (word, phase_k) where:
- word is the input vector, mutated and resized in-place
- phase_k ∈ 0:3 encodes the accumulated phase as (im)^phase_k

Phase encoding: 0 → 1, 1 → i, 2 → -1, 3 → -i

# Algorithm
1. Stable sort by site (using `_pauli_site`)
2. Two-pointer reduction: read_idx scans groups, write_idx writes results
3. Resize to truncate trailing elements
"""
function simplify!(::Type{PauliAlgebra}, word::Vector{T}) where {T<:Integer}
    phase_k = UInt8(0)  # Start with phase = 1 = (im)^0

    # Empty or single: nothing to simplify
    length(word) <= 1 && return (word, phase_k)

    # Stage 1: Sort by site (stable sort preserves within-site order for determinism)
    sort!(word, alg=InsertionSort, by=_pauli_site)

    # Stage 2: Reduce same-site groups using two-pointer in-place modification
    # read_idx scans through sorted word, write_idx tracks output position
    write_idx = 0
    read_idx = 1
    while read_idx <= length(word)
        current_site = _pauli_site(word[read_idx])

        # Find extent of this site's group: word[read_idx:j-1]
        j = read_idx + 1
        while j <= length(word) && _pauli_site(word[j]) == current_site
            j += 1
        end

        # Reduce all operators in word[read_idx:j-1] to at most one operator
        # current_type: -1 = identity accumulated, 0/1/2 = X/Y/Z accumulated
        current_type = _pauli_type(word[read_idx])
        for k in (read_idx + 1):(j - 1)
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

        # Write result for this site in-place (skip if identity)
        if current_type != -1
            write_idx += 1
            word[write_idx] = T(_pauli_index(current_site, current_type))
        end

        read_idx = j  # Move to next site group
    end

    # Truncate trailing elements
    resize!(word, write_idx)

    return (word, phase_k)
end

"""
    simplify(m::NormalMonomial{PauliAlgebra,T}) where T

Simplify a Pauli algebra monomial.

Returns a `(word, phase)` tuple containing the simplified word and accumulated phase.
The original monomial is unchanged.

# Algebraic Rules
- σᵢ² = I (involution)
- Different sites commute
- Same site cyclic products: σₓσᵧ = iσz, σᵧσz = iσₓ, σzσₓ = iσᵧ
"""

"""
    simplify(::Type{PauliAlgebra}, word::Vector{T}) where {T<:Integer}

Simplify a raw Pauli word into canonical form.

This is the primary entry point for Pauli simplification. Takes a raw word vector
and returns a `(word, phase)` tuple.
"""
function simplify(::Type{PauliAlgebra}, word::Vector{T}) where {T<:Integer}
    return simplify!(PauliAlgebra, copy(word))
end
