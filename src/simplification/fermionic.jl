"""
    Fermionic Algebra Simplification using Generalized Wick's Theorem

Normal orders fermionic creation/annihilation operators using the systematic
Wick's theorem algorithm from the Generalized Time-Independent Wick Theorem.

# Algebraic Rules
- {aᵢ, aⱼ†} = δᵢⱼ  (anticommutation gives delta correction)
- {aᵢ, aⱼ} = 0      (annihilation-annihilation)
- {aᵢ†, aⱼ†} = 0    (creation-creation)

# Variable Encoding
- Annihilation aᵢ: positive index `i`
- Creation aᵢ†: negative index `-i`

# Normal Ordering
- All creators (negative) on LEFT
- All annihilators (positive) on RIGHT
- Within groups: creators sorted by mode descending, annihilators by mode ascending

# Wick's Theorem Algorithm (4 steps)
1. Identify valid contractions (pairs with non-zero contraction values)
2. Generate all non-overlapping combinations of contractions
3. Form N-products for each combination (normal-ordered remaining operators)
4. Evaluate signs using permutation parity
"""

"""
    has_even_parity(m::NormalMonomial{FermionicAlgebra,T}) where T -> Bool

Check if a fermionic monomial has even parity (even number of operators).

In fermionic systems, parity superselection forbids observables with odd fermion
number from having non-zero expectation values. Only operators with even parity
(even number of creation/annihilation operators) can have non-zero physical
expectation values.

# Arguments
- `m::NormalMonomial{FermionicAlgebra,T}`: A fermionic monomial

# Returns
- `true` if the monomial has an even number of operators (including 0 for identity)
- `false` if the monomial has an odd number of operators

# Examples
```jldoctest
julia> using NCTSSoS

julia> registry, (a, a_dag) = create_fermionic_variables(1:2);

julia> has_even_parity(one(typeof(a[1])))  # Identity has 0 operators
true

julia> has_even_parity(a[1])  # Single annihilation: 1 operator (odd)
false

julia> has_even_parity(a_dag[1] * a[1])  # Number operator: 2 operators (even)
true
```

See also: [`FermionicAlgebra`](@ref), [`simplify`](@ref)
"""
function has_even_parity(m::NormalMonomial{FermionicAlgebra,T}) where {T}
    return iseven(length(m.word))
end

"""
    Base.iszero(m::NormalMonomial{FermionicAlgebra,T}) -> Bool

Check if a fermionic monomial is zero due to nilpotency.

# Algorithm
A fermionic monomial is zero if for any mode, the net surplus of operators
(annihilation minus creation, or vice versa) is 2 or more. This is because:
- `aᵢ aᵢ = 0` (two annihilations with no creation to contract)
- `aᵢ† aᵢ† = 0` (two creations with no annihilation to contract)

The key insight is that anticommutation `{aᵢ, aⱼ} = 0` allows reordering,
but if there's a net surplus of ≥2 same-type operators for any mode,
they cannot all be contracted away and the result is zero.

# Examples
```jldoctest
julia> m = NormalMonomial{FermionicAlgebra}(Int32[1, 1]);  # a₁ a₁ = 0 (surplus 2)

julia> iszero(m)
true

julia> m2 = NormalMonomial{FermionicAlgebra}(Int32[1, 2, 1]);  # a₁ a₂ a₁ = -a₁ a₁ a₂ = 0

julia> iszero(m2)
true

julia> m3 = NormalMonomial{FermionicAlgebra}(Int32[1, 1, -1]);  # a₁ a₁ a₁† (surplus 1, not zero yet)

julia> iszero(m3)
false

julia> m4 = NormalMonomial{FermionicAlgebra}(Int32[1, -1, 1, -1]);  # a₁ a₁† a₁ a₁† ≠ 0

julia> iszero(m4)
false
```
"""
function Base.iszero(m::NormalMonomial{FermionicAlgebra,T}) where {T}
    word = m.word
    isempty(word) && return false

    # Find max mode to size the tracking array
    max_mode = maximum(_operator_mode, word)

    # Track net flux (annihilation - creation) for each mode
    # +1 for annihilation, -1 for creation
    flux = zeros(Int, max_mode)

    @inbounds for op in word
        mode = _operator_mode(op)
        if _is_creation(op)
            flux[mode] -= 1
        else
            flux[mode] += 1
        end
    end

    # Check if any mode has a surplus of 2 or more operators
    # |flux| >= 2 implies the operator is nilpotent
    @inbounds for f in flux
        abs(f) >= 2 && return true
    end

    return false
end

"""
    _is_nilpotent_by_flux(word::Vector{T}) where {T<:Signed} -> Bool

Check if a raw fermionic word is nilpotent using flux-based detection.
A fermionic word is zero if any mode has net surplus >= 2 (same-type operators).

This is the Vector{T} version of iszero(::NormalMonomial{FermionicAlgebra}) for use
before a NormalMonomial is constructed.
"""
function _is_nilpotent_by_flux(word::Vector{T}) where {T<:Signed}
    # Find max mode to size the tracking array
    max_mode = maximum(_operator_mode, word)

    # Track net flux (annihilation - creation) for each mode
    # +1 for annihilation, -1 for creation
    flux = zeros(Int, max_mode)

    @inbounds for op in word
        mode = _operator_mode(op)
        if _is_creation(op)
            flux[mode] -= 1
        else
            flux[mode] += 1
        end
    end

    # Check if any mode has a surplus of 2 or more operators
    # |flux| >= 2 implies the operator is nilpotent
    @inbounds for f in flux
        abs(f) >= 2 && return true
    end

    return false
end

"""
    _permutation_parity(perm::AbstractVector{<:Integer}) -> Int

Compute the parity (sign) of a permutation using the cycle counting algorithm.

# Algorithm
The sign of a permutation equals `(-1)^(n - c)` where:
- `n` is the length of the permutation
- `c` is the number of disjoint cycles

# Returns
- `1` for even permutations (even number of swaps)
- `-1` for odd permutations (odd number of swaps)
"""
function _permutation_parity(perm::VT) where {VT<:AbstractVector{<:Integer}}
    n = length(perm)
    n == 0 && return 1

    visited = falses(n)
    num_cycles = 0

    for i in 1:n
        if !visited[i]
            num_cycles += 1
            j = i
            while !visited[j]
                visited[j] = true
                j = perm[j]
            end
        end
    end

    return iseven(n - num_cycles) ? 1 : -1
end

"""
    _find_valid_contractions(word::Vector{T}) -> Vector{Tuple{Int,Int}}

Step 1 of Wick's algorithm: Identify all pairs with non-zero contractions.

For fermions, a contraction is non-zero only when:
- Position i has annihilation operator (positive index)
- Position j > i has creation operator (negative index)
- Both have the same mode (|word[i]| == |word[j]|)

This represents aᵢ followed by aᵢ† where {aᵢ, aᵢ†} = 1.
"""
function _find_valid_contractions(word::Vector{T}) where {T}
    contractions = Tuple{Int,Int}[]
    n = length(word)

    @inbounds for i in 1:n-1
        for j in i+1:n
            op_i = word[i]
            op_j = word[j]

            # Non-zero contraction: annihilation at i, creation at j, same mode
            if !_is_creation(op_i) && _is_creation(op_j) && _operator_mode(op_i) == _operator_mode(op_j)
                push!(contractions, (i, j))
            end
        end
    end

    return contractions
end

"""
    _generate_nonoverlapping_combinations(contractions::Vector{Tuple{Int,Int}}) -> Vector{Vector{Tuple{Int,Int}}}

Step 2 of Wick's algorithm: Generate all non-overlapping combinations.

Returns all subsets of contractions where no two pairs share an index.
Includes the empty set (no contractions).
"""
function _generate_nonoverlapping_combinations(contractions::Vector{Tuple{Int,Int}})
    combinations = Vector{Vector{Tuple{Int,Int}}}()
    push!(combinations, Tuple{Int,Int}[])  # Empty set (no contractions)

    _generate_combinations_recursive!(combinations, contractions, Tuple{Int,Int}[], Set{Int}(), 1)

    return combinations
end

function _generate_combinations_recursive!(
    result::Vector{Vector{Tuple{Int,Int}}},
    contractions::Vector{Tuple{Int,Int}},
    current::Vector{Tuple{Int,Int}},
    used_indices::Set{Int},
    start_idx::Int
)
    for i in start_idx:length(contractions)
        (idx1, idx2) = contractions[i]

        # Skip if either index is already used
        if idx1 ∈ used_indices || idx2 ∈ used_indices
            continue
        end

        # Add this contraction
        new_current = push!(copy(current), contractions[i])
        new_used = union(used_indices, Set([idx1, idx2]))

        push!(result, new_current)

        # Recurse for additional contractions
        _generate_combinations_recursive!(result, contractions, new_current, new_used, i + 1)
    end
end

"""
    _contraction_sign(word::Vector{T}, contraction::Vector{Tuple{Int,Int}}) where T

Compute the sign from bringing contracted pairs together.

For each contraction (i, j), we count the number of uncontracted operators
between positions i and j. Each such operator requires a swap to bring the
pair together, contributing a factor of -1.
"""
function _contraction_sign(word::Vector{T}, contraction::Vector{Tuple{Int,Int}}) where {T}
    isempty(contraction) && return 1

    # Sort contractions by first index for consistent processing
    sorted_contraction = sort(contraction, by=x -> x[1])

    # Track which positions are contracted
    contracted = Set{Int}()
    for (i, j) in contraction
        push!(contracted, i)
        push!(contracted, j)
    end

    sign = 1

    # For each contraction, count operators between i and j
    for (i, j) in sorted_contraction
        # Count uncontracted positions between i and j (exclusive)
        swaps = 0
        for k in (i+1):(j-1)
            if k ∉ contracted
                swaps += 1
            end
        end

        sign *= iseven(swaps) ? 1 : -1
    end

    return sign
end

"""
    _compute_normal_ordered_term(word::Vector{T}, contraction::Vector{Tuple{Int,Int}}) where T

Steps 3-4 of Wick's algorithm: Form N-product and evaluate sign.

For a given contraction combination:
1. Remove contracted operator pairs from the word
2. Normal-order remaining operators (creators left, annihilators right, sorted by mode)
3. Compute sign from the permutation required
4. Check for nilpotency (duplicate operators → zero)

Returns (coefficient, normal_ordered_word).
"""
function _compute_normal_ordered_term(word::Vector{T}, contraction::Vector{Tuple{Int,Int}}) where {T}
    # Collect indices to remove
    contracted_indices = Set{Int}()
    for (i, j) in contraction
        push!(contracted_indices, i)
        push!(contracted_indices, j)
    end

    # Build remaining operators with their original positions
    remaining = [(word[i], i) for i in 1:length(word) if i ∉ contracted_indices]

    if isempty(remaining)
        # All operators contracted - result is a scalar
        sign = _contraction_sign(word, contraction)
        return (sign * 1.0, T[])
    end

    # Check for duplicate operators (nilpotency: aᵢ² = 0)
    ops = [op for (op, _) in remaining]
    if length(ops) != length(Set(ops))
        return (0.0, T[])
    end

    # Sort to normal order: creators (negative) first by mode, then annihilators (positive) by mode
    sorted_remaining = sort(remaining, by=x -> normal_order_key(x[1]))

    # Extract the normal-ordered operators
    normal_word = T[op for (op, _) in sorted_remaining]

    # Compute sign from reordering the remaining operators
    n = length(remaining)
    remaining_positions = [pos for (_, pos) in remaining]

    # Build permutation: for each position in remaining, where does it end up in sorted?
    perm = zeros(Int, n)
    for (new_idx, (_, orig_pos)) in enumerate(sorted_remaining)
        old_idx = findfirst(==(orig_pos), remaining_positions)
        perm[old_idx] = new_idx
    end

    reorder_sign = _permutation_parity(perm)

    # Sign from contractions (bringing pairs together)
    contract_sign = _contraction_sign(word, contraction)

    total_sign = reorder_sign * contract_sign

    return (Float64(total_sign), normal_word)
end

# Note: Like-term combination happens via `_combine_physics_terms` below.

"""
    _simplify_fermionic_word!(word::Vector{T}) where {T<:Integer} -> Vector{Tuple{Int,Vector{T}}}

Normal-order a fermionic word using Generalized Wick's Theorem.

Returns a vector of (coefficient, normal_ordered_word) pairs representing the sum:
  Σ coeffs[i] * word[i]

This is the low-level function used by `simplify(FermionicAlgebra, word)`.
Integer coefficients are exact (from anticommutation signs).

# Algorithm
1. Find all valid contractions (aᵢ...aᵢ† pairs with same mode)
2. Generate all non-overlapping contraction combinations
3. For each combination, compute the normal-ordered term with sign
4. Combine like terms
"""
function _simplify_fermionic_word!(word::Vector{T}) where {T<:Integer}
    # Handle empty word
    if isempty(word)
        return Tuple{Int,Vector{T}}[(1, T[])]
    end

    # Check for nilpotency using flux-based detection.
    # A fermionic word is zero if any mode has net surplus >= 2 (same-type operators).
    if _is_nilpotent_by_flux(word)
        return Tuple{Int,Vector{T}}[(0, T[])]
    end

    # Step 1: Find valid contractions
    contractions = _find_valid_contractions(word)

    # Step 2: Generate non-overlapping combinations
    combinations = _generate_nonoverlapping_combinations(contractions)

    # Steps 3-4: Compute each term
    result_pairs = Tuple{Int,Vector{T}}[]

    for combo in combinations
        (coef, normal_word) = _compute_normal_ordered_term(word, combo)
        coef_int = Int(coef)
        if coef_int != 0
            push!(result_pairs, (coef_int, normal_word))
        end
    end

    # Combine like terms
    grouped = Dict{Vector{T},Int}()
    for (coef, w) in result_pairs
        grouped[w] = get(grouped, w, 0) + coef
    end

    result = Tuple{Int,Vector{T}}[]
    for (w, coef) in grouped
        if coef != 0
            push!(result, (coef, w))
        end
    end

    if isempty(result)
        push!(result, (0, T[]))
    end

    return result
end

"""
    simplify(::Type{FermionicAlgebra}, word::Vector{T}) where {T<:Signed} -> Vector{Tuple{Int,Vector{T}}}

Simplify a raw fermionic word into normal-ordered form.

This is the primary entry point for fermionic simplification. Takes a raw word vector
and returns a vector of `(coefficient, normal_ordered_word)` pairs representing the
canonical normal-ordered representation.
"""
simplify(::Type{FermionicAlgebra}, word::Vector{T}) where {T<:Signed} = simplify!(FermionicAlgebra, copy(word))

"""
    simplify!(::Type{FermionicAlgebra}, word::Vector{T}) where {T<:Signed} -> Vector{Tuple{Int,Vector{T}}}

Low-level in-place normal ordering of a fermionic word.

Returns a vector of `(coefficient, normal_ordered_word)` pairs representing the
PBW expansion as a sum. Each word in the result is unique (coefficients are
accumulated for identical words). Modifies `word` in-place (filters zeros).

This is the internal workhorse; most users should call `simplify(FermionicAlgebra, word)`.

# Algorithm
Uses Generalized Wick's Theorem:
1. Find all valid contractions (aᵢ...aᵢ† pairs with same mode)
2. Generate all non-overlapping contraction combinations
3. For each combination, compute the normal-ordered term with sign
4. Combine like terms

# Throws
- `ErrorException` if simplification produces an empty result (should not happen for valid input)
"""
function simplify!(::Type{FermionicAlgebra}, word::Vector{T}) where {T<:Signed}
    filter!(!iszero, word)

    # Handle empty word (identity)
    if isempty(word)
        return [(1, T[])]
    end

    # Early exit for nilpotent monomials using flux-based detection
    if _is_nilpotent_by_flux(word)
        return [(0, T[])]
    end

    # Step 1: Find valid contractions
    contractions = _find_valid_contractions(word)

    # Step 2: Generate non-overlapping combinations
    combinations = _generate_nonoverlapping_combinations(contractions)

    # Steps 3-4: Compute each term and combine like terms
    word_coeffs = Dict{Vector{T},Int}()

    for combo in combinations
        (coef, normal_word) = _compute_normal_ordered_term(word, combo)
        int_coef = round(Int, coef)  # fermionic coefficients are always ±1
        if int_coef != 0
            word_coeffs[normal_word] = get(word_coeffs, normal_word, 0) + int_coef
        end
    end

    # Construct final (coef, word) pairs, filtering zero coefficients
    result = Tuple{Int,Vector{T}}[]
    for (w, coef) in word_coeffs
        coef == 0 && continue
        push!(result, (coef, w))
    end

    isempty(result) && return [(0, T[])]

    return result
end

"""
    _is_fermionic_nilpotent(word::Vector{T}) where {T<:Signed} -> Bool

Check if a fermionic word contains repeated identical operators (nilpotent: aᵢ² = 0 or (a†ᵢ)² = 0).

Note: aᵢ and a†ᵢ are DIFFERENT operators (indices 1 and -1 for mode 1), so
[1, -1] is NOT nilpotent, but [1, 1] or [-1, -1] are.
"""
function _is_fermionic_nilpotent(word::Vector{T}) where {T<:Signed}
    seen = Set{T}()
    for idx in word
        idx ∈ seen && return true
        push!(seen, idx)
    end
    return false
end

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_fermionic_word!(word::Vector{T}) where {T<:Signed}

Check that a fermionic word is in normal-ordered form. Throws `ArgumentError` if invalid.

Normal order requirements:
- Creation operators (negative indices) come before annihilation operators (positive)
- Creators sorted by mode descending; annihilators sorted by mode ascending

This is used by `NormalMonomial{FermionicAlgebra,T}` constructor to enforce invariants.
"""
function _validate_word(::Type{FermionicAlgebra}, word::Vector{T}) where {T<:Signed}
    isempty(word) && return nothing

    if !is_normal_ordered(word)
        throw(ArgumentError(
            "Fermionic word is not normal-ordered. " *
            "Use simplify(FermionicAlgebra, word) for auto-normal-ordering."
        ))
    end

    # Normal-form monomials must be non-nilpotent: aᵢ² = 0 and (aᵢ†)² = 0.
    # In normal order (creators then annihilators, each sorted by mode), duplicate operators
    # can only appear as adjacent equal entries.
    @inbounds for i in 2:length(word)
        if word[i] == word[i - 1]
            throw(ArgumentError(
                "Fermionic word is nilpotent (duplicate operator). " *
                "Use simplify(FermionicAlgebra, word) for canonicalization."
            ))
        end
    end
    return nothing
end
