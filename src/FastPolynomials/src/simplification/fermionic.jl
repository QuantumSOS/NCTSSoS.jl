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
- Within groups: sorted by mode (absolute value)

# Wick's Theorem Algorithm (4 steps)
1. Identify valid contractions (pairs with non-zero contraction values)
2. Generate all non-overlapping combinations of contractions
3. Form N-products for each combination (normal-ordered remaining operators)
4. Evaluate signs using permutation parity
"""

# Helper functions
@inline _is_creation(op::Integer) = op < 0
@inline _fermi_mode(op::Integer) = abs(op)

"""
    Base.iszero(m::Monomial{FermionicAlgebra,T}) -> Bool

Check if a fermionic monomial is zero due to nilpotency.

# Algorithm
Scan through operators tracking the last seen operator type for each mode.
If we see the same operator type (creation or annihilation) for the same mode
twice without an intervening opposite type, the monomial is zero.

This works because consecutive same-type operators for a mode cannot be
contracted away and will produce nilpotent terms (aᵢ² = 0, (aᵢ†)² = 0).

# Examples
```jldoctest
julia> m = Monomial{FermionicAlgebra}(Int32[1, 1]);  # a₁ a₁ = 0

julia> iszero(m)
true

julia> m2 = Monomial{FermionicAlgebra}(Int32[1, -1, 1, -1]);  # a₁ a₁† a₁ a₁† ≠ 0

julia> iszero(m2)
false

julia> m3 = Monomial{FermionicAlgebra}(Int32[1, -1, 1]);  # a₁ a₁† a₁ = a₁ ≠ 0

julia> iszero(m3)
false
```
"""
function Base.iszero(m::Monomial{FermionicAlgebra,T}) where {T}
    word = m.word
    isempty(word) && return false

    # Find max mode to size the tracking array
    max_mode = maximum(_fermi_mode, word)

    # Track last seen operator type: 0 = unseen, 1 = annihilation, -1 = creation
    last_seen = zeros(Int8, max_mode)

    for op in word
        mode = _fermi_mode(op)
        op_type = _is_creation(op) ? Int8(-1) : Int8(1)

        if last_seen[mode] == op_type
            return true  # Same type twice without intervening opposite
        end
        last_seen[mode] = op_type
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

    for i in 1:n-1
        for j in i+1:n
            op_i = word[i]
            op_j = word[j]

            # Non-zero contraction: annihilation at i, creation at j, same mode
            if !_is_creation(op_i) && _is_creation(op_j) && _fermi_mode(op_i) == _fermi_mode(op_j)
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
    sorted_remaining = sort(remaining, by=x -> (_is_creation(x[1]) ? 0 : 1, _fermi_mode(x[1])))

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

"""
    _combine_like_terms_fermi(terms::Vector{Term}) -> Vector{Term}

Combine terms with identical monomials by summing coefficients.
"""
function _combine_like_terms_fermi(terms::Vector{Term{Monomial{FermionicAlgebra,T},Float64}}) where {T}
    isempty(terms) && return [Term(0.0, Monomial{FermionicAlgebra}(T[]))]

    grouped = Dict{Vector{T},Float64}()
    for t in terms
        key = t.monomial.word
        grouped[key] = get(grouped, key, 0.0) + t.coefficient
    end

    result = Term{Monomial{FermionicAlgebra,T},Float64}[]
    for (word, coef) in grouped
        if coef != 0.0
            push!(result, Term(coef, Monomial{FermionicAlgebra}(word)))
        end
    end

    if isempty(result)
        push!(result, Term(0.0, Monomial{FermionicAlgebra}(T[])))
    end

    return result
end

"""
    simplify!(m::Monomial{FermionicAlgebra,T}) where T -> Polynomial{FermionicAlgebra,T,Float64}

Simplify fermionic monomial using Generalized Wick's Theorem.

# Algorithm
1. Find all valid contractions (aᵢ...aᵢ† pairs with same mode)
2. Generate all non-overlapping contraction combinations
3. For each combination, compute the normal-ordered term with sign
4. Combine like terms

# Returns
Polynomial representing the normal-ordered expansion (potentially with multiple terms due to anticommutation).
"""
function simplify!(m::Monomial{FermionicAlgebra,T}) where {T}
    word = m.word

    # Handle empty word
    if isempty(word)
        return Polynomial([Term(1.0, Monomial{FermionicAlgebra}(T[]))])
    end

    # Early exit for nilpotent monomials (aᵢ² = 0)
    if iszero(m)
        return Polynomial([Term(0.0, Monomial{FermionicAlgebra}(T[]))])
    end

    # Step 1: Find valid contractions
    contractions = _find_valid_contractions(word)

    # Step 2: Generate non-overlapping combinations
    combinations = _generate_nonoverlapping_combinations(contractions)

    # Steps 3-4: Compute each term
    result_terms = Term{Monomial{FermionicAlgebra,T},Float64}[]

    for combo in combinations
        (coef, normal_word) = _compute_normal_ordered_term(word, combo)
        if coef != 0.0
            push!(result_terms, Term(coef, Monomial{FermionicAlgebra}(normal_word)))
        end
    end

    # Combine like terms
    result_terms = _combine_like_terms_fermi(result_terms)

    return Polynomial(result_terms)
end

"""
    simplify(m::Monomial{FermionicAlgebra,T}) where T -> Polynomial{FermionicAlgebra,T,Float64}

Non-mutating version of simplify!.
Returns a Polynomial representing the normal-ordered expansion (potentially with multiple terms due to anticommutation).
"""
function simplify(m::Monomial{FermionicAlgebra,T}) where {T}
    word_copy = copy(m.word)
    m_copy = Monomial{FermionicAlgebra,T}(word_copy, m.hash)
    simplify!(m_copy)
end
