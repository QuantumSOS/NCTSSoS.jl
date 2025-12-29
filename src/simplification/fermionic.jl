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

# Note: Helper functions (_is_creation, _operator_mode, normal_order_key,
# combine_like_terms) are in utils.jl

"""
    has_even_parity(m::Monomial{FermionicAlgebra,T}) where T -> Bool

Check if a fermionic monomial has even parity (even number of operators).

In fermionic systems, parity superselection forbids observables with odd fermion
number from having non-zero expectation values. Only operators with even parity
(even number of creation/annihilation operators) can have non-zero physical
expectation values.

# Arguments
- `m::Monomial{FermionicAlgebra,T}`: A fermionic monomial

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
function has_even_parity(m::Monomial{FermionicAlgebra,T}) where {T}
    return iseven(length(m.word))
end

"""
    Base.iszero(m::Monomial{FermionicAlgebra,T}) -> Bool

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
julia> m = Monomial{FermionicAlgebra}(Int32[1, 1]);  # a₁ a₁ = 0 (surplus 2)

julia> iszero(m)
true

julia> m2 = Monomial{FermionicAlgebra}(Int32[1, 2, 1]);  # a₁ a₂ a₁ = -a₁ a₁ a₂ = 0

julia> iszero(m2)
true

julia> m3 = Monomial{FermionicAlgebra}(Int32[1, 1, -1]);  # a₁ a₁ a₁† (surplus 1, not zero yet)

julia> iszero(m3)
false

julia> m4 = Monomial{FermionicAlgebra}(Int32[1, -1, 1, -1]);  # a₁ a₁† a₁ a₁† ≠ 0

julia> iszero(m4)
false
```
"""
function Base.iszero(m::Monomial{FermionicAlgebra,T}) where {T}
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

# Note: combine_like_terms is now in utils.jl as a shared helper

"""
    simplify(m::Monomial{FermionicAlgebra,T}) where T -> Polynomial{FermionicAlgebra,T,Float64}

Simplify and normal-order a fermionic algebra monomial using Generalized Wick's Theorem.

Returns a Polynomial representing the normal-ordered expansion (potentially with multiple
terms due to anticommutation). The original monomial is unchanged.

# Algorithm
1. Find all valid contractions (aᵢ...aᵢ† pairs with same mode)
2. Generate all non-overlapping contraction combinations
3. For each combination, compute the normal-ordered term with sign
4. Combine like terms

# Algebraic Rules
- {aᵢ, aⱼ†} = δᵢⱼ (anticommutation gives delta correction)
- {aᵢ, aⱼ} = 0 (annihilation-annihilation anticommute)
- {aᵢ†, aⱼ†} = 0 (creation-creation anticommute)
- aᵢ² = 0, (aᵢ†)² = 0 (nilpotency)

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{FermionicAlgebra}(Int32[1, -1]);  # a₁ a₁†

julia> poly = simplify(m);

julia> length(terms(poly))  # Two terms: 1 - a₁† a₁
2
```

# Note
Unlike other algebra types, fermionic simplification returns a `Polynomial` (not a `Monomial`)
because anticommutation creates multiple terms. There is no `simplify!` variant since the
return type differs from the input type.
"""
function simplify(m::Monomial{FermionicAlgebra,T}) where {T}
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
    result_terms = combine_like_terms(result_terms)

    return Polynomial(result_terms)
end
