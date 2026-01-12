"""
    Bosonic Algebra Simplification

Implements normal ordering for bosonic creation/annihilation operators satisfying:
- [cᵢ, cⱼ†] = δᵢⱼ (creation-annihilation commute with delta correction)
- [cᵢ, cⱼ] = 0 (annihilation-annihilation commute)
- [cᵢ†, cⱼ†] = 0 (creation-creation commute)

# Variable Encoding Convention
Same as fermionic:
- Annihilation operator cᵢ: positive index `i`
- Creation operator cᵢ†: negative index `-i`

# Normal Ordering
All creation operators (c†, negative indices) are placed to the LEFT of
all annihilation operators (c, positive indices). Within each group,
creators are sorted by mode descending and annihilators by mode ascending.

# Key Difference from Fermionic
Unlike fermionic operators:
1. Bosonic operators COMMUTE (no sign change when swapping)
2. But cᵢ cⱼ† = cⱼ† cᵢ + δᵢⱼ creates an additional term when modes match
3. Bosonic operators are NOT nilpotent (cᵢ cᵢ ≠ 0)

# Algorithm
Uses a group-based algorithm with closed-form formulas based on:
**"Rook numbers and the normal ordering problem"**

1. **Group by mode**: Partition operators by mode, preserving relative order
2. **Single-mode normal form**: For each mode, compute normal form using
   rook numbers on Ferrers boards (Eq. 1.40 from the paper)
3. **Product expansion**: Expand product across all mode groups
4. **Construct terms**: Build final monomial words from the expansion

# Mathematical Foundation
For a word ω in a single mode with m creations and n annihilations:
  ω = Σₖ₌₀^{min(m,n)} rₖ(Bω) (a†)^{m-k} a^{n-k}

where rₖ(Bω) is the k-th rook number of the Ferrers board Bω.
The board cell (i,j) is included if creation i appears AFTER annihilation j,
representing a contraction that produces a delta term.

# Complexity Analysis
- **Grouping**: O(n) single pass to partition by mode
- **Single-mode**: O(K²) per mode where K = operators per mode
- **Product expansion**: O(∏ᵢ |terms_i|) where |terms_i| = terms per mode
- **Total**: O(n + M×K² + R) where M = modes, R = result size

This is polynomial time for most practical cases, compared to O(2^n) for
naive iterative expansion.

**TODO**: Verify complexity bounds against literature.

# Returns
A vector of `(coefficient, word)` pairs representing the PBW expansion.

# Examples
```jldoctest
julia> using NCTSSoS

julia> word = Int32[1, -1];  # c₁ c₁† (raw word; not in Bosonic normal form)

julia> m = simplify(BosonicAlgebra, word);

julia> length(terms(m))  # Two terms: c₁† c₁ and identity
2
```
"""

"""
    simplify!(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed} -> Vector{Tuple{Int,Vector{T}}}

Low-level in-place normal ordering of a bosonic word.

Returns a vector of `(coefficient, normal_ordered_word)` pairs representing the
PBW expansion as a sum. Each word in the result is unique (coefficients are
accumulated for identical words). Modifies `word` in-place (filters zeros, sorts by mode).

This is the internal workhorse; most users should call `simplify(BosonicAlgebra, word)`.

# Algorithm
Uses rook numbers on Ferrers boards (arXiv:quant-ph/0507206):
1. Sort by mode (preserving relative order within each mode)
2. Find mode boundaries
3. Compute single-mode normal forms using rook number formula
4. Expand Cartesian product across all modes
5. Combine like terms by final word and construct result pairs

# Throws
- `ErrorException` if simplification produces an empty result (should not happen for valid input)
"""
function simplify!(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed}
    filter!(!iszero, word)

    # Handle empty word (identity)
    if isempty(word)
        return [(1, T[])]
    end

    _stable_sort_by_site!(word)

    # Step 1: Group by mode (stable sort)
    sep = find_site_separation(word)
    n_modes = length(sep) - 1

    # Step 2: Compute single-mode normal forms
    modes = Vector{Int}(undef, n_modes)
    mode_results = Vector{Tuple{Int,Vector{Tuple{Int,Int,Int}}}}(undef, n_modes)
    for i in 1:n_modes
        modes[i] = _operator_mode(word[sep[i]+1])
        ops = @view word[sep[i]+1:sep[i+1]]
        mode_results[i] = (i, single_mode_normal_form(ops))
    end

    # Step 3: Expand product - but return (Int, Vector{T}) pairs instead of Terms
    term_lists = [terms for (_, terms) in mode_results]

    # Cartesian product of all mode terms, combine like terms by final word
    word_coeffs = Dict{Vector{T},Int}()
    for combo in Iterators.product(term_lists...)
        coef = prod(t[1] for t in combo)
        coef == 0 && continue

        # Build the normal-ordered word
        new_word = T[]
        for i in reverse(eachindex(modes))
            m = modes[i]
            append!(new_word, fill(T(-m), combo[i][2]))
        end
        for (i, m) in enumerate(modes)
            append!(new_word, fill(T(m), combo[i][3]))
        end

        word_coeffs[new_word] = get(word_coeffs, new_word, 0) + coef
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
    simplify(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed}

Simplify a raw bosonic word into normal-ordered form.

This is the primary entry point for bosonic simplification. Takes a raw word vector
and returns a vector of `(coefficient, word)` pairs representing the PBW expansion.
"""
simplify(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed} = simplify!(BosonicAlgebra, copy(word))

# =============================================================================
# Group-Based Algorithm with Closed-Form Formulas (Rook Numbers)
# Based on arXiv:quant-ph/0507206 - Combinatorics of boson normal ordering
# =============================================================================

"""
    find_site_separation(word::Vector{T}) -> Vector{Int}

Find mode boundaries in a pre-sorted word.

Assumes `word` is already sorted by mode (via `_stable_sort_by_site!`).
Returns separation indices `[0, end_mode1, end_mode2, ...]` where
operators for mode i are at indices `sep[i]+1:sep[i+1]`.

# Example
```julia
word = Int32[1, -1, -2, 2]  # Already sorted: mode 1 ops, then mode 2 ops
sep = find_site_separation(word)
# sep = [0, 2, 4]
```
"""
function find_site_separation(word::Vector{T})::Vector{Int} where {T<:Signed}
    isempty(word) && return [0]

    # Single sweep to find partition boundaries
    sep = [0]
    current_mode = _operator_mode(word[1])
    @inbounds for i in 2:length(word)
        m = _operator_mode(word[i])
        if m != current_mode
            push!(sep, i - 1)
            current_mode = m
        end
    end
    push!(sep, length(word))

    return sep
end

"""
    build_ferrers_board(ops::AbstractVector{T}) -> Vector{Int}

Build a Ferrers board from an operator sequence (single mode).
Returns column heights where board[j] = number of creations AFTER annihilation j.

The board has n columns (one per annihilation) and m rows (one per creation).
Cell (i, j) is included if creation i appears AFTER annihilation j in the word.
This corresponds to pairs that need contraction (swapping a with a†).

# Example
```julia
ops = Int32[1, -1]  # a a† (annihilation at 1, creation at 2)
board = build_ferrers_board(ops)
# Creation is AFTER annihilation → contraction possible → height 1
# board = [1]

ops = Int32[-1, 1]  # a† a (creation at 1, annihilation at 2)
board = build_ferrers_board(ops)
# Creation is BEFORE annihilation → no contraction → height 0
# board = [0]
```
"""
function build_ferrers_board(ops::V)::Vector{Int} where {T, V<:AbstractVector{T}}
    # Count annihilations (positive = annihilation)
    n_annihilations = count(op -> !_is_creation(op), ops)
    n_annihilations == 0 && return Int[]

    total_creations = length(ops) - n_annihilations

    # Build board: for each annihilation, count creations that appear AFTER it
    board = Vector{Int}(undef, n_annihilations)
    ann_idx = 0
    creations_seen = 0

    for op in ops
        if _is_creation(op)
            creations_seen += 1
        else
            ann_idx += 1
            # Creations after this annihilation = total - creations seen so far
            board[ann_idx] = total_creations - creations_seen
        end
    end

    return board
end

"""
    compute_rook_numbers(board::Vector{Int}) -> Vector{Int}

Compute rook numbers for a Ferrers board using dynamic programming.
Returns rook_numbers where rook_numbers[k+1] = rₖ(B) (1-indexed).

rₖ(B) = number of ways to place k non-attacking rooks on board B.

# Example
```julia
board = [1, 1]  # Two columns, each with height 1 (one row)
rook = compute_rook_numbers(board)
# r₀ = 1 (one way to place 0 rooks)
# r₁ = 2 (two cells to place one rook)
# r₂ = 0 (can't place 2 rooks in 1 row, but we only track up to n_rows)
# rook = [1, 2]
```
"""
function compute_rook_numbers(board::Vector{Int})::Vector{Int}
    n_cols = length(board)
    n_cols == 0 && return [1]  # r₀ = 1

    n_rows = maximum(board; init=0)
    n_rows == 0 && return [1]  # No cells, only r₀ = 1

    max_rooks = min(n_cols, n_rows)

    # DP: rook[k+1] = number of ways to place k rooks
    rook = zeros(Int, max_rooks + 1)
    rook[1] = 1  # r₀ = 1

    for j in 1:n_cols
        h = board[j]
        # Process in reverse to avoid using updated values
        for k in min(j, max_rooks):-1:1
            # When placing k-th rook in column j:
            # - h cells available in column
            # - (k-1) rows already occupied by previous rooks
            # - available = max(0, h - (k-1))
            available = max(0, h - (k - 1))
            rook[k+1] += rook[k] * available
        end
    end

    # Trim trailing zeros
    last_nonzero = findlast(!=(0), rook)
    return isnothing(last_nonzero) ? [1] : rook[1:last_nonzero]
end

"""
    single_mode_normal_form(ops::AbstractVector{T}) -> Vector{Tuple{Int,Int,Int}}

Compute the normal-ordered form of a single-mode operator sequence.
Uses the rook number formula from Eq. 1.40 of arXiv:quant-ph/0507206.

Returns vector of (coefficient, num_creations, num_annihilations) tuples.
Each tuple represents a term in the normal-ordered expansion.

# Formula
For word ω with m creations and n annihilations:
  ω = Σₖ rₖ(Bω) (a†)^{m-k} a^{n-k}
where Bω is the Ferrers board and rₖ is the k-th rook number.

# Example
```julia
ops = Int32[1, -1]  # a a†
result = single_mode_normal_form(ops)
# result = [(1, 1, 1), (1, 0, 0)]  # a† a + 1
```
"""
function single_mode_normal_form(ops::V)::Vector{Tuple{Int,Int,Int}} where {T<:Signed, V<:AbstractVector{T}}
    # Count creations and annihilations
    m = count(_is_creation, ops)
    n = length(ops) - m  # annihilations

    # Edge cases: no contractions possible
    if m == 0
        return [(1, 0, n)]  # Just annihilations
    end
    if n == 0
        return [(1, m, 0)]  # Just creations
    end

    # Build Ferrers board and compute rook numbers
    board = build_ferrers_board(ops)
    rook = compute_rook_numbers(board)

    # Build result: for each k, add term rₖ × (a†)^{m-k} a^{n-k}
    result = Tuple{Int,Int,Int}[]
    for k in 0:length(rook)-1
        coef = rook[k+1]
        if coef != 0
            push!(result, (coef, m - k, n - k))
        end
    end

    return result
end

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_word(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed}

Validate that a bosonic word is in normal-ordered form.

Throws `ArgumentError` if the word violates normal ordering. Called by the
`NormalMonomial{BosonicAlgebra,T}` constructor to enforce invariants.

Normal order requirements:
- Creation operators (negative indices) before annihilation operators (positive)
- Creators sorted by mode descending; annihilators sorted by mode ascending
"""
function _validate_word(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed}
    isempty(word) && return nothing

    is_normal_ordered(word) ||
        throw(ArgumentError(
            "Bosonic word is not normal-ordered. " *
            "Use simplify(BosonicAlgebra, word) for auto-normal-ordering."
        ))
    return nothing
end
