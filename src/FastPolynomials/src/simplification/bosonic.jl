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
operators are sorted by mode (absolute value of index).

# Key Difference from Fermionic
Unlike fermionic operators:
1. Bosonic operators COMMUTE (no sign change when swapping)
2. But cᵢ cⱼ† = cⱼ† cᵢ + δᵢⱼ creates an additional term when modes match
3. Bosonic operators are NOT nilpotent (cᵢ cᵢ ≠ 0)
4. **Returns `Vector{Term}` instead of single `Term`** due to term generation

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
`Vector{Term{Monomial{BosonicAlgebra,T},Float64}}` - multiple terms due to delta corrections

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{BosonicAlgebra}(Int32[1, -1]);  # c₁ c₁†

julia> terms = simplify(m);

julia> length(terms)  # Two terms: c₁† c₁ and identity
2
```
"""

# Encoding helper functions (exported for testing)

"""
    bosonic_mode(op::Integer) -> Int

Extract the mode (site) index from a bosonic operator.
Same as fermionic.
"""
@inline bosonic_mode(op::Integer) = abs(op)

"""
    bosonic_sort_key(op::T) where T

Compute sort key for bosonic operators.
Creations (negative) come first, then annihilations (positive).
Within each type, sort by mode in ascending order.
"""
@inline function bosonic_sort_key(op)
    is_cr = _is_creation(op)
    mode = bosonic_mode(op)
    return (is_cr ? 0 : 1, mode)
end

"""
    is_normal_ordered(word::Vector{T}) where T -> Bool

Check if a bosonic word is already in normal order.
"""
function is_normal_ordered(word::Vector{T}) where {T}
    length(word) <= 1 && return true

    for i in 1:length(word)-1
        key_i = bosonic_sort_key(word[i])
        key_i1 = bosonic_sort_key(word[i+1])
        if key_i > key_i1
            return false
        end
    end
    return true
end

"""
    find_first_out_of_order_bosonic(word::Vector{T}) where T -> Int

Find the first position where operators are out of normal order.
Returns the index, or 0 if already in normal order (type-stable).
"""
function find_first_out_of_order_bosonic(word::Vector{T}) where {T}
    for i in 1:length(word)-1
        key_i = bosonic_sort_key(word[i])
        key_i1 = bosonic_sort_key(word[i+1])
        if key_i > key_i1
            return i
        end
    end
    return 0  # Type-stable: returns Int, not Nothing
end

"""
    simplify!(m::Monomial{BosonicAlgebra,T}) where T -> Polynomial{BosonicAlgebra,T,Float64}

In-place simplification and normal ordering of a bosonic algebra monomial.
Returns a Polynomial (potentially with multiple terms) due to delta corrections from commutation.

Uses a group-based algorithm with closed-form formulas (rook numbers on Ferrers boards)
based on . This is more efficient than iterative expansion,
especially for expressions with multiple modes.

**Algorithm:**
1. Group operators by mode (preserving relative order within each mode)
2. For each mode, compute normal form using rook number formula
3. Expand product across all modes
4. Construct final terms

**Note:** Returns a Polynomial because [c, c†] = 1 creates delta corrections:
c c† = c† c + 1

# Examples
```jldoctest
julia> m = Monomial{BosonicAlgebra}(Int32[1, -1]);  # c₁ c₁†

julia> poly = simplify!(m);

julia> length(terms(poly))
2
```
"""
function simplify!(m::Monomial{BosonicAlgebra,T}) where {T}
    term_vec = simplify_bosonic_grouped!(m)
    return Polynomial(term_vec)
end

"""
    simplify(m::Monomial{BosonicAlgebra,T}) where T -> Polynomial{BosonicAlgebra,T,Float64}

Simplify and normal order a bosonic algebra monomial.

Non-mutating version - creates a copy and simplifies it.
Returns a Polynomial (potentially with multiple terms) due to delta corrections from commutation.

# Algebraic Rules
- Commutation: [cᵢ, cⱼ†] = δᵢⱼ (swapping creates delta term when modes match)
- Bosons commute: [cᵢ, cⱼ] = 0, [cᵢ†, cⱼ†] = 0 (no sign change)
- Normal ordering: c† on left, c on right
- **NOT nilpotent**: cᵢ cᵢ ≠ 0 (unlike fermions)

# Examples
```jldoctest
julia> m = Monomial{BosonicAlgebra}(Int32[2, -1]);  # c₂ c₁†

julia> poly = simplify(m);

julia> length(terms(poly))  # Just one term (different modes, no delta)
1

julia> monomials(poly)[1].word
2-element Vector{Int32}:
 -1
  2
```
"""
function simplify(m::Monomial{BosonicAlgebra,T}) where {T}
    # Create a copy of the word to work with
    word_copy = copy(m.word)
    m_copy = Monomial{BosonicAlgebra,T}(word_copy)
    simplify!(m_copy)
end

# =============================================================================
# Group-Based Algorithm with Closed-Form Formulas (Rook Numbers)
# Based on arXiv:quant-ph/0507206 - Combinatorics of boson normal ordering
# =============================================================================

"""
    group_by_mode!(word::Vector{T}) -> Vector{Int}

Group operators by mode in-place using stable partitioning.
Preserves the relative order of operators within each mode.

Returns separation indices [0, end_mode1, end_mode2, ...] where
operators for mode m are at indices sep[m]+1:sep[m+1].

# Example
```julia
word = Int32[1, -2, -1, 2]  # a₁ a₂† a₁† a₂
sep = group_by_mode!(word)
# word becomes [1, -1, -2, 2] (mode 1 ops, then mode 2 ops)
# sep = [0, 2, 4]
```
"""
function group_by_mode!(word::Vector{T})::Vector{Int} where {T}
    isempty(word) && return [0]

    # Stable sort by mode - groups operators while preserving relative order
    sort!(word; by=bosonic_mode, alg=Base.InsertionSort)

    # Single sweep to find partition boundaries
    sep = [0]
    current_mode = bosonic_mode(word[1])
    for i in 2:length(word)
        m = bosonic_mode(word[i])
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
function single_mode_normal_form(ops::V)::Vector{Tuple{Int,Int,Int}} where {T, V<:AbstractVector{T}}
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

"""
    expand_and_construct(mode_results, modes, ::Type{T}) -> Vector{Term}

Expand the product of normal-ordered terms across all modes and construct final Terms.

# Arguments
- `mode_results`: Vector of (mode_idx, terms) where terms = [(coef, num_creations, num_annihilations), ...]
- `modes`: Sorted vector of actual mode indices
- `T`: Integer type for the monomial word

# Example
```julia
mode_results = [(1, [(1, 1, 1), (1, 0, 0)]), (2, [(1, 1, 0)])]  # (a₁†a₁ + 1) × a₂†
modes = [1, 2]
result = expand_and_construct(mode_results, modes, Int32)
# Returns terms for: a₁†a₂†a₁ + a₂†
```
"""
function expand_and_construct(
    mode_results::Vector{Tuple{Int,Vector{Tuple{Int,Int,Int}}}},
    modes::Vector{Int},
    ::Type{T}
)::Vector{Term{Monomial{BosonicAlgebra,T},Float64}} where {T}
    @assert !isempty(mode_results) "mode_results must be non-empty (caller handles empty case)"
    @assert length(modes) == length(mode_results) "modes and mode_results must have same length"

    n_modes = length(modes)
    term_lists = [terms for (_, terms) in mode_results]

    # Cartesian product of all mode terms, combine like terms
    combined = Dict{Tuple{Vector{Int},Vector{Int}},Int}()
    for combo in Iterators.product(term_lists...)
        coef = prod(t[1] for t in combo)
        coef == 0 && continue
        cre = [combo[i][2] for i in 1:n_modes]
        ann = [combo[i][3] for i in 1:n_modes]
        key = (cre, ann)
        combined[key] = get(combined, key, 0) + coef
    end

    # Construct final Term vector
    result = Term{Monomial{BosonicAlgebra,T},Float64}[]
    for ((cre, ann), coef) in combined
        coef == 0 && continue
        word = T[]
        for (i, m) in enumerate(modes)
            append!(word, fill(T(-m), cre[i]))
        end
        for (i, m) in enumerate(modes)
            append!(word, fill(T(m), ann[i]))
        end
        push!(result, Term(Float64(coef), Monomial{BosonicAlgebra}(word)))
    end

    return result
end

"""
    simplify_bosonic_grouped!(m::Monomial{BosonicAlgebra,T}) -> Vector{Term}

Simplify a bosonic monomial using the group-based algorithm with rook numbers.

Algorithm:
1. Group operators by mode (preserving relative order within each mode)
2. For each mode, compute normal form using rook number formula
3. Expand product across all modes
4. Construct final terms

This is more efficient than iterative expansion for multi-mode expressions.
"""
function simplify_bosonic_grouped!(
    m::Monomial{BosonicAlgebra,T}
)::Vector{Term{Monomial{BosonicAlgebra,T},Float64}} where {T}
    word = m.word

    # Handle empty monomial
    if isempty(word)
        return [Term(1.0, m)]
    end

    # Step 1: Group by mode (stable sort)
    sep = group_by_mode!(word)
    n_modes = length(sep) - 1

    # Step 2: Compute single-mode normal forms
    modes = Vector{Int}(undef, n_modes)
    mode_results = Vector{Tuple{Int,Vector{Tuple{Int,Int,Int}}}}(undef, n_modes)
    for i in 1:n_modes
        modes[i] = bosonic_mode(word[sep[i]+1])
        ops = @view word[sep[i]+1:sep[i+1]]
        mode_results[i] = (i, single_mode_normal_form(ops))
    end

    # Step 3: Expand product and construct terms
    return expand_and_construct(mode_results, modes, T)
end
