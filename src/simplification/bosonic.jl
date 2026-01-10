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
`Monomial{BosonicAlgebra}` - a canonical PBW expansion (integer coefficients internally).
Convert to numeric coefficients with `Polynomial(m)` if needed.

# Examples
```jldoctest
julia> using NCTSSoS

julia> word = Int32[1, -1];  # c₁ c₁† (raw word; not in Bosonic normal form)

julia> m = simplify(BosonicAlgebra, word);

julia> length(terms(m))  # Two terms: c₁† c₁ and identity
2
```
"""

# Note: Shared helpers (_is_creation, _operator_mode, normal_order_key,
# find_first_out_of_order, is_normal_ordered) are in `src/util/helpers.jl`.

"""
    simplify(m::NormalMonomial{BosonicAlgebra,T}) where T -> Monomial{BosonicAlgebra}

Simplify and normal-order a bosonic algebra monomial.

Returns a `Monomial` (canonical PBW expansion). Bosonic commutation may introduce
additional terms; coefficients are stored internally as exact integers.

Uses a group-based algorithm with closed-form formulas (rook numbers on Ferrers boards)
based on arXiv:quant-ph/0507206. This is more efficient than iterative expansion,
especially for expressions with multiple modes.

# Algebraic Rules
- Commutation: [cᵢ, cⱼ†] = δᵢⱼ (swapping creates delta term when modes match)
- Bosons commute: [cᵢ, cⱼ] = 0, [cᵢ†, cⱼ†] = 0 (no sign change)
- Normal ordering: c† on left, c on right
- **NOT nilpotent**: cᵢ cᵢ ≠ 0 (unlike fermions)

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{BosonicAlgebra}(Int32[-1, 2]);  # c₁† c₂ (already normal ordered)

julia> sm = simplify(m);

julia> length(terms(sm))  # One term
1

julia> monomials(sm)[1].word
2-element Vector{Int32}:
 -1
  2
```

# Note
Unlike monoid-like algebras, bosonic simplification may return multiple terms because
commutation creates delta corrections. There is no `simplify!` variant since the return
type differs from the input type.
"""
function simplify(m::NormalMonomial{BosonicAlgebra,T}) where {T}
    pairs = _simplify_bosonic_word!(copy(m.word))

    coeffs = Int[]
    monos = NormalMonomial{BosonicAlgebra,T}[]
    for (c, w) in pairs
        c == 0 && continue
        push!(coeffs, c)
        push!(monos, NormalMonomial{BosonicAlgebra,T}(w, _OWNED_NORMAL_MONOMIAL))
    end

    if isempty(coeffs)
        return zero(Monomial{BosonicAlgebra,T})
    end

    terms = _combine_physics_terms(coeffs, monos)
    return Monomial(terms)
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
    sort!(word; by=_operator_mode, alg=Base.InsertionSort)

    # Single sweep to find partition boundaries
    sep = [0]
    current_mode = _operator_mode(word[1])
    for i in 2:length(word)
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
    _simplify_bosonic_word!(word::Vector{T}) where {T<:Integer} -> Vector{Tuple{Int,Vector{T}}}

Normal-order a bosonic word using the group-based algorithm with rook numbers.

Returns a vector of (coefficient, normal_ordered_word) pairs representing the sum:
  Σ coeffs[i] * word[i]

This is the low-level function used by `simplify(BosonicAlgebra, word)`.
Integer coefficients are exact (from commutation relations).

# Algorithm
1. Group operators by mode (preserving relative order within each mode)
2. For each mode, compute normal form using rook number formula
3. Expand product across all modes
4. Construct final (coef, word) pairs
"""
function _simplify_bosonic_word!(word::Vector{T}) where {T<:Integer}
    # Handle empty word
    if isempty(word)
        return Tuple{Int,Vector{T}}[(1, T[])]
    end

    # Make a copy since group_by_mode! modifies in place
    word_copy = copy(word)

    # Step 1: Group by mode (stable sort)
    sep = group_by_mode!(word_copy)
    n_modes = length(sep) - 1

    # Step 2: Compute single-mode normal forms
    modes = Vector{Int}(undef, n_modes)
    mode_results = Vector{Tuple{Int,Vector{Tuple{Int,Int,Int}}}}(undef, n_modes)
    for i in 1:n_modes
        modes[i] = _operator_mode(word_copy[sep[i]+1])
        ops = @view word_copy[sep[i]+1:sep[i+1]]
        mode_results[i] = (i, single_mode_normal_form(ops))
    end

    # Step 3: Expand product - but return (Int, Vector{T}) pairs instead of Terms
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

    # Construct final (coef, word) pairs
    result = Tuple{Int,Vector{T}}[]
    for ((cre, ann), coef) in combined
        coef == 0 && continue
        new_word = T[]
        for i in reverse(eachindex(modes))
            m = modes[i]
            append!(new_word, fill(T(-m), cre[i]))
        end
        for (i, m) in enumerate(modes)
            append!(new_word, fill(T(m), ann[i]))
        end
        push!(result, (coef, new_word))
    end

    if isempty(result)
        push!(result, (0, T[]))
    end

    return result
end

# =============================================================================
# Validation (for NormalMonomial constructor)
# =============================================================================

"""
    _validate_bosonic_word!(word::Vector{T}) where {T<:Signed}

Check that a bosonic word is in normal-ordered form. Throws `ArgumentError` if invalid.

Normal order requirements:
- Creation operators (negative indices) come before annihilation operators (positive)
- Creators sorted by mode descending; annihilators sorted by mode ascending

This is used by `NormalMonomial{BosonicAlgebra,T}` constructor to enforce invariants.
"""
function _validate_bosonic_word!(word::Vector{T}) where {T<:Signed}
    isempty(word) && return nothing

    if !is_normal_ordered(word)
        throw(ArgumentError(
            "Bosonic word is not normal-ordered. " *
            "Use simplify(BosonicAlgebra, word) for auto-normal-ordering."
        ))
    end
    return nothing
end

# Connect validation hook used by `NormalMonomial{A,T}` inner constructor.
_validate_word!(::Type{BosonicAlgebra}, word::Vector{T}) where {T<:Signed} =
    _validate_bosonic_word!(word)

# =============================================================================
# Specialized Outer Constructor (validates, rejects non-normal-ordered)
# =============================================================================

"""
    NormalMonomial{BosonicAlgebra}(word::Vector{T}) where {T<:Signed}

Construct a Bosonic monomial, validating that the input is in normal-ordered form.

Throws `ArgumentError` if the word is not normal-ordered. For non-normal-ordered words,
use `simplify(BosonicAlgebra, word)` which auto-normal-orders and returns a `Monomial`
(iterable as `(c_internal, NormalMonomial)` pairs).

Normal order requirements:
- Creation operators (negative indices) come before annihilation operators (positive)
- Creators sorted by mode descending; annihilators sorted by mode ascending

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{BosonicAlgebra}(Int32[-1, 1]);  # c₁† c₁ - normal ordered

julia> m.word
2-element Vector{Int32}:
 -1
  1

julia> NormalMonomial{BosonicAlgebra}(Int32[1, -1])  # c₁ c₁† - NOT normal ordered
ERROR: ArgumentError: Bosonic word is not normal-ordered. Use simplify(BosonicAlgebra, word) for auto-normal-ordering.
```
"""
function NormalMonomial{BosonicAlgebra}(word::Vector{T}) where {T<:Signed}
    word_filtered = filter(!iszero, word)
    return NormalMonomial{BosonicAlgebra,T}(word_filtered, _OWNED_NORMAL_MONOMIAL)
end

function Base.:*(m1::NormalMonomial{BosonicAlgebra,T}, m2::NormalMonomial{BosonicAlgebra,T}) where {T<:Integer}
    return simplify(BosonicAlgebra, vcat(m1.word, m2.word))
end
