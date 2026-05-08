# =============================================================================
# MomentProblem linear-form cache data types
# =============================================================================

"""
    key_lt(a, b) -> Bool

Deterministic ordering for canonical moment keys. Vector keys use lexicographic
ordering by contents; generic keys fall back to `isless`.
"""
key_lt(a, b) = isless(a, b)

function key_lt(a::AbstractVector, b::AbstractVector)
    for (ai, bi) in zip(a, b)
        isequal(ai, bi) || return isless(ai, bi)
    end
    return length(a) < length(b)
end

"""
    key_isequal(a, b) -> Bool

Equality predicate matching `key_lt`'s value semantics. This avoids relying on
mutable array identity for the current `Vector{T}` canonical-key representation.
"""
key_isequal(a, b) = isequal(a, b)

function key_isequal(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return false
    return all(isequal(ai, bi) for (ai, bi) in zip(a, b))
end

# ── linear forms ──────────────────────────────────────────────────────────────

"""
    LinearMomentForm{K,C}

Sorted, deduplicated, zero-pruned linear combination of canonical moments.
`terms[i] = (key, coefficient)` with keys in ascending `key_lt` order.
"""
struct LinearMomentForm{K,C}
    terms::Vector{Pair{K,C}}

    function LinearMomentForm{K,C}(pairs) where {K,C}
        terms = Pair{K,C}[]
        for (key, coef) in pairs
            converted = convert(C, coef)
            iszero(converted) && continue
            push!(terms, convert(K, key) => converted)
        end

        sort!(terms; by=x -> x.first, lt=key_lt)

        out = Pair{K,C}[]
        for (key, coef) in terms
            if !isempty(out) && key_isequal(out[end].first, key)
                combined = out[end].second + coef
                if iszero(combined)
                    pop!(out)
                else
                    out[end] = out[end].first => combined
                end
            else
                push!(out, key => coef)
            end
        end

        return new{K,C}(out)
    end
end

LinearMomentForm(pairs::AbstractVector{<:Pair{K,C}}) where {K,C} =
    LinearMomentForm{K,C}(pairs)

Base.length(form::LinearMomentForm) = length(form.terms)
Base.isempty(form::LinearMomentForm) = isempty(form.terms)
Base.iterate(form::LinearMomentForm, state...) = iterate(form.terms, state...)

function _assert_linear_moment_form_invariants(form::LinearMomentForm)
    prev = nothing
    for (idx, (key, coef)) in enumerate(form.terms)
        iszero(coef) && throw(ArgumentError("LinearMomentForm contains a zero coefficient at term $idx"))
        if prev !== nothing
            key_lt(prev, key) || throw(ArgumentError(
                "LinearMomentForm terms must be strictly sorted and deduplicated; violation at term $idx"
            ))
        end
        prev = key
    end
    return nothing
end

# ── pivots ───────────────────────────────────────────────────────────────────

"""
    Pivot{C}

Records that a physical canonical moment is represented by one PSD/HPSD block
entry with a unit phase. `adjoint` disambiguates shared Hermitian upper-triangle
positions for a key and its adjoint.
"""
struct Pivot{C}
    block::Int
    row::Int
    col::Int
    phase::C
    adjoint::Bool

    function Pivot{C}(block::Integer, row::Integer, col::Integer, phase, adjoint::Bool=false) where {C}
        _assert_positive_index("pivot block", block)
        _assert_positive_index("pivot row", row)
        _assert_positive_index("pivot col", col)
        return new{C}(Int(block), Int(row), Int(col), convert(C, phase), adjoint)
    end
end

Pivot(block::Integer, row::Integer, col::Integer, phase, adjoint::Bool=false) =
    Pivot{typeof(phase)}(block, row, col, phase, adjoint)

# ── block provenance ─────────────────────────────────────────────────────────

abstract type BlockOrigin end

struct MomentMatrixOrigin <: BlockOrigin
    clique::Int
    ts_block::Int

    function MomentMatrixOrigin(clique::Integer, ts_block::Integer)
        _assert_positive_index("moment-matrix clique", clique)
        _assert_positive_index("moment-matrix term-sparsity block", ts_block)
        return new(Int(clique), Int(ts_block))
    end
end

struct LocalizingOrigin <: BlockOrigin
    clique::Int
    cons_idx::Int
    ts_block::Int

    function LocalizingOrigin(clique::Integer, cons_idx::Integer, ts_block::Integer)
        _assert_positive_index("localizing clique", clique)
        _assert_positive_index("localizing constraint index", cons_idx)
        _assert_positive_index("localizing term-sparsity block", ts_block)
        return new(Int(clique), Int(cons_idx), Int(ts_block))
    end
end

struct GlobalOrigin <: BlockOrigin
    cons_idx::Int

    function GlobalOrigin(cons_idx::Integer)
        _assert_positive_index("global constraint index", cons_idx)
        return new(Int(cons_idx))
    end
end

struct AuxOrigin{K} <: BlockOrigin
    reason::Symbol
    key_group::Vector{K}

    function AuxOrigin{K}(reason::Symbol, key_group::AbstractVector{K}) where {K}
        reason == :orphan_packing || throw(ArgumentError(
            "Unsupported auxiliary block reason $(repr(reason)); expected :orphan_packing"
        ))
        return new{K}(reason, collect(key_group))
    end
end

AuxOrigin(reason::Symbol, key_group::AbstractVector{K}) where {K} = AuxOrigin{K}(reason, key_group)

struct BlockMeta{M}
    cone::Symbol
    origin::BlockOrigin
    row_labels::Vector{M}

    function BlockMeta{M}(cone::Symbol, origin::BlockOrigin, row_labels::AbstractVector{M}) where {M}
        _assert_psd_block_cone(cone)
        return new{M}(cone, origin, collect(row_labels))
    end
end

BlockMeta(cone::Symbol, origin::BlockOrigin, row_labels::AbstractVector{M}) where {M} =
    BlockMeta{M}(cone, origin, row_labels)

struct PSDBlockLin{K,C,M}
    size::Int
    entries::Matrix{LinearMomentForm{K,C}}
    meta::BlockMeta{M}

    function PSDBlockLin{K,C,M}(
        block_size::Integer,
        entries::Matrix{LinearMomentForm{K,C}},
        meta::BlockMeta{M},
    ) where {K,C,M}
        block_size >= 0 || throw(ArgumentError("PSD block size must be nonnegative, got $block_size"))
        n = Int(block_size)
        Base.size(entries) == (n, n) || throw(DimensionMismatch(
            "PSD linear block entries must have size ($n, $n), got $(Base.size(entries))"
        ))
        length(meta.row_labels) == n || throw(DimensionMismatch(
            "PSD block row label count $(length(meta.row_labels)) does not match block size $n"
        ))
        for form in entries
            _assert_linear_moment_form_invariants(form)
        end
        return new{K,C,M}(n, entries, meta)
    end
end

PSDBlockLin(entries::Matrix{LinearMomentForm{K,C}}, meta::BlockMeta{M}) where {K,C,M} =
    PSDBlockLin{K,C,M}(Base.size(entries, 1), entries, meta)

# ── scalar constraints ───────────────────────────────────────────────────────

abstract type ConstraintOrigin end

struct ZeroMatrixOrigin <: ConstraintOrigin
    constraint_idx::Int
    row::Int
    col::Int
    part::Symbol

    function ZeroMatrixOrigin(constraint_idx::Integer, row::Integer, col::Integer, part::Symbol)
        _assert_positive_index("zero constraint index", constraint_idx)
        _assert_positive_index("zero constraint row", row)
        _assert_positive_index("zero constraint col", col)
        part in (:real, :imag, :scalar) || throw(ArgumentError(
            "Unsupported zero-constraint part $(repr(part)); expected :real, :imag, or :scalar"
        ))
        return new(Int(constraint_idx), Int(row), Int(col), part)
    end
end

struct ParityOrigin <: ConstraintOrigin
    cons_idx::Int

    function ParityOrigin(cons_idx::Integer)
        _assert_positive_index("parity constraint index", cons_idx)
        return new(Int(cons_idx))
    end
end

struct MomentEqOrigin <: ConstraintOrigin
    cons_idx::Int
    row::Int

    function MomentEqOrigin(cons_idx::Integer, row::Integer)
        _assert_positive_index("moment equality constraint index", cons_idx)
        _assert_positive_index("moment equality row", row)
        return new(Int(cons_idx), Int(row))
    end
end

struct IdentityOrigin <: ConstraintOrigin end

struct ScalarLinearConstraint{K,C}
    form::LinearMomentForm{K,C}
    kind::Symbol
    origin::ConstraintOrigin

    function ScalarLinearConstraint{K,C}(
        form::LinearMomentForm{K,C},
        kind::Symbol,
        origin::ConstraintOrigin,
    ) where {K,C}
        kind in (:zero, :identity_norm) || throw(ArgumentError(
            "Unsupported scalar linear constraint kind $(repr(kind)); expected :zero or :identity_norm"
        ))
        _assert_linear_moment_form_invariants(form)
        return new{K,C}(form, kind, origin)
    end
end

ScalarLinearConstraint(form::LinearMomentForm{K,C}, kind::Symbol, origin::ConstraintOrigin) where {K,C} =
    ScalarLinearConstraint{K,C}(form, kind, origin)

# ── linear data cache ────────────────────────────────────────────────────────

"""
    MomentLinearData{K,C,M}

Derived linear-form view of a `MomentProblem`. Populated by `moment_relax`
after all symbolic constraint mutations have completed.
"""
struct MomentLinearData{K,C,M}
    moments::Vector{K}
    moment_index::Dict{K,Int}
    identity::K
    key_to_monomial::Dict{K,M}
    adjoint_key::Dict{K,K}

    psd_blocks_lin::Vector{PSDBlockLin{K,C,M}}
    psd_block_constraint_idx::Vector{Int}
    zero_constraints::Vector{ScalarLinearConstraint{K,C}}
    objective_lin::LinearMomentForm{K,C}

    pivots::Dict{K,Pivot{C}}
    pivot_at::Dict{Tuple{Int,Int,Int},Vector{K}}
    free_keys::Vector{K}

    function MomentLinearData{K,C,M}(
        moments::Vector{K},
        moment_index::Dict{K,Int},
        identity::K,
        key_to_monomial::Dict{K,M},
        adjoint_key::Dict{K,K},
        psd_blocks_lin::Vector{PSDBlockLin{K,C,M}},
        psd_block_constraint_idx::Vector{Int},
        zero_constraints::Vector{ScalarLinearConstraint{K,C}},
        objective_lin::LinearMomentForm{K,C},
        pivots::Dict{K,Pivot{C}},
        pivot_at::Dict{Tuple{Int,Int,Int},Vector{K}},
        free_keys::Vector{K},
    ) where {K,C,M}
        _assert_moment_linear_data_invariants(
            moments,
            moment_index,
            identity,
            key_to_monomial,
            adjoint_key,
            psd_blocks_lin,
            psd_block_constraint_idx,
            zero_constraints,
            objective_lin,
            pivots,
            pivot_at,
            free_keys,
        )
        return new{K,C,M}(
            moments,
            moment_index,
            identity,
            key_to_monomial,
            adjoint_key,
            psd_blocks_lin,
            psd_block_constraint_idx,
            zero_constraints,
            objective_lin,
            pivots,
            pivot_at,
            free_keys,
        )
    end
end

function assert_moment_linear_data_invariants(linear::MomentLinearData)
    return _assert_moment_linear_data_invariants(
        linear.moments,
        linear.moment_index,
        linear.identity,
        linear.key_to_monomial,
        linear.adjoint_key,
        linear.psd_blocks_lin,
        linear.psd_block_constraint_idx,
        linear.zero_constraints,
        linear.objective_lin,
        linear.pivots,
        linear.pivot_at,
        linear.free_keys,
    )
end

function assert_moment_linear_data_invariants(linear::MomentLinearData, constraints::AbstractVector)
    assert_moment_linear_data_invariants(linear)
    return _assert_moment_linear_data_constraint_invariants(linear, constraints)
end

# =============================================================================
# Invariant helpers
# =============================================================================

function _assert_positive_index(name::AbstractString, value::Integer)
    value >= 1 || throw(ArgumentError("$name must be positive, got $value"))
    return nothing
end

function _assert_psd_block_cone(cone::Symbol)
    cone in (:PSD, :HPSD, :AuxHPSD) || throw(ArgumentError(
        "Unsupported PSD linear block cone $(repr(cone)); expected :PSD, :HPSD, or :AuxHPSD"
    ))
    return nothing
end

function _keys_match(xs::AbstractVector, ys::AbstractVector)
    length(xs) == length(ys) || return false
    sx = sort!(collect(xs); lt=key_lt)
    sy = sort!(collect(ys); lt=key_lt)
    return all(key_isequal(x, y) for (x, y) in zip(sx, sy))
end

function _has_key_value(keys, key)
    return any(candidate -> key_isequal(candidate, key), keys)
end

function _has_key_value(keys::AbstractSet, key)
    return key in keys
end

function _get_key_value(dict::Dict{K,V}, key::K, what::AbstractString) where {K,V}
    haskey(dict, key) && return dict[key]
    throw(ArgumentError("Missing $what for canonical key $(repr(key))"))
end

function _assert_form_keys_in_moment_universe(form::LinearMomentForm{K}, moment_keys) where {K}
    _assert_linear_moment_form_invariants(form)
    for (key, _) in form
        _has_key_value(moment_keys, key) || throw(ArgumentError(
            "Linear form references key outside MomentLinearData.moments: $(repr(key))"
        ))
    end
    return nothing
end

function _assert_moment_index(moments::Vector{K}, moment_index::Dict{K,Int}) where {K}
    length(moment_index) == length(moments) || throw(ArgumentError(
        "moment_index has $(length(moment_index)) entries for $(length(moments)) moments"
    ))
    for (idx, key) in enumerate(moments)
        stored_idx = _get_key_value(moment_index, key, "moment index")
        stored_idx == idx || throw(ArgumentError(
            "moment_index maps key $(repr(key)) to $stored_idx, expected $idx"
        ))
    end
    return nothing
end

function _assert_key_to_monomial(moments::Vector{K}, key_to_monomial::Dict{K,M}) where {K,M}
    _keys_match(moments, collect(keys(key_to_monomial))) || throw(ArgumentError(
        "MomentLinearData.moments must match keys(key_to_monomial)"
    ))

    for (key, mono) in key_to_monomial
        derived = symmetric_canon(expval(mono))
        key_isequal(derived, key) || throw(ArgumentError(
            "key_to_monomial representative does not canonicalize back to its key: " *
            "got $(repr(derived)), expected $(repr(key))"
        ))
    end
    return nothing
end

function _assert_identity_key(::Type{M}, identity) where {M}
    derived = symmetric_canon(expval(one(M)))
    key_isequal(identity, derived) || throw(ArgumentError(
        "MomentLinearData.identity is $(repr(identity)); expected $(repr(derived))"
    ))
    return nothing
end

function _assert_adjoint_keys(moments::Vector{K}, moment_set, adjoint_key::Dict{K,K}, key_to_monomial::Dict{K,M}) where {K,M}
    length(adjoint_key) == length(moments) || throw(ArgumentError(
        "adjoint_key has $(length(adjoint_key)) entries for $(length(moments)) moments"
    ))

    is_complex = _moment_linear_is_complex_monomial(M)
    for key in moments
        adj = _get_key_value(adjoint_key, key, "adjoint key")
        _has_key_value(moment_set, adj) || throw(ArgumentError(
            "adjoint_key maps $(repr(key)) outside the moment universe to $(repr(adj))"
        ))

        if is_complex
            adj_adj = _get_key_value(adjoint_key, adj, "double adjoint key")
            key_isequal(adj_adj, key) || throw(ArgumentError(
                "adjoint_key is not involutive at $(repr(key)): got $(repr(adj_adj))"
            ))
        else
            key_isequal(adj, key) || throw(ArgumentError(
                "real-algebra adjoint_key must be the identity map at $(repr(key))"
            ))
        end
    end
    return nothing
end

_moment_linear_is_complex_monomial(::Type{<:NormalMonomial{A,T}}) where {A<:AlgebraType,T<:Integer} =
    _is_complex_problem(A)
_moment_linear_is_complex_monomial(::Type) = true

function _assert_psd_blocks(
    psd_blocks_lin::Vector{PSDBlockLin{K,C,M}},
    psd_block_constraint_idx::Vector{Int},
    moments,
) where {K,C,M}
    length(psd_blocks_lin) == length(psd_block_constraint_idx) || throw(ArgumentError(
        "psd_blocks_lin has $(length(psd_blocks_lin)) blocks but psd_block_constraint_idx has $(length(psd_block_constraint_idx)) entries"
    ))
    for constraint_idx in psd_block_constraint_idx
        _assert_positive_index("PSD block constraint index", constraint_idx)
    end
    for block in psd_blocks_lin
        block.size == Base.size(block.entries, 1) == Base.size(block.entries, 2) || throw(DimensionMismatch(
            "PSD linear block size $(block.size) does not match entries size $(Base.size(block.entries))"
        ))
        length(block.meta.row_labels) == block.size || throw(DimensionMismatch(
            "PSD linear block row label count $(length(block.meta.row_labels)) does not match size $(block.size)"
        ))
        for form in block.entries
            _assert_form_keys_in_moment_universe(form, moments)
        end
    end
    return nothing
end

function _form_coefficient(form::LinearMomentForm{K,C}, key::K) where {K,C}
    for (candidate, coef) in form
        key_isequal(candidate, key) && return coef
    end
    return zero(C)
end

function _assert_self_adjoint_form(form::LinearMomentForm{K,C}, adjoint_key::Dict{K,K}) where {K,C}
    C <: Complex || return nothing
    rtol = sqrt(eps(float(real(one(C)))))
    for (key, coef) in form
        adj = _get_key_value(adjoint_key, key, "adjoint key")
        adj_coef = _form_coefficient(form, adj)
        isapprox(adj_coef, conj(coef); atol=0, rtol=rtol) || throw(ArgumentError(
            "zero constraint form is not self-adjoint at key $(repr(key)): " *
            "coefficient $(repr(coef)), adjoint coefficient $(repr(adj_coef))"
        ))
    end
    return nothing
end

function _assert_zero_constraints(
    zero_constraints::Vector{ScalarLinearConstraint{K,C}},
    moments,
    adjoint_key::Dict{K,K},
) where {K,C}
    for zc in zero_constraints
        zc.kind == :zero || throw(ArgumentError(
            "zero_constraints contains constraint of kind $(repr(zc.kind)); expected :zero"
        ))
        _assert_form_keys_in_moment_universe(zc.form, moments)
        _assert_self_adjoint_form(zc.form, adjoint_key)
    end
    return nothing
end

function _phase_is_unit(phase)
    r = abs(phase)
    return isapprox(r, one(r); atol=0, rtol=sqrt(eps(float(one(r)))))
end

function _assert_pivots(
    moments::Vector{K},
    moment_set,
    pivots::Dict{K,Pivot{C}},
    pivot_at::Dict{Tuple{Int,Int,Int},Vector{K}},
    free_keys::Vector{K},
    adjoint_key::Dict{K,K},
    psd_blocks_lin::Vector{PSDBlockLin{K,C,M}},
) where {K,C,M}
    pivot_keys = collect(keys(pivots))
    union_keys = vcat(pivot_keys, free_keys)
    _keys_match(union_keys, moments) || throw(ArgumentError(
        "pivots and free_keys must cover exactly MomentLinearData.moments"
    ))

    free_key_set = Set(free_keys)
    for key in pivot_keys
        _has_key_value(free_key_set, key) && throw(ArgumentError(
            "key $(repr(key)) appears both in pivots and free_keys"
        ))
    end

    expected_pivot_at = Dict{Tuple{Int,Int,Int},Vector{K}}()
    is_complex = _moment_linear_is_complex_monomial(M)
    for (key, pivot) in pivots
        _has_key_value(moment_set, key) || throw(ArgumentError(
            "pivot key outside moment universe: $(repr(key))"
        ))
        1 <= pivot.block <= length(psd_blocks_lin) || throw(ArgumentError(
            "pivot for key $(repr(key)) references block $(pivot.block), but there are $(length(psd_blocks_lin)) PSD blocks"
        ))
        block = psd_blocks_lin[pivot.block]
        1 <= pivot.row <= block.size || throw(ArgumentError("pivot row $(pivot.row) outside block size $(block.size)"))
        1 <= pivot.col <= block.size || throw(ArgumentError("pivot col $(pivot.col) outside block size $(block.size)"))
        _phase_is_unit(pivot.phase) || throw(ArgumentError(
            "pivot phase for key $(repr(key)) is not unit: $(repr(pivot.phase))"
        ))
        if !is_complex
            pivot.phase == one(C) || throw(ArgumentError(
                "real-algebra pivot phase must be 1, got $(repr(pivot.phase))"
            ))
        end
        push!(get!(expected_pivot_at, (pivot.block, pivot.row, pivot.col), K[]), key)
    end

    _keys_match(collect(keys(expected_pivot_at)), collect(keys(pivot_at))) || throw(ArgumentError(
        "pivot_at keys do not match pivot positions"
    ))

    for (position, keys_at_position) in pivot_at
        expected_keys = expected_pivot_at[position]
        _keys_match(keys_at_position, expected_keys) || throw(ArgumentError(
            "pivot_at[$position] does not list exactly the keys pivoting at that entry"
        ))
        length(keys_at_position) <= 2 || throw(ArgumentError(
            "pivot_at[$position] has $(length(keys_at_position)) keys; expected at most 2"
        ))
        if length(keys_at_position) == 2
            k1, k2 = keys_at_position
            adj = _get_key_value(adjoint_key, k1, "adjoint key")
            key_isequal(adj, k2) || throw(ArgumentError(
                "two keys sharing pivot position $position must be adjoints"
            ))
        end
    end

    return nothing
end

function _assert_moment_linear_data_invariants(
    moments::Vector{K},
    moment_index::Dict{K,Int},
    identity::K,
    key_to_monomial::Dict{K,M},
    adjoint_key::Dict{K,K},
    psd_blocks_lin::Vector{PSDBlockLin{K,C,M}},
    psd_block_constraint_idx::Vector{Int},
    zero_constraints::Vector{ScalarLinearConstraint{K,C}},
    objective_lin::LinearMomentForm{K,C},
    pivots::Dict{K,Pivot{C}},
    pivot_at::Dict{Tuple{Int,Int,Int},Vector{K}},
    free_keys::Vector{K},
) where {K,C,M}
    issorted(moments; lt=key_lt) || throw(ArgumentError("MomentLinearData.moments must be sorted by key_lt"))
    for idx in 2:length(moments)
        key_isequal(moments[idx - 1], moments[idx]) && throw(ArgumentError(
            "MomentLinearData.moments contains duplicate key $(repr(moments[idx]))"
        ))
    end

    moment_set = Set(moments)

    _assert_moment_index(moments, moment_index)
    _assert_key_to_monomial(moments, key_to_monomial)
    _assert_identity_key(M, identity)
    _assert_adjoint_keys(moments, moment_set, adjoint_key, key_to_monomial)
    _assert_form_keys_in_moment_universe(objective_lin, moment_set)
    _assert_psd_blocks(psd_blocks_lin, psd_block_constraint_idx, moment_set)
    _assert_zero_constraints(zero_constraints, moment_set, adjoint_key)
    _assert_pivots(moments, moment_set, pivots, pivot_at, free_keys, adjoint_key, psd_blocks_lin)
    return nothing
end

function _assert_moment_linear_data_constraint_invariants(linear::MomentLinearData, constraints::AbstractVector)
    for (block_idx, constraint_idx) in enumerate(linear.psd_block_constraint_idx)
        1 <= constraint_idx <= length(constraints) || throw(ArgumentError(
            "PSD block $block_idx points to constraint index $constraint_idx, but there are $(length(constraints)) constraints"
        ))
        cone, mat = constraints[constraint_idx]
        cone == linear.psd_blocks_lin[block_idx].meta.cone || throw(ArgumentError(
            "PSD block $block_idx cone $(repr(linear.psd_blocks_lin[block_idx].meta.cone)) does not match constraint cone $(repr(cone))"
        ))
        linear.psd_blocks_lin[block_idx].size == Base.size(mat, 1) || throw(DimensionMismatch(
            "PSD block $block_idx size $(linear.psd_blocks_lin[block_idx].size) does not match constraint matrix size $(Base.size(mat))"
        ))
    end
    return nothing
end
