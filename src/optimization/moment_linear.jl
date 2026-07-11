# =============================================================================
# MomentProblem linear-form cache data types
# =============================================================================

"""
    key_lt(a, b) -> Bool

Deterministic ordering for canonical moment keys. Vector keys use lexicographic
ordering by contents; generic keys fall back to `isless`.
"""
key_lt(a, b) = isless(a, b)

function key_lt(a::Vector, b::Vector)
    la = length(a)
    lb = length(b)
    @inbounds for idx in 1:min(la, lb)
        ai = a[idx]
        bi = b[idx]
        isequal(ai, bi) || return isless(ai, bi)
    end
    return la < lb
end

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

function key_isequal(a::Vector, b::Vector)
    length(a) == length(b) || return false
    @inbounds for idx in eachindex(a, b)
        isequal(a[idx], b[idx]) || return false
    end
    return true
end

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

        return _linear_moment_form_from_owned_pairs!(terms)
    end

    function LinearMomentForm{K,C}(terms::Vector{Pair{K,C}}, ::Val{:trusted}) where {K,C}
        return new{K,C}(terms)
    end
end

function _linear_moment_form_from_owned_pairs!(terms::Vector{Pair{K,C}}) where {K,C}
    isempty(terms) && return LinearMomentForm{K,C}(terms, Val(:trusted))
    if length(terms) == 1
        if iszero(terms[1].second)
            resize!(terms, 0)
        end
        return LinearMomentForm{K,C}(terms, Val(:trusted))
    elseif length(terms) == 2
        first_term = terms[1]
        second_term = terms[2]
        if key_isequal(first_term.first, second_term.first)
            coef = first_term.second + second_term.second
            if iszero(coef)
                resize!(terms, 0)
            else
                terms[1] = first_term.first => coef
                resize!(terms, 1)
            end
            return LinearMomentForm{K,C}(terms, Val(:trusted))
        end

        write_idx = 0
        if !iszero(first_term.second)
            write_idx += 1
            terms[write_idx] = first_term
        end
        if !iszero(second_term.second)
            write_idx += 1
            terms[write_idx] = second_term
        end
        resize!(terms, write_idx)
        if write_idx == 2 && key_lt(terms[2].first, terms[1].first)
            terms[1], terms[2] = terms[2], terms[1]
        end
        return LinearMomentForm{K,C}(terms, Val(:trusted))
    end

    sort!(terms; by=x -> x.first, lt=key_lt)

    write_idx = 0
    current_key = terms[1].first
    current_coef = terms[1].second

    @inbounds for i in 2:length(terms)
        key = terms[i].first
        coef = terms[i].second
        if key_isequal(current_key, key)
            current_coef += coef
        else
            if !iszero(current_coef)
                write_idx += 1
                terms[write_idx] = current_key => current_coef
            end
            current_key = key
            current_coef = coef
        end
    end

    if !iszero(current_coef)
        write_idx += 1
        terms[write_idx] = current_key => current_coef
    end

    resize!(terms, write_idx)
    return LinearMomentForm{K,C}(terms, Val(:trusted))
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
    trusted_self_adjoint::Bool

    function ScalarLinearConstraint{K,C}(
        form::LinearMomentForm{K,C},
        kind::Symbol,
        origin::ConstraintOrigin,
        trusted_self_adjoint::Bool=false,
    ) where {K,C}
        kind in (:zero, :identity_norm) || throw(ArgumentError(
            "Unsupported scalar linear constraint kind $(repr(kind)); expected :zero or :identity_norm"
        ))
        _assert_linear_moment_form_invariants(form)
        return new{K,C}(form, kind, origin, trusted_self_adjoint)
    end
end

ScalarLinearConstraint(
    form::LinearMomentForm{K,C},
    kind::Symbol,
    origin::ConstraintOrigin;
    trusted_self_adjoint::Bool=false,
) where {K,C} =
    ScalarLinearConstraint{K,C}(form, kind, origin, trusted_self_adjoint)

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
        ;
        stage_times_ns=nothing,
        stage_prefix::Symbol=:moment_linear_data,
        zero_constraint_keys_registered::Bool=false,
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
            stage_times_ns=stage_times_ns,
            stage_prefix=stage_prefix,
            zero_constraint_keys_registered=zero_constraint_keys_registered,
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
    terms = form.terms
    lo = firstindex(terms)
    hi = lastindex(terms)
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        candidate = terms[mid].first
        if key_isequal(candidate, key)
            return terms[mid].second
        elseif key_lt(candidate, key)
            lo = mid + 1
        else
            hi = mid - 1
        end
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
    ;
    stage_times_ns=nothing,
    stage_prefix::Symbol=:moment_linear_data_zero_constraints,
    keys_registered::Bool=false,
) where {K,C}
    if stage_times_ns === nothing
        for zc in zero_constraints
            zc.kind == :zero || throw(ArgumentError(
                "zero_constraints contains constraint of kind $(repr(zc.kind)); expected :zero"
            ))
            if !keys_registered
                _assert_form_keys_in_moment_universe(zc.form, moments)
            end
            zc.trusted_self_adjoint || _assert_self_adjoint_form(zc.form, adjoint_key)
        end
        return nothing
    end

    for zc in zero_constraints
        stage_start_ns = time_ns()
        zc.kind == :zero || throw(ArgumentError(
            "zero_constraints contains constraint of kind $(repr(zc.kind)); expected :zero"
        ))
        stage_key = Symbol(stage_prefix, "_kind")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)

        stage_start_ns = time_ns()
        if keys_registered
            stage_key = Symbol(stage_prefix, "_trusted_keys")
        else
            _assert_form_keys_in_moment_universe(zc.form, moments)
            stage_key = Symbol(stage_prefix, "_moment_keys")
        end
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)

        stage_start_ns = time_ns()
        if zc.trusted_self_adjoint
            stage_key = Symbol(stage_prefix, "_trusted_self_adjoint")
        else
            _assert_self_adjoint_form(zc.form, adjoint_key)
            stage_key = Symbol(stage_prefix, "_self_adjoint")
        end
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
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
    free_key_set = Set(free_keys)
    length(pivots) + length(free_keys) == length(moments) || throw(ArgumentError(
        "pivots and free_keys must cover exactly MomentLinearData.moments"
    ))
    for key in moments
        (haskey(pivots, key) || _has_key_value(free_key_set, key)) || throw(ArgumentError(
            "pivots and free_keys must cover moment key $(repr(key))"
        ))
    end

    for key in keys(pivots)
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

    length(expected_pivot_at) == length(pivot_at) || throw(ArgumentError(
        "pivot_at keys do not match pivot positions"
    ))
    for position in keys(expected_pivot_at)
        haskey(pivot_at, position) || throw(ArgumentError(
            "pivot_at keys do not match pivot positions"
        ))
    end

    for (position, keys_at_position) in pivot_at
        expected_keys = expected_pivot_at[position]
        (length(keys_at_position) == length(expected_keys) &&
            all(key -> _has_key_value(expected_keys, key), keys_at_position)) || throw(ArgumentError(
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
    ;
    stage_times_ns=nothing,
    stage_prefix::Symbol=:moment_linear_data,
    zero_constraint_keys_registered::Bool=false,
) where {K,C,M}
    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    issorted(moments; lt=key_lt) || throw(ArgumentError("MomentLinearData.moments must be sorted by key_lt"))
    for idx in 2:length(moments)
        key_isequal(moments[idx - 1], moments[idx]) && throw(ArgumentError(
            "MomentLinearData.moments contains duplicate key $(repr(moments[idx]))"
        ))
    end
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_moments_sorted")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    moment_set = Set(moments)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_moment_set")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_moment_index(moments, moment_index)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_moment_index")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_key_to_monomial(moments, key_to_monomial)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_key_to_monomial")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_identity_key(M, identity)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_identity")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_adjoint_keys(moments, moment_set, adjoint_key, key_to_monomial)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_adjoint_keys")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_form_keys_in_moment_universe(objective_lin, moment_set)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_objective")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_psd_blocks(psd_blocks_lin, psd_block_constraint_idx, moment_set)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_psd_blocks")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_zero_constraints(
        zero_constraints,
        moment_set,
        adjoint_key;
        stage_times_ns=stage_times_ns,
        stage_prefix=Symbol(stage_prefix, "_zero_constraints"),
        keys_registered=zero_constraint_keys_registered,
    )
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_zero_constraints")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = stage_times_ns === nothing ? 0 : time_ns()
    _assert_pivots(moments, moment_set, pivots, pivot_at, free_keys, adjoint_key, psd_blocks_lin)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_pivots")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end
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

# ── staged construction ──────────────────────────────────────────────────────

"""
    MomentLinearBuilder(K, C, M)

Staging buffer for constructing `MomentLinearData` without first materializing
full polynomial matrices. Inputs are copied on insertion so caller-owned scratch
keys can be reused or mutated before finalization.
"""
mutable struct MomentLinearBuilder{K,C,M<:NormalMonomial}
    identity::K
    key_to_monomial::Dict{K,M}
    objective_terms::Vector{Pair{K,C}}
    psd_blocks_lin::Vector{PSDBlockLin{K,C,M}}
    psd_block_constraint_idx::Vector{Int}
    zero_constraints::Vector{ScalarLinearConstraint{K,C}}
    zero_constraint_keys_registered::Bool
    finalized::Bool
end

function MomentLinearBuilder(::Type{K}, ::Type{C}, ::Type{M}) where {K,C,M<:NormalMonomial}
    identity_mono = one(M)
    identity = _owned_moment_key(K, symmetric_canon(expval(identity_mono)))
    key_to_monomial = Dict{K,M}(identity => identity_mono)
    return MomentLinearBuilder{K,C,M}(
        identity,
        key_to_monomial,
        Pair{K,C}[],
        PSDBlockLin{K,C,M}[],
        Int[],
        ScalarLinearConstraint{K,C}[],
        true,
        false,
    )
end

function _owned_moment_key(::Type{K}, key) where {K}
    converted = convert(K, key)
    return converted isa AbstractVector ? copy(converted) : converted
end

function _check_not_finalized!(builder::MomentLinearBuilder)
    builder.finalized && throw(ArgumentError("MomentLinearBuilder has already been finalized"))
    return nothing
end

function _builder_stored_key(key_to_monomial::Dict{K,M}, key::K) where {K,M}
    stored = getkey(key_to_monomial, key, nothing)
    stored === nothing || return stored

    # Fallback for mutable key types whose stored hash no longer matches their
    # current value. Canonical moment keys should not mutate, so this is cold.
    for stored_key in keys(key_to_monomial)
        key_isequal(stored_key, key) && return stored_key
    end
    return nothing
end

function _register_builder_moment_key!(
    key_to_monomial::Dict{K,M},
    ::Type{K},
    key,
    mono::M,
) where {K,M}
    owned = _owned_moment_key(K, key)
    stored = _builder_stored_key(key_to_monomial, owned)
    stored === nothing || return stored
    key_to_monomial[owned] = mono
    return owned
end

function register_moment!(builder::MomentLinearBuilder{K,C,M}, mono::M) where {K,C,M}
    _check_not_finalized!(builder)
    key = _register_builder_moment_key!(
        builder.key_to_monomial,
        K,
        symmetric_canon(expval(mono)),
        mono,
    )
    return _owned_moment_key(K, key)
end

function register_moment!(builder::MomentLinearBuilder{K,C,M}, key, mono::M) where {K,C,M}
    _check_not_finalized!(builder)
    stored = _register_builder_moment_key!(builder.key_to_monomial, K, key, mono)
    return _owned_moment_key(K, stored)
end

function _owned_builder_pairs(::Type{K}, ::Type{C}, pairs) where {K,C}
    owned = Pair{K,C}[]
    for (key, coef) in pairs
        converted = convert(C, coef)
        iszero(converted) && continue
        push!(owned, _owned_moment_key(K, key) => converted)
    end
    return owned
end

function _builder_linear_form(::Type{K}, ::Type{C}, pairs) where {K,C}
    return _linear_moment_form_from_owned_pairs!(_owned_builder_pairs(K, C, pairs))
end

function _owned_builder_form(::Type{K}, ::Type{C}, form::LinearMomentForm{K,C}) where {K,C}
    terms = Pair{K,C}[]
    sizehint!(terms, length(form))
    for (key, coef) in form
        push!(terms, _owned_moment_key(K, key) => coef)
    end
    return LinearMomentForm{K,C}(terms, Val(:trusted))
end

function add_objective_terms!(builder::MomentLinearBuilder{K,C}, pairs) where {K,C}
    _check_not_finalized!(builder)
    append!(builder.objective_terms, _owned_builder_pairs(K, C, pairs))
    return builder
end

function add_psd_block!(
    builder::MomentLinearBuilder{K,C,M},
    cone::Symbol,
    entry_pairs::AbstractMatrix,
    meta::BlockMeta{M};
    constraint_idx::Integer,
) where {K,C,M}
    _check_not_finalized!(builder)
    _assert_positive_index("PSD block constraint index", constraint_idx)
    _assert_psd_block_cone(cone)
    cone == meta.cone || throw(ArgumentError(
        "PSD block cone $(repr(cone)) does not match metadata cone $(repr(meta.cone))"
    ))

    nrows, ncols = size(entry_pairs)
    nrows == ncols || throw(DimensionMismatch("$cone staged block must be square, got $(size(entry_pairs))"))

    entries = Matrix{LinearMomentForm{K,C}}(undef, nrows, ncols)
    for idx in eachindex(entry_pairs)
        entries[idx] = _builder_linear_form(K, C, entry_pairs[idx])
    end

    push!(builder.psd_blocks_lin, PSDBlockLin{K,C,M}(nrows, entries, meta))
    push!(builder.psd_block_constraint_idx, Int(constraint_idx))
    return builder
end

function add_psd_block!(
    builder::MomentLinearBuilder{K,C,M},
    cone::Symbol,
    entry_forms::AbstractMatrix{<:LinearMomentForm{K,C}},
    meta::BlockMeta{M};
    constraint_idx::Integer,
) where {K,C,M}
    _check_not_finalized!(builder)
    _assert_positive_index("PSD block constraint index", constraint_idx)
    _assert_psd_block_cone(cone)
    cone == meta.cone || throw(ArgumentError(
        "PSD block cone $(repr(cone)) does not match metadata cone $(repr(meta.cone))"
    ))

    nrows, ncols = size(entry_forms)
    nrows == ncols || throw(DimensionMismatch("$cone staged block must be square, got $(size(entry_forms))"))

    entries = Matrix{LinearMomentForm{K,C}}(undef, nrows, ncols)
    for idx in eachindex(entry_forms)
        entries[idx] = _owned_builder_form(K, C, entry_forms[idx])
    end

    push!(builder.psd_blocks_lin, PSDBlockLin{K,C,M}(nrows, entries, meta))
    push!(builder.psd_block_constraint_idx, Int(constraint_idx))
    return builder
end

function add_zero_constraint!(
    builder::MomentLinearBuilder{K,C},
    pairs,
    origin::ConstraintOrigin;
    kind::Symbol=:zero,
) where {K,C}
    _check_not_finalized!(builder)
    form = _builder_linear_form(K, C, pairs)
    push!(builder.zero_constraints, ScalarLinearConstraint(form, kind, origin))
    builder.zero_constraint_keys_registered = false
    return builder
end

function add_zero_constraint!(
    builder::MomentLinearBuilder{K,C},
    form::LinearMomentForm{K,C},
    origin::ConstraintOrigin;
    kind::Symbol=:zero,
) where {K,C}
    _check_not_finalized!(builder)
    push!(
        builder.zero_constraints,
        ScalarLinearConstraint(_owned_builder_form(K, C, form), kind, origin),
    )
    builder.zero_constraint_keys_registered = false
    return builder
end

"""
Internal append for freshly generated, already-owned linear forms.
Public insertion paths copy caller-owned keys before staging constraints.
"""
function _add_zero_constraint_trusted!(
    builder::MomentLinearBuilder{K,C},
    form::LinearMomentForm{K,C},
    origin::ConstraintOrigin;
    kind::Symbol=:zero,
    trusted_self_adjoint::Bool=false,
) where {K,C}
    _check_not_finalized!(builder)
    push!(
        builder.zero_constraints,
        ScalarLinearConstraint(form, kind, origin; trusted_self_adjoint),
    )
    return builder
end

function _owned_key_to_monomial(builder::MomentLinearBuilder{K,C,M}) where {K,C,M}
    key_to_monomial = Dict{K,M}()
    for (key, mono) in builder.key_to_monomial
        _register_builder_moment_key!(key_to_monomial, K, key, mono)
    end
    return key_to_monomial
end

function _take_builder_psd_blocks!(builder::MomentLinearBuilder{K,C,M}) where {K,C,M}
    blocks = builder.psd_blocks_lin
    builder.psd_blocks_lin = PSDBlockLin{K,C,M}[]
    return blocks
end

function _take_builder_zero_constraints!(builder::MomentLinearBuilder{K,C}) where {K,C}
    constraints = builder.zero_constraints
    builder.zero_constraints = ScalarLinearConstraint{K,C}[]
    return constraints
end

function finalize!(
    builder::MomentLinearBuilder{K,C,M},
    ;
    stage_times_ns=nothing,
    stage_prefix::Symbol=:moment_linear_finalize,
) where {K,C,A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T}}
    _check_not_finalized!(builder)
    builder.finalized = true

    stage_start_ns = time_ns()
    key_to_monomial = _owned_key_to_monomial(builder)
    _close_adjoint_keys!(key_to_monomial, K, A)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_key_copy")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    moments = sort!(collect(keys(key_to_monomial)); lt=key_lt)
    moment_index = Dict{K,Int}(key => idx for (idx, key) in enumerate(moments))
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_moment_index")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    adjoint_key = Dict{K,K}()
    if _is_complex_problem(A)
        for key in moments
            mono = _get_key_value(key_to_monomial, key, "representative monomial")
            adjoint_key[key] = _owned_moment_key(K, _moment_key(K, _moment_linear_adjoint_monomial(mono)))
        end
    else
        for key in moments
            adjoint_key[key] = key
        end
    end
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_adjoint_key")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    objective_lin = _builder_linear_form(K, C, builder.objective_terms)
    psd_blocks_lin = _take_builder_psd_blocks!(builder)
    psd_block_constraint_idx = builder.psd_block_constraint_idx
    builder.psd_block_constraint_idx = Int[]
    zero_constraints = _take_builder_zero_constraints!(builder)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_take_constraints")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    pivots = _discover_linear_pivots(psd_blocks_lin, adjoint_key)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_pivots")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    free_keys = K[key for key in moments if !haskey(pivots, key)]
    pivot_at = _build_pivot_at(pivots)
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_free_keys")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    linear = MomentLinearData{K,C,M}(
        moments,
        moment_index,
        _owned_moment_key(K, builder.identity),
        key_to_monomial,
        adjoint_key,
        psd_blocks_lin,
        psd_block_constraint_idx,
        zero_constraints,
        objective_lin,
        pivots,
        pivot_at,
        free_keys,
        stage_times_ns=stage_times_ns,
        stage_prefix=Symbol(stage_prefix, "_construct_data"),
        zero_constraint_keys_registered=builder.zero_constraint_keys_registered,
    )
    if stage_times_ns !== nothing
        stage_key = Symbol(stage_prefix, "_construct_data")
        stage_times_ns[stage_key] =
            get(stage_times_ns, stage_key, 0) + Int(time_ns() - stage_start_ns)
    end
    return linear
end
