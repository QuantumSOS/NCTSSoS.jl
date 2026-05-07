# =============================================================================
# Structured periodic PQG V2RDM assembler
# =============================================================================

"""
Structured V2RDM block provenance for PQG moment data built without the generic
non-commutative polynomial path.
"""
struct V2RDMBlockOrigin <: BlockOrigin
    family::Symbol
    block_key::Tuple
end

const V2RDM_SPINS = (1, -1) # α, β as two_Sz signs

@inline _v2rdm_mode_index(k::Integer, spin::Integer, orb::Integer, norb::Integer) =
    2 * Int(k) * Int(norb) + (spin == 1 ? 0 : Int(norb)) + Int(orb)

@inline function _v2rdm_mode_k(mode::Integer, norb::Integer)
    return (Int(abs(mode)) - 1) ÷ (2 * Int(norb))
end

@inline function _v2rdm_mode_spin(mode::Integer, norb::Integer)
    family = (Int(abs(mode)) - 1) ÷ Int(norb)
    return iseven(family) ? 1 : -1
end

@inline function _v2rdm_particle_sign(op::Integer)
    return op < 0 ? 1 : -1
end

@inline function _v2rdm_normal_key(op::Integer)
    return (op > 0, op)
end

function _v2rdm_is_normal_ordered(word::AbstractVector{<:Signed})
    @inbounds for i in 1:(length(word) - 1)
        _v2rdm_normal_key(word[i]) > _v2rdm_normal_key(word[i + 1]) && return false
    end
    return true
end

function _v2rdm_first_out_of_order(word::AbstractVector{<:Signed})
    @inbounds for i in 1:(length(word) - 1)
        _v2rdm_normal_key(word[i]) > _v2rdm_normal_key(word[i + 1]) && return i
    end
    return 0
end

function _v2rdm_has_duplicate_operator(word::AbstractVector)
    @inbounds for i in 2:length(word)
        word[i] == word[i - 1] && return true
    end
    return false
end

function _v2rdm_accumulate_word!(out::Vector{Pair{Vector{T},Int}}, word::Vector{T}, coeff::Integer) where {T<:Signed}
    iszero(coeff) && return out
    for idx in eachindex(out)
        if key_isequal(out[idx].first, word)
            combined = out[idx].second + Int(coeff)
            if iszero(combined)
                deleteat!(out, idx)
            else
                out[idx] = out[idx].first => combined
            end
            return out
        end
    end
    push!(out, copy(word) => Int(coeff))
    return out
end

"""
    normal_order_terms(raw) -> Vector{Pair{Vector{T},Int}}

Small, index-level fermionic normal-orderer for the PQG fast path.  It handles
the degree-≤4 words produced by D/Q/G blocks, the Hamiltonian, and the closed
form equality rows.  It deliberately avoids `Polynomial`, `simplify`, and the
PBW canonicalization stack.
"""
function normal_order_terms(raw::Vector{T}) where {T<:Signed}
    pending = Pair{Vector{T},Int}[copy(raw) => 1]
    out = Pair{Vector{T},Int}[]

    while !isempty(pending)
        word, coeff = pop!(pending)
        iszero(coeff) && continue

        idx = _v2rdm_first_out_of_order(word)
        if idx == 0
            _v2rdm_has_duplicate_operator(word) && continue
            _v2rdm_accumulate_word!(out, word, coeff)
            continue
        end

        a = word[idx]
        b = word[idx + 1]

        # Swap distinct fermionic operators: xy = -yx, plus the contraction
        # a_i a_i† = 1 - a_i† a_i when an annihilator crosses its own creator.
        swapped = copy(word)
        swapped[idx], swapped[idx + 1] = swapped[idx + 1], swapped[idx]
        push!(pending, swapped => -coeff)

        if a > 0 && b < 0 && a == -b
            contracted = Vector{T}(undef, length(word) - 2)
            out_idx = 0
            @inbounds for src_idx in eachindex(word)
                (src_idx == idx || src_idx == idx + 1) && continue
                out_idx += 1
                contracted[out_idx] = word[src_idx]
            end
            push!(pending, contracted => coeff)
        end
    end

    sort!(out; by = first, lt = key_lt)
    return out
end

@inline function _v2rdm_adjoint_word(word::Vector{T}) where {T<:Signed}
    out = Vector{T}(undef, length(word))
    @inbounds for i in eachindex(word)
        out[i] = -word[length(word) - i + 1]
    end
    return out
end

@inline function _v2rdm_raw_dot3(a::Vector{T}, m::Vector{T}, b::Vector{T}) where {T<:Signed}
    out = Vector{T}(undef, length(a) + length(m) + length(b))
    @inbounds begin
        for i in eachindex(a)
            out[i] = -a[length(a) - i + 1]
        end
        o = length(a)
        for i in eachindex(m)
            out[o + i] = m[i]
        end
        o += length(m)
        for i in eachindex(b)
            out[o + i] = b[i]
        end
    end
    return out
end

function _v2rdm_block_label(word::Vector{T}; nk::Int, norb::Int) where {T<:Signed}
    K = 0
    two_Sz = 0
    @inbounds for op in word
        ps = _v2rdm_particle_sign(op)
        K = mod(K + ps * _v2rdm_mode_k(op, norb), nk)
        two_Sz += ps * _v2rdm_mode_spin(op, norb)
    end
    return (; K, two_Sz)
end

function _v2rdm_block_key(word::Vector{T}, blocking::Symbol; nk::Int, norb::Int) where {T<:Signed}
    label = _v2rdm_block_label(word; nk, norb)
    if blocking == :none
        return ()
    elseif blocking == :momentum
        return (label.K,)
    elseif blocking == :spin
        return (label.K, label.two_Sz)
    else
        throw(ArgumentError("unknown blocking $(repr(blocking)); expected :momentum, :spin, or :none"))
    end
end

function _v2rdm_group_rows(rows::Vector{M}, blocking::Symbol; nk::Int, norb::Int) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    blocking == :none && return Tuple{Tuple,Vector{M}}[((), rows)]

    buckets = Dict{Tuple,Vector{M}}()
    for row in rows
        key = _v2rdm_block_key(row.word, blocking; nk, norb)
        push!(get!(buckets, key, M[]), row)
    end
    return Tuple{Tuple,Vector{M}}[(key, buckets[key]) for key in sort!(collect(keys(buckets)))]
end

function _v2rdm_group_rows_with_identity(rows::Vector{M}, identity::M, blocking::Symbol; nk::Int, norb::Int) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    if blocking == :none
        return Tuple{Tuple,Vector{M}}[((), M[identity; rows])]
    end

    buckets = Dict{Tuple,Vector{M}}()
    for row in rows
        key = _v2rdm_block_key(row.word, blocking; nk, norb)
        push!(get!(buckets, key, M[]), row)
    end
    identity_key = blocking == :momentum ? (0,) : (0, 0)
    pushfirst!(get!(buckets, identity_key, M[]), identity)
    return Tuple{Tuple,Vector{M}}[(key, buckets[key]) for key in sort!(collect(keys(buckets)))]
end

function _v2rdm_index_type(n_modes::Integer)
    n_modes <= typemax(Int8) && return Int8
    n_modes <= typemax(Int16) && return Int16
    n_modes <= typemax(Int32) && return Int32
    return Int64
end

function _v2rdm_monomial(::Type{T}, word) where {T<:Signed}
    return NormalMonomial{FermionicAlgebra,T}(T[x for x in word])
end

function _v2rdm_row_bases(::Type{T}; nk::Int, norb::Int, include_one_d::Bool) where {T<:Signed}
    n_modes = 2 * nk * norb
    M = NormalMonomial{FermionicAlgebra,T}

    oneD = M[]
    D = M[]
    Q = M[]
    G = M[]

    if include_one_d
        for p in 1:n_modes
            push!(oneD, _v2rdm_monomial(T, (p,)))
        end
    end

    for p in 1:(n_modes - 1), q in (p + 1):n_modes
        push!(D, _v2rdm_monomial(T, (p, q)))
        # This matches the symbolic basis: only_monomial(a†_p a†_q) stores the
        # canonical monomial a†_q a†_p and discards the construction sign.
        push!(Q, _v2rdm_monomial(T, (-q, -p)))
    end

    for p in 1:n_modes, q in 1:n_modes
        push!(G, _v2rdm_monomial(T, (-p, q)))
    end

    return (; oneD, D, Q, G)
end

# Public-by-name helper requested by the benchmark task.  The assembler itself
# uses the lower-level word routines to avoid the symbolic polynomial stack.
function canonical_moment_id(kind::Symbol, K::Integer, two_Sz::Integer, multi_index)
    if kind == :identity
        return Int[]
    elseif kind == :word
        return collect(Int, multi_index)
    elseif kind == :one
        p, q = multi_index
        terms = normal_order_terms(Int[-p, q])
        isempty(terms) && return Int[]
        return first(only(terms))
    elseif kind == :raw
        terms = normal_order_terms(collect(Int, multi_index))
        length(terms) == 1 || throw(ArgumentError("raw canonical id expanded to $(length(terms)) terms for K=$K two_Sz=$two_Sz"))
        return first(only(terms))
    else
        throw(ArgumentError("unknown canonical moment kind $(repr(kind))"))
    end
end

mutable struct _V2RDMAccumulator{T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    key_to_monomial::Dict{Vector{T},M}
end

function _v2rdm_get_key(dict::Dict{Vector{T},V}, key::Vector{T}) where {T,V}
    haskey(dict, key) && return getkey(dict, key, key)
    for candidate in keys(dict)
        key_isequal(candidate, key) && return candidate
    end
    return nothing
end

function _v2rdm_register_key!(acc::_V2RDMAccumulator{T,M}, key::Vector{T}) where {T,M}
    existing = _v2rdm_get_key(acc.key_to_monomial, key)
    existing !== nothing && return existing
    stored = copy(key)
    acc.key_to_monomial[stored] = M(stored)
    return stored
end

function _v2rdm_add_term!(pairs::Vector{Pair{Vector{T},ComplexF64}}, acc::_V2RDMAccumulator{T}, key::Vector{T}, coeff) where {T}
    c = ComplexF64(coeff)
    iszero(c) && return pairs
    interned = _v2rdm_register_key!(acc, key)
    push!(pairs, interned => c)
    return pairs
end

function _v2rdm_add_normal_ordered!(pairs::Vector{Pair{Vector{T},ComplexF64}}, acc::_V2RDMAccumulator{T}, raw::Vector{T}, coeff) where {T}
    c = ComplexF64(coeff)
    iszero(c) && return pairs
    for key_coeff in normal_order_terms(raw)
        _v2rdm_add_term!(pairs, acc, key_coeff.first, c * key_coeff.second)
    end
    return pairs
end

function _v2rdm_add_hermitian_raw!(pairs::Vector{Pair{Vector{T},ComplexF64}}, acc::_V2RDMAccumulator{T}, raw::Vector{T}, coeff) where {T}
    c = ComplexF64(coeff)
    iszero(c) && return pairs
    _v2rdm_add_normal_ordered!(pairs, acc, raw, 0.5 * c)
    _v2rdm_add_normal_ordered!(pairs, acc, _v2rdm_adjoint_word(raw), 0.5 * conj(c))
    return pairs
end

function _v2rdm_form(acc::_V2RDMAccumulator{T}, pairs::Vector{Pair{Vector{T},ComplexF64}}) where {T}
    return LinearMomentForm{Vector{T},ComplexF64}(pairs)
end

function _v2rdm_empty_form(::Type{T}) where {T<:Signed}
    return LinearMomentForm{Vector{T},ComplexF64}(Pair{Vector{T},ComplexF64}[])
end

function _v2rdm_adjoint_key(key::Vector{T}) where {T<:Signed}
    # Adjoint of a normal fermionic word is normal in this convention.
    return _v2rdm_adjoint_word(key)
end

function _v2rdm_close_adjoint_keys!(acc::_V2RDMAccumulator{T}) where {T}
    idx = 1
    keys_to_visit = collect(keys(acc.key_to_monomial))
    while idx <= length(keys_to_visit)
        key = keys_to_visit[idx]
        adj = _v2rdm_adjoint_key(key)
        existing = _v2rdm_get_key(acc.key_to_monomial, adj)
        if existing === nothing
            stored = _v2rdm_register_key!(acc, adj)
            push!(keys_to_visit, stored)
        end
        idx += 1
    end
    return nothing
end

function _v2rdm_raw_real_part_form(raw::LinearMomentForm{Vector{T},ComplexF64}, acc::_V2RDMAccumulator{T}, adjoint_key::Dict{Vector{T},Vector{T}}) where {T}
    pairs = Pair{Vector{T},ComplexF64}[]
    for (key, coef) in raw
        adj = _get_key_value(adjoint_key, key, "adjoint key")
        _v2rdm_add_term!(pairs, acc, key, 0.5 * coef)
        _v2rdm_add_term!(pairs, acc, adj, 0.5 * conj(coef))
    end
    return _v2rdm_form(acc, pairs)
end

function _v2rdm_raw_imag_part_form(raw::LinearMomentForm{Vector{T},ComplexF64}, acc::_V2RDMAccumulator{T}, adjoint_key::Dict{Vector{T},Vector{T}}) where {T}
    pairs = Pair{Vector{T},ComplexF64}[]
    for (key, coef) in raw
        adj = _get_key_value(adjoint_key, key, "adjoint key")
        _v2rdm_add_term!(pairs, acc, key, -0.5im * coef)
        _v2rdm_add_term!(pairs, acc, adj, 0.5im * conj(coef))
    end
    return _v2rdm_form(acc, pairs)
end

function _v2rdm_number_operator_terms(::Type{T}; nk::Int, norb::Int) where {T<:Signed}
    n_modes = 2 * nk * norb
    return [T[-mode, mode] for mode in 1:n_modes]
end

function _v2rdm_identity_key!(acc::_V2RDMAccumulator{T,M}) where {T,M}
    return _v2rdm_register_key!(acc, T[])
end

"""
    objective_from_integrals(h1e, eri; nk, norb)

Build the active electronic Hamiltonian objective directly as a linear moment
form from PySCF integrals.  The integral normalization matches the symbolic H4
exporter: h1e/Nk and eri/Nk², with the usual 1/2 two-electron prefactor.
"""
function objective_from_integrals(h1e, eri, acc::_V2RDMAccumulator{T}; nk::Int, norb::Int) where {T<:Signed}
    pairs = Pair{Vector{T},ComplexF64}[]

    for k in 0:(nk - 1)
        haskey(h1e, k) || throw(ArgumentError("missing h1e block for k=$k"))
        size(h1e[k], 1) >= norb && size(h1e[k], 2) >= norb || throw(DimensionMismatch(
            "h1e[$k] must have at least size ($norb,$norb), got $(size(h1e[k]))"
        ))
        hk = h1e[k] / nk
        for p in 1:norb, q in 1:norb, spin in V2RDM_SPINS
            coeff = hk[p, q]
            iszero(coeff) && continue
            mp = _v2rdm_mode_index(k, spin, p, norb)
            mq = _v2rdm_mode_index(k, spin, q, norb)
            _v2rdm_add_hermitian_raw!(pairs, acc, T[-mp, mq], coeff)
        end
    end

    spin_channels = ((1, 1), (1, -1), (-1, 1), (-1, -1))
    for key in sort!(collect(keys(eri)))
        k1, k2, k3, k4 = key
        all(k -> 0 <= k < nk, key) || throw(ArgumentError(
            "ERI block $key has k index outside 0:$(nk - 1)"
        ))
        size(eri[key]) == (norb, norb, norb, norb) || throw(DimensionMismatch(
            "eri[$key] must have size ($norb,$norb,$norb,$norb), got $(size(eri[key]))"
        ))
        mod(k1 + k2 - k3 - k4, nk) == 0 || throw(ArgumentError(
            "ERI block $key violates momentum conservation modulo nk=$nk"
        ))
        V = eri[key] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff0 = V[p, r, q, s]
            iszero(coeff0) && continue
            coeff = 0.5 * coeff0
            for (σ, τ) in spin_channels
                mp = _v2rdm_mode_index(k1, σ, p, norb)
                mq = _v2rdm_mode_index(k2, τ, q, norb)
                mr = _v2rdm_mode_index(k3, σ, r, norb)
                ms = _v2rdm_mode_index(k4, τ, s, norb)
                _v2rdm_add_hermitian_raw!(pairs, acc, T[-mp, -mq, ms, mr], coeff)
            end
        end
    end

    return _v2rdm_form(acc, pairs)
end

function objective_from_integrals(h1e, eri; nk::Int, norb::Int)
    T = _v2rdm_index_type(2 * nk * norb)
    M = NormalMonomial{FermionicAlgebra,T}
    acc = _V2RDMAccumulator{T,M}(Dict{Vector{T},M}())
    _v2rdm_identity_key!(acc)
    return objective_from_integrals(h1e, eri, acc; nk, norb)
end

"""
    pqg_block_entries(family, K, two_Sz, rows, acc)

Construct one PQG block as closed-form linear moment entries.  `family` is one
of `:D`, `:Q`, `:G`, or `:oneD`; `K`/`two_Sz` are metadata for the residual
symmetry block.
"""
function pqg_block_entries(family::Symbol, K, two_Sz, rows::Vector{M}, acc::_V2RDMAccumulator{T}) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    n = length(rows)
    entries = Matrix{LinearMomentForm{Vector{T},ComplexF64}}(undef, n, n)
    empty = _v2rdm_empty_form(T)

    for j in 1:n, i in 1:n
        raw = _v2rdm_raw_dot3(rows[i].word, T[], rows[j].word)
        pairs = Pair{Vector{T},ComplexF64}[]
        _v2rdm_add_normal_ordered!(pairs, acc, raw, 1.0)
        entries[i, j] = isempty(pairs) ? empty : _v2rdm_form(acc, pairs)
    end
    return entries
end

function _v2rdm_raw_scalar_constraint_form(acc::_V2RDMAccumulator{T}, terms::Vector{Pair{Vector{T},ComplexF64}}, rhs) where {T}
    pairs = copy(terms)
    _v2rdm_add_term!(pairs, acc, T[], -ComplexF64(rhs))
    return _v2rdm_form(acc, pairs)
end

function _v2rdm_number_row_form(acc::_V2RDMAccumulator{T}, row::M; nk::Int, norb::Int, total_electrons::Int) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    pairs = Pair{Vector{T},ComplexF64}[]
    for nword in _v2rdm_number_operator_terms(T; nk, norb)
        _v2rdm_add_normal_ordered!(pairs, acc, _v2rdm_raw_dot3(row.word, nword, T[]), 1.0)
    end
    # -N * b†
    _v2rdm_add_normal_ordered!(pairs, acc, _v2rdm_raw_dot3(row.word, T[], T[]), -Float64(total_electrons))
    return _v2rdm_form(acc, pairs)
end

function _v2rdm_trace_q_form(acc::_V2RDMAccumulator{T}, Q_rows::Vector{M}) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    pairs = Pair{Vector{T},ComplexF64}[]
    for row in Q_rows
        _v2rdm_add_normal_ordered!(pairs, acc, _v2rdm_raw_dot3(row.word, T[], row.word), 1.0)
    end
    return pairs
end

function _v2rdm_trace_g_form(acc::_V2RDMAccumulator{T}, G_rows::Vector{M}) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    pairs = Pair{Vector{T},ComplexF64}[]
    for row in G_rows
        _v2rdm_add_normal_ordered!(pairs, acc, _v2rdm_raw_dot3(row.word, T[], row.word), 1.0)
    end
    return pairs
end

function _v2rdm_trace_d_spin_forms(acc::_V2RDMAccumulator{T}, D_rows::Vector{M}; norb::Int) where {T<:Signed,M<:NormalMonomial{FermionicAlgebra,T}}
    aa = Pair{Vector{T},ComplexF64}[]
    ab = Pair{Vector{T},ComplexF64}[]
    bb = Pair{Vector{T},ComplexF64}[]
    for row in D_rows
        p, q = row.word
        target = if _v2rdm_mode_spin(p, norb) == 1 && _v2rdm_mode_spin(q, norb) == 1
            aa
        elseif _v2rdm_mode_spin(p, norb) == -1 && _v2rdm_mode_spin(q, norb) == -1
            bb
        else
            ab
        end
        _v2rdm_add_normal_ordered!(target, acc, _v2rdm_raw_dot3(row.word, T[], row.word), 1.0)
    end
    return (; aa, ab, bb)
end

function _v2rdm_singlet_s2_form(acc::_V2RDMAccumulator{T}; nk::Int, norb::Int, total_electrons::Int) where {T<:Signed}
    pairs = Pair{Vector{T},ComplexF64}[]
    for ki in 0:(nk - 1), pi in 1:norb, kj in 0:(nk - 1), pj in 1:norb
        up_i = _v2rdm_mode_index(ki, 1, pi, norb)
        dn_j = _v2rdm_mode_index(kj, -1, pj, norb)
        dn_i = _v2rdm_mode_index(ki, -1, pi, norb)
        up_j = _v2rdm_mode_index(kj, 1, pj, norb)
        _v2rdm_add_hermitian_raw!(pairs, acc, T[-up_i, -dn_j, dn_i, up_j], 1.0)
    end
    _v2rdm_add_term!(pairs, acc, T[], -0.5 * total_electrons)
    return _v2rdm_form(acc, pairs)
end

function _v2rdm_ensure_adjoint_keys!(
    acc::_V2RDMAccumulator{T},
    adjoint_key::Dict{Vector{T},Vector{T}},
    raw::LinearMomentForm{Vector{T},ComplexF64},
) where {T}
    for (key, _) in raw
        stored = _v2rdm_register_key!(acc, key)
        adj = _v2rdm_register_key!(acc, _v2rdm_adjoint_key(stored))
        adjoint_key[stored] = adj
        adjoint_key[adj] = stored
    end
    return adjoint_key
end

function _v2rdm_append_complex_zero_rows!(rows::Vector{ScalarLinearConstraint{Vector{T},ComplexF64}}, raw::LinearMomentForm{Vector{T},ComplexF64}, acc::_V2RDMAccumulator{T}, adjoint_key::Dict{Vector{T},Vector{T}}, cons_idx::Int, row_idx::Int) where {T}
    _v2rdm_ensure_adjoint_keys!(acc, adjoint_key, raw)

    real_form = _v2rdm_raw_real_part_form(raw, acc, adjoint_key)
    isempty(real_form) || push!(rows, ScalarLinearConstraint(real_form, :zero, MomentEqOrigin(cons_idx, row_idx)))

    imag_form = _v2rdm_raw_imag_part_form(raw, acc, adjoint_key)
    isempty(imag_form) || push!(rows, ScalarLinearConstraint(imag_form, :zero, MomentEqOrigin(cons_idx, row_idx)))
    return rows
end

"""
    pqg_equality_rows(...)

Closed-form scalar equality generator for the structured PQG relaxation.  It
matches the symbolic benchmarker's moment-equality convention: the particle
number polynomial is applied one-sided to every degree-≤2 PQG row, while the
quartic trace/S² constraints only use the identity row.
"""
function pqg_equality_rows(
    acc::_V2RDMAccumulator{T},
    adjoint_key::Dict{Vector{T},Vector{T}},
    row_bases::NamedTuple;
    nk::Int,
    norb::Int,
    nelec_per_cell::Int,
    spin_resolved_trace::Bool,
    singlet_s2::Bool,
) where {T<:Signed}
    M = NormalMonomial{FermionicAlgebra,T}
    total_electrons = nk * nelec_per_cell
    if spin_resolved_trace && isodd(total_electrons)
        throw(ArgumentError("spin-resolved D traces assume M_s = 0 and an even total electron count"))
    end
    n_modes = 2 * nk * norb
    n_holes = n_modes - total_electrons
    rows = ScalarLinearConstraint{Vector{T},ComplexF64}[]

    moment_eq_rows = M[]
    append!(moment_eq_rows, row_bases.oneD)
    append!(moment_eq_rows, row_bases.D)
    append!(moment_eq_rows, row_bases.Q)
    append!(moment_eq_rows, row_bases.G)
    push!(moment_eq_rows, one(M))
    sort!(moment_eq_rows)
    moment_eq_rows = sorted_unique!(moment_eq_rows)

    cons_idx = 1
    row_idx = 1
    for row in moment_eq_rows
        raw = _v2rdm_number_row_form(acc, row; nk, norb, total_electrons)
        _v2rdm_append_complex_zero_rows!(rows, raw, acc, adjoint_key, cons_idx, row_idx)
        row_idx += 1
    end

    # The remaining constraints are quartic, so the symbolic truncation keeps
    # only the identity row.
    cons_idx += 1
    q_pairs = _v2rdm_trace_q_form(acc, row_bases.Q)
    raw_q = _v2rdm_raw_scalar_constraint_form(acc, q_pairs, n_holes * (n_holes - 1) ÷ 2)
    _v2rdm_append_complex_zero_rows!(rows, raw_q, acc, adjoint_key, cons_idx, 1)

    cons_idx += 1
    g_pairs = _v2rdm_trace_g_form(acc, row_bases.G)
    raw_g = _v2rdm_raw_scalar_constraint_form(acc, g_pairs, total_electrons * (n_holes + 1))
    _v2rdm_append_complex_zero_rows!(rows, raw_g, acc, adjoint_key, cons_idx, 1)

    if spin_resolved_trace
        n_alpha = total_electrons ÷ 2
        n_beta = total_electrons ÷ 2
        d_forms = _v2rdm_trace_d_spin_forms(acc, row_bases.D; norb)

        cons_idx += 1
        raw_aa = _v2rdm_raw_scalar_constraint_form(acc, d_forms.aa, n_alpha * (n_alpha - 1) ÷ 2)
        _v2rdm_append_complex_zero_rows!(rows, raw_aa, acc, adjoint_key, cons_idx, 1)

        cons_idx += 1
        raw_ab = _v2rdm_raw_scalar_constraint_form(acc, d_forms.ab, n_alpha * n_beta)
        _v2rdm_append_complex_zero_rows!(rows, raw_ab, acc, adjoint_key, cons_idx, 1)

        cons_idx += 1
        raw_bb = _v2rdm_raw_scalar_constraint_form(acc, d_forms.bb, n_beta * (n_beta - 1) ÷ 2)
        _v2rdm_append_complex_zero_rows!(rows, raw_bb, acc, adjoint_key, cons_idx, 1)
    else
        cons_idx += 1
        d_pairs = Pair{Vector{T},ComplexF64}[]
        for row in row_bases.D
            _v2rdm_add_normal_ordered!(d_pairs, acc, _v2rdm_raw_dot3(row.word, T[], row.word), 1.0)
        end
        raw_d = _v2rdm_raw_scalar_constraint_form(acc, d_pairs, total_electrons * (total_electrons - 1) ÷ 2)
        _v2rdm_append_complex_zero_rows!(rows, raw_d, acc, adjoint_key, cons_idx, 1)
    end

    if singlet_s2
        cons_idx += 1
        raw_s2 = _v2rdm_singlet_s2_form(acc; nk, norb, total_electrons)
        _v2rdm_append_complex_zero_rows!(rows, raw_s2, acc, adjoint_key, cons_idx, 1)
    end

    return rows
end

function _v2rdm_discover_pivots(psd_blocks_lin::Vector{PSDBlockLin{Vector{T},ComplexF64,M}}, adjoint_key::Dict{Vector{T},Vector{T}}) where {T<:Signed,M}
    pivots = Dict{Vector{T},Pivot{ComplexF64}}()
    for (block_idx, block) in enumerate(psd_blocks_lin)
        for i in 1:block.size, j in i:block.size
            form = block.entries[i, j]
            length(form.terms) == 1 || continue
            key, coef = only(form.terms)
            phase = if coef == 1
                1.0 + 0.0im
            elseif coef == -1
                -1.0 + 0.0im
            elseif coef == im
                0.0 + 1.0im
            elseif coef == -im
                0.0 - 1.0im
            else
                nothing
            end
            phase === nothing && continue
            _v2rdm_get_key(pivots, key) === nothing && (pivots[key] = Pivot{ComplexF64}(block_idx, i, j, phase, false))
            adj = _get_key_value(adjoint_key, key, "adjoint key")
            if i != j && !key_isequal(adj, key) && _v2rdm_get_key(pivots, adj) === nothing
                pivots[adj] = Pivot{ComplexF64}(block_idx, i, j, phase, true)
            end
        end
    end
    return pivots
end

function _v2rdm_build_pivot_at(pivots::Dict{Vector{T},Pivot{ComplexF64}}) where {T}
    pivot_at = Dict{Tuple{Int,Int,Int},Vector{Vector{T}}}()
    for (key, pivot) in pivots
        push!(get!(pivot_at, (pivot.block, pivot.row, pivot.col), Vector{T}[]), key)
    end
    for keys_at_position in values(pivot_at)
        sort!(keys_at_position; lt=key_lt)
    end
    return pivot_at
end

function _v2rdm_adjoint_dict(acc::_V2RDMAccumulator{T}) where {T}
    adjoint_key = Dict{Vector{T},Vector{T}}()
    for key in keys(acc.key_to_monomial)
        adj = _v2rdm_register_key!(acc, _v2rdm_adjoint_key(key))
        adjoint_key[key] = adj
    end
    return adjoint_key
end

"""
    build_pqg_moment_data(h1e, eri; nk, norb, nelec_per_cell, blocking,
                          spin_resolved_trace, singlet_s2) -> MomentLinearData

Build the periodic PQG V2RDM moment problem directly as `MomentLinearData` from
integrals and fixed PQG combinatorics.  This bypasses `Polynomial`, `simplify`,
and `moment_relax` on the hot path.
"""
function build_pqg_moment_data(
    h1e,
    eri;
    nk::Int,
    norb::Int,
    nelec_per_cell::Int,
    blocking::Symbol = :momentum,
    spin_resolved_trace::Bool = true,
    singlet_s2::Bool = true,
    include_one_d::Bool = false,
)
    nk >= 1 || throw(ArgumentError("nk must be positive"))
    norb >= 1 || throw(ArgumentError("norb must be positive"))
    nelec_per_cell >= 1 || throw(ArgumentError("nelec_per_cell must be positive"))
    total_electrons = nk * nelec_per_cell
    n_modes = 2 * nk * norb
    total_electrons <= n_modes || throw(ArgumentError(
        "total electrons ($total_electrons) exceed spin-orbital modes ($n_modes)"
    ))
    spin_resolved_trace && isodd(total_electrons) && throw(ArgumentError(
        "spin-resolved D traces assume M_s = 0 and an even total electron count"
    ))
    blocking in (:momentum, :spin, :none) || throw(ArgumentError("blocking must be :momentum, :spin, or :none"))

    T = _v2rdm_index_type(n_modes)
    M = NormalMonomial{FermionicAlgebra,T}
    acc = _V2RDMAccumulator{T,M}(Dict{Vector{T},M}())
    identity = _v2rdm_identity_key!(acc)
    identity_mono = one(M)

    rows = _v2rdm_row_bases(T; nk, norb, include_one_d)

    psd_blocks = PSDBlockLin{Vector{T},ComplexF64,M}[]
    constraint_idx = Int[]

    function append_family!(family::Symbol, grouped)
        for (key, block_rows) in grouped
            isempty(block_rows) && continue
            K = isempty(key) ? nothing : key[1]
            two_Sz = length(key) >= 2 ? key[2] : nothing
            entries = pqg_block_entries(family, K, two_Sz, block_rows, acc)
            meta = BlockMeta{M}(:HPSD, V2RDMBlockOrigin(family, key), block_rows)
            push!(psd_blocks, PSDBlockLin{Vector{T},ComplexF64,M}(length(block_rows), entries, meta))
            push!(constraint_idx, length(constraint_idx) + 1)
        end
    end

    if include_one_d
        append_family!(:oneD, _v2rdm_group_rows(rows.oneD, blocking; nk, norb))
    end
    append_family!(:D, _v2rdm_group_rows(rows.D, blocking; nk, norb))
    append_family!(:Q, _v2rdm_group_rows(rows.Q, blocking; nk, norb))
    append_family!(:G, _v2rdm_group_rows_with_identity(rows.G, identity_mono, blocking; nk, norb))

    objective = objective_from_integrals(h1e, eri, acc; nk, norb)

    # Close under adjoints before splitting complex equality rows.
    _v2rdm_close_adjoint_keys!(acc)
    adjoint_key = _v2rdm_adjoint_dict(acc)

    zero_constraints = pqg_equality_rows(
        acc,
        adjoint_key,
        rows;
        nk,
        norb,
        nelec_per_cell,
        spin_resolved_trace,
        singlet_s2,
    )

    # Equality rows can introduce adjoint partners while being split.
    _v2rdm_close_adjoint_keys!(acc)
    adjoint_key = _v2rdm_adjoint_dict(acc)

    moments = sort!(collect(keys(acc.key_to_monomial)); lt = key_lt)
    moment_index = Dict{Vector{T},Int}(key => idx for (idx, key) in enumerate(moments))
    pivots = _v2rdm_discover_pivots(psd_blocks, adjoint_key)
    free_keys = Vector{T}[key for key in moments if _v2rdm_get_key(pivots, key) === nothing]
    sort!(free_keys; lt = key_lt)
    pivot_at = _v2rdm_build_pivot_at(pivots)

    return MomentLinearData{Vector{T},ComplexF64,M}(
        moments,
        moment_index,
        identity,
        acc.key_to_monomial,
        adjoint_key,
        psd_blocks,
        constraint_idx,
        zero_constraints,
        objective,
        pivots,
        pivot_at,
        free_keys,
    )
end
