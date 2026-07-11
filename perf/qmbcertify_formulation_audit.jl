#!/usr/bin/env julia

# Compare QMBCertify's source-level Heisenberg A1 formulation pieces against
# the current NCTSSoS source-like finite-axis construction.  This is a
# no-solver diagnostic: it builds row/basis sets, but does not lower or solve
# an SDP.
#
# Large L values require NCTS_QMB_ALLOW_LARGE=true plus explicit wall/RSS
# estimates and safe load/memory telemetry before QMBCertify or NCTSSoS load.

Base.include(@__MODULE__, joinpath(@__DIR__, "shared_load_guard.jl"))

function _formulation_audit_preimport_large_run_pressure_guard()
    _ncts_load_guard_parse_bool("NCTS_QMB_ALLOW_LARGE", false) ||
        return nothing
    lengths_label = strip(get(ENV, "NCTS_FORMULATION_AUDIT_NS", "20"))
    _ncts_check_large_run_pressure_guard(
        env_prefix="NCTS_QMB",
        label="QMBCertify formulation audit Ls=$lengths_label",
    )
    return nothing
end

_formulation_audit_preimport_large_run_pressure_guard()

using LinearAlgebra
using Pkg
using Printf
using TOML

using NCTSSoS

include(joinpath(@__DIR__, "qmbcertify_reference_runs.jl"))

function _parse_int_env(name::AbstractString, default::Int)
    raw = strip(get(ENV, name, string(default)))
    isempty(raw) && return default
    return parse(Int, raw)
end

function _audit_lengths()
    return _parse_int_list_env("NCTS_FORMULATION_AUDIT_NS", "20")
end

function _qmb_word_key(word)
    isempty(word) && return "I"
    return join(string.(Int.(word)), ",")
end

function _qmb_rdm_words!(qmb::Module, tsupp::Vector{Vector{UInt16}}, len::Int, rdm)
    rdm == false && return nothing
    kk = Int(rdm)
    kk in (8, 9, 10) || return nothing

    ranges = ntuple(_ -> 0:3, kk)
    sites = collect(1:kk)
    for ind_tuple in Iterators.product(ranges...)
        ind = collect(ind_tuple)
        all(axis -> iseven(count(==(axis), ind)), 1:3) || continue
        word = UInt16[]
        for (site, axis) in zip(sites, ind)
            axis == 0 && continue
            push!(word, UInt16(3 * (site - 1) + axis))
        end
        reduced = Base.invokelatest(qmb.reduce4, word, len; lattice="chain")
        push!(tsupp, Vector{UInt16}(reduced))
    end
    return nothing
end

function _qmb_chain_basis_and_tsupp(qmb::Module, len::Int, args)
    d = Int(args["d"])
    extra = Int(get(args, "extra", 0))
    pso = Int(get(args, "pso", 0))
    rdm = _rdm_value(get(args, "rdm", false))

    basis = Vector{Vector{Vector{UInt16}}}(undef, 2)
    tsupp = Vector{UInt16}[UInt16[]]

    for i in 1:2
        basis[i] = Base.invokelatest(
            qmb.get_basis,
            len,
            i - 1,
            d;
            lattice="chain",
            extra,
            three_type=[1; 1],
        )
        k = length(basis[i]) ÷ len
        for j in 1:k
            for r in 1:Int(len / 2)
                word = UInt16[basis[i][len * (j - 1) + 1]; basis[i][len * (j - 1) + r + 1]]
                reduced, coef = Base.invokelatest(
                    qmb.reduce!,
                    word;
                    L=len,
                    lattice="chain",
                    realify=true,
                )
                coef != 0 && push!(tsupp, Vector{UInt16}(reduced))
            end
        end
        for j1 in 1:(k - 1), j2 in (j1 + 1):k
            for r in 1:len
                word = UInt16[basis[i][len * (j1 - 1) + 1]; basis[i][len * (j2 - 1) + r]]
                reduced, coef = Base.invokelatest(
                    qmb.reduce!,
                    word;
                    L=len,
                    lattice="chain",
                    realify=true,
                )
                coef != 0 && push!(tsupp, Vector{UInt16}(reduced))
            end
        end

        if pso > 0
            # PSD state-opt support is not needed for A1, but include the
            # pso-local basis count in the audit output when requested.
            count(word -> length(word) <= pso, basis[i])
        end
    end

    _qmb_rdm_words!(qmb, tsupp, len, rdm)
    sort!(tsupp)
    unique!(tsupp)
    return basis, tsupp
end

function _qmb_lso_row_key(qmb::Module, mon::Vector{UInt16}, len::Int)
    row = Dict{String,Float64}()
    for i in 1:len, axis in 1:3
        word = UInt16[3 * (i - 1) + axis; 3 * mod(i, len) + axis; mon]
        reduced, coef = Base.invokelatest(qmb.reduce!, word; L=len, lattice="chain")
        imag_coef = imag(coef)
        iszero(imag_coef) && continue
        key = _qmb_word_key(reduced)
        row[key] = get(row, key, 0.0) + 0.75 * Float64(imag_coef)
    end
    return _row_key(row)
end

function _row_key(row::Dict{String,Float64})
    pairs = [(key, value) for (key, value) in row if !iszero(value)]
    sort!(pairs; by=first)
    isempty(pairs) && return ""
    return join(("$key=$(@sprintf("%.12g", value))" for (key, value) in pairs), ";")
end

function _nctssos_lso_row_key(qmb::Module, test_mono, len::Int)
    mon = UInt16.(test_mono.word)
    return _qmb_lso_row_key(qmb, mon, len)
end

function _set_sample(values::Set{String}; limit::Int=3)
    sample = collect(values)
    sort!(sample)
    length(sample) > limit && resize!(sample, limit)
    return sample
end

function _audit_histogram_pairs(values::AbstractVector{<:Integer})
    counts = Dict{Int,Int}()
    for value in values
        counts[Int(value)] = get(counts, Int(value), 0) + 1
    end
    return sort!(collect(counts); by=first)
end

function _print_pair_hist(label::AbstractString, values::AbstractVector{<:Integer})
    println("$label=$(_audit_histogram_pairs(values))")
end

function _qmb_word_to_nctssos(::Type{M}, word::Vector{UInt16}) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    simplified, _ = NCTSSoS.simplify(PauliAlgebra, T.(word))
    return M(copy(simplified))
end

function _qmb_basis_family_records(
    basis::Vector{Vector{Vector{UInt16}}},
    len::Int,
    ::Type{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    records = NamedTuple[]
    for parity in 1:2
        k = length(basis[parity]) ÷ len
        for family in 1:k
            word = basis[parity][len * (family - 1) + 1]
            mono = _qmb_word_to_nctssos(M, word)
            rep = NCTSSoS._translation_orbit_representative(mono, len)
            orbit_length = NCTSSoS._translation_orbit_length(rep, len)
            push!(
                records,
                (
                    parity=parity - 1,
                    family,
                    word=Int.(word),
                    word_length=length(word),
                    orbit_length,
                    representative=rep,
                ),
            )
        end
    end
    return records
end

function _qmb_basis_as_nctssos(
    basis::Vector{Vector{Vector{UInt16}}},
    len::Int,
    ::Type{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    monos = M[one(M)]
    for parity in 1:2
        for word in basis[parity]
            push!(monos, _qmb_word_to_nctssos(M, word))
        end
    end
    return sort!(unique!(monos))
end

function _audit_records_sample(records; limit::Int=5)
    sample = collect(records)
    length(sample) > limit && resize!(sample, limit)
    return [
        (
            parity=record.parity,
            family=record.family,
            word=record.word,
            word_length=record.word_length,
            orbit_length=record.orbit_length,
        )
        for record in sample
    ]
end

function _objective_coefficients_by_key(poly)
    coeffs = Dict{String,Float64}()
    for (coef, mono) in poly
        key = _qmb_word_key(UInt16.(mono.word))
        coeffs[key] = get(coeffs, key, 0.0) + Float64(real(coef))
    end
    return filter!(pair -> !iszero(pair.second), coeffs)
end

function _qmb_total_objective_coefficients_by_key(
    args,
    len::Int,
    ::Type{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    supp, coe = _heisenberg_input(args)
    coeffs = Dict{String,Float64}()
    for (word, coef) in zip(supp, coe)
        mono = NCTSSoS._qmbcertify_chain_support_monomial(M, word, len)
        mono === nothing && continue
        key = _qmb_word_key(UInt16.(mono.word))
        coeffs[key] = get(coeffs, key, 0.0) + len * Float64(coef)
    end
    return filter!(pair -> !iszero(pair.second), coeffs)
end

function _objective_coefficient_delta_metrics(lhs::Dict{String,Float64}, rhs::Dict{String,Float64})
    keys_union = union(keys(lhs), keys(rhs))
    max_abs_delta = isempty(keys_union) ? 0.0 : maximum(
        abs(get(lhs, key, 0.0) - get(rhs, key, 0.0)) for key in keys_union
    )
    lhs_only = setdiff(Set(keys(lhs)), Set(keys(rhs)))
    rhs_only = setdiff(Set(keys(rhs)), Set(keys(lhs)))
    return (
        max_abs_delta=max_abs_delta,
        lhs_only_count=length(lhs_only),
        rhs_only_count=length(rhs_only),
        lhs_only_sample=_set_sample(lhs_only),
        rhs_only_sample=_set_sample(rhs_only),
    )
end

function _nctssos_qmb_objective_metrics(
    pop,
    args,
    len::Int,
    ::Type{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    MP = Polynomial{PauliAlgebra,T,Float64}
    reduced = NCTSSoS._qmbcertify_reduce_chain_polynomial(pop.objective, len, MP)
    nctssos_coeffs = _objective_coefficients_by_key(reduced)
    qmb_coeffs = _qmb_total_objective_coefficients_by_key(args, len, M)
    delta = _objective_coefficient_delta_metrics(nctssos_coeffs, qmb_coeffs)
    return (
        nctssos_term_count=length(nctssos_coeffs),
        qmb_total_term_count=length(qmb_coeffs),
        max_abs_delta=delta.max_abs_delta,
        nctssos_only_count=delta.lhs_only_count,
        qmb_only_count=delta.rhs_only_count,
        nctssos_only_sample=delta.lhs_only_sample,
        qmb_only_sample=delta.rhs_only_sample,
        match=delta.max_abs_delta <= 1e-12 &&
            delta.lhs_only_count == 0 &&
            delta.rhs_only_count == 0,
    )
end

function _qmb_source_rdm_blocks(qmb_path::AbstractString, k::Int)
    source_path = joinpath(qmb_path, "src", "rdm_positivity.jl")
    source = read(source_path, String)
    pattern = Regex("function posepsd$(k)!.*?blocks = (\\[.*?\\])\\n    pos", "s")
    found = match(pattern, source)
    found === nothing && throw(ArgumentError(
        "Could not find posepsd$(k)! block rows in $(repr(source_path))."
    ))
    parsed = eval(Meta.parse(only(found.captures)))
    return [Int.(block) for block in parsed]
end

function _qmb_source_rdm_entry(row::Int, code0::Int, k::Int)
    row0 = row - 1
    col0 = 0
    coeff = ComplexF64(1)
    code = code0
    for site in 1:k
        axis = mod(code, 4)
        code = div(code, 4)
        bitpos = k - site
        row_bit = (row0 >> bitpos) & 1
        col_bit = row_bit
        if axis == 1
            col_bit = 1 - row_bit
        elseif axis == 2
            col_bit = 1 - row_bit
            coeff *= iszero(row_bit) ? -im : im
        elseif axis == 3
            coeff *= iszero(row_bit) ? 1 : -1
        elseif !iszero(axis)
            throw(ArgumentError("Unsupported local Pauli axis $axis."))
        end
        col0 |= col_bit << bitpos
    end
    abs(imag(coeff)) <= 1e-12 || throw(ArgumentError(
        "QMBCertify even-axis RDM entry produced non-real coefficient."
    ))
    return col0 + 1, real(coeff)
end

function _qmbcertify_rdm_source_entry_metrics(
    qmb_path::AbstractString,
    len::Int,
    rdm,
    ::Type{M};
    max_k::Int,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    rdm === false && return (
        skipped=true,
        skip_reason="no_rdm",
        k=0,
        row_block_exact_match=false,
        source_block_sizes=Int[],
        nctssos_block_sizes=Int[],
        entry_match=false,
        entry_mismatch_count=0,
        entry_max_abs_delta=0.0,
        entry_compared_count=0,
        scale=0.0,
    )

    kk = Int(rdm)
    source_blocks = _qmb_source_rdm_blocks(qmb_path, kk)
    nctssos_blocks = pauli_qmbcertify_rdm_blocks(kk; ambient_sites=len)
    row_block_exact_match = source_blocks == nctssos_blocks
    if kk > max_k || !row_block_exact_match
        return (
            skipped=true,
            skip_reason=kk > max_k ? "k_above_entry_max" : "row_block_mismatch",
            k=kk,
            row_block_exact_match,
            source_block_sizes=length.(source_blocks),
            nctssos_block_sizes=length.(nctssos_blocks),
            entry_match=false,
            entry_mismatch_count=0,
            entry_max_abs_delta=0.0,
            entry_compared_count=0,
            scale=inv(float(2^kk)),
        )
    end

    K = Vector{T}
    C = Float64
    nctssos_entries = NCTSSoS._qmbcertify_chain_rdm_linear_blocks(
        len,
        kk,
        nctssos_blocks,
        K,
        C,
    )
    row_to_block = zeros(Int, 1 << kk)
    row_to_local = zeros(Int, 1 << kk)
    for (block_idx, rows) in pairs(source_blocks), (local_idx, row) in pairs(rows)
        row_to_block[row] = block_idx
        row_to_local[row] = local_idx
    end

    scale = inv(float(2^kk))
    source_terms = [
        [Pair{K,C}[] for _ in 1:length(rows), _ in 1:length(rows)]
        for rows in source_blocks
    ]
    for code0 in 0:(4^kk - 1)
        word = NCTSSoS._qmbcertify_even_axis_word(kk, code0)
        word === nothing && continue
        reduced = NCTSSoS._qmbcertify_chain_support_monomial(M, word, len)
        reduced === nothing && continue
        key = NCTSSoS._moment_key(K, reduced)
        for (block_idx, rows) in pairs(source_blocks), (local_idx, row) in pairs(rows)
            col, coeff = _qmb_source_rdm_entry(row, code0, kk)
            row_to_block[col] == block_idx || continue
            push!(
                source_terms[block_idx][local_idx, row_to_local[col]],
                key => scale * coeff,
            )
        end
    end

    mismatch_count = 0
    max_abs_delta = 0.0
    compared_count = 0
    for block_idx in eachindex(source_blocks), idx in eachindex(nctssos_entries[block_idx])
        source_form = NCTSSoS._linear_moment_form_from_owned_pairs!(
            source_terms[block_idx][idx],
        )
        nctssos_form = nctssos_entries[block_idx][idx]
        compared_count += 1
        if length(source_form) != length(nctssos_form)
            mismatch_count += 1
            max_abs_delta = Inf
            continue
        end
        for (source_term, nctssos_term) in zip(source_form, nctssos_form)
            if source_term.first != nctssos_term.first
                mismatch_count += 1
                max_abs_delta = Inf
                break
            end
            max_abs_delta = max(max_abs_delta, abs(source_term.second - nctssos_term.second))
        end
    end

    return (
        skipped=false,
        skip_reason="",
        k=kk,
        row_block_exact_match,
        source_block_sizes=length.(source_blocks),
        nctssos_block_sizes=length.(nctssos_blocks),
        entry_match=mismatch_count == 0 && max_abs_delta <= 1e-12,
        entry_mismatch_count=mismatch_count,
        entry_max_abs_delta=max_abs_delta,
        entry_compared_count=compared_count,
        scale,
    )
end

function _qmb_source_reduce_key(
    qmb::Module,
    ::Type{K},
    ::Type{M},
    word::AbstractVector{UInt16},
    len::Int;
    realify::Bool=false,
) where {K,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    reduced, coef = Base.invokelatest(
        qmb.reduce!,
        copy(Vector{UInt16}(word));
        L=len,
        lattice="chain",
        realify,
    )
    iszero(coef) && return nothing, 0.0
    mono = _qmb_word_to_nctssos(M, Vector{UInt16}(reduced))
    return NCTSSoS._moment_key(K, mono), Float64(real(coef))
end

function _push_qmb_source_entry!(
    entries::Matrix{Vector{Pair{K,C}}},
    row::Int,
    col::Int,
    key,
    coef,
    ::Type{C};
    atol::Real=1e-12,
) where {K,C}
    key === nothing && return nothing
    value = convert(C, coef)
    abs(value) <= atol && return nothing
    i, j = row <= col ? (row, col) : (col, row)
    push!(entries[i, j], key => value)
    return nothing
end

function _qmb_source_form_from_pairs!(
    pairs::Vector{Pair{K,C}};
    atol::Real=1e-12,
) where {K,C}
    form = NCTSSoS._linear_moment_form_from_owned_pairs!(pairs)
    cleaned = Pair{K,C}[]
    sizehint!(cleaned, length(form))
    for (key, coef) in form
        abs(coef) <= atol && continue
        push!(cleaned, key => coef)
    end
    return NCTSSoS._linear_moment_form_from_owned_pairs!(cleaned)
end

function _qmb_source_identity_key(::Type{K}, ::Type{M}) where {K,M}
    return NCTSSoS._moment_key(K, one(M))
end

function _qmb_source_base_tables(
    qmb::Module,
    basis::Vector{Vector{Vector{UInt16}}},
    len::Int,
    ::Type{K},
    ::Type{M},
) where {K,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    tables = Any[]
    half = div(len, 2)
    for family_words in basis
        family_count = length(family_words) ÷ len
        identity_cross_keys = Union{Nothing,K}[nothing for _ in 1:family_count]
        for family in 1:family_count
            word = family_words[len * (family - 1) + 1]
            key, _ = _qmb_source_reduce_key(qmb, K, M, word, len)
            identity_cross_keys[family] = key
        end

        diagonal_terms = [
            Vector{Tuple{Union{Nothing,K},Float64}}(undef, half)
            for _ in 1:family_count
        ]
        for family in 1:family_count
            left = family_words[len * (family - 1) + 1]
            for shift in 1:half
                right = family_words[len * (family - 1) + shift + 1]
                diagonal_terms[family][shift] =
                    _qmb_source_reduce_key(qmb, K, M, UInt16[left; right], len; realify=true)
            end
        end

        offdiagonal_terms = Dict{Tuple{Int,Int},Vector{Tuple{Union{Nothing,K},Float64}}}()
        for left_family in 1:(family_count - 1), right_family in (left_family + 1):family_count
            left = family_words[len * (left_family - 1) + 1]
            terms = Vector{Tuple{Union{Nothing,K},Float64}}(undef, len)
            for shift in 1:len
                right = family_words[len * (right_family - 1) + shift]
                terms[shift] =
                    _qmb_source_reduce_key(qmb, K, M, UInt16[left; right], len; realify=true)
            end
            offdiagonal_terms[(left_family, right_family)] = terms
        end

        push!(
            tables,
            (;
                family_count,
                identity_cross_keys,
                diagonal_terms,
                offdiagonal_terms,
            ),
        )
    end
    return tables
end

function _qmb_source_base_diagonal_terms!(
    entries::Matrix{Vector{Pair{K,C}}},
    diagonal_terms::Vector{Tuple{Union{Nothing,K},Float64}},
    momentum::Int,
    row::Int,
    ::Type{K},
    ::Type{M},
    ::Type{C},
) where {K,C,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    identity_key = _qmb_source_identity_key(K, M)
    _push_qmb_source_entry!(entries, row, row, identity_key, one(C), C)

    half = length(diagonal_terms)
    for shift in 1:half
        key, coef = diagonal_terms[shift]
        key === nothing && continue
        factor = shift == half ?
            coef * (-1.0)^momentum :
            2 * coef * cos(pi * shift * momentum / half)
        _push_qmb_source_entry!(entries, row, row, key, factor, C)
    end
    return nothing
end

function _qmb_source_base_upper_forms(
    record,
    source_table,
    ::Type{K},
    ::Type{M},
    ::Type{C},
) where {K,C<:Real,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    family_count = source_table.family_count
    entries = [Pair{K,C}[] for _ in 1:record.block_size, _ in 1:record.block_size]
    momentum = Int(record.momentum)

    if record.identity_row
        identity_key = _qmb_source_identity_key(K, M)
        _push_qmb_source_entry!(entries, 1, 1, identity_key, one(C), C)
        for family in 1:family_count
            key = source_table.identity_cross_keys[family]
            key === nothing && continue
            _push_qmb_source_entry!(
                entries,
                1,
                family + 1,
                key,
                2 * sqrt(2 * length(source_table.diagonal_terms[family])),
                C,
            )
        end
    end

    offset = record.identity_row ? 1 : 0
    for family in 1:family_count
        row = offset + family
        if record.realified
            _qmb_source_base_diagonal_terms!(
                entries,
                source_table.diagonal_terms[family],
                momentum,
                row,
                K,
                M,
                C,
            )
            _qmb_source_base_diagonal_terms!(
                entries,
                source_table.diagonal_terms[family],
                momentum,
                family_count + row,
                K,
                M,
                C,
            )
        else
            _qmb_source_base_diagonal_terms!(
                entries,
                source_table.diagonal_terms[family],
                momentum,
                row,
                K,
                M,
                C,
            )
        end
    end

    for left_family in 1:(family_count - 1), right_family in (left_family + 1):family_count
        terms = source_table.offdiagonal_terms[(left_family, right_family)]
        len = length(terms)
        for shift in 1:len
            key, coef = terms[shift]
            key === nothing && continue
            angle = 2pi * (shift - 1) * momentum / len
            real_coef = 2 * coef * cos(angle)
            if record.realified
                _push_qmb_source_entry!(
                    entries,
                    left_family,
                    right_family,
                    key,
                    real_coef,
                    C,
                )
                _push_qmb_source_entry!(
                    entries,
                    family_count + left_family,
                    family_count + right_family,
                    key,
                    real_coef,
                    C,
                )
                skew_coef = -2 * coef * sin(angle)
                _push_qmb_source_entry!(
                    entries,
                    family_count + left_family,
                    right_family,
                    key,
                    skew_coef,
                    C,
                )
                _push_qmb_source_entry!(
                    entries,
                    family_count + right_family,
                    left_family,
                    key,
                    -skew_coef,
                    C,
                )
            else
                _push_qmb_source_entry!(
                    entries,
                    offset + left_family,
                    offset + right_family,
                    key,
                    real_coef,
                    C,
                )
            end
        end
    end

    forms = Matrix{NCTSSoS.LinearMomentForm{K,C}}(undef, size(entries))
    for idx in eachindex(entries)
        forms[idx] = _qmb_source_form_from_pairs!(entries[idx])
    end
    return forms
end

function _qmb_source_block_label(block)
    origin = block.meta.origin
    return hasproperty(origin, :label) ? origin.label : nothing
end

function _qmb_source_base_blocks(linear)
    return [
        block for block in linear.psd_blocks_lin
        if (_qmb_source_block_label(block) isa NamedTuple) &&
            haskey(_qmb_source_block_label(block), :feature) &&
            _qmb_source_block_label(block).feature == :qmbcertify_base
    ]
end

function _qmb_source_label_matches_record(label, record)
    label isa NamedTuple || return false
    return haskey(label, :parity) &&
        haskey(label, :momentum) &&
        haskey(label, :identity_row) &&
        haskey(label, :realified) &&
        label.parity == record.parity &&
        label.momentum == record.momentum &&
        label.identity_row == record.identity_row &&
        label.realified == record.realified
end

function _qmb_source_symmetric_entry_form(
    left::NCTSSoS.LinearMomentForm{K,C},
    right::NCTSSoS.LinearMomentForm{K,C},
) where {K,C}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(left) + length(right))
    append!(pairs, left.terms)
    append!(pairs, right.terms)
    return NCTSSoS._linear_moment_form_from_owned_pairs!(pairs)
end

function _qmb_source_symmetric_entry_form(block, row::Int, col::Int)
    row == col && return block.entries[row, col]
    return _qmb_source_symmetric_entry_form(
        block.entries[row, col],
        block.entries[col, row],
    )
end

function _qmb_source_form_delta(source_form, nctssos_form; atol::Real)
    if length(source_form) != length(nctssos_form)
        return false, Inf
    end
    max_abs_delta = 0.0
    for (source_term, nctssos_term) in zip(source_form, nctssos_form)
        if !NCTSSoS.key_isequal(source_term.first, nctssos_term.first)
            return false, Inf
        end
        max_abs_delta = max(max_abs_delta, abs(source_term.second - nctssos_term.second))
    end
    return max_abs_delta <= atol, max_abs_delta
end

function _qmb_source_base_entry_skipped(reason::AbstractString)
    return (
        skipped=true,
        skip_reason=String(reason),
        block_count=0,
        block_label_match=false,
        entry_match=false,
        entry_mismatch_count=0,
        entry_max_abs_delta=0.0,
        entry_compared_count=0,
        mismatch_sample=Any[],
    )
end

function _qmbcertify_base_source_entry_metrics(
    qmb::Module,
    basis::Vector{Vector{Vector{UInt16}}},
    len::Int,
    linear::NCTSSoS.MomentLinearData{K,C,M};
    atol::Real=1e-10,
) where {K,C,M}
    isempty(linear.psd_blocks_lin) && return _qmb_source_base_entry_skipped("no_psd_blocks")

    C <: Real || return _qmb_source_base_entry_skipped("nonreal_linear_coefficients")

    family_counts = [length(words) ÷ len for words in basis]
    records = NCTSSoS._qmbcertify_chain_base_block_records(family_counts, len)
    source_tables = _qmb_source_base_tables(qmb, basis, len, K, M)
    blocks = _qmb_source_base_blocks(linear)
    label_match = length(blocks) == length(records)
    mismatch_count = 0
    compared_count = 0
    max_abs_delta = 0.0
    sample = Any[]

    for (block_idx, record) in enumerate(records)
        if block_idx > length(blocks)
            mismatch_count += 1
            length(sample) < 5 && push!(sample, (block=block_idx, reason=:missing_block))
            continue
        end
        block = blocks[block_idx]
        label = _qmb_source_block_label(block)
        if !_qmb_source_label_matches_record(label, record)
            label_match = false
            mismatch_count += 1
            length(sample) < 5 && push!(
                sample,
                (block=block_idx, reason=:label_mismatch, label=label, record=record),
            )
            continue
        end
        if block.size != record.block_size
            mismatch_count += 1
            length(sample) < 5 && push!(
                sample,
                (block=block_idx, reason=:size_mismatch, source=record.block_size, nctssos=block.size),
            )
            continue
        end

        source_forms = _qmb_source_base_upper_forms(
            record,
            source_tables[record.parity + 1],
            K,
            M,
            C,
        )
        for col in 1:block.size, row in 1:col
            source_form = source_forms[row, col]
            nctssos_form = _qmb_source_symmetric_entry_form(block, row, col)
            ok, delta = _qmb_source_form_delta(source_form, nctssos_form; atol)
            compared_count += 1
            max_abs_delta = max(max_abs_delta, delta)
            ok && continue
            mismatch_count += 1
            if length(sample) < 5
                push!(
                    sample,
                    (
                        block=block_idx,
                        row,
                        col,
                        source_terms=length(source_form),
                        nctssos_terms=length(nctssos_form),
                        delta,
                    ),
                )
            end
        end
    end

    return (
        skipped=false,
        skip_reason="",
        block_count=length(blocks),
        block_label_match=label_match,
        entry_match=label_match && mismatch_count == 0 && max_abs_delta <= atol,
        entry_mismatch_count=mismatch_count,
        entry_max_abs_delta=max_abs_delta,
        entry_compared_count=compared_count,
        mismatch_sample=sample,
    )
end

function _try_nctssos_qmb_basis_formulation(pop, ops, d::Int, rdm, qmb_basis)
    try
        linear, report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            d;
            basis=qmb_basis,
            reflection_symmetry=true,
            sign_symmetry=true,
            real_moment_matrix=true,
            axis_rotation_symmetry=true,
            axis_rotation_equalities=true,
            axis_rotation_quotient=true,
            contiguous_rdm_k=rdm,
            contiguous_rdm_decomposition=:qmbcertify,
            contiguous_rdm_support=:extend,
            linear_state_opt_width=7,
            linear_state_opt_mode=:qmbcertify,
        )
        metrics = translation_report_metrics(report)
        return (
            ok=true,
            moments=length(linear.moments),
            zero_constraints=length(linear.zero_constraints),
            lso_rows=metrics.linear_state_opt_row_count,
            psd_block_histogram=_audit_histogram_pairs(report.psd_block_sizes),
            max_block=maximum(report.psd_block_sizes),
        )
    catch err
        return (
            ok=false,
            error_type=string(typeof(err)),
            error_message=sprint(showerror, err),
        )
    end
end

function _profile_bool_arg(args, key::AbstractString, default::Bool)
    raw = get(args, key, default)
    raw isa Bool && return raw
    raw isa AbstractString && lowercase(raw) in ("true", "1", "yes") && return true
    raw isa AbstractString && lowercase(raw) in ("false", "0", "no") && return false
    throw(ArgumentError("Unsupported boolean profile argument $key=$(repr(raw))."))
end

function _try_nctssos_qmb_source_base_formulation(
    pop,
    ops,
    d::Int,
    rdm,
    args;
    qmb=nothing,
    qmb_basis_source=nothing,
)
    try
        extra = Int(get(args, "extra", 0))
        pso = Int(get(args, "pso", 0))
        lso = _profile_bool_arg(args, "lso", true)
        linear_state_width = lso ? 7 : nothing
        psd_state_width = pso > 0 ? pso : nothing
        linear, report = NCTSSoS._pauli_qmbcertify_chain_base_linear_relaxation(
            pop,
            ops,
            d;
            extra,
            three_type=(1, 1),
            real_moment_matrix=true,
            contiguous_rdm_k=rdm,
            contiguous_rdm_decomposition=:qmbcertify,
            contiguous_rdm_support=:extend,
            linear_state_opt_width=linear_state_width,
            linear_state_opt_mode=:qmbcertify,
            psd_state_opt_width=psd_state_width,
        )
        n = length(first(ops))
        base_entry_metrics = if qmb === nothing || qmb_basis_source === nothing
            _qmb_source_base_entry_skipped("missing_qmb_source")
        else
            _qmbcertify_base_source_entry_metrics(qmb, qmb_basis_source, n, linear)
        end
        metrics = translation_report_metrics(report)
        source_base_blocks = Int[
            size for (label, size) in zip(report.block_labels, report.psd_block_sizes)
            if label isa NamedTuple && haskey(label, :feature) &&
                label.feature == :qmbcertify_base
        ]
        return (
            ok=true,
            moments=length(linear.moments),
            free_keys=length(linear.free_keys),
            zero_constraints=length(linear.zero_constraints),
            lso_rows=metrics.linear_state_opt_row_count,
            psd_state_opt_blocks=metrics.psd_state_opt_block_count,
            rdm_blocks=metrics.contiguous_rdm_psd_block_sizes,
            source_base_block_histogram=_audit_histogram_pairs(source_base_blocks),
            psd_block_histogram=_audit_histogram_pairs(report.psd_block_sizes),
            max_block=maximum(report.psd_block_sizes),
            solve_supported=metrics.solve_supported,
            solve_blocker=metrics.solve_blocker,
            base_entry_metrics,
        )
    catch err
        return (
            ok=false,
            error_type=string(typeof(err)),
            error_message=sprint(showerror, err),
        )
    end
end

function _print_nctssos_source_base_metrics(metrics)
    println("nctssos_source_base_ok=$(metrics.ok)")
    if metrics.ok
        println("nctssos_source_base_linear_moments=$(metrics.moments)")
        println("nctssos_source_base_free_keys=$(metrics.free_keys)")
        println("nctssos_source_base_zero_constraints=$(metrics.zero_constraints)")
        println("nctssos_source_base_lso_row_count=$(metrics.lso_rows)")
        println("nctssos_source_base_psd_state_opt_block_count=$(metrics.psd_state_opt_blocks)")
        println("nctssos_source_base_rdm_blocks=$(metrics.rdm_blocks)")
        println("nctssos_source_base_block_histogram=$(metrics.source_base_block_histogram)")
        println("nctssos_source_base_full_psd_block_histogram=$(metrics.psd_block_histogram)")
        println("nctssos_source_base_psd_max_block=$(metrics.max_block)")
        println("nctssos_source_base_solve_supported=$(metrics.solve_supported)")
        println("nctssos_source_base_solve_blocker=$(metrics.solve_blocker)")
        entry = metrics.base_entry_metrics
        println("nctssos_source_base_entry_skipped=$(entry.skipped)")
        println("nctssos_source_base_entry_skip_reason=$(entry.skip_reason)")
        println("nctssos_source_base_entry_block_count=$(entry.block_count)")
        println("nctssos_source_base_entry_block_label_match=$(entry.block_label_match)")
        println("nctssos_source_base_entry_compared_count=$(entry.entry_compared_count)")
        println("nctssos_source_base_entry_match=$(entry.entry_match)")
        println("nctssos_source_base_entry_mismatch_count=$(entry.entry_mismatch_count)")
        println("nctssos_source_base_entry_max_abs_delta=$(entry.entry_max_abs_delta)")
        println("nctssos_source_base_entry_mismatch_sample=$(entry.mismatch_sample)")
    else
        println("nctssos_source_base_error_type=$(metrics.error_type)")
        println("nctssos_source_base_error_message=$(metrics.error_message)")
    end
    return nothing
end

function _audit_one(
    len::Int,
    args,
    qmb::Module;
    qmb_path::AbstractString,
    try_qmb_basis::Bool,
    qmb_basis_probe_max_l::Int,
    compare_direct::Bool,
    rdm_entry_max_k::Int,
)
    d = Int(args["d"])
    rdm = _rdm_value(get(args, "rdm", false))
    lol = _lol_value(get(args, "lol", "Int(L/2)"), len)

    basis, tsupp = _qmb_chain_basis_and_tsupp(qmb, len, args)
    qmb_base_blocks = _qmb_base_block_sizes(basis, len)
    qmb_rdm_blocks = _rdm_block_sizes(rdm)
    generated = Base.invokelatest(qmb.generate_mons, len, lol, Int(rdm) - 1)
    filtered = Base.invokelatest(
        qmb.filter_mons,
        generated,
        tsupp,
        [1],
        [0.75],
        len;
        lattice="chain",
    )
    qmb_lso_keys = Set{String}()
    for mon in filtered
        key = _qmb_lso_row_key(qmb, Vector{UInt16}(mon), len)
        isempty(key) || push!(qmb_lso_keys, key)
    end

    registry, ops = create_pauli_variables(1:len)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    M = eltype(first(ops))
    qmb_family_records = _qmb_basis_family_records(basis, len, M)
    qmb_family_reps = Set(record.representative for record in qmb_family_records)
    qmb_basis = _qmb_basis_as_nctssos(basis, len, M)
    contiguous_basis = pauli_contiguous_chain_basis(ops, d)
    contiguous_reps = Set(NCTSSoS._translation_orbit_representatives(contiguous_basis, len))
    qmb_not_contiguous = [
        record for record in qmb_family_records if !(record.representative in contiguous_reps)
    ]
    contiguous_not_qmb = [
        rep for rep in sort!(collect(contiguous_reps)) if !(rep in qmb_family_reps)
    ]
    qmb_short_orbit = [
        record for record in qmb_family_records if record.orbit_length != len
    ]
    ncts_tests = NCTSSoS._contiguous_state_opt_tests(M, len, 7; sign_symmetry=true)
    ncts_qmbcertify_lso_tests = NCTSSoS._linear_state_opt_tests(
        M,
        len,
        7;
        sign_symmetry=true,
        mode=:qmbcertify,
    )
    ncts_lso_keys = Set{String}()
    for test_mono in ncts_tests
        key = _nctssos_lso_row_key(qmb, test_mono, len)
        isempty(key) || push!(ncts_lso_keys, key)
    end
    ncts_qmbcertify_lso_keys = Set{String}()
    for test_mono in ncts_qmbcertify_lso_tests
        key = _nctssos_lso_row_key(qmb, test_mono, len)
        isempty(key) || push!(ncts_qmbcertify_lso_keys, key)
    end

    linear = nothing
    report = nothing
    metrics = nothing
    ncts_base_block_sizes = Int[]
    if compare_direct
        linear, report = NCTSSoS._pauli_translation_base_linear_relaxation(
            pop,
            ops,
            d;
            reflection_symmetry=true,
            sign_symmetry=true,
            real_moment_matrix=true,
            axis_rotation_symmetry=true,
            axis_rotation_equalities=true,
            axis_rotation_quotient=true,
            contiguous_rdm_k=rdm,
            contiguous_rdm_decomposition=:qmbcertify,
            contiguous_rdm_support=:extend,
            linear_state_opt_width=7,
        )
        metrics = translation_report_metrics(report)
        ncts_base_block_sizes = Int[
            size for (label, size) in zip(report.block_labels, report.psd_block_sizes)
            if !(label isa NamedTuple && haskey(label, :feature))
        ]
    end
    source_base_metrics = _try_nctssos_qmb_source_base_formulation(
        pop,
        ops,
        d,
        rdm,
        args,
        qmb=qmb,
        qmb_basis_source=basis,
    )
    objective_metrics = _nctssos_qmb_objective_metrics(pop, args, len, M)
    rdm_source_metrics = _qmbcertify_rdm_source_entry_metrics(
        qmb_path,
        len,
        rdm,
        M;
        max_k=rdm_entry_max_k,
    )

    intersection = intersect(qmb_lso_keys, ncts_lso_keys)
    qmb_only = setdiff(qmb_lso_keys, ncts_lso_keys)
    ncts_only = setdiff(ncts_lso_keys, qmb_lso_keys)
    qmbcertify_lso_intersection = intersect(qmb_lso_keys, ncts_qmbcertify_lso_keys)
    qmbcertify_lso_qmb_only = setdiff(qmb_lso_keys, ncts_qmbcertify_lso_keys)
    qmbcertify_lso_ncts_only = setdiff(ncts_qmbcertify_lso_keys, qmb_lso_keys)

    println("\n## L=$len")
    println("qmb_basis_lengths=$(length.(basis))")
    println("qmb_basis_k=$([length(basis[i]) ÷ len for i in 1:2])")
    println("qmb_basis_family_count=$(length(qmb_family_records))")
    println("qmb_basis_unique_nctssos_count=$(length(qmb_basis))")
    _print_pair_hist("qmb_basis_family_word_length_histogram", Int[record.word_length for record in qmb_family_records])
    _print_pair_hist("qmb_basis_family_translation_orbit_length_histogram", Int[record.orbit_length for record in qmb_family_records])
    println("qmb_basis_short_translation_orbit_family_count=$(length(qmb_short_orbit))")
    println("qmb_basis_short_translation_orbit_sample=$(_audit_records_sample(qmb_short_orbit))")
    println("nctssos_contiguous_order_representative_count=$(length(contiguous_reps))")
    println("qmb_family_reps_in_nctssos_contiguous_order_count=$(length(qmb_family_records) - length(qmb_not_contiguous))")
    println("qmb_family_reps_not_in_nctssos_contiguous_order_count=$(length(qmb_not_contiguous))")
    println("qmb_family_reps_not_in_nctssos_contiguous_order_sample=$(_audit_records_sample(qmb_not_contiguous))")
    println("nctssos_contiguous_order_reps_not_in_qmb_family_count=$(length(contiguous_not_qmb))")
    println("nctssos_contiguous_order_reps_not_in_qmb_family_sample=$(_set_sample(Set(sprint(show, rep) for rep in contiguous_not_qmb)))")
    println("nctssos_qmb_objective_total_match=$(objective_metrics.match)")
    println("nctssos_qmb_objective_nctssos_term_count=$(objective_metrics.nctssos_term_count)")
    println("nctssos_qmb_objective_qmb_total_term_count=$(objective_metrics.qmb_total_term_count)")
    println("nctssos_qmb_objective_max_abs_delta=$(objective_metrics.max_abs_delta)")
    println("nctssos_qmb_objective_nctssos_only_count=$(objective_metrics.nctssos_only_count)")
    println("nctssos_qmb_objective_qmb_only_count=$(objective_metrics.qmb_only_count)")
    println("nctssos_qmb_objective_nctssos_only_sample=$(objective_metrics.nctssos_only_sample)")
    println("nctssos_qmb_objective_qmb_only_sample=$(objective_metrics.qmb_only_sample)")
    println("qmb_tsupp_count=$(length(tsupp))")
    _print_pair_hist("qmb_base_block_histogram", qmb_base_blocks)
    _print_pair_hist("qmb_rdm_block_histogram", qmb_rdm_blocks)
    println("qmb_lso_generated_count=$(length(generated))")
    println("qmb_lso_filtered_count=$(length(filtered))")
    println("qmb_lso_distinct_row_count=$(length(qmb_lso_keys))")
    println("qmb_rdm_source_k=$(rdm_source_metrics.k)")
    println("qmb_rdm_source_row_block_exact_match=$(rdm_source_metrics.row_block_exact_match)")
    println("qmb_rdm_source_block_sizes=$(rdm_source_metrics.source_block_sizes)")
    println("nctssos_qmb_rdm_block_sizes=$(rdm_source_metrics.nctssos_block_sizes)")
    println("qmb_rdm_source_entry_skipped=$(rdm_source_metrics.skipped)")
    println("qmb_rdm_source_entry_skip_reason=$(rdm_source_metrics.skip_reason)")
    println("qmb_rdm_source_entry_scale=$(rdm_source_metrics.scale)")
    println("qmb_rdm_source_entry_compared_count=$(rdm_source_metrics.entry_compared_count)")
    println("qmb_rdm_source_entry_match=$(rdm_source_metrics.entry_match)")
    println("qmb_rdm_source_entry_mismatch_count=$(rdm_source_metrics.entry_mismatch_count)")
    println("qmb_rdm_source_entry_max_abs_delta=$(rdm_source_metrics.entry_max_abs_delta)")
    println("nctssos_lso_test_count=$(length(ncts_tests))")
    println("nctssos_lso_distinct_qmb_reduced_row_count=$(length(ncts_lso_keys))")
    println("lso_row_intersection_count=$(length(intersection))")
    println("lso_row_qmb_only_count=$(length(qmb_only))")
    println("lso_row_nctssos_only_count=$(length(ncts_only))")
    println("lso_row_qmb_only_sample=$(_set_sample(qmb_only))")
    println("lso_row_nctssos_only_sample=$(_set_sample(ncts_only))")
    println("nctssos_qmbcertify_lso_test_count=$(length(ncts_qmbcertify_lso_tests))")
    println("nctssos_qmbcertify_lso_distinct_row_count=$(length(ncts_qmbcertify_lso_keys))")
    println("qmbcertify_lso_row_intersection_count=$(length(qmbcertify_lso_intersection))")
    println("qmbcertify_lso_row_qmb_only_count=$(length(qmbcertify_lso_qmb_only))")
    println("qmbcertify_lso_row_nctssos_only_count=$(length(qmbcertify_lso_ncts_only))")
    println("qmbcertify_lso_row_qmb_only_sample=$(_set_sample(qmbcertify_lso_qmb_only))")
    println("qmbcertify_lso_row_nctssos_only_sample=$(_set_sample(qmbcertify_lso_ncts_only))")
    println("nctssos_direct_compare=$(compare_direct)")
    if compare_direct
        println("nctssos_linear_moments=$(length(linear.moments))")
        println("nctssos_free_keys=$(length(linear.free_keys))")
        println("nctssos_zero_constraints=$(length(linear.zero_constraints))")
        println("nctssos_lso_report_row_count=$(metrics.linear_state_opt_row_count)")
        _print_pair_hist("nctssos_base_psd_block_histogram", ncts_base_block_sizes)
        _print_pair_hist("nctssos_full_psd_block_histogram", report.psd_block_sizes)
        println("nctssos_psd_max_block=$(maximum(report.psd_block_sizes))")
        println("nctssos_report_rdm_blocks=$(metrics.contiguous_rdm_psd_block_sizes)")
        println("nctssos_report_axis_quotient=$(metrics.axis_rotation_quotient)")
    else
        println("nctssos_linear_moments=skipped")
        println("nctssos_free_keys=skipped")
        println("nctssos_zero_constraints=skipped")
        println("nctssos_lso_report_row_count=skipped")
        println("nctssos_base_psd_block_histogram=skipped")
        println("nctssos_full_psd_block_histogram=skipped")
        println("nctssos_psd_max_block=skipped")
        println("nctssos_report_rdm_blocks=skipped")
        println("nctssos_report_axis_quotient=skipped")
    end
    _print_nctssos_source_base_metrics(source_base_metrics)
    if compare_direct && try_qmb_basis && len <= qmb_basis_probe_max_l
        println("nctssos_qmb_basis_probe=$(_try_nctssos_qmb_basis_formulation(pop, ops, d, rdm, qmb_basis))")
    else
        println("nctssos_qmb_basis_probe=(skipped=true, max_l=$qmb_basis_probe_max_l)")
    end
    return nothing
end

function main()
    profile = String(get(ENV, "NCTS_FORMULATION_AUDIT_PROFILE", "A1"))
    args = _profile_args(profile)
    source = _profile_source(profile)
    qmb_path = String(get(ENV, "NCTS_QMBCERTIFY_PATH", _default_qmbcertify_path(source)))
    lengths = _audit_lengths()
    try_qmb_basis = _parse_bool_env("NCTS_FORMULATION_AUDIT_TRY_QMB_BASIS", true)
    qmb_basis_probe_max_l = _parse_int_env("NCTS_FORMULATION_AUDIT_QMB_BASIS_MAX_L", 20)
    compare_direct = _parse_bool_env("NCTS_FORMULATION_AUDIT_DIRECT_COMPARE", true)
    rdm_entry_max_k = _parse_int_env("NCTS_FORMULATION_AUDIT_RDM_ENTRY_MAX_K", 8)
    _check_size_guard(lengths)

    println("# QMBCertify/NCTSSoS formulation audit")
    println("profile=$profile")
    println("lengths=$(join(lengths, ","))")
    println("qmbcertify_path=$qmb_path")
    println("solve=false")
    println("dualize=false")
    println("try_qmb_basis=$try_qmb_basis")
    println("qmb_basis_probe_max_l=$qmb_basis_probe_max_l")
    println("compare_direct=$compare_direct")
    println("rdm_entry_max_k=$rdm_entry_max_k")

    qmb = _load_qmbcertify(qmb_path, source)
    for len in lengths
        _audit_one(
            len,
            args,
            qmb;
            qmb_path,
            try_qmb_basis,
            qmb_basis_probe_max_l,
            compare_direct,
            rdm_entry_max_k,
        )
        flush(stdout)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
