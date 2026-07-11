# SU(2)-invariant Pauli moment-coordinate quotient.

const _PAULI_SU2_SINGLET_COUNTS = (1, 0, 1, 1, 3, 6, 15, 36, 91)
const _PAULI_SU2_SPIN1_CG_CACHE = Dict{NTuple{4,Int},Float64}()
const _PAULI_SU2_SPIN1_CG_CACHE_LOCK = ReentrantLock()

function _pauli_su2_factorial(n::Int)
    n >= 0 || throw(DomainError(n, "factorial argument must be non-negative"))
    return factorial(big(n))
end

function _pauli_su2_spin1_cg_uncached(
    parent_spin2::Int,
    parent_m2::Int,
    local_m2::Int,
    child_spin2::Int,
)
    all(iseven, (parent_spin2, parent_m2, local_m2, child_spin2)) ||
        throw(ArgumentError("spin-1 coupling requires even doubled quantum numbers"))
    local_m2 in (-2, 0, 2) ||
        throw(DomainError(local_m2, "local m2 must be -2, 0, or 2"))

    j1 = parent_spin2 ÷ 2
    m1 = parent_m2 ÷ 2
    j2 = 1
    m2 = local_m2 ÷ 2
    j = child_spin2 ÷ 2
    m = m1 + m2
    abs(m1) <= j1 && abs(m2) <= j2 && abs(m) <= j || return 0.0
    abs(j1 - j2) <= j <= j1 + j2 || return 0.0

    triangle = (
        BigInt(2j + 1) * _pauli_su2_factorial(j + j1 - j2) *
        _pauli_su2_factorial(j - j1 + j2) *
        _pauli_su2_factorial(j1 + j2 - j)
    ) // _pauli_su2_factorial(j1 + j2 + j + 1)
    magnetic = (
        _pauli_su2_factorial(j + m) * _pauli_su2_factorial(j - m) *
        _pauli_su2_factorial(j1 - m1) * _pauli_su2_factorial(j1 + m1) *
        _pauli_su2_factorial(j2 - m2) * _pauli_su2_factorial(j2 + m2)
    ) // one(BigInt)

    kmin = maximum((0, j2 - j - m1, j1 + m2 - j))
    kmax = minimum((j1 + j2 - j, j1 - m1, j2 + m2))
    kmin <= kmax || return 0.0
    series = zero(Rational{BigInt})
    for k in kmin:kmax
        denominator = _pauli_su2_factorial(k) *
            _pauli_su2_factorial(j1 + j2 - j - k) *
            _pauli_su2_factorial(j1 - m1 - k) *
            _pauli_su2_factorial(j2 + m2 - k) *
            _pauli_su2_factorial(j - j2 + m1 + k) *
            _pauli_su2_factorial(j - j1 - m2 + k)
        numerator = isodd(k) ? -one(BigInt) : one(BigInt)
        series += numerator // denominator
    end
    value = sqrt(Float64(triangle * magnetic)) * Float64(series)
    return iszero(value) ? 0.0 : value
end

function _pauli_su2_spin1_cg(
    parent_spin2::Int,
    parent_m2::Int,
    local_m2::Int,
    child_spin2::Int,
)
    key = (parent_spin2, parent_m2, local_m2, child_spin2)
    return lock(_PAULI_SU2_SPIN1_CG_CACHE_LOCK) do
        get!(_PAULI_SU2_SPIN1_CG_CACHE, key) do
            _pauli_su2_spin1_cg_uncached(key...)
        end
    end
end

function _pauli_su2_singlet_spin_paths(support_size::Int)
    support_size >= 0 ||
        throw(DomainError(support_size, "support size must be non-negative"))
    paths = Vector{Vector{Int}}()

    function visit!(path::Vector{Int}, remaining::Int)
        current = last(path)
        if iszero(remaining)
            iszero(current) && push!(paths, copy(path))
            return
        end
        for child in abs(current - 2):2:(current + 2)
            child <= 2 * (remaining - 1) || continue
            push!(path, child)
            visit!(path, remaining - 1)
            pop!(path)
        end
    end

    visit!(Int[0], support_size)
    return paths
end

function _pauli_su2_singlet_cartesian_row!(
    row::Vector{Tuple{Int,ComplexF64}},
    path::Vector{Int},
    support_size::Int,
    local_transform::Matrix{ComplexF64},
    atol::Float64,
)
    max_state_count = 2 * support_size + 1
    amplitudes = zeros(ComplexF64, max_state_count)
    next_amplitudes = similar(amplitudes)

    for code0 in 0:(3^support_size - 1)
        Base.fill!(amplitudes, 0)
        amplitudes[1] = 1
        code = code0
        for site in 1:support_size
            Base.fill!(next_amplitudes, 0)
            cartesian = mod(code, 3) + 1
            code = div(code, 3)
            parent_spin2 = path[site]
            child_spin2 = path[site + 1]
            for parent_idx in 1:(parent_spin2 + 1)
                amplitude = amplitudes[parent_idx]
                iszero(amplitude) && continue
                parent_m2 = -parent_spin2 + 2 * (parent_idx - 1)
                for (local_idx, local_m2) in enumerate((-2, 0, 2))
                    child_m2 = parent_m2 + local_m2
                    abs(child_m2) <= child_spin2 || continue
                    cg = _pauli_su2_spin1_cg(
                        parent_spin2,
                        parent_m2,
                        local_m2,
                        child_spin2,
                    )
                    iszero(cg) && continue
                    child_idx = (child_m2 + child_spin2) ÷ 2 + 1
                    next_amplitudes[child_idx] +=
                        amplitude * cg * local_transform[local_idx, cartesian]
                end
            end
            amplitudes, next_amplitudes = next_amplitudes, amplitudes
        end
        value = amplitudes[1]
        abs(value) <= atol || push!(row, (code0 + 1, value))
    end
    if !isempty(row)
        first_value = first(row)[2]
        phase = conj(first_value) / abs(first_value)
        for idx in eachindex(row)
            column, coefficient = row[idx]
            row[idx] = (column, coefficient * phase)
        end
    end
    return row
end

function _pauli_su2_word_singlet_rows(
    support_size::Integer;
    atol::Real=1e-13,
)
    s = Int(support_size)
    s >= 0 || throw(DomainError(support_size, "support size must be non-negative"))
    s <= 8 || throw(ArgumentError(
        "moment quotient supports active support through 8; got $s",
    ))
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "atol must be non-negative"))

    local_transform = _pauli_su2_word_local_spherical_transform()
    rows = Vector{Vector{Tuple{Int,ComplexF64}}}()
    for path in _pauli_su2_singlet_spin_paths(s)
        row = Tuple{Int,ComplexF64}[]
        _pauli_su2_singlet_cartesian_row!(
            row,
            path,
            s,
            local_transform,
            tolerance,
        )
        push!(rows, row)
    end
    length(rows) == _PAULI_SU2_SINGLET_COUNTS[s + 1] || throw(ArgumentError(
        "singlet path count mismatch for support $s: got $(length(rows))",
    ))
    return _sparse_transform_rows(rows, 3^s)
end

function _dense_sparse_transform_rows(transform::_SparseTransformRows{C}) where {C}
    dense = zeros(C, size(transform))
    for (row_idx, row) in pairs(transform.rows)
        for (column, coefficient) in row
            dense[row_idx, column] = coefficient
        end
    end
    return dense
end

const _PAULI_SU2_INTEGER_INVARIANT_ROWS_CACHE =
    Dict{Int,Vector{Dict{Int,Int}}}()
const _PAULI_SU2_INTEGER_INVARIANT_ROWS_CACHE_LOCK = ReentrantLock()
const _PAULI_SU2_PROJECTED_RANK_CACHE = Dict{Any,Int}()
const _PAULI_SU2_PROJECTED_RANK_CACHE_LOCK = ReentrantLock()

function _pauli_su2_index_pairings(indices::Vector{Int})
    isempty(indices) && return [NTuple{2,Int}[]]
    iseven(length(indices)) || throw(ArgumentError(
        "pairing index count must be even, got $(length(indices))",
    ))
    first_index = first(indices)
    pairings = Vector{Vector{NTuple{2,Int}}}()
    for partner_position in 2:length(indices)
        partner = indices[partner_position]
        remaining = Int[
            indices[idx] for idx in 2:length(indices) if idx != partner_position
        ]
        for suffix in _pauli_su2_index_pairings(remaining)
            push!(pairings, vcat(NTuple{2,Int}[(first_index, partner)], suffix))
        end
    end
    return pairings
end

const _PAULI_SU2_EPSILON_ASSIGNMENTS = (
    ((0, 1, 2), 1),
    ((1, 2, 0), 1),
    ((2, 0, 1), 1),
    ((0, 2, 1), -1),
    ((2, 1, 0), -1),
    ((1, 0, 2), -1),
)

function _pauli_su2_integer_invariant_row(
    support_size::Int,
    pairing::Vector{NTuple{2,Int}};
    epsilon_positions::Union{Nothing,NTuple{3,Int}}=nothing,
)
    axes = zeros(Int, support_size)
    row = Dict{Int,Int}()

    function emit!()
        if epsilon_positions === nothing
            code0 = _pauli_su2_axes_code(Tuple(axes))
            row[code0] = get(row, code0, 0) + 1
            return
        end
        for (epsilon_axes, sign) in _PAULI_SU2_EPSILON_ASSIGNMENTS
            for idx in 1:3
                axes[epsilon_positions[idx]] = epsilon_axes[idx]
            end
            code0 = _pauli_su2_axes_code(Tuple(axes))
            row[code0] = get(row, code0, 0) + sign
        end
    end

    function assign_pair!(pair_idx::Int)
        if pair_idx > length(pairing)
            emit!()
            return
        end
        left, right = pairing[pair_idx]
        for axis in 0:2
            axes[left] = axis
            axes[right] = axis
            assign_pair!(pair_idx + 1)
        end
    end

    assign_pair!(1)
    filter!(pair -> !iszero(last(pair)), row)
    return row
end

function _pauli_su2_integer_invariant_rows_uncached(support_size::Int)
    support_size >= 0 || throw(DomainError(
        support_size,
        "support size must be non-negative",
    ))
    rows = Dict{Int,Int}[]
    indices = collect(1:support_size)
    if iseven(support_size)
        for pairing in _pauli_su2_index_pairings(indices)
            push!(rows, _pauli_su2_integer_invariant_row(support_size, pairing))
        end
    elseif support_size >= 3
        for first_idx in 1:(support_size - 2)
            for second_idx in (first_idx + 1):(support_size - 1)
                for third_idx in (second_idx + 1):support_size
                    epsilon_positions = (first_idx, second_idx, third_idx)
                    remaining = Int[
                        idx for idx in indices if idx ∉ epsilon_positions
                    ]
                    for pairing in _pauli_su2_index_pairings(remaining)
                        push!(
                            rows,
                            _pauli_su2_integer_invariant_row(
                                support_size,
                                pairing;
                                epsilon_positions,
                            ),
                        )
                    end
                end
            end
        end
    end
    return rows
end

function _pauli_su2_integer_invariant_rows(support_size::Int)
    return lock(_PAULI_SU2_INTEGER_INVARIANT_ROWS_CACHE_LOCK) do
        get!(_PAULI_SU2_INTEGER_INVARIANT_ROWS_CACHE, support_size) do
            _pauli_su2_integer_invariant_rows_uncached(support_size)
        end
    end
end

function _pauli_su2_exact_sparse_integer_rank(rows::Vector{Dict{Int,Int}})
    pivots = Dict{Int,Dict{Int,Rational{BigInt}}}()
    for integer_row in rows
        row = Dict{Int,Rational{BigInt}}(
            column => (BigInt(value) // one(BigInt))
            for (column, value) in integer_row if !iszero(value)
        )
        while !isempty(row)
            pivot_column = minimum(keys(row))
            pivot_row = get(pivots, pivot_column, nothing)
            if pivot_row === nothing
                pivot_value = row[pivot_column]
                for column in collect(keys(row))
                    row[column] /= pivot_value
                end
                pivots[pivot_column] = row
                break
            end
            factor = row[pivot_column]
            for (column, value) in pivot_row
                updated = get(row, column, zero(Rational{BigInt})) - factor * value
                if iszero(updated)
                    delete!(row, column)
                else
                    row[column] = updated
                end
            end
        end
    end
    return length(pivots)
end

function _pauli_su2_exact_projected_singlet_rank(
    support_orbit::Tuple,
    n_sites::Int,
    observed_patterns::Set{Tuple{Vararg{UInt8}}},
)
    _, permutations = _pauli_su2_support_stabilizer_permutations(
        support_orbit,
        n_sites,
    )
    observed_codes = sort!(Int[_pauli_su2_axes_code(axes) for axes in observed_patterns])
    observed_columns = Dict(code0 => idx for (idx, code0) in pairs(observed_codes))
    cache_key = (
        length(support_orbit),
        Tuple(Tuple(permutation) for permutation in permutations),
        Tuple(observed_codes),
    )
    return lock(_PAULI_SU2_PROJECTED_RANK_CACHE_LOCK) do
        get!(_PAULI_SU2_PROJECTED_RANK_CACHE, cache_key) do
            restricted_rows = Dict{Int,Int}[]
            for invariant_row in _pauli_su2_integer_invariant_rows(length(support_orbit))
                symmetrized = Dict{Int,Int}()
                for (code0, coefficient) in invariant_row
                    for permutation in permutations
                        permuted_code0 = _pauli_su2_permute_axes_code(code0, permutation)
                        symmetrized[permuted_code0] =
                            get(symmetrized, permuted_code0, 0) + coefficient
                    end
                end
                restricted = Dict{Int,Int}()
                for (code0, coefficient) in symmetrized
                    column = get(observed_columns, code0, 0)
                    iszero(column) && continue
                    restricted[column] = coefficient
                end
                filter!(pair -> !iszero(last(pair)), restricted)
                isempty(restricted) || push!(restricted_rows, restricted)
            end
            _pauli_su2_exact_sparse_integer_rank(restricted_rows)
        end
    end
end

function _pauli_su2_select_expected_singlet_rows(
    singlet::AbstractMatrix{<:Complex},
    expected_rank::Int;
    atol::Real,
)
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "atol must be non-negative"))
    matrix = ComplexF64.(singlet)
    row_count = size(matrix, 1)
    0 <= expected_rank <= row_count || throw(DimensionMismatch(
        "exact projected singlet rank $expected_rank is outside 0:$row_count",
    ))
    expected_rank == row_count && return matrix
    iszero(expected_rank) && return zeros(ComplexF64, 0, size(matrix, 2))

    factor = qr(adjoint(matrix), ColumnNorm())
    selected_rows = sort!(Int[factor.p[idx] for idx in 1:expected_rank])
    independent = matrix[selected_rows, :]
    row_coefficients = adjoint(independent) \ adjoint(matrix)
    scale = max(maximum(abs, matrix; init=0.0), 1.0)
    residual = maximum(
        abs,
        adjoint(row_coefficients) * independent - matrix;
        init=0.0,
    )
    residual <= 10 * tolerance * scale || throw(ArgumentError(
        "SU(2) projected singlet row selection for exact rank $expected_rank " *
        "has residual $residual above tolerance $(10 * tolerance * scale)",
    ))
    return independent
end

function _pauli_su2_select_pivot_columns(
    singlet::AbstractMatrix{<:Complex};
    condition_limit::Real,
)
    limit = Float64(condition_limit)
    isfinite(limit) && limit > 0 || throw(DomainError(
        condition_limit,
        "condition_limit must be positive and finite",
    ))
    channel_count, source_count = size(singlet)
    channel_count <= source_count || throw(DimensionMismatch(
        "singlet transform has $channel_count rows but only $source_count source columns",
    ))
    if iszero(channel_count)
        return Int[], zeros(ComplexF64, source_count, 0), 1.0
    end

    matrix = ComplexF64.(singlet)
    factor = qr(matrix, ColumnNorm())
    pivots = sort!(Int[factor.p[idx] for idx in 1:channel_count])
    source_from_channels = adjoint(matrix)
    chart = source_from_channels[pivots, :]
    chart_condition = Float64(cond(chart))
    isfinite(chart_condition) && chart_condition <= limit || throw(ArgumentError(
        "SU(2) moment quotient coordinate-chart condition " *
        "$chart_condition exceeds limit $limit",
    ))
    reconstruction = adjoint(adjoint(chart) \ adjoint(source_from_channels))
    return pivots, Matrix{ComplexF64}(reconstruction), chart_condition
end

struct PauliSU2MomentOrbitQuotient{K,C}
    support_orbit::Tuple{Vararg{Int}}
    source_coordinate_mask::Symbol
    source_keys::Vector{K}
    pivot_keys::Vector{K}
    singlet_labels::Vector{Any}
    reconstruction::Matrix{C}
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    pivot_residual::Float64
    invariant_residual::Float64
    reconstruction_residual::Float64
    condition::Float64
end

struct PauliSU2MomentQuotientDescriptor{K,C}
    n_sites::Int
    sign_symmetry::Bool
    coordinate_provenance::Symbol
    orbits::Vector{PauliSU2MomentOrbitQuotient{K,C}}
    source_to_orbit_row::Dict{K,Tuple{Int,Int}}
    raw_moment_count::Int
    quotient_moment_count::Int
    support_orbit_count::Int
    singlet_channel_support_counts::Vector{Pair{Int,Int}}
    max_pivot_residual::Float64
    max_invariant_residual::Float64
    max_reconstruction_residual::Float64
    max_condition::Float64
end

function _pauli_su2_axes_code(axes::Tuple)
    code = 0
    stride = 1
    for axis in axes
        code += Int(axis) * stride
        stride *= 3
    end
    return code
end

function _pauli_su2_support_stabilizer_permutations(
    support_orbit::Tuple,
    n_sites::Int,
)
    positions = Dict{Int,Int}(Int(site) => idx for (idx, site) in pairs(support_orbit))
    permutations = Vector{Vector{Int}}()
    shifts = Int[]
    for shift in 0:(n_sites - 1)
        permutation = Vector{Int}(undef, length(support_orbit))
        preserves_support = true
        for (source_idx, site) in pairs(support_orbit)
            target_site = mod1(Int(site) + shift, n_sites)
            target_idx = get(positions, target_site, 0)
            if iszero(target_idx)
                preserves_support = false
                break
            end
            permutation[source_idx] = target_idx
        end
        preserves_support || continue
        push!(shifts, shift)
        push!(permutations, permutation)
    end
    isempty(permutations) && throw(ArgumentError(
        "support orbit $support_orbit has no identity stabilizer",
    ))
    return shifts, permutations
end

function _pauli_su2_permute_axes_code(code0::Int, permutation::Vector{Int})
    width = length(permutation)
    axes = digits(code0; base=3, pad=width)
    permuted = zeros(Int, width)
    for source_idx in 1:width
        permuted[permutation[source_idx]] = axes[source_idx]
    end
    return sum(permuted[idx] * 3^(idx - 1) for idx in 1:width; init=0)
end

function _pauli_su2_stabilizer_axis_patterns(
    support_orbit::Tuple,
    n_sites::Int,
)
    _, permutations = _pauli_su2_support_stabilizer_permutations(
        support_orbit,
        n_sites,
    )
    patterns = Set{Tuple{Vararg{UInt8}}}()
    for code0 in 0:(3^length(support_orbit) - 1)
        representative = minimum(
            _pauli_type_word_from_code(
                length(support_orbit),
                _pauli_su2_permute_axes_code(code0, permutation),
            )
            for permutation in permutations
        )
        push!(patterns, representative)
    end
    return patterns
end

function _pauli_su2_stabilizer_singlet_rows(
    support_orbit::Tuple,
    n_sites::Int,
    atol::Float64,
)
    singlet = _dense_sparse_transform_rows(
        _pauli_su2_word_singlet_rows(length(support_orbit); atol),
    )
    channel_count = size(singlet, 1)
    shifts, permutations = _pauli_su2_support_stabilizer_permutations(
        support_orbit,
        n_sites,
    )
    iszero(channel_count) && return singlet, shifts

    projector = zeros(ComplexF64, channel_count, channel_count)
    for permutation in permutations
        permuted_columns = Int[
            _pauli_su2_permute_axes_code(code0, permutation) + 1
            for code0 in 0:(size(singlet, 2) - 1)
        ]
        projector .+= singlet[:, permuted_columns] * adjoint(singlet)
    end
    projector ./= length(permutations)
    projector = Matrix{ComplexF64}(Hermitian((projector + adjoint(projector)) / 2))
    idempotence_residual = maximum(abs, projector * projector - projector; init=0.0)
    idempotence_residual <= 10 * atol || throw(ArgumentError(
        "support-stabilizer projector residual $idempotence_residual exceeds " *
        "tolerance on support orbit $support_orbit",
    ))
    projector_trace = LinearAlgebra.tr(projector)
    invariant_channel_count = round(Int, real(projector_trace))
    abs(real(projector_trace) - invariant_channel_count) <= 10 * atol ||
        throw(ArgumentError(
            "support-stabilizer projector has non-integral trace $projector_trace " *
            "on support orbit $support_orbit",
        ))
    iszero(invariant_channel_count) &&
        return zeros(ComplexF64, 0, size(singlet, 2)), shifts

    factor = qr(adjoint(projector), ColumnNorm())
    selected_rows = sort!(Int[factor.p[idx] for idx in 1:invariant_channel_count])
    invariant_rows = projector[selected_rows, :] * singlet
    return Matrix{ComplexF64}(invariant_rows), shifts
end

function _pauli_su2_orbit_quotient(
    source_keys::Vector{K},
    support_orbit::Tuple,
    source_coordinate_mask::Symbol,
    singlet::Matrix{ComplexF64},
    stabilizer_shifts::Vector{Int},
    atol::Float64,
    condition_limit::Float64,
) where {K}
    pivots, reconstruction, chart_condition = _pauli_su2_select_pivot_columns(
        singlet;
        condition_limit,
    )
    pivot_keys = K[_owned_moment_key(K, source_keys[idx]) for idx in pivots]
    permutation = sortperm(pivot_keys; lt=key_lt)
    pivots = pivots[permutation]
    pivot_keys = pivot_keys[permutation]
    reconstruction = reconstruction[:, permutation]

    source_from_channels = adjoint(singlet)
    chart = source_from_channels[pivots, :]
    pivot_residual = maximum(
        abs,
        reconstruction[pivots, :] - I;
        init=0.0,
    )
    invariant_coordinates = source_from_channels \ reconstruction
    invariant_residual = maximum(
        abs,
        reconstruction - source_from_channels * invariant_coordinates;
        init=0.0,
    )
    reconstruction_residual = maximum(
        abs,
        reconstruction * chart - source_from_channels;
        init=0.0,
    )
    maximum((pivot_residual, invariant_residual, reconstruction_residual)) <= atol ||
        throw(ArgumentError(
            "SU(2) moment quotient residual exceeds atol=$atol on support orbit " *
            "$support_orbit: pivot=$pivot_residual, invariant=$invariant_residual, " *
            "reconstruction=$reconstruction_residual",
        ))

    labels = Any[
        (
            feature=:pauli_su2_moment_singlet_channel,
            support_orbit=support_orbit,
            support_size=length(support_orbit),
            source_coordinate_mask,
            stabilizer_shifts=Tuple(stabilizer_shifts),
            invariant_channel=channel,
        )
        for channel in 1:size(singlet, 1)
    ]
    return PauliSU2MomentOrbitQuotient{K,ComplexF64}(
        support_orbit,
        source_coordinate_mask,
        K[_owned_moment_key(K, key) for key in source_keys],
        pivot_keys,
        labels,
        reconstruction,
        :complex_algebraic_float64,
        :complex_sqrt_rational,
        pivot_residual,
        invariant_residual,
        reconstruction_residual,
        chart_condition,
    )
end

function _pauli_su2_moment_quotient_descriptor(
    linear::MomentLinearData{K,C,M},
    n_sites::Integer;
    sign_symmetry::Bool=false,
    allow_registered_projection::Bool=false,
    atol::Real=1e-11,
    condition_limit::Real=1e10,
) where {K,C,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    n = Int(n_sites)
    n > 0 || throw(DomainError(n_sites, "n_sites must be positive"))
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "atol must be non-negative"))
    limit = Float64(condition_limit)
    isfinite(limit) && limit > 0 || throw(DomainError(
        condition_limit,
        "condition_limit must be positive and finite",
    ))

    basis = M[linear.key_to_monomial[key] for key in linear.moments]
    support_word_counts, support_axis_patterns =
        _pauli_translation_support_orbit_word_patterns(basis, n)
    supports = sort!(
        collect(keys(support_word_counts));
        by=support -> (length(support), support),
    )
    any(support -> length(support) > 8, supports) && throw(ArgumentError(
        "SU(2) moment quotient supports active support through 8",
    ))
    basis_lookup = _pauli_su2_translation_orbit_basis_lookup(basis, n)

    orbits = PauliSU2MomentOrbitQuotient{K,ComplexF64}[]
    source_to_orbit_row = Dict{K,Tuple{Int,Int}}()
    support_counts = Dict{Int,Int}()
    for support in supports
        complete_patterns = _pauli_su2_stabilizer_axis_patterns(support, n)
        sign_trivial_patterns = if sign_symmetry
            filter(
                axes -> _pauli_axis_signature(Tuple(axis + 0x01 for axis in axes)) == 0x00,
                complete_patterns,
            )
        else
            Set{Tuple{Vararg{UInt8}}}()
        end
        observed_patterns = support_axis_patterns[support]
        source_coordinate_mask = if observed_patterns == complete_patterns
            :complete
        elseif sign_symmetry && observed_patterns == sign_trivial_patterns
            :sign_trivial
        elseif sign_symmetry &&
                issubset(sign_trivial_patterns, observed_patterns) &&
                issubset(observed_patterns, complete_patterns)
            :sign_trivial_superset
        elseif allow_registered_projection && issubset(observed_patterns, complete_patterns)
            :generated_projection
        else
            throw(ArgumentError(
                "Pauli SU(2) moment quotient requires either the complete " *
                "support-stabilizer axis-pattern set or a registered superset of " *
                "the explicit trivial-sign mask. Projected coordinate masks are " *
                "accepted only with generated-translation provenance; support=$support, " *
                "observed=$(repr(observed_patterns)), " *
                "complete=$(repr(complete_patterns)), " *
                "sign_trivial=$(repr(sign_trivial_patterns))",
            ))
        end
        patterns = sort!(collect(observed_patterns))
        columns = Int[
            basis_lookup[(support=support, axes=axes)] for axes in patterns
        ]
        source_keys = K[
            _owned_moment_key(K, linear.moments[column]) for column in columns
        ]
        full_singlet, stabilizer_shifts = _pauli_su2_stabilizer_singlet_rows(
            support,
            n,
            tolerance,
        )
        if source_coordinate_mask in (:sign_trivial, :sign_trivial_superset)
            nontrivial_patterns = setdiff(complete_patterns, sign_trivial_patterns)
            nontrivial_columns = Int[
                _pauli_su2_axes_code(axes) + 1 for axes in nontrivial_patterns
            ]
            nontrivial_residual = maximum(
                abs,
                full_singlet[:, nontrivial_columns];
                init=0.0,
            )
            nontrivial_residual <= tolerance || throw(ArgumentError(
                "Pauli SU(2) singlet rows violate the declared sign-symmetry " *
                "coordinate mask on support $support: maximum nontrivial-sign " *
                "coefficient=$nontrivial_residual, atol=$tolerance",
            ))
        end
        representative_columns = Int[_pauli_su2_axes_code(axes) + 1 for axes in patterns]
        projected_singlet = full_singlet[:, representative_columns]
        observed_singlet = if source_coordinate_mask == :generated_projection
            exact_rank = _pauli_su2_exact_projected_singlet_rank(
                support,
                n,
                observed_patterns,
            )
            _pauli_su2_select_expected_singlet_rows(
                projected_singlet,
                exact_rank;
                atol=tolerance,
            )
        else
            projected_singlet
        end
        orbit = _pauli_su2_orbit_quotient(
            source_keys,
            support,
            source_coordinate_mask,
            observed_singlet,
            stabilizer_shifts,
            tolerance,
            limit,
        )
        push!(orbits, orbit)
        orbit_idx = length(orbits)
        for (source_row, key) in pairs(orbit.source_keys)
            haskey(source_to_orbit_row, key) && throw(ArgumentError(
                "duplicate SU(2) quotient source key $(repr(key))",
            ))
            source_to_orbit_row[_owned_moment_key(K, key)] = (orbit_idx, source_row)
        end
        support_size = length(support)
        support_counts[support_size] = get(support_counts, support_size, 0) +
            length(orbit.pivot_keys)
    end

    length(source_to_orbit_row) == length(linear.moments) || throw(ArgumentError(
        "SU(2) quotient mapped $(length(source_to_orbit_row)) source moments, " *
        "expected $(length(linear.moments))",
    ))
    quotient_moment_count = sum(length(orbit.pivot_keys) for orbit in orbits; init=0)
    return PauliSU2MomentQuotientDescriptor{K,ComplexF64}(
        n,
        sign_symmetry,
        allow_registered_projection ? :generated_translation_relaxation : :strict_input,
        orbits,
        source_to_orbit_row,
        length(linear.moments),
        quotient_moment_count,
        length(orbits),
        Pair{Int,Int}[
            support_size => support_counts[support_size]
            for support_size in sort!(collect(keys(support_counts)))
        ],
        maximum((orbit.pivot_residual for orbit in orbits); init=0.0),
        maximum((orbit.invariant_residual for orbit in orbits); init=0.0),
        maximum((orbit.reconstruction_residual for orbit in orbits); init=0.0),
        maximum((orbit.condition for orbit in orbits); init=1.0),
    )
end

struct PauliSU2QuotientTransform
    family::Symbol
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    base::Any
    descriptor::Any

    function PauliSU2QuotientTransform(base, descriptor)
        return new(
            :pauli_su2_moment_quotient,
            :complex_algebraic_float64,
            :complex_sqrt_rational,
            base,
            descriptor,
        )
    end
end

struct PauliSU2QuotientBlockOrigin <: BlockOrigin
    label::Any
    logical_row_labels::Vector{Any}
    transform::PauliSU2QuotientTransform
    source_origin::BlockOrigin

    function PauliSU2QuotientBlockOrigin(
        source_origin::BlockOrigin,
        descriptor::PauliSU2MomentQuotientDescriptor,
    )
        label = hasproperty(source_origin, :label) ?
            getproperty(source_origin, :label) : nothing
        logical_row_labels = hasproperty(source_origin, :logical_row_labels) ?
            Any[label for label in getproperty(source_origin, :logical_row_labels)] :
            Any[]
        base_transform = hasproperty(source_origin, :transform) ?
            getproperty(source_origin, :transform) : nothing
        return new(
            label,
            logical_row_labels,
            PauliSU2QuotientTransform(base_transform, descriptor),
            source_origin,
        )
    end
end

function _pauli_su2_convert_quotient_coefficient(
    ::Type{C},
    coefficient,
    atol::Float64,
    key,
) where {C<:Number}
    if C <: Real
        imaginary_residual = abs(imag(coefficient))
        imaginary_residual <= atol || throw(ArgumentError(
            "SU(2) moment quotient produced imaginary coefficient residual " *
            "$imaginary_residual for real moment key $(repr(key)); atol=$atol",
        ))
        return convert(C, real(coefficient))
    end
    return convert(C, coefficient)
end

function _pauli_su2_rewrite_form(
    form::LinearMomentForm{K,SC},
    descriptor::PauliSU2MomentQuotientDescriptor{K,QC},
    ::Type{C};
    atol::Real,
) where {K,SC<:Number,QC<:Number,C<:Number}
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "atol must be non-negative"))
    W = promote_type(SC, QC)
    wide_terms = Pair{K,W}[]
    for (source_key, source_coefficient) in form
        location = get(descriptor.source_to_orbit_row, source_key, nothing)
        location === nothing && throw(ArgumentError(
            "SU(2) moment quotient has no source coordinate for key " *
            "$(repr(source_key))",
        ))
        orbit_idx, source_row = location
        orbit = descriptor.orbits[orbit_idx]
        for pivot_col in eachindex(orbit.pivot_keys)
            coefficient = source_coefficient *
                orbit.reconstruction[source_row, pivot_col]
            iszero(coefficient) && continue
            push!(
                wide_terms,
                _owned_moment_key(K, orbit.pivot_keys[pivot_col]) =>
                    convert(W, coefficient),
            )
        end
    end

    wide_form = _linear_moment_form_from_owned_pairs!(wide_terms)
    terms = Pair{K,C}[]
    sizehint!(terms, length(wide_form))
    for (key, coefficient) in wide_form
        abs(coefficient) <= tolerance && continue
        converted = _pauli_su2_convert_quotient_coefficient(
            C,
            coefficient,
            tolerance,
            key,
        )
        iszero(converted) && continue
        push!(terms, _owned_moment_key(K, key) => converted)
    end
    return _linear_moment_form_from_owned_pairs!(terms)
end

struct PauliSU2MomentRewritePlan{K,C}
    pivot_keys::Vector{K}
    source_rows::Dict{K,Vector{Pair{Int,C}}}
end

struct PauliSU2MomentRewriteWorkspace{C}
    accumulator::Vector{C}
    touched::Vector{Int}
    seen::BitVector
end

function _pauli_su2_moment_rewrite_plan(
    descriptor::PauliSU2MomentQuotientDescriptor{K,C},
) where {K,C<:Number}
    pivot_keys = K[
        _owned_moment_key(K, key)
        for orbit in descriptor.orbits for key in orbit.pivot_keys
    ]
    sort!(pivot_keys; lt=key_lt)
    for idx in 2:length(pivot_keys)
        key_isequal(pivot_keys[idx - 1], pivot_keys[idx]) && throw(ArgumentError(
            "duplicate global SU(2) quotient pivot key $(repr(pivot_keys[idx]))",
        ))
    end
    pivot_index = Dict{K,Int}(
        _owned_moment_key(K, key) => idx for (idx, key) in pairs(pivot_keys)
    )

    source_rows = Dict{K,Vector{Pair{Int,C}}}()
    for orbit in descriptor.orbits
        orbit_pivot_indices = Int[
            _get_key_value(pivot_index, key, "global SU(2) quotient pivot index")
            for key in orbit.pivot_keys
        ]
        for (source_row, source_key) in pairs(orbit.source_keys)
            row = Pair{Int,C}[]
            sizehint!(row, length(orbit_pivot_indices))
            for (pivot_col, global_idx) in pairs(orbit_pivot_indices)
                coefficient = orbit.reconstruction[source_row, pivot_col]
                iszero(coefficient) && continue
                push!(row, global_idx => coefficient)
            end
            sort!(row; by=first)
            source_rows[_owned_moment_key(K, source_key)] = row
        end
    end
    length(source_rows) == descriptor.raw_moment_count || throw(ArgumentError(
        "SU(2) rewrite plan mapped $(length(source_rows)) source moments, " *
        "expected $(descriptor.raw_moment_count)",
    ))
    return PauliSU2MomentRewritePlan{K,C}(pivot_keys, source_rows)
end

function _pauli_su2_moment_rewrite_workspace(
    plan::PauliSU2MomentRewritePlan{K,PC},
    ::Type{SC},
) where {K,PC<:Number,SC<:Number}
    C = promote_type(PC, SC)
    return PauliSU2MomentRewriteWorkspace{C}(
        zeros(C, length(plan.pivot_keys)),
        Int[],
        falses(length(plan.pivot_keys)),
    )
end

function _pauli_su2_rewrite_form_with_plan!(
    workspace::PauliSU2MomentRewriteWorkspace{W},
    form::LinearMomentForm{K,SC},
    plan::PauliSU2MomentRewritePlan{K,PC},
    ::Type{C};
    atol::Real,
) where {W<:Number,K,SC<:Number,PC<:Number,C<:Number}
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "atol must be non-negative"))
    isempty(workspace.touched) || throw(ArgumentError(
        "SU(2) rewrite workspace was not cleared after its previous form",
    ))

    for (source_key, source_coefficient) in form
        row = get(plan.source_rows, source_key, nothing)
        row === nothing && throw(ArgumentError(
            "SU(2) rewrite plan has no source coordinate for key " *
            "$(repr(source_key))",
        ))
        for (pivot_idx, reconstruction_coefficient) in row
            if !workspace.seen[pivot_idx]
                workspace.seen[pivot_idx] = true
                push!(workspace.touched, pivot_idx)
            end
            workspace.accumulator[pivot_idx] += convert(
                W,
                source_coefficient * reconstruction_coefficient,
            )
        end
    end

    sort!(workspace.touched)
    terms = Pair{K,C}[]
    sizehint!(terms, length(workspace.touched))
    for pivot_idx in workspace.touched
        coefficient = workspace.accumulator[pivot_idx]
        workspace.accumulator[pivot_idx] = zero(W)
        workspace.seen[pivot_idx] = false
        abs(coefficient) <= tolerance && continue
        key = plan.pivot_keys[pivot_idx]
        converted = _pauli_su2_convert_quotient_coefficient(
            C,
            coefficient,
            tolerance,
            key,
        )
        iszero(converted) && continue
        push!(terms, _owned_moment_key(K, key) => converted)
    end
    empty!(workspace.touched)
    return LinearMomentForm{K,C}(terms, Val(:trusted))
end

function _pauli_su2_quotient_block_meta(
    meta::BlockMeta{M},
    descriptor::PauliSU2MomentQuotientDescriptor,
) where {M}
    return BlockMeta(
        meta.cone,
        PauliSU2QuotientBlockOrigin(meta.origin, descriptor),
        meta.row_labels,
    )
end

function _pauli_su2_quotient_redundant_zero_origin(origin::ConstraintOrigin)
    feature = _translation_zero_origin_histogram_key(origin).feature
    return feature in (
        :pauli_su2_base_zero,
        :pauli_su2_translation_orbit_singlet_channel_equality,
    )
end

function _pauli_su2_quotient_linear_data(
    linear::MomentLinearData{K,C,M},
    n_sites::Integer;
    sign_symmetry::Bool=false,
    allow_registered_projection::Bool=false,
    atol::Real=1e-11,
    condition_limit::Real=1e10,
    stage_times_ns=nothing,
) where {K,C<:Number,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "atol must be non-negative"))
    descriptor = _pauli_su2_moment_quotient_descriptor(
        linear,
        n_sites;
        sign_symmetry,
        allow_registered_projection,
        atol=tolerance,
        condition_limit,
    )
    rewrite_plan = _pauli_su2_moment_rewrite_plan(descriptor)
    rewrite_workspace = _pauli_su2_moment_rewrite_workspace(rewrite_plan, C)
    rewrite(form) = _pauli_su2_rewrite_form_with_plan!(
        rewrite_workspace,
        form,
        rewrite_plan,
        C;
        atol=tolerance,
    )
    builder = MomentLinearBuilder(K, C, M)
    for orbit in descriptor.orbits, pivot_key in orbit.pivot_keys
        monomial = _get_key_value(
            linear.key_to_monomial,
            pivot_key,
            "SU(2) quotient pivot monomial",
        )
        register_moment!(builder, pivot_key, monomial)
    end

    objective = rewrite(linear.objective_lin)
    add_objective_terms!(builder, objective)

    for (block_idx, block) in pairs(linear.psd_blocks_lin)
        entries = Matrix{LinearMomentForm{K,C}}(undef, size(block.entries))
        for entry_idx in eachindex(block.entries)
            entries[entry_idx] = rewrite(block.entries[entry_idx])
        end
        add_psd_block!(
            builder,
            block.meta.cone,
            entries,
            _pauli_su2_quotient_block_meta(block.meta, descriptor);
            constraint_idx=linear.psd_block_constraint_idx[block_idx],
        )
    end

    eliminated_count = 0
    eliminated_histogram = Dict{Any,Int}()
    for constraint in linear.zero_constraints
        histogram_key = _translation_zero_origin_histogram_key(constraint.origin)
        form = rewrite(constraint.form)
        if _pauli_su2_quotient_redundant_zero_origin(constraint.origin)
            constraint.kind == :zero || throw(ArgumentError(
                "SU(2) moment quotient cannot eliminate non-zero-kind invariant " *
                "constraint $(repr(constraint.kind)) from $(typeof(constraint.origin))",
            ))
            if !isempty(form)
                residual = maximum(abs(coefficient) for (_, coefficient) in form)
                throw(ArgumentError(
                    "SU(2) moment quotient labeled invariant constraint did not " *
                    "rewrite to zero: origin=$(repr(constraint.origin)), " *
                    "maximum residual=$residual, atol=$tolerance",
                ))
            end
            eliminated_count += 1
            eliminated_histogram[histogram_key] =
                get(eliminated_histogram, histogram_key, 0) + 1
            continue
        end
        if isempty(form)
            constraint.kind == :zero || throw(ArgumentError(
                "SU(2) moment quotient eliminated non-zero-kind scalar constraint " *
                "$(repr(constraint.kind)) from $(typeof(constraint.origin))",
            ))
            eliminated_count += 1
            eliminated_histogram[histogram_key] =
                get(eliminated_histogram, histogram_key, 0) + 1
            continue
        end
        _add_zero_constraint_trusted!(
            builder,
            form,
            constraint.origin;
            kind=constraint.kind,
            trusted_self_adjoint=constraint.trusted_self_adjoint,
        )
    end

    quotient = finalize!(
        builder;
        stage_times_ns,
        stage_prefix=:su2_moment_quotient,
    )
    histogram = Pair{Any,Int}[
        Pair{Any,Int}(entry.first, entry.second)
        for entry in sort!(collect(eliminated_histogram); by=entry -> repr(entry.first))
    ]
    return (
        linear=quotient,
        descriptor=descriptor,
        eliminated_zero_row_count=eliminated_count,
        eliminated_zero_feature_histogram=histogram,
    )
end
