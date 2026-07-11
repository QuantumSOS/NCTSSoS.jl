"""
    pauli_contiguous_chain_basis(registry, degree; periodic=true)

Build the sparse Pauli half-basis with contiguous support on a one-dimensional
chain.  The basis contains the identity and every Pauli word supported on a
contiguous window of length `1:degree`.  Periodic wrapping is enabled by default.

This is the basis used for large translation-invariant spin-chain relaxations:
for `degree < nqubits` and `periodic=true` its size is
`1 + nqubits * sum(3^k for k in 1:degree)` instead of the full NPA basis size.
"""
function pauli_contiguous_chain_basis(
    registry::VariableRegistry{PauliAlgebra,T},
    degree::Integer;
    periodic::Bool=true,
) where {T<:Integer}
    d = Int(degree)
    d >= 0 || throw(ArgumentError("`degree` must be nonnegative; got $d."))

    nqubits = _pauli_chain_nqubits(registry)
    basis = Set{NormalMonomial{PauliAlgebra,T}}()
    push!(basis, one(NormalMonomial{PauliAlgebra,T}))
    d == 0 && return sort!(collect(basis))

    max_width = min(d, nqubits)
    sizehint!(basis, _pauli_contiguous_chain_basis_size_hint(nqubits, max_width; periodic))

    for width in 1:max_width
        starts = periodic ? (1:nqubits) : (1:(nqubits - width + 1))
        for start in starts
            sites = if periodic
                [mod1(start + offset, nqubits) for offset in 0:(width - 1)]
            else
                collect(start:(start + width - 1))
            end
            for code in 0:(3^width - 1)
                push!(basis, _pauli_chain_word(T, sites, code))
            end
        end
    end

    return sort!(collect(basis))
end

function _pauli_chain_nqubits(registry::VariableRegistry{PauliAlgebra,T}) where {T<:Integer}
    !isempty(registry.idx_to_variables) || throw(ArgumentError("Pauli chain basis needs a non-empty registry."))
    nqubits = maximum(_pauli_site(idx) for idx in keys(registry.idx_to_variables))
    for site in 1:nqubits, pauli_type in (_PAULI_X_TYPE, _PAULI_Y_TYPE, _PAULI_Z_TYPE)
        idx = convert(T, _pauli_index(site, pauli_type))
        haskey(registry.idx_to_variables, idx) || throw(ArgumentError(
            "Pauli chain basis needs a complete site-contiguous Pauli registry; missing index $idx for site $site."
        ))
    end
    length(registry.idx_to_variables) == 3 * nqubits || throw(ArgumentError(
        "Pauli chain basis needs a complete site-contiguous Pauli registry with 3 variables per site."
    ))
    return nqubits
end

function _pauli_chain_ops_from_registry(registry::VariableRegistry{PauliAlgebra,T}) where {T<:Unsigned}
    nqubits = _pauli_chain_nqubits(registry)
    M = NormalMonomial{PauliAlgebra,T}
    make_ops(pauli_type) = M[M(T[convert(T, _pauli_index(site, pauli_type))]) for site in 1:nqubits]
    return (
        make_ops(_PAULI_X_TYPE),
        make_ops(_PAULI_Y_TYPE),
        make_ops(_PAULI_Z_TYPE),
    )
end

function _pauli_contiguous_chain_basis_size_hint(nqubits::Integer, degree::Integer; periodic::Bool)
    if periodic
        return 1 + Int(nqubits) * sum(3^width for width in 1:Int(degree); init=0)
    end
    return 1 + sum((Int(nqubits) - width + 1) * 3^width for width in 1:Int(degree); init=0)
end

function _pauli_chain_word(::Type{T}, sites::AbstractVector{<:Integer}, code::Integer) where {T<:Integer}
    raw = Vector{T}(undef, length(sites))
    value = Int(code)
    @inbounds for i in eachindex(sites)
        pauli_type = value % 3
        value ÷= 3
        raw[i] = convert(T, _pauli_index(sites[i], pauli_type))
    end

    word, phase = simplify(PauliAlgebra, raw)
    phase == UInt8(0) || throw(ArgumentError(
        "Internal Pauli chain basis construction produced a phased word; this should not happen for distinct sites."
    ))
    return NormalMonomial{PauliAlgebra,T}(word)
end

"""
    pauli_sign_symmetry(nqubits; integer_type=Int)

Construct the global Pauli sign symmetry induced by conjugation with
`∏ᵢ σᵢᶻ`: `Xᵢ ↦ -Xᵢ`, `Yᵢ ↦ -Yᵢ`, and `Zᵢ ↦ Zᵢ`.
"""
function pauli_sign_symmetry(nqubits::Integer; integer_type::Type{T}=Int) where {T<:Integer}
    n = _clifford_validate_nqubits(nqubits)
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    sizehint!(images, 2 * n)
    for site in 1:n
        x = _pauli_letter(T, site, _PAULI_X_TYPE)
        y = _pauli_letter(T, site, _PAULI_Y_TYPE)
        images[x] = (-1, x)
        images[y] = (-1, y)
    end
    return CliffordSymmetry(images; nqubits=n)
end


# =============================================================================
# Pauli spin-chain helpers and translation-invariant relaxations
# =============================================================================

"""
    pauli_contiguous_chain_basis((σx, σy, σz), order; periodic=true)

Build the contiguous-support Pauli half-basis used by large 1D spin-chain
relaxations:

```math
B_d = {1} ∪ {σ_i^{a_1} σ_{i+1}^{a_2} ⋯ σ_{i+ℓ-1}^{a_ℓ} : 1 ≤ ℓ ≤ d}.
```

For a periodic chain with `order < N`, this basis has
`1 + N * sum(3^ℓ for ℓ in 1:order)` elements.  This is intentionally not the
full NPA basis; it is a local-basis relaxation, so bounds remain rigorous but
may be weaker.
"""
function pauli_contiguous_chain_basis(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    periodic::Bool=true,
)
    σx, σy, σz, n = _validate_pauli_chain_ops(ops)
    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    if periodic && d > n
        throw(ArgumentError("Periodic contiguous Pauli chain basis currently requires `order <= n_sites`; got order=$d, n_sites=$n."))
    end

    M = eltype(σx)
    T = eltype(first(σx).word)
    basis = M[one(first(σx))]
    d == 0 && return basis

    # Exact size hint for the common non-wrapping case.
    if periodic
        sizehint!(basis, 1 + n * sum(3^ℓ for ℓ in 1:d))
    else
        sizehint!(basis, 1 + sum(max(n - ℓ + 1, 0) * 3^ℓ for ℓ in 1:d))
    end

    axes = (σx, σy, σz)
    word = T[]
    for ℓ in 1:d
        nstarts = periodic ? n : max(n - ℓ + 1, 0)
        for start in 1:nstarts
            for code0 in 0:(3^ℓ - 1)
                empty!(word)
                code = code0
                for offset in 0:(ℓ - 1)
                    axis = mod(code, 3) + 1
                    code = div(code, 3)
                    site = periodic ? mod1(start + offset, n) : start + offset
                    push!(word, only(axes[axis][site].word))
                end
                simplified, phase = simplify(PauliAlgebra, word)
                phase == 0x00 || throw(ArgumentError(
                    "Internal error: contiguous Pauli basis construction produced non-real phase $phase."
                ))
                push!(basis, M(copy(simplified)))
            end
        end
    end

    return sorted_unique!(basis)
end

function _pauli_contiguous_chain_orbit_representatives(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    periodic::Bool=true,
)
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))

    if !periodic || d >= n
        basis = pauli_contiguous_chain_basis(ops, d; periodic)
        return _translation_orbit_representatives(basis, n)
    end

    M = eltype(σx)
    T = eltype(first(σx).word)
    reps = M[one(first(σx))]
    d == 0 && return reps
    sizehint!(reps, 1 + sum(3^ℓ for ℓ in 1:d; init=0))

    axes = (ops[1], ops[2], ops[3])
    word = T[]
    for ℓ in 1:d
        for code0 in 0:(3^ℓ - 1)
            empty!(word)
            code = code0
            for offset in 0:(ℓ - 1)
                axis = mod(code, 3) + 1
                code = div(code, 3)
                push!(word, only(axes[axis][offset + 1].word))
            end
            simplified, phase = simplify(PauliAlgebra, word)
            phase == 0x00 || throw(ArgumentError(
                "Internal error: contiguous Pauli orbit representative construction produced non-real phase $phase."
            ))
            mono = M(copy(simplified))
            push!(reps, _translation_orbit_representative(mono, n))
        end
    end

    return sorted_unique!(reps)
end

function _histogram_pairs(values)
    counts = Dict{Int,Int}()
    for value in values
        key = Int(value)
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by=first)
end

function _value_histogram_pairs(values)
    counts = Dict{Any,Int}()
    for value in values
        counts[value] = get(counts, value, 0) + 1
    end
    return sort!(collect(counts); by=pair -> string(first(pair)))
end

function _label_size_histogram(labels::AbstractVector, sizes::AbstractVector)
    length(labels) == length(sizes) || throw(ArgumentError(
        "Cannot build label-size histogram from $(length(labels)) labels and $(length(sizes)) sizes."
    ))

    counts = Dict{Any,Int}()
    for (label, size) in zip(labels, sizes)
        key = label isa NamedTuple ?
            merge(label, (size=Int(size),)) :
            (label=label, size=Int(size))
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by=pair -> repr(pair.first))
end

const _PAULI_AXIS_PERMUTATIONS = (
    (0x01, 0x02, 0x03),
    (0x01, 0x03, 0x02),
    (0x02, 0x01, 0x03),
    (0x02, 0x03, 0x01),
    (0x03, 0x01, 0x02),
    (0x03, 0x02, 0x01),
)

function _pauli_axis_word_from_code(width::Int, code0::Int)
    word = Vector{UInt8}(undef, width)
    code = code0
    @inbounds for idx in 1:width
        word[idx] = UInt8(mod(code, 3) + 1)
        code = div(code, 3)
    end
    return Tuple(word)
end

function _pauli_axis_signature(word::Tuple{Vararg{UInt8}})
    px = false
    py = false
    pz = false
    for axis in word
        if axis == 0x01
            px = !px
        elseif axis == 0x02
            py = !py
        else
            pz = !pz
        end
    end
    return UInt8((xor(px, py) ? 0x01 : 0x00) | (xor(py, pz) ? 0x02 : 0x00))
end

function _pauli_axis_signature_counts(order::Int)
    counts = Dict{UInt8,Int}()
    for width in 1:order
        for code in 0:(3^width - 1)
            signature = _pauli_axis_signature(_pauli_axis_word_from_code(width, code))
            counts[signature] = get(counts, signature, 0) + 1
        end
    end
    return counts
end

function _pauli_axis_palindrome_counts(order::Int)
    even_counts = Dict{UInt8,Int}()
    odd_counts = Dict{UInt8,Int}()
    for width in 1:order
        for code in 0:(3^width - 1)
            word = _pauli_axis_word_from_code(width, code)
            reverse(collect(word)) == collect(word) || continue
            signature = _pauli_axis_signature(word)
            counts = iseven(width) ? even_counts : odd_counts
            counts[signature] = get(counts, signature, 0) + 1
        end
    end
    return even_counts, odd_counts
end

function _pauli_axis_word_orbit(word::Tuple{Vararg{UInt8}})
    orbit = Set{Tuple{Vararg{UInt8}}}()
    for perm in _PAULI_AXIS_PERMUTATIONS
        push!(orbit, ntuple(idx -> perm[Int(word[idx])], length(word)))
    end
    return orbit
end

function _pauli_axis_word_orbit_summary(order::Int)
    representatives = Set{Tuple{Vararg{UInt8}}}()
    orbit_sizes = Int[]

    push!(representatives, ())
    push!(orbit_sizes, 1)

    for width in 1:order
        for code in 0:(3^width - 1)
            word = _pauli_axis_word_from_code(width, code)
            orbit = _pauli_axis_word_orbit(word)
            representative = minimum(collect(orbit))
            representative in representatives && continue
            push!(representatives, representative)
            push!(orbit_sizes, length(orbit))
        end
    end

    return (
        axis_orbit_count=length(representatives),
        axis_orbit_size_histogram=_histogram_pairs(orbit_sizes),
        max_axis_orbit_size=isempty(orbit_sizes) ? 0 : maximum(orbit_sizes),
    )
end

function _pauli_axis_orbit_summary(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    orbit_reps::Vector{M},
    order::Int,
    ;
    strict::Bool=true,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    eltype(σx) == M || throw(ArgumentError(
        "Pauli axis orbit diagnostics require operators and orbit representatives to use the same monomial type."
    ))

    axis_group = CliffordSymmetryGroup(
        pauli_global_axis_rotation_generators(ops);
        nqubits=n,
        integer_type=T,
    )
    rep_set = Set(orbit_reps)
    remaining = Set(orbit_reps)
    axis_orbits = Vector{Vector{M}}()

    open_summary() = (
        n_sites=n,
        order=order,
        translation_orbit_count=length(orbit_reps),
        axis_orbit_closed=false,
        axis_orbit_count=0,
        axis_group_order=length(axis_group),
        axis_orbit_sizes=Int[],
        axis_orbit_size_histogram=Pair{Int,Int}[],
        max_axis_orbit_size=0,
        reduction_ratio=0.0,
        translation_orbit_representatives=orbit_reps,
        axis_orbit_representatives=M[],
        axis_orbits=Vector{M}[],
    )

    for rep in orbit_reps
        rep in remaining || continue
        orbit = Set{M}()
        for g in axis_group
            _, image = _act_monomial(_clifford_group_value(g), rep)
            reduced = _translation_orbit_representative(image, n)
            if !(reduced in rep_set)
                strict || return open_summary()
                throw(ArgumentError(
                    "Global Pauli-axis rotation maps translation representative $rep outside the supplied representative set: $reduced."
                ))
            end
            push!(orbit, reduced)
        end
        orbit_vec = sort!(collect(orbit))
        push!(axis_orbits, orbit_vec)
        for mono in orbit_vec
            delete!(remaining, mono)
        end
    end

    orbit_sizes = length.(axis_orbits)
    axis_reps = M[orbit[1] for orbit in axis_orbits]
    return (
        n_sites=n,
        order=order,
        translation_orbit_count=length(orbit_reps),
        axis_orbit_closed=true,
        axis_orbit_count=length(axis_orbits),
        axis_group_order=length(axis_group),
        axis_orbit_sizes=orbit_sizes,
        axis_orbit_size_histogram=_histogram_pairs(orbit_sizes),
        max_axis_orbit_size=isempty(orbit_sizes) ? 0 : maximum(orbit_sizes),
        reduction_ratio=length(axis_orbits) == 0 ? Inf : length(orbit_reps) / length(axis_orbits),
        translation_orbit_representatives=orbit_reps,
        axis_orbit_representatives=axis_reps,
        axis_orbits=axis_orbits,
    )
end

function _pauli_axis_translation_orbit_representatives(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis=nothing,
    periodic::Bool=true,
    context::AbstractString="Pauli axis orbit diagnostics",
    basis_context::AbstractString="axis-orbit diagnostic basis",
    require_identity::Bool=true,
    basis_is_orbit_representatives::Bool=false,
)
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    M = eltype(σx)
    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))

    orbit_reps = if isnothing(basis)
        _pauli_contiguous_chain_orbit_representatives(ops, d; periodic)
    else
        local_basis = M[mono for mono in basis]
        !require_identity || one(first(σx)) in local_basis || throw(ArgumentError(
            "$context require a basis containing the identity."
        ))
        _check_pauli_chain_support(local_basis, n; context=basis_context)
        if basis_is_orbit_representatives
            length(Set(local_basis)) == length(local_basis) || throw(ArgumentError(
                "$context require unique translation-orbit representatives."
            ))
            local_basis
        else
            _check_translation_basis_closure(local_basis, n)
            _translation_orbit_representatives(local_basis, n)
        end
    end

    return n, d, orbit_reps
end

"""
    pauli_axis_orbit_diagnostics((σx, σy, σz), order; basis=nothing, periodic=true)

Group translation orbit representatives of the contiguous Pauli-chain basis
under global Pauli-axis Clifford rotations.  This is a structural diagnostic for
the missing axis-rotation part of the specialized translation backend.
"""
function pauli_axis_orbit_diagnostics(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis=nothing,
    periodic::Bool=true,
    basis_is_orbit_representatives::Bool=false,
)
    _, d, orbit_reps = _pauli_axis_translation_orbit_representatives(
        ops,
        order;
        basis,
        periodic,
        context="Pauli axis orbit diagnostics",
        basis_context="axis-orbit diagnostic basis",
        require_identity=true,
    )
    return _pauli_axis_orbit_summary(ops, orbit_reps, d; strict=true)
end

function _pauli_axis_rotation_generator_action_matrices(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    orbit_reps::Vector{M},
    n_sites::Integer,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    rep_indices = Dict{M,Int}(rep => idx for (idx, rep) in pairs(orbit_reps))
    matrices = SparseMatrixCSC{Int,Int}[]
    sizehint!(matrices, 2)
    dim = length(orbit_reps)

    for (generator_idx, rotation) in pairs(pauli_global_axis_rotation_generators(ops))
        action = spzeros(Int, dim, dim)
        for (col, source) in pairs(orbit_reps)
            coefficient, image = _act_monomial(rotation, source)
            target = _translation_orbit_representative(image, n_sites)
            row = get(rep_indices, target, 0)
            row == 0 && throw(ArgumentError(
                "Pauli axis rotation action matrices require the translation orbit basis " *
                "to be closed under global Pauli-axis rotations; generator $generator_idx " *
                "maps $source to $target, which is missing."
            ))
            coefficient ∈ (-1, 1) || throw(ArgumentError(
                "Pauli axis rotation generator $generator_idx maps $source with non-sign coefficient $coefficient."
            ))
            action[row, col] = Int(coefficient)
        end
        push!(matrices, action)
    end

    return matrices
end

function _sparse_signed_matrix_key(matrix::SparseMatrixCSC{Int,Int})
    rows, cols, values = findnz(matrix)
    return (Tuple(rows), Tuple(cols), Tuple(values))
end

function _pauli_axis_rotation_group_action_matrices(
    generator_matrices::Vector{SparseMatrixCSC{Int,Int}},
)
    isempty(generator_matrices) && return SparseMatrixCSC{Int,Int}[]
    dim = size(first(generator_matrices), 1)
    all(size(matrix) == (dim, dim) for matrix in generator_matrices) || throw(ArgumentError(
        "Pauli axis rotation generator action matrices must all be square matrices with the same size."
    ))

    identity = sparse(1:dim, 1:dim, ones(Int, dim), dim, dim)
    group_matrices = SparseMatrixCSC{Int,Int}[identity]
    seen = Set{Any}([_sparse_signed_matrix_key(identity)])
    idx = 1
    while idx <= length(group_matrices)
        current = group_matrices[idx]
        for generator in generator_matrices
            candidate = generator * current
            key = _sparse_signed_matrix_key(candidate)
            key in seen && continue
            push!(seen, key)
            push!(group_matrices, candidate)
        end
        idx += 1
    end
    return group_matrices
end

function _pauli_axis_rotation_group_element_orders(
    group_matrices::Vector{SparseMatrixCSC{Int,Int}},
)
    isempty(group_matrices) && return Int[]
    identity_key = _sparse_signed_matrix_key(first(group_matrices))
    orders = Int[]
    sizehint!(orders, length(group_matrices))
    for element in group_matrices
        product = element
        order = 1
        while _sparse_signed_matrix_key(product) != identity_key
            product = product * element
            order += 1
            order <= length(group_matrices) || throw(ArgumentError(
                "Pauli axis rotation group action matrix did not close within the enumerated group size."
            ))
        end
        push!(orders, order)
    end
    return orders
end

const _PAULI_AXIS_CLASS_ORDER = (:identity, :edge_180, :face_180, :body_120, :face_90)

const _PAULI_AXIS_IRREP_SPECS = (
    (
        label=:trivial,
        dimension=1,
        characters=(identity=1, edge_180=1, face_180=1, body_120=1, face_90=1),
    ),
    (
        label=:orientation_sign,
        dimension=1,
        characters=(identity=1, edge_180=-1, face_180=1, body_120=1, face_90=-1),
    ),
    (
        label=:two_dimensional,
        dimension=2,
        characters=(identity=2, edge_180=0, face_180=2, body_120=-1, face_90=0),
    ),
    (
        label=:standard,
        dimension=3,
        characters=(identity=3, edge_180=1, face_180=-1, body_120=0, face_90=-1),
    ),
    (
        label=:axis_vector,
        dimension=3,
        characters=(identity=3, edge_180=-1, face_180=-1, body_120=0, face_90=1),
    ),
)

function _pauli_axis_rotation_class_label(
    vector_action::AbstractMatrix{<:Integer},
    element_order::Integer,
)
    order = Int(element_order)
    order == 1 && return :identity
    order == 3 && return :body_120
    order == 4 && return :face_90
    if order == 2
        diagonal_nonzeros = count(!iszero, diag(vector_action))
        diagonal_nonzeros == 3 && return :face_180
        diagonal_nonzeros == 1 && return :edge_180
    end
    throw(ArgumentError(
        "Could not classify Pauli axis-rotation group element with order $order and action $vector_action."
    ))
end

function _pauli_axis_class_size_histogram(class_labels::AbstractVector{Symbol})
    counts = Dict{Symbol,Int}()
    for label in class_labels
        counts[label] = get(counts, label, 0) + 1
    end
    return Pair{Symbol,Int}[label => counts[label] for label in _PAULI_AXIS_CLASS_ORDER if haskey(counts, label)]
end

function _pauli_axis_rotation_class_labels(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
)
    vector_actions = pauli_axis_rotation_action_matrices(ops, 1)
    vector_actions.axis_group_order == 24 || throw(ArgumentError(
        "Pauli axis irrep projectors require the faithful order-1 axis representation; " *
        "got group order $(vector_actions.axis_group_order)."
    ))
    group_actions = Matrix.(vector_actions.group_action_matrices)
    orders = vector_actions.group_element_orders
    indices = Dict(
        mono => idx for (idx, mono) in pairs(vector_actions.translation_orbit_representatives)
    )
    σx, σy, σz, _ = _validate_pauli_chain_ops(ops)
    vector_indices = [indices[σx[1]], indices[σy[1]], indices[σz[1]]]
    return Symbol[
        _pauli_axis_rotation_class_label(action[vector_indices, vector_indices], order)
        for (action, order) in zip(group_actions, orders)
    ]
end

function _matrix_rank_from_projector(projector::AbstractMatrix{<:Real}, atol::Real)
    return count(value -> value > atol, svdvals(Matrix(projector)))
end

"""
    pauli_axis_rotation_action_matrices((σx, σy, σz), order; basis=nothing, periodic=true)

Return sparse signed-permutation matrices for the two global Pauli-axis
Clifford generators on the translation-orbit representative basis.  Matrix
columns are source representatives and rows are translation-reduced images.
The returned group matrices enumerate the finite closure generated by those
two matrices, with identity first.
"""
function pauli_axis_rotation_action_matrices(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis=nothing,
    periodic::Bool=true,
    basis_is_orbit_representatives::Bool=false,
)
    n, d, orbit_reps = _pauli_axis_translation_orbit_representatives(
        ops,
        order;
        basis,
        periodic,
        context="Pauli axis rotation action matrices",
        basis_context="axis-rotation action matrix basis",
        require_identity=false,
        basis_is_orbit_representatives,
    )
    generator_action_matrices = _pauli_axis_rotation_generator_action_matrices(ops, orbit_reps, n)
    group_action_matrices = _pauli_axis_rotation_group_action_matrices(generator_action_matrices)
    return (
        n_sites=n,
        order=d,
        translation_orbit_count=length(orbit_reps),
        generator_labels=[:global_h, :global_s],
        axis_group_order=length(group_action_matrices),
        translation_orbit_representatives=orbit_reps,
        generator_action_matrices=generator_action_matrices,
        group_action_matrices=group_action_matrices,
        group_element_orders=_pauli_axis_rotation_group_element_orders(group_action_matrices),
    )
end

"""
    pauli_axis_rotation_irrep_projectors((σx, σy, σz), order; basis=nothing, periodic=true, atol=1e-10)

Return central isotypic projectors for the finite global Pauli-axis rotation
group acting on the translation-orbit representative basis.  These dense
projectors are structural data for the standalone axis-rotation PSD block reducer.
"""
function pauli_axis_rotation_irrep_projectors(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis=nothing,
    periodic::Bool=true,
    atol::Real=1e-10,
    basis_is_orbit_representatives::Bool=false,
)
    actions = pauli_axis_rotation_action_matrices(
        ops,
        order;
        basis,
        periodic,
        basis_is_orbit_representatives,
    )
    if actions.axis_group_order == 1
        identity = Matrix{Float64}(I, actions.translation_orbit_count, actions.translation_orbit_count)
        return (
            n_sites=actions.n_sites,
            order=actions.order,
            translation_orbit_count=actions.translation_orbit_count,
            axis_group_order=actions.axis_group_order,
            class_size_histogram=[:identity => 1],
            irrep_labels=[:trivial],
            irrep_dimensions=[1],
            projector_ranks=[actions.translation_orbit_count],
            irrep_multiplicities=[actions.translation_orbit_count],
            projector_matrices=[identity],
        )
    end
    actions.axis_group_order == 24 || throw(ArgumentError(
        "Pauli axis irrep projectors require a faithful 24-element axis action; " *
        "got group order $(actions.axis_group_order)."
    ))

    class_labels = _pauli_axis_rotation_class_labels(ops)
    length(class_labels) == length(actions.group_action_matrices) || throw(ArgumentError(
        "Pauli axis irrep projector classification produced $(length(class_labels)) labels " *
        "for $(length(actions.group_action_matrices)) group matrices."
    ))

    group_matrices = Matrix{Float64}.(actions.group_action_matrices)
    group_order = length(group_matrices)
    projector_matrices = Matrix{Float64}[]
    projector_ranks = Int[]
    irrep_multiplicities = Int[]

    for spec in _PAULI_AXIS_IRREP_SPECS
        projector = zeros(Float64, actions.translation_orbit_count, actions.translation_orbit_count)
        for (class_label, matrix) in zip(class_labels, group_matrices)
            projector .+= getproperty(spec.characters, class_label) .* matrix
        end
        projector .*= spec.dimension / group_order
        rank = _matrix_rank_from_projector(projector, atol)
        rank % spec.dimension == 0 || throw(ArgumentError(
            "Pauli axis irrep $(spec.label) projector rank $rank is not divisible " *
            "by irrep dimension $(spec.dimension)."
        ))
        push!(projector_matrices, projector)
        push!(projector_ranks, rank)
        push!(irrep_multiplicities, div(rank, spec.dimension))
    end

    return (
        n_sites=actions.n_sites,
        order=actions.order,
        translation_orbit_count=actions.translation_orbit_count,
        axis_group_order=actions.axis_group_order,
        class_size_histogram=_pauli_axis_class_size_histogram(class_labels),
        irrep_labels=[spec.label for spec in _PAULI_AXIS_IRREP_SPECS],
        irrep_dimensions=[spec.dimension for spec in _PAULI_AXIS_IRREP_SPECS],
        projector_ranks=projector_ranks,
        irrep_multiplicities=irrep_multiplicities,
        projector_matrices=projector_matrices,
    )
end

function _pauli_projector_range_basis(
    projector::AbstractMatrix{<:Number},
    atol::Real,
)
    sym_projector = Hermitian(0.5 .* (Matrix(projector) .+ Matrix(projector)'))
    eig = eigen(sym_projector)
    keep = findall(value -> value > atol, eig.values)
    isempty(keep) && return zeros(Float64, size(projector, 1), 0)
    return Matrix(eig.vectors[:, keep])
end

"""
    pauli_axis_rotation_isotypic_transform((σx, σy, σz), order; basis=nothing, periodic=true, atol=1e-10)

Return an orthonormal transform whose column groups span the isotypic
components of the finite global Pauli-axis rotation action on the
translation-orbit representative basis.
"""
function pauli_axis_rotation_isotypic_transform(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis=nothing,
    periodic::Bool=true,
    atol::Real=1e-10,
    basis_is_orbit_representatives::Bool=false,
)
    projectors = pauli_axis_rotation_irrep_projectors(
        ops,
        order;
        basis,
        periodic,
        atol,
        basis_is_orbit_representatives,
    )
    transform_blocks = Matrix{Float64}[]
    sizehint!(transform_blocks, length(projectors.projector_matrices))
    for projector in projectors.projector_matrices
        push!(transform_blocks, _pauli_projector_range_basis(projector, atol))
    end

    transform_matrix = zeros(Float64, projectors.translation_orbit_count, 0)
    for block in transform_blocks
        transform_matrix = hcat(transform_matrix, block)
    end

    isotypic_block_sizes = size.(transform_blocks, 2)
    sum(isotypic_block_sizes) == projectors.translation_orbit_count || throw(ArgumentError(
        "Pauli axis isotypic transform recovered $(sum(isotypic_block_sizes)) columns " *
        "for a $(projectors.translation_orbit_count)-dimensional representation."
    ))

    return (
        n_sites=projectors.n_sites,
        order=projectors.order,
        translation_orbit_count=projectors.translation_orbit_count,
        axis_group_order=projectors.axis_group_order,
        class_size_histogram=projectors.class_size_histogram,
        irrep_labels=projectors.irrep_labels,
        irrep_dimensions=projectors.irrep_dimensions,
        isotypic_block_sizes=isotypic_block_sizes,
        irrep_multiplicities=projectors.irrep_multiplicities,
        transform_blocks=transform_blocks,
        transform_matrix=transform_matrix,
    )
end

"""
    closure_check(:contiguous_rdm, basis; n_sites, k)
    closure_check(:linear_state_opt, basis; n_sites, hamiltonian, test_width, sign_symmetry=false)
    closure_check(:psd_state_opt, basis; n_sites, hamiltonian, test_width, sign_symmetry=false)

Verify that the Pauli half-basis generates every translation-reduced moment
needed by an additional translation-invariant constraint family before model
construction.
"""
function closure_check(
    feature::Symbol,
    basis::Vector{M};
    n_sites=nothing,
    k=nothing,
    hamiltonian=nothing,
    test_width=nothing,
    sign_symmetry::Bool=false,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    n_sites === nothing && throw(ArgumentError("closure_check($feature, ...) requires keyword `n_sites`."))
    if feature == :contiguous_rdm
        k === nothing && throw(ArgumentError("closure_check(:contiguous_rdm, ...) requires keyword `k`."))
        return _check_contiguous_rdm_closure(basis, Int(n_sites), Int(k))
    elseif feature == :linear_state_opt
        hamiltonian === nothing && throw(ArgumentError("closure_check(:linear_state_opt, ...) requires keyword `hamiltonian`."))
        test_width === nothing && throw(ArgumentError("closure_check(:linear_state_opt, ...) requires keyword `test_width`."))
        return _check_linear_state_opt_closure(
            basis,
            Int(n_sites),
            hamiltonian,
            Int(test_width);
            sign_symmetry,
        )
    elseif feature == :psd_state_opt
        hamiltonian === nothing && throw(ArgumentError("closure_check(:psd_state_opt, ...) requires keyword `hamiltonian`."))
        test_width === nothing && throw(ArgumentError("closure_check(:psd_state_opt, ...) requires keyword `test_width`."))
        return _check_psd_state_opt_closure(
            basis,
            Int(n_sites),
            hamiltonian,
            Int(test_width);
            sign_symmetry,
        )
    end
    throw(ArgumentError(
        "Unsupported Pauli closure_check feature $feature. Supported features: :contiguous_rdm, :linear_state_opt, :psd_state_opt."
    ))
end

function _check_contiguous_rdm_closure(
    basis::Vector{M},
    n_sites::Int,
    k::Int,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    n_sites > 0 || throw(DomainError(n_sites, "`n_sites` must be positive."))
    0 <= k <= n_sites || throw(DomainError(k, "`k` must satisfy 0 <= k <= n_sites."))
    isempty(basis) && throw(ArgumentError("closure_check(:contiguous_rdm, ...) requires a non-empty basis."))
    _check_pauli_chain_support(basis, n_sites; context="closure-check basis")

    generated = _translation_reduced_pair_moment_set(basis, n_sites)
    required = _contiguous_rdm_required_moment_set(M, n_sites, k)
    missing = M[m for m in required if !(m in generated)]
    isempty(missing) && return nothing

    sort!(missing)
    shown = join((sprint(show, mono) for mono in Iterators.take(missing, 5)), ", ")
    length(missing) > 5 && (shown *= ", ...")
    throw(ArgumentError(
        "closure_check(:contiguous_rdm, ...) missing $(length(missing)) translation-reduced moment(s) for k=$k: [$shown]. Increase the half-basis order before adding this RDM block."
    ))
end

function _translation_reduced_pair_moment_set(
    basis::Vector{M},
    n_sites::Int,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    moments = Set{M}()
    one_mono = one(first(basis))
    buf = T[]
    for row in basis, col in basis
        word, phase_k = simplify!(
            PauliAlgebra,
            _neat_dot3!(buf, row, one_mono, col),
        )
        phase_k == 0x04 && continue
        mono = M(copy(word))
        push!(moments, _translation_orbit_representative(mono, n_sites))
    end
    return moments
end

function _contiguous_rdm_required_moment_set(
    ::Type{M},
    n_sites::Int,
    k::Int,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    required = Set{M}()
    word = T[]
    for code0 in 0:(4^k - 1)
        empty!(word)
        code = code0
        for offset in 0:(k - 1)
            pauli_code = mod(code, 4)
            code = div(code, 4)
            pauli_code == 0 && continue
            push!(word, _pauli_index_from_site_type(T, offset + 1, pauli_code - 1))
        end
        simplified, phase = simplify(PauliAlgebra, word)
        phase == 0x00 || throw(ArgumentError(
            "Internal error: contiguous RDM moment enumeration produced non-real phase $phase."
        ))
        mono = M(copy(simplified))
        push!(required, _translation_orbit_representative(mono, n_sites))
    end
    return required
end

function _check_linear_state_opt_closure(
    basis::Vector{M},
    n_sites::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    test_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,C<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    n_sites > 0 || throw(DomainError(n_sites, "`n_sites` must be positive."))
    test_width >= 0 || throw(DomainError(test_width, "`test_width` must be non-negative."))
    isempty(basis) && throw(ArgumentError("closure_check(:linear_state_opt, ...) requires a non-empty basis."))
    _check_pauli_chain_support(basis, n_sites; context="closure-check basis")
    _check_pauli_chain_support(hamiltonian, n_sites; context="linear state-opt Hamiltonian")

    generated = _translation_reduced_pair_moment_set(basis, n_sites)
    required = _linear_state_opt_required_moment_set(
        M,
        n_sites,
        hamiltonian,
        test_width;
        sign_symmetry,
    )
    missing = M[m for m in required if !(m in generated)]
    isempty(missing) && return nothing

    sort!(missing)
    shown = join((sprint(show, mono) for mono in Iterators.take(missing, 5)), ", ")
    length(missing) > 5 && (shown *= ", ...")
    throw(ArgumentError(
        "closure_check(:linear_state_opt, ...) missing $(length(missing)) translation-reduced commutator moment(s) for test_width=$test_width: [$shown]. Increase the half-basis order before adding linear state-opt rows."
    ))
end

function _contiguous_state_opt_tests(
    ::Type{M},
    n_sites::Int,
    test_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    tests = Set{M}()
    word = T[]
    for width in 1:test_width
        for code0 in 0:(3^width - 1)
            empty!(word)
            code = code0
            for offset in 0:(width - 1)
                pauli_type = mod(code, 3)
                code = div(code, 3)
                push!(word, _pauli_index_from_site_type(T, offset + 1, pauli_type))
            end
            simplified, phase = simplify(PauliAlgebra, word)
            phase == 0x00 || throw(ArgumentError(
                "Internal error: state-opt test enumeration produced non-real phase $phase."
            ))
            mono = M(copy(simplified))
            sign_symmetry && _pauli_sign_signature(mono) != 0x00 && continue
            push!(tests, _translation_orbit_representative(mono, n_sites))
        end
    end
    return sort!(collect(tests))
end

_closed_int_range(lo::Int, hi::Int) = lo <= hi ? (lo:hi) : (1:0)

function _push_qmbcertify_linear_state_opt_test!(
    tests::Set{M},
    raw_word::Vector{Int},
    sign_symmetry::Bool,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    word = T.(sort(raw_word))
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError(
        "Internal error: QMBCertify linear state-opt test produced non-real phase $phase."
    ))
    mono = M(copy(simplified))
    sign_symmetry && _pauli_sign_signature(mono) != 0x00 && return nothing
    push!(tests, mono)
    return nothing
end

function _qmbcertify_linear_state_opt_tests(
    ::Type{M},
    n_sites::Int,
    test_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    tests = Set{M}()
    b1 = fld(n_sites, 2)
    b2 = test_width

    for i in _closed_int_range(2, n_sites - 1), j in _closed_int_range(i + 1, n_sites)
        _push_qmbcertify_linear_state_opt_test!(
            tests,
            [1, 3 * (i - 1) + 2, 3 * (j - 1) + 3],
            sign_symmetry,
        )
    end

    for i in _closed_int_range(2, b1 - 1),
        j in _closed_int_range(i + 1, b1),
        k in _closed_int_range(2, b1 - 1),
        l in _closed_int_range(k + 1, b1)

        length(unique([i, j, k, l])) == 4 || continue
        for axis in 1:3
            _push_qmbcertify_linear_state_opt_test!(
                tests,
                [
                    1,
                    3 * (i - 1) + 2,
                    3 * (j - 1) + 3,
                    3 * (k - 1) + axis,
                    3 * (l - 1) + axis,
                ],
                sign_symmetry,
            )
        end
    end

    for i in _closed_int_range(2, b2 - 1),
        j in _closed_int_range(i + 1, b2),
        k in _closed_int_range(2, b2 - 3),
        l in _closed_int_range(k + 1, b2 - 2),
        u in _closed_int_range(l + 1, b2 - 1),
        v in _closed_int_range(u + 1, b2)

        length(unique([i, j, k, l, u, v])) == 6 || continue
        for axis in 1:3
            _push_qmbcertify_linear_state_opt_test!(
                tests,
                [
                    1,
                    3 * (i - 1) + 2,
                    3 * (j - 1) + 3,
                    3 * (k - 1) + axis,
                    3 * (l - 1) + axis,
                    3 * (u - 1) + axis,
                    3 * (v - 1) + axis,
                ],
                sign_symmetry,
            )
        end
    end

    for i in _closed_int_range(2, b2 - 1),
        j in _closed_int_range(i + 1, b2),
        k in _closed_int_range(2, b2 - 1),
        l in _closed_int_range(k + 1, b2),
        u in _closed_int_range(2, b2 - 1),
        v in _closed_int_range(u + 1, b2)

        length(unique([i, j, k, l, u, v])) == 6 || continue
        for axis1 in 1:2, axis2 in (axis1 + 1):3
            _push_qmbcertify_linear_state_opt_test!(
                tests,
                [
                    1,
                    3 * (i - 1) + 2,
                    3 * (j - 1) + 3,
                    3 * (k - 1) + axis1,
                    3 * (l - 1) + axis1,
                    3 * (u - 1) + axis2,
                    3 * (v - 1) + axis2,
                ],
                sign_symmetry,
            )
        end
    end

    return sort!(collect(tests))
end

function _linear_state_opt_tests(
    ::Type{M},
    n_sites::Int,
    test_width::Int;
    sign_symmetry::Bool,
    mode::Symbol,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    mode == :contiguous && return _contiguous_state_opt_tests(
        M,
        n_sites,
        test_width;
        sign_symmetry,
    )
    mode == :qmbcertify && return _qmbcertify_linear_state_opt_tests(
        M,
        n_sites,
        test_width;
        sign_symmetry,
    )
    throw(ArgumentError("Unsupported linear_state_opt_mode=$mode."))
end

function _linear_state_opt_required_moment_set(
    ::Type{M},
    n_sites::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    test_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,C<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    required = Set{M}()
    for test_mono in _contiguous_state_opt_tests(M, n_sites, test_width; sign_symmetry)
        commutator = hamiltonian * test_mono - test_mono * hamiltonian
        for (coef, mono) in commutator
            iszero(coef) && continue
            push!(required, _translation_orbit_representative(mono, n_sites))
        end
    end
    return required
end

function _check_psd_state_opt_closure(
    basis::Vector{M},
    n_sites::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    test_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,C<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    n_sites > 0 || throw(DomainError(n_sites, "`n_sites` must be positive."))
    test_width >= 0 || throw(DomainError(test_width, "`test_width` must be non-negative."))
    isempty(basis) && throw(ArgumentError("closure_check(:psd_state_opt, ...) requires a non-empty basis."))
    _check_pauli_chain_support(basis, n_sites; context="closure-check basis")
    _check_pauli_chain_support(hamiltonian, n_sites; context="PSD state-opt Hamiltonian")

    generated = _translation_reduced_pair_moment_set(basis, n_sites)
    required = _psd_state_opt_required_moment_set(
        M,
        n_sites,
        hamiltonian,
        test_width;
        sign_symmetry,
    )
    missing = M[m for m in required if !(m in generated)]
    isempty(missing) && return nothing

    sort!(missing)
    shown = join((sprint(show, mono) for mono in Iterators.take(missing, 5)), ", ")
    length(missing) > 5 && (shown *= ", ...")
    throw(ArgumentError(
        "closure_check(:psd_state_opt, ...) missing $(length(missing)) translation-reduced PSD state-opt moment(s) for test_width=$test_width: [$shown]. Increase the half-basis order before adding PSD state-opt blocks."
    ))
end

function _psd_state_opt_required_moment_set(
    ::Type{M},
    n_sites::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    test_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,C<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    required = Set{M}()
    tests = _contiguous_state_opt_tests(M, n_sites, test_width; sign_symmetry)
    translated = Dict{M,Vector{M}}()
    for w in tests
        translated[w] = [_translate_pauli_monomial(w, r, n_sites) for r in 0:(n_sites - 1)]
    end
    term_cache = _translation_psd_state_opt_term_cache(
        tests,
        n_sites,
        hamiltonian,
        translated,
    )
    union!(required, _translation_psd_state_opt_required_moments(term_cache))
    return required
end

function _pauli_chain_support_width(
    mono::NormalMonomial{PauliAlgebra},
    n_sites::Integer,
)
    isempty(mono.word) && return 0
    n = Int(n_sites)
    sites = sort!(unique(Int[_pauli_site(idx) for idx in mono.word]))
    length(sites) == 1 && return 1

    max_gap = 0
    for idx in 1:length(sites)
        current = sites[idx]
        next = idx == length(sites) ? sites[1] + n : sites[idx + 1]
        max_gap = max(max_gap, next - current)
    end
    return n - max_gap + 1
end

function _pauli_chain_polynomial_support_width(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer,
) where {T<:Unsigned,C<:Number}
    width = 0
    for (_, mono) in poly.terms
        width = max(width, _pauli_chain_support_width(mono, n_sites))
    end
    return width
end

function _check_pauli_translation_degree_closure_limits(
    order::Int,
    n_sites::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    rdm_ks::Vector{Int},
    rdm_support::Symbol,
    linear_state_width::Int,
    psd_state_width::Int,
) where {T<:Unsigned,C<:Number}
    max_rdm_k = min(n_sites, 2order)
    if rdm_support == :closed
        oversized = [k for k in rdm_ks if k > max_rdm_k]
        isempty(oversized) || throw(ArgumentError(
            "`contiguous_rdm_k` with `contiguous_rdm_support=:closed` requires k <= 2order " *
            "(max $max_rdm_k for order=$order, n_sites=$n_sites); got $(oversized). " *
            "Increase `order` or use `contiguous_rdm_support=:extend` intentionally."
        ))
    end

    hamiltonian_width = _pauli_chain_polynomial_support_width(hamiltonian, n_sites)
    max_linear_width = max(0, 2order - hamiltonian_width + 1)
    linear_state_width <= max_linear_width || throw(ArgumentError(
        "`linear_state_opt_width=$linear_state_width` is not closed by order=$order " *
        "for Hamiltonian support width $hamiltonian_width; maximum closed width is " *
        "$max_linear_width. Increase `order` before adding linear state-opt rows."
    ))

    max_psd_width = max(0, fld(2order - hamiltonian_width, 2))
    psd_state_width <= max_psd_width || throw(ArgumentError(
        "`psd_state_opt_width=$psd_state_width` is not closed by order=$order " *
        "for Hamiltonian support width $hamiltonian_width; maximum closed width is " *
        "$max_psd_width. Increase `order` before adding PSD state-opt blocks."
    ))
    return nothing
end

function _check_pauli_translation_constraint_closure(
    basis::Union{Nothing,Vector{NormalMonomial{PauliAlgebra,T}}},
    order::Int,
    n_sites::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    rdm_ks::Vector{Int},
    rdm_support::Symbol,
    linear_state_width::Int,
    psd_state_width::Int;
    sign_symmetry::Bool,
) where {T<:Unsigned,C<:Number}
    if basis === nothing
        return _check_pauli_translation_degree_closure_limits(
            order,
            n_sites,
            hamiltonian,
            rdm_ks,
            rdm_support,
            linear_state_width,
            psd_state_width,
        )
    end

    if rdm_support == :closed
        for k in rdm_ks
            closure_check(:contiguous_rdm, basis; n_sites, k)
        end
    end
    linear_state_width > 0 && closure_check(
        :linear_state_opt,
        basis;
        n_sites,
        hamiltonian,
        test_width=linear_state_width,
        sign_symmetry,
    )
    psd_state_width > 0 && closure_check(
        :psd_state_opt,
        basis;
        n_sites,
        hamiltonian,
        test_width=psd_state_width,
        sign_symmetry,
    )
    return nothing
end

"""
    heisenberg_chain_hamiltonian((σx, σy, σz); coupling=1/4, periodic=true)

Construct the spin-1/2 XXX Heisenberg chain Hamiltonian
`coupling * Σᵢ (σxᵢσxᵢ₊₁ + σyᵢσyᵢ₊₁ + σzᵢσzᵢ₊₁)`.

Use `coupling=1/4` for the usual `S_i · S_j` normalization.
"""
function heisenberg_chain_hamiltonian(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector};
    coupling::Number=1 / 4,
    periodic::Bool=true,
)
    σx, σy, σz, n = _validate_pauli_chain_ops(ops)
    n > 1 || throw(ArgumentError("Heisenberg chain needs at least two sites; got $n."))
    nbonds = periodic ? n : n - 1
    h = zero(coupling * one(first(σx)))
    axes = (σx, σy, σz)
    for i in 1:nbonds
        j = periodic ? mod1(i + 1, n) : i + 1
        for op in axes
            h += coupling * op[i] * op[j]
        end
    end
    return h
end

"""
    pauli_chain_translation((σx, σy, σz); shift=1)

Return the Clifford symmetry translating a periodic Pauli chain by `shift` sites.
"""
function pauli_chain_translation(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector};
    shift::Integer=1,
)
    σx, σy, σz, n = _validate_pauli_chain_ops(ops)
    M = eltype(σx)
    s = mod(Int(shift), n)
    images = Dict{M,Tuple{Int,M}}()
    for op in (σx, σy, σz), i in 1:n
        images[op[i]] = (1, op[mod1(i + s, n)])
    end
    return CliffordSymmetry(images; nqubits=n)
end

"""
    pauli_chain_reflection((σx, σy, σz))

Return the Clifford symmetry reversing a Pauli chain: site `i ↦ N + 1 - i`.
"""
function pauli_chain_reflection(ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector})
    σx, σy, σz, n = _validate_pauli_chain_ops(ops)
    M = eltype(σx)
    images = Dict{M,Tuple{Int,M}}()
    for op in (σx, σy, σz), i in 1:n
        images[op[i]] = (1, op[n + 1 - i])
    end
    return CliffordSymmetry(images; nqubits=n)
end

"""
    pauli_global_axis_rotation_generators((σx, σy, σz))

Return two global Clifford generators for the octahedral Pauli-axis rotation
group: a global Hadamard (`x ↔ z`, `y ↦ -y`) and a global phase rotation
(`x ↦ y`, `y ↦ -x`).
"""
function pauli_global_axis_rotation_generators(ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector})
    σx, σy, σz, n = _validate_pauli_chain_ops(ops)
    M = eltype(σx)

    global_h = Dict{M,Tuple{Int,M}}()
    global_s = Dict{M,Tuple{Int,M}}()
    for i in 1:n
        global_h[σx[i]] = (1, σz[i])
        global_h[σz[i]] = (1, σx[i])
        global_h[σy[i]] = (-1, σy[i])

        global_s[σx[i]] = (1, σy[i])
        global_s[σy[i]] = (-1, σx[i])
    end

    return CliffordSymmetry[
        CliffordSymmetry(global_h; nqubits=n),
        CliffordSymmetry(global_s; nqubits=n),
    ]
end

"""
    heisenberg_chain_symmetry_spec((σx, σy, σz); translation=true, reflection=true, sign=false, axis_rotations=true, ...)

Build a `SymmetrySpec` for a periodic Heisenberg XXX chain from the common
structural generators: lattice translation, reflection, optional global
Heisenberg sign flips, and global Pauli-axis rotations.  This is a convenience
wrapper around `CliffordSymmetry`.
"""
function heisenberg_chain_symmetry_spec(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector};
    translation::Bool=true,
    reflection::Bool=true,
    sign::Bool=false,
    axis_rotations::Bool=true,
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    generators = CliffordSymmetry[]
    translation && push!(generators, pauli_chain_translation(ops))
    reflection && push!(generators, pauli_chain_reflection(ops))
    if sign
        σx, _, _, n = _validate_pauli_chain_ops(ops)
        T = eltype(first(σx).word)
        push!(generators, pauli_sign_symmetry(n; integer_type=T))
    end
    axis_rotations && append!(generators, pauli_global_axis_rotation_generators(ops))
    isempty(generators) && throw(ArgumentError("At least one Heisenberg chain symmetry generator must be enabled."))
    return SymmetrySpec(generators; check_invariance, offblock_check)
end

function _count_clifford_generator(
    generators::AbstractVector{<:CliffordSymmetry},
    target::CliffordSymmetry,
)
    return count(generator -> generator == target, generators)
end

"""
    pauli_chain_fast_path_profile((σx, σy, σz), symmetry)

Conservatively inspect a Pauli-chain `SymmetrySpec` and report whether its
finite Clifford generators match the specialized translation-invariant backend.
This is a recognition helper only; it does not route `cs_nctssos` to the fast
path.
"""
function pauli_chain_fast_path_profile(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    symmetry::SymmetrySpec,
)
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    T = eltype(first(σx).word)
    generators = symmetry.clifford_generators

    translation_count = _count_clifford_generator(generators, pauli_chain_translation(ops))
    reflection_count = _count_clifford_generator(generators, pauli_chain_reflection(ops))
    sign_count = _count_clifford_generator(generators, pauli_sign_symmetry(n; integer_type=T))
    axis_generators = pauli_global_axis_rotation_generators(ops)
    axis_generator_counts = [
        _count_clifford_generator(generators, generator) for generator in axis_generators
    ]
    axis_count = sum(axis_generator_counts)
    axis_rotation_symmetry = all(count -> count > 0, axis_generator_counts)

    recognized = Set{CliffordSymmetry}()
    push!(recognized, pauli_chain_translation(ops))
    push!(recognized, pauli_chain_reflection(ops))
    push!(recognized, pauli_sign_symmetry(n; integer_type=T))
    union!(recognized, axis_generators)
    unrecognized_count = count(generator -> !(generator in recognized), generators)

    unsupported = Symbol[]
    isempty(symmetry.generators) || push!(unsupported, :signed_permutation)
    isempty(symmetry.fermionic_generators) || push!(unsupported, :fermionic_permutation)
    isnothing(symmetry.sector) || push!(unsupported, :fermionic_sector)
    isnothing(symmetry.spin_adaptation) || push!(unsupported, :fermionic_spin_adaptation)
    isnothing(symmetry.pauli_charge) || push!(unsupported, :pauli_charge)
    isnothing(symmetry.pauli_singlet) || push!(unsupported, :pauli_singlet)
    axis_supported = axis_rotation_symmetry
    (axis_count == 0 || axis_supported) || push!(unsupported, :axis_rotation)
    unrecognized_count == 0 || push!(unsupported, :unrecognized_clifford)

    missing = Symbol[]
    translation_count > 0 || push!(missing, :translation)

    return (
        n_sites=n,
        translation_symmetry=translation_count > 0,
        reflection_symmetry=reflection_count > 0,
        sign_symmetry=sign_count > 0,
        axis_rotation_symmetry=axis_rotation_symmetry,
        translation_generator_count=translation_count,
        reflection_generator_count=reflection_count,
        sign_generator_count=sign_count,
        axis_rotation_generator_count=axis_count,
        unrecognized_clifford_generator_count=unrecognized_count,
        supported_by_translation_fast_path=isempty(unsupported) && isempty(missing),
        unsupported_features=unsupported,
        missing_required_features=missing,
        check_invariance=symmetry.check_invariance,
        offblock_check=symmetry.offblock_check,
    )
end

"""
    PauliSU2RDMBlock

Local SU(2) block metadata for a `k`-qubit RDM.  `spin2` stores `2j`,
`multiplicity` is the Schur-Weyl multiplicity of the spin-`j` irrep, and
`irrep_dimension == spin2 + 1`.
"""
struct PauliSU2RDMBlock
    spin2::Int
    multiplicity::Int
    irrep_dimension::Int

    function PauliSU2RDMBlock(spin2::Integer, multiplicity::Integer)
        s2 = Int(spin2)
        mult = Int(multiplicity)
        s2 >= 0 || throw(DomainError(spin2, "`spin2` must be non-negative."))
        mult > 0 || throw(DomainError(multiplicity, "`multiplicity` must be positive."))
        return new(s2, mult, s2 + 1)
    end
end

"""
    PauliSU2OperatorBlock

SU(2) block metadata for the Pauli operator space on a local support.  Each
site contributes `spin-0 ⊕ spin-1`, corresponding to `I ⊕ span(X,Y,Z)`.
"""
struct PauliSU2OperatorBlock
    spin2::Int
    multiplicity::Int
    irrep_dimension::Int

    function PauliSU2OperatorBlock(spin2::Integer, multiplicity::Integer)
        s2 = Int(spin2)
        mult = Int(multiplicity)
        s2 >= 0 || throw(DomainError(spin2, "`spin2` must be non-negative."))
        mult > 0 || throw(DomainError(multiplicity, "`multiplicity` must be positive."))
        return new(s2, mult, s2 + 1)
    end
end

"""
    PauliSU2WordBlock

SU(2) block metadata for Pauli words on a fixed active support.  Every active
letter is non-identity and transforms as spin 1.
"""
struct PauliSU2WordBlock
    spin2::Int
    multiplicity::Int
    irrep_dimension::Int

    function PauliSU2WordBlock(spin2::Integer, multiplicity::Integer)
        s2 = Int(spin2)
        mult = Int(multiplicity)
        s2 >= 0 || throw(DomainError(spin2, "`spin2` must be non-negative."))
        mult > 0 || throw(DomainError(multiplicity, "`multiplicity` must be positive."))
        return new(s2, mult, s2 + 1)
    end
end

"""
    PauliSU2BasisBlock

SU(2) block metadata for a support-complete Pauli word basis.  `logical_row_labels`
stores one provenance label per multiplicity row that survives the
irreducible-tensor reduction.  The coefficient-domain fields describe the
Clebsch-Gordan transform family needed to realize the block.
"""
struct PauliSU2BasisBlock
    spin2::Int
    multiplicity::Int
    irrep_dimension::Int
    logical_row_labels::Vector{Any}
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol

    function PauliSU2BasisBlock(
        spin2::Integer,
        logical_row_labels::AbstractVector;
        coefficient_domain::Symbol=:complex_algebraic_float64,
        exact_coefficient_domain::Symbol=:complex_sqrt_rational,
    )
        s2 = Int(spin2)
        labels = Any[label for label in logical_row_labels]
        s2 >= 0 || throw(DomainError(spin2, "`spin2` must be non-negative."))
        !isempty(labels) || throw(ArgumentError("Pauli SU(2) basis blocks need at least one logical row label."))
        return new(
            s2,
            length(labels),
            s2 + 1,
            labels,
            coefficient_domain,
            exact_coefficient_domain,
        )
    end
end

"""
    PauliSU2BasisTransformBlock

Concrete SU(2) row transform for one total-spin sector of a
support-complete Pauli word basis.  `row_labels` stores one provenance label
per transformed row, including the magnetic sublevel.
"""
struct PauliSU2BasisTransformBlock
    spin2::Int
    multiplicity::Int
    irrep_dimension::Int
    row_labels::Vector{Any}
    transform::Matrix{ComplexF64}
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol

    function PauliSU2BasisTransformBlock(
        spin2::Integer,
        row_labels::AbstractVector,
        transform::AbstractMatrix{<:Complex};
        coefficient_domain::Symbol=:complex_algebraic_float64,
        exact_coefficient_domain::Symbol=:complex_sqrt_rational,
    )
        s2 = Int(spin2)
        labels = Any[label for label in row_labels]
        mat = ComplexF64.(transform)
        s2 >= 0 || throw(DomainError(spin2, "`spin2` must be non-negative."))
        !isempty(labels) || throw(ArgumentError("Pauli SU(2) basis transform blocks need at least one row label."))
        length(labels) == size(mat, 1) || throw(DimensionMismatch(
            "Pauli SU(2) basis transform block has $(length(labels)) row labels but $(size(mat, 1)) transform rows."
        ))
        irrep_dimension = s2 + 1
        length(labels) % irrep_dimension == 0 || throw(ArgumentError(
            "Pauli SU(2) basis transform block for spin2=$s2 has a non-integral multiplicity."
        ))
        return new(
            s2,
            length(labels) ÷ irrep_dimension,
            irrep_dimension,
            labels,
            mat,
            coefficient_domain,
            exact_coefficient_domain,
        )
    end
end

struct PauliSU2BasisMomentOrigin <: BlockOrigin
    label::Any
    logical_row_labels::Vector{Any}
    transform::Any

    function PauliSU2BasisMomentOrigin(
        label,
        logical_row_labels::AbstractVector;
        transform=nothing,
    )
        return new(label, Any[row_label for row_label in logical_row_labels], transform)
    end
end

struct PauliSU2SingletChannelEqualityOriginSeed
    label::Any
end

struct PauliSU2SingletChannelEqualityOrigin <: ConstraintOrigin
    label::Any
    row::Int
    col::Int
    part::Symbol

    function PauliSU2SingletChannelEqualityOrigin(
        label,
        row::Integer,
        col::Integer,
        part::Symbol,
    )
        _assert_positive_index("Pauli SU(2) singlet-channel equality row", row)
        _assert_positive_index("Pauli SU(2) singlet-channel equality col", col)
        part in (:real, :imag, :scalar) || throw(ArgumentError(
            "Unsupported Pauli SU(2) singlet-channel equality part $(repr(part)); " *
            "expected :real, :imag, or :scalar"
        ))
        return new(label, Int(row), Int(col), part)
    end
end

function _zero_constraint_origin(
    seed::PauliSU2SingletChannelEqualityOriginSeed,
    constraint_idx::Integer,
    row::Integer,
    col::Integer,
    part::Symbol,
)
    return PauliSU2SingletChannelEqualityOrigin(seed.label, row, col, part)
end

"""
    pauli_su2_rdm_blocks(k)

Return the SU(2) Schur-Weyl block summary for a `k`-qubit local RDM.  The
result is ordered by increasing `spin2`; PSD reduction keeps one
`multiplicity × multiplicity` block per spin sector.
"""
function pauli_su2_rdm_blocks(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    blocks = PauliSU2RDMBlock[]
    for spin2 in (kk % 2):2:kk
        rank = (kk - spin2) ÷ 2
        previous = rank == 0 ? 0 : binomial(kk, rank - 1)
        multiplicity = binomial(kk, rank) - previous
        push!(blocks, PauliSU2RDMBlock(spin2, multiplicity))
    end
    return blocks
end

_pauli_su2_real_psd_block_size(block) = block.multiplicity

_pauli_su2_real_psd_block_size(block::PauliSU2RDMBlock) =
    block.multiplicity > 1 ? 2 * block.multiplicity : block.multiplicity

_pauli_su2_block_coefficient_domain(block) =
    hasproperty(block, :coefficient_domain) ? block.coefficient_domain : nothing

_pauli_su2_block_coefficient_domain(::PauliSU2RDMBlock) =
    :algebraic_float64

_pauli_su2_block_coefficient_domain(::PauliSU2WordBlock) =
    :complex_algebraic_float64

_pauli_su2_block_exact_coefficient_domain(block) =
    hasproperty(block, :exact_coefficient_domain) ? block.exact_coefficient_domain : nothing

_pauli_su2_block_exact_coefficient_domain(::PauliSU2RDMBlock) =
    :sqrt_rational

_pauli_su2_block_exact_coefficient_domain(::PauliSU2WordBlock) =
    :complex_sqrt_rational

function _pauli_su2_storage_metrics(blocks, full_side::Integer, scalar_bytes::Integer)
    side = Int(full_side)
    side >= 0 || throw(DomainError(full_side, "`full_side` must be non-negative."))
    bytes = Int(scalar_bytes)
    bytes > 0 || throw(ArgumentError("`scalar_bytes` must be positive; got $scalar_bytes."))

    full_dense_entries = side * side
    full_symmetric_entries = side * (side + 1) ÷ 2
    reduced_dense_entries = sum(block.multiplicity^2 for block in blocks; init=0)
    reduced_symmetric_entries = sum(
        block.multiplicity * (block.multiplicity + 1) ÷ 2 for block in blocks; init=0
    )
    real_psd_block_sizes = [_pauli_su2_real_psd_block_size(block) for block in blocks]
    real_psd_dense_entries = sum(size * size for size in real_psd_block_sizes; init=0)
    real_psd_symmetric_entries = sum(
        size * (size + 1) ÷ 2 for size in real_psd_block_sizes; init=0
    )
    transformed_block_sizes = [
        block.irrep_dimension * block.multiplicity for block in blocks
    ]
    block_coefficient_domains = [_pauli_su2_block_coefficient_domain(block) for block in blocks]
    block_exact_coefficient_domains = [
        _pauli_su2_block_exact_coefficient_domain(block) for block in blocks
    ]

    return (
        full_side=side,
        reduced_block_sizes=[block.multiplicity for block in blocks],
        transformed_block_sizes=transformed_block_sizes,
        block_coefficient_domains=block_coefficient_domains,
        block_exact_coefficient_domains=block_exact_coefficient_domains,
        block_coefficient_domain_histogram=_value_histogram_pairs(block_coefficient_domains),
        block_exact_coefficient_domain_histogram=_value_histogram_pairs(block_exact_coefficient_domains),
        reduced_max_block=isempty(blocks) ? 0 : maximum(block.multiplicity for block in blocks),
        transformed_max_block=isempty(transformed_block_sizes) ? 0 : maximum(transformed_block_sizes),
        reduced_total_block_side=sum(block.multiplicity for block in blocks; init=0),
        transformed_total_block_side=sum(transformed_block_sizes; init=0),
        full_dense_entries=full_dense_entries,
        full_symmetric_entries=full_symmetric_entries,
        reduced_dense_entries=reduced_dense_entries,
        reduced_symmetric_entries=reduced_symmetric_entries,
        real_psd_block_sizes=real_psd_block_sizes,
        real_psd_max_block=isempty(real_psd_block_sizes) ? 0 : maximum(real_psd_block_sizes),
        real_psd_total_block_side=sum(real_psd_block_sizes; init=0),
        real_psd_dense_entries=real_psd_dense_entries,
        real_psd_symmetric_entries=real_psd_symmetric_entries,
        full_dense_bytes=full_dense_entries * bytes,
        full_symmetric_bytes=full_symmetric_entries * bytes,
        reduced_dense_bytes=reduced_dense_entries * bytes,
        reduced_symmetric_bytes=reduced_symmetric_entries * bytes,
        real_psd_dense_bytes=real_psd_dense_entries * bytes,
        real_psd_symmetric_bytes=real_psd_symmetric_entries * bytes,
    )
end

function _pauli_su2_reduction_accounting(metrics)
    active_dense_entries = sum(
        block.irrep_dimension * block.multiplicity^2 for block in metrics.blocks; init=0
    )
    offblock_entry_count = metrics.full_dense_entries - active_dense_entries
    copy_entry_count = active_dense_entries - metrics.reduced_dense_entries
    accounted_entry_count =
        offblock_entry_count + copy_entry_count + metrics.reduced_dense_entries

    return (
        active_dense_entries=active_dense_entries,
        offblock_entry_count=offblock_entry_count,
        copy_entry_count=copy_entry_count,
        accounted_entry_count=accounted_entry_count,
    )
end

"""
    pauli_su2_rdm_metrics(k; scalar_bytes=sizeof(Float64))

Return storage-size estimates for a `k`-qubit local RDM before and after the
Schur-Weyl SU(2) split.  The reduced estimates count one multiplicity PSD block
per total-spin sector; `real_psd_*` fields report the solver-facing real PSD
storage convention used by the translation/RDM target path.
"""
function pauli_su2_rdm_metrics(k::Integer; scalar_bytes::Integer=sizeof(Float64))
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    blocks = pauli_su2_rdm_blocks(kk)
    return (blocks=blocks, _pauli_su2_storage_metrics(blocks, 1 << kk, scalar_bytes)...)
end

"""
    pauli_su2_rdm_reduction_diagnostics(k; scalar_bytes=sizeof(Float64))

Return entry-count diagnostics for the Schur-Weyl RDM reduction.  Off-block
entries are forbidden by SU(2) selection rules, copy entries equate equivalent
magnetic-sublevel multiplicity blocks, and reduced entries remain in the PSD
multiplicity matrices.
"""
function pauli_su2_rdm_reduction_diagnostics(
    k::Integer;
    scalar_bytes::Integer=sizeof(Float64),
)
    metrics = pauli_su2_rdm_metrics(k; scalar_bytes)
    return (; metrics..., _pauli_su2_reduction_accounting(metrics)...)
end

@inline _qmbcertify_smod(i::Integer, s::Integer) = iszero(mod(Int(i), Int(s))) ? Int(s) : mod(Int(i), Int(s))
@inline _qmbcertify_site(label::Integer) = cld(Int(label), 3)
@inline _qmbcertify_axis(label::Integer) = _qmbcertify_smod(label, 3)

function _qmbcertify_rot(label::Integer)
    lbl = Int(label)
    lbl == 1 && return (2, 3)
    lbl == 2 && return (3, 1)
    lbl == 3 && return (1, 2)
    throw(ArgumentError("QMBCertify axis label must be 1, 2, or 3; got $label."))
end

function _qmbcertify_chain_source_basis_words(
    n_sites::Int,
    label::Int,
    order::Int;
    extra::Int,
    three_type::Tuple{Int,Int},
)
    words = Vector{UInt16}[]
    n3 = 3 * n_sites
    sorted_word(values...) = sort!(UInt16[Int(value) for value in values])

    if label > 0
        rot1, rot2 = _qmbcertify_rot(label)
        a1 = ((rot1, rot2), (rot2, rot1))
        a2 = (
            (label, label, rot1, rot2),
            (rot2, rot1, label, label),
            (label, rot1, label, rot2),
            (rot2, label, rot1, label),
            (label, rot1, rot2, label),
            (label, rot2, rot1, label),
            (rot1, label, label, rot2),
            (rot2, label, label, rot1),
            (rot1, label, rot2, label),
            (label, rot2, label, rot1),
            (rot1, rot2, label, label),
            (label, label, rot2, rot1),
            (rot1, rot1, rot1, rot2),
            (rot2, rot1, rot1, rot1),
            (rot1, rot1, rot2, rot1),
            (rot1, rot2, rot1, rot1),
            (rot2, rot2, rot2, rot1),
            (rot1, rot2, rot2, rot2),
            (rot2, rot2, rot1, rot2),
            (rot2, rot1, rot2, rot2),
        )

        for site in 1:n_sites
            push!(words, UInt16[3 * (site - 1) + label])
        end
        if order > 2
            for site in 1:n_sites
                push!(
                    words,
                    sorted_word(
                        3 * (site - 1) + label,
                        _qmbcertify_smod(3 * (site - 1 + three_type[1]) + label, n3),
                        _qmbcertify_smod(3 * (site - 1 + sum(three_type)) + label, n3),
                    ),
                )
            end
            for base_axis in (rot1, rot2), changed in 1:3
                axes = [base_axis, base_axis, base_axis]
                axes[changed] = label
                for site in 1:n_sites
                    push!(
                        words,
                        sorted_word(
                            3 * (site - 1) + axes[1],
                            _qmbcertify_smod(3 * (site - 1 + three_type[1]) + axes[2], n3),
                            _qmbcertify_smod(3 * (site - 1 + sum(three_type)) + axes[3], n3),
                        ),
                    )
                end
            end
        end
        for offset in 0:extra, axes in a1, site in 1:n_sites
            push!(
                words,
                sorted_word(
                    3 * (site - 1) + axes[1],
                    _qmbcertify_smod(3 * (site + offset) + axes[2], n3),
                ),
            )
        end
        if order > 3
            for axes in a2, site in 1:n_sites
                push!(
                    words,
                    sorted_word(
                        3 * (site - 1) + axes[1],
                        _qmbcertify_smod(3 * site + axes[2], n3),
                        _qmbcertify_smod(3 * (site + 1) + axes[3], n3),
                        _qmbcertify_smod(3 * (site + 2) + axes[4], n3),
                    ),
                )
            end
        end
        return words
    end

    label == 0 || throw(ArgumentError("QMBCertify chain basis label must be 0, 1, 2, or 3; got $label."))
    a1 = (
        (1, 2, 3),
        (3, 2, 1),
        (1, 3, 2),
        (2, 3, 1),
        (2, 1, 3),
        (3, 1, 2),
    )
    a2 = (
        (1, 1, 1, 1),
        (2, 2, 2, 2),
        (3, 3, 3, 3),
        (1, 2, 2, 1),
        (2, 1, 1, 2),
        (1, 3, 3, 1),
        (3, 1, 1, 3),
        (3, 2, 2, 3),
        (2, 3, 3, 2),
        (1, 1, 2, 2),
        (2, 2, 1, 1),
        (1, 2, 1, 2),
        (2, 1, 2, 1),
        (1, 1, 3, 3),
        (3, 3, 1, 1),
        (1, 3, 1, 3),
        (3, 1, 3, 1),
        (3, 3, 2, 2),
        (2, 2, 3, 3),
        (3, 2, 3, 2),
        (2, 3, 2, 3),
    )
    for offset in 0:extra, axis in 1:3, site in 1:n_sites
        push!(
            words,
            sorted_word(
                3 * (site - 1) + axis,
                _qmbcertify_smod(3 * (site + offset) + axis, n3),
            ),
        )
    end
    if order > 3
        for axes in a2, site in 1:n_sites
            push!(
                words,
                sorted_word(
                    3 * (site - 1) + axes[1],
                    _qmbcertify_smod(3 * site + axes[2], n3),
                    _qmbcertify_smod(3 * (site + 1) + axes[3], n3),
                    _qmbcertify_smod(3 * (site + 2) + axes[4], n3),
                ),
            )
        end
    end
    if order > 2
        for axes in a1, site in 1:n_sites
            push!(
                words,
                sorted_word(
                    3 * (site - 1) + axes[1],
                    _qmbcertify_smod(3 * (site - 1 + three_type[1]) + axes[2], n3),
                    _qmbcertify_smod(3 * (site - 1 + sum(three_type)) + axes[3], n3),
                ),
            )
        end
    end
    return words
end

function _qmbcertify_source_word_monomial(
    ::Type{M},
    source_word::Vector{UInt16},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    word = T.(source_word)
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError(
        "QMBCertify chain source word $source_word simplified with non-real phase $phase."
    ))
    return M(copy(simplified))
end

function _qmbcertify_chain_base_block_sizes(family_counts::AbstractVector{<:Integer}, n_sites::Int)
    return getproperty.(_qmbcertify_chain_base_block_records(family_counts, n_sites), :block_size)
end

function _qmbcertify_chain_base_block_records(
    family_counts::AbstractVector{<:Integer},
    n_sites::Int,
)
    iseven(n_sites) || throw(ArgumentError(
        "QMBCertify chain base block accounting currently requires even n_sites; got $n_sites."
    ))
    length(family_counts) == 2 || throw(ArgumentError(
        "QMBCertify chain base block accounting expects two parity family counts."
    ))
    half = div(n_sites, 2)
    records = NamedTuple[]
    for parity_index in 1:2
        parity = parity_index - 1
        k = Int(family_counts[parity_index])
        for momentum in 0:half
            identity_row = parity == 0 && momentum == 0
            realified = !(momentum == 0 || momentum == half)
            block_size = if identity_row
                k + 1
            elseif realified
                2k
            else
                k
            end
            push!(
                records,
                (
                    parity=parity,
                    momentum=momentum,
                    block_size=block_size,
                    family_count=k,
                    identity_row=identity_row,
                    realified=realified,
                ),
            )
        end
    end
    return records
end

function _qmbcertify_sign_zero_word(word::AbstractVector{<:Integer})
    axis_counts = zeros(Int, 3)
    for label in word
        axis_counts[mod(Int(label), 3) + 1] += 1
    end
    return any(isodd, axis_counts)
end

function _qmbcertify_reduce_chain_support_word(
    raw_word::AbstractVector{<:Integer},
    n_sites::Int,
)
    simplified, _ = simplify(PauliAlgebra, Int.(raw_word))
    return _qmbcertify_reduce_simplified_chain_support_word(simplified, n_sites)
end

function _qmbcertify_reduce_simplified_chain_support_word(
    simplified::AbstractVector{<:Integer},
    n_sites::Int,
)
    _qmbcertify_sign_zero_word(simplified) && return nothing
    return _qmbcertify_reduce4_chain(Int.(simplified), n_sites)
end

function _qmbcertify_chain_support_term(
    ::Type{M},
    simplified::AbstractVector{<:Integer},
    phase_k::UInt8,
    n_sites::Int;
    realify::Bool=false,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    phase_k == 0x04 && return nothing, 0
    coef = realify ?
        _qmbcertify_realified_phase_number(phase_k) :
        _coeff_to_number(PauliAlgebra, phase_k)
    iszero(coef) && return nothing, coef

    reduced = _qmbcertify_reduce_simplified_chain_support_word(simplified, n_sites)
    reduced === nothing && return nothing, zero(coef)
    mono = isempty(reduced) ? one(M) : _unchecked_monomial(PauliAlgebra, T.(collect(reduced)))
    return mono, coef
end

function _qmbcertify_chain_support_term(
    ::Type{M},
    raw_word::AbstractVector{<:Integer},
    n_sites::Int;
    realify::Bool=false,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    simplified, phase_k = simplify(PauliAlgebra, Int.(raw_word))
    return _qmbcertify_chain_support_term(M, simplified, phase_k, n_sites; realify)
end

function _qmbcertify_chain_support_monomial(
    ::Type{M},
    raw_word::AbstractVector{<:Integer},
    n_sites::Int,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    mono, _ = _qmbcertify_chain_support_term(M, raw_word, n_sites)
    return mono
end

function _qmbcertify_reduce_chain_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Int,
    ::Type{P},
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T}}
    PC = _coefficient_type(P)
    M = NormalMonomial{PauliAlgebra,T}
    terms = Tuple{PC,M}[]
    sizehint!(terms, length(poly.terms))
    for (coef, mono) in poly
        reduced, phase_coef = _qmbcertify_chain_support_term(M, mono.word, n_sites)
        reduced === nothing && continue
        value = convert(PC, coef * phase_coef)
        iszero(value) || push!(terms, (value, reduced))
    end
    return convert(P, Polynomial(terms))
end

function _qmbcertify_chain_base_support_metrics(
    source_words_by_parity::AbstractVector,
    family_count_by_parity::AbstractVector{<:Integer},
    n_sites::Int,
)
    support = Set{Tuple{Vararg{Int}}}([()])
    diagonal_support = Set{Tuple{Vararg{Int}}}()
    offdiagonal_support = Set{Tuple{Vararg{Int}}}()
    nonzero_rows = 0
    zero_rows = 0
    diagonal_nonzero_rows = 0
    offdiagonal_nonzero_rows = 0

    for parity in eachindex(source_words_by_parity)
        source_words = source_words_by_parity[parity]
        family_count = Int(family_count_by_parity[parity])
        for family in 1:family_count
            left = source_words[n_sites * (family - 1) + 1]
            for shift in 1:div(n_sites, 2)
                reduced = _qmbcertify_reduce_chain_support_word(
                    vcat(left, source_words[n_sites * (family - 1) + shift + 1]),
                    n_sites,
                )
                if reduced === nothing
                    zero_rows += 1
                else
                    nonzero_rows += 1
                    diagonal_nonzero_rows += 1
                    push!(support, reduced)
                    push!(diagonal_support, reduced)
                end
            end
        end

        for left_family in 1:(family_count - 1), right_family in (left_family + 1):family_count
            left = source_words[n_sites * (left_family - 1) + 1]
            for shift in 1:n_sites
                reduced = _qmbcertify_reduce_chain_support_word(
                    vcat(left, source_words[n_sites * (right_family - 1) + shift]),
                    n_sites,
                )
                if reduced === nothing
                    zero_rows += 1
                else
                    nonzero_rows += 1
                    offdiagonal_nonzero_rows += 1
                    push!(support, reduced)
                    push!(offdiagonal_support, reduced)
                end
            end
        end
    end

    support_words = sort!(collect(support))
    return (
        nonzero_row_count=nonzero_rows,
        zero_row_count=zero_rows,
        diagonal_nonzero_row_count=diagonal_nonzero_rows,
        offdiagonal_nonzero_row_count=offdiagonal_nonzero_rows,
        unique_count=length(support_words),
        diagonal_unique_count=length(diagonal_support),
        offdiagonal_unique_count=length(offdiagonal_support),
        words=support_words,
        word_length_histogram=_histogram_pairs(length.(support_words)),
    )
end

"""
    pauli_qmbcertify_chain_basis((σx, σy, σz), order; extra=0, three_type=(1, 1))

Return source-pinned structural data for QMBCertify's 1D-chain `get_basis`
half-basis.  This mirrors the chain branch used by QMBCertify for labels `0`
and `1`; it is diagnostic/target data, not a promise that the current
translation DFT backend can solve with this basis.  In particular, `extra`
can introduce short translation orbits.
"""
function pauli_qmbcertify_chain_basis(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    extra::Integer=0,
    three_type::Tuple{<:Integer,<:Integer}=(1, 1),
)
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    ext = Int(extra)
    ext >= 0 || throw(ArgumentError("`extra` must be non-negative; got $extra."))
    tt = (Int(three_type[1]), Int(three_type[2]))
    all(>(0), tt) || throw(ArgumentError("`three_type` entries must be positive; got $three_type."))

    M = eltype(σx)
    source_words_by_parity = [
        _qmbcertify_chain_source_basis_words(n, label, d; extra=ext, three_type=tt)
        for label in (0, 1)
    ]
    basis_by_parity = [
        [_qmbcertify_source_word_monomial(M, word) for word in source_words]
        for source_words in source_words_by_parity
    ]
    family_count_by_parity = Int[length(source_words) ÷ n for source_words in source_words_by_parity]
    all(
        length(source_words) == n * family_count
        for (source_words, family_count) in zip(source_words_by_parity, family_count_by_parity)
    ) || throw(ArgumentError("Internal QMBCertify chain basis family count is not translation-regular."))

    family_records = NamedTuple[]
    for parity in 1:2
        for family in 1:family_count_by_parity[parity]
            word = source_words_by_parity[parity][n * (family - 1) + 1]
            mono = basis_by_parity[parity][n * (family - 1) + 1]
            rep = _translation_orbit_representative(mono, n)
            push!(
                family_records,
                (
                    parity=parity - 1,
                    family,
                    word=copy(word),
                    word_length=length(word),
                    monomial=mono,
                    representative=rep,
                    translation_orbit_length=_translation_orbit_length(rep, n),
                ),
            )
        end
    end

    unique_basis = M[one(first(σx))]
    for basis in basis_by_parity
        append!(unique_basis, basis)
    end
    sorted_unique!(unique_basis)
    source_row_count = sum(length, basis_by_parity)
    unique_row_count_by_parity = length.(unique.(basis_by_parity))
    duplicate_row_count_by_parity = length.(basis_by_parity) .- unique_row_count_by_parity
    short_orbit_family_count = count(
        record -> record.translation_orbit_length < n,
        family_records,
    )
    short_orbit_overcomplete_row_count = sum(
        record -> n - record.translation_orbit_length,
        family_records;
        init=0,
    )
    base_block_sizes = _qmbcertify_chain_base_block_sizes(family_count_by_parity, n)
    base_block_records = _qmbcertify_chain_base_block_records(family_count_by_parity, n)
    base_support_metrics = _qmbcertify_chain_base_support_metrics(
        source_words_by_parity,
        family_count_by_parity,
        n,
    )

    return (
        n_sites=n,
        order=d,
        extra=ext,
        three_type=tt,
        source_words_by_parity=source_words_by_parity,
        basis_by_parity=basis_by_parity,
        unique_basis=unique_basis,
        source_row_count=source_row_count,
        unique_nonidentity_row_count=length(unique_basis) - 1,
        duplicate_source_row_count=source_row_count - (length(unique_basis) - 1),
        unique_row_count_by_parity=unique_row_count_by_parity,
        duplicate_row_count_by_parity=duplicate_row_count_by_parity,
        family_count_by_parity=family_count_by_parity,
        family_records=family_records,
        short_orbit_family_count=short_orbit_family_count,
        short_orbit_overcomplete_row_count=short_orbit_overcomplete_row_count,
        base_block_records=base_block_records,
        base_block_sizes=base_block_sizes,
        base_block_histogram=_histogram_pairs(base_block_sizes),
        base_block_count=length(base_block_sizes),
        base_max_block=maximum(base_block_sizes),
        base_total_block_side=sum(base_block_sizes),
        base_dense_entries=sum(size -> size^2, base_block_sizes),
        base_symmetric_entries=sum(size -> size * (size + 1) ÷ 2, base_block_sizes),
        base_support_words=base_support_metrics.words,
        base_support_nonzero_row_count=base_support_metrics.nonzero_row_count,
        base_support_zero_row_count=base_support_metrics.zero_row_count,
        base_support_diagonal_nonzero_row_count=
            base_support_metrics.diagonal_nonzero_row_count,
        base_support_offdiagonal_nonzero_row_count=
            base_support_metrics.offdiagonal_nonzero_row_count,
        base_support_unique_count=base_support_metrics.unique_count,
        base_support_diagonal_unique_count=base_support_metrics.diagonal_unique_count,
        base_support_offdiagonal_unique_count=base_support_metrics.offdiagonal_unique_count,
        base_support_word_length_histogram=base_support_metrics.word_length_histogram,
        family_word_length_histogram=_histogram_pairs(getproperty.(family_records, :word_length)),
        family_translation_orbit_length_histogram=
            _histogram_pairs(getproperty.(family_records, :translation_orbit_length)),
    )
end

function _qmbcertify_axis_permutations(word::Vector{Int})
    return Tuple{Vararg{Int}}[
        Tuple(3 * (_qmbcertify_site(label) - 1) + perm[_qmbcertify_axis(label)] for label in word)
        for perm in _PAULI_AXIS_PERMUTATIONS
    ]
end

function _qmbcertify_reduce4_chain(word::Vector{Int}, ambient_sites::Int)
    isempty(word) && return ()

    candidates = Tuple{Vararg{Int}}[]
    l = length(word)
    for start in 1:l
        shifted = Vector{Int}(undef, l)
        shift = 3 * (_qmbcertify_site(word[start]) - 1)
        idx = 1
        for source in start:l
            shifted[idx] = word[source] - shift
            idx += 1
        end
        for source in 1:(start - 1)
            shifted[idx] = word[source] + 3 * ambient_sites - shift
            idx += 1
        end
        append!(candidates, _qmbcertify_axis_permutations(shifted))

        reflected = reverse(shifted)
        last_site = _qmbcertify_site(shifted[end])
        reflected_shifted = [
            3 * (last_site - _qmbcertify_site(label)) + _qmbcertify_axis(label)
            for label in reflected
        ]
        append!(candidates, _qmbcertify_axis_permutations(reflected_shifted))
    end
    return minimum(candidates)
end

function _qmbcertify_even_axis_word(k::Int, code0::Int)
    counts = zeros(Int, 3)
    word = Int[]
    code = code0
    for site in 1:k
        axis = mod(code, 4)
        code = div(code, 4)
        iszero(axis) && continue
        counts[axis] += 1
        push!(word, 3 * (site - 1) + axis)
    end
    all(iseven, counts) || return nothing
    return word
end

function _qmbcertify_rdm_entry_sign(row0::Int, code0::Int, k::Int)
    sign = 1
    phase = 0
    col0 = row0
    code = code0
    for offset in 0:(k - 1)
        axis = mod(code, 4)
        code = div(code, 4)
        iszero(axis) && continue
        row_bit = (row0 >> offset) & 1
        if axis == 1
            col0 ⊻= 1 << offset
        elseif axis == 2
            col0 ⊻= 1 << offset
            phase = mod(phase + (iszero(row_bit) ? 3 : 1), 4)
        else
            iszero(row_bit) || (sign = -sign)
        end
    end
    phase == 0 || phase == 2 || throw(ArgumentError(
        "Internal error: QMBCertify RDM even-axis word produced non-real matrix entry."
    ))
    phase == 2 && (sign = -sign)
    return col0, sign
end

function _union_find_root!(parent::Vector{Int}, idx::Int)
    root = idx
    while parent[root] != root
        root = parent[root]
    end
    while parent[idx] != idx
        next = parent[idx]
        parent[idx] = root
        idx = next
    end
    return root
end

function _union_find_union!(parent::Vector{Int}, rank::Vector{Int}, a::Int, b::Int)
    root_a = _union_find_root!(parent, a)
    root_b = _union_find_root!(parent, b)
    root_a == root_b && return nothing
    if rank[root_a] < rank[root_b]
        parent[root_a] = root_b
    elseif rank[root_a] > rank[root_b]
        parent[root_b] = root_a
    else
        parent[root_b] = root_a
        rank[root_a] += 1
    end
    return nothing
end

function _int_vector_lexless(left::Vector{Int}, right::Vector{Int})
    for idx in 1:min(length(left), length(right))
        left[idx] == right[idx] && continue
        return left[idx] < right[idx]
    end
    return length(left) < length(right)
end

function _sort_qmbcertify_rdm_blocks!(blocks::Vector{Vector{Int}})
    return sort!(blocks; lt=(left, right) -> begin
        length(left) == length(right) && return _int_vector_lexless(left, right)
        return length(left) > length(right)
    end)
end

function _unique_qmbcertify_rdm_block_representatives(
    blocks::Vector{Vector{Int}},
)
    representatives = Vector{Int}[]
    seen_sizes = Set{Int}()
    for block in blocks
        size = length(block)
        size in seen_sizes && continue
        push!(seen_sizes, size)
        push!(representatives, block)
    end
    return representatives
end

const _QMBCERTIFY_RDM_BLOCK_TEMPLATE_CACHE =
    Dict{Tuple{Int,Int},Vector{Vector{Int}}}()

function _qmbcertify_rdm_component_blocks(
    k::Integer;
    ambient_sites::Integer=100,
)
    kk = Int(k)
    kk > 0 || throw(DomainError(k, "`k` must be positive."))
    ambient = Int(ambient_sites)
    ambient >= kk || throw(ArgumentError(
        "`ambient_sites` must be at least k for QMBCertify RDM diagnostics; got ambient_sites=$ambient, k=$kk."
    ))

    dim = 1 << kk
    edge_coeffs = Dict{Int,Dict{Tuple{Vararg{Int}},Int}}()
    for code0 in 0:(4^kk - 1)
        word = _qmbcertify_even_axis_word(kk, code0)
        word === nothing && continue
        reduced = _qmbcertify_reduce4_chain(word, ambient)
        for row0 in 0:(dim - 1)
            col0, sign = _qmbcertify_rdm_entry_sign(row0, code0, kk)
            row0 < col0 || continue
            edge_key = row0 * dim + col0
            coeffs = get!(edge_coeffs, edge_key) do
                Dict{Tuple{Vararg{Int}},Int}()
            end
            coeffs[reduced] = get(coeffs, reduced, 0) + sign
        end
    end

    parent = collect(1:dim)
    rank = zeros(Int, dim)
    for (edge_key, coeffs) in edge_coeffs
        any(!iszero, values(coeffs)) || continue
        row0 = div(edge_key, dim)
        col0 = mod(edge_key, dim)
        _union_find_union!(parent, rank, row0 + 1, col0 + 1)
    end

    components = Dict{Int,Vector{Int}}()
    for idx in 1:dim
        root = _union_find_root!(parent, idx)
        push!(get!(components, root) do
            Int[]
        end, idx)
    end
    blocks = [sort!(rows) for rows in values(components) if length(rows) > 1]
    _sort_qmbcertify_rdm_blocks!(blocks)
    return _unique_qmbcertify_rdm_block_representatives(blocks)
end

function _qmbcertify_rdm_block_templates(
    k::Integer;
    ambient_sites::Integer=100,
)
    kk = Int(k)
    ambient = Int(ambient_sites)
    return get!(_QMBCERTIFY_RDM_BLOCK_TEMPLATE_CACHE, (kk, ambient)) do
        _qmbcertify_rdm_component_blocks(kk; ambient_sites=ambient)
    end
end

function _copy_qmbcertify_rdm_blocks(blocks::Vector{Vector{Int}})
    return [copy(block) for block in blocks]
end

"""
    pauli_qmbcertify_rdm_blocks(k; ambient_sites=100)

Return the one-based local computational-basis row blocks used by
QMBCertify's `posepsd8!`-style RDM storage for a `k`-site contiguous RDM.

The blocks are derived from QMBCertify's even-axis Pauli filter and chain
`reduce4` equivalence under translations, reflection, and global Pauli-axis
permutations.  Singleton components are omitted, matching QMBCertify's
hardcoded reusable PSD block variables.  Model construction uses these blocks
only on the direct-linear finite-axis quotient path, where the surrounding
moments carry the matching translation/reflection/global-axis identifications.
"""
function pauli_qmbcertify_rdm_blocks(
    k::Integer;
    ambient_sites::Integer=100,
)
    return _copy_qmbcertify_rdm_blocks(
        _qmbcertify_rdm_block_templates(k; ambient_sites),
    )
end

"""
    pauli_qmbcertify_rdm_block_sizes(k; ambient_sites=100)

Derive the local RDM PSD block sizes used by QMBCertify's `posepsd8!`-style
construction.  The diagnostic follows QMBCertify's even-axis Pauli filter and
its chain `reduce4` equivalence under translations, reflection, and global
Pauli-axis permutations, then returns the unique non-singleton connected
component sizes of the generic local RDM sparsity graph.  Repeated components
are represented once, matching QMBCertify's reusable PSD block variables.

The direct-linear finite-axis quotient path can emit these row groups as RDM
PSD blocks; other constructors still use this as structural/reference data.
"""
function pauli_qmbcertify_rdm_block_sizes(
    k::Integer;
    ambient_sites::Integer=100,
)
    kk = Int(k)
    kk > 0 || throw(DomainError(k, "`k` must be positive."))
    ambient = Int(ambient_sites)
    ambient >= kk || throw(ArgumentError(
        "`ambient_sites` must be at least k for QMBCertify RDM diagnostics; got ambient_sites=$ambient, k=$kk."
    ))
    if ambient == 100
        source_pinned = get(_QMBCERTIFY_RDM_REFERENCE_BLOCK_SIZES, kk, nothing)
        source_pinned === nothing || return copy(source_pinned)
    end
    return length.(_qmbcertify_rdm_block_templates(kk; ambient_sites=ambient))
end

"""
    pauli_qmbcertify_rdm_block_metrics(k; ambient_sites=100, scalar_bytes=sizeof(Float64))

Return storage metrics for [`pauli_qmbcertify_rdm_block_sizes`](@ref).
"""
function pauli_qmbcertify_rdm_block_metrics(
    k::Integer;
    ambient_sites::Integer=100,
    scalar_bytes::Integer=sizeof(Float64),
)
    blocks = pauli_qmbcertify_rdm_block_sizes(k; ambient_sites)
    bytes = Int(scalar_bytes)
    bytes > 0 || throw(ArgumentError("`scalar_bytes` must be positive; got $scalar_bytes."))
    dense_entries = sum(size * size for size in blocks; init=0)
    symmetric_entries = sum(size * (size + 1) ÷ 2 for size in blocks; init=0)
    return (
        k=Int(k),
        block_sizes=blocks,
        n_blocks=length(blocks),
        max_block=isempty(blocks) ? 0 : maximum(blocks),
        total_block_side=sum(blocks; init=0),
        dense_entries=dense_entries,
        symmetric_entries=symmetric_entries,
        dense_bytes=dense_entries * bytes,
        symmetric_bytes=symmetric_entries * bytes,
        requires_construction=false,
    )
end

function _su2_coupled_spin2(parent_spin2::Int, local_spin2::Int)
    local_spin2 == 0 && return Int[parent_spin2]
    first_spin2 = abs(parent_spin2 - local_spin2)
    return collect(first_spin2:2:(parent_spin2 + local_spin2))
end

function _su2_tensor_power_multiplicities(local_spins2, k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    counts = Dict{Int,Int}(0 => 1)
    for _ in 1:kk
        next_counts = Dict{Int,Int}()
        for (parent_spin2, multiplicity) in counts
            for local_spin2 in local_spins2
                for child_spin2 in _su2_coupled_spin2(parent_spin2, Int(local_spin2))
                    next_counts[child_spin2] = get(next_counts, child_spin2, 0) + multiplicity
                end
            end
        end
        counts = next_counts
    end
    return counts
end

"""
    pauli_su2_operator_blocks(k)

Return the SU(2) irrep summary for the Pauli operator space on `k` local sites.
The decomposition is `(spin-0 ⊕ spin-1)^{⊗k}`: identities are scalar and each
non-identity Pauli letter transforms as spin 1.
"""
function pauli_su2_operator_blocks(k::Integer)
    counts = _su2_tensor_power_multiplicities((0, 2), k)
    return [
        PauliSU2OperatorBlock(spin2, counts[spin2])
        for spin2 in sort!(collect(keys(counts)))
    ]
end

"""
    pauli_su2_operator_metrics(k; scalar_bytes=sizeof(Float64))

Return storage-size estimates for the Pauli operator space on `k` local sites
before and after the SU(2) irreducible-tensor split.
"""
function pauli_su2_operator_metrics(k::Integer; scalar_bytes::Integer=sizeof(Float64))
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    blocks = pauli_su2_operator_blocks(kk)
    return (blocks=blocks, _pauli_su2_storage_metrics(blocks, 4^kk, scalar_bytes)...)
end

"""
    pauli_su2_operator_reduction_diagnostics(k; scalar_bytes=sizeof(Float64))

Return entry-count diagnostics for the SU(2) irreducible-tensor reduction of
the full Pauli operator space `(spin-0 ⊕ spin-1)^{⊗k}`.
"""
function pauli_su2_operator_reduction_diagnostics(
    k::Integer;
    scalar_bytes::Integer=sizeof(Float64),
)
    metrics = pauli_su2_operator_metrics(k; scalar_bytes)
    return (; metrics..., _pauli_su2_reduction_accounting(metrics)...)
end

const _PAULI_SU2_OPERATOR_LOCAL_STATES = (
    (spin2=0, m2=0),
    (spin2=2, m2=-2),
    (spin2=2, m2=0),
    (spin2=2, m2=2),
)

const _PAULI_SU2_WORD_LOCAL_STATES = (
    (spin2=2, m2=-2),
    (spin2=2, m2=0),
    (spin2=2, m2=2),
)

function _pauli_su2_tensor_dimension(radix::Int, k::Int, max_dimension::Int)
    max_dimension > 0 || throw(ArgumentError(
        "`max_dimension` must be positive; got $max_dimension."
    ))

    dim = 1
    for _ in 1:k
        dim <= max_dimension ÷ radix || throw(ArgumentError(
            "SU(2) tensor diagnostic dimension exceeds max_dimension=$max_dimension."
        ))
        dim *= radix
    end
    return dim
end

function _pauli_su2_tensor_strides(radix::Int, k::Int)
    strides = Vector{Int}(undef, k)
    stride = 1
    for site in 1:k
        strides[site] = stride
        stride *= radix
    end
    return strides
end

function _pauli_su2_raise_table(local_states)
    radix = length(local_states)
    raised_code = fill(0, radix)
    raised_coeff = zeros(Float64, radix)
    for code in 1:radix
        local_state = local_states[code]
        local_state.m2 < local_state.spin2 || continue
        target_m2 = local_state.m2 + 2
        for target in 1:radix
            target_local = local_states[target]
            if target_local.spin2 == local_state.spin2 && target_local.m2 == target_m2
                raised_code[code] = target
                raised_coeff[code] = 0.5 * sqrt(
                    (local_state.spin2 - local_state.m2) *
                    (local_state.spin2 + local_state.m2 + 2),
                )
                break
            end
        end
    end
    return raised_code, raised_coeff
end

function _pauli_su2_tensor_generators(local_states, k::Int, max_dimension::Int)
    radix = length(local_states)
    dim = _pauli_su2_tensor_dimension(radix, k, max_dimension)
    strides = _pauli_su2_tensor_strides(radix, k)
    raised_code, raised_coeff = _pauli_su2_raise_table(local_states)

    m2_by_state = Vector{Int}(undef, dim)
    jplus = zeros(Float64, dim, dim)
    for state in 0:(dim - 1)
        m2 = 0
        for site in 1:k
            stride = strides[site]
            code = (state ÷ stride) % radix + 1
            local_state = local_states[code]
            m2 += local_state.m2

            target_code = raised_code[code]
            iszero(target_code) && continue
            target_state = state + (target_code - code) * stride
            jplus[target_state + 1, state + 1] += raised_coeff[code]
        end
        m2_by_state[state + 1] = m2
    end

    jminus_jplus = transpose(jplus) * jplus
    jz_diag = 0.5 .* m2_by_state
    casimir = Matrix{Float64}(jminus_jplus)
    for idx in 1:dim
        m = jz_diag[idx]
        casimir[idx, idx] += m^2 + m
    end
    return jz_diag, casimir, m2_by_state
end

function _pauli_su2_spin2_from_casimir(value::Real)
    spin2_float = sqrt(max(0.0, 1 + 4 * Float64(value))) - 1
    spin2 = round(Int, spin2_float)
    eigenvalue = 0.25 * spin2 * (spin2 + 2)
    return spin2, abs(Float64(value) - eigenvalue)
end

function _pauli_su2_casimir_eigenbasis(casimir::Matrix{Float64}, m2_by_state::Vector{Int})
    dim = length(m2_by_state)
    states_by_m2 = Dict{Int,Vector{Int}}()
    for (idx, m2) in pairs(m2_by_state)
        push!(get!(states_by_m2, m2, Int[]), idx)
    end

    transform = Matrix{Float64}(undef, dim, dim)
    labels = Vector{Tuple{Int,Int}}(undef, dim)
    spin_identification_residual = 0.0
    col = 1
    for m2 in sort!(collect(keys(states_by_m2)))
        rows = states_by_m2[m2]
        eigen = LinearAlgebra.eigen(LinearAlgebra.Symmetric(casimir[rows, rows]))
        for local_col in axes(eigen.vectors, 2)
            vector = zeros(Float64, dim)
            vector[rows] = eigen.vectors[:, local_col]
            spin2, residual = _pauli_su2_spin2_from_casimir(eigen.values[local_col])
            spin_identification_residual = max(spin_identification_residual, residual)
            transform[:, col] = vector
            labels[col] = (spin2, m2)
            col += 1
        end
    end
    return transform, labels, spin_identification_residual
end

function _pauli_su2_offdiag_residual(mat::Matrix{Float64})
    residual = 0.0
    for col in axes(mat, 2), row in axes(mat, 1)
        row == col && continue
        residual = max(residual, abs(mat[row, col]))
    end
    return residual
end

function _pauli_su2_tensor_spin_diagnostics(local_states, k::Int, max_dimension::Int)
    jz_diag, casimir, m2_by_state =
        _pauli_su2_tensor_generators(local_states, k, max_dimension)
    transform, labels, spin_identification_residual =
        _pauli_su2_casimir_eigenbasis(casimir, m2_by_state)

    gram_residual = transpose(transform) * transform - I
    jz_transform = transpose(transform) * (jz_diag .* transform)
    casimir_transform = transpose(transform) * casimir * transform
    max_sz_residual = 0.0
    max_casimir_residual = spin_identification_residual
    spin_state_counts = Dict{Int,Int}()

    for (col, (spin2, m2)) in pairs(labels)
        vector = transform[:, col]
        m = 0.5 * m2
        max_sz_residual = max(
            max_sz_residual,
            LinearAlgebra.norm(jz_diag .* vector .- m .* vector),
        )
        eigenvalue = 0.25 * spin2 * (spin2 + 2)
        max_casimir_residual = max(
            max_casimir_residual,
            LinearAlgebra.norm(casimir * vector .- eigenvalue .* vector),
        )
        spin_state_counts[spin2] = get(spin_state_counts, spin2, 0) + 1
    end

    spin_multiplicities = Pair{Int,Int}[]
    for spin2 in sort!(collect(keys(spin_state_counts)))
        state_count = spin_state_counts[spin2]
        irrep_dimension = spin2 + 1
        state_count % irrep_dimension == 0 || throw(ArgumentError(
            "Internal SU(2) diagnostic found non-integral multiplicity for spin2=$spin2."
        ))
        push!(spin_multiplicities, spin2 => state_count ÷ irrep_dimension)
    end

    return (
        dimension=size(transform, 1),
        state_count=length(labels),
        spin_multiplicities=spin_multiplicities,
        unitarity_residual=maximum(abs, gram_residual),
        sz_residual=max_sz_residual,
        casimir_residual=max_casimir_residual,
        offblock_residual=max(
            _pauli_su2_offdiag_residual(jz_transform),
            _pauli_su2_offdiag_residual(casimir_transform),
        ),
    )
end

function _pauli_su2_word_local_spherical_transform()
    invsqrt2 = inv(sqrt(2.0))
    return ComplexF64[
        invsqrt2 -im * invsqrt2 0
        0 0 1
        -invsqrt2 -im * invsqrt2 0
    ]
end

function _pauli_su2_tensor_power_matrix(local_transform::Matrix{ComplexF64}, k::Int)
    k == 0 && return ones(ComplexF64, 1, 1)

    result = copy(local_transform)
    for _ in 2:k
        result = kron(result, local_transform)
    end
    return result
end

function _pauli_su2_spin_multiplicities_from_labels(labels::Vector{Tuple{Int,Int}})
    spin_state_counts = Dict{Int,Int}()
    for (spin2, _) in labels
        spin_state_counts[spin2] = get(spin_state_counts, spin2, 0) + 1
    end

    spin_multiplicities = Pair{Int,Int}[]
    for spin2 in sort!(collect(keys(spin_state_counts)))
        irrep_dimension = spin2 + 1
        state_count = spin_state_counts[spin2]
        state_count % irrep_dimension == 0 || throw(ArgumentError(
            "Internal SU(2) transform found non-integral multiplicity for spin2=$spin2."
        ))
        push!(spin_multiplicities, spin2 => state_count ÷ irrep_dimension)
    end
    return spin_multiplicities
end

"""
    pauli_su2_word_transform(support_size; max_dimension=4096)

Return the dense fixed-support transform from Cartesian Pauli-word rows to
SU(2)-coupled tensor rows for `(spin-1)^{⊗support_size}`.  Columns are ordered
by local Cartesian letters `(X, Y, Z)` with the same base-3 convention as
`pauli_contiguous_chain_basis`; rows are ordered by the internally generated
coupled spin basis.
"""
function pauli_su2_word_transform(
    support_size::Integer;
    max_dimension::Integer=4096,
)
    s = Int(support_size)
    s >= 0 || throw(DomainError(support_size, "`support_size` must be non-negative."))

    jz_diag, casimir, m2_by_state =
        _pauli_su2_tensor_generators(_PAULI_SU2_WORD_LOCAL_STATES, s, Int(max_dimension))
    magnetic_to_cartesian = _pauli_su2_tensor_power_matrix(
        _pauli_su2_word_local_spherical_transform(),
        s,
    )
    coupled_to_magnetic, labels, spin_identification_residual =
        _pauli_su2_casimir_eigenbasis(casimir, m2_by_state)
    transform = ComplexF64.(transpose(coupled_to_magnetic)) * magnetic_to_cartesian
    gram_residual = transform * transform' - I

    return (
        support_size=s,
        dimension=size(transform, 1),
        state_count=length(labels),
        labels=[
            (spin2=spin2, m2=m2)
            for (spin2, m2) in labels
        ],
        spin_multiplicities=_pauli_su2_spin_multiplicities_from_labels(labels),
        coefficient_domain=:complex_algebraic_float64,
        exact_coefficient_domain=:complex_sqrt_rational,
        spin_identification_residual=spin_identification_residual,
        unitarity_residual=maximum(abs, gram_residual),
        transform=transform,
    )
end

"""
    pauli_su2_word_singlet_channels(support_size; max_dimension=4096)

Return the spin-0 channel rows of the fixed-support Pauli-word SU(2)
transform.  This exposes the irreducible-tensor singlet channels that underlie
higher-degree SU(2)-invariant scalar equalities.
"""
function pauli_su2_word_singlet_channels(
    support_size::Integer;
    max_dimension::Integer=4096,
)
    word_transform = pauli_su2_word_transform(support_size; max_dimension)
    singlet_rows = findall(
        label -> label.spin2 == 0 && label.m2 == 0,
        word_transform.labels,
    )
    singlet_transform = word_transform.transform[singlet_rows, :]
    singlet_orthonormality_residual = isempty(singlet_rows) ? 0.0 :
        maximum(abs, singlet_transform * adjoint(singlet_transform) - I)

    return (
        support_size=word_transform.support_size,
        dimension=word_transform.dimension,
        channel_count=length(singlet_rows),
        rows=singlet_rows,
        row_labels=word_transform.labels[singlet_rows],
        coefficient_domain=word_transform.coefficient_domain,
        exact_coefficient_domain=word_transform.exact_coefficient_domain,
        unitarity_residual=word_transform.unitarity_residual,
        singlet_orthonormality_residual=singlet_orthonormality_residual,
        transform=singlet_transform,
    )
end

function _check_pauli_su2_atol(atol::Real)
    tolerance = Float64(atol)
    tolerance >= 0 || throw(DomainError(atol, "`atol` must be non-negative."))
    return tolerance
end

function _pauli_su2_transform_column_forms(
    transform::AbstractMatrix{<:Number};
    atol::Real=1e-12,
)
    tolerance = _check_pauli_su2_atol(atol)

    forms = Vector{Vector{Tuple{Int,ComplexF64}}}(undef, size(transform, 1))
    for row_idx in axes(transform, 1)
        form = Tuple{Int,ComplexF64}[]
        for col_idx in axes(transform, 2)
            coefficient = ComplexF64(transform[row_idx, col_idx])
            abs(coefficient) <= tolerance && continue
            push!(form, (Int(col_idx), coefficient))
        end
        forms[Int(row_idx)] = form
    end
    return forms
end

function _pauli_su2_row_orthonormality_residual(transform::AbstractMatrix{<:Number})
    size(transform, 1) == 0 && return 0.0
    return maximum(abs, transform * adjoint(transform) - I)
end

function _pauli_su2_row_cross_residual(
    left::AbstractMatrix{<:Number},
    right::AbstractMatrix{<:Number},
)
    (size(left, 1) == 0 || size(right, 1) == 0) && return 0.0
    return maximum(abs, left * adjoint(right))
end

"""
    pauli_su2_word_singlet_channel_equalities(support_size; max_dimension=4096, atol=1e-12)

Return the non-singlet rows of the fixed-support Pauli-word SU(2) transform.
These rows are the scalar linear forms that must vanish when a Pauli
expectation vector is constrained to lie in the singlet channel.
"""
function pauli_su2_word_singlet_channel_equalities(
    support_size::Integer;
    max_dimension::Integer=4096,
    atol::Real=1e-12,
)
    word_transform = pauli_su2_word_transform(support_size; max_dimension)
    _check_pauli_su2_atol(atol)
    equality_rows = findall(label -> label.spin2 != 0, word_transform.labels)
    singlet_rows = findall(label -> label.spin2 == 0 && label.m2 == 0, word_transform.labels)
    equality_transform = word_transform.transform[equality_rows, :]
    singlet_transform = word_transform.transform[singlet_rows, :]

    return (
        support_size=word_transform.support_size,
        dimension=word_transform.dimension,
        equality_count=length(equality_rows),
        rows=equality_rows,
        row_labels=word_transform.labels[equality_rows],
        coefficient_domain=word_transform.coefficient_domain,
        exact_coefficient_domain=word_transform.exact_coefficient_domain,
        equality_orthonormality_residual=_pauli_su2_row_orthonormality_residual(
            equality_transform,
        ),
        singlet_orthogonality_residual=_pauli_su2_row_cross_residual(
            equality_transform,
            singlet_transform,
        ),
        transform=equality_transform,
        column_forms=_pauli_su2_transform_column_forms(equality_transform; atol),
    )
end

function _pauli_su2_reverse_base3_code(code0::Int, support_size::Int)
    code = code0
    reversed = 0
    for _ in 1:support_size
        digit = mod(code, 3)
        code = div(code, 3)
        reversed = 3 * reversed + digit
    end
    return reversed
end

function _pauli_su2_word_reversal_matrix(support_size::Int)
    dim = 3^support_size
    mat = zeros(ComplexF64, dim, dim)
    for code in 0:(dim - 1)
        reflected = _pauli_su2_reverse_base3_code(code, support_size)
        mat[reflected + 1, code + 1] = 1.0 + 0.0im
    end
    return mat
end

function _pauli_su2_label_rows(labels)
    rows_by_spin_m = Dict{Tuple{Int,Int},Vector{Int}}()
    for (row, label) in pairs(labels)
        push!(get!(rows_by_spin_m, (label.spin2, label.m2)) do
            Int[]
        end, row)
    end
    return rows_by_spin_m
end

"""
    pauli_su2_word_reflection_diagnostics(support_size; max_dimension=4096)

Return diagnostics for the reversal action on fixed-support Pauli words after
the SU(2) word transform.  Reflection should commute with global SU(2), so it
must be block-diagonal by `(spin2, m2)`.  The returned parity multiplicities are
computed from one magnetic sublevel per spin sector.
"""
function pauli_su2_word_reflection_diagnostics(
    support_size::Integer;
    max_dimension::Integer=4096,
)
    word_transform = pauli_su2_word_transform(support_size; max_dimension)
    reflected_cartesian = _pauli_su2_word_reversal_matrix(word_transform.support_size)
    reflected_coupled =
        word_transform.transform * reflected_cartesian * adjoint(word_transform.transform)
    rows_by_spin_m = _pauli_su2_label_rows(word_transform.labels)
    spin_m_keys = sort!(collect(keys(rows_by_spin_m)))

    spin_magnetic_offblock_residual = 0.0
    for (idx, key_i) in pairs(spin_m_keys)
        rows_i = rows_by_spin_m[key_i]
        for key_j in spin_m_keys[(idx + 1):end]
            rows_j = rows_by_spin_m[key_j]
            spin_magnetic_offblock_residual = max(
                spin_magnetic_offblock_residual,
                _max_abs_entry(reflected_coupled[rows_i, rows_j]),
                _max_abs_entry(reflected_coupled[rows_j, rows_i]),
            )
        end
    end

    spin_values = sort!(unique(first.(spin_m_keys)))
    spin_parities = NamedTuple[]
    magnetic_trace_copy_residual = 0.0
    trace_round_residual = 0.0
    for spin2 in spin_values
        m2_values = sort!(Int[m2 for (s2, m2) in spin_m_keys if s2 == spin2])
        reference_trace = nothing
        reference_multiplicity = nothing
        trace_values = Pair{Int,Float64}[]
        for m2 in m2_values
            rows = rows_by_spin_m[(spin2, m2)]
            block_trace = LinearAlgebra.tr(reflected_coupled[rows, rows])
            trace_real = Float64(real(block_trace))
            trace_round_residual = max(
                trace_round_residual,
                abs(trace_real - round(trace_real)),
                abs(Float64(imag(block_trace))),
            )
            reference_trace === nothing || (
                magnetic_trace_copy_residual =
                    max(magnetic_trace_copy_residual, abs(trace_real - reference_trace))
            )
            reference_trace = reference_trace === nothing ? trace_real : reference_trace
            reference_multiplicity =
                reference_multiplicity === nothing ? length(rows) : reference_multiplicity
            push!(trace_values, m2 => trace_real)
        end

        trace_int = round(Int, reference_trace)
        multiplicity = Int(reference_multiplicity)
        plus_multiplicity = (multiplicity + trace_int) ÷ 2
        minus_multiplicity = (multiplicity - trace_int) ÷ 2
        push!(
            spin_parities,
            (
                spin2=spin2,
                multiplicity=multiplicity,
                plus_multiplicity=plus_multiplicity,
                minus_multiplicity=minus_multiplicity,
                trace=trace_int,
                trace_values=trace_values,
            ),
        )
    end

    hermitian_residual = _max_abs_entry(reflected_coupled - adjoint(reflected_coupled))
    involution_residual = _max_abs_entry(reflected_coupled * reflected_coupled - I)
    return (
        support_size=word_transform.support_size,
        dimension=word_transform.dimension,
        spin_parities=spin_parities,
        unitarity_residual=word_transform.unitarity_residual,
        hermitian_residual=hermitian_residual,
        involution_residual=involution_residual,
        spin_magnetic_offblock_residual=spin_magnetic_offblock_residual,
        magnetic_trace_copy_residual=magnetic_trace_copy_residual,
        trace_round_residual=trace_round_residual,
        max_reflection_residual=max(
            hermitian_residual,
            involution_residual,
            spin_magnetic_offblock_residual,
            magnetic_trace_copy_residual,
            trace_round_residual,
        ),
    )
end

"""
    pauli_su2_operator_spin_diagnostics(k; max_dimension=4096)

Return numerical diagnostics showing that the Pauli operator space
`(spin-0 ⊕ spin-1)^{⊗k}` decomposes into the same total-spin sectors reported
by [`pauli_su2_operator_blocks`](@ref).  This dense check is meant for small
local supports.
"""
function pauli_su2_operator_spin_diagnostics(
    k::Integer;
    max_dimension::Integer=4096,
)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    return (
        k=kk,
        _pauli_su2_tensor_spin_diagnostics(
            _PAULI_SU2_OPERATOR_LOCAL_STATES,
            kk,
            Int(max_dimension),
        )...,
    )
end

"""
    pauli_su2_word_blocks(support_size)

Return the SU(2) irrep summary for Pauli words on a fixed active support of
length `support_size`.  This is `(spin-1)^{⊗support_size}` and excludes local
identity choices.
"""
function pauli_su2_word_blocks(support_size::Integer)
    counts = _su2_tensor_power_multiplicities((2,), support_size)
    return [
        PauliSU2WordBlock(spin2, counts[spin2])
        for spin2 in sort!(collect(keys(counts)))
    ]
end

"""
    pauli_su2_word_metrics(support_size; scalar_bytes=sizeof(Float64))

Return storage-size estimates for Pauli words on a fixed active support before
and after the SU(2) irreducible-tensor split.
"""
function pauli_su2_word_metrics(support_size::Integer; scalar_bytes::Integer=sizeof(Float64))
    s = Int(support_size)
    s >= 0 || throw(DomainError(support_size, "`support_size` must be non-negative."))

    blocks = pauli_su2_word_blocks(s)
    return (blocks=blocks, _pauli_su2_storage_metrics(blocks, 3^s, scalar_bytes)...)
end

"""
    pauli_su2_word_reduction_diagnostics(support_size; scalar_bytes=sizeof(Float64))

Return entry-count diagnostics for the SU(2) reduction of fixed-support Pauli
words `(spin-1)^{⊗support_size}`.
"""
function pauli_su2_word_reduction_diagnostics(
    support_size::Integer;
    scalar_bytes::Integer=sizeof(Float64),
)
    metrics = pauli_su2_word_metrics(support_size; scalar_bytes)
    return (; metrics..., _pauli_su2_reduction_accounting(metrics)...)
end

"""
    pauli_su2_word_spin_diagnostics(support_size; max_dimension=4096)

Return numerical diagnostics showing that fixed-support Pauli words decompose as
`(spin-1)^{⊗support_size}`, matching [`pauli_su2_word_blocks`](@ref).  This
dense check is intended for small support sizes.
"""
function pauli_su2_word_spin_diagnostics(
    support_size::Integer;
    max_dimension::Integer=4096,
)
    s = Int(support_size)
    s >= 0 || throw(DomainError(support_size, "`support_size` must be non-negative."))

    return (
        support_size=s,
        _pauli_su2_tensor_spin_diagnostics(
            _PAULI_SU2_WORD_LOCAL_STATES,
            s,
            Int(max_dimension),
        )...,
    )
end

function _pauli_support_tuple(mono::NormalMonomial{PauliAlgebra})
    isempty(mono.word) && return ()

    sites = Int[_pauli_site(idx) for idx in mono.word]
    sort!(sites)
    unique!(sites)
    return Tuple(sites)
end

function _pauli_support_word_counts(basis::AbstractVector{<:NormalMonomial{PauliAlgebra}})
    !isempty(basis) || throw(ArgumentError("Pauli SU(2) basis summary requires a non-empty basis."))

    support_word_counts = Dict{Tuple{Vararg{Int}},Int}()
    for mono in basis
        support = _pauli_support_tuple(mono)
        support_word_counts[support] = get(support_word_counts, support, 0) + 1
    end
    return support_word_counts
end

function _check_pauli_su2_support_complete(support_word_counts::Dict{Tuple{Vararg{Int}},Int})
    for (support, observed) in support_word_counts
        expected = 3^length(support)
        observed == expected || throw(ArgumentError(
            "Pauli SU(2) basis summary requires complete active supports; support $(support) has $observed words, expected $expected."
        ))
    end
    return nothing
end

function _translation_support_axis_orbit_representative(
    mono::NormalMonomial{PauliAlgebra},
    n_sites::Integer,
)
    n = Int(n_sites)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    support = _pauli_support_tuple(mono)
    isempty(support) && return (support=(), axes=())
    all(site -> 1 <= site <= n, support) || throw(ArgumentError(
        "Support $support is not contained in a Pauli chain with n_sites=$n."
    ))

    shifted = Vector{Tuple{Int,UInt8}}(undef, length(mono.word))
    best = nothing
    for shift in 0:(n - 1)
        for (idx, pauli_idx) in pairs(mono.word)
            shifted[idx] = (
                mod1(_pauli_site(pauli_idx) + shift, n),
                UInt8(_pauli_type(pauli_idx)),
            )
        end
        sort!(shifted; by=first)
        candidate = (
            support=Tuple(site_axis[1] for site_axis in shifted),
            axes=Tuple(site_axis[2] for site_axis in shifted),
        )
        if best === nothing || isless(candidate, best)
            best = candidate
        end
    end
    return best
end

function _pauli_translation_support_orbit_word_patterns(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer,
)
    !isempty(basis) || throw(ArgumentError(
        "Pauli translation-orbit SU(2) basis summary requires a non-empty basis."
    ))

    support_word_counts = Dict{Tuple{Vararg{Int}},Int}()
    support_axis_patterns = Dict{Tuple{Vararg{Int}},Set{Tuple{Vararg{UInt8}}}}()
    for mono in basis
        orbit_key = _translation_support_axis_orbit_representative(mono, n_sites)
        support = orbit_key.support
        support_word_counts[support] = get(support_word_counts, support, 0) + 1
        patterns = get!(support_axis_patterns, support) do
            Set{Tuple{Vararg{UInt8}}}()
        end
        push!(patterns, orbit_key.axes)
    end
    return support_word_counts, support_axis_patterns
end

function _check_pauli_su2_translation_orbit_support_complete(
    support_word_counts::Dict{Tuple{Vararg{Int}},Int},
    support_axis_patterns::Dict{Tuple{Vararg{Int}},Set{Tuple{Vararg{UInt8}}}},
)
    for (support, observed) in support_word_counts
        expected = 3^length(support)
        observed == expected || throw(ArgumentError(
            "Pauli translation-orbit SU(2) basis summary requires complete active " *
            "support orbits; support orbit $(support) has $observed words, " *
            "expected $expected."
        ))
        unique_patterns = length(support_axis_patterns[support])
        unique_patterns == expected || throw(ArgumentError(
            "Pauli translation-orbit SU(2) basis summary requires complete active " *
            "support orbits; support orbit $(support) has $unique_patterns unique " *
            "axis patterns, expected $expected."
        ))
    end
    return nothing
end

"""
    pauli_su2_basis_summary(basis)

Summarize the SU(2) content of a Pauli word basis whose active supports are
complete.  For every support of size `s`, the basis must contain all `3^s`
non-identity Pauli words on that support.
"""
function pauli_su2_basis_summary(basis::AbstractVector{<:NormalMonomial{PauliAlgebra}})
    support_word_counts = _pauli_support_word_counts(basis)
    _check_pauli_su2_support_complete(support_word_counts)

    support_size_counts = Dict{Int,Int}()
    spin_multiplicities = Dict{Int,Int}()
    for support in keys(support_word_counts)
        support_size = length(support)
        support_size_counts[support_size] = get(support_size_counts, support_size, 0) + 1
        for block in pauli_su2_word_blocks(support_size)
            spin_multiplicities[block.spin2] =
                get(spin_multiplicities, block.spin2, 0) + block.multiplicity
        end
    end

    support_counts = [
        support_size => support_size_counts[support_size]
        for support_size in sort!(collect(keys(support_size_counts)))
    ]
    blocks = [
        PauliSU2WordBlock(spin2, spin_multiplicities[spin2])
        for spin2 in sort!(collect(keys(spin_multiplicities)))
    ]
    return (support_counts=support_counts, blocks=blocks, word_count=length(basis))
end

"""
    pauli_su2_translation_orbit_basis_summary(basis, n_sites)

Summarize the SU(2) content of a translation-orbit Pauli word basis.  For every
active support orbit of size `s`, the basis must contain one representative for
each of the `3^s` non-identity Pauli words modulo chain translations.  This is
a structural summary for the future translation/SU(2) reducer; it does not
build SU(2)-reduced PSD blocks.
"""
function pauli_su2_translation_orbit_basis_summary(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer,
)
    support_word_counts, support_axis_patterns =
        _pauli_translation_support_orbit_word_patterns(basis, n_sites)
    _check_pauli_su2_translation_orbit_support_complete(
        support_word_counts,
        support_axis_patterns,
    )

    support_size_counts = Dict{Int,Int}()
    spin_multiplicities = Dict{Int,Int}()
    for support in keys(support_word_counts)
        support_size = length(support)
        support_size_counts[support_size] = get(support_size_counts, support_size, 0) + 1
        for block in pauli_su2_word_blocks(support_size)
            spin_multiplicities[block.spin2] =
                get(spin_multiplicities, block.spin2, 0) + block.multiplicity
        end
    end

    support_counts = [
        support_size => support_size_counts[support_size]
        for support_size in sort!(collect(keys(support_size_counts)))
    ]
    blocks = [
        PauliSU2WordBlock(spin2, spin_multiplicities[spin2])
        for spin2 in sort!(collect(keys(spin_multiplicities)))
    ]
    return (support_counts=support_counts, blocks=blocks, word_count=length(basis))
end

function _pauli_su2_translation_orbit_basis_multiplicity_label(
    support_orbit::Tuple,
    spin2::Int,
    local_multiplicity::Int,
)
    return (
        feature=:pauli_su2_translation_orbit_multiplicity,
        support_orbit=support_orbit,
        support_size=length(support_orbit),
        spin2=spin2,
        local_multiplicity=local_multiplicity,
    )
end

function _pauli_su2_basis_multiplicity_label(
    support::Tuple,
    spin2::Int,
    local_multiplicity::Int,
)
    return (
        feature=:pauli_su2_basis_multiplicity,
        support=support,
        support_size=length(support),
        spin2=spin2,
        local_multiplicity=local_multiplicity,
    )
end

"""
    pauli_su2_translation_orbit_basis_reduction_plan(basis, n_sites; scalar_bytes=sizeof(Float64))

Return the explicit SU(2) multiplicity-block plan for a translation-orbit
Pauli word basis.  The plan is structural only: it labels the multiplicity
rows that a future translation/SU(2) reducer must build, but does not construct
PSD matrices or Clebsch-Gordan transforms.
"""
function pauli_su2_translation_orbit_basis_reduction_plan(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer;
    scalar_bytes::Integer=sizeof(Float64),
)
    summary = pauli_su2_translation_orbit_basis_summary(basis, n_sites)
    support_word_counts, support_axis_patterns =
        _pauli_translation_support_orbit_word_patterns(basis, n_sites)
    _check_pauli_su2_translation_orbit_support_complete(
        support_word_counts,
        support_axis_patterns,
    )

    labels_by_spin = Dict{Int,Vector{Any}}()
    supports = sort!(collect(keys(support_word_counts)); by=support -> (length(support), support))
    for support in supports
        support_size = length(support)
        for word_block in pauli_su2_word_blocks(support_size)
            labels = get!(labels_by_spin, word_block.spin2) do
                Any[]
            end
            for local_multiplicity in 1:word_block.multiplicity
                push!(
                    labels,
                    _pauli_su2_translation_orbit_basis_multiplicity_label(
                        support,
                        word_block.spin2,
                        local_multiplicity,
                    ),
                )
            end
        end
    end

    blocks = [
        PauliSU2BasisBlock(spin2, labels_by_spin[spin2])
        for spin2 in sort!(collect(keys(labels_by_spin)))
    ]
    return (
        support_counts=summary.support_counts,
        blocks=blocks,
        _pauli_su2_storage_metrics(blocks, summary.word_count, scalar_bytes)...,
    )
end

"""
    pauli_su2_translation_orbit_basis_reduction_diagnostics(basis, n_sites; scalar_bytes=sizeof(Float64))

Return entry-count diagnostics for the structural SU(2) reduction plan of a
translation-orbit Pauli word basis.
"""
function pauli_su2_translation_orbit_basis_reduction_diagnostics(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer;
    scalar_bytes::Integer=sizeof(Float64),
)
    metrics = pauli_su2_translation_orbit_basis_reduction_plan(
        basis,
        n_sites;
        scalar_bytes,
    )
    return (; metrics..., _pauli_su2_reduction_accounting(metrics)...)
end

"""
    pauli_su2_basis_reduction_plan(basis; scalar_bytes=sizeof(Float64))

Return the explicit SU(2) multiplicity-block plan for a support-complete Pauli
word basis.  The plan contains one [`PauliSU2BasisBlock`](@ref) per total spin
sector, with stable logical row labels for the surviving multiplicity rows.
"""
function pauli_su2_basis_reduction_plan(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    scalar_bytes::Integer=sizeof(Float64),
)
    summary = pauli_su2_basis_summary(basis)
    support_word_counts = _pauli_support_word_counts(basis)
    _check_pauli_su2_support_complete(support_word_counts)

    labels_by_spin = Dict{Int,Vector{Any}}()
    supports = sort!(collect(keys(support_word_counts)); by=support -> (length(support), support))
    for support in supports
        support_size = length(support)
        for word_block in pauli_su2_word_blocks(support_size)
            labels = get!(labels_by_spin, word_block.spin2) do
                Any[]
            end
            for local_multiplicity in 1:word_block.multiplicity
                push!(
                    labels,
                    _pauli_su2_basis_multiplicity_label(
                        support,
                        word_block.spin2,
                        local_multiplicity,
                    ),
                )
            end
        end
    end

    blocks = [
        PauliSU2BasisBlock(spin2, labels_by_spin[spin2])
        for spin2 in sort!(collect(keys(labels_by_spin)))
    ]
    return (
        support_counts=summary.support_counts,
        blocks=blocks,
        _pauli_su2_storage_metrics(blocks, summary.word_count, scalar_bytes)...,
    )
end

function _pauli_su2_basis_state_label(
    support::Tuple,
    spin2::Int,
    m2::Int,
    local_multiplicity::Int,
)
    return (
        feature=:pauli_su2_basis_state,
        support=support,
        support_size=length(support),
        spin2=spin2,
        m2=m2,
        local_multiplicity=local_multiplicity,
    )
end

function _pauli_su2_translation_orbit_basis_state_label(
    support_orbit::Tuple,
    spin2::Int,
    m2::Int,
    local_multiplicity::Int,
)
    return (
        feature=:pauli_su2_translation_orbit_basis_state,
        support_orbit=support_orbit,
        support_size=length(support_orbit),
        spin2=spin2,
        m2=m2,
        local_multiplicity=local_multiplicity,
    )
end

function _pauli_su2_basis_lookup(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
) where {T<:Integer}
    lookup = Dict{NormalMonomial{PauliAlgebra,T},Int}()
    for idx in eachindex(basis)
        mono = basis[idx]
        haskey(lookup, mono) && throw(ArgumentError(
            "Pauli SU(2) basis transform requires unique basis words; duplicate word $mono."
        ))
        lookup[mono] = Int(idx)
    end
    return lookup
end

function _pauli_su2_support_columns(
    ::Type{T},
    support::Tuple,
    basis_lookup::Dict{NormalMonomial{PauliAlgebra,T},Int},
) where {T<:Integer}
    sites = Int[site for site in support]
    columns = Vector{Int}(undef, 3^length(sites))
    for code in 0:(length(columns) - 1)
        mono = _pauli_chain_word(T, sites, code)
        column = get(basis_lookup, mono, nothing)
        column === nothing && throw(ArgumentError(
            "Pauli SU(2) basis transform requires complete active supports; missing word $mono on support $support."
        ))
        columns[code + 1] = column
    end
    return columns
end

function _pauli_su2_translation_orbit_basis_lookup(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer,
)
    lookup = Dict{Any,Int}()
    for idx in eachindex(basis)
        orbit_key = _translation_support_axis_orbit_representative(basis[idx], n_sites)
        key = (support=orbit_key.support, axes=orbit_key.axes)
        haskey(lookup, key) && throw(ArgumentError(
            "Pauli translation-orbit SU(2) transform requires unique support/axis orbit representatives; duplicate key $key."
        ))
        lookup[key] = Int(idx)
    end
    return lookup
end

function _pauli_type_word_from_code(width::Int, code0::Int)
    word = Vector{UInt8}(undef, width)
    code = code0
    @inbounds for idx in 1:width
        word[idx] = UInt8(mod(code, 3))
        code = div(code, 3)
    end
    return Tuple(word)
end

function _pauli_su2_translation_orbit_support_columns(
    support_orbit::Tuple,
    basis_lookup::Dict{Any,Int},
)
    columns = Vector{Int}(undef, 3^length(support_orbit))
    for code in 0:(length(columns) - 1)
        axes = _pauli_type_word_from_code(length(support_orbit), code)
        key = (support=support_orbit, axes=axes)
        column = get(basis_lookup, key, nothing)
        column === nothing && throw(ArgumentError(
            "Pauli translation-orbit SU(2) transform requires complete active support orbits; missing axis pattern $axes on support orbit $support_orbit."
        ))
        columns[code + 1] = column
    end
    return columns
end

function _pauli_monomial_from_support_axes(
    ::Type{T},
    support::Tuple,
    axes::Tuple,
) where {T<:Integer}
    length(support) == length(axes) || throw(DimensionMismatch(
        "Support/axis pattern length mismatch: $(length(support)) sites, $(length(axes)) axes."
    ))
    word = T[
        _pauli_index_from_site_type(T, Int(site), UInt8(axis))
        for (site, axis) in zip(support, axes)
    ]
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError(
        "Internal error: support-axis Pauli monomial construction produced phase $phase."
    ))
    return NormalMonomial{PauliAlgebra,T}(copy(simplified))
end

function _pauli_translation_shift_between(
    source::NormalMonomial{PauliAlgebra,T},
    target::NormalMonomial{PauliAlgebra,T},
    n_sites::Integer,
) where {T<:Unsigned}
    n = Int(n_sites)
    for shift in 0:(n - 1)
        _translate_pauli_monomial(source, shift, n) == target && return shift
    end
    throw(ArgumentError(
        "Internal error: $target is not in the translation orbit of $source."
    ))
end

function _pauli_su2_translation_orbit_momentum_rephase_factors(
    basis::AbstractVector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    momentum::Integer,
) where {T<:Unsigned}
    n = Int(n_sites)
    k = mod(Int(momentum), n)
    phases = Vector{ComplexF64}(undef, length(basis))
    for (idx, mono) in pairs(basis)
        if isone(mono)
            phases[idx] = 1.0 + 0.0im
            continue
        end
        orbit_key = _translation_support_axis_orbit_representative(mono, n)
        canonical = _pauli_monomial_from_support_axes(
            T,
            orbit_key.support,
            orbit_key.axes,
        )
        shift = _pauli_translation_shift_between(canonical, mono, n)
        phases[idx] = conj(ComplexF64(_momentum_phase(Float64, k, shift, n)))
    end
    return phases
end

function _pauli_su2_rephase_transform_blocks(
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock},
    phases::AbstractVector{ComplexF64},
)
    phase_matrix = Diagonal(phases)
    return PauliSU2BasisTransformBlock[
        PauliSU2BasisTransformBlock(
            block.spin2,
            block.row_labels,
            block.transform * phase_matrix;
            coefficient_domain=block.coefficient_domain,
            exact_coefficient_domain=block.exact_coefficient_domain,
        )
        for block in transform_blocks
    ]
end

function _push_pauli_su2_basis_transform_row!(
    rows_by_spin::Dict{Int,Vector{Vector{ComplexF64}}},
    labels_by_spin::Dict{Int,Vector{Any}},
    spin2::Int,
    label,
    row::Vector{ComplexF64},
)
    rows = get!(rows_by_spin, spin2) do
        Vector{ComplexF64}[]
    end
    labels = get!(labels_by_spin, spin2) do
        Any[]
    end
    push!(rows, row)
    push!(labels, label)
    return nothing
end

function _pauli_su2_basis_transform_matrix(rows::Vector{Vector{ComplexF64}}, width::Int)
    mat = Matrix{ComplexF64}(undef, length(rows), width)
    for (idx, row) in pairs(rows)
        length(row) == width || throw(DimensionMismatch(
            "Pauli SU(2) basis transform row has width $(length(row)), expected $width."
        ))
        mat[idx, :] = row
    end
    return mat
end

function _pauli_su2_basis_transform_blocks_from_columns(
    supports,
    basis_length::Integer,
    columns_for_support,
    label_for_row;
    max_dimension::Integer=4096,
)
    rows_by_spin = Dict{Int,Vector{Vector{ComplexF64}}}()
    labels_by_spin = Dict{Int,Vector{Any}}()
    word_transforms_by_support_size = Dict{Int,Any}()
    for support in supports
        support_size = length(support)
        columns = columns_for_support(support)
        word_transform = get!(word_transforms_by_support_size, support_size) do
            pauli_su2_word_transform(support_size; max_dimension)
        end
        size(word_transform.transform, 2) == length(columns) || throw(DimensionMismatch(
            "Pauli SU(2) word transform has $(size(word_transform.transform, 2)) columns for support $support, expected $(length(columns))."
        ))

        local_multiplicity_counts = Dict{Tuple{Int,Int},Int}()
        for row_idx in axes(word_transform.transform, 1)
            word_label = word_transform.labels[row_idx]
            spin2 = Int(word_label.spin2)
            m2 = Int(word_label.m2)
            local_key = (spin2, m2)
            local_multiplicity = get(local_multiplicity_counts, local_key, 0) + 1
            local_multiplicity_counts[local_key] = local_multiplicity

            row = zeros(ComplexF64, Int(basis_length))
            row[columns] = word_transform.transform[row_idx, :]
            _push_pauli_su2_basis_transform_row!(
                rows_by_spin,
                labels_by_spin,
                spin2,
                label_for_row(support, spin2, m2, local_multiplicity),
                row,
            )
        end
    end

    return [
        PauliSU2BasisTransformBlock(
            spin2,
            labels_by_spin[spin2],
            _pauli_su2_basis_transform_matrix(rows_by_spin[spin2], Int(basis_length)),
        )
        for spin2 in sort!(collect(keys(rows_by_spin)))
    ]
end

function _pauli_su2_singlet_channel_data(
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock},
    basis_size::Integer,
)
    singlet_block = findfirst(block -> block.spin2 == 0, transform_blocks)
    if singlet_block === nothing
        transform = zeros(ComplexF64, 0, Int(basis_size))
        row_labels = Any[]
        coefficient_domain = :complex_algebraic_float64
        exact_coefficient_domain = :complex_sqrt_rational
    else
        block = transform_blocks[singlet_block]
        transform = block.transform
        row_labels = block.row_labels
        coefficient_domain = block.coefficient_domain
        exact_coefficient_domain = block.exact_coefficient_domain
    end
    residual = isempty(row_labels) ? 0.0 :
        maximum(abs, transform * adjoint(transform) - I)

    return (
        channel_count=length(row_labels),
        row_labels=row_labels,
        coefficient_domain=coefficient_domain,
        exact_coefficient_domain=exact_coefficient_domain,
        singlet_orthonormality_residual=residual,
        transform=transform,
    )
end

function _pauli_su2_singlet_channel_equality_data(
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock},
    basis_size::Integer;
    atol::Real=1e-12,
)
    selected_blocks = [
        block
        for block in transform_blocks
        if block.spin2 != 0
    ]
    singlet_data = _pauli_su2_singlet_channel_data(transform_blocks, basis_size)

    if isempty(selected_blocks)
        transform = zeros(ComplexF64, 0, Int(basis_size))
        row_labels = Any[]
        coefficient_domain = :complex_algebraic_float64
        exact_coefficient_domain = :complex_sqrt_rational
    else
        transform = reduce(vcat, [block.transform for block in selected_blocks])
        row_labels = Any[]
        for block in selected_blocks
            append!(row_labels, block.row_labels)
        end
        coefficient_domain = first(selected_blocks).coefficient_domain
        exact_coefficient_domain = first(selected_blocks).exact_coefficient_domain
    end

    return (
        equality_count=length(row_labels),
        row_labels=row_labels,
        coefficient_domain=coefficient_domain,
        exact_coefficient_domain=exact_coefficient_domain,
        equality_orthonormality_residual=_pauli_su2_row_orthonormality_residual(
            transform,
        ),
        singlet_orthogonality_residual=_pauli_su2_row_cross_residual(
            transform,
            singlet_data.transform,
        ),
        transform=transform,
        column_forms=_pauli_su2_transform_column_forms(transform; atol),
    )
end

function _pauli_su2_column_forms_to_basis_forms(
    column_forms,
    basis::AbstractVector,
)
    return [
        [
            (basis[column], coefficient)
            for (column, coefficient) in form
        ]
        for form in column_forms
    ]
end

"""
    pauli_su2_basis_transform_blocks(basis; max_dimension=4096)

Build concrete SU(2) row-transform blocks for a support-complete Pauli word
basis.  Rows are grouped by total spin and columns follow the input `basis`
order, so stacking all returned block transforms gives a unitary transform of
the full basis.
"""
function pauli_su2_basis_transform_blocks(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}};
    max_dimension::Integer=4096,
) where {T<:Integer}
    support_word_counts = _pauli_support_word_counts(basis)
    _check_pauli_su2_support_complete(support_word_counts)
    basis_lookup = _pauli_su2_basis_lookup(basis)
    supports = sort!(collect(keys(support_word_counts)); by=support -> (length(support), support))

    return _pauli_su2_basis_transform_blocks_from_columns(
        supports,
        length(basis),
        support -> _pauli_su2_support_columns(T, support, basis_lookup),
        _pauli_su2_basis_state_label;
        max_dimension,
    )
end

"""
    pauli_su2_basis_singlet_channels(basis; max_dimension=4096)

Return the spin-0 rows of the support-complete Pauli-word SU(2) transform.
The input basis must satisfy the same support-completeness contract as
[`pauli_su2_basis_transform_blocks`](@ref).
"""
function pauli_su2_basis_singlet_channels(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    max_dimension::Integer=4096,
)
    transform_blocks = pauli_su2_basis_transform_blocks(basis; max_dimension)
    data = _pauli_su2_singlet_channel_data(transform_blocks, length(basis))
    return (
        basis_size=length(basis),
        data...,
        channel_supports=[
            label.support
            for label in data.row_labels
            if label isa NamedTuple && haskey(label, :support)
        ],
    )
end

"""
    pauli_su2_basis_singlet_channel_equalities(basis; max_dimension=4096, atol=1e-12)

Return the non-singlet SU(2) transform rows for a support-complete Pauli word
basis.  The returned `column_forms` are sparse forms over the input basis
columns and are the candidate scalar equalities for singlet-channel emission.
"""
function pauli_su2_basis_singlet_channel_equalities(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    max_dimension::Integer=4096,
    atol::Real=1e-12,
)
    transform_blocks = pauli_su2_basis_transform_blocks(basis; max_dimension)
    data = _pauli_su2_singlet_channel_equality_data(
        transform_blocks,
        length(basis);
        atol,
    )
    return (
        basis_size=length(basis),
        data...,
        basis_forms=_pauli_su2_column_forms_to_basis_forms(data.column_forms, basis),
        equality_supports=[
            label.support
            for label in data.row_labels
            if label isa NamedTuple && haskey(label, :support)
        ],
    )
end

"""
    pauli_su2_translation_orbit_basis_transform_blocks(basis, n_sites; max_dimension=4096)

Build concrete SU(2) row-transform blocks for a translation-orbit Pauli word
basis.  The input must contain one representative for every non-identity Pauli
axis pattern on each active support orbit.  Rows are grouped by total spin and
columns follow the input `basis` order.
"""
function pauli_su2_translation_orbit_basis_transform_blocks(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    max_dimension::Integer=4096,
)
    throw(ArgumentError(
        "`pauli_su2_translation_orbit_basis_transform_blocks` requires `n_sites`."
    ))
end

function pauli_su2_translation_orbit_basis_transform_blocks(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    max_dimension::Integer=4096,
    momentum::Union{Nothing,Integer}=nothing,
) where {T<:Integer}
    support_word_counts, support_axis_patterns =
        _pauli_translation_support_orbit_word_patterns(basis, n_sites)
    _check_pauli_su2_translation_orbit_support_complete(
        support_word_counts,
        support_axis_patterns,
    )
    basis_lookup = _pauli_su2_translation_orbit_basis_lookup(basis, n_sites)
    supports = sort!(collect(keys(support_word_counts)); by=support -> (length(support), support))

    blocks = _pauli_su2_basis_transform_blocks_from_columns(
        supports,
        length(basis),
        support -> _pauli_su2_translation_orbit_support_columns(support, basis_lookup),
        _pauli_su2_translation_orbit_basis_state_label;
        max_dimension,
    )
    momentum === nothing && return blocks
    T <: Unsigned || throw(ArgumentError(
        "Momentum-rephased translation/SU(2) transforms require unsigned Pauli indices."
    ))
    phases = _pauli_su2_translation_orbit_momentum_rephase_factors(
        NormalMonomial{PauliAlgebra,T}[mono for mono in basis],
        n_sites,
        momentum,
    )
    return _pauli_su2_rephase_transform_blocks(blocks, phases)
end

"""
    pauli_su2_translation_orbit_singlet_channels(basis, n_sites; max_dimension=4096, momentum=nothing)

Return the spin-0 rows of the translation-orbit Pauli-word SU(2) transform.
With `momentum` set, the same singlet rows are returned after the translation
momentum rephasing used by the translation/SU(2) reducer.
"""
function pauli_su2_translation_orbit_singlet_channels(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    max_dimension::Integer=4096,
    momentum::Union{Nothing,Integer}=nothing,
)
    throw(ArgumentError(
        "`pauli_su2_translation_orbit_singlet_channels` requires `n_sites`."
    ))
end

function pauli_su2_translation_orbit_singlet_channels(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer;
    max_dimension::Integer=4096,
    momentum::Union{Nothing,Integer}=nothing,
)
    transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
        basis,
        n_sites;
        max_dimension,
        momentum,
    )
    data = _pauli_su2_singlet_channel_data(transform_blocks, length(basis))

    return (
        n_sites=Int(n_sites),
        basis_size=length(basis),
        momentum=isnothing(momentum) ? nothing : mod(Int(momentum), Int(n_sites)),
        data...,
        channel_support_orbits=[
            label.support_orbit
            for label in data.row_labels
            if label isa NamedTuple && haskey(label, :support_orbit)
        ],
    )
end

"""
    pauli_su2_translation_orbit_singlet_channel_equalities(basis, n_sites; max_dimension=4096, momentum=nothing, atol=1e-12)

Return the non-singlet SU(2) rows for a translation-orbit Pauli word basis,
optionally after the same momentum rephasing used by the translation/SU(2)
reducer.
"""
function pauli_su2_translation_orbit_singlet_channel_equalities(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    max_dimension::Integer=4096,
    momentum::Union{Nothing,Integer}=nothing,
    atol::Real=1e-12,
)
    throw(ArgumentError(
        "`pauli_su2_translation_orbit_singlet_channel_equalities` requires `n_sites`."
    ))
end

function pauli_su2_translation_orbit_singlet_channel_equalities(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer;
    max_dimension::Integer=4096,
    momentum::Union{Nothing,Integer}=nothing,
    atol::Real=1e-12,
)
    transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
        basis,
        n_sites;
        max_dimension,
        momentum,
    )
    data = _pauli_su2_singlet_channel_equality_data(
        transform_blocks,
        length(basis);
        atol,
    )

    return (
        n_sites=Int(n_sites),
        basis_size=length(basis),
        momentum=isnothing(momentum) ? nothing : mod(Int(momentum), Int(n_sites)),
        data...,
        basis_forms=_pauli_su2_column_forms_to_basis_forms(data.column_forms, basis),
        equality_support_orbits=[
            label.support_orbit
            for label in data.row_labels
            if label isa NamedTuple && haskey(label, :support_orbit)
        ],
    )
end

function _pauli_su2_reference_m2(block::PauliSU2BasisTransformBlock, reference_m2)
    ref = if reference_m2 === :lowest
        -block.spin2
    elseif reference_m2 === :highest
        block.spin2
    elseif reference_m2 === :zero
        0
    elseif reference_m2 isa Integer
        Int(reference_m2)
    else
        throw(ArgumentError(
            "Unsupported SU(2) reference_m2 $(repr(reference_m2)); expected :lowest, :highest, :zero, or an integer."
        ))
    end

    (-block.spin2 <= ref <= block.spin2 && iseven(ref + block.spin2)) || throw(ArgumentError(
        "Invalid reference_m2=$ref for spin2=$(block.spin2)."
    ))
    return ref
end

function _pauli_su2_basis_reference_rows(
    block::PauliSU2BasisTransformBlock,
    reference_m2,
)
    ref = _pauli_su2_reference_m2(block, reference_m2)
    rows = _pauli_su2_basis_m2_rows(block, ref)
    return ref, rows, Any[block.row_labels[row] for row in rows]
end

function _pauli_su2_basis_reference_rows(
    block::PauliSU2BasisTransformBlock,
    reference_m2,
    rows_by_m2::Dict{Int,Vector{Int}},
)
    ref = _pauli_su2_reference_m2(block, reference_m2)
    rows = _pauli_su2_basis_m2_rows(block, ref, rows_by_m2)
    return ref, rows, Any[block.row_labels[row] for row in rows]
end

function _pauli_su2_basis_m2_rows(block::PauliSU2BasisTransformBlock, m2::Integer)
    mm2 = Int(m2)
    rows = findall(label -> label.m2 == mm2, block.row_labels)
    length(rows) == block.multiplicity || throw(ArgumentError(
        "Pauli SU(2) transform block for spin2=$(block.spin2) has $(length(rows)) rows for m2=$mm2, expected $(block.multiplicity)."
    ))
    return rows
end

function _pauli_su2_basis_m2_rows(
    block::PauliSU2BasisTransformBlock,
    m2::Integer,
    rows_by_m2::Dict{Int,Vector{Int}},
)
    mm2 = Int(m2)
    rows = get(rows_by_m2, mm2, Int[])
    length(rows) == block.multiplicity || throw(ArgumentError(
        "Pauli SU(2) transform block for spin2=$(block.spin2) has $(length(rows)) rows for m2=$mm2, expected $(block.multiplicity)."
    ))
    return rows
end

function _pauli_su2_basis_m2_row_lookup(block::PauliSU2BasisTransformBlock)
    rows_by_m2 = Dict{Int,Vector{Int}}()
    for (row, label) in pairs(block.row_labels)
        push!(get!(rows_by_m2, Int(label.m2)) do
            Int[]
        end, row)
    end
    for m2 in (-block.spin2):2:block.spin2
        _pauli_su2_basis_m2_rows(block, m2, rows_by_m2)
    end
    return rows_by_m2
end

function _check_pauli_su2_basis_moment_matrix(mat::AbstractMatrix, basis_length::Integer)
    size(mat, 1) == size(mat, 2) || throw(DimensionMismatch(
        "Pauli SU(2) basis moment block must be square, got $(size(mat))."
    ))
    size(mat, 1) == Int(basis_length) || throw(DimensionMismatch(
        "Pauli SU(2) basis moment block has side $(size(mat, 1)), expected $(Int(basis_length))."
    ))
    return nothing
end

function _pauli_su2_reduced_moment_blocks_from_transform_blocks(
    mat::AbstractMatrix,
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock};
    reference_m2,
)
    return [
        let ref_rows = _pauli_su2_basis_reference_rows(block, reference_m2)
            ref_m2, rows, row_labels = ref_rows
            transform = block.transform[rows, :]
            (
                spin2=block.spin2,
                multiplicity=block.multiplicity,
                irrep_dimension=block.irrep_dimension,
                reference_m2=ref_m2,
                row_labels=row_labels,
                matrix=transform * mat * adjoint(transform),
                transform=block,
                coefficient_domain=block.coefficient_domain,
                exact_coefficient_domain=block.exact_coefficient_domain,
            )
        end
        for block in transform_blocks
    ]
end

"""
    pauli_su2_basis_reduced_moment_blocks(mat, basis; max_dimension=4096, reference_m2=:lowest)

Conjugate a basis-indexed moment block by the SU(2) basis transform and return
one reduced multiplicity block per total-spin sector.  The default keeps the
lowest magnetic sublevel (`m2 = -spin2`) as the representative copy.
"""
function pauli_su2_basis_reduced_moment_blocks(
    mat::AbstractMatrix,
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}};
    max_dimension::Integer=4096,
    reference_m2=:lowest,
) where {T<:Integer}
    _check_pauli_su2_basis_moment_matrix(mat, length(basis))

    transform_blocks = pauli_su2_basis_transform_blocks(basis; max_dimension)
    return _pauli_su2_reduced_moment_blocks_from_transform_blocks(
        mat,
        transform_blocks;
        reference_m2,
    )
end

"""
    pauli_su2_translation_orbit_basis_reduced_moment_blocks(mat, basis, n_sites; max_dimension=4096, reference_m2=:lowest)

Conjugate a translation-orbit-basis moment block by the SU(2) transform on
orbit representatives and return one reduced multiplicity block per total-spin
sector.
"""
function pauli_su2_translation_orbit_basis_reduced_moment_blocks(
    mat::AbstractMatrix,
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer;
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    momentum::Union{Nothing,Integer}=nothing,
)
    _check_pauli_su2_basis_moment_matrix(mat, length(basis))

    transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
        basis,
        n_sites;
        max_dimension,
        momentum,
    )
    return _pauli_su2_reduced_moment_blocks_from_transform_blocks(
        mat,
        transform_blocks;
        reference_m2,
    )
end

_max_abs_entry(mat::AbstractArray{<:Number}) = isempty(mat) ? 0.0 : Float64(maximum(abs, mat))

function _pauli_su2_moment_reduction_diagnostics_from_transform_blocks(
    mat::AbstractMatrix{<:Number},
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock};
    reference_m2,
)
    stacked_transform = reduce(vcat, [block.transform for block in transform_blocks])
    unitarity_residual = _max_abs_entry(stacked_transform * adjoint(stacked_transform) - I)

    spin_offblock_residual = 0.0
    for block_i in eachindex(transform_blocks)
        transform_i = transform_blocks[block_i].transform
        for block_j in (block_i + 1):lastindex(transform_blocks)
            transform_j = transform_blocks[block_j].transform
            cross = transform_i * mat * adjoint(transform_j)
            spin_offblock_residual = max(spin_offblock_residual, _max_abs_entry(cross))
        end
    end

    magnetic_offdiag_residual = 0.0
    magnetic_copy_residual = 0.0
    for block in transform_blocks
        transformed = block.transform * mat * adjoint(block.transform)
        ref_m2, ref_rows, _ = _pauli_su2_basis_reference_rows(block, reference_m2)
        reference_matrix = transformed[ref_rows, ref_rows]
        m2_values = collect((-block.spin2):2:block.spin2)

        for m2 in m2_values
            rows = m2 == ref_m2 ? ref_rows : _pauli_su2_basis_m2_rows(block, m2)
            magnetic_copy_residual = max(
                magnetic_copy_residual,
                _max_abs_entry(transformed[rows, rows] - reference_matrix),
            )
        end
        for (idx, m2_a) in pairs(m2_values)
            rows_a = m2_a == ref_m2 ? ref_rows : _pauli_su2_basis_m2_rows(block, m2_a)
            for m2_b in m2_values[(idx + 1):end]
                rows_b = m2_b == ref_m2 ? ref_rows : _pauli_su2_basis_m2_rows(block, m2_b)
                magnetic_offdiag_residual = max(
                    magnetic_offdiag_residual,
                    _max_abs_entry(transformed[rows_a, rows_b]),
                    _max_abs_entry(transformed[rows_b, rows_a]),
                )
            end
        end
    end

    return (
        full_side=size(mat, 1),
        reduced_block_sizes=[block.multiplicity for block in transform_blocks],
        transformed_block_sizes=[size(block.transform, 1) for block in transform_blocks],
        reference_m2s=[_pauli_su2_reference_m2(block, reference_m2) for block in transform_blocks],
        unitarity_residual=unitarity_residual,
        spin_offblock_residual=spin_offblock_residual,
        magnetic_offdiag_residual=magnetic_offdiag_residual,
        magnetic_copy_residual=magnetic_copy_residual,
        max_reduction_residual=max(
            spin_offblock_residual,
            magnetic_offdiag_residual,
            magnetic_copy_residual,
        ),
    )
end

"""
    pauli_su2_basis_moment_reduction_diagnostics(mat, basis; max_dimension=4096, reference_m2=:lowest)

Return numeric residuals checking whether a basis-indexed moment block has the
SU(2) form `⊕_J (I_J ⊗ X_J)` after the Pauli-word Clebsch-Gordan transform.
"""
function pauli_su2_basis_moment_reduction_diagnostics(
    mat::AbstractMatrix{<:Number},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}};
    max_dimension::Integer=4096,
    reference_m2=:lowest,
) where {T<:Integer}
    _check_pauli_su2_basis_moment_matrix(mat, length(basis))

    transform_blocks = pauli_su2_basis_transform_blocks(basis; max_dimension)
    return _pauli_su2_moment_reduction_diagnostics_from_transform_blocks(
        mat,
        transform_blocks;
        reference_m2,
    )
end

"""
    pauli_su2_translation_orbit_basis_moment_reduction_diagnostics(mat, basis, n_sites; max_dimension=4096, reference_m2=:lowest)

Return numeric residuals checking the SU(2) block form of a moment block indexed
by translation-orbit Pauli representatives.
"""
function pauli_su2_translation_orbit_basis_moment_reduction_diagnostics(
    mat::AbstractMatrix{<:Number},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}},
    n_sites::Integer;
    max_dimension::Integer=4096,
    reference_m2=:lowest,
)
    _check_pauli_su2_basis_moment_matrix(mat, length(basis))

    transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
        basis,
        n_sites;
        max_dimension,
    )
    return _pauli_su2_moment_reduction_diagnostics_from_transform_blocks(
        mat,
        transform_blocks;
        reference_m2,
    )
end

function _collect_pauli_su2_moment_basis(
    objective::Polynomial{PauliAlgebra,T,C},
    constraints::Vector{Tuple{Symbol, Matrix{Polynomial{PauliAlgebra,T,C}}}},
) where {T<:Integer,C<:Number}
    basis = NormalMonomial{PauliAlgebra,T}[]
    append!(basis, monomials(objective))
    for (_, mat) in constraints, poly in mat
        append!(basis, monomials(poly))
    end
    return sorted_unique!(basis)
end

function _pauli_su2_basis_form_polynomial(
    form,
    ::Type{MP},
) where {T<:Integer,C<:Number,MP<:Polynomial{PauliAlgebra,T,C}}
    M = NormalMonomial{PauliAlgebra,T}
    terms = Tuple{C,M}[]
    sizehint!(terms, length(form))
    for (mono, coefficient) in form
        converted = convert(C, coefficient)
        iszero(converted) && continue
        push!(terms, (converted, convert(M, mono)))
    end
    return _polynomial_from_owned_terms!(terms)
end

function _pauli_su2_basis_form_linear_moment_form(
    form,
    ::Type{K},
    ::Type{C},
    ::Type{M};
    atol::Real,
) where {K,C,M<:NormalMonomial{PauliAlgebra}}
    FC = C <: Real ? Complex{C} : C
    tolerance = real(one(FC)) isa Real ? typeof(real(one(FC)))(atol) : atol
    pairs = Pair{K,FC}[]
    sizehint!(pairs, length(form))
    for (mono, coefficient) in form
        converted = convert(FC, coefficient)
        abs(converted) <= tolerance && continue
        push!(pairs, _axis_moment_key(K, convert(M, mono)) => converted)
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _append_pauli_su2_translation_singlet_channel_equality_linear_constraints!(
    builder::MomentLinearBuilder{K,C,M},
    basis::Vector{M},
    n_sites::Integer,
    momentum::Integer;
    max_dimension::Integer,
    phase_atol::Real,
    atol::Real,
    registered_key_tokens=nothing,
) where {K,C,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    equality_data = pauli_su2_translation_orbit_singlet_channel_equalities(
        basis,
        n_sites;
        max_dimension,
        momentum,
        atol,
    )
    for (row_idx, form) in pairs(equality_data.basis_forms)
        linear_form = _pauli_su2_basis_form_linear_moment_form(
            form,
            K,
            C,
            M;
            atol,
        )
        isempty(linear_form) && continue
        label = (
            feature=:pauli_su2_translation_orbit_singlet_channel_equality,
            decomposition=:translation_su2,
            momentum=equality_data.momentum,
            row=Int(row_idx),
            su2_label=equality_data.row_labels[row_idx],
            coefficient_domain=equality_data.coefficient_domain,
            exact_coefficient_domain=equality_data.exact_coefficient_domain,
        )
        _add_translation_zero_linear_form!(
            builder,
            linear_form;
            phase_atol,
            label,
            registered_key_tokens,
            origin_factory=_pauli_su2_singlet_channel_linear_origin,
        )
    end
    return equality_data
end

function _append_pauli_su2_singlet_channel_equality_constraints!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    zero_origin_by_constraint::Dict{Int,Any};
    max_dimension::Integer,
    atol::Real,
) where {T<:Integer,C<:Number,MP<:Polynomial{PauliAlgebra,T,C}}
    equality_data = pauli_su2_basis_singlet_channel_equalities(
        basis;
        max_dimension,
        atol,
    )
    for (row_idx, form) in pairs(equality_data.basis_forms)
        poly = _pauli_su2_basis_form_polynomial(form, MP)
        iszero(poly) && continue
        label = (
            feature=:pauli_su2_singlet_channel_equality,
            row=Int(row_idx),
            su2_label=equality_data.row_labels[row_idx],
            coefficient_domain=equality_data.coefficient_domain,
            exact_coefficient_domain=equality_data.exact_coefficient_domain,
        )
        before = length(constraints)
        _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        for constraint_idx in (before + 1):length(constraints)
            zero_origin_by_constraint[constraint_idx] =
                PauliSU2SingletChannelEqualityOriginSeed(label)
        end
    end
    return equality_data
end

function _append_pauli_su2_translation_orbit_singlet_channel_equality_constraints!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    ::Type{MP},
    zero_origin_by_constraint::Dict{Int,Any};
    max_dimension::Integer,
    atol::Real,
) where {T<:Integer,C<:Number,MP<:Polynomial{PauliAlgebra,T,C}}
    equality_data = pauli_su2_translation_orbit_singlet_channel_equalities(
        basis,
        n_sites;
        max_dimension,
        atol,
    )
    for (row_idx, form) in pairs(equality_data.basis_forms)
        poly = _pauli_su2_basis_form_polynomial(form, MP)
        iszero(poly) && continue
        label = (
            feature=:pauli_su2_translation_orbit_singlet_channel_equality,
            row=Int(row_idx),
            n_sites=equality_data.n_sites,
            momentum=equality_data.momentum,
            su2_label=equality_data.row_labels[row_idx],
            coefficient_domain=equality_data.coefficient_domain,
            exact_coefficient_domain=equality_data.exact_coefficient_domain,
        )
        before = length(constraints)
        _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        for constraint_idx in (before + 1):length(constraints)
            zero_origin_by_constraint[constraint_idx] =
                PauliSU2SingletChannelEqualityOriginSeed(label)
        end
    end
    return equality_data
end

function _append_pauli_su2_translation_singlet_channel_equality_constraints!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    zero_origin_by_constraint::Dict{Int,Any},
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    momentum::Integer,
    ::Type{MP};
    max_dimension::Integer,
    atol::Real,
) where {T<:Unsigned,C<:Number,MP<:Polynomial{PauliAlgebra,T,C}}
    equality_data = pauli_su2_translation_orbit_singlet_channel_equalities(
        basis,
        n_sites;
        max_dimension,
        momentum,
        atol,
    )
    for (row_idx, form) in pairs(equality_data.basis_forms)
        poly = _pauli_su2_basis_form_polynomial(form, MP)
        iszero(poly) && continue
        label = (
            feature=:pauli_su2_translation_orbit_singlet_channel_equality,
            decomposition=:translation_su2,
            momentum=equality_data.momentum,
            row=Int(row_idx),
            su2_label=equality_data.row_labels[row_idx],
            coefficient_domain=equality_data.coefficient_domain,
            exact_coefficient_domain=equality_data.exact_coefficient_domain,
        )
        before = length(constraints)
        _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        for constraint_idx in (before + 1):length(constraints)
            zero_origin_by_constraint[constraint_idx] =
                PauliSU2SingletChannelEqualityOriginSeed(label)
        end
    end
    return equality_data
end

function _pauli_su2_basis_moment_block_meta(
    ::Type{M},
    cone::Symbol,
    block,
) where {M<:NormalMonomial{PauliAlgebra}}
    label = (
        feature=:pauli_su2_basis_moment,
        spin2=block.spin2,
        reference_m2=block.reference_m2,
    )
    return BlockMeta{M}(
        cone,
        PauliSU2BasisMomentOrigin(
            label,
            block.row_labels;
            transform=block.transform,
        ),
        fill(one(M), size(block.matrix, 1)),
    )
end

function _pauli_su2_translation_orbit_basis_moment_block_meta(
    ::Type{M},
    cone::Symbol,
    block,
) where {M<:NormalMonomial{PauliAlgebra}}
    label = (
        feature=:pauli_su2_translation_orbit_basis_moment,
        spin2=block.spin2,
        reference_m2=block.reference_m2,
    )
    return BlockMeta{M}(
        cone,
        PauliSU2BasisMomentOrigin(
            label,
            block.row_labels;
            transform=block.transform,
        ),
        fill(one(M), size(block.matrix, 1)),
    )
end

"""
    pauli_su2_basis_moment_problem(objective, basis; assume_su2_invariant=false, ops=nothing, max_dimension=4096, reference_m2=:lowest, singlet_channel_equalities=false, singlet_channel_atol=1e-12)

Build a low-level Pauli `MomentProblem` whose base moment matrix is reduced by
the support-complete SU(2) Pauli-word decomposition.  This helper is intended
for small structural tests and as the integration point for future fast paths;
it does not apply translation, sign, reflection symmetry, or an automatic
translation/reflection proof.  Pass `ops=(σx, σy, σz)` to verify global
Pauli-axis-rotation invariance of the objective, or pass
`assume_su2_invariant=true` only after doing that proof elsewhere.
"""
function pauli_su2_basis_moment_problem(
    objective::Polynomial{PauliAlgebra,T,C},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}};
    assume_su2_invariant::Bool=false,
    ops=nothing,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Integer,C<:Number}
    verified_su2_invariant = assume_su2_invariant
    if ops !== nothing
        _check_pauli_axis_rotation_invariance(objective, ops; context="SU(2) basis objective")
        verified_su2_invariant = true
    end
    verified_su2_invariant || throw(ArgumentError(
        "`pauli_su2_basis_moment_problem` drops SU(2) off-sector and magnetic-copy rows; " *
        "pass `ops=(σx, σy, σz)` to verify global Pauli-axis-rotation invariance or " *
        "`assume_su2_invariant=true` after verifying it elsewhere."
    ))

    M = NormalMonomial{PauliAlgebra,T}
    MP_P = Polynomial{PauliAlgebra,T,ComplexF64}
    basis_vec = M[mono for mono in basis]
    objective_mp = convert(MP_P, objective)
    _, basis_moment_block = _build_constraint_matrix(one(MP_P), basis_vec, :HPSD)

    constraints = Tuple{Symbol,Matrix{MP_P}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()
    for block in pauli_su2_basis_reduced_moment_blocks(
        basis_moment_block,
        basis_vec;
        max_dimension,
        reference_m2,
    )
        before = length(constraints)
        _append_constraint!(constraints, :HPSD, block.matrix, MP_P)
        block_meta_by_constraint[before + 1] =
            _pauli_su2_basis_moment_block_meta(M, :HPSD, block)
    end
    if singlet_channel_equalities
        _append_pauli_su2_singlet_channel_equality_constraints!(
            constraints,
            basis_vec,
            MP_P,
            zero_origin_by_constraint;
            max_dimension,
            atol=singlet_channel_atol,
        )
    end

    constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
    _check_polynomial_moments_covered(
        objective_mp,
        constraint_basis,
        "SU(2) basis objective",
    )
    total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
    return MomentProblem{PauliAlgebra,T,M,MP_P}(
        objective_mp,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=false,
    )
end

"""
    pauli_su2_translation_orbit_basis_moment_problem(objective, moment_block, basis, n_sites; assume_su2_invariant=false, ops=nothing, max_dimension=4096, reference_m2=:lowest, singlet_channel_equalities=false, singlet_channel_atol=1e-12)

Build a low-level Pauli `MomentProblem` from an already constructed
translation-orbit-basis moment block, reducing that block by the SU(2)
Clebsch-Gordan transform on orbit representatives.  This helper does not build
translation momentum products itself; callers must supply the polynomial block
to reduce.
"""
function pauli_su2_translation_orbit_basis_moment_problem(
    objective::Polynomial{PauliAlgebra,T,C},
    moment_block::AbstractMatrix{<:Polynomial{PauliAlgebra,T,C2}},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    assume_su2_invariant::Bool=false,
    ops=nothing,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Integer,C<:Number,C2<:Number}
    verified_su2_invariant = assume_su2_invariant
    if ops !== nothing
        _check_pauli_axis_rotation_invariance(
            objective,
            ops;
            context="translation-orbit SU(2) basis objective",
        )
        verified_su2_invariant = true
    end
    verified_su2_invariant || throw(ArgumentError(
        "`pauli_su2_translation_orbit_basis_moment_problem` drops SU(2) " *
        "off-sector and magnetic-copy rows; pass `ops=(σx, σy, σz)` to " *
        "verify global Pauli-axis-rotation invariance or " *
        "`assume_su2_invariant=true` after verifying it elsewhere."
    ))

    M = NormalMonomial{PauliAlgebra,T}
    MP_P = Polynomial{PauliAlgebra,T,ComplexF64}
    basis_vec = M[mono for mono in basis]
    objective_mp = convert(MP_P, objective)
    block_mp = Matrix{MP_P}(undef, size(moment_block))
    for idx in eachindex(moment_block)
        block_mp[idx] = convert(MP_P, moment_block[idx])
    end

    constraints = Tuple{Symbol,Matrix{MP_P}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()
    for block in pauli_su2_translation_orbit_basis_reduced_moment_blocks(
        block_mp,
        basis_vec,
        n_sites;
        max_dimension,
        reference_m2,
    )
        before = length(constraints)
        _append_constraint!(constraints, :HPSD, block.matrix, MP_P)
        block_meta_by_constraint[before + 1] =
            _pauli_su2_translation_orbit_basis_moment_block_meta(M, :HPSD, block)
    end
    if singlet_channel_equalities
        _append_pauli_su2_translation_orbit_singlet_channel_equality_constraints!(
            constraints,
            basis_vec,
            n_sites,
            MP_P,
            zero_origin_by_constraint;
            max_dimension,
            atol=singlet_channel_atol,
        )
    end

    constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
    _check_polynomial_moments_covered(
        objective_mp,
        constraint_basis,
        "translation-orbit SU(2) basis objective",
    )
    total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
    return MomentProblem{PauliAlgebra,T,M,MP_P}(
        objective_mp,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=false,
    )
end

function _translation_orbit_representative_translates(
    basis::Vector{M},
    n_sites::Integer,
) where {M<:NormalMonomial{PauliAlgebra}}
    n = Int(n_sites)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    isempty(basis) && throw(ArgumentError(
        "Translation-orbit representative translates require a non-empty basis."
    ))
    one_mono = one(first(basis))
    translated = Dict{M,Vector{M}}()
    for rep in basis
        translated[rep] = isone(rep) ?
            fill(one_mono, n) :
            [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
    end
    return translated
end

"""
    pauli_su2_translation_orbit_momentum_blocks(basis, n_sites, momentum; reflection_symmetry=false, max_dimension=4096, reference_m2=:lowest)

Build one translation momentum polynomial block on a translation-orbit Pauli
word basis, then reduce it by the SU(2) Clebsch-Gordan transform on orbit
representatives.  Nonzero momentum sectors drop the identity row, matching the
translation fast path.  This is a structural helper for the future
translation/SU(2) reducer; it does not apply sign, realification, or solver
construction.  With `reflection_symmetry=true`, only reflection-fixed momentum
sectors are accepted and each SU(2) multiplicity block is split by reflection
parity.
"""
function pauli_su2_translation_orbit_momentum_blocks(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    momentum::Integer;
    reflection_symmetry::Bool=false,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
) where {T<:Integer}
    n = Int(n_sites)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    k = mod(Int(momentum), n)
    M = NormalMonomial{PauliAlgebra,T}
    basis_vec = M[mono for mono in basis]
    sector_basis = iszero(k) ? basis_vec : M[mono for mono in basis_vec if !isone(mono)]
    translated = _translation_orbit_representative_translates(sector_basis, n)
    rep_cache = Dict{M,M}()
    MP_P = Polynomial{PauliAlgebra,T,ComplexF64}
    moment_block = _translation_momentum_block(
        sector_basis,
        k,
        n,
        translated,
        rep_cache,
        MP_P,
    )
    blocks = pauli_su2_translation_orbit_basis_reduced_moment_blocks(
        moment_block,
        sector_basis,
        n;
        max_dimension,
        reference_m2,
        momentum=k,
    )
    if reflection_symmetry
        blocks = _pauli_su2_translation_reflection_fixed_momentum_blocks(
            blocks,
            sector_basis,
            k,
            n,
        )
    end
    return [
        merge(
            block,
            (
                momentum=k,
                label=reflection_symmetry ? block.label : (
                    feature=:moment_matrix,
                    decomposition=:translation_su2,
                    momentum=k,
                    spin2=block.spin2,
                ),
            ),
        )
        for block in blocks
    ]
end

function _pauli_su2_translation_zero_label(momentum::Int, reason::Symbol; fields...)
    return (
        feature=:pauli_su2_base_zero,
        decomposition=:translation_su2,
        momentum=momentum,
        reason=reason,
        fields...,
    )
end

function _pauli_su2_m2_rows_or_reference(
    block,
    m2::Int,
    reference_m2::Int,
    ref_rows,
    rows_by_m2=nothing,
)
    m2 == reference_m2 && return ref_rows
    rows_by_m2 === nothing && return _pauli_su2_basis_m2_rows(block, m2)
    return _pauli_su2_basis_m2_rows(block, m2, rows_by_m2)
end

function _append_pauli_su2_translation_zero_matrix!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    mat::AbstractMatrix,
    ::Type{MP};
    phase_atol,
    zero_origin_by_constraint::Dict{Int,Any},
    label,
) where {MP<:Polynomial{PauliAlgebra}}
    for col in axes(mat, 2), row in axes(mat, 1)
        _append_translation_zero_polynomial!(
            constraints,
            mat[row, col],
            MP;
            phase_atol,
            zero_origin_by_constraint,
            label=merge(label, (row=row, col=col)),
        )
    end
    return nothing
end

function _append_pauli_su2_translation_wigner_eckart_constraints!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    block_meta_by_constraint::Dict{Int,BlockMeta{M}},
    zero_origin_by_constraint::Dict{Int,Any},
    moment_block::AbstractMatrix{MP},
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock},
    momentum::Int,
    ::Type{MP},
    ::Type{M};
    reference_m2,
    phase_atol,
    append_psd::Bool=true,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},MP<:Polynomial{PauliAlgebra,T}}
    reduced_blocks = Any[]
    for block_i in eachindex(transform_blocks)
        transform_i = transform_blocks[block_i].transform
        for block_j in eachindex(transform_blocks)
            block_i == block_j && continue
            block_j < block_i && continue
            transform_j = transform_blocks[block_j].transform
            cross = transform_i * moment_block * adjoint(transform_j)
            _append_pauli_su2_translation_zero_matrix!(
                constraints,
                cross,
                MP;
                phase_atol,
                zero_origin_by_constraint,
                label=_pauli_su2_translation_zero_label(
                    momentum,
                    :spin_offblock;
                    row_spin2=transform_blocks[block_i].spin2,
                    col_spin2=transform_blocks[block_j].spin2,
                ),
            )
            _append_pauli_su2_translation_zero_matrix!(
                constraints,
                adjoint(cross),
                MP;
                phase_atol,
                zero_origin_by_constraint,
                label=_pauli_su2_translation_zero_label(
                    momentum,
                    :spin_offblock;
                    row_spin2=transform_blocks[block_j].spin2,
                    col_spin2=transform_blocks[block_i].spin2,
                ),
            )
        end
    end

    for block in transform_blocks
        rows_by_m2 = _pauli_su2_basis_m2_row_lookup(block)
        transformed = block.transform * moment_block * adjoint(block.transform)
        ref_m2, ref_rows, row_labels = _pauli_su2_basis_reference_rows(
            block,
            reference_m2,
            rows_by_m2,
        )
        reference_matrix = transformed[ref_rows, ref_rows]
        m2_values = collect((-block.spin2):2:block.spin2)

        for row_m2 in m2_values, col_m2 in m2_values
            row_m2 == col_m2 && continue
            row_rows = _pauli_su2_m2_rows_or_reference(
                block,
                row_m2,
                ref_m2,
                ref_rows,
                rows_by_m2,
            )
            col_rows = _pauli_su2_m2_rows_or_reference(
                block,
                col_m2,
                ref_m2,
                ref_rows,
                rows_by_m2,
            )
            _append_pauli_su2_translation_zero_matrix!(
                constraints,
                transformed[row_rows, col_rows],
                MP;
                phase_atol,
                zero_origin_by_constraint,
                label=_pauli_su2_translation_zero_label(
                    momentum,
                    :magnetic_offdiag;
                    spin2=block.spin2,
                    row_m2=row_m2,
                    col_m2=col_m2,
                ),
            )
        end

        for m2 in m2_values
            m2 == ref_m2 && continue
            rows = _pauli_su2_basis_m2_rows(block, m2, rows_by_m2)
            for col in 1:block.multiplicity, row in 1:block.multiplicity
                _append_translation_zero_polynomial!(
                    constraints,
                    transformed[rows[row], rows[col]] - reference_matrix[row, col],
                    MP;
                    phase_atol,
                    zero_origin_by_constraint,
                    label=_pauli_su2_translation_zero_label(
                        momentum,
                        :magnetic_copy;
                        spin2=block.spin2,
                        m2=m2,
                        reference_m2=ref_m2,
                        row_multiplicity=row,
                        col_multiplicity=col,
                    ),
                )
            end
        end

        reduced_block = (
            spin2=block.spin2,
            multiplicity=block.multiplicity,
            irrep_dimension=block.irrep_dimension,
            reference_m2=ref_m2,
            row_labels=row_labels,
            matrix=reference_matrix,
            transform=block,
            coefficient_domain=block.coefficient_domain,
            exact_coefficient_domain=block.exact_coefficient_domain,
            label=(
                feature=:moment_matrix,
                decomposition=:translation_su2,
                momentum=momentum,
                spin2=block.spin2,
            ),
        )
        push!(reduced_blocks, reduced_block)
        if append_psd
            _append_constraint!(constraints, :HPSD, reference_matrix, MP)
            block_meta_by_constraint[length(constraints)] =
                _pauli_su2_translation_orbit_momentum_block_meta(M, :HPSD, reduced_block)
        end
    end
    return reduced_blocks
end

function _append_pauli_su2_translation_zero_linear_matrix!(
    builder::MomentLinearBuilder{K,LC,M},
    mat::AbstractMatrix{<:LinearMomentForm};
    phase_atol,
    label,
    registered_key_tokens=nothing,
) where {K,LC,M}
    for col in axes(mat, 2), row in axes(mat, 1)
        _add_translation_zero_linear_form!(
            builder,
            mat[row, col];
            phase_atol,
            label=merge(label, (row=row, col=col)),
            registered_key_tokens,
        )
    end
    return nothing
end

function _append_transformed_pauli_su2_zero_linear_block!(
    builder::MomentLinearBuilder{K,LC,M},
    entries::Matrix{LinearMomentForm{K,C}},
    U_left::AbstractMatrix{<:Number},
    U_right::AbstractMatrix{<:Number};
    phase_atol,
    label,
    registered_key_tokens=nothing,
) where {K,LC,C,M}
    source_size = size(entries)
    size(U_left, 2) == source_size[1] || throw(DimensionMismatch(
        "left transform has $(size(U_left, 2)) columns but linear block has $(source_size[1]) rows"
    ))
    size(U_right, 2) == source_size[2] || throw(DimensionMismatch(
        "right transform has $(size(U_right, 2)) columns but linear block has $(source_size[2]) columns"
    ))

    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(phase_atol) : phase_atol
    for j in axes(U_right, 1), i in axes(U_left, 1)
        pairs = Pair{K,C}[]
        for b in axes(entries, 2)
            right_coeff = conj(U_right[j, b])
            abs(right_coeff) <= tolerance && continue
            for a in axes(entries, 1)
                left_coeff = U_left[i, a]
                abs(left_coeff) <= tolerance && continue
                _add_scaled_linear_form_terms!(
                    pairs,
                    left_coeff * right_coeff,
                    entries[a, b],
                    C;
                    atol=tolerance,
                )
            end
        end
        _add_translation_zero_linear_form!(
            builder,
            _linear_moment_form_from_owned_pairs!(pairs);
            phase_atol,
            label=merge(label, (row=i, col=j)),
            registered_key_tokens,
        )
    end
    return nothing
end

function _sparse_transformed_linear_form(
    entries::Matrix{LinearMomentForm{K,C}},
    left_row,
    right_row,
    ::Type{C};
    atol::Real,
    scratch_acc::Union{Nothing,Dict{K,C}}=nothing,
    scratch_acc_sizehint=nothing,
) where {K,C}
    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(atol) : atol
    pair_count_hint = 0
    for (a, _) in left_row
        for (b, _) in right_row
            pair_count_hint += length(entries[a, b])
        end
    end
    pair_count_hint == 0 && return LinearMomentForm{K,C}(Pair{K,C}[], Val(:trusted))

    if pair_count_hint >= _PAULI_SU2_WIGNER_SPARSE_FORM_ACCUMULATOR_THRESHOLD
        acc = scratch_acc === nothing ? Dict{K,C}() : scratch_acc
        empty!(acc)
        hint = min(pair_count_hint >>> 1, 8192)
        if scratch_acc_sizehint === nothing
            sizehint!(acc, hint)
        elseif hint > scratch_acc_sizehint[]
            sizehint!(acc, hint)
            scratch_acc_sizehint[] = hint
        end
        for (a, left_coeff) in left_row
            for (b, right_coeff) in right_row
                entry_form = entries[a, b]
                isempty(entry_form) && continue
                _add_scaled_linear_form_terms!(
                    acc,
                    left_coeff * conj(right_coeff),
                    entry_form,
                    C;
                    atol=tolerance,
                )
            end
        end
        return _linear_moment_form_from_accumulator!(acc; atol=tolerance)
    end

    pairs = Pair{K,C}[]
    sizehint!(pairs, pair_count_hint)
    for (a, left_coeff) in left_row
        for (b, right_coeff) in right_row
            entry_form = entries[a, b]
            isempty(entry_form) && continue
            _add_scaled_linear_form_terms!(
                pairs,
                left_coeff * conj(right_coeff),
                entry_form,
                C;
                atol=tolerance,
            )
        end
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _sparse_transform_linear_block_entries(
    entries::Matrix{LinearMomentForm{K,C}},
    transform,
    row_indices::AbstractVector{<:Integer},
    col_indices::AbstractVector{<:Integer};
    atol::Real,
) where {K,C}
    transformed = Matrix{LinearMomentForm{K,C}}(undef, length(row_indices), length(col_indices))
    scratch_acc = Dict{K,C}()
    scratch_acc_sizehint = Ref(0)
    for (col_out, col_idx) in pairs(col_indices)
        right_row = transform.rows[Int(col_idx)]
        for (row_out, row_idx) in pairs(row_indices)
            transformed[row_out, col_out] = _sparse_transformed_linear_form(
                entries,
                transform.rows[Int(row_idx)],
                right_row,
                C;
                atol,
                scratch_acc,
                scratch_acc_sizehint,
            )
        end
    end
    return transformed
end

function _indexed_linear_block_entries(
    entries::AbstractMatrix{LinearMomentForm{K,C}},
) where {K,C}
    index_to_key = K[]
    seen = Set{K}()
    for form in entries
        for (key, _) in form
            key in seen && continue
            push!(seen, key)
            push!(index_to_key, key)
        end
    end
    sort!(index_to_key; lt=key_lt)

    key_to_index = Dict{K,Int}()
    sizehint!(key_to_index, length(index_to_key))
    for (idx, key) in pairs(index_to_key)
        key_to_index[key] = Int(idx)
    end

    indexed = Matrix{LinearMomentForm{Int,C}}(undef, size(entries))
    for idx in eachindex(entries)
        form = entries[idx]
        pairs = Pair{Int,C}[]
        sizehint!(pairs, length(form))
        for (key, coef) in form
            push!(pairs, key_to_index[key] => coef)
        end
        indexed[idx] = LinearMomentForm{Int,C}(pairs, Val(:trusted))
    end
    return indexed, index_to_key
end

function _rekey_indexed_linear_form(
    form::LinearMomentForm{Int,C},
    index_to_key::AbstractVector{K},
) where {K,C}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(form))
    for (idx, coef) in form
        1 <= idx <= length(index_to_key) || throw(BoundsError(index_to_key, idx))
        push!(pairs, index_to_key[idx] => coef)
    end
    return LinearMomentForm{K,C}(pairs, Val(:trusted))
end

function _rekey_indexed_linear_block_entries(
    entries::AbstractMatrix{LinearMomentForm{Int,C}},
    index_to_key::AbstractVector{K},
) where {K,C}
    rekeyed = Matrix{LinearMomentForm{K,C}}(undef, size(entries))
    for idx in eachindex(entries)
        rekeyed[idx] = _rekey_indexed_linear_form(entries[idx], index_to_key)
    end
    return rekeyed
end

const _PAULI_SU2_WIGNER_SPARSE_FORM_ACCUMULATOR_THRESHOLD = 1024

function _append_sparse_transformed_pauli_su2_zero_linear_block!(
    builder::MomentLinearBuilder{K,LC,M},
    entries::Matrix{LinearMomentForm{EK,C}},
    U_left,
    U_right;
    phase_atol,
    label,
    registered_key_tokens=nothing,
    stage_times_ns=nothing,
    stage_prefix::Symbol=:spin,
    register_keys::Bool=true,
    adjoint_key_cache=nothing,
    entry_index_to_key=nothing,
) where {K,LC,EK,C,M}
    source_size = size(entries)
    size(U_left, 2) == source_size[1] || throw(DimensionMismatch(
        "left transform has $(size(U_left, 2)) columns but linear block has $(source_size[1]) rows"
    ))
    size(U_right, 2) == source_size[2] || throw(DimensionMismatch(
        "right transform has $(size(U_right, 2)) columns but linear block has $(source_size[2]) columns"
    ))

    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(phase_atol) : phase_atol
    form_build_ns = 0
    form_hint_ns = 0
    form_accumulate_ns = 0
    form_reduce_ns = 0
    zero_append_ns = 0
    scratch_acc = Dict{EK,C}()
    scratch_acc_sizehint = 0
    for j in eachindex(U_right.rows), i in eachindex(U_left.rows)
        stage_start_ns = time_ns()
        pair_count_hint = 0
        for (a, _) in U_left.rows[i]
            for (b, _) in U_right.rows[j]
                pair_count_hint += length(entries[a, b])
            end
        end
        form_hint_ns += Int(time_ns() - stage_start_ns)
        pair_count_hint == 0 && continue

        stage_start_ns = time_ns()
        form = if pair_count_hint >= _PAULI_SU2_WIGNER_SPARSE_FORM_ACCUMULATOR_THRESHOLD
            empty!(scratch_acc)
            hint = min(pair_count_hint >>> 1, 8192)
            if hint > scratch_acc_sizehint
                sizehint!(scratch_acc, hint)
                scratch_acc_sizehint = hint
            end
            for (a, left_coeff) in U_left.rows[i]
                for (b, right_coeff) in U_right.rows[j]
                    entry_form = entries[a, b]
                    isempty(entry_form) && continue
                    _add_scaled_linear_form_terms!(
                        scratch_acc,
                        left_coeff * conj(right_coeff),
                        entry_form,
                        C;
                        atol=tolerance,
                    )
                end
            end
            form_accumulate_ns += Int(time_ns() - stage_start_ns)

            stage_start_ns = time_ns()
            reduced = _linear_moment_form_from_accumulator!(scratch_acc; atol=tolerance)
            form_reduce_ns += Int(time_ns() - stage_start_ns)
            reduced
        else
            pairs = Pair{EK,C}[]
            sizehint!(pairs, pair_count_hint)
            for (a, left_coeff) in U_left.rows[i]
                for (b, right_coeff) in U_right.rows[j]
                    entry_form = entries[a, b]
                    isempty(entry_form) && continue
                    _add_scaled_linear_form_terms!(
                        pairs,
                        left_coeff * conj(right_coeff),
                        entry_form,
                        C;
                        atol=tolerance,
                    )
                end
            end
            form_accumulate_ns += Int(time_ns() - stage_start_ns)

            stage_start_ns = time_ns()
            reduced = _linear_moment_form_from_owned_pairs!(pairs)
            form_reduce_ns += Int(time_ns() - stage_start_ns)
            reduced
        end

        stage_start_ns = time_ns()
        if !isempty(form)
            if entry_index_to_key === nothing
                _add_translation_zero_linear_form!(
                    builder,
                    form;
                    phase_atol,
                    label=merge(label, (row=i, col=j)),
                    registered_key_tokens,
                    register_keys,
                    adjoint_key_cache,
                )
            else
                _add_translation_zero_indexed_linear_form!(
                    builder,
                    form,
                    entry_index_to_key;
                    phase_atol,
                    label=merge(label, (row=i, col=j)),
                    registered_key_tokens,
                    register_keys,
                    adjoint_key_cache,
                )
            end
        end
        zero_append_ns += Int(time_ns() - stage_start_ns)
    end
    form_build_ns = form_hint_ns + form_accumulate_ns + form_reduce_ns
    if stage_times_ns !== nothing
        form_build_key = Symbol("su2_wigner_", stage_prefix, "_form_build")
        form_hint_key = Symbol("su2_wigner_", stage_prefix, "_form_hint")
        form_accumulate_key = Symbol("su2_wigner_", stage_prefix, "_form_accumulate")
        form_reduce_key = Symbol("su2_wigner_", stage_prefix, "_form_reduce")
        zero_append_key = Symbol("su2_wigner_", stage_prefix, "_zero_append")
        stage_times_ns[form_build_key] =
            get(stage_times_ns, form_build_key, 0) + form_build_ns
        stage_times_ns[form_hint_key] =
            get(stage_times_ns, form_hint_key, 0) + form_hint_ns
        stage_times_ns[form_accumulate_key] =
            get(stage_times_ns, form_accumulate_key, 0) + form_accumulate_ns
        stage_times_ns[form_reduce_key] =
            get(stage_times_ns, form_reduce_key, 0) + form_reduce_ns
        stage_times_ns[zero_append_key] =
            get(stage_times_ns, zero_append_key, 0) + zero_append_ns
    end
    return (form_build_ns=form_build_ns, zero_append_ns=zero_append_ns)
end

function _append_pauli_su2_translation_wigner_eckart_linear!(
    builder::MomentLinearBuilder{K,LC,M},
    moment_entries::Matrix{LinearMomentForm{K,C}},
    transform_blocks::AbstractVector{PauliSU2BasisTransformBlock},
    momentum::Int;
    reference_m2,
    phase_atol,
    stage_times_ns::Union{Nothing,Dict{Symbol,Int}}=nothing,
    registered_key_tokens=nothing,
    register_keys::Bool=true,
    emit_invariance_rows::Bool=true,
) where {K,LC,C,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    reduced_blocks = Any[]
    adjoint_key_cache = Dict{K,K}()
    stage_start_ns = time_ns()
    indexed_moment_entries, entry_index_to_key = _indexed_linear_block_entries(moment_entries)
    if stage_times_ns !== nothing
        stage_times_ns[:su2_wigner_index_entries] =
            get(stage_times_ns, :su2_wigner_index_entries, 0) +
            Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    sparse_transforms = [
        _pauli_sparse_transform_rows(block.transform; atol=phase_atol)
        for block in transform_blocks
    ]
    if stage_times_ns !== nothing
        stage_times_ns[:su2_wigner_sparse_transforms] =
            get(stage_times_ns, :su2_wigner_sparse_transforms, 0) +
            Int(time_ns() - stage_start_ns)
    end

    if emit_invariance_rows
        for block_i in eachindex(transform_blocks)
            transform_i = sparse_transforms[block_i]
            for block_j in eachindex(transform_blocks)
                block_i == block_j && continue
                block_j < block_i && continue
                transform_j = sparse_transforms[block_j]
                stage_start_ns = time_ns()
                spin_pair_stage_times = _append_sparse_transformed_pauli_su2_zero_linear_block!(
                    builder,
                    indexed_moment_entries,
                    transform_i,
                    transform_j;
                    phase_atol,
                    label=_pauli_su2_translation_zero_label(
                        momentum,
                        :spin_offblock;
                        row_spin2=transform_blocks[block_i].spin2,
                        col_spin2=transform_blocks[block_j].spin2,
                    ),
                    registered_key_tokens,
                    stage_times_ns,
                    register_keys,
                    adjoint_key_cache,
                    entry_index_to_key,
                )
                if stage_times_ns !== nothing
                    elapsed_ns = Int(time_ns() - stage_start_ns)
                    stage_times_ns[:su2_wigner_spin_stream] =
                        get(stage_times_ns, :su2_wigner_spin_stream, 0) + elapsed_ns
                    spin_pair_key = Symbol(
                        "su2_wigner_spin_pair_",
                        transform_blocks[block_i].spin2,
                        "_",
                        transform_blocks[block_j].spin2,
                    )
                    spin_pair_form_build_key = Symbol(spin_pair_key, "_form_build")
                    spin_pair_zero_append_key = Symbol(spin_pair_key, "_zero_append")
                    stage_times_ns[spin_pair_key] =
                        get(stage_times_ns, spin_pair_key, 0) + elapsed_ns
                    stage_times_ns[spin_pair_form_build_key] =
                        get(stage_times_ns, spin_pair_form_build_key, 0) +
                        spin_pair_stage_times.form_build_ns
                    stage_times_ns[spin_pair_zero_append_key] =
                        get(stage_times_ns, spin_pair_zero_append_key, 0) +
                        spin_pair_stage_times.zero_append_ns
                end
            end
        end
    end

    for block_idx in eachindex(transform_blocks)
        block = transform_blocks[block_idx]
        transform = sparse_transforms[block_idx]
        rows_by_m2 = _pauli_su2_basis_m2_row_lookup(block)
        stage_start_ns = time_ns()
        ref_m2, ref_rows, row_labels = _pauli_su2_basis_reference_rows(
            block,
            reference_m2,
            rows_by_m2,
        )
        indexed_reference_matrix = _sparse_transform_linear_block_entries(
            indexed_moment_entries,
            transform,
            ref_rows,
            ref_rows;
            atol=phase_atol,
        )
        reference_matrix = _rekey_indexed_linear_block_entries(
            indexed_reference_matrix,
            entry_index_to_key,
        )
        if stage_times_ns !== nothing
            stage_times_ns[:su2_wigner_diag_transform] =
                get(stage_times_ns, :su2_wigner_diag_transform, 0) +
                Int(time_ns() - stage_start_ns)
        end
        if emit_invariance_rows
            m2_values = collect((-block.spin2):2:block.spin2)
            sparse_transform_by_m2 = Dict{Int,typeof(transform)}()
            sizehint!(sparse_transform_by_m2, length(m2_values))
            for m2 in m2_values
                rows = _pauli_su2_m2_rows_or_reference(block, m2, ref_m2, ref_rows, rows_by_m2)
                sparse_transform_by_m2[m2] = _select_sparse_transform_rows(transform, rows)
            end

            for row_m2 in m2_values, col_m2 in m2_values
                row_m2 == col_m2 && continue
                col_m2 < row_m2 && continue
                row_transform = sparse_transform_by_m2[row_m2]
                col_transform = sparse_transform_by_m2[col_m2]
                stage_start_ns = time_ns()
                _append_sparse_transformed_pauli_su2_zero_linear_block!(
                    builder,
                    indexed_moment_entries,
                    row_transform,
                    col_transform;
                    phase_atol,
                    label=_pauli_su2_translation_zero_label(
                        momentum,
                        :magnetic_offdiag;
                        spin2=block.spin2,
                        row_m2=row_m2,
                        col_m2=col_m2,
                    ),
                    registered_key_tokens,
                    stage_times_ns,
                    stage_prefix=:magnetic_offdiag,
                    register_keys,
                    adjoint_key_cache,
                    entry_index_to_key,
                )
                if stage_times_ns !== nothing
                    stage_times_ns[:su2_wigner_magnetic_offdiag_append] =
                        get(stage_times_ns, :su2_wigner_magnetic_offdiag_append, 0) +
                        Int(time_ns() - stage_start_ns)
                end
            end

            magnetic_copy_scratch_acc = Dict{Int,C}()
            magnetic_copy_scratch_sizehint = Ref(0)
            for m2 in m2_values
                m2 == ref_m2 && continue
                row_transform = sparse_transform_by_m2[m2]
                for col in 1:block.multiplicity, row in 1:block.multiplicity
                    row > col && continue
                    stage_total_start_ns = time_ns()
                    stage_start_ns = stage_total_start_ns
                    form = _sparse_transformed_linear_form(
                        indexed_moment_entries,
                        row_transform.rows[row],
                        row_transform.rows[col],
                        C;
                        atol=phase_atol,
                        scratch_acc=magnetic_copy_scratch_acc,
                        scratch_acc_sizehint=magnetic_copy_scratch_sizehint,
                    )
                    if stage_times_ns !== nothing
                        stage_times_ns[:su2_wigner_magnetic_copy_form_build] =
                            get(stage_times_ns, :su2_wigner_magnetic_copy_form_build, 0) +
                            Int(time_ns() - stage_start_ns)
                    end
                    stage_start_ns = time_ns()
                    copy_form = _subtract_linear_forms(
                        form,
                        indexed_reference_matrix[row, col];
                        atol=phase_atol,
                    )
                    if stage_times_ns !== nothing
                        stage_times_ns[:su2_wigner_magnetic_copy_form_subtract] =
                            get(stage_times_ns, :su2_wigner_magnetic_copy_form_subtract, 0) +
                            Int(time_ns() - stage_start_ns)
                    end
                    stage_start_ns = time_ns()
                    if !isempty(copy_form)
                        _add_translation_zero_indexed_linear_form!(
                            builder,
                            copy_form,
                            entry_index_to_key;
                            phase_atol,
                            label=_pauli_su2_translation_zero_label(
                                momentum,
                                :magnetic_copy;
                                spin2=block.spin2,
                                m2=m2,
                                reference_m2=ref_m2,
                                row_multiplicity=row,
                                col_multiplicity=col,
                            ),
                            registered_key_tokens,
                            register_keys,
                            adjoint_key_cache,
                        )
                    end
                    if stage_times_ns !== nothing
                        stage_times_ns[:su2_wigner_magnetic_copy_zero_append] =
                            get(stage_times_ns, :su2_wigner_magnetic_copy_zero_append, 0) +
                            Int(time_ns() - stage_start_ns)
                        stage_times_ns[:su2_wigner_magnetic_copy_append] =
                            get(stage_times_ns, :su2_wigner_magnetic_copy_append, 0) +
                            Int(time_ns() - stage_total_start_ns)
                    end
                end
            end
        end

        push!(
            reduced_blocks,
            (
                spin2=block.spin2,
                multiplicity=block.multiplicity,
                irrep_dimension=block.irrep_dimension,
                reference_m2=ref_m2,
                row_labels=row_labels,
                matrix=reference_matrix,
                transform=block,
                coefficient_domain=block.coefficient_domain,
                exact_coefficient_domain=block.exact_coefficient_domain,
                label=(
                    feature=:moment_matrix,
                    decomposition=:translation_su2,
                    momentum=momentum,
                    spin2=block.spin2,
                ),
            ),
        )
    end
    return reduced_blocks
end

function _pauli_su2_translation_wigner_sector_blocks!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    block_meta_by_constraint::Dict{Int,BlockMeta{M}},
    zero_origin_by_constraint::Dict{Int,Any},
    basis_vec::Vector{M},
    n::Int,
    k::Int,
    ::Type{MP},
    ::Type{M};
    max_dimension::Integer,
    reference_m2,
    phase_atol,
    append_psd::Bool=true,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},MP<:Polynomial{PauliAlgebra,T}}
    sector_basis = iszero(k) ? basis_vec : M[mono for mono in basis_vec if !isone(mono)]
    translated = _translation_orbit_representative_translates(sector_basis, n)
    moment_block = _translation_momentum_block(
        sector_basis,
        k,
        n,
        translated,
        Dict{M,M}(),
        MP,
    )
    transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
        sector_basis,
        n;
        max_dimension,
        momentum=k,
    )
    reduced_blocks = _append_pauli_su2_translation_wigner_eckart_constraints!(
        constraints,
        block_meta_by_constraint,
        zero_origin_by_constraint,
        moment_block,
        transform_blocks,
        k,
        MP,
        M;
        reference_m2,
        phase_atol,
        append_psd,
    )
    if singlet_channel_equalities
        _append_pauli_su2_translation_singlet_channel_equality_constraints!(
            constraints,
            zero_origin_by_constraint,
            sector_basis,
            n,
            k,
            MP;
            max_dimension,
            atol=singlet_channel_atol,
        )
    end
    return sector_basis, reduced_blocks
end

"""
    pauli_su2_translation_orbit_wigner_eckart_moment_problem(objective, basis, n_sites, momentum; assume_su2_invariant=false, ops=nothing, max_dimension=4096, reference_m2=:lowest, phase_atol=1e-8)

Build one complex translation/SU(2) momentum sector while explicitly emitting
the Wigner-Eckart zero/copy rows that the reduced base-moment helper currently
drops.  This is a low-level structural/debug constructor; it does not yet handle
realification, reflection pairing, RDM/state-opt blocks, or solver dispatch.
"""
function pauli_su2_translation_orbit_wigner_eckart_moment_problem(
    objective::Polynomial{PauliAlgebra,T,C},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    momentum::Integer;
    assume_su2_invariant::Bool=false,
    ops=nothing,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    phase_atol::Real=1.0e-8,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Unsigned,C<:Number}
    verified_su2_invariant = assume_su2_invariant
    if ops !== nothing
        assume_su2_invariant || _check_pauli_axis_rotation_invariance(
            objective,
            ops;
            context="translation-orbit SU(2) Wigner-Eckart objective",
        )
        verified_su2_invariant = true
    end
    verified_su2_invariant || throw(ArgumentError(
        "`pauli_su2_translation_orbit_wigner_eckart_moment_problem` requires " *
        "a verified SU(2)-invariant objective."
    ))

    n = Int(n_sites)
    _check_translation_invariance(
        objective,
        n;
        context="translation-orbit SU(2) Wigner-Eckart objective",
    )
    k = mod(Int(momentum), n)

    M = NormalMonomial{PauliAlgebra,T}
    MP = Polynomial{PauliAlgebra,T,ComplexF64}
    basis_vec = M[mono for mono in basis]
    sector_basis = iszero(k) ? basis_vec : M[mono for mono in basis_vec if !isone(mono)]
    translated = _translation_orbit_representative_translates(sector_basis, n)
    moment_block = _translation_momentum_block(
        sector_basis,
        k,
        n,
        translated,
        Dict{M,M}(),
        MP,
    )
    transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
        sector_basis,
        n;
        max_dimension,
        momentum=k,
    )

    constraints = Tuple{Symbol,Matrix{MP}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()
    _append_pauli_su2_translation_wigner_eckart_constraints!(
        constraints,
        block_meta_by_constraint,
        zero_origin_by_constraint,
        moment_block,
        transform_blocks,
        k,
        MP,
        M;
        reference_m2,
        phase_atol,
    )
    if singlet_channel_equalities
        _append_pauli_su2_translation_singlet_channel_equality_constraints!(
            constraints,
            zero_origin_by_constraint,
            sector_basis,
            n,
            k,
            MP;
            max_dimension,
            atol=singlet_channel_atol,
        )
    end

    objective_mp = convert(MP, _translation_orbit_reduce_polynomial(objective, n))
    constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
    _check_polynomial_moments_covered(
        objective_mp,
        constraint_basis,
        "translation-orbit SU(2) Wigner-Eckart objective",
    )
    total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
    return MomentProblem{PauliAlgebra,T,M,MP}(
        objective_mp,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=false,
    )
end

"""
    pauli_su2_translation_orbit_wigner_eckart_moment_problem(objective, basis, n_sites; assume_su2_invariant=false, ops=nothing, momenta=nothing, max_dimension=4096, reference_m2=:lowest, phase_atol=1e-8)

Build all selected complex translation/SU(2) momentum sectors while explicitly
emitting the Wigner-Eckart zero/copy rows.  This is the all-sector counterpart
of the single-momentum debug constructor; it still does not perform
realification, reflection pairing, RDM/state-opt blocks, or solver dispatch.
"""
function pauli_su2_translation_orbit_wigner_eckart_moment_problem(
    objective::Polynomial{PauliAlgebra,T,C},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    assume_su2_invariant::Bool=false,
    ops=nothing,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    phase_atol::Real=1.0e-8,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Unsigned,C<:Number}
    verified_su2_invariant = assume_su2_invariant
    if ops !== nothing
        assume_su2_invariant || _check_pauli_axis_rotation_invariance(
            objective,
            ops;
            context="translation-orbit SU(2) Wigner-Eckart objective",
        )
        verified_su2_invariant = true
    end
    verified_su2_invariant || throw(ArgumentError(
        "`pauli_su2_translation_orbit_wigner_eckart_moment_problem` requires " *
        "a verified SU(2)-invariant objective."
    ))

    n = Int(n_sites)
    _check_translation_invariance(
        objective,
        n;
        context="translation-orbit SU(2) Wigner-Eckart objective",
    )
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix=false)
    0 in sectors || throw(ArgumentError(
        "All-sector translation/SU(2) Wigner-Eckart construction requires " *
        "momentum sector 0 because it carries the objective moments."
    ))

    M = NormalMonomial{PauliAlgebra,T}
    MP = Polynomial{PauliAlgebra,T,ComplexF64}
    basis_vec = M[mono for mono in basis]
    constraints = Tuple{Symbol,Matrix{MP}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()

    for k in sectors
        sector_basis = iszero(k) ? basis_vec : M[mono for mono in basis_vec if !isone(mono)]
        translated = _translation_orbit_representative_translates(sector_basis, n)
        moment_block = _translation_momentum_block(
            sector_basis,
            k,
            n,
            translated,
            Dict{M,M}(),
            MP,
        )
        transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
            sector_basis,
            n;
            max_dimension,
            momentum=k,
        )
        _append_pauli_su2_translation_wigner_eckart_constraints!(
            constraints,
            block_meta_by_constraint,
            zero_origin_by_constraint,
            moment_block,
            transform_blocks,
            k,
            MP,
            M;
            reference_m2,
            phase_atol,
        )
        if singlet_channel_equalities
            _append_pauli_su2_translation_singlet_channel_equality_constraints!(
                constraints,
                zero_origin_by_constraint,
                sector_basis,
                n,
                k,
                MP;
                max_dimension,
                atol=singlet_channel_atol,
            )
        end
    end

    objective_mp = convert(MP, _translation_orbit_reduce_polynomial(objective, n))
    constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
    _check_polynomial_moments_covered(
        objective_mp,
        constraint_basis,
        "translation-orbit SU(2) Wigner-Eckart objective",
    )
    total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
    return MomentProblem{PauliAlgebra,T,M,MP}(
        objective_mp,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=false,
    )
end

"""
    pauli_su2_translation_orbit_zero_momentum_blocks(basis, n_sites; max_dimension=4096, reference_m2=:lowest)

Build the zero-momentum translation polynomial block on a translation-orbit
Pauli word basis, then reduce it by the SU(2) Clebsch-Gordan transform on orbit
representatives.
"""
function pauli_su2_translation_orbit_zero_momentum_blocks(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    max_dimension::Integer=4096,
    reference_m2=:lowest,
) where {T<:Integer}
    return pauli_su2_translation_orbit_momentum_blocks(
        basis,
        n_sites,
        0;
        max_dimension,
        reference_m2,
    )
end

"""
    pauli_su2_translation_orbit_momentum_block_bundle(basis, n_sites; real_moment_matrix=true, momenta=nothing, reflection_symmetry=false, max_dimension=4096, reference_m2=:lowest)

Build SU(2)-reduced translation momentum blocks for the same canonical
momentum-sector list used by the translation fast path.  The result is a
structural bundle with stable momentum/spin labels and logical block-size
histograms.  It deliberately does not apply sign symmetry or solver
construction.  Reflection splitting keeps fixed sectors as complex HPSD blocks
and, when `real_moment_matrix=true`, emits non-fixed conjugate-pair sectors as
real reflection-parity PSD blocks.
"""
function pauli_su2_translation_orbit_momentum_block_bundle(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    real_moment_matrix::Bool=true,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    reflection_symmetry::Bool=false,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
) where {T<:Integer}
    n = Int(n_sites)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    real_moment_matrix || 0 in sectors || throw(ArgumentError(
        "Momentum sector 0 is required because it carries the normalized identity moment."
    ))
    if reflection_symmetry && !real_moment_matrix
        nonfixed = Int[k for k in sectors if !_translation_reflection_fixed_momentum(k, n)]
        isempty(nonfixed) || throw(ArgumentError(
            "Translation/SU(2) complex reflection reduction requires " *
            "reflection-fixed momentum sectors; got non-fixed sector(s) $nonfixed."
        ))
    end

    blocks = Any[]
    for k in sectors
        sector_blocks = if reflection_symmetry && !_translation_reflection_fixed_momentum(k, n)
            T <: Unsigned || throw(ArgumentError(
                "Translation/SU(2) real reflection sectors require unsigned Pauli indices."
            ))
            pauli_su2_translation_orbit_real_reflection_momentum_blocks(
                basis,
                n,
                k;
                max_dimension,
            )
        else
            pauli_su2_translation_orbit_momentum_blocks(
                basis,
                n,
                k;
                reflection_symmetry,
                max_dimension,
                reference_m2,
            )
        end
        append!(blocks, sector_blocks)
    end

    logical_sizes = Int[block.multiplicity for block in blocks]
    block_labels = Any[block.label for block in blocks]
    block_coefficient_domains = [block.coefficient_domain for block in blocks]
    block_exact_coefficient_domains = [block.exact_coefficient_domain for block in blocks]
    return (
        n_sites=n,
        momentum_sectors=sectors,
        momentum_sector_count=length(sectors),
        real_moment_matrix=real_moment_matrix,
        blocks=blocks,
        block_labels=block_labels,
        logical_block_sizes=logical_sizes,
        n_blocks=length(blocks),
        logical_max_block=isempty(logical_sizes) ? 0 : maximum(logical_sizes),
        logical_total_block_side=sum(logical_sizes; init=0),
        logical_block_histogram=_histogram_pairs(logical_sizes),
        logical_block_feature_histogram=_label_size_histogram(block_labels, logical_sizes),
        block_coefficient_domains=block_coefficient_domains,
        block_exact_coefficient_domains=block_exact_coefficient_domains,
        block_coefficient_domain_histogram=_value_histogram_pairs(block_coefficient_domains),
        block_exact_coefficient_domain_histogram=_value_histogram_pairs(block_exact_coefficient_domains),
        requires_realification=real_moment_matrix,
        reflection_combined=reflection_symmetry,
    )
end

function _pauli_su2_translation_reflection_matrix(
    basis::Vector{M},
    k::Int,
    n::Int,
    ;
    require_fixed::Bool=false,
) where {M<:NormalMonomial{PauliAlgebra}}
    require_fixed && !_translation_reflection_fixed_momentum(k, n) && throw(ArgumentError(
        "Translation/SU(2) reflection reduction currently requires " *
        "reflection-fixed momentum sectors; got k=$k for n_sites=$n."
    ))
    basis_index = Dict{M,Int}(mono => idx for (idx, mono) in pairs(basis))
    mat = zeros(ComplexF64, length(basis), length(basis))
    for (idx, mono) in pairs(basis)
        phase, image = _translation_reflection_image(mono, k, n)
        image_idx = get(basis_index, image, 0)
        image_idx == 0 && throw(ArgumentError(
            "Reflection image $image of translation/SU(2) orbit row $mono " *
            "is missing from the momentum basis."
        ))
        mat[image_idx, idx] = phase
    end
    return mat
end

function _pauli_su2_reflection_row_label(parity::Symbol, base_labels, coeffs; atol::Float64)
    terms = Any[]
    for (label, coeff) in zip(base_labels, coeffs)
        abs(coeff) <= atol && continue
        push!(terms, (coefficient=ComplexF64(coeff), label=label))
    end
    return (
        feature=:pauli_su2_translation_reflection_adapted_row,
        reflection_parity=parity,
        terms=terms,
    )
end

function _pauli_su2_translation_reflection_transform(
    block,
    row_basis::Matrix{ComplexF64},
    momentum::Int,
    parity::Symbol,
)
    return (
        family=:pauli_translation_su2_reflection,
        coefficient_domain=:complex_algebraic_float64,
        exact_coefficient_domain=:complex_sqrt_rational,
        momentum=momentum,
        spin2=block.spin2,
        reference_m2=block.reference_m2,
        reflection_parity=parity,
        base=block.transform,
        row_basis=row_basis,
    )
end

function _pauli_su2_translation_reflection_fixed_momentum_blocks(
    blocks,
    sector_basis::Vector{M},
    k::Int,
    n::Int;
    atol::Float64=1.0e-8,
) where {M<:NormalMonomial{PauliAlgebra}}
    reflection_matrix = _pauli_su2_translation_reflection_matrix(
        sector_basis,
        k,
        n;
        require_fixed=true,
    )
    adapted_blocks = Any[]
    for block in blocks
        transform_block = block.transform
        ref_rows = _pauli_su2_basis_m2_rows(transform_block, block.reference_m2)
        reference_transform = transform_block.transform[ref_rows, :]
        reduced_reflection =
            reference_transform * reflection_matrix * adjoint(reference_transform)
        hermitian_reflection =
            Hermitian((reduced_reflection + adjoint(reduced_reflection)) / 2)
        reflection_residual = max(
            _max_abs_entry(reduced_reflection - adjoint(reduced_reflection)),
            _max_abs_entry(reduced_reflection * reduced_reflection - I),
        )
        reflection_residual <= 100 * atol || throw(ArgumentError(
            "Translation/SU(2) reflection action for momentum $k, " *
            "spin2=$(block.spin2) is not a stable Hermitian involution; " *
            "residual=$reflection_residual."
        ))

        eig = LinearAlgebra.eigen(hermitian_reflection)
        for (parity, sign) in ((:plus, 1.0), (:minus, -1.0))
            cols = findall(value -> abs(value - sign) <= 100 * atol, eig.values)
            isempty(cols) && continue
            vectors = eig.vectors[:, cols]
            row_basis = Matrix{ComplexF64}(adjoint(vectors))
            adapted_matrix = row_basis * block.matrix * adjoint(row_basis)
            row_labels = Any[
                _pauli_su2_reflection_row_label(
                    parity,
                    block.row_labels,
                    row_basis[row_idx, :];
                    atol,
                )
                for row_idx in axes(row_basis, 1)
            ]
            push!(
                adapted_blocks,
                (
                    spin2=block.spin2,
                    multiplicity=size(adapted_matrix, 1),
                    irrep_dimension=block.irrep_dimension,
                    reference_m2=block.reference_m2,
                    row_labels=row_labels,
                    matrix=adapted_matrix,
                    transform=_pauli_su2_translation_reflection_transform(
                        block,
                        row_basis,
                        k,
                        parity,
                    ),
                    coefficient_domain=:complex_algebraic_float64,
                    exact_coefficient_domain=:complex_sqrt_rational,
                    label=(
                        feature=:moment_matrix,
                        decomposition=:translation_su2_reflection,
                        momentum=k,
                        reflection_fixed=true,
                        spin2=block.spin2,
                        reflection_parity=parity,
                    ),
                ),
            )
        end
    end
    return adapted_blocks
end

function _pauli_su2_real_reflection_row_label(
    parity::Symbol,
    base_labels,
    coeffs;
    atol::Float64,
)
    n_complex = length(base_labels)
    terms = Any[]
    for (idx, coeff) in pairs(coeffs)
        abs(coeff) <= atol && continue
        label_idx = idx <= n_complex ? idx : idx - n_complex
        part = idx <= n_complex ? :real : :imag
        push!(
            terms,
            (
                coefficient=Float64(coeff),
                part=part,
                label=base_labels[label_idx],
            ),
        )
    end
    return (
        feature=:pauli_su2_translation_real_reflection_adapted_row,
        reflection_parity=parity,
        terms=terms,
    )
end

function _pauli_su2_translation_real_reflection_transform(
    block,
    row_basis::Matrix{Float64},
    momentum::Int,
    parity::Symbol,
)
    return (
        family=:pauli_translation_su2_real_reflection,
        coefficient_domain=:complex_algebraic_float64,
        exact_coefficient_domain=:complex_sqrt_rational,
        momentum=momentum,
        spin2=block.spin2,
        reference_m2=block.reference_m2,
        reflection_parity=parity,
        base=block.transform,
        row_basis=row_basis,
    )
end

function _pauli_su2_translation_real_reflection_row_bases(
    block,
    reflection_matrix::Matrix{ComplexF64};
    atol::Float64,
)
    transform_block = block.transform
    ref_rows = _pauli_su2_basis_m2_rows(transform_block, block.reference_m2)
    reference_transform = transform_block.transform[ref_rows, :]
    antiunitary_reflection =
        reference_transform * reflection_matrix * transpose(reference_transform)
    re = real.(antiunitary_reflection)
    im = imag.(antiunitary_reflection)
    real_action = [re im; im -re]
    action_residual = max(
        _max_abs_entry(real_action - transpose(real_action)),
        _max_abs_entry(real_action * real_action - I),
    )
    action_residual <= 100 * atol || throw(ArgumentError(
        "Translation/SU(2) real reflection action for spin2=$(block.spin2) " *
        "is not a stable real involution; residual=$action_residual."
    ))

    eig = LinearAlgebra.eigen(Symmetric((real_action + transpose(real_action)) / 2))
    row_blocks = NamedTuple[]
    for (parity, sign) in ((:plus, 1.0), (:minus, -1.0))
        cols = findall(value -> abs(value - sign) <= 100 * atol, eig.values)
        isempty(cols) && continue
        vectors = eig.vectors[:, cols]
        row_basis = Matrix{Float64}(transpose(vectors))
        row_labels = Any[
            _pauli_su2_real_reflection_row_label(
                parity,
                block.row_labels,
                row_basis[row_idx, :];
                atol,
            )
            for row_idx in axes(row_basis, 1)
        ]
        push!(
            row_blocks,
            (
                reflection_parity=parity,
                row_basis=row_basis,
                row_labels=row_labels,
            ),
        )
    end
    return row_blocks
end

"""
    pauli_su2_translation_orbit_real_reflection_momentum_blocks(basis, n_sites, momentum; max_dimension=4096, reference_m2=:zero, phase_atol=1e-8)

Build real PSD blocks for one non-reflection-fixed translation momentum sector
after SU(2) reduction and reflection pairing.  This is the low-level
counterpart of the ordinary translation real-reflection adapter; it constructs
the antiunitary reflection action on each SU(2) multiplicity space and splits
the realified block into reflection parity sectors.  The base moment-matrix
constructors use this helper for all-sector real SU(2)+reflection reduction;
RDM and state-opt SU(2) combinations are separate reducers.
"""
function pauli_su2_translation_orbit_real_reflection_momentum_blocks(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    momentum::Integer;
    max_dimension::Integer=4096,
    reference_m2=:zero,
    phase_atol::Real=1.0e-8,
) where {T<:Unsigned}
    n = Int(n_sites)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    k = mod(Int(momentum), n)
    _translation_reflection_fixed_momentum(k, n) && throw(ArgumentError(
        "Real translation/SU(2) reflection pairing requires a non-fixed " *
        "momentum sector; got k=$k for n_sites=$n."
    ))

    M = NormalMonomial{PauliAlgebra,T}
    basis_vec = M[mono for mono in basis]
    sector_basis = M[mono for mono in basis_vec if !isone(mono)]
    complex_blocks = pauli_su2_translation_orbit_momentum_blocks(
        sector_basis,
        n,
        k;
        max_dimension,
        reference_m2,
    )
    return _pauli_su2_translation_real_reflection_momentum_blocks_from_reduced(
        complex_blocks,
        sector_basis,
        k,
        n;
        phase_atol,
    )
end

function _pauli_su2_translation_real_reflection_momentum_blocks_from_reduced(
    complex_blocks,
    sector_basis::Vector{M},
    k::Int,
    n::Int;
    phase_atol::Real=1.0e-8,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    reflection_matrix = _pauli_su2_translation_reflection_matrix(sector_basis, k, n)
    MP_R = Polynomial{PauliAlgebra,T,Float64}
    atol = Float64(phase_atol)
    real_blocks = Any[]
    for block in complex_blocks
        realified_matrix = _realify_hermitian_block_full(block.matrix, MP_R; atol)
        for reflected in _pauli_su2_translation_real_reflection_row_bases(
            block,
            reflection_matrix;
            atol,
        )
            adapted_matrix = _symmetrize_real_polynomial_block(
                _transform_polynomial_block(
                    realified_matrix,
                    reflected.row_basis,
                    reflected.row_basis;
                    atol,
                ),
                MP_R,
            )
            push!(
                real_blocks,
                (
                    spin2=block.spin2,
                    multiplicity=size(adapted_matrix, 1),
                    irrep_dimension=block.irrep_dimension,
                    reference_m2=block.reference_m2,
                    row_labels=reflected.row_labels,
                    matrix=adapted_matrix,
                    transform=_pauli_su2_translation_real_reflection_transform(
                        block,
                        reflected.row_basis,
                        k,
                        reflected.reflection_parity,
                    ),
                    coefficient_domain=:complex_algebraic_float64,
                    exact_coefficient_domain=:complex_sqrt_rational,
                    label=(
                        feature=:moment_matrix,
                        decomposition=:translation_su2_reflection,
                        momentum=k,
                        reflection_fixed=false,
                        spin2=block.spin2,
                        reflection_parity=reflected.reflection_parity,
                    ),
                ),
            )
        end
    end
    return real_blocks
end

function _pauli_su2_translation_reflection_fixed_linear_blocks(
    blocks,
    sector_basis::Vector{M},
    k::Int,
    n::Int;
    atol::Float64=1.0e-8,
) where {M<:NormalMonomial{PauliAlgebra}}
    reflection_matrix = _pauli_su2_translation_reflection_matrix(
        sector_basis,
        k,
        n;
        require_fixed=true,
    )
    adapted_blocks = Any[]
    for block in blocks
        transform_block = block.transform
        ref_rows = _pauli_su2_basis_m2_rows(transform_block, block.reference_m2)
        reference_transform = transform_block.transform[ref_rows, :]
        reduced_reflection =
            reference_transform * reflection_matrix * adjoint(reference_transform)
        hermitian_reflection =
            Hermitian((reduced_reflection + adjoint(reduced_reflection)) / 2)
        reflection_residual = max(
            _max_abs_entry(reduced_reflection - adjoint(reduced_reflection)),
            _max_abs_entry(reduced_reflection * reduced_reflection - I),
        )
        reflection_residual <= 100 * atol || throw(ArgumentError(
            "Translation/SU(2) reflection action for momentum $k, " *
            "spin2=$(block.spin2) is not a stable Hermitian involution; " *
            "residual=$reflection_residual."
        ))

        eig = LinearAlgebra.eigen(hermitian_reflection)
        indexed_matrix, index_to_key = _indexed_linear_block_entries(block.matrix)
        for (parity, sign) in ((:plus, 1.0), (:minus, -1.0))
            cols = findall(value -> abs(value - sign) <= 100 * atol, eig.values)
            isempty(cols) && continue
            vectors = eig.vectors[:, cols]
            row_basis = Matrix{ComplexF64}(adjoint(vectors))
            row_basis_sparse = _pauli_sparse_transform_rows(row_basis; atol)
            indexed_adapted_matrix = _transform_linear_block(
                indexed_matrix,
                row_basis_sparse,
                row_basis_sparse;
                atol,
            )
            adapted_matrix = _rekey_indexed_linear_block_entries(
                indexed_adapted_matrix,
                index_to_key,
            )
            row_labels = Any[
                _pauli_su2_reflection_row_label(
                    parity,
                    block.row_labels,
                    row_basis[row_idx, :];
                    atol,
                )
                for row_idx in axes(row_basis, 1)
            ]
            push!(
                adapted_blocks,
                (
                    spin2=block.spin2,
                    multiplicity=size(adapted_matrix, 1),
                    irrep_dimension=block.irrep_dimension,
                    reference_m2=block.reference_m2,
                    row_labels=row_labels,
                    matrix=adapted_matrix,
                    transform=_pauli_su2_translation_reflection_transform(
                        block,
                        row_basis,
                        k,
                        parity,
                    ),
                    coefficient_domain=:complex_algebraic_float64,
                    exact_coefficient_domain=:complex_sqrt_rational,
                    label=(
                        feature=:moment_matrix,
                        decomposition=:translation_su2_reflection,
                        momentum=k,
                        reflection_fixed=true,
                        spin2=block.spin2,
                        reflection_parity=parity,
                    ),
                ),
            )
        end
    end
    return adapted_blocks
end

function _pauli_su2_translation_real_reflection_linear_blocks_from_reduced(
    complex_blocks,
    sector_basis::Vector{M},
    k::Int,
    n::Int,
    ::Type{R};
    phase_atol::Real=1.0e-8,
) where {R<:Real,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    reflection_matrix = _pauli_su2_translation_reflection_matrix(sector_basis, k, n)
    atol = Float64(phase_atol)
    real_blocks = Any[]
    for block in complex_blocks
        realified_matrix = _symmetrize_real_linear_block(
            _realify_hermitian_linear_block_full(block.matrix, R; atol=R(phase_atol)),
        )
        indexed_realified_matrix, index_to_key =
            _indexed_linear_block_entries(realified_matrix)
        for reflected in _pauli_su2_translation_real_reflection_row_bases(
            block,
            reflection_matrix;
            atol,
        )
            row_basis_sparse = _pauli_sparse_transform_rows(reflected.row_basis; atol)
            indexed_adapted_matrix = _transform_real_symmetric_linear_block(
                indexed_realified_matrix,
                row_basis_sparse;
                atol,
            )
            adapted_matrix = _rekey_indexed_linear_block_entries(
                indexed_adapted_matrix,
                index_to_key,
            )
            push!(
                real_blocks,
                (
                    spin2=block.spin2,
                    multiplicity=size(adapted_matrix, 1),
                    irrep_dimension=block.irrep_dimension,
                    reference_m2=block.reference_m2,
                    row_labels=reflected.row_labels,
                    matrix=adapted_matrix,
                    transform=_pauli_su2_translation_real_reflection_transform(
                        block,
                        reflected.row_basis,
                        k,
                        reflected.reflection_parity,
                    ),
                    coefficient_domain=:complex_algebraic_float64,
                    exact_coefficient_domain=:complex_sqrt_rational,
                    label=(
                        feature=:moment_matrix,
                        decomposition=:translation_su2_reflection,
                        momentum=k,
                        reflection_fixed=false,
                        spin2=block.spin2,
                        reflection_parity=reflected.reflection_parity,
                    ),
                ),
            )
        end
    end
    return real_blocks
end

function _pauli_su2_translation_orbit_momentum_block_meta(
    ::Type{M},
    cone::Symbol,
    block,
) where {M<:NormalMonomial{PauliAlgebra}}
    return BlockMeta{M}(
        cone,
        PauliSU2BasisMomentOrigin(
            block.label,
            block.row_labels;
            transform=block.transform,
        ),
        fill(one(M), size(block.matrix, 1)),
    )
end

function _append_translation_complex_zero_constraints!(
    constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    zero_origin_by_constraint::Dict{Int,Any},
    source_constraints::Vector{Tuple{Symbol,Matrix{MP}}},
    source_zero_origin_by_constraint::Dict{Int,Any},
) where {MP<:Polynomial}
    for (source_idx, (cone, mat)) in pairs(source_constraints)
        cone == :Zero || throw(ArgumentError(
            "Expected Wigner-Eckart zero constraints only; source constraint " *
            "$source_idx has cone $cone."
        ))
        before = length(constraints)
        _append_constraint!(constraints, :Zero, mat, MP)
        source_seed = get(source_zero_origin_by_constraint, Int(source_idx), nothing)
        if source_seed !== nothing
            zero_origin_by_constraint[before + 1] = source_seed
        end
    end
    return nothing
end

function _append_realified_translation_zero_constraints!(
    constraints::Vector{Tuple{Symbol,Matrix{PReal}}},
    zero_origin_by_constraint::Dict{Int,Any},
    source_constraints::Vector{Tuple{Symbol,Matrix{PComplex}}},
    source_zero_origin_by_constraint::Dict{Int,Any},
    ::Type{PReal};
    phase_atol::Real,
    fallback_label,
) where {
    T<:Unsigned,
    PReal<:Polynomial{PauliAlgebra,T},
    PComplex<:Polynomial{PauliAlgebra,T},
}
    for (source_idx, (cone, mat)) in pairs(source_constraints)
        cone == :Zero || throw(ArgumentError(
            "Expected Wigner-Eckart zero constraints only; source constraint " *
            "$source_idx has cone $cone."
        ))
        source_seed = get(source_zero_origin_by_constraint, Int(source_idx), nothing)
        source_label = source_seed === nothing ?
            merge(fallback_label, (source_constraint=Int(source_idx),)) :
            _translation_zero_seed_label(source_seed)
        origin_seed_factory = source_seed === nothing ?
            TranslationZeroOriginSeed :
            _translation_zero_seed_factory(source_seed)
        for col in axes(mat, 2), row in axes(mat, 1)
            _append_translation_zero_polynomial!(
                constraints,
                mat[row, col],
                PReal;
                phase_atol,
                zero_origin_by_constraint,
                label=source_label,
                origin_seed_factory,
            )
        end
    end
    return nothing
end

function _realify_pauli_su2_translation_orbit_moment_problem(
    complex_mp::MomentProblem{PauliAlgebra,T,M,PComplex},
    ::Type{PReal};
    phase_atol::Real,
) where {
    T<:Unsigned,
    M<:NormalMonomial{PauliAlgebra,T},
    PComplex<:Polynomial{PauliAlgebra,T},
    PReal<:Polynomial{PauliAlgebra,T},
}
    constraints = Tuple{Symbol,Matrix{PReal}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()
    source_zero_origin_by_constraint = Dict{Int,Any}()
    _seed_translation_zero_origins_from_linear!(
        source_zero_origin_by_constraint,
        complex_mp.constraints,
        complex_mp.linear,
        PauliAlgebra,
    )
    psd_idx = 0
    for (constraint_idx, (cone, mat)) in pairs(complex_mp.constraints)
        if cone == :HPSD
            psd_idx += 1
            block = complex_mp.linear.psd_blocks_lin[psd_idx]
            real_mat = _symmetrize_real_polynomial_block(
                _realify_hermitian_block_full(mat, PReal; atol=phase_atol),
                PReal;
                atol=phase_atol,
            )
            push!(constraints, (:PSD, real_mat))
            block_meta_by_constraint[length(constraints)] = BlockMeta{M}(
                :PSD,
                block.meta.origin,
                fill(one(M), size(real_mat, 1)),
            )
        elseif cone == :Zero
            _append_realified_translation_zero_constraints!(
                constraints,
                zero_origin_by_constraint,
                Tuple{Symbol,Matrix{PComplex}}[(cone, mat)],
                Dict{Int,Any}(
                    1 => get(
                        source_zero_origin_by_constraint,
                        Int(constraint_idx),
                        TranslationZeroOriginSeed((
                            feature=:pauli_su2_base_zero,
                            decomposition=:translation_su2,
                            reason=:realified_zero,
                            source_constraint=constraint_idx,
                        )),
                    ),
                ),
                PReal;
                phase_atol,
                fallback_label=(
                    feature=:pauli_su2_base_zero,
                    decomposition=:translation_su2,
                    reason=:realified_zero,
                ),
            )
        else
            throw(ArgumentError(
                "Translation/SU(2) realification expects HPSD or Zero constraints; " *
                "constraint $constraint_idx has cone $cone."
            ))
        end
    end
    psd_idx == length(complex_mp.linear.psd_blocks_lin) || throw(ArgumentError(
        "Translation/SU(2) realification saw $psd_idx constraints but " *
        "$(length(complex_mp.linear.psd_blocks_lin)) linear PSD blocks."
    ))

    objective = _real_part_polynomial(complex_mp.objective, PReal; atol=phase_atol)
    total_basis = _collect_pauli_su2_moment_basis(objective, constraints)
    return MomentProblem{PauliAlgebra,T,M,PReal}(
        objective,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=true,
    )
end

function _pauli_su2_translation_orbit_reflected_real_moment_problem(
    objective::Polynomial{PauliAlgebra,T,C},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
    ::Type{PReal};
    assume_su2_invariant::Bool=false,
    ops=nothing,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    phase_atol::Real=1.0e-8,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {
    T<:Unsigned,
    C<:Number,
    PReal<:Polynomial{PauliAlgebra,T},
}
    verified_su2_invariant = assume_su2_invariant
    if ops !== nothing
        assume_su2_invariant || _check_pauli_axis_rotation_invariance(
            objective,
            ops;
            context="translation-orbit SU(2) reflected objective",
        )
        verified_su2_invariant = true
    end
    verified_su2_invariant || throw(ArgumentError(
        "`pauli_su2_translation_orbit_moment_problem` drops SU(2) " *
        "off-sector and magnetic-copy rows; pass `ops=(σx, σy, σz)` to " *
        "verify global Pauli-axis-rotation invariance or " *
        "`assume_su2_invariant=true` after verifying it elsewhere."
    ))

    n = Int(n_sites)
    _check_translation_invariance(
        objective,
        n;
        context="translation-orbit SU(2) reflected objective",
    )

    R = _coefficient_type(PReal)
    atol = R(phase_atol)
    M = NormalMonomial{PauliAlgebra,T}
    MP_C = Polynomial{PauliAlgebra,T,ComplexF64}
    objective_reduced = _translation_orbit_reduce_polynomial(objective, n)
    objective_mp = _real_part_polynomial(convert(MP_C, objective_reduced), PReal; atol)
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix=true)

    constraints = Tuple{Symbol,Matrix{PReal}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()
    wigner_zero_constraints = Tuple{Symbol,Matrix{MP_C}}[]
    wigner_zero_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    wigner_zero_origin_by_constraint = Dict{Int,Any}()
    basis_vec = M[mono for mono in basis]
    for k in sectors
        sector_reference_m2 =
            _translation_reflection_fixed_momentum(k, n) ? reference_m2 : :zero
        sector_basis, reduced_blocks = _pauli_su2_translation_wigner_sector_blocks!(
            wigner_zero_constraints,
            wigner_zero_meta_by_constraint,
            wigner_zero_origin_by_constraint,
            basis_vec,
            n,
            k,
            MP_C,
            M;
            max_dimension,
            reference_m2=sector_reference_m2,
            phase_atol,
            append_psd=false,
            singlet_channel_equalities,
            singlet_channel_atol,
        )
        sector_blocks = if _translation_reflection_fixed_momentum(k, n)
            _pauli_su2_translation_reflection_fixed_momentum_blocks(
                reduced_blocks,
                sector_basis,
                k,
                n,
            )
        else
            _pauli_su2_translation_real_reflection_momentum_blocks_from_reduced(
                reduced_blocks,
                sector_basis,
                k,
                n;
                phase_atol=Float64(phase_atol),
            )
        end

        for block in sector_blocks
            mat = if _translation_reflection_fixed_momentum(k, n)
                _symmetrize_real_polynomial_block(
                    _realify_hermitian_block_full(block.matrix, PReal; atol),
                    PReal;
                    atol,
                )
            else
                map(entry -> convert(PReal, entry), block.matrix)
            end
            push!(constraints, (:PSD, mat))
            block_meta_by_constraint[length(constraints)] = BlockMeta{M}(
                :PSD,
                PauliSU2BasisMomentOrigin(
                    block.label,
                    block.row_labels;
                    transform=block.transform,
                ),
                fill(one(M), size(mat, 1)),
            )
        end
    end
    _append_realified_translation_zero_constraints!(
        constraints,
        zero_origin_by_constraint,
        wigner_zero_constraints,
        wigner_zero_origin_by_constraint,
        PReal;
        phase_atol=atol,
        fallback_label=(
            feature=:pauli_su2_base_zero,
            decomposition=:translation_su2_reflection,
            reason=:wigner_eckart,
        ),
    )

    constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
    _check_polynomial_moments_covered(
        objective_mp,
        constraint_basis,
        "translation-orbit SU(2) reflected objective",
    )
    total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
    return MomentProblem{PauliAlgebra,T,M,PReal}(
        objective_mp,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=true,
    )
end

"""
    pauli_su2_translation_orbit_moment_problem(objective, basis, n_sites; assume_su2_invariant=false, ops=nothing, real_moment_matrix=true, momenta=nothing, reflection_symmetry=false, max_dimension=4096, reference_m2=:lowest, phase_atol=1e-8, singlet_channel_equalities=false, singlet_channel_atol=1e-12)

Build a low-level Pauli `MomentProblem` by constructing the canonical
translation momentum blocks on a translation-orbit basis and reducing each one
by the SU(2) Clebsch-Gordan transform.  The emitted constraints are complex
HPSD blocks with translation/SU(2) momentum-spin labels; this helper still
does not perform sign reduction, direct JuMP emission, or an SDP solve.
Reflection reduction is supported only for reflection-fixed requested momentum
sectors.
"""
function pauli_su2_translation_orbit_moment_problem(
    objective::Polynomial{PauliAlgebra,T,C},
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    assume_su2_invariant::Bool=false,
    ops=nothing,
    real_moment_matrix::Bool=true,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    reflection_symmetry::Bool=false,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    phase_atol::Real=1.0e-8,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
) where {T<:Integer,C<:Number}
    verified_su2_invariant = assume_su2_invariant
    if ops !== nothing
        assume_su2_invariant || _check_pauli_axis_rotation_invariance(
            objective,
            ops;
            context="translation-orbit SU(2) objective",
        )
        verified_su2_invariant = true
    end
    verified_su2_invariant || throw(ArgumentError(
        "`pauli_su2_translation_orbit_moment_problem` drops SU(2) " *
        "off-sector and magnetic-copy rows; pass `ops=(σx, σy, σz)` to " *
        "verify global Pauli-axis-rotation invariance or " *
        "`assume_su2_invariant=true` after verifying it elsewhere."
    ))
    if T <: Unsigned
        _check_translation_invariance(
            objective,
            n_sites;
            context="translation-orbit SU(2) objective",
        )
        if reflection_symmetry && real_moment_matrix
            n = Int(n_sites)
            sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix=true)
            if any(k -> !_translation_reflection_fixed_momentum(k, n), sectors)
                MP_R = Polynomial{PauliAlgebra,T,Float64}
                return _pauli_su2_translation_orbit_reflected_real_moment_problem(
                    objective,
                    basis,
                    n,
                    MP_R;
                    assume_su2_invariant,
                    ops,
                    momenta,
                    max_dimension,
                    reference_m2,
                    phase_atol,
                    singlet_channel_equalities,
                    singlet_channel_atol,
                )
            end
        end
    end

    M = NormalMonomial{PauliAlgebra,T}
    MP_P = Polynomial{PauliAlgebra,T,ComplexF64}
    objective_reduced = T <: Unsigned ?
        _translation_orbit_reduce_polynomial(objective, n_sites) :
        objective
    objective_mp = convert(MP_P, objective_reduced)
    if reflection_symmetry && T <: Unsigned
        n = Int(n_sites)
        sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
        nonfixed = Int[k for k in sectors if !_translation_reflection_fixed_momentum(k, n)]
        isempty(nonfixed) || throw(ArgumentError(
            "Translation/SU(2) complex reflection reduction requires " *
            "reflection-fixed momentum sectors; got non-fixed sector(s) $nonfixed."
        ))

        constraints = Tuple{Symbol,Matrix{MP_P}}[]
        block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
        zero_origin_by_constraint = Dict{Int,Any}()
        wigner_zero_constraints = Tuple{Symbol,Matrix{MP_P}}[]
        wigner_zero_meta_by_constraint = Dict{Int,BlockMeta{M}}()
        wigner_zero_origin_by_constraint = Dict{Int,Any}()
        basis_vec = M[mono for mono in basis]
        for k in sectors
            sector_basis, reduced_blocks = _pauli_su2_translation_wigner_sector_blocks!(
                wigner_zero_constraints,
                wigner_zero_meta_by_constraint,
                wigner_zero_origin_by_constraint,
                basis_vec,
                n,
                k,
                MP_P,
                M;
                max_dimension,
                reference_m2,
                phase_atol,
                append_psd=false,
                singlet_channel_equalities,
                singlet_channel_atol,
            )
            for block in _pauli_su2_translation_reflection_fixed_momentum_blocks(
                reduced_blocks,
                sector_basis,
                k,
                n,
            )
                before = length(constraints)
                _append_constraint!(constraints, :HPSD, block.matrix, MP_P)
                block_meta_by_constraint[before + 1] =
                    _pauli_su2_translation_orbit_momentum_block_meta(M, :HPSD, block)
            end
        end
        _append_translation_complex_zero_constraints!(
            constraints,
            zero_origin_by_constraint,
            wigner_zero_constraints,
            wigner_zero_origin_by_constraint,
        )

        constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
        _check_polynomial_moments_covered(
            objective_mp,
            constraint_basis,
            "translation-orbit SU(2) objective",
        )
        total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
        return MomentProblem{PauliAlgebra,T,M,MP_P}(
            objective_mp,
            constraints,
            total_basis,
            length(total_basis);
            block_meta_by_constraint,
            zero_origin_by_constraint,
            real_moments=false,
        )
    end

    bundle = pauli_su2_translation_orbit_momentum_block_bundle(
        basis,
        n_sites;
        real_moment_matrix,
        momenta,
        reflection_symmetry,
        max_dimension,
        reference_m2,
    )

    constraints = Tuple{Symbol,Matrix{MP_P}}[]
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    zero_origin_by_constraint = Dict{Int,Any}()
    for block in bundle.blocks
        before = length(constraints)
        _append_constraint!(constraints, :HPSD, block.matrix, MP_P)
        block_meta_by_constraint[before + 1] =
            _pauli_su2_translation_orbit_momentum_block_meta(M, :HPSD, block)
    end
    if singlet_channel_equalities
        sectors = _pauli_chain_momentum_sectors(n_sites, momenta; real_moment_matrix)
        basis_vec = M[mono for mono in basis]
        for k in sectors
            sector_basis = iszero(k) ? basis_vec : M[mono for mono in basis_vec if !isone(mono)]
            _append_pauli_su2_translation_singlet_channel_equality_constraints!(
                constraints,
                zero_origin_by_constraint,
                sector_basis,
                n_sites,
                k,
                MP_P;
                max_dimension,
                atol=singlet_channel_atol,
            )
        end
    end

    constraint_basis = _collect_pauli_su2_moment_basis(zero(objective_mp), constraints)
    _check_polynomial_moments_covered(
        objective_mp,
        constraint_basis,
        "translation-orbit SU(2) objective",
    )
    total_basis = sorted_unique!(vcat(monomials(objective_mp), constraint_basis))
    return MomentProblem{PauliAlgebra,T,M,MP_P}(
        objective_mp,
        constraints,
        total_basis,
        length(total_basis);
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=false,
    )
end

"""
    pauli_su2_basis_metrics(basis; scalar_bytes=sizeof(Float64))

Return storage-size estimates for a support-complete Pauli word basis before and
after an SU(2) irreducible-tensor moment-matrix split.  The reduced estimates
count one multiplicity PSD block per total-spin sector.
"""
function pauli_su2_basis_metrics(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    scalar_bytes::Integer=sizeof(Float64),
)
    summary = pauli_su2_basis_summary(basis)
    return (
        support_counts=summary.support_counts,
        blocks=summary.blocks,
        _pauli_su2_storage_metrics(summary.blocks, summary.word_count, scalar_bytes)...,
    )
end

function _pauli_su2_blocks_from_support_counts(support_counts)
    spin_multiplicities = Dict{Int,Int}()
    for (support_size, support_count) in support_counts
        for block in pauli_su2_word_blocks(support_size)
            spin_multiplicities[block.spin2] =
                get(spin_multiplicities, block.spin2, 0) +
                Int(support_count) * block.multiplicity
        end
    end
    return [
        PauliSU2WordBlock(spin2, spin_multiplicities[spin2])
        for spin2 in sort!(collect(keys(spin_multiplicities)))
    ]
end

function _pauli_su2_word_singlet_channel_count(support_size::Integer)
    for block in pauli_su2_word_blocks(support_size)
        block.spin2 == 0 && return block.multiplicity
    end
    return 0
end

function _pauli_su2_singlet_channel_support_counts(support_counts)
    channel_counts = Pair{Int,Int}[]
    for (support_size, support_count) in support_counts
        channel_count =
            Int(support_count) * _pauli_su2_word_singlet_channel_count(support_size)
        channel_count > 0 && push!(channel_counts, Int(support_size) => channel_count)
    end
    return channel_counts
end

function _pauli_su2_translation_reflection_fixed_blocks(
    order::Int;
    include_identity::Bool,
    flip_odd_supports::Bool,
    max_dimension::Int,
)
    plus_counts = Dict{Int,Int}()
    minus_counts = Dict{Int,Int}()
    include_identity && (plus_counts[0] = 1)

    for support_size in 1:order
        diag = pauli_su2_word_reflection_diagnostics(
            support_size;
            max_dimension,
        )
        flip_parity = flip_odd_supports && isodd(support_size)
        for parity in diag.spin_parities
            plus = flip_parity ? parity.minus_multiplicity : parity.plus_multiplicity
            minus = flip_parity ? parity.plus_multiplicity : parity.minus_multiplicity
            plus_counts[parity.spin2] = get(plus_counts, parity.spin2, 0) + plus
            minus_counts[parity.spin2] = get(minus_counts, parity.spin2, 0) + minus
        end
    end

    blocks = NamedTuple[]
    for spin2 in sort!(unique!(vcat(collect(keys(plus_counts)), collect(keys(minus_counts)))))
        plus = get(plus_counts, spin2, 0)
        minus = get(minus_counts, spin2, 0)
        plus > 0 && push!(
            blocks,
            (size=plus, spin2=spin2, reflection_parity=:plus),
        )
        minus > 0 && push!(
            blocks,
            (size=minus, spin2=spin2, reflection_parity=:minus),
        )
    end
    return blocks
end

function _pauli_su2_translation_reflection_fixed_block_sizes(
    order::Int;
    include_identity::Bool,
    flip_odd_supports::Bool,
    max_dimension::Int,
)
    blocks = _pauli_su2_translation_reflection_fixed_blocks(
        order;
        include_identity,
        flip_odd_supports,
        max_dimension,
    )
    sizes = Int[block.size for block in blocks]
    return sizes
end

function _pauli_su2_word_block_accounting(blocks)
    full_side = sum(block.irrep_dimension * block.multiplicity for block in blocks; init=0)
    active_dense_entries = sum(
        block.irrep_dimension * block.multiplicity^2 for block in blocks; init=0
    )
    reduced_dense_entries = sum(block.multiplicity^2 for block in blocks; init=0)
    return (
        full_dense_entries=full_side^2,
        active_dense_entries=active_dense_entries,
        reduced_dense_entries=reduced_dense_entries,
    )
end

function _add_pauli_su2_accounting(left, right)
    return (
        full_dense_entries=left.full_dense_entries + right.full_dense_entries,
        active_dense_entries=left.active_dense_entries + right.active_dense_entries,
        reduced_dense_entries=left.reduced_dense_entries + right.reduced_dense_entries,
    )
end

function _pauli_su2_reflection_fixed_accounting(blocks)
    total = (full_dense_entries=0, active_dense_entries=0, reduced_dense_entries=0)
    for reflection_parity in (:plus, :minus)
        parity_blocks = [block for block in blocks if block.reflection_parity == reflection_parity]
        isempty(parity_blocks) && continue
        full_side = sum((block.spin2 + 1) * block.size for block in parity_blocks; init=0)
        active_dense_entries = sum(
            (block.spin2 + 1) * block.size^2 for block in parity_blocks; init=0
        )
        reduced_dense_entries = sum(block.size^2 for block in parity_blocks; init=0)
        total = _add_pauli_su2_accounting(
            total,
            (
                full_dense_entries=full_side^2,
                active_dense_entries=active_dense_entries,
                reduced_dense_entries=reduced_dense_entries,
            ),
        )
    end
    return total
end

function _pauli_su2_accounting_summary(accounting)
    offblock_entry_count = accounting.full_dense_entries - accounting.active_dense_entries
    copy_entry_count = accounting.active_dense_entries - accounting.reduced_dense_entries
    return (
        su2_full_dense_entries=accounting.full_dense_entries,
        su2_active_dense_entries=accounting.active_dense_entries,
        su2_reduced_dense_entries=accounting.reduced_dense_entries,
        offblock_entry_count=offblock_entry_count,
        copy_entry_count=copy_entry_count,
        accounted_entry_count=offblock_entry_count + copy_entry_count +
            accounting.reduced_dense_entries,
    )
end

function _pauli_su2_accounting_record(label, accounting)
    return (label=label, _pauli_su2_accounting_summary(accounting)...)
end

function _pauli_su2_reflection_fixed_accounting_records(momentum::Int, blocks)
    records = Any[]
    for reflection_parity in (:plus, :minus)
        parity_blocks = [block for block in blocks if block.reflection_parity == reflection_parity]
        isempty(parity_blocks) && continue
        full_side = sum((block.spin2 + 1) * block.size for block in parity_blocks; init=0)
        accounting = (
            full_dense_entries=full_side^2,
            active_dense_entries=sum(
                (block.spin2 + 1) * block.size^2 for block in parity_blocks; init=0
            ),
            reduced_dense_entries=sum(block.size^2 for block in parity_blocks; init=0),
        )
        push!(
            records,
            _pauli_su2_accounting_record(
                (
                    feature=:moment_matrix,
                    decomposition=:translation_su2_reflection,
                    momentum=momentum,
                    reflection_fixed=true,
                    reflection_parity=reflection_parity,
                ),
                accounting,
            ),
        )
    end
    return records
end

"""
    pauli_su2_contiguous_chain_structural_targets(n_sites, order; scalar_bytes=sizeof(Float64))

Return analytic SU(2) irreducible-tensor structural targets for the periodic
contiguous Pauli-chain basis without constructing the site-space basis.  The
formula is valid in the no-short-orbit regime `n_sites > 2order`, where each
active support width `1:order` appears once per site and contains all `3^width`
non-identity Pauli words.  These are full-basis SU(2) targets, not the future
translation/reflection-combined block sizes.
"""
function pauli_su2_contiguous_chain_structural_targets(
    n_sites::Integer,
    order::Integer;
    scalar_bytes::Integer=sizeof(Float64),
)
    n = Int(n_sites)
    d = Int(order)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    2d < n || throw(ArgumentError(
        "Analytic Pauli SU(2) contiguous-chain targets require `n_sites > 2order`; got n_sites=$n, order=$d."
    ))

    support_counts = Pair{Int,Int}[0 => 1]
    for support_size in 1:d
        push!(support_counts, support_size => n)
    end

    blocks = _pauli_su2_blocks_from_support_counts(support_counts)
    basis_size = 1 + n * sum(3^support_size for support_size in 1:d; init=0)
    metrics = (
        support_counts=support_counts,
        blocks=blocks,
        _pauli_su2_storage_metrics(blocks, basis_size, Int(scalar_bytes))...,
    )
    singlet_channel_support_counts =
        _pauli_su2_singlet_channel_support_counts(support_counts)

    return (
        n_sites=n,
        order=d,
        basis_size=basis_size,
        word_count=basis_size,
        metrics...,
        singlet_channel_count=sum(last, singlet_channel_support_counts; init=0),
        singlet_channel_support_counts=singlet_channel_support_counts,
        singlet_channel_equality_row_count=basis_size -
            sum(last, singlet_channel_support_counts; init=0),
        _pauli_su2_reduction_accounting(metrics)...,
        translation_combined=false,
        solve_supported=false,
        solve_blocker=:structural_target_only,
        solve_blocker_reason=_TRANSLATION_SOLVE_STRUCTURAL_TARGET_REASON,
        estimated_model_size_gate_status=:blocked_missing_scalar_equality_estimate,
        estimated_model_size_gate_reason=_TRANSLATION_STRUCTURAL_MODEL_SIZE_GATE_REASON,
        solve_unsupported_block_features=Symbol[],
        solve_unsupported_zero_features=Symbol[],
        requires_construction=false,
        assumptions=(
            periodic_contiguous_basis=true,
            support_complete=true,
            no_short_orbits=true,
            n_sites_gt_2order=true,
            translation_combined=false,
        ),
    )
end

"""
    pauli_su2_translation_orbit_structural_targets(n_sites, order; real_moment_matrix=true, momenta=nothing, reflection_symmetry=false, scalar_bytes=sizeof(Float64), max_dimension=4096)

Return analytic SU(2) structural targets after translation has reduced a
periodic contiguous Pauli-chain basis to one representative per active support
width and Pauli word.  The target is valid for `n_sites > 2order`.  It is a
shape target for the future translation/SU(2) reducer: it does not include
solver construction or any numerical SDP solve.  With
`reflection_symmetry=true`, the target uses the fixed-support reflection
diagnostics to split SU(2) multiplicity spaces by reversal parity.
"""
function pauli_su2_translation_orbit_structural_targets(
    n_sites::Integer,
    order::Integer;
    real_moment_matrix::Bool=true,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    reflection_symmetry::Bool=false,
    scalar_bytes::Integer=sizeof(Float64),
    max_dimension::Integer=4096,
)
    n = Int(n_sites)
    d = Int(order)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    2d < n || throw(ArgumentError(
        "Analytic Pauli translation/SU(2) orbit targets require `n_sites > 2order`; got n_sites=$n, order=$d."
    ))
    bytes = Int(scalar_bytes)
    bytes > 0 || throw(ArgumentError("`scalar_bytes` must be positive; got $scalar_bytes."))

    orbit_support_counts = Pair{Int,Int}[0 => 1]
    nonidentity_support_counts = Pair{Int,Int}[]
    for support_size in 1:d
        push!(orbit_support_counts, support_size => 1)
        push!(nonidentity_support_counts, support_size => 1)
    end
    zero_momentum_blocks = _pauli_su2_blocks_from_support_counts(orbit_support_counts)
    nonzero_momentum_blocks = isempty(nonidentity_support_counts) ?
        PauliSU2WordBlock[] :
        _pauli_su2_blocks_from_support_counts(nonidentity_support_counts)

    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    real_moment_matrix || 0 in sectors || throw(ArgumentError(
        "Momentum sector 0 is required because it carries the normalized identity moment."
    ))
    if reflection_symmetry && !real_moment_matrix
        nonfixed = Int[k for k in sectors if !_translation_reflection_fixed_momentum(k, n)]
        isempty(nonfixed) || throw(ArgumentError(
            "`reflection_symmetry=true, real_moment_matrix=false` SU(2) orbit targets " *
            "are supported only for reflection-fixed momentum sectors; got non-fixed " *
            "sector(s) $nonfixed."
        ))
    end

    logical_sizes = Int[]
    psd_sizes = Int[]
    block_labels = Any[]
    block_coefficient_domains = Symbol[]
    block_exact_coefficient_domains = Symbol[]
    accounting = (full_dense_entries=0, active_dense_entries=0, reduced_dense_entries=0)
    accounting_records = Any[]
    if reflection_symmetry
        for k in sectors
            if _translation_reflection_fixed_momentum(k, n)
                fixed_blocks = _pauli_su2_translation_reflection_fixed_blocks(
                    d;
                    include_identity=iszero(k),
                    flip_odd_supports=!iszero(k),
                    max_dimension=Int(max_dimension),
                )
                fixed_sizes = Int[block.size for block in fixed_blocks]
                append!(logical_sizes, fixed_sizes)
                append!(psd_sizes, real_moment_matrix ? 2 .* fixed_sizes : fixed_sizes)
                accounting = _add_pauli_su2_accounting(
                    accounting,
                    _pauli_su2_reflection_fixed_accounting(fixed_blocks),
                )
                append!(
                    accounting_records,
                    _pauli_su2_reflection_fixed_accounting_records(k, fixed_blocks),
                )
                for block in fixed_blocks
                    push!(
                        block_labels,
                        (
                            feature=:moment_matrix,
                            decomposition=:translation_su2_reflection,
                            momentum=k,
                            reflection_fixed=true,
                            spin2=block.spin2,
                            reflection_parity=block.reflection_parity,
                        ),
                    )
                    push!(block_coefficient_domains, :complex_algebraic_float64)
                    push!(block_exact_coefficient_domains, :complex_sqrt_rational)
                end
            else
                blocks = nonzero_momentum_blocks
                sector_accounting = _pauli_su2_word_block_accounting(blocks)
                for reflection_parity in (:plus, :minus)
                    accounting = _add_pauli_su2_accounting(accounting, sector_accounting)
                    push!(
                        accounting_records,
                        _pauli_su2_accounting_record(
                            (
                                feature=:moment_matrix,
                                decomposition=:translation_su2_reflection,
                                momentum=k,
                                reflection_fixed=false,
                                reflection_parity=reflection_parity,
                            ),
                            sector_accounting,
                        ),
                    )
                end
                for block in blocks
                    push!(logical_sizes, block.multiplicity)
                    push!(logical_sizes, block.multiplicity)
                    push!(psd_sizes, block.multiplicity)
                    push!(psd_sizes, block.multiplicity)
                    push!(
                        block_labels,
                        (
                            feature=:moment_matrix,
                            decomposition=:translation_su2_reflection,
                            momentum=k,
                            reflection_fixed=false,
                            spin2=block.spin2,
                            reflection_parity=:plus,
                        ),
                    )
                    push!(block_coefficient_domains, :complex_algebraic_float64)
                    push!(block_exact_coefficient_domains, :complex_sqrt_rational)
                    push!(
                        block_labels,
                        (
                            feature=:moment_matrix,
                            decomposition=:translation_su2_reflection,
                            momentum=k,
                            reflection_fixed=false,
                            spin2=block.spin2,
                            reflection_parity=:minus,
                        ),
                    )
                    push!(block_coefficient_domains, :complex_algebraic_float64)
                    push!(block_exact_coefficient_domains, :complex_sqrt_rational)
                end
            end
        end
    else
        for k in sectors
            blocks = k == 0 ? zero_momentum_blocks : nonzero_momentum_blocks
            sector_accounting = _pauli_su2_word_block_accounting(blocks)
            accounting = _add_pauli_su2_accounting(accounting, sector_accounting)
            push!(
                accounting_records,
                _pauli_su2_accounting_record(
                    (
                        feature=:moment_matrix,
                        decomposition=:translation_su2,
                        momentum=k,
                    ),
                    sector_accounting,
                ),
            )
            append!(logical_sizes, (block.multiplicity for block in blocks))
            for block in blocks
                push!(
                    block_labels,
                    (
                        feature=:moment_matrix,
                        decomposition=:translation_su2,
                        momentum=k,
                        spin2=block.spin2,
                    ),
                )
                push!(block_coefficient_domains, :complex_algebraic_float64)
                push!(block_exact_coefficient_domains, :complex_sqrt_rational)
            end
        end
        psd_sizes = real_moment_matrix ? [2 * size for size in logical_sizes] : copy(logical_sizes)
    end
    psd_dense_entries = sum(size * size for size in psd_sizes; init=0)
    psd_symmetric_entries = sum(
        real_moment_matrix ? size * (size + 1) ÷ 2 : size * size for size in psd_sizes;
        init=0,
    )
    su2_accounting = _pauli_su2_accounting_summary(accounting)
    wigner_offblock_zero_row_budget = 2 * su2_accounting.offblock_entry_count
    wigner_magnetic_copy_zero_row_budget = su2_accounting.copy_entry_count
    wigner_zero_row_budget =
        wigner_offblock_zero_row_budget + wigner_magnetic_copy_zero_row_budget
    singlet_channel_support_counts =
        _pauli_su2_singlet_channel_support_counts(orbit_support_counts)
    orbit_basis_size = 1 + sum(3^support_size for support_size in 1:d; init=0)
    singlet_channel_count = sum(last, singlet_channel_support_counts; init=0)

    return (
        n_sites=n,
        order=d,
        orbit_basis_size=orbit_basis_size,
        support_counts=orbit_support_counts,
        singlet_channel_count=singlet_channel_count,
        singlet_channel_support_counts=singlet_channel_support_counts,
        singlet_channel_equality_row_count=orbit_basis_size - singlet_channel_count,
        zero_momentum_blocks=zero_momentum_blocks,
        nonzero_momentum_blocks=nonzero_momentum_blocks,
        momentum_sectors=sectors,
        momentum_sector_count=length(sectors),
        real_moment_matrix=real_moment_matrix,
        block_labels=block_labels,
        block_coefficient_domains=block_coefficient_domains,
        block_exact_coefficient_domains=block_exact_coefficient_domains,
        block_coefficient_domain_histogram=_value_histogram_pairs(block_coefficient_domains),
        block_exact_coefficient_domain_histogram=_value_histogram_pairs(block_exact_coefficient_domains),
        logical_block_sizes=logical_sizes,
        psd_block_sizes=psd_sizes,
        n_blocks=length(psd_sizes),
        logical_max_block=isempty(logical_sizes) ? 0 : maximum(logical_sizes),
        psd_max_block=isempty(psd_sizes) ? 0 : maximum(psd_sizes),
        logical_total_block_side=sum(logical_sizes; init=0),
        psd_total_block_side=sum(psd_sizes; init=0),
        logical_block_histogram=_histogram_pairs(logical_sizes),
        psd_block_histogram=_histogram_pairs(psd_sizes),
        logical_block_feature_histogram=_label_size_histogram(block_labels, logical_sizes),
        psd_block_feature_histogram=_label_size_histogram(block_labels, psd_sizes),
        psd_dense_entries=psd_dense_entries,
        psd_symmetric_entries=psd_symmetric_entries,
        psd_dense_bytes=psd_dense_entries * bytes,
        psd_symmetric_bytes=psd_symmetric_entries * bytes,
        su2_accounting...,
        wigner_offblock_zero_row_budget=wigner_offblock_zero_row_budget,
        wigner_magnetic_copy_zero_row_budget=wigner_magnetic_copy_zero_row_budget,
        wigner_zero_row_budget=wigner_zero_row_budget,
        su2_accounting_records=accounting_records,
        su2_accounting_record_count=length(accounting_records),
        translation_combined=true,
        reflection_combined=reflection_symmetry,
        sign_symmetry_subsumed=true,
        solve_supported=false,
        solve_blocker=:structural_target_only,
        solve_blocker_reason=_TRANSLATION_SOLVE_STRUCTURAL_TARGET_REASON,
        estimated_model_size_gate_status=:blocked_missing_scalar_equality_estimate,
        estimated_model_size_gate_reason=_TRANSLATION_STRUCTURAL_MODEL_SIZE_GATE_REASON,
        solve_unsupported_block_features=Symbol[],
        solve_unsupported_zero_features=Symbol[],
        requires_construction=false,
        assumptions=(
            periodic_contiguous_basis=true,
            support_complete=true,
            no_short_orbits=true,
            n_sites_gt_2order=true,
            translation_combined=true,
            reflection_combined=reflection_symmetry,
        ),
    )
end

"""
    pauli_su2_basis_reduction_diagnostics(basis; scalar_bytes=sizeof(Float64))

Return entry-count diagnostics for the SU(2) irreducible-tensor reduction of a
support-complete Pauli word basis.
"""
function pauli_su2_basis_reduction_diagnostics(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    scalar_bytes::Integer=sizeof(Float64),
)
    metrics = pauli_su2_basis_metrics(basis; scalar_bytes)
    return (; metrics..., _pauli_su2_reduction_accounting(metrics)...)
end

"""
    pauli_su2_basis_spin_diagnostics(basis; max_dimension=4096)

Return direct-sum SU(2) spin diagnostics for a support-complete Pauli word
basis.  Each active support is checked as a fixed-support spin-1 tensor product,
then multiplicities are accumulated across supports.
"""
function pauli_su2_basis_spin_diagnostics(
    basis::AbstractVector{<:NormalMonomial{PauliAlgebra}};
    max_dimension::Integer=4096,
)
    summary = pauli_su2_basis_summary(basis)
    spin_multiplicity_counts = Dict{Int,Int}()
    dimension = 0
    unitarity_residual = 0.0
    sz_residual = 0.0
    casimir_residual = 0.0
    offblock_residual = 0.0

    for (support_size, support_count) in summary.support_counts
        diag = pauli_su2_word_spin_diagnostics(support_size; max_dimension)
        dimension += support_count * diag.dimension
        unitarity_residual = max(unitarity_residual, diag.unitarity_residual)
        sz_residual = max(sz_residual, diag.sz_residual)
        casimir_residual = max(casimir_residual, diag.casimir_residual)
        offblock_residual = max(offblock_residual, diag.offblock_residual)
        for (spin2, multiplicity) in diag.spin_multiplicities
            spin_multiplicity_counts[spin2] =
                get(spin_multiplicity_counts, spin2, 0) +
                support_count * multiplicity
        end
    end
    dimension == summary.word_count || throw(ArgumentError(
        "Internal Pauli SU(2) basis diagnostic dimension mismatch: got $dimension, expected $(summary.word_count)."
    ))

    spin_multiplicities = [
        spin2 => spin_multiplicity_counts[spin2]
        for spin2 in sort!(collect(keys(spin_multiplicity_counts)))
    ]
    return (
        support_counts=summary.support_counts,
        dimension=dimension,
        state_count=dimension,
        spin_multiplicities=spin_multiplicities,
        unitarity_residual=unitarity_residual,
        sz_residual=sz_residual,
        casimir_residual=casimir_residual,
        offblock_residual=offblock_residual,
    )
end

"""
    TranslationInvariantReport

Diagnostic summary for [`pauli_translation_invariant_moment_relaxation`](@ref).
"""
struct TranslationInvariantReport
    n_sites::Int
    order::Int
    basis_size::Int
    orbit_basis_size::Int
    axis_orbit_closed::Bool
    axis_orbit_basis_size::Int
    axis_orbit_size_histogram::Vector{Pair{Int,Int}}
    axis_reduction_ratio::Float64
    momentum_sectors::Vector{Int}
    sign_symmetry::Bool
    logical_block_sizes::Vector{Int}
    psd_block_sizes::Vector{Int}
    block_labels::Vector{Any}
    block_coefficient_domains::Vector{Union{Nothing,Symbol}}
    block_exact_coefficient_domains::Vector{Union{Nothing,Symbol}}
    n_unique_moment_matrix_elements::Int
    linear_moment_count::Int
    zero_constraint_count::Int
    construction_time_ns::Int
    construction_stage_times_ns::Dict{Symbol,Int}
    real_moment_matrix::Bool
    product_cache_hits::Int
    product_cache_misses::Int
    product_cache_entries::Int
    product_cache_hit_rate::Float64
    zero_constraint_feature_histogram::Vector{Pair{Any,Int}}
    axis_rotation_quotient::Bool
    axis_rotation_raw_moment_key_count::Int
    axis_rotation_moment_class_count::Int
    axis_rotation_quotient_moment_key_count::Int
    axis_rotation_forced_zero_moment_class_count::Int
    axis_rotation_moment_quotient_reduction_ratio::Float64
    su2_moment_quotient::Bool
    su2_moment_raw_count::Int
    su2_moment_quotient_count::Int
    su2_moment_quotient_reduction_ratio::Float64
    su2_moment_support_orbit_count::Int
    su2_moment_singlet_channel_support_counts::Vector{Pair{Int,Int}}
    su2_moment_max_pivot_residual::Float64
    su2_moment_max_invariant_residual::Float64
    su2_moment_max_reconstruction_residual::Float64
    su2_moment_max_condition::Float64
    su2_moment_eliminated_zero_row_count::Int
    su2_moment_eliminated_zero_feature_histogram::Vector{Pair{Any,Int}}
end

function TranslationInvariantReport(args::Vararg{Any,32})
    linear_moment_count = Int(args[17])
    return TranslationInvariantReport(
        args...,
        false,
        linear_moment_count,
        linear_moment_count,
        1.0,
        0,
        Pair{Int,Int}[],
        0.0,
        0.0,
        0.0,
        1.0,
        0,
        Pair{Any,Int}[],
    )
end

struct TranslationDFTTransform
    family::Symbol
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    n_sites::Int
    momentum::Int
    real_moment_matrix::Bool

    function TranslationDFTTransform(n_sites::Integer, momentum::Integer, real_moment_matrix::Bool)
        _assert_positive_index("translation DFT site count", n_sites)
        k = Int(momentum)
        0 <= k < Int(n_sites) || throw(ArgumentError(
            "Translation DFT momentum must satisfy 0 <= k < n_sites; got $k for n_sites=$n_sites."
        ))
        return new(:translation_dft, :cyclotomic_float64, :cyclotomic, Int(n_sites), k, real_moment_matrix)
    end
end

struct TranslationReflectionTransform
    family::Symbol
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    n_sites::Int
    momentum::Int
    reflection::Int
    real_moment_matrix::Bool
    antiunitary::Bool
    base::TranslationDFTTransform

    function TranslationReflectionTransform(
        n_sites::Integer,
        momentum::Integer,
        reflection::Integer,
        real_moment_matrix::Bool,
        antiunitary::Bool,
    )
        parity = Int(reflection)
        parity in (-1, 1) || throw(ArgumentError(
            "Translation reflection parity must be -1 or 1; got $reflection."
        ))
        base = TranslationDFTTransform(n_sites, momentum, real_moment_matrix)
        return new(
            :translation_dft_reflection,
            :cyclotomic_float64,
            :cyclotomic_sqrt_rational,
            base.n_sites,
            base.momentum,
            parity,
            real_moment_matrix,
            Bool(antiunitary),
            base,
        )
    end
end

struct TranslationAxisIrrepTransform
    family::Symbol
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    base::Union{TranslationDFTTransform,TranslationReflectionTransform}
    irrep_label::Symbol
    irrep_dimension::Int
    irrep_multiplicity::Int
    axis_group_order::Int

    function TranslationAxisIrrepTransform(
        base::Union{TranslationDFTTransform,TranslationReflectionTransform},
        irrep_label::Symbol,
        irrep_dimension::Integer,
        irrep_multiplicity::Integer,
        axis_group_order::Integer,
    )
        family = base isa TranslationReflectionTransform ?
            :translation_dft_reflection_axis_irrep :
            :translation_dft_axis_irrep
        return new(
            family,
            :real_algebraic_float64,
            :axis_character_projector,
            base,
            irrep_label,
            Int(irrep_dimension),
            Int(irrep_multiplicity),
            Int(axis_group_order),
        )
    end
end

struct TranslationSU2RDMTransform
    family::Symbol
    coefficient_domain::Symbol
    exact_coefficient_domain::Symbol
    k::Int
    spin2::Int
    multiplicity::Int
    reference_m2::Int
    schur_matrix::Matrix{Float64}

    function TranslationSU2RDMTransform(
        k::Integer,
        block::PauliSU2RDMBlock,
        schur_matrix::Matrix{Float64},
    )
        kk = Int(k)
        kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))
        return new(
            :pauli_su2_rdm,
            :algebraic_float64,
            :sqrt_rational,
            kk,
            block.spin2,
            block.multiplicity,
            -block.spin2,
            schur_matrix,
        )
    end
end

struct TranslationInvariantBlockOrigin <: BlockOrigin
    label::Any
    logical_row_labels::Vector{Any}
    transform::Any

    function TranslationInvariantBlockOrigin(
        label,
        logical_row_labels::AbstractVector=Any[];
        transform=nothing,
    )
        return new(label, Any[row_label for row_label in logical_row_labels], transform)
    end
end

struct TranslationZeroOriginSeed
    label::Any
end

struct _SparseTransformRows{C<:Number}
    rows::Vector{Vector{Tuple{Int,C}}}
    ncols::Int
end

Base.size(transform::_SparseTransformRows) = (length(transform.rows), transform.ncols)
Base.size(transform::_SparseTransformRows, dim::Integer) = size(transform)[Int(dim)]

function _sparse_transform_rows(
    rows::Vector{Vector{Tuple{Int,C}}},
    ncols::Integer,
) where {C<:Number}
    cols = Int(ncols)
    cols >= 0 || throw(DomainError(ncols, "sparse transform column count must be non-negative."))
    for row in rows, (idx, _) in row
        1 <= idx <= cols || throw(BoundsError(1:cols, idx))
    end
    return _SparseTransformRows{C}(rows, cols)
end

function _pauli_sparse_transform_rows(
    matrix::AbstractMatrix{C};
    atol::Real=0,
) where {C<:Number}
    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(atol) : atol
    rows = Vector{Vector{Tuple{Int,C}}}(undef, size(matrix, 1))
    for (row_idx, i) in enumerate(axes(matrix, 1))
        row = Tuple{Int,C}[]
        for j in axes(matrix, 2)
            value = matrix[i, j]
            abs(value) <= tolerance && continue
            push!(row, (Int(j), value))
        end
        rows[row_idx] = row
    end
    return _sparse_transform_rows(rows, size(matrix, 2))
end

function _select_sparse_transform_rows(
    transform::_SparseTransformRows{C},
    rows::AbstractVector{<:Integer},
) where {C<:Number}
    selected = Vector{Vector{Tuple{Int,C}}}(undef, length(rows))
    for (out_idx, row) in pairs(rows)
        src_idx = Int(row)
        1 <= src_idx <= length(transform.rows) || throw(BoundsError(transform.rows, src_idx))
        selected[out_idx] = transform.rows[src_idx]
    end
    return _SparseTransformRows{C}(selected, transform.ncols)
end

struct TranslationZeroOrigin <: ConstraintOrigin
    label::Any
    row::Int
    col::Int
    part::Symbol

    function TranslationZeroOrigin(label, row::Integer, col::Integer, part::Symbol)
        _assert_positive_index("translation zero row", row)
        _assert_positive_index("translation zero col", col)
        part in (:real, :imag, :scalar) || throw(ArgumentError(
            "Unsupported translation zero part $(repr(part)); expected :real, :imag, or :scalar"
        ))
        return new(label, Int(row), Int(col), part)
    end
end

function _zero_constraint_origin(
    seed::TranslationZeroOriginSeed,
    constraint_idx::Integer,
    row::Integer,
    col::Integer,
    part::Symbol,
)
    return TranslationZeroOrigin(seed.label, row, col, part)
end

_translation_zero_seed_from_origin(origin::TranslationZeroOrigin) =
    TranslationZeroOriginSeed(origin.label)

_translation_zero_seed_from_origin(origin::PauliSU2SingletChannelEqualityOrigin) =
    PauliSU2SingletChannelEqualityOriginSeed(origin.label)

function _translation_zero_seed_from_origin(origin)
    throw(ArgumentError(
        "Translation zero provenance expected TranslationZeroOrigin or " *
        "PauliSU2SingletChannelEqualityOrigin, got $(typeof(origin)).",
    ))
end

_translation_zero_seed_label(seed::TranslationZeroOriginSeed) = seed.label
_translation_zero_seed_label(seed::PauliSU2SingletChannelEqualityOriginSeed) = seed.label

_translation_zero_seed_factory(::TranslationZeroOriginSeed) = TranslationZeroOriginSeed
_translation_zero_seed_factory(::PauliSU2SingletChannelEqualityOriginSeed) =
    PauliSU2SingletChannelEqualityOriginSeed

function _seed_translation_zero_origins_from_linear!(
    zero_origin_by_constraint::Dict{Int,Any},
    constraints,
    linear::MomentLinearData{K,C,M},
    ::Type{A},
) where {K,C,M,A<:AlgebraType}
    zero_constraints = linear.zero_constraints
    zero_idx = 1
    for (constraint_idx, (cone, mat)) in pairs(constraints)
        cone == :Zero || continue
        tmp = ScalarLinearConstraint{K,C}[]
        _append_zero_linear_constraints!(
            tmp,
            A,
            K,
            C,
            linear.adjoint_key,
            Int(constraint_idx),
            mat,
            Dict{Int,Any}(),
        )
        row_count = length(tmp)
        iszero(row_count) && continue
        zero_idx + row_count - 1 <= length(zero_constraints) || throw(ArgumentError(
            "Translation zero provenance has fewer linear zero rows than symbolic Zero constraints."
        ))
        origin = zero_constraints[zero_idx].origin
        zero_origin_by_constraint[Int(constraint_idx)] =
            _translation_zero_seed_from_origin(origin)
        zero_idx += row_count
    end
    zero_idx == length(zero_constraints) + 1 || throw(ArgumentError(
        "Translation zero provenance has extra linear zero rows after symbolic Zero constraints."
    ))
    return zero_origin_by_constraint
end

function _translation_label_field(label, key::Symbol)
    return label isa NamedTuple && haskey(label, key) ? getproperty(label, key) : nothing
end

function _translation_zero_origin_histogram_key(origin::TranslationZeroOrigin)
    label = origin.label
    return (
        feature=_translation_label_field(label, :feature),
        decomposition=_translation_label_field(label, :decomposition),
        reason=_translation_label_field(label, :reason),
    )
end

function _translation_zero_origin_histogram_key(
    origin::PauliSU2SingletChannelEqualityOrigin,
)
    label = origin.label
    return (
        feature=_translation_label_field(label, :feature),
        decomposition=_translation_label_field(label, :decomposition),
        reason=_translation_label_field(label, :reason),
    )
end

function _translation_zero_origin_histogram_key(origin::ConstraintOrigin)
    return (feature=Symbol(nameof(typeof(origin))), decomposition=nothing, reason=nothing)
end

"""
    translation_zero_origin_histogram(mp_or_linear)

Return sorted counts of zero/equality rows by translation provenance.  Histogram
keys are `(feature, decomposition, reason)`; missing label fields are reported as
`nothing`.
"""
function _translation_zero_origin_histogram(zero_constraints)
    counts = Dict{Any,Int}()
    for zero_constraint in zero_constraints
        key = _translation_zero_origin_histogram_key(zero_constraint.origin)
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by=pair -> repr(pair.first))
end

function _translation_zero_report_histogram(zero_constraints)
    return Pair{Any,Int}[
        Pair{Any,Int}(pair.first, pair.second) for
            pair in _translation_zero_origin_histogram(zero_constraints)
    ]
end

translation_zero_origin_histogram(mp::MomentProblem) =
    _translation_zero_origin_histogram(mp.linear.zero_constraints)

translation_zero_origin_histogram(linear::MomentLinearData) =
    _translation_zero_origin_histogram(linear.zero_constraints)

function _translation_origin_field(origin, field::Symbol)
    return hasproperty(origin, field) ? getproperty(origin, field) : nothing
end

function _translation_origin_vector(origin, field::Symbol)
    value = _translation_origin_field(origin, field)
    value === nothing && return nothing
    return Any[item for item in value]
end

function _translation_transform_field(transform, field::Symbol)
    transform === nothing && return nothing
    return hasproperty(transform, field) ? getproperty(transform, field) : nothing
end

function _translation_zero_origin_record(index::Integer, zero_constraint)
    origin = zero_constraint.origin
    label = _translation_origin_field(origin, :label)
    return (
        index=Int(index),
        kind=zero_constraint.kind,
        origin=origin,
        label=label,
        feature=_translation_label_field(label, :feature),
        decomposition=_translation_label_field(label, :decomposition),
        reason=_translation_label_field(label, :reason),
        row=_translation_origin_field(origin, :row),
        col=_translation_origin_field(origin, :col),
        part=_translation_origin_field(origin, :part),
        term_count=length(zero_constraint.form),
    )
end

"""
    translation_linear_provenance(mp_or_linear)

Return deterministic provenance records for translation-path linear data.
`psd_blocks` exposes each solver block's label, logical row labels, raw row
labels, transform descriptor, and coefficient domains.  `zero_constraints`
exposes equality/zero-row origins and linear-form term counts.
"""
function translation_linear_provenance(linear::MomentLinearData)
    psd_blocks = [
        begin
            origin = block.meta.origin
            transform = _translation_origin_field(origin, :transform)
            (
                index=Int(block_idx),
                cone=block.meta.cone,
                size=block.size,
                origin=origin,
                label=_translation_origin_field(origin, :label),
                logical_row_labels=_translation_origin_vector(origin, :logical_row_labels),
                row_labels=copy(block.meta.row_labels),
                transform=transform,
                transform_family=_translation_transform_field(transform, :family),
                coefficient_domain=_translation_transform_field(transform, :coefficient_domain),
                exact_coefficient_domain=_translation_transform_field(
                    transform,
                    :exact_coefficient_domain,
                ),
            )
        end
        for (block_idx, block) in pairs(linear.psd_blocks_lin)
    ]
    zero_constraints = [
        _translation_zero_origin_record(zero_idx, zero_constraint)
        for (zero_idx, zero_constraint) in pairs(linear.zero_constraints)
    ]
    return (psd_blocks=psd_blocks, zero_constraints=zero_constraints)
end

translation_linear_provenance(mp::MomentProblem) =
    translation_linear_provenance(mp.linear)

function TranslationInvariantReport(
    n_sites::Int,
    order::Int,
    basis_size::Int,
    orbit_basis_size::Int,
    momentum_sectors::Vector{Int},
    sign_symmetry::Bool,
    psd_block_sizes::Vector{Int},
    block_labels::Vector{Any},
    n_unique_moment_matrix_elements::Int,
    real_moment_matrix::Bool,
)
    return TranslationInvariantReport(
        n_sites,
        order,
        basis_size,
        orbit_basis_size,
        false,
        0,
        Pair{Int,Int}[],
        0.0,
        momentum_sectors,
        sign_symmetry,
        copy(psd_block_sizes),
        psd_block_sizes,
        block_labels,
        fill(nothing, length(block_labels)),
        fill(nothing, length(block_labels)),
        n_unique_moment_matrix_elements,
        n_unique_moment_matrix_elements,
        0,
        0,
        Dict{Symbol,Int}(),
        real_moment_matrix,
        0,
        0,
        0,
        0.0,
        Pair{Any,Int}[],
        false,
        0,
        0,
        0,
        0,
        0.0,
    )
end

function _translation_report_with_su2_moment_quotient(
    report::TranslationInvariantReport,
    quotient_result,
    ;
    construction_time_ns::Integer=report.construction_time_ns,
)
    descriptor = quotient_result.descriptor
    raw_count = descriptor.raw_moment_count
    quotient_count = descriptor.quotient_moment_count
    raw_count > 0 || throw(ArgumentError(
        "SU(2) moment quotient report requires a positive raw moment count",
    ))
    base_fields = ntuple(32) do idx
        idx == 17 && return quotient_count
        idx == 18 && return length(quotient_result.linear.zero_constraints)
        idx == 19 && return Int(construction_time_ns)
        idx == 26 && return _translation_zero_report_histogram(
            quotient_result.linear.zero_constraints,
        )
        return getfield(report, idx)
    end
    return TranslationInvariantReport(
        base_fields...,
        true,
        raw_count,
        quotient_count,
        quotient_count / raw_count,
        descriptor.support_orbit_count,
        copy(descriptor.singlet_channel_support_counts),
        descriptor.max_pivot_residual,
        descriptor.max_invariant_residual,
        descriptor.max_reconstruction_residual,
        descriptor.max_condition,
        quotient_result.eliminated_zero_row_count,
        copy(quotient_result.eliminated_zero_feature_histogram),
    )
end

function _apply_translation_su2_moment_quotient(
    linear::MomentLinearData,
    report::TranslationInvariantReport,
    n_sites::Integer;
    enabled::Bool,
    allow_registered_projection::Bool,
    atol::Real,
    condition_limit::Real,
    construction_start_ns::Integer,
    stage_times_ns::Dict{Symbol,Int},
)
    enabled || return (linear=linear, report=report)
    quotient_start_ns = time_ns()
    quotient_result = _pauli_su2_quotient_linear_data(
        linear,
        n_sites;
        sign_symmetry=report.sign_symmetry,
        allow_registered_projection,
        atol,
        condition_limit,
        stage_times_ns,
    )
    quotient_done_ns = time_ns()
    stage_times_ns[:su2_moment_quotient] = Int(quotient_done_ns - quotient_start_ns)
    quotient_report = _translation_report_with_su2_moment_quotient(
        report,
        quotient_result;
        construction_time_ns=Int(quotient_done_ns - construction_start_ns),
    )
    return (linear=quotient_result.linear, report=quotient_report)
end

const _EMPTY_AXIS_ROTATION_QUOTIENT_STATS = (
    equality_row_count=0,
    raw_moment_key_count=0,
    moment_class_count=0,
    quotient_moment_key_count=0,
    forced_zero_moment_class_count=0,
    reduction_ratio=0.0,
)

function _translation_block_meta_by_constraint(
    ::Type{M},
    constraints::Vector{Tuple{Symbol,Matrix{P}}},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
) where {A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P<:Polynomial{A,T}}
    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    block_idx = 0
    for (constraint_idx, (cone, mat)) in pairs(constraints)
        cone in (:PSD, :HPSD) || continue
        block_idx += 1
        block_idx <= length(block_labels) || throw(ArgumentError(
            "Translation block metadata is missing a label for PSD/HPSD constraint $constraint_idx."
        ))
        block_idx <= length(block_logical_row_labels) || throw(ArgumentError(
            "Translation block metadata is missing logical row labels for PSD/HPSD constraint $constraint_idx."
        ))
        block_idx <= length(block_transforms) || throw(ArgumentError(
            "Translation block metadata is missing transform provenance for PSD/HPSD constraint $constraint_idx."
        ))
        logical_rows = block_logical_row_labels[block_idx]
        (length(logical_rows) == size(mat, 1) || 2 * length(logical_rows) == size(mat, 1)) || throw(
            DimensionMismatch(
                "Translation block $(block_labels[block_idx]) has $(length(logical_rows)) logical row labels but solver-facing size $(size(mat, 1))."
            )
        )
        block_meta_by_constraint[constraint_idx] = BlockMeta{M}(
            cone,
            TranslationInvariantBlockOrigin(
                block_labels[block_idx],
                logical_rows;
                transform=block_transforms[block_idx],
            ),
            fill(one(M), size(mat, 1)),
        )
    end
    block_idx == length(block_labels) || throw(ArgumentError(
        "Translation report has $(length(block_labels)) block labels but constraints contain $block_idx PSD/HPSD blocks."
    ))
    block_idx == length(block_logical_row_labels) || throw(ArgumentError(
        "Translation report has $(length(block_logical_row_labels)) block row-label groups but constraints contain $block_idx PSD/HPSD blocks."
    ))
    block_idx == length(block_transforms) || throw(ArgumentError(
        "Translation report has $(length(block_transforms)) block transform descriptors but constraints contain $block_idx PSD/HPSD blocks."
    ))
    return block_meta_by_constraint
end

function _translation_block_coefficient_domains(block_transforms::Vector{Any})
    return Union{Nothing,Symbol}[
        transform === nothing ? nothing : transform.coefficient_domain
        for transform in block_transforms
    ]
end

function _translation_block_exact_coefficient_domains(block_transforms::Vector{Any})
    return Union{Nothing,Symbol}[
        transform === nothing ? nothing : transform.exact_coefficient_domain
        for transform in block_transforms
    ]
end

function _translation_su2_base_report(
    su2_mp::MomentProblem,
    n::Int,
    d::Int,
    local_basis_size::Int,
    orbit_reps::AbstractVector,
    axis_orbit_summary,
    sectors::Vector{Int},
    sign_symmetry::Bool,
    construction_time_ns::Integer,
    construction_stage_times_ns::Dict{Symbol,Int},
    real_moment_matrix::Bool,
)
    block_labels = Any[block.meta.origin.label for block in su2_mp.linear.psd_blocks_lin]
    logical_block_sizes = Int[
        length(block.meta.origin.logical_row_labels)
        for block in su2_mp.linear.psd_blocks_lin
    ]
    psd_block_sizes = Int[block.size for block in su2_mp.linear.psd_blocks_lin]
    block_transforms = Any[
        block.meta.origin.transform for block in su2_mp.linear.psd_blocks_lin
    ]
    return TranslationInvariantReport(
        n,
        d,
        local_basis_size,
        length(orbit_reps),
        axis_orbit_summary.axis_orbit_closed,
        axis_orbit_summary.axis_orbit_count,
        axis_orbit_summary.axis_orbit_size_histogram,
        axis_orbit_summary.reduction_ratio,
        sectors,
        sign_symmetry,
        logical_block_sizes,
        psd_block_sizes,
        block_labels,
        _translation_block_coefficient_domains(block_transforms),
        _translation_block_exact_coefficient_domains(block_transforms),
        su2_mp.n_unique_moment_matrix_elements,
        length(su2_mp.linear.moments),
        length(su2_mp.linear.zero_constraints),
        Int(construction_time_ns),
        construction_stage_times_ns,
        real_moment_matrix,
        0,
        0,
        0,
        0.0,
        _translation_zero_report_histogram(su2_mp.linear.zero_constraints),
        false,
        0,
        0,
        0,
        0,
        0.0,
    )
end

function _translation_su2_base_report(
    linear::MomentLinearData,
    n::Int,
    d::Int,
    local_basis_size::Int,
    orbit_reps::AbstractVector,
    axis_orbit_summary,
    sectors::Vector{Int},
    sign_symmetry::Bool,
    construction_time_ns::Integer,
    construction_stage_times_ns::Dict{Symbol,Int},
    real_moment_matrix::Bool,
    product_cache_hits::Integer=0,
    product_cache_misses::Integer=0,
    product_cache_entries::Integer=0,
    product_cache_hit_rate::Real=0.0,
)
    block_labels = Any[block.meta.origin.label for block in linear.psd_blocks_lin]
    logical_block_sizes = Int[
        length(block.meta.origin.logical_row_labels)
        for block in linear.psd_blocks_lin
    ]
    psd_block_sizes = Int[block.size for block in linear.psd_blocks_lin]
    block_transforms = Any[
        block.meta.origin.transform for block in linear.psd_blocks_lin
    ]
    return TranslationInvariantReport(
        n,
        d,
        local_basis_size,
        length(orbit_reps),
        axis_orbit_summary.axis_orbit_closed,
        axis_orbit_summary.axis_orbit_count,
        axis_orbit_summary.axis_orbit_size_histogram,
        axis_orbit_summary.reduction_ratio,
        sectors,
        sign_symmetry,
        logical_block_sizes,
        psd_block_sizes,
        block_labels,
        _translation_block_coefficient_domains(block_transforms),
        _translation_block_exact_coefficient_domains(block_transforms),
        length(linear.moments),
        length(linear.moments),
        length(linear.zero_constraints),
        Int(construction_time_ns),
        construction_stage_times_ns,
        real_moment_matrix,
        Int(product_cache_hits),
        Int(product_cache_misses),
        Int(product_cache_entries),
        Float64(product_cache_hit_rate),
        _translation_zero_report_histogram(linear.zero_constraints),
        false,
        0,
        0,
        0,
        0,
        0.0,
    )
end

function _translation_report_with_linear(
    base_report::TranslationInvariantReport,
    linear::MomentLinearData,
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_transforms::Vector{Any},
    construction_time_ns::Integer,
    construction_stage_times_ns::Dict{Symbol,Int},
)
    return TranslationInvariantReport(
        base_report.n_sites,
        base_report.order,
        base_report.basis_size,
        base_report.orbit_basis_size,
        base_report.axis_orbit_closed,
        base_report.axis_orbit_basis_size,
        base_report.axis_orbit_size_histogram,
        base_report.axis_reduction_ratio,
        base_report.momentum_sectors,
        base_report.sign_symmetry,
        logical_block_sizes,
        block_sizes,
        block_labels,
        _translation_block_coefficient_domains(block_transforms),
        _translation_block_exact_coefficient_domains(block_transforms),
        base_report.n_unique_moment_matrix_elements,
        length(linear.moments),
        length(linear.zero_constraints),
        Int(construction_time_ns),
        construction_stage_times_ns,
        base_report.real_moment_matrix,
        base_report.product_cache_hits,
        base_report.product_cache_misses,
        base_report.product_cache_entries,
        base_report.product_cache_hit_rate,
        _translation_zero_report_histogram(linear.zero_constraints),
        base_report.axis_rotation_quotient,
        base_report.axis_rotation_raw_moment_key_count,
        base_report.axis_rotation_moment_class_count,
        base_report.axis_rotation_quotient_moment_key_count,
        base_report.axis_rotation_forced_zero_moment_class_count,
        base_report.axis_rotation_moment_quotient_reduction_ratio,
    )
end

function _pauli_su2_linear_block_meta(
    ::Type{M},
    cone::Symbol,
    block,
    block_size::Integer,
) where {M<:NormalMonomial{PauliAlgebra}}
    return BlockMeta{M}(
        cone,
        PauliSU2BasisMomentOrigin(
            block.label,
            block.row_labels;
            transform=block.transform,
        ),
        fill(one(M), Int(block_size)),
    )
end

function _append_pauli_su2_linear_psd_block!(
    builder::MomentLinearBuilder{K,LC,M},
    block,
    entries::Matrix{LinearMomentForm{K,LC}},
    cone::Symbol,
    constraint_idx::Int,
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_transforms::Vector{Any},
) where {K,LC,M<:NormalMonomial{PauliAlgebra}}
    add_psd_block!(
        builder,
        cone,
        entries,
        _pauli_su2_linear_block_meta(M, cone, block, size(entries, 1));
        constraint_idx,
    )
    push!(logical_block_sizes, length(block.row_labels))
    push!(block_sizes, size(entries, 1))
    push!(block_labels, block.label)
    push!(block_transforms, block.transform)
    return nothing
end

function _pauli_su2_translation_base_linear_relaxation(
    objective::Polynomial{PauliAlgebra,T,C0},
    orbit_reps::Vector{M},
    n::Int,
    sectors::Vector{Int};
    real_moment_matrix::Bool,
    reflection_symmetry::Bool,
    phase_atol::R,
    max_dimension::Integer=4096,
    reference_m2=:lowest,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
    emit_invariance_rows::Bool=true,
) where {T<:Unsigned,C0<:Number,M<:NormalMonomial{PauliAlgebra,T},R<:Real}
    !real_moment_matrix && reflection_symmetry && throw(ArgumentError(
        "Complex direct-linear translation/SU(2) base construction does not " *
        "yet support reflection_symmetry=true."
    ))
    CComplex = Complex{R}
    LC = real_moment_matrix ? R : CComplex
    K = typeof(symmetric_canon(expval(one(M))))
    MP_R = Polynomial{PauliAlgebra,T,LC}
    objective_mp = convert(MP_R, _translation_orbit_reduce_polynomial(objective, n))

    builder = MomentLinearBuilder(K, LC, M)
    _add_translation_objective_terms!(builder, objective_mp)
    registered_key_tokens = _translation_registered_key_tokens(builder)

    translated = Dict{M,Vector{M}}()
    for rep in orbit_reps
        translated[rep] = isone(rep) ?
            fill(rep, n) :
            [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
    end
    nontrivial_reps = M[rep for rep in orbit_reps if !isone(rep)]
    rep_cache = Dict{M,M}()
    product_cache = _TranslationProductCache(
        Dict{Tuple{M,M},Vector{Tuple{Int,CComplex,M}}}(),
        0,
        0,
    )

    logical_block_sizes = Int[]
    block_sizes = Int[]
    block_labels = Any[]
    block_transforms = Any[]
    constraint_idx = 0
    cone = real_moment_matrix ? :PSD : :HPSD
    atol = R(phase_atol)
    base_stage_times_ns = Dict{Symbol,Int}()
    for k in sectors
        sector_basis = iszero(k) ? orbit_reps : nontrivial_reps
        isempty(sector_basis) && continue
        stage_start_ns = time_ns()
        complex_entries = _translation_momentum_block_linear_entries(
            sector_basis,
            k,
            n,
            translated,
            rep_cache,
            K,
            CComplex,
            product_cache;
            real_moment_matrix=false,
            atol,
        )
        base_stage_times_ns[:su2_base_entries] =
            get(base_stage_times_ns, :su2_base_entries, 0) + Int(time_ns() - stage_start_ns)
        if emit_invariance_rows
            stage_start_ns = time_ns()
            _register_translation_linear_entries!(
                builder,
                complex_entries,
                registered_key_tokens,
            )
            base_stage_times_ns[:su2_base_register_entries] =
                get(base_stage_times_ns, :su2_base_register_entries, 0) +
                Int(time_ns() - stage_start_ns)
        end
        stage_start_ns = time_ns()
        transform_blocks = pauli_su2_translation_orbit_basis_transform_blocks(
            sector_basis,
            n;
            max_dimension,
            momentum=k,
        )
        base_stage_times_ns[:su2_base_transforms] =
            get(base_stage_times_ns, :su2_base_transforms, 0) + Int(time_ns() - stage_start_ns)
        sector_reference_m2 = reflection_symmetry && !_translation_reflection_fixed_momentum(k, n) ?
            :zero :
            reference_m2
        stage_start_ns = time_ns()
        reduced_blocks = _append_pauli_su2_translation_wigner_eckart_linear!(
            builder,
            complex_entries,
            transform_blocks,
            k;
            reference_m2=sector_reference_m2,
            phase_atol=atol,
            stage_times_ns=base_stage_times_ns,
            registered_key_tokens,
            register_keys=false,
            emit_invariance_rows,
        )
        base_stage_times_ns[:su2_base_wigner_rows] =
            get(base_stage_times_ns, :su2_base_wigner_rows, 0) + Int(time_ns() - stage_start_ns)

        if singlet_channel_equalities && emit_invariance_rows
            stage_start_ns = time_ns()
            _append_pauli_su2_translation_singlet_channel_equality_linear_constraints!(
                builder,
                sector_basis,
                n,
                k;
                max_dimension,
                phase_atol=atol,
                atol=singlet_channel_atol,
                registered_key_tokens,
            )
            base_stage_times_ns[:su2_base_singlet_channel_rows] =
                get(base_stage_times_ns, :su2_base_singlet_channel_rows, 0) +
                Int(time_ns() - stage_start_ns)
        end

        stage_start_ns = time_ns()
        sector_blocks = if reflection_symmetry
            if _translation_reflection_fixed_momentum(k, n)
                _pauli_su2_translation_reflection_fixed_linear_blocks(
                    reduced_blocks,
                    sector_basis,
                    k,
                    n;
                    atol=Float64(phase_atol),
                )
            else
                _pauli_su2_translation_real_reflection_linear_blocks_from_reduced(
                    reduced_blocks,
                    sector_basis,
                    k,
                    n,
                    LC;
                    phase_atol=atol,
                )
            end
        else
            reduced_blocks
        end
        base_stage_times_ns[:su2_base_reflection] =
            get(base_stage_times_ns, :su2_base_reflection, 0) + Int(time_ns() - stage_start_ns)

        stage_start_ns = time_ns()
        for block in sector_blocks
            entries = if reflection_symmetry && !_translation_reflection_fixed_momentum(k, n)
                block.matrix
            elseif real_moment_matrix
                _symmetrize_real_linear_block(
                    _realify_hermitian_linear_block_full(block.matrix, LC; atol),
                )
            else
                _hermitianize_pauli_linear_block(block.matrix)
            end
            entries = _prepare_translation_linear_entries!(
                builder,
                entries,
                registered_key_tokens,
                nothing,
            )
            constraint_idx += 1
            _append_pauli_su2_linear_psd_block!(
                builder,
                block,
                entries,
                cone,
                constraint_idx,
                logical_block_sizes,
                block_sizes,
                block_labels,
                block_transforms,
            )
        end
        base_stage_times_ns[:su2_base_psd_append] =
            get(base_stage_times_ns, :su2_base_psd_append, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    if emit_invariance_rows || singlet_channel_equalities
        _deduplicate_translation_zero_constraints!(builder)
        base_stage_times_ns[:su2_base_finalize_deduplicate_zero_constraints] =
            get(base_stage_times_ns, :su2_base_finalize_deduplicate_zero_constraints, 0) +
            Int(time_ns() - stage_start_ns)
    end
    dedup_done_ns = time_ns()
    linear = finalize!(
        builder;
        stage_times_ns=base_stage_times_ns,
        stage_prefix=:su2_base_finalize,
    )
    finalize_done_ns = time_ns()
    base_stage_times_ns[:su2_base_finalize_linear] =
        get(base_stage_times_ns, :su2_base_finalize_linear, 0) +
        Int(finalize_done_ns - dedup_done_ns)
    base_stage_times_ns[:su2_base_finalize] =
        get(base_stage_times_ns, :su2_base_finalize, 0) + Int(finalize_done_ns - stage_start_ns)
    return (
        linear=linear,
        logical_block_sizes=logical_block_sizes,
        block_sizes=block_sizes,
        block_labels=block_labels,
        block_transforms=block_transforms,
        stage_times_ns=base_stage_times_ns,
        product_cache_hits=product_cache.hits,
        product_cache_misses=product_cache.misses,
        product_cache_entries=length(product_cache.terms),
        product_cache_hit_rate=_translation_product_cache_hit_rate(product_cache),
    )
end

function _append_translation_direct_extend_addons_to_su2_linear(
    pop::PolyOpt{PauliAlgebra,T,P},
    base_linear::MomentLinearData{K,LC,M},
    base_report::TranslationInvariantReport,
    n_sites::Integer;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    linear_state_width::Int,
    linear_state_opt_mode::Symbol,
    psd_state_width::Int,
    rdm_ks::Vector{Int},
    sectors::Vector{Int},
    phase_atol::R,
    stage_times_ns::Union{Nothing,Dict{Symbol,Int}}=nothing,
) where {
    K,
    LC,
    R<:Real,
    T<:Unsigned,
    C0<:Number,
    P<:Polynomial{PauliAlgebra,T,C0},
    M<:NormalMonomial{PauliAlgebra,T},
}
    n = Int(n_sites)
    stage_start_ns = time_ns()
    builder = _translation_append_builder_from_linear(base_linear)
    registered_key_tokens = _translation_registered_key_tokens(builder)
    constraint_moment_basis = _translation_builder_psd_moment_basis(builder)
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_clone] =
            get(stage_times_ns, :su2_extend_clone, 0) + Int(time_ns() - stage_start_ns)
    end
    MP_R = _pauli_chain_real_coeff_type(C0)
    MP_C = Complex{MP_R}
    MP_P = Polynomial{PauliAlgebra,T,LC}
    cone = real_moment_matrix ? :PSD : :HPSD

    logical_block_sizes = copy(base_report.logical_block_sizes)
    block_sizes = copy(base_report.psd_block_sizes)
    block_labels = copy(base_report.block_labels)
    block_logical_row_labels = Vector{Any}[
        Any[row for row in block.meta.origin.logical_row_labels]
        for block in base_linear.psd_blocks_lin
    ]
    block_transforms = Any[
        block.meta.origin.transform for block in base_linear.psd_blocks_lin
    ]
    constraint_idx = isempty(base_linear.psd_block_constraint_idx) ?
        0 :
        maximum(base_linear.psd_block_constraint_idx)

    stage_start_ns = time_ns()
    _add_translation_linear_state_opt_linear_constraints!(
        builder,
        pop.objective,
        n,
        linear_state_width,
        sign_symmetry,
        constraint_moment_basis,
        MP_P,
        true;
        registered_key_tokens,
        linear_state_opt_mode,
    )
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_linear_state_opt] =
            get(stage_times_ns, :su2_extend_linear_state_opt, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    if psd_state_width > 0
        state_opt_basis = _contiguous_state_opt_tests(M, n, psd_state_width; sign_symmetry)
        if !isempty(state_opt_basis)
            psd_state_opt_translated = Dict{M,Vector{M}}()
            for rep in state_opt_basis
                psd_state_opt_translated[rep] = [
                    _translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)
                ]
            end
            state_opt_term_cache = _translation_psd_state_opt_term_cache(
                state_opt_basis,
                n,
                pop.objective,
                psd_state_opt_translated,
            )
            for k in sectors
                complex_entries = _hermitianize_pauli_linear_block(
                    _translation_psd_state_opt_linear_block(
                        state_opt_basis,
                        k,
                        n,
                        state_opt_term_cache,
                        K,
                        MP_C,
                    ),
                )
                entries = real_moment_matrix ?
                    _realify_hermitian_linear_block(complex_entries, LC; atol=MP_R(phase_atol)) :
                    complex_entries
                label = (feature=:psd_state_opt, width=psd_state_width, momentum=k)
                logical_rows = Any[mono for mono in state_opt_basis]
                transform = TranslationDFTTransform(n, k, real_moment_matrix)
                entries = _prepare_translation_linear_entries!(
                    builder,
                    entries,
                    registered_key_tokens,
                    nothing,
                )
                constraint_idx += 1
                add_psd_block!(
                    builder,
                    cone,
                    entries,
                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                    constraint_idx,
                )
                push!(logical_block_sizes, size(complex_entries, 1))
                push!(block_sizes, size(entries, 1))
                push!(block_labels, label)
                push!(block_logical_row_labels, logical_rows)
                push!(block_transforms, transform)
            end
        end
    end
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_psd_state_opt] =
            get(stage_times_ns, :su2_extend_psd_state_opt, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    for rdm_k in rdm_ks
        constraint_idx = _add_translation_su2_extend_rdm_linear_blocks!(
            builder,
            n,
            rdm_k,
            MP_C,
            real_moment_matrix,
            MP_R(phase_atol),
            registered_key_tokens,
            nothing,
            logical_block_sizes,
            block_sizes,
            block_labels,
            block_logical_row_labels,
            block_transforms,
            constraint_idx,
            stage_times_ns,
        )
    end
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_rdm] =
            get(stage_times_ns, :su2_extend_rdm, 0) + Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    linear = finalize!(
        builder;
        stage_times_ns,
        stage_prefix=:su2_extend_finalize,
    )
    finalize_ns = Int(time_ns() - stage_start_ns)
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_finalize_linear] =
            get(stage_times_ns, :su2_extend_finalize_linear, 0) + finalize_ns
        stage_times_ns[:su2_extend_finalize] =
            get(stage_times_ns, :su2_extend_finalize, 0) + finalize_ns
    end
    return (
        linear=linear,
        logical_block_sizes=logical_block_sizes,
        block_sizes=block_sizes,
        block_labels=block_labels,
        block_transforms=block_transforms,
    )
end

function _append_translation_extra_constraints_to_su2_base_moment_problem(
    pop::PolyOpt{PauliAlgebra,T,P},
    su2_mp::MomentProblem{PauliAlgebra,T,M,MP},
    orbit_reps::Vector{M},
    n_sites::Integer;
    real_moment_matrix::Bool,
    sign_symmetry::Bool,
    linear_state_width::Int,
    axis_rotation_equalities::Bool,
    rdm_ks::Vector{Int}=Int[],
    rdm_decomposition::Symbol=:full,
    rdm_support::Symbol=:closed,
    psd_state_width::Int=0,
    sectors::Vector{Int}=Int[],
    phase_atol::Real=1e-12,
    ops=nothing,
) where {
    T<:Unsigned,
    P<:Polynomial{PauliAlgebra,T},
    M<:NormalMonomial{PauliAlgebra,T},
    CMP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
}
    (
        isempty(pop.eq_constraints) &&
        isempty(pop.ineq_constraints) &&
        isempty(pop.moment_eq_constraints) &&
        isempty(rdm_ks) &&
        !axis_rotation_equalities &&
        iszero(linear_state_width) &&
        iszero(psd_state_width)
    ) && return su2_mp
    ops === nothing && throw(ArgumentError(
        "Translation/SU(2) extra constraints require Pauli chain operators."
    ))

    constraints = copy(su2_mp.constraints)
    existing_blocks = su2_mp.linear.psd_blocks_lin
    logical_block_sizes = Int[
        length(block.meta.origin.logical_row_labels) for block in existing_blocks
    ]
    block_sizes = Int[block.size for block in existing_blocks]
    block_labels = Any[block.meta.origin.label for block in existing_blocks]
    block_logical_row_labels = Vector{Any}[
        Any[row for row in block.meta.origin.logical_row_labels] for block in existing_blocks
    ]
    block_transforms = Any[block.meta.origin.transform for block in existing_blocks]
    zero_origin_by_constraint = Dict{Int,Any}()
    _seed_translation_zero_origins_from_linear!(
        zero_origin_by_constraint,
        constraints,
        su2_mp.linear,
        PauliAlgebra,
    )
    moment_basis = _collect_pauli_su2_moment_basis(
        zero(su2_mp.objective),
        su2_mp.constraints,
    )

    if axis_rotation_equalities
        _append_translation_axis_rotation_constraints!(
            constraints,
            ops,
            n_sites,
            moment_basis,
            MP,
            zero_origin_by_constraint,
        )
    end
    _append_translation_scalar_constraints!(
        constraints,
        pop,
        n_sites,
        moment_basis,
        MP,
        logical_block_sizes,
        block_sizes,
        block_labels,
        block_logical_row_labels,
        block_transforms,
        zero_origin_by_constraint,
    )
    _append_translation_moment_eq_constraints!(
        constraints,
        pop,
        n_sites,
        orbit_reps,
        sign_symmetry,
        moment_basis,
        MP,
        zero_origin_by_constraint,
    )
    _append_translation_linear_state_opt_constraints!(
        constraints,
        pop.objective,
        n_sites,
        linear_state_width,
        sign_symmetry,
        moment_basis,
        MP,
        zero_origin_by_constraint,
    )
    if psd_state_width > 0 || !isempty(rdm_ks)
        MP_R = _pauli_chain_real_coeff_type(CMP)
        BLOCK_P = Polynomial{PauliAlgebra,T,_pauli_chain_complex_coeff_type(CMP)}
        if psd_state_width > 0
            _append_translation_psd_state_opt_constraints!(
                constraints,
                pop.objective,
                n_sites,
                psd_state_width,
                sectors,
                sign_symmetry,
                moment_basis,
                BLOCK_P,
                MP,
                real_moment_matrix,
                MP_R(phase_atol),
                logical_block_sizes,
                block_sizes,
                block_labels,
                block_logical_row_labels,
                block_transforms,
            )
        end
    end
    if !isempty(rdm_ks)
        rdm_support == :closed || throw(ArgumentError(
            "Translation/SU(2) base moment reduction currently combines with RDM only for `contiguous_rdm_support=:closed`; got $(repr(rdm_support))."
        ))
        MP_R = _pauli_chain_real_coeff_type(CMP)
        BLOCK_P = Polynomial{PauliAlgebra,T,_pauli_chain_complex_coeff_type(CMP)}
        _append_translation_contiguous_rdm_constraints!(
            constraints,
            n_sites,
            rdm_ks,
            rdm_decomposition,
            rdm_support,
            moment_basis,
            BLOCK_P,
            MP,
            real_moment_matrix,
            MP_R(phase_atol),
            zero_origin_by_constraint,
            logical_block_sizes,
            block_sizes,
            block_labels,
            block_logical_row_labels,
            block_transforms,
        )
    end

    block_meta_by_constraint = Dict{Int,BlockMeta{M}}()
    psd_idx = 0
    for (constraint_idx, (cone, mat)) in pairs(constraints)
        cone in (:PSD, :HPSD) || continue
        psd_idx += 1
        if psd_idx <= length(existing_blocks)
            block_meta_by_constraint[constraint_idx] = existing_blocks[psd_idx].meta
        else
            block_meta_by_constraint[constraint_idx] = _translation_block_meta(
                M,
                cone,
                size(mat, 1),
                block_labels[psd_idx],
                block_logical_row_labels[psd_idx],
                block_transforms[psd_idx],
            )
        end
    end
    psd_idx == length(block_labels) || throw(ArgumentError(
        "Translation/SU(2) scalar constraint metadata has $(length(block_labels)) block labels but constraints contain $psd_idx PSD/HPSD blocks."
    ))

    total_basis = _collect_pauli_su2_moment_basis(su2_mp.objective, constraints)
    return MomentProblem{PauliAlgebra,T,M,MP}(
        su2_mp.objective,
        constraints,
        total_basis,
        su2_mp.n_unique_moment_matrix_elements;
        block_meta_by_constraint,
        zero_origin_by_constraint,
        real_moments=real_moment_matrix,
    )
end

function _translation_transform_descriptor(
    n_sites::Integer,
    momentum::Integer,
    real_moment_matrix::Bool,
    label,
)
    if haskey(label, :reflection)
        return TranslationReflectionTransform(
            n_sites,
            momentum,
            label.reflection,
            real_moment_matrix,
            !_translation_reflection_fixed_momentum(Int(momentum), Int(n_sites)),
        )
    end
    return TranslationDFTTransform(n_sites, momentum, real_moment_matrix)
end

function Base.show(io::IO, report::TranslationInvariantReport)
    print(
        io,
        "TranslationInvariantReport(n_sites=$(report.n_sites), order=$(report.order), " *
        "basis_size=$(report.basis_size), orbit_basis_size=$(report.orbit_basis_size), " *
        "axis_orbit_closed=$(report.axis_orbit_closed), " *
        "axis_orbit_basis_size=$(report.axis_orbit_basis_size), " *
        "axis_reduction_ratio=$(report.axis_reduction_ratio), " *
        "momentum_sectors=$(report.momentum_sectors), sign_symmetry=$(report.sign_symmetry), " *
        "logical_block_sizes=$(report.logical_block_sizes), " *
        "psd_block_sizes=$(report.psd_block_sizes), " *
        "block_coefficient_domains=$(report.block_coefficient_domains), " *
        "block_exact_coefficient_domains=$(report.block_exact_coefficient_domains), " *
        "n_unique_moment_matrix_elements=$(report.n_unique_moment_matrix_elements), " *
        "linear_moment_count=$(report.linear_moment_count), " *
        "zero_constraint_count=$(report.zero_constraint_count), " *
        "construction_time_ns=$(report.construction_time_ns), " *
        "construction_stage_times_ns=$(report.construction_stage_times_ns), " *
        "real_moment_matrix=$(report.real_moment_matrix), " *
        "product_cache_entries=$(report.product_cache_entries), " *
        "product_cache_hit_rate=$(report.product_cache_hit_rate))"
    )
end

"""
    translation_block_histogram(report; kind=:logical)

Return a sorted block-size histogram for a [`TranslationInvariantReport`](@ref).
Use `kind=:logical` for complex logical block sizes, or `kind=:psd` for the
solver-facing PSD block sizes.
"""
function translation_block_histogram(report::TranslationInvariantReport; kind::Symbol=:logical)
    sizes = if kind == :logical
        report.logical_block_sizes
    elseif kind == :psd
        report.psd_block_sizes
    else
        throw(ArgumentError("Unsupported translation block histogram kind $(repr(kind)); expected :logical or :psd."))
    end

    counts = Dict{Int,Int}()
    for size in sizes
        counts[size] = get(counts, size, 0) + 1
    end
    return sort!(collect(counts); by=first)
end

function _translation_block_feature_histogram_key(label, size::Integer)
    feature = _translation_label_field(label, :feature)
    decomposition = _translation_label_field(label, :decomposition)
    if feature === nothing
        feature = :moment_matrix
        decomposition = if _translation_label_field(label, :axis_irrep) !== nothing
            _translation_label_field(label, :reflection) === nothing ?
                :translation_axis_irrep :
                :translation_reflection_axis_irrep
        elseif _translation_label_field(label, :reflection) === nothing
            :translation
        else
            :translation_reflection
        end
    end
    return (feature=feature, decomposition=decomposition, size=Int(size))
end

"""
    translation_block_feature_histogram(report; kind=:logical)

Return a sorted histogram of translation-invariant PSD/HPSD block sizes grouped
by block feature and decomposition.  Use `kind=:logical` for complex logical
block sizes, or `kind=:psd` for the solver-facing PSD block sizes.
"""
function translation_block_feature_histogram(
    report::TranslationInvariantReport;
    kind::Symbol=:logical,
)
    sizes = if kind == :logical
        report.logical_block_sizes
    elseif kind == :psd
        report.psd_block_sizes
    else
        throw(ArgumentError("Unsupported translation block feature histogram kind $(repr(kind)); expected :logical or :psd."))
    end

    length(sizes) == length(report.block_labels) || throw(ArgumentError(
        "Translation report has $(length(sizes)) block sizes but $(length(report.block_labels)) labels."
    ))

    counts = Dict{Any,Int}()
    for (label, size) in zip(report.block_labels, sizes)
        key = _translation_block_feature_histogram_key(label, size)
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by=pair -> repr(pair.first))
end

function _translation_report_has_label_field(report::TranslationInvariantReport, key::Symbol)
    return any(
        label -> label isa NamedTuple && haskey(label, key),
        report.block_labels,
    )
end

function _translation_report_has_feature(
    report::TranslationInvariantReport,
    feature::Symbol;
    decomposition=nothing,
)
    return any(report.block_labels) do label
        label isa NamedTuple || return false
        haskey(label, :feature) && label.feature == feature || return false
        decomposition === nothing && return true
        return haskey(label, :decomposition) && label.decomposition == decomposition
    end
end

function _translation_report_has_zero_feature(
    report::TranslationInvariantReport,
    feature::Symbol;
    decomposition=nothing,
    reason=nothing,
)
    return any(report.zero_constraint_feature_histogram) do pair
        label = first(pair)
        label isa NamedTuple || return false
        haskey(label, :feature) && label.feature == feature || return false
        if decomposition !== nothing
            haskey(label, :decomposition) && label.decomposition == decomposition || return false
        end
        if reason !== nothing
            haskey(label, :reason) && label.reason == reason || return false
        end
        return true
    end
end

function _translation_report_zero_feature_count(
    report::TranslationInvariantReport,
    feature::Symbol,
    ;
    decomposition=nothing,
    reason=nothing,
)
    count = 0
    for pair in report.zero_constraint_feature_histogram
        label = first(pair)
        label isa NamedTuple || continue
        haskey(label, :feature) && label.feature == feature || continue
        if decomposition !== nothing
            haskey(label, :decomposition) && label.decomposition == decomposition || continue
        end
        if reason !== nothing
            haskey(label, :reason) && label.reason == reason || continue
        end
        count += last(pair)
    end
    return count
end

function _translation_report_feature_block_sizes(
    report::TranslationInvariantReport,
    feature::Symbol;
    kind::Symbol,
)
    sizes = if kind == :logical
        report.logical_block_sizes
    elseif kind == :psd
        report.psd_block_sizes
    else
        throw(ArgumentError("Unsupported translation feature block size kind $(repr(kind)); expected :logical or :psd."))
    end
    length(sizes) == length(report.block_labels) || throw(ArgumentError(
        "Translation report has $(length(sizes)) block sizes but $(length(report.block_labels)) labels."
    ))
    return Int[
        Int(size) for (label, size) in zip(report.block_labels, sizes) if
            label isa NamedTuple && haskey(label, :feature) && label.feature == feature
    ]
end

"""
    translation_symmetry_profile(report)

Summarize which symmetry reductions are actually represented in a
[`TranslationInvariantReport`](@ref).  The specialized Pauli-chain backend always
uses translation; sign and reflection are optional.  Global Pauli-axis rotations
are reported when the momentum blocks carry axis-irrep labels.
"""
function translation_symmetry_profile(report::TranslationInvariantReport)
    reflection_symmetry = _translation_report_has_label_field(report, :reflection) ||
        _translation_report_has_label_field(report, :reflection_parity)
    axis_rotation_symmetry = _translation_report_has_label_field(report, :axis_irrep)
    missing = Symbol[]
    report.sign_symmetry || push!(missing, :sign)
    reflection_symmetry || push!(missing, :reflection)
    axis_rotation_symmetry || push!(missing, :axis_rotation)

    return (
        translation_symmetry=true,
        sign_symmetry=report.sign_symmetry,
        reflection_symmetry=reflection_symmetry,
        axis_rotation_symmetry=axis_rotation_symmetry,
        axis_orbit_closed=report.axis_orbit_closed,
        axis_orbit_basis_size=report.axis_orbit_basis_size,
        axis_orbit_size_histogram=copy(report.axis_orbit_size_histogram),
        axis_reduction_ratio=report.axis_reduction_ratio,
        full_discrete_profile=isempty(missing),
        missing_discrete_symmetries=missing,
        real_moment_matrix=report.real_moment_matrix,
        momentum_sector_count=length(report.momentum_sectors),
        axis_rotation_equalities=_translation_report_has_zero_feature(
            report,
            :axis_rotation_equality,
        ),
        axis_rotation_quotient=report.axis_rotation_quotient,
        su2_moment_symmetry=_translation_report_has_feature(
            report,
            :moment_matrix;
            decomposition=:translation_su2,
        ) || _translation_report_has_feature(
            report,
            :moment_matrix;
            decomposition=:translation_su2_reflection,
        ),
        su2_reflection_moment_symmetry=_translation_report_has_feature(
            report,
            :moment_matrix;
            decomposition=:translation_su2_reflection,
        ),
        contiguous_rdm=_translation_report_has_feature(report, :contiguous_rdm),
        u1_rdm_symmetry=_translation_report_has_feature(
            report,
            :contiguous_rdm;
            decomposition=:u1,
        ),
        su2_rdm_symmetry=_translation_report_has_feature(
            report,
            :contiguous_rdm;
            decomposition=:su2,
        ),
        linear_state_opt=_translation_report_has_zero_feature(report, :linear_state_opt),
        psd_state_opt=_translation_report_has_feature(report, :psd_state_opt),
    )
end

const _TRANSLATION_SOLVE_BASE_SU2_BLOCKED_REASON =
    "This translation/SU(2) moment reducer is not yet solve-supported: " *
    "the public SDP solve is currently validated only for base SU(2) " *
    "moment blocks with Wigner zero/copy rows, SU(2) singlet-channel equality " *
    "rows, invariant scalar constraints, moment-equality rows, axis-rotation " *
    "equality rows, linear state-opt rows, PSD state-opt blocks, and " *
    "closed-support full/U(1)/SU(2) RDM blocks. " *
    "Use the translation path without " *
    "`su2_symmetry=true`, drop the extra SU(2) combination, or use the " *
    "contiguous RDM SU(2) reducer."

const _TRANSLATION_SOLVE_STRUCTURAL_TARGET_REASON =
    "Structural target only: this analytic payload reports block and " *
    "constraint-accounting estimates but does not construct or validate an " *
    "SDP solve. Build a translation/SU(2) moment report to check public " *
    "solve support."

const _TRANSLATION_STRUCTURAL_MODEL_SIZE_GATE_REASON =
    "Structural target only: no scalar-equality upper bound or dense-Schur " *
    "proxy has been constructed. Use a solver-probe/model-sizing run before " *
    "launching lowering or solve."

function _translation_label_is_base_su2_moment(label)
    label isa NamedTuple || return false
    haskey(label, :feature) && label.feature == :moment_matrix || return false
    haskey(label, :decomposition) || return false
    return label.decomposition in (:translation_su2, :translation_su2_reflection)
end

function _translation_label_is_validated_base_su2_block(label)
    _translation_label_is_base_su2_moment(label) && return true
    label isa NamedTuple || return false
    haskey(label, :feature) || return false
    label.feature == :scalar_inequality && return true
    label.feature == :psd_state_opt && return true
    return label.feature == :contiguous_rdm &&
        (
            !haskey(label, :decomposition) ||
            label.decomposition in (:u1, :su2)
        )
end

function _translation_zero_label_is_base_su2_wigner(label)
    label isa NamedTuple || return false
    haskey(label, :feature) && label.feature == :pauli_su2_base_zero || return false
    haskey(label, :decomposition) && label.decomposition == :translation_su2 || return false
    haskey(label, :reason) || return false
    return label.reason in (:spin_offblock, :magnetic_offdiag, :magnetic_copy)
end

function _translation_zero_label_is_su2_singlet_channel_equality(label)
    label isa NamedTuple || return false
    haskey(label, :feature) &&
        label.feature == :pauli_su2_translation_orbit_singlet_channel_equality ||
        return false
    return haskey(label, :decomposition) && label.decomposition == :translation_su2
end

function _translation_zero_label_is_validated_base_su2(label)
    _translation_zero_label_is_base_su2_wigner(label) && return true
    _translation_zero_label_is_su2_singlet_channel_equality(label) && return true
    label isa NamedTuple || return false
    haskey(label, :feature) || return false
    label.feature in (
        :scalar_equality,
        :moment_equality,
        :axis_rotation_equality,
        :linear_state_opt,
    ) && return true
    return label.feature == :contiguous_rdm_zero &&
        haskey(label, :decomposition) &&
        label.decomposition in (:u1, :su2) &&
        haskey(label, :reason) &&
        label.reason in (:magnetization_offblock, :spin_magnetic_offblock, :magnetic_copy)
end

_translation_label_feature(label) =
    label isa NamedTuple && haskey(label, :feature) ? label.feature : :unknown

function _translation_invalid_base_su2_block_features(report::TranslationInvariantReport)
    features = Symbol[]
    for label in report.block_labels
        _translation_label_is_validated_base_su2_block(label) && continue
        push!(features, _translation_label_feature(label))
    end
    return sort!(unique(features))
end

function _translation_base_su2_source_zero_histogram(
    report::TranslationInvariantReport,
)
    report.su2_moment_quotient || return report.zero_constraint_feature_histogram
    return vcat(
        report.zero_constraint_feature_histogram,
        report.su2_moment_eliminated_zero_feature_histogram,
    )
end

function _translation_invalid_base_su2_zero_features(report::TranslationInvariantReport)
    features = Symbol[]
    for pair in _translation_base_su2_source_zero_histogram(report)
        label = first(pair)
        _translation_zero_label_is_validated_base_su2(label) && continue
        push!(features, _translation_label_feature(label))
    end
    return sort!(unique(features))
end

function _translation_report_has_base_su2_wigner_rows(report::TranslationInvariantReport)
    histogram = _translation_base_su2_source_zero_histogram(report)
    has_reason(reason) = any(histogram) do pair
        label = first(pair)
        label isa NamedTuple || return false
        return haskey(label, :feature) && label.feature == :pauli_su2_base_zero &&
            haskey(label, :reason) && label.reason == reason
    end
    return has_reason(:spin_offblock) && has_reason(:magnetic_offdiag) &&
        has_reason(:magnetic_copy)
end

function _translation_report_is_validated_base_su2_solve(report::TranslationInvariantReport)
    isempty(report.block_labels) && return false
    all(_translation_label_is_validated_base_su2_block, report.block_labels) ||
        return false
    (report.su2_moment_quotient || _translation_report_has_base_su2_wigner_rows(report)) ||
        return false
    return all(
        pair -> _translation_zero_label_is_validated_base_su2(first(pair)),
        _translation_base_su2_source_zero_histogram(report),
    )
end

function _translation_su2_moment_quotient_validation_error(
    report::TranslationInvariantReport,
)
    report.su2_moment_quotient || return nothing
    raw_count = report.su2_moment_raw_count
    quotient_count = report.su2_moment_quotient_count
    0 < quotient_count <= raw_count || return (
        "invalid moment counts raw=$raw_count, quotient=$quotient_count"
    )
    quotient_count == report.linear_moment_count || return (
        "quotient count $quotient_count does not match linear moment count " *
        "$(report.linear_moment_count)"
    )
    expected_ratio = quotient_count / raw_count
    isfinite(report.su2_moment_quotient_reduction_ratio) &&
        isapprox(
            report.su2_moment_quotient_reduction_ratio,
            expected_ratio;
            atol=1e-12,
            rtol=1e-12,
        ) || return "invalid quotient reduction ratio"
    report.su2_moment_support_orbit_count > 0 || return "no support orbits"
    sum(last, report.su2_moment_singlet_channel_support_counts; init=0) ==
        quotient_count || return "singlet-channel support counts do not sum to the quotient count"
    residuals = (
        report.su2_moment_max_pivot_residual,
        report.su2_moment_max_invariant_residual,
        report.su2_moment_max_reconstruction_residual,
    )
    all(value -> isfinite(value) && 0 <= value <= 1e-8, residuals) || return (
        "quotient residuals are non-finite or exceed 1e-8: $(repr(residuals))"
    )
    isfinite(report.su2_moment_max_condition) &&
        report.su2_moment_max_condition > 0 || return "invalid coordinate-chart condition"
    report.su2_moment_eliminated_zero_row_count >= 0 || return (
        "negative eliminated zero-row count"
    )
    sum(last, report.su2_moment_eliminated_zero_feature_histogram; init=0) ==
        report.su2_moment_eliminated_zero_row_count || return (
            "eliminated zero-row histogram does not match its count"
        )
    return nothing
end

"""
    translation_solve_support(report)

Report whether a [`TranslationInvariantReport`](@ref) is supported by the
public solve wrapper.  Unvalidated reducer combinations return
`supported=false` with a human-readable `reason`.
"""
function translation_solve_support(report::TranslationInvariantReport)
    quotient_error = _translation_su2_moment_quotient_validation_error(report)
    if quotient_error !== nothing
        return (
            supported=false,
            reason="Invalid SU(2) moment quotient report: $quotient_error.",
            blocker=:su2_moment_quotient,
            unsupported_block_features=Symbol[],
            unsupported_zero_features=Symbol[],
        )
    end
    profile = translation_symmetry_profile(report)
    if profile.su2_moment_symmetry
        if _translation_report_is_validated_base_su2_solve(report)
            return (
                supported=true,
                reason="",
                blocker=nothing,
                unsupported_block_features=Symbol[],
                unsupported_zero_features=Symbol[],
            )
        end
        unsupported_block_features =
            _translation_invalid_base_su2_block_features(report)
        unsupported_zero_features =
            _translation_invalid_base_su2_zero_features(report)
        detail = isempty(unsupported_block_features) &&
            isempty(unsupported_zero_features) ? "" : " Unsupported block features: " *
            repr(unsupported_block_features) * "; unsupported zero-row features: " *
            repr(unsupported_zero_features) * "."
        return (
            supported=false,
            reason=_TRANSLATION_SOLVE_BASE_SU2_BLOCKED_REASON * detail,
            blocker=:base_translation_su2_moment_reducer,
            unsupported_block_features=unsupported_block_features,
            unsupported_zero_features=unsupported_zero_features,
        )
    end
    return (
        supported=true,
        reason="",
        blocker=nothing,
        unsupported_block_features=Symbol[],
        unsupported_zero_features=Symbol[],
    )
end

"""
    translation_report_metrics(report; scalar_bytes=sizeof(Float64))

Return deterministic structural metrics for a [`TranslationInvariantReport`](@ref):
block histograms, moment count, aggregate PSD dimensions, and dense/symmetric
scalar-storage byte estimates.  These byte counts are shape estimates, not solver
RSS measurements.
"""
function translation_report_metrics(
    report::TranslationInvariantReport;
    scalar_bytes::Integer=sizeof(Float64),
)
    bytes = Int(scalar_bytes)
    bytes > 0 || throw(ArgumentError("`scalar_bytes` must be positive; got $scalar_bytes."))

    psd_dense_entries = 0
    psd_symmetric_entries = 0
    for size in report.psd_block_sizes
        psd_dense_entries += size * size
        psd_symmetric_entries += report.real_moment_matrix ?
            size * (size + 1) ÷ 2 :
            size * size
    end
    profile = translation_symmetry_profile(report)
    solve_support = translation_solve_support(report)
    sos_dual_scalar_equalities = report.real_moment_matrix ?
        max(report.linear_moment_count - 1, 0) :
        max(2 * report.linear_moment_count - 1, 0)
    rdm_logical_sizes = _translation_report_feature_block_sizes(
        report,
        :contiguous_rdm;
        kind=:logical,
    )
    rdm_psd_sizes = _translation_report_feature_block_sizes(
        report,
        :contiguous_rdm;
        kind=:psd,
    )
    psd_state_logical_sizes = _translation_report_feature_block_sizes(
        report,
        :psd_state_opt;
        kind=:logical,
    )
    psd_state_psd_sizes = _translation_report_feature_block_sizes(
        report,
        :psd_state_opt;
        kind=:psd,
    )

    return (
        n_sites=report.n_sites,
        order=report.order,
        basis_size=report.basis_size,
        orbit_basis_size=report.orbit_basis_size,
        axis_orbit_closed=report.axis_orbit_closed,
        axis_orbit_basis_size=report.axis_orbit_basis_size,
        axis_orbit_size_histogram=copy(report.axis_orbit_size_histogram),
        axis_reduction_ratio=report.axis_reduction_ratio,
        translation_symmetry=profile.translation_symmetry,
        sign_symmetry=profile.sign_symmetry,
        reflection_symmetry=profile.reflection_symmetry,
        axis_rotation_symmetry=profile.axis_rotation_symmetry,
        axis_rotation_equalities=profile.axis_rotation_equalities,
        axis_rotation_quotient=profile.axis_rotation_quotient,
        full_discrete_profile=profile.full_discrete_profile,
        missing_discrete_symmetries=copy(profile.missing_discrete_symmetries),
        su2_moment_symmetry=profile.su2_moment_symmetry,
        su2_reflection_moment_symmetry=profile.su2_reflection_moment_symmetry,
        contiguous_rdm=profile.contiguous_rdm,
        u1_rdm_symmetry=profile.u1_rdm_symmetry,
        su2_rdm_symmetry=profile.su2_rdm_symmetry,
        linear_state_opt=profile.linear_state_opt,
        psd_state_opt=profile.psd_state_opt,
        momentum_sector_count=length(report.momentum_sectors),
        n_blocks=length(report.psd_block_sizes),
        moment_count=report.n_unique_moment_matrix_elements,
        base_moment_count=report.n_unique_moment_matrix_elements,
        linear_moment_count=report.linear_moment_count,
        su2_moment_quotient=report.su2_moment_quotient,
        su2_moment_raw_count=report.su2_moment_raw_count,
        su2_moment_quotient_count=report.su2_moment_quotient_count,
        su2_moment_quotient_reduction_ratio=
            report.su2_moment_quotient_reduction_ratio,
        su2_moment_support_orbit_count=report.su2_moment_support_orbit_count,
        su2_moment_singlet_channel_support_counts=
            copy(report.su2_moment_singlet_channel_support_counts),
        su2_moment_max_pivot_residual=report.su2_moment_max_pivot_residual,
        su2_moment_max_invariant_residual=report.su2_moment_max_invariant_residual,
        su2_moment_max_reconstruction_residual=
            report.su2_moment_max_reconstruction_residual,
        su2_moment_max_condition=report.su2_moment_max_condition,
        su2_moment_eliminated_zero_row_count=
            report.su2_moment_eliminated_zero_row_count,
        su2_moment_eliminated_zero_feature_histogram=
            copy(report.su2_moment_eliminated_zero_feature_histogram),
        zero_constraint_count=report.zero_constraint_count,
        estimated_sos_dual_scalar_equalities_upper_bound=sos_dual_scalar_equalities,
        estimated_sos_dual_dense_schur_bytes=bytes * sos_dual_scalar_equalities^2,
        zero_constraint_feature_histogram=copy(report.zero_constraint_feature_histogram),
        contiguous_rdm_block_count=length(rdm_psd_sizes),
        contiguous_rdm_logical_block_sizes=rdm_logical_sizes,
        contiguous_rdm_psd_block_sizes=rdm_psd_sizes,
        contiguous_rdm_zero_row_count=_translation_report_zero_feature_count(
            report,
            :contiguous_rdm_zero,
        ),
        scalar_equality_row_count=_translation_report_zero_feature_count(
            report,
            :scalar_equality,
        ),
        axis_rotation_equality_row_count=_translation_report_zero_feature_count(
            report,
            :axis_rotation_equality,
        ),
        axis_rotation_raw_moment_key_count=report.axis_rotation_raw_moment_key_count,
        axis_rotation_moment_class_count=report.axis_rotation_moment_class_count,
        axis_rotation_quotient_moment_key_count=report.axis_rotation_quotient_moment_key_count,
        axis_rotation_forced_zero_moment_class_count=report.axis_rotation_forced_zero_moment_class_count,
        axis_rotation_moment_quotient_reduction_ratio=report.axis_rotation_moment_quotient_reduction_ratio,
        moment_equality_row_count=_translation_report_zero_feature_count(
            report,
            :moment_equality,
        ),
        linear_state_opt_row_count=_translation_report_zero_feature_count(
            report,
            :linear_state_opt,
        ),
        su2_singlet_channel_equality_row_count=_translation_report_zero_feature_count(
            report,
            :pauli_su2_translation_orbit_singlet_channel_equality,
        ),
        su2_base_zero_row_count=_translation_report_zero_feature_count(
            report,
            :pauli_su2_base_zero,
        ),
        su2_base_spin_offblock_row_count=_translation_report_zero_feature_count(
            report,
            :pauli_su2_base_zero;
            reason=:spin_offblock,
        ),
        su2_base_magnetic_offdiag_row_count=_translation_report_zero_feature_count(
            report,
            :pauli_su2_base_zero;
            reason=:magnetic_offdiag,
        ),
        su2_base_magnetic_copy_row_count=_translation_report_zero_feature_count(
            report,
            :pauli_su2_base_zero;
            reason=:magnetic_copy,
        ),
        psd_state_opt_block_count=length(psd_state_psd_sizes),
        psd_state_opt_logical_block_sizes=psd_state_logical_sizes,
        psd_state_opt_psd_block_sizes=psd_state_psd_sizes,
        solve_supported=solve_support.supported,
        solve_blocker=solve_support.blocker,
        solve_blocker_reason=solve_support.reason,
        solve_unsupported_block_features=solve_support.unsupported_block_features,
        solve_unsupported_zero_features=solve_support.unsupported_zero_features,
        construction_time_ns=report.construction_time_ns,
        construction_time_seconds=report.construction_time_ns / 1e9,
        construction_stage_times_ns=copy(report.construction_stage_times_ns),
        construction_stage_time_seconds=Dict(
            key => value / 1e9 for (key, value) in report.construction_stage_times_ns
        ),
        logical_max_block=isempty(report.logical_block_sizes) ? 0 : maximum(report.logical_block_sizes),
        psd_max_block=isempty(report.psd_block_sizes) ? 0 : maximum(report.psd_block_sizes),
        logical_total_block_side=sum(report.logical_block_sizes; init=0),
        psd_total_block_side=sum(report.psd_block_sizes; init=0),
        psd_dense_entries=psd_dense_entries,
        psd_symmetric_entries=psd_symmetric_entries,
        psd_dense_bytes=psd_dense_entries * bytes,
        psd_symmetric_bytes=psd_symmetric_entries * bytes,
        block_coefficient_domains=copy(report.block_coefficient_domains),
        block_exact_coefficient_domains=copy(report.block_exact_coefficient_domains),
        block_coefficient_domain_histogram=_value_histogram_pairs(report.block_coefficient_domains),
        block_exact_coefficient_domain_histogram=_value_histogram_pairs(report.block_exact_coefficient_domains),
        product_cache_hits=report.product_cache_hits,
        product_cache_misses=report.product_cache_misses,
        product_cache_lookups=report.product_cache_hits + report.product_cache_misses,
        product_cache_entries=report.product_cache_entries,
        product_cache_hit_rate=report.product_cache_hit_rate,
        logical_block_histogram=translation_block_histogram(report; kind=:logical),
        psd_block_histogram=translation_block_histogram(report; kind=:psd),
        logical_block_feature_histogram=translation_block_feature_histogram(
            report; kind=:logical
        ),
        psd_block_feature_histogram=translation_block_feature_histogram(report; kind=:psd),
    )
end

function _pauli_translation_target_block_sizes(
    orbit_basis_size::Int,
    signature_counts::Dict{UInt8,Int},
    sectors::Vector{Int};
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
)
    nontrivial_size = orbit_basis_size - 1
    logical_sizes = Int[]

    if sign_symmetry
        signatures = sort!(collect(keys(signature_counts)))
        0x00 in signatures || push!(signatures, 0x00)
        sort!(signatures)
        for k in sectors
            for signature in signatures
                block_size = get(signature_counts, signature, 0)
                if k == 0 && signature == 0x00
                    block_size += 1
                end
                block_size > 0 && push!(logical_sizes, block_size)
            end
        end
    else
        for k in sectors
            block_size = k == 0 ? orbit_basis_size : nontrivial_size
            block_size > 0 && push!(logical_sizes, block_size)
        end
    end

    psd_sizes = real_moment_matrix ? [2 * size for size in logical_sizes] : copy(logical_sizes)
    return logical_sizes, psd_sizes
end

function _pauli_reflection_fixed_block_sizes(
    signature_counts::Dict{UInt8,Int},
    palindrome_even_counts::Dict{UInt8,Int},
    palindrome_odd_counts::Dict{UInt8,Int},
    k::Int,
    n::Int;
    sign_symmetry::Bool,
)
    total_count = sum(values(signature_counts); init=0)
    even_count = sum(values(palindrome_even_counts); init=0)
    odd_count = sum(values(palindrome_odd_counts); init=0)
    signatures = sign_symmetry ? sort!(collect(keys(signature_counts))) : UInt8[0x00]
    sign_symmetry && !(0x00 in signatures) && push!(signatures, 0x00)
    sort!(signatures)

    sizes = Int[]
    for signature in signatures
        nontrivial_count = sign_symmetry ? get(signature_counts, signature, 0) : total_count
        pal_even = sign_symmetry ? get(palindrome_even_counts, signature, 0) : even_count
        pal_odd = sign_symmetry ? get(palindrome_odd_counts, signature, 0) : odd_count
        pair_count = (nontrivial_count - pal_even - pal_odd) ÷ 2

        fixed_plus = if k == 0
            pal_even + pal_odd + (signature == 0x00 ? 1 : 0)
        else
            pal_even
        end
        fixed_minus = k == 0 ? 0 : pal_odd

        plus_size = fixed_plus + pair_count
        minus_size = fixed_minus + pair_count
        plus_size > 0 && push!(sizes, plus_size)
        minus_size > 0 && push!(sizes, minus_size)
    end
    return sizes
end

function _pauli_translation_reflection_target_block_sizes(
    orbit_basis_size::Int,
    signature_counts::Dict{UInt8,Int},
    sectors::Vector{Int},
    n::Int;
    order::Int,
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
)
    if !real_moment_matrix
        nonfixed = Int[k for k in sectors if !_translation_reflection_fixed_momentum(k, n)]
        isempty(nonfixed) || throw(ArgumentError(
            "`reflection_symmetry=true, real_moment_matrix=false` structural targets " *
            "are supported only for reflection-fixed momentum sectors; got non-fixed " *
            "sector(s) $nonfixed."
        ))
    end

    palindrome_even_counts, palindrome_odd_counts = _pauli_axis_palindrome_counts(
        order,
    )
    total_count = orbit_basis_size - 1
    signatures = sign_symmetry ? sort!(collect(keys(signature_counts))) : UInt8[0x00]
    sign_symmetry && !(0x00 in signatures) && push!(signatures, 0x00)
    sort!(signatures)

    logical_sizes = Int[]
    psd_sizes = Int[]
    for k in sectors
        if _translation_reflection_fixed_momentum(k, n)
            fixed_sizes = _pauli_reflection_fixed_block_sizes(
                signature_counts,
                palindrome_even_counts,
                palindrome_odd_counts,
                k,
                n;
                sign_symmetry,
            )
            append!(logical_sizes, fixed_sizes)
            append!(psd_sizes, real_moment_matrix ? 2 .* fixed_sizes : fixed_sizes)
        else
            for signature in signatures
                block_size = sign_symmetry ? get(signature_counts, signature, 0) : total_count
                block_size > 0 || continue
                push!(logical_sizes, block_size)
                push!(logical_sizes, block_size)
                push!(psd_sizes, block_size)
                push!(psd_sizes, block_size)
            end
        end
    end
    return logical_sizes, psd_sizes
end

function _pauli_axis_target_solver_size(
    logical_size::Integer,
    axis_irrep::Symbol;
    real_moment_matrix::Bool,
    realified_reflection::Bool=false,
)
    size = Int(logical_size)
    real_moment_matrix || return size
    realified_reflection && return size
    return axis_irrep == :trivial ? size : 2 * size
end

function _pauli_translation_axis_target_block_sizes(
    n::Int,
    d::Int,
    sectors::Vector{Int};
    reflection_symmetry::Bool,
    real_moment_matrix::Bool,
)
    T = select_uint_type(3, n)
    ops = _pauli_target_chain_ops(T, n)
    orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    nontrivial_reps = [mono for mono in orbit_reps if !isone(mono)]

    logical_sizes = Int[]
    psd_sizes = Int[]
    for k in sectors
        block_basis = iszero(k) ? orbit_reps : nontrivial_reps
        if reflection_symmetry
            if _translation_reflection_fixed_momentum(k, n)
                for reflected in _translation_reflection_adapted_blocks(
                    block_basis,
                    k,
                    n;
                    atol=1e-10,
                )
                    for adapted in _translation_reflection_axis_isotypic_adapted_blocks(
                        ops,
                        block_basis,
                        d,
                        reflected.reflection,
                        reflected.row_basis;
                        realified=false,
                        atol=1e-10,
                    )
                        logical_size = length(adapted.row_labels)
                        push!(logical_sizes, logical_size)
                        push!(
                            psd_sizes,
                            _pauli_axis_target_solver_size(
                                logical_size,
                                adapted.axis_irrep;
                                real_moment_matrix,
                            ),
                        )
                    end
                end
            else
                real_moment_matrix || throw(ArgumentError(
                    "`axis_rotation_symmetry=true, reflection_symmetry=true, real_moment_matrix=false` structural targets " *
                    "are supported only for reflection-fixed momentum sectors; got non-fixed sector $k."
                ))
                for reflected in _translation_real_reflection_adapted_blocks(
                    block_basis,
                    k,
                    n;
                    atol=1e-10,
                )
                    for adapted in _translation_reflection_axis_isotypic_adapted_blocks(
                        ops,
                        block_basis,
                        d,
                        reflected.reflection,
                        reflected.row_basis;
                        realified=true,
                        atol=1e-10,
                    )
                        logical_size = length(adapted.row_labels)
                        push!(logical_sizes, logical_size)
                        push!(
                            psd_sizes,
                            _pauli_axis_target_solver_size(
                                logical_size,
                                adapted.axis_irrep;
                                real_moment_matrix,
                                realified_reflection=true,
                            ),
                        )
                    end
                end
            end
        else
            for adapted in _translation_axis_isotypic_adapted_blocks(
                ops,
                block_basis,
                d;
                atol=1e-10,
            )
                logical_size = length(adapted.row_labels)
                push!(logical_sizes, logical_size)
                push!(
                    psd_sizes,
                    _pauli_axis_target_solver_size(
                        logical_size,
                        adapted.axis_irrep;
                        real_moment_matrix,
                    ),
                )
            end
        end
    end
    return logical_sizes, psd_sizes
end

@inline _triangular_entry_count(size::Int) = size * (size + 1) ÷ 2

function _pauli_translation_target_product_cache_stats(
    orbit_basis_size::Int,
    signature_counts::Dict{UInt8,Int},
    sectors::Vector{Int};
    sign_symmetry::Bool,
)
    nontrivial_size = orbit_basis_size - 1
    hits = 0
    misses = 0

    if sign_symmetry
        signatures = sort!(collect(keys(signature_counts)))
        0x00 in signatures || push!(signatures, 0x00)
        sort!(signatures)
        seen_nontrivial = Set{UInt8}()
        for k in sectors
            for signature in signatures
                nontrivial_count = get(signature_counts, signature, 0)
                identity_count = k == 0 && signature == 0x00 ? 1 : 0
                block_size = nontrivial_count + identity_count
                block_size > 0 || continue

                nontrivial_pairs = _triangular_entry_count(nontrivial_count)
                if nontrivial_count > 0
                    if signature in seen_nontrivial
                        hits += nontrivial_pairs
                    else
                        misses += nontrivial_pairs
                        push!(seen_nontrivial, signature)
                    end
                end
                identity_count == 1 && (misses += nontrivial_count + 1)
            end
        end
    else
        seen_nontrivial = false
        for k in sectors
            identity_count = k == 0 ? 1 : 0
            block_size = nontrivial_size + identity_count
            block_size > 0 || continue

            nontrivial_pairs = _triangular_entry_count(nontrivial_size)
            if nontrivial_size > 0
                if seen_nontrivial
                    hits += nontrivial_pairs
                else
                    misses += nontrivial_pairs
                    seen_nontrivial = true
                end
            end
            identity_count == 1 && (misses += nontrivial_size + 1)
        end
    end

    lookups = hits + misses
    return (
        product_cache_hits=hits,
        product_cache_misses=misses,
        product_cache_lookups=lookups,
        product_cache_entries=misses,
        product_cache_hit_rate=iszero(lookups) ? 0.0 : hits / lookups,
    )
end

function _pauli_target_domain_vector(domain::Union{Nothing,Symbol}, count::Integer)
    return Union{Nothing,Symbol}[domain for _ in 1:Int(count)]
end

function _pauli_translation_rdm_target_block_domains(
    rdm_ks::Vector{Int},
    decomposition::Symbol,
)
    coefficient_domains = Union{Nothing,Symbol}[]
    exact_coefficient_domains = Union{Nothing,Symbol}[]
    for rdm_k in rdm_ks
        if decomposition == :full
            push!(coefficient_domains, nothing)
            push!(exact_coefficient_domains, nothing)
        elseif decomposition == :u1
            for _ in _contiguous_rdm_u1_sectors(rdm_k)
                push!(coefficient_domains, nothing)
                push!(exact_coefficient_domains, nothing)
            end
        elseif decomposition == :su2
            for _ in pauli_su2_rdm_blocks(rdm_k)
                push!(coefficient_domains, :algebraic_float64)
                push!(exact_coefficient_domains, :sqrt_rational)
            end
        elseif decomposition == :qmbcertify
            for _ in pauli_qmbcertify_rdm_block_sizes(rdm_k)
                push!(coefficient_domains, nothing)
                push!(exact_coefficient_domains, nothing)
            end
        else
            throw(ArgumentError("Unsupported contiguous RDM decomposition $(repr(decomposition))."))
        end
    end
    return coefficient_domains, exact_coefficient_domains
end

function _pauli_translation_rdm_target_block_sizes(
    rdm_ks::Vector{Int},
    decomposition::Symbol,
    real_moment_matrix::Bool,
)
    logical_sizes = Int[]
    psd_sizes = Int[]
    for rdm_k in rdm_ks
        if decomposition == :full
            dim = 1 << rdm_k
            push!(logical_sizes, dim)
            push!(psd_sizes, real_moment_matrix ? 2 * dim : dim)
        elseif decomposition == :u1
            for sector in _contiguous_rdm_u1_sectors(rdm_k)
                dim = length(sector)
                push!(logical_sizes, dim)
                push!(psd_sizes, real_moment_matrix && dim > 1 ? 2 * dim : dim)
            end
        elseif decomposition == :su2
            for block in pauli_su2_rdm_blocks(rdm_k)
                push!(logical_sizes, block.multiplicity)
                push!(
                    psd_sizes,
                    real_moment_matrix ?
                        _pauli_su2_real_psd_block_size(block) :
                        block.multiplicity,
                )
            end
        elseif decomposition == :qmbcertify
            real_moment_matrix || throw(ArgumentError(
                "`contiguous_rdm_decomposition=:qmbcertify` structural targets are real PSD targets only."
            ))
            blocks = get(_QMBCERTIFY_RDM_REFERENCE_BLOCK_SIZES, rdm_k, nothing)
            blocks === nothing && throw(ArgumentError(
                "`contiguous_rdm_decomposition=:qmbcertify` has source-pinned targets only for k=8, 9, or 10; got k=$rdm_k."
            ))
            append!(logical_sizes, blocks)
            append!(psd_sizes, blocks)
        else
            throw(ArgumentError("Unsupported contiguous RDM decomposition $(repr(decomposition))."))
        end
    end
    return logical_sizes, psd_sizes
end

function _pauli_translation_u1_rdm_target_zero_rows(k::Int; real_moment_matrix::Bool)
    state_count = 1 << k
    active_entries = sum(length(sector)^2 for sector in _contiguous_rdm_u1_sectors(k); init=0)
    offblock_entries = state_count * state_count - active_entries
    return 2 * offblock_entries
end

function _pauli_translation_rdm_target_zero_row_count(
    rdm_ks::Vector{Int},
    decomposition::Symbol,
    real_moment_matrix::Bool,
)
    rows = 0
    for rdm_k in rdm_ks
        decomposition == :u1 || continue
        rows += _pauli_translation_u1_rdm_target_zero_rows(
            rdm_k;
            real_moment_matrix,
        )
    end
    return rows
end

const _QMBCERTIFY_RDM_REFERENCE_BLOCK_SIZES = Dict(
    8 => [72, 64, 56],
    9 => [136, 120],
    10 => [256, 272, 240],
)

function _pauli_qmbcertify_rdm_reference_metrics(k::Int; scalar_bytes::Int)
    blocks = get(_QMBCERTIFY_RDM_REFERENCE_BLOCK_SIZES, k, nothing)
    blocks === nothing && return nothing
    dense_entries = sum(size * size for size in blocks; init=0)
    symmetric_entries = sum(size * (size + 1) ÷ 2 for size in blocks; init=0)
    return (
        k=k,
        block_sizes=copy(blocks),
        n_blocks=length(blocks),
        max_block=maximum(blocks),
        total_block_side=sum(blocks; init=0),
        dense_entries=dense_entries,
        symmetric_entries=symmetric_entries,
        dense_bytes=dense_entries * scalar_bytes,
        symmetric_bytes=symmetric_entries * scalar_bytes,
        source=:qmbcertify,
        requires_construction=false,
    )
end

function _pauli_qmbcertify_rdm_reference_metrics(
    rdm_ks::Vector{Int};
    scalar_bytes::Int,
)
    metrics = Any[]
    for k in rdm_ks
        metric = _pauli_qmbcertify_rdm_reference_metrics(k; scalar_bytes)
        metric === nothing || push!(metrics, metric)
    end
    return metrics
end

function _pauli_qmbcertify_rdm_reference_summary(metrics::AbstractVector)
    block_sizes = Int[]
    dense_entries = 0
    symmetric_entries = 0
    dense_bytes = 0
    symmetric_bytes = 0
    requires_construction = false

    for metric in metrics
        append!(block_sizes, metric.block_sizes)
        dense_entries += metric.dense_entries
        symmetric_entries += metric.symmetric_entries
        dense_bytes += metric.dense_bytes
        symmetric_bytes += metric.symmetric_bytes
        requires_construction |= metric.requires_construction
    end

    return (
        block_sizes=block_sizes,
        n_blocks=length(block_sizes),
        max_block=isempty(block_sizes) ? 0 : maximum(block_sizes),
        total_block_side=sum(block_sizes; init=0),
        dense_entries=dense_entries,
        symmetric_entries=symmetric_entries,
        dense_bytes=dense_bytes,
        symmetric_bytes=symmetric_bytes,
        requires_construction=requires_construction,
    )
end

function _pauli_qmbcertify_base_reference_metrics(
    n::Int,
    d::Int,
    extra;
    three_type::Tuple{Int,Int}=(1, 1),
    scalar_bytes::Int,
    active::Bool=false,
)
    if isnothing(extra)
        return (
            extra=nothing,
            source_row_count=0,
            orbit_family_count=0,
            family_count_by_parity=Int[],
            block_records=NamedTuple[],
            block_sizes=Int[],
            block_histogram=Pair{Int,Int}[],
            n_blocks=0,
            max_block=0,
            total_block_side=0,
            dense_entries=0,
            symmetric_entries=0,
            dense_bytes=0,
            symmetric_bytes=0,
            support_nonzero_row_count=0,
            support_zero_row_count=0,
            support_diagonal_nonzero_row_count=0,
            support_offdiagonal_nonzero_row_count=0,
            support_unique_count=0,
            support_diagonal_unique_count=0,
            support_offdiagonal_unique_count=0,
            support_word_length_histogram=Pair{Int,Int}[],
            active=false,
            requires_construction=false,
        )
    end

    ext = Int(extra)
    ext >= 0 || throw(ArgumentError("`qmbcertify_base_extra` must be non-negative; got $extra."))
    all(>(0), three_type) || throw(ArgumentError(
        "`qmbcertify_base_three_type` entries must be positive; got $three_type."
    ))
    source_words_by_parity = [
        _qmbcertify_chain_source_basis_words(n, label, d; extra=ext, three_type)
        for label in (0, 1)
    ]
    family_counts = Int[length(source_words) ÷ n for source_words in source_words_by_parity]
    source_row_count = sum(length, source_words_by_parity; init=0)
    orbit_family_count = sum(family_counts; init=0)
    records = _qmbcertify_chain_base_block_records(family_counts, n)
    block_sizes = getproperty.(records, :block_size)
    dense_entries = sum(size * size for size in block_sizes; init=0)
    symmetric_entries = sum(size * (size + 1) ÷ 2 for size in block_sizes; init=0)
    support_metrics = _qmbcertify_chain_base_support_metrics(
        source_words_by_parity,
        family_counts,
        n,
    )
    return (
        extra=ext,
        source_row_count=source_row_count,
        orbit_family_count=orbit_family_count,
        family_count_by_parity=family_counts,
        block_records=records,
        block_sizes=block_sizes,
        block_histogram=_histogram_pairs(block_sizes),
        n_blocks=length(block_sizes),
        max_block=isempty(block_sizes) ? 0 : maximum(block_sizes),
        total_block_side=sum(block_sizes; init=0),
        dense_entries=dense_entries,
        symmetric_entries=symmetric_entries,
        dense_bytes=dense_entries * scalar_bytes,
        symmetric_bytes=symmetric_entries * scalar_bytes,
        support_nonzero_row_count=support_metrics.nonzero_row_count,
        support_zero_row_count=support_metrics.zero_row_count,
        support_diagonal_nonzero_row_count=support_metrics.diagonal_nonzero_row_count,
        support_offdiagonal_nonzero_row_count=support_metrics.offdiagonal_nonzero_row_count,
        support_unique_count=support_metrics.unique_count,
        support_diagonal_unique_count=support_metrics.diagonal_unique_count,
        support_offdiagonal_unique_count=support_metrics.offdiagonal_unique_count,
        support_word_length_histogram=support_metrics.word_length_histogram,
        active=active,
        requires_construction=true,
    )
end

function _pauli_qmbcertify_psd_state_opt_target_summary(
    n::Int,
    d::Int,
    extra,
    width::Int;
    three_type::Tuple{Int,Int}=(1, 1),
)
    iszero(width) && return (
        logical_sizes=Int[],
        psd_sizes=Int[],
        candidate_count=0,
    )
    ext = isnothing(extra) ? 0 : Int(extra)
    ops = _pauli_target_chain_ops(select_uint_type(3, n), n)
    qmb_basis = pauli_qmbcertify_chain_basis(ops, d; extra=ext, three_type)
    psd_rows = _qmbcertify_chain_psd_state_opt_rows(
        qmb_basis,
        _qmbcertify_chain_family_representatives(qmb_basis),
        width,
    )
    records = _qmbcertify_chain_psd_state_opt_block_records(
        psd_rows.family_count_by_parity,
        n,
    )
    block_sizes = getproperty.(records, :block_size)
    return (
        logical_sizes=copy(block_sizes),
        psd_sizes=block_sizes,
        candidate_count=sum(psd_rows.family_count_by_parity; init=0),
    )
end

function _qmbcertify_support_moment_set(
    ::Type{M},
    support_words::AbstractVector,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    moments = Set{M}()
    sizehint!(moments, length(support_words))
    for word in support_words
        mono = isempty(word) ?
            one(M) :
            NormalMonomial{PauliAlgebra,T}(T.(collect(word)))
        push!(moments, mono)
    end
    return moments
end

function _push_qmbcertify_rdm_target_moments!(
    moments::Set{M},
    n::Int,
    rdm_ks::AbstractVector{<:Integer},
    ::Type{K},
) where {M<:NormalMonomial{PauliAlgebra},K}
    isempty(rdm_ks) && return moments
    for rdm_k in rdm_ks
        row_blocks = _qmbcertify_rdm_block_templates(Int(rdm_k); ambient_sites=n)
        linear_blocks = _qmbcertify_chain_rdm_linear_blocks(
            n,
            Int(rdm_k),
            row_blocks,
            K,
            ComplexF64,
        )
        for entries in linear_blocks
            union!(moments, _linear_entries_moment_basis(M, entries))
        end
    end
    return moments
end

function _push_qmbcertify_psd_state_opt_target_moments!(
    moments::Set{M},
    qmb_basis,
    family_rows_by_parity,
    n::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    width::Int,
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    iszero(width) && return moments
    psd_rows = _qmbcertify_chain_psd_state_opt_rows(
        qmb_basis,
        family_rows_by_parity,
        width,
    )
    for block_basis in psd_rows.rows_by_parity
        isempty(block_basis) && continue
        translated = Dict(
            rep => [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
            for rep in block_basis
        )
        term_cache = _qmbcertify_chain_psd_state_opt_term_cache(
            block_basis,
            n,
            hamiltonian,
            translated,
        )
        union!(moments, _translation_psd_state_opt_required_moments(term_cache))
    end
    return moments
end

function _pauli_qmbcertify_linear_state_opt_target_summary(
    n::Int,
    d::Int,
    extra,
    width::Int;
    three_type::Tuple{Int,Int}=(1, 1),
    rdm_ks::AbstractVector{<:Integer}=Int[],
    psd_state_opt_width::Int=0,
)
    iszero(width) && return (row_count=0, candidate_count=0)
    ext = isnothing(extra) ? 0 : Int(extra)
    T = select_uint_type(3, n)
    M = NormalMonomial{PauliAlgebra,T}
    MP = Polynomial{PauliAlgebra,T,Float64}
    source_words_by_parity = [
        _qmbcertify_chain_source_basis_words(n, label, d; extra=ext, three_type)
        for label in (0, 1)
    ]
    family_counts = Int[length(source_words) ÷ n for source_words in source_words_by_parity]
    support_metrics = _qmbcertify_chain_base_support_metrics(
        source_words_by_parity,
        family_counts,
        n,
    )
    moment_set = _qmbcertify_support_moment_set(M, support_metrics.words)
    ops = _pauli_target_chain_ops(T, n)
    hamiltonian = heisenberg_chain_hamiltonian(ops)
    if !isempty(rdm_ks) || psd_state_opt_width > 0
        K = typeof(symmetric_canon(expval(one(M))))
        _push_qmbcertify_rdm_target_moments!(moment_set, n, rdm_ks, K)
        if psd_state_opt_width > 0
            qmb_basis = pauli_qmbcertify_chain_basis(
                ops,
                d;
                extra=ext,
                three_type,
            )
            _push_qmbcertify_psd_state_opt_target_moments!(
                moment_set,
                qmb_basis,
                _qmbcertify_chain_family_representatives(qmb_basis),
                n,
                hamiltonian,
                psd_state_opt_width,
            )
        end
    end
    tests = _qmbcertify_linear_state_opt_tests(M, n, width; sign_symmetry=true)
    seen_rows = Set{MP}()
    for test_mono in tests
        row = im * (hamiltonian * test_mono - test_mono * hamiltonian)
        iszero(row) && continue
        reduced = _qmbcertify_reduce_chain_polynomial(row, n, MP)
        iszero(reduced) && continue
        all(mono -> mono in moment_set, monomials(reduced)) || continue
        push!(seen_rows, _qmbcertify_linear_state_opt_row_key(reduced))
    end
    return (row_count=length(seen_rows), candidate_count=length(tests))
end

function _pauli_state_opt_test_count(test_width::Int; sign_symmetry::Bool)
    test_width <= 0 && return 0
    if sign_symmetry
        return get(_pauli_axis_signature_counts(test_width), 0x00, 0)
    end
    return sum(3^width for width in 1:test_width; init=0)
end

function _pauli_translation_linear_state_opt_target_row_count(
    n::Int,
    test_width::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
)
    test_width <= 0 && return 0
    T = select_uint_type(3, n)
    return _pauli_translation_linear_state_opt_target_row_count(
        T,
        n,
        test_width;
        sign_symmetry,
        real_moment_matrix,
    )
end

function _pauli_translation_linear_state_opt_target_row_count(
    ::Type{T},
    n::Int,
    test_width::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
) where {T<:Unsigned}
    M = NormalMonomial{PauliAlgebra,T}
    MP = real_moment_matrix ?
        Polynomial{PauliAlgebra,T,Float64} :
        Polynomial{PauliAlgebra,T,ComplexF64}
    ops = _pauli_target_chain_ops(T, n)
    hamiltonian = heisenberg_chain_hamiltonian(ops)
    rows = 0
    for test_mono in _contiguous_state_opt_tests(M, n, test_width; sign_symmetry)
        row = im * (hamiltonian * test_mono - test_mono * hamiltonian)
        iszero(row) && continue
        reduced = _translation_reduced_constraint_polynomial(
            row,
            n,
            MP;
            context="Linear state-opt structural target row",
        )
        iszero(reduced) || (rows += 1)
    end
    return rows
end

function _pauli_translation_moment_eq_h2_target_row_count(
    n::Int,
    d::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
)
    T = select_uint_type(3, n)
    return _pauli_translation_moment_eq_h2_target_row_count(
        T,
        n,
        d;
        sign_symmetry,
        real_moment_matrix,
        momenta,
    )
end

function _moment_eq_h2_target_closure_error(
    n::Int,
    d::Int,
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
)
    return ArgumentError(
        "`moment_eq_h2=true` structural target is not closed for n_sites=$n, " *
        "order=$d, sign_symmetry=$sign_symmetry, real_moment_matrix=$real_moment_matrix. " *
        "Some translated H^2 row multipliers require moments or coefficient domains " *
        "outside the current contiguous translation moment representation; use the " *
        "order-2 H^2 smoke target or construct the relaxation explicitly."
    )
end

function _pauli_translation_moment_eq_h2_target_row_count(
    ::Type{T},
    n::Int,
    d::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
) where {T<:Unsigned}
    M = NormalMonomial{PauliAlgebra,T}
    MP = real_moment_matrix ?
        Polynomial{PauliAlgebra,T,Float64} :
        Polynomial{PauliAlgebra,T,ComplexF64}
    ops = _pauli_target_chain_ops(T, n)
    orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    moment_basis = _pauli_translation_target_base_moment_basis(
        orbit_reps,
        n;
        sign_symmetry,
        real_moment_matrix,
        momenta,
    )
    hamiltonian = heisenberg_chain_hamiltonian(ops)
    moment_eq = hamiltonian * hamiltonian
    row_bases, row_degrees = _translation_moment_eq_row_bases(
        orbit_reps;
        sign_symmetry,
    )
    rows = _truncate_moment_eq_row_bases(row_bases, row_degrees, moment_eq)
    isempty(rows) && return 0

    one_mono = one(M)
    buf = T[]
    row_count = 0
    for row_mono in rows
        terms = Tuple{ComplexF64,M}[]
        sizehint!(terms, length(row_mono) * length(moment_eq.terms))
        for (c_row, row_word) in row_mono
            conj_row = _conj_coef(PauliAlgebra, c_row)
            for (coef, mono) in moment_eq.terms
                _push_scaled_buffered_terms!(
                    terms,
                    conj_row * coef,
                    PauliAlgebra,
                    simplify!(PauliAlgebra, _neat_dot3!(buf, row_word, mono, one_mono)),
                    T,
                    ComplexF64,
                )
            end
        end
        try
            poly = _polynomial_from_owned_terms!(terms)
            iszero(poly) && continue
            reduced = _translation_reduced_constraint_polynomial(
                poly,
                n,
                MP;
                context="Moment equality structural target",
            )
            iszero(reduced) && continue
            _check_polynomial_moments_covered(
                reduced,
                moment_basis,
                "Moment equality structural target",
            )
            reduced
        catch err
            err isa ArgumentError || rethrow()
            throw(_moment_eq_h2_target_closure_error(
                n,
                d,
                sign_symmetry,
                real_moment_matrix,
            ))
        end
        row_count += 1
    end
    return row_count
end

function _pauli_psd_state_opt_target_block_sizes(
    test_width::Int,
    sectors::Vector{Int};
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    full_realification::Bool=false,
)
    block_size = _pauli_state_opt_test_count(test_width; sign_symmetry)
    block_size == 0 && return Int[], Int[]
    logical_sizes = fill(block_size, length(sectors))
    psd_sizes = real_moment_matrix ?
        [full_realification || !iszero(sector) ? 2 * block_size : block_size for sector in sectors] :
        copy(logical_sizes)
    return logical_sizes, psd_sizes
end

function _pauli_translation_psd_state_opt_target_support_summary(
    ::Type{T},
    n::Int,
    d::Int,
    h::Int,
    test_width::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
) where {T<:Unsigned}
    M = NormalMonomial{PauliAlgebra,T}
    test_width == 0 && return (
        checked=true,
        base_support_closed=true,
        missing_moment_count=0,
        missing_moment_sample=String[],
    )
    h == 2 || return (
        checked=false,
        base_support_closed=false,
        missing_moment_count=0,
        missing_moment_sample=String[],
    )

    ops = _pauli_target_chain_ops(T, n)
    hamiltonian = heisenberg_chain_hamiltonian(ops)
    orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    moment_basis = _pauli_translation_target_base_moment_basis(
        orbit_reps,
        n;
        sign_symmetry,
        real_moment_matrix,
        momenta,
    )
    moment_set = Set(moment_basis)
    required = _psd_state_opt_required_moment_set(
        M,
        n,
        hamiltonian,
        test_width;
        sign_symmetry,
    )
    missing = M[moment for moment in required if !(moment in moment_set)]
    sort!(missing)
    return (
        checked=true,
        base_support_closed=isempty(missing),
        missing_moment_count=length(missing),
        missing_moment_sample=String[sprint(show, mono) for mono in Iterators.take(missing, 5)],
    )
end

function _pauli_target_chain_ops(::Type{T}, n::Int) where {T<:Unsigned}
    M = NormalMonomial{PauliAlgebra,T}
    make_ops(pauli_type) = M[
        M(T[convert(T, _pauli_index(site, pauli_type))]) for site in 1:n
    ]
    return (
        make_ops(_PAULI_X_TYPE),
        make_ops(_PAULI_Y_TYPE),
        make_ops(_PAULI_Z_TYPE),
    )
end

function _push_translation_target_moments!(
    moments::Set{M},
    terms::AbstractVector{<:Tuple{<:Number,M}},
    ::Type{C},
) where {M<:NormalMonomial{PauliAlgebra},C<:Number}
    accum = Dict{M,C}()
    for (coef, mono) in terms
        converted = convert(C, coef)
        iszero(converted) && continue
        accum[mono] = get(accum, mono, zero(C)) + converted
    end
    for (mono, coef) in accum
        iszero(coef) || push!(moments, mono)
    end
    return moments
end

function _pauli_translation_target_base_moment_basis(
    orbit_reps::Vector{M},
    n::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    nontrivial_reps = M[mono for mono in orbit_reps if !isone(mono)]
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    translated = _translation_orbit_representative_translates(orbit_reps, n)
    rep_cache = Dict{M,M}()
    product_cache = _TranslationProductCache(
        Dict{Tuple{M,M},Vector{Tuple{Int,ComplexF64,M}}}(),
        0,
        0,
    )
    moments = Set{M}()

    for k in sectors
        sector_basis = iszero(k) ? orbit_reps : nontrivial_reps
        blocks = sign_symmetry ? _pauli_signature_blocks(sector_basis) : [(:all, sector_basis)]
        for (_, block_basis) in blocks
            isempty(block_basis) && continue
            for j in eachindex(block_basis), i in 1:j
                terms = _translation_momentum_entry_terms(
                    block_basis[i],
                    block_basis[j],
                    k,
                    n,
                    translated,
                    rep_cache,
                    ComplexF64,
                    product_cache,
                )
                if i == j
                    _push_translation_target_moments!(
                        moments,
                        _translation_hermitian_diagonal_terms(ComplexF64, terms),
                        ComplexF64,
                    )
                else
                    _push_translation_target_moments!(moments, terms, ComplexF64)
                    _push_translation_target_moments!(
                        moments,
                        _translation_adjoint_entry_terms(ComplexF64, terms),
                        ComplexF64,
                    )
                end
            end
        end
    end

    return sorted_unique!(collect(moments))
end

function _pauli_translation_target_axis_rotation_equality_row_count(
    n::Int,
    d::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
)
    T = select_uint_type(3, n)
    return _pauli_translation_target_axis_rotation_equality_row_count(
        T,
        n,
        d;
        sign_symmetry,
        real_moment_matrix,
        momenta,
    )
end

function _pauli_translation_target_axis_rotation_equality_row_count(
    ::Type{T},
    n::Int,
    d::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
) where {T<:Unsigned}
    ops = _pauli_target_chain_ops(T, n)
    orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    moment_basis = _pauli_translation_target_base_moment_basis(
        orbit_reps,
        n;
        sign_symmetry,
        real_moment_matrix,
        momenta,
    )
    return length(_translation_axis_rotation_equality_rows(ops, n, moment_basis))
end

function _axis_quotient_find!(
    parent::Vector{Int},
    sign_to_parent::Vector{Int},
    idx::Int,
)
    p = parent[idx]
    p == idx && return idx, 1
    root, root_sign = _axis_quotient_find!(parent, sign_to_parent, p)
    sign_to_parent[idx] *= root_sign
    parent[idx] = root
    return root, sign_to_parent[idx]
end

function _axis_quotient_union!(
    parent::Vector{Int},
    sign_to_parent::Vector{Int},
    sizes::Vector{Int},
    forced_zero::AbstractVector{Bool},
    source_idx::Int,
    target_idx::Int,
    coefficient::Integer,
)
    coef = Int(coefficient)
    abs(coef) == 1 || throw(ArgumentError(
        "Axis-rotation quotient statistics only support ±1 equality coefficients; got $coefficient."
    ))

    source_root, source_sign = _axis_quotient_find!(parent, sign_to_parent, source_idx)
    target_root, target_sign = _axis_quotient_find!(parent, sign_to_parent, target_idx)
    if source_root == target_root
        if source_sign != coef * target_sign
            forced_zero[source_root] = true
        end
        return nothing
    end

    link_sign = coef * target_sign * source_sign
    if sizes[source_root] < sizes[target_root]
        parent[source_root] = target_root
        sign_to_parent[source_root] = link_sign
        sizes[target_root] += sizes[source_root]
        forced_zero[target_root] |= forced_zero[source_root]
    else
        parent[target_root] = source_root
        sign_to_parent[target_root] = link_sign
        sizes[source_root] += sizes[target_root]
        forced_zero[source_root] |= forced_zero[target_root]
    end
    return nothing
end

function _axis_moment_key_token(mono::NormalMonomial{PauliAlgebra})
    return _translation_registration_key_token(symmetric_canon(expval(mono)))
end

function _axis_moment_key(::Type{K}, mono::NormalMonomial{PauliAlgebra}) where {K}
    return _owned_moment_key(K, symmetric_canon(expval(mono)))
end

function _translation_axis_rotation_quotient_stats(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    n_sites::Integer,
    moment_basis::Vector{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    token_to_index = Dict{Any,Int}()
    for mono in moment_basis
        token = _axis_moment_key_token(mono)
        haskey(token_to_index, token) && continue
        token_to_index[token] = length(token_to_index) + 1
    end

    n_keys = length(token_to_index)
    parent = collect(1:n_keys)
    sign_to_parent = ones(Int, n_keys)
    sizes = ones(Int, n_keys)
    forced_zero = falses(n_keys)

    rows = _translation_axis_rotation_equality_rows(ops, n_sites, moment_basis)
    for row in rows
        source_idx = token_to_index[_axis_moment_key_token(row.source)]
        target_idx = token_to_index[_axis_moment_key_token(row.target)]
        _axis_quotient_union!(
            parent,
            sign_to_parent,
            sizes,
            forced_zero,
            source_idx,
            target_idx,
            row.coefficient,
        )
    end

    roots = Set{Int}()
    forced_roots = Set{Int}()
    for idx in 1:n_keys
        root, _ = _axis_quotient_find!(parent, sign_to_parent, idx)
        push!(roots, root)
        forced_zero[root] && push!(forced_roots, root)
    end
    component_count = length(roots)
    forced_zero_class_count = length(forced_roots)
    quotient_count = component_count - forced_zero_class_count
    return (
        equality_row_count=length(rows),
        raw_moment_key_count=n_keys,
        moment_class_count=component_count,
        quotient_moment_key_count=quotient_count,
        forced_zero_moment_class_count=forced_zero_class_count,
        reduction_ratio=quotient_count == 0 ? Inf : n_keys / quotient_count,
    )
end

function _translation_axis_rotation_quotient_map(
    ::Type{K},
    ::Type{M},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    n_sites::Integer,
    moment_basis::Vector{M},
) where {K,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    token_to_index = Dict{Any,Int}()
    token_to_key = Dict{Any,K}()
    token_to_mono = Dict{Any,M}()
    for mono in moment_basis
        key = _axis_moment_key(K, mono)
        token = _translation_registration_key_token(key)
        haskey(token_to_index, token) && continue
        token_to_index[token] = length(token_to_index) + 1
        token_to_key[token] = key
        token_to_mono[token] = _pauli_monomial_from_moment_key(M, key)
    end

    n_keys = length(token_to_index)
    parent = collect(1:n_keys)
    sign_to_parent = ones(Int, n_keys)
    sizes = ones(Int, n_keys)
    forced_zero = falses(n_keys)

    rows = _translation_axis_rotation_equality_rows(ops, n_sites, moment_basis)
    for row in rows
        source_idx = token_to_index[_axis_moment_key_token(row.source)]
        target_idx = token_to_index[_axis_moment_key_token(row.target)]
        _axis_quotient_union!(
            parent,
            sign_to_parent,
            sizes,
            forced_zero,
            source_idx,
            target_idx,
            row.coefficient,
        )
    end

    index_to_token = Vector{Any}(undef, n_keys)
    for (token, idx) in token_to_index
        index_to_token[idx] = token
    end

    members_by_root = Dict{Int,Vector{Any}}()
    sign_by_token = Dict{Any,Int}()
    for idx in 1:n_keys
        root, sign = _axis_quotient_find!(parent, sign_to_parent, idx)
        token = index_to_token[idx]
        push!(get!(members_by_root, root, Any[]), token)
        sign_by_token[token] = sign
    end

    mapping = Dict{Any,Any}()
    key_mapping = Dict{K,Any}()
    representative_pairs = Pair{K,M}[]
    forced_zero_class_count = 0
    for (root, tokens) in members_by_root
        if forced_zero[root]
            forced_zero_class_count += 1
            for token in tokens
                info = (
                    forced_zero=true,
                    sign=0,
                    key=token_to_key[token],
                    mono=token_to_mono[token],
                )
                mapping[token] = info
                key_mapping[token_to_key[token]] = info
            end
            continue
        end

        rep_token = first(tokens)
        for token in Iterators.drop(tokens, 1)
            if key_lt(token_to_key[token], token_to_key[rep_token])
                rep_token = token
            end
        end
        rep_sign = sign_by_token[rep_token]
        rep_key = token_to_key[rep_token]
        rep_mono = token_to_mono[rep_token]
        push!(representative_pairs, rep_key => rep_mono)
        for token in tokens
            info = (
                forced_zero=false,
                sign=sign_by_token[token] * rep_sign,
                key=rep_key,
                mono=rep_mono,
            )
            mapping[token] = info
            key_mapping[token_to_key[token]] = info
        end
    end

    component_count = length(members_by_root)
    quotient_count = component_count - forced_zero_class_count
    return (
        map=mapping,
        key_map=key_mapping,
        representative_pairs=representative_pairs,
        equality_row_count=length(rows),
        raw_moment_key_count=n_keys,
        moment_class_count=component_count,
        quotient_moment_key_count=quotient_count,
        forced_zero_moment_class_count=forced_zero_class_count,
        reduction_ratio=quotient_count == 0 ? Inf : n_keys / quotient_count,
    )
end

function _axis_quotient_representative_moment_map(
    ::Type{K},
    ::Type{M},
    quotient,
) where {K,M}
    key_to_monomial = Dict{K,M}()
    sizehint!(key_to_monomial, length(quotient.representative_pairs))
    for (key, mono) in quotient.representative_pairs
        key_to_monomial[key] = mono
    end
    return key_to_monomial
end

function _axis_quotient_info(quotient, key)
    info = get(quotient.key_map, key, nothing)
    info === nothing && throw(ArgumentError(
        "Axis-rotation quotient map is missing moment key $(repr(key))."
    ))
    return info
end

@inline function _axis_quotient_rewrite_key!(
    key_to_monomial::Dict{K,M},
    ::Type{K},
    info,
    ::Val{false},
) where {K,M}
    return _register_builder_moment_key!(
        key_to_monomial,
        K,
        info.key,
        info.mono,
    )
end

@inline function _axis_quotient_rewrite_key!(
    ::Dict{K,M},
    ::Type{K},
    info,
    ::Val{true},
) where {K,M}
    return info.key
end

function _axis_quotient_rewrite_form!(
    key_to_monomial::Dict{K,M},
    ::Type{K},
    ::Type{C},
    form::LinearMomentForm{K,FC},
    quotient,
    representatives_registered::Val=Val(false),
) where {K,C,FC,M}
    if length(form) > 8
        coefficients = Dict{K,C}()
        for (key, coef) in form
            info = _axis_quotient_info(quotient, key)
            info.forced_zero && continue
            stored = _axis_quotient_rewrite_key!(
                key_to_monomial,
                K,
                info,
                representatives_registered,
            )
            value = convert(C, coef) * convert(C, info.sign)
            iszero(value) && continue
            if haskey(coefficients, stored)
                updated = coefficients[stored] + value
                if iszero(updated)
                    delete!(coefficients, stored)
                else
                    coefficients[stored] = updated
                end
            else
                coefficients[stored] = value
            end
        end
        pairs = Pair{K,C}[]
        sizehint!(pairs, length(coefficients))
        for (key, coef) in coefficients
            push!(pairs, key => coef)
        end
        return _linear_moment_form_from_owned_pairs!(pairs)
    end

    pairs = Pair{K,C}[]
    sizehint!(pairs, length(form))
    for (key, coef) in form
        info = _axis_quotient_info(quotient, key)
        info.forced_zero && continue
        stored = _axis_quotient_rewrite_key!(
            key_to_monomial,
            K,
            info,
            representatives_registered,
        )
        push!(pairs, stored => convert(C, coef) * convert(C, info.sign))
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

struct _LinearMomentFormCacheKey{K,C}
    terms::Vector{Pair{K,C}}
    hash::UInt
end

function _linear_moment_form_cache_key(form::LinearMomentForm{K,C}) where {K,C}
    h = hash(length(form), UInt(0x35b840f4c91d3021))
    for (key, coef) in form
        h = hash(coef, hash(key, h))
    end
    return _LinearMomentFormCacheKey{K,C}(form.terms, h)
end

Base.hash(key::_LinearMomentFormCacheKey, h::UInt) = hash(key.hash, h)

function Base.isequal(
    lhs::_LinearMomentFormCacheKey{K,C},
    rhs::_LinearMomentFormCacheKey{K,C},
) where {K,C}
    length(lhs.terms) == length(rhs.terms) || return false
    for idx in eachindex(lhs.terms)
        lterm = lhs.terms[idx]
        rterm = rhs.terms[idx]
        key_isequal(lterm.first, rterm.first) || return false
        isequal(lterm.second, rterm.second) || return false
    end
    return true
end

function _axis_quotient_rewrite_pairs!(
    key_to_monomial::Dict{K,M},
    ::Type{K},
    ::Type{C},
    pairs,
    quotient,
    representatives_registered::Val=Val(false),
) where {K,C,M}
    return _axis_quotient_rewrite_form!(
        key_to_monomial,
        K,
        C,
        _linear_moment_form_from_owned_pairs!(_owned_builder_pairs(K, C, pairs)),
        quotient,
        representatives_registered,
    ).terms
end

function _apply_translation_axis_rotation_quotient!(
    builder::MomentLinearBuilder{K,C,M},
    quotient,
) where {K,C,M}
    _check_not_finalized!(builder)

    key_to_monomial = _axis_quotient_representative_moment_map(K, M, quotient)
    identity_info = _axis_quotient_info(quotient, builder.identity)
    identity_info.forced_zero && throw(ArgumentError(
        "Axis-rotation quotient forced the identity moment to zero."
    ))
    builder.identity = _owned_moment_key(K, identity_info.key)

    builder.objective_terms = _axis_quotient_rewrite_pairs!(
        key_to_monomial,
        K,
        C,
        builder.objective_terms,
        quotient,
        Val(true),
    )

    rewritten_blocks = PSDBlockLin{K,C,M}[]
    sizehint!(rewritten_blocks, length(builder.psd_blocks_lin))
    for block in builder.psd_blocks_lin
        entries = Matrix{LinearMomentForm{K,C}}(undef, size(block.entries))
        for idx in eachindex(block.entries)
            entries[idx] = _axis_quotient_rewrite_form!(
                key_to_monomial,
                K,
                C,
                block.entries[idx],
                quotient,
                Val(true),
            )
        end
        push!(rewritten_blocks, PSDBlockLin{K,C,M}(block.size, entries, block.meta))
    end
    builder.psd_blocks_lin = rewritten_blocks

    rewritten_zero_constraints = ScalarLinearConstraint{K,C}[]
    for zc in builder.zero_constraints
        form = _axis_quotient_rewrite_form!(
            key_to_monomial,
            K,
            C,
            zc.form,
            quotient,
            Val(true),
        )
        isempty(form) && continue
        push!(
            rewritten_zero_constraints,
            ScalarLinearConstraint(
                form,
                zc.kind,
                zc.origin;
                trusted_self_adjoint=zc.trusted_self_adjoint,
            ),
        )
    end
    builder.zero_constraints = rewritten_zero_constraints
    builder.key_to_monomial = key_to_monomial
    return nothing
end

function _pauli_translation_target_axis_rotation_quotient_stats(
    ::Type{T},
    n::Int,
    d::Int;
    sign_symmetry::Bool,
    real_moment_matrix::Bool,
    momenta,
) where {T<:Unsigned}
    ops = _pauli_target_chain_ops(T, n)
    orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    moment_basis = _pauli_translation_target_base_moment_basis(
        orbit_reps,
        n;
        sign_symmetry,
        real_moment_matrix,
        momenta,
    )
    return _translation_axis_rotation_quotient_stats(ops, n, moment_basis)
end

function _pauli_target_block_feature_histogram(
    base_sizes::Vector{Int},
    psd_state_sizes::Vector{Int},
    rdm_sizes::Vector{Int};
    base_decomposition::Symbol,
    rdm_decomposition::Symbol,
)
    counts = Dict{Any,Int}()
    rdm_decomp = rdm_decomposition == :full ? nothing : rdm_decomposition
    for (feature, decomposition, sizes) in (
        (:moment_matrix, base_decomposition, base_sizes),
        (:psd_state_opt, nothing, psd_state_sizes),
        (:contiguous_rdm, rdm_decomp, rdm_sizes),
    )
        for size in sizes
            key = (feature=feature, decomposition=decomposition, size=size)
            counts[key] = get(counts, key, 0) + 1
        end
    end
    return sort!(collect(counts); by=pair -> repr(pair.first))
end

function _pauli_target_zero_feature_histogram(
    axis_rotation_equality_row_count::Int,
    rdm_zero_row_count::Int,
    rdm_decomposition::Symbol,
    linear_state_opt_row_count::Int,
    moment_equality_row_count::Int,
)
    counts = Dict{Any,Int}()
    if axis_rotation_equality_row_count > 0
        key = (feature=:axis_rotation_equality, decomposition=nothing, reason=nothing)
        counts[key] = get(counts, key, 0) + axis_rotation_equality_row_count
    end
    if rdm_zero_row_count > 0 && rdm_decomposition == :u1
        key = (
            feature=:contiguous_rdm_zero,
            decomposition=:u1,
            reason=:magnetization_offblock,
        )
        counts[key] = get(counts, key, 0) + rdm_zero_row_count
    end
    if linear_state_opt_row_count > 0
        key = (feature=:linear_state_opt, decomposition=nothing, reason=nothing)
        counts[key] = get(counts, key, 0) + linear_state_opt_row_count
    end
    if moment_equality_row_count > 0
        key = (feature=:moment_equality, decomposition=nothing, reason=nothing)
        counts[key] = get(counts, key, 0) + moment_equality_row_count
    end
    return sort!(collect(counts); by=pair -> repr(pair.first))
end

"""
    pauli_translation_structural_targets(n_sites, order; sign_symmetry=true, real_moment_matrix=true, momenta=nothing, reflection_symmetry=false, hamiltonian_width=2, contiguous_rdm_k=nothing, contiguous_rdm_decomposition=:full, contiguous_rdm_support=:closed, u1_symmetry=false, su2_symmetry=false, base_su2_extend_rdm=false, axis_rotation_symmetry=false, axis_rotation_equalities=false, axis_rotation_quotient=false, moment_eq_h2=false, linear_state_opt_width=nothing, psd_state_opt_width=nothing, qmbcertify_base_construct=false, qmbcertify_base_extra=nothing, qmbcertify_base_three_type=(1, 1), scalar_bytes=sizeof(Float64))

Return analytic structural targets for the sparse periodic Pauli-chain
translation backend without constructing the site-space basis or PSD blocks.
The formula is valid for the contiguous basis in the no-short-orbit regime
`n_sites > 2order`.

The returned named tuple mirrors the construction-facing metrics that can be
known from `(n_sites, order)`: site-space and translation-orbit basis sizes,
momentum-sector count, block-size histograms, product-cache accounting, degree
closure limits, and the global Pauli-axis orbit diagnostic.
With `reflection_symmetry=true`, the target uses the no-short-orbit reversal
formula matching the realified reflection adaptation used by the constructor.
When RDM or PSD state-opt options are provided, their PSD block shapes are added
to the aggregate block metrics.  PSD state-opt has two support notions: the
degree bound `2a + h <= 2d`, and actual coverage by the base contiguous moment
PSD support.  Targets reject uncovered PSD-state profiles unless an extending
RDM support policy is present.  For the default nearest-neighbor Heisenberg
Hamiltonian width, linear state-opt also reports the exact nonzero row count;
for other Hamiltonian widths, the candidate count remains available but the
add-on zero-row count is marked incomplete.
With `axis_rotation_symmetry=true`, the target computes the finite global
axis-isotypic base PSD block shapes from the small translation-orbit
representative action and includes the required axis equality rows.
With `axis_rotation_quotient=true`, the target follows the direct-linear
finite-axis quotient storage convention for add-on blocks; real PSD-state-opt
blocks are fully realified, including the zero-momentum sector.
With `base_su2_extend_rdm=true`, the base moment-matrix target uses the
direct-linear translation/SU(2) reduction and appends extending SU(2) RDM
blocks.  This mode is shape/accounting-only and can account for requested
finite-axis equality rows and Heisenberg `H^2` moment-equality rows.
The `moment_eq_h2=true` option similarly counts exact nonzero translated rows
for the Heisenberg `H^2` moment-equality smoke constraint.
`contiguous_rdm_decomposition=:qmbcertify` is supported here as a shape-only
target for QMBCertify's source-pinned RDM block sizes.  `qmbcertify_base_extra`
adds source-pinned base-block reference metrics for comparison with the active
contiguous base.  `qmbcertify_base_construct=true` switches the active base
target to the QMBCertify source-family blocks, matching the direct linear
constructor's base PSD block shapes without building the model.  Actual
QMBCertify RDM construction is limited to the direct-linear finite-axis
quotient path.  Source-base linear state-opt rows are counted exactly with the
QMBCertify support quotient and sign-canonical row deduplication.
"""
function pauli_translation_structural_targets(
    n_sites::Integer,
    order::Integer;
    sign_symmetry::Bool=true,
    real_moment_matrix::Bool=true,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    reflection_symmetry::Bool=false,
    hamiltonian_width::Integer=2,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition::Symbol=:full,
    contiguous_rdm_support::Symbol=:closed,
    u1_symmetry::Bool=false,
    su2_symmetry::Bool=false,
    base_su2_extend_rdm::Bool=false,
    axis_rotation_symmetry::Bool=false,
    axis_rotation_equalities::Bool=false,
    axis_rotation_quotient::Bool=false,
    moment_eq_h2::Bool=false,
    linear_state_opt_width=nothing,
    psd_state_opt_width=nothing,
    qmbcertify_base_construct::Bool=false,
    qmbcertify_base_extra=nothing,
    qmbcertify_base_three_type::Tuple{<:Integer,<:Integer}=(1, 1),
    scalar_bytes::Integer=sizeof(Float64),
)
    n = Int(n_sites)
    d = Int(order)
    h = Int(hamiltonian_width)
    n > 0 || throw(ArgumentError("`n_sites` must be positive; got $n_sites."))
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    h > 0 || throw(ArgumentError("`hamiltonian_width` must be positive; got $hamiltonian_width."))
    moment_eq_h2 && h != 2 && throw(ArgumentError(
        "`moment_eq_h2=true` structural targets are defined for the nearest-neighbor Heisenberg Hamiltonian with `hamiltonian_width=2`."
    ))
    bytes = Int(scalar_bytes)
    bytes > 0 || throw(ArgumentError("`scalar_bytes` must be positive; got $scalar_bytes."))
    2d < n || throw(ArgumentError(
        "Analytic Pauli translation structural targets require `n_sites > 2order`; got n_sites=$n, order=$d."
    ))

    nontrivial_orbit_size = sum(3^width for width in 1:d; init=0)
    contiguous_orbit_basis_size = 1 + nontrivial_orbit_size
    contiguous_basis_size = 1 + n * nontrivial_orbit_size
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    real_moment_matrix || 0 in sectors || throw(ArgumentError(
        "Momentum sector 0 is required because it carries the normalized identity moment."
    ))
    rdm_ks = _normalize_contiguous_rdm_k(contiguous_rdm_k)
    rdm_decomposition = _normalize_contiguous_rdm_decomposition(
        contiguous_rdm_decomposition,
        u1_symmetry,
        su2_symmetry,
    )
    rdm_support = _normalize_contiguous_rdm_support(contiguous_rdm_support)
    qmbcertify_base_construct && base_su2_extend_rdm && throw(ArgumentError(
        "`qmbcertify_base_construct=true` cannot be combined with `base_su2_extend_rdm=true`."
    ))
    qmbcertify_base_construct && reflection_symmetry && throw(ArgumentError(
        "`qmbcertify_base_construct=true` structural targets currently expose the source-family translation blocks directly; do not also set `reflection_symmetry=true`."
    ))
    qmbcertify_base_construct && !sign_symmetry && throw(ArgumentError(
        "`qmbcertify_base_construct=true` structural targets currently require `sign_symmetry=true`, matching the source-family block split."
    ))
    if base_su2_extend_rdm
        su2_symmetry && rdm_decomposition == :su2 && rdm_support == :extend ||
            throw(ArgumentError(
                "`base_su2_extend_rdm=true` structural targets require " *
                "`su2_symmetry=true`, `contiguous_rdm_decomposition=:su2`, " *
                "and `contiguous_rdm_support=:extend`."
            ))
        real_moment_matrix || throw(ArgumentError(
            "`base_su2_extend_rdm=true` structural targets require `real_moment_matrix=true`."
        ))
    end
    linear_state_width = _normalize_linear_state_opt_width(linear_state_opt_width)
    psd_state_width = _normalize_psd_state_opt_width(psd_state_opt_width)
    extend_addon_support = rdm_support == :extend && !isempty(rdm_ks)
    if qmbcertify_base_construct && !isempty(rdm_ks)
        rdm_decomposition == :qmbcertify || throw(ArgumentError(
            "`qmbcertify_base_construct=true` structural targets support RDM add-ons only with `contiguous_rdm_decomposition=:qmbcertify`."
        ))
        rdm_support == :extend || throw(ArgumentError(
            "`qmbcertify_base_construct=true` structural targets support RDM add-ons only with `contiguous_rdm_support=:extend`."
        ))
        real_moment_matrix || throw(ArgumentError(
            "`qmbcertify_base_construct=true` structural targets support RDM add-ons only with `real_moment_matrix=true`."
        ))
    end
    qmbcertify_base_construct && psd_state_width > 0 && !real_moment_matrix && throw(ArgumentError(
        "`qmbcertify_base_construct=true` structural targets support PSD state-opt add-ons only with `real_moment_matrix=true`."
    ))

    signature_counts = _pauli_axis_signature_counts(d)
    standalone_axis_target = axis_rotation_symmetry && !su2_symmetry
    axis_rotation_quotient && !standalone_axis_target && throw(ArgumentError(
        "`axis_rotation_quotient=true` structural targets require the finite-axis " *
        "target path (`axis_rotation_symmetry=true`, without `su2_symmetry=true`)."
    ))
    qmbcertify_base_construct && standalone_axis_target && throw(ArgumentError(
        "`qmbcertify_base_construct=true` structural targets cannot be combined with finite-axis isotypic targets."
    ))
    qmbcertify_base_construct && axis_rotation_equalities && throw(ArgumentError(
        "`qmbcertify_base_construct=true` structural targets cannot be combined with `axis_rotation_equalities=true`."
    ))
    effective_base_sign_symmetry = sign_symmetry && !standalone_axis_target
    effective_qmbcertify_base_extra =
        qmbcertify_base_construct && isnothing(qmbcertify_base_extra) ?
        0 :
        qmbcertify_base_extra
    qmbcertify_base_reference = _pauli_qmbcertify_base_reference_metrics(
        n,
        d,
        effective_qmbcertify_base_extra;
        three_type=(
            Int(qmbcertify_base_three_type[1]),
            Int(qmbcertify_base_three_type[2]),
        ),
        scalar_bytes=bytes,
        active=qmbcertify_base_construct,
    )
    basis_size = qmbcertify_base_construct ?
        qmbcertify_base_reference.source_row_count + 1 :
        contiguous_basis_size
    orbit_basis_size = qmbcertify_base_construct ?
        qmbcertify_base_reference.orbit_family_count + 1 :
        contiguous_orbit_basis_size
    base_su2_target = base_su2_extend_rdm ?
        pauli_su2_translation_orbit_structural_targets(
            n,
            d;
            real_moment_matrix,
            momenta,
            reflection_symmetry,
            scalar_bytes=bytes,
        ) :
        nothing
    base_logical_sizes, base_psd_sizes = if base_su2_extend_rdm
        copy(something(base_su2_target).logical_block_sizes),
        copy(something(base_su2_target).psd_block_sizes)
    elseif qmbcertify_base_construct
        copy(qmbcertify_base_reference.block_sizes),
        copy(qmbcertify_base_reference.block_sizes)
    elseif standalone_axis_target
        _pauli_translation_axis_target_block_sizes(
            n,
            d,
            sectors;
            reflection_symmetry,
            real_moment_matrix,
        )
    elseif reflection_symmetry
        _pauli_translation_reflection_target_block_sizes(
            orbit_basis_size,
            signature_counts,
            sectors,
            n;
            order=d,
            sign_symmetry=effective_base_sign_symmetry,
            real_moment_matrix,
        )
    else
        _pauli_translation_target_block_sizes(
            orbit_basis_size,
            signature_counts,
            sectors;
            sign_symmetry=effective_base_sign_symmetry,
            real_moment_matrix,
        )
    end
    logical_sizes = copy(base_logical_sizes)
    psd_sizes = copy(base_psd_sizes)

    axis_summary = _pauli_axis_word_orbit_summary(d)
    max_contiguous_rdm_k = 2d
    max_linear_state_opt_width = max(0, 2d - h + 1)
    max_psd_state_opt_width = max(0, fld(2d - h, 2))

    for rdm_k in rdm_ks
        1 <= rdm_k <= n || throw(DomainError(
            rdm_k,
            "`contiguous_rdm_k` must satisfy 1 <= k <= n_sites for structural targets.",
        ))
    end
    if rdm_support == :closed
        for rdm_k in rdm_ks
            rdm_k <= max_contiguous_rdm_k || throw(ArgumentError(
                "`contiguous_rdm_k=$rdm_k` with `contiguous_rdm_support=:closed` requires k <= 2order for structural targets; got order=$d."
            ))
        end
    end
    linear_state_width <= max_linear_state_opt_width || throw(ArgumentError(
        "`linear_state_opt_width=$linear_state_width` is not closed by order=$d and hamiltonian_width=$h."
    ))
    psd_state_width <= max_psd_state_opt_width || throw(ArgumentError(
        "`psd_state_opt_width=$psd_state_width` is not closed by order=$d and hamiltonian_width=$h."
    ))
    psd_state_support = qmbcertify_base_construct ?
        (
            checked=false,
            base_support_closed=true,
            missing_moment_count=0,
            missing_moment_sample=String[],
        ) :
        _pauli_translation_psd_state_opt_target_support_summary(
            select_uint_type(3, n),
            n,
            d,
            h,
            psd_state_width;
            sign_symmetry,
            real_moment_matrix,
            momenta,
        )
    if psd_state_width > 0 &&
            psd_state_support.checked &&
            !psd_state_support.base_support_closed &&
            !extend_addon_support
        shown = join(psd_state_support.missing_moment_sample, ", ")
        throw(ArgumentError(
            "`psd_state_opt_width=$psd_state_width` is degree-closed but not covered " *
            "by the closed base contiguous moment support for structural targets. " *
            "Use `contiguous_rdm_support=:extend` with an extending RDM add-on " *
            "when the combined profile intentionally introduces PSD-state-opt " *
            "moments outside the base support. Missing $(psd_state_support.missing_moment_count) " *
            "moment key(s); first missing: [$shown]."
        ))
    end

    qmbcertify_psd_state_target = qmbcertify_base_construct ?
        _pauli_qmbcertify_psd_state_opt_target_summary(
            n,
            d,
            effective_qmbcertify_base_extra,
            psd_state_width;
            three_type=(
                Int(qmbcertify_base_three_type[1]),
                Int(qmbcertify_base_three_type[2]),
            ),
        ) :
        (
            logical_sizes=Int[],
            psd_sizes=Int[],
            candidate_count=0,
        )
    qmbcertify_linear_state_target = qmbcertify_base_construct ?
        _pauli_qmbcertify_linear_state_opt_target_summary(
            n,
            d,
            effective_qmbcertify_base_extra,
            linear_state_width;
            three_type=(
                Int(qmbcertify_base_three_type[1]),
                Int(qmbcertify_base_three_type[2]),
            ),
            rdm_ks=rdm_ks,
            psd_state_opt_width=psd_state_width,
        ) :
        (row_count=0, candidate_count=0)
    psd_state_logical_sizes, psd_state_psd_sizes = qmbcertify_base_construct ?
        (
            qmbcertify_psd_state_target.logical_sizes,
            qmbcertify_psd_state_target.psd_sizes,
        ) :
        _pauli_psd_state_opt_target_block_sizes(
            psd_state_width,
            sectors;
            sign_symmetry,
            real_moment_matrix,
            full_realification=axis_rotation_quotient,
        )
    append!(logical_sizes, psd_state_logical_sizes)
    append!(psd_sizes, psd_state_psd_sizes)
    rdm_logical_sizes, rdm_psd_sizes = _pauli_translation_rdm_target_block_sizes(
        rdm_ks,
        rdm_decomposition,
        real_moment_matrix,
    )
    append!(logical_sizes, rdm_logical_sizes)
    append!(psd_sizes, rdm_psd_sizes)
    rdm_coefficient_domains, rdm_exact_coefficient_domains =
        _pauli_translation_rdm_target_block_domains(
            rdm_ks,
            rdm_decomposition,
        )
    block_coefficient_domains = vcat(
        base_su2_extend_rdm ?
            Union{Nothing,Symbol}[something(base_su2_target).block_coefficient_domains...] :
            qmbcertify_base_construct ?
            _pauli_target_domain_vector(
                :cyclotomic_float64,
                length(base_psd_sizes),
            ) :
            _pauli_target_domain_vector(
                standalone_axis_target ? :real_algebraic_float64 : :cyclotomic_float64,
                length(base_psd_sizes),
            ),
        _pauli_target_domain_vector(
            :cyclotomic_float64,
            length(psd_state_psd_sizes),
        ),
        rdm_coefficient_domains,
    )
    block_exact_coefficient_domains = vcat(
        base_su2_extend_rdm ?
            Union{Nothing,Symbol}[something(base_su2_target).block_exact_coefficient_domains...] :
            qmbcertify_base_construct ?
            _pauli_target_domain_vector(:cyclotomic, length(base_psd_sizes)) :
            _pauli_target_domain_vector(
                standalone_axis_target ? :axis_character_projector :
                (reflection_symmetry ? :cyclotomic_sqrt_rational : :cyclotomic),
                length(base_psd_sizes),
            ),
        _pauli_target_domain_vector(:cyclotomic, length(psd_state_psd_sizes)),
        rdm_exact_coefficient_domains,
    )
    rdm_zero_row_count = _pauli_translation_rdm_target_zero_row_count(
        rdm_ks,
        rdm_decomposition,
        real_moment_matrix,
    )

    psd_dense_entries = sum(size * size for size in psd_sizes; init=0)
    psd_symmetric_entries = sum(
        real_moment_matrix ? size * (size + 1) ÷ 2 : size * size for size in psd_sizes;
        init=0,
    )
    product_cache_stats = qmbcertify_base_construct ?
        (
            product_cache_hits=0,
            product_cache_misses=0,
            product_cache_lookups=0,
            product_cache_entries=0,
            product_cache_hit_rate=0.0,
        ) :
        _pauli_translation_target_product_cache_stats(
            orbit_basis_size,
            signature_counts,
            sectors;
            sign_symmetry=effective_base_sign_symmetry,
        )
    effective_axis_rotation_equalities =
        axis_rotation_equalities || standalone_axis_target || axis_rotation_quotient
    axis_rotation_quotient_stats = effective_axis_rotation_equalities ?
        _pauli_translation_target_axis_rotation_quotient_stats(
            select_uint_type(3, n),
            n,
            d;
            sign_symmetry=standalone_axis_target ? false : sign_symmetry,
            real_moment_matrix,
            momenta,
        ) :
        (
            equality_row_count=0,
            raw_moment_key_count=0,
            moment_class_count=0,
            quotient_moment_key_count=0,
            forced_zero_moment_class_count=0,
            reduction_ratio=0.0,
        )
    axis_rotation_equality_row_count = axis_rotation_quotient_stats.equality_row_count
    linear_state_opt_row_count = h == 2 && qmbcertify_base_construct ?
        qmbcertify_linear_state_target.row_count :
        h == 2 ?
        _pauli_translation_linear_state_opt_target_row_count(
            n,
            linear_state_width;
            sign_symmetry,
            real_moment_matrix,
        ) :
        0
    moment_equality_row_count = moment_eq_h2 ?
        _pauli_translation_moment_eq_h2_target_row_count(
            n,
            d;
            sign_symmetry,
            real_moment_matrix,
            momenta,
        ) :
        0
    add_on_zero_row_count =
        axis_rotation_equality_row_count +
        rdm_zero_row_count +
        linear_state_opt_row_count +
        moment_equality_row_count
    known_zero_histogram = _pauli_target_zero_feature_histogram(
        axis_rotation_equality_row_count,
        rdm_zero_row_count,
        rdm_decomposition,
        linear_state_opt_row_count,
        moment_equality_row_count,
    )
    base_decomposition = base_su2_extend_rdm ?
        (reflection_symmetry ? :translation_su2_reflection : :translation_su2) :
        qmbcertify_base_construct ?
        :qmbcertify_base :
        standalone_axis_target ?
        (reflection_symmetry ? :translation_reflection_axis_irrep : :translation_axis_irrep) :
        (reflection_symmetry ? :translation_reflection : :translation)
    su2_base_accounting = base_su2_extend_rdm ?
        (
            su2_base_full_dense_entries=something(base_su2_target).su2_full_dense_entries,
            su2_base_active_dense_entries=something(base_su2_target).su2_active_dense_entries,
            su2_base_reduced_dense_entries=something(base_su2_target).su2_reduced_dense_entries,
            su2_base_offblock_entry_count=something(base_su2_target).offblock_entry_count,
            su2_base_copy_entry_count=something(base_su2_target).copy_entry_count,
            su2_base_accounted_entry_count=something(base_su2_target).accounted_entry_count,
            su2_base_accounting_records=something(base_su2_target).su2_accounting_records,
            su2_base_accounting_record_count=something(base_su2_target).su2_accounting_record_count,
            su2_base_singlet_channel_count=something(base_su2_target).singlet_channel_count,
            su2_base_singlet_channel_support_counts=something(base_su2_target).singlet_channel_support_counts,
            su2_base_singlet_channel_equality_row_count=something(base_su2_target).singlet_channel_equality_row_count,
            su2_base_offblock_zero_row_budget=something(base_su2_target).wigner_offblock_zero_row_budget,
            su2_base_magnetic_copy_zero_row_budget=something(base_su2_target).wigner_magnetic_copy_zero_row_budget,
            su2_base_zero_row_budget=something(base_su2_target).wigner_zero_row_budget,
        ) :
        (
            su2_base_full_dense_entries=0,
            su2_base_active_dense_entries=0,
            su2_base_reduced_dense_entries=0,
            su2_base_offblock_entry_count=0,
            su2_base_copy_entry_count=0,
            su2_base_accounted_entry_count=0,
            su2_base_accounting_records=Any[],
            su2_base_accounting_record_count=0,
            su2_base_singlet_channel_count=0,
            su2_base_singlet_channel_support_counts=Pair{Int,Int}[],
            su2_base_singlet_channel_equality_row_count=0,
            su2_base_offblock_zero_row_budget=0,
            su2_base_magnetic_copy_zero_row_budget=0,
            su2_base_zero_row_budget=0,
        )
    qmbcertify_rdm_references = _pauli_qmbcertify_rdm_reference_metrics(
        rdm_ks;
        scalar_bytes=bytes,
    )
    qmbcertify_rdm_reference_summary = _pauli_qmbcertify_rdm_reference_summary(
        qmbcertify_rdm_references,
    )
    return (
        n_sites=n,
        order=d,
        hamiltonian_width=h,
        basis_size=basis_size,
        orbit_basis_size=orbit_basis_size,
        translation_orbit_count=orbit_basis_size,
        momentum_sectors=sectors,
        momentum_sector_count=length(sectors),
        sign_symmetry=sign_symmetry,
        reflection_symmetry=reflection_symmetry,
        axis_rotation_symmetry=standalone_axis_target,
        axis_rotation_quotient=axis_rotation_quotient,
        real_moment_matrix=real_moment_matrix,
        sign_signature_sizes=sort!(collect(signature_counts); by=first),
        logical_block_sizes=logical_sizes,
        psd_block_sizes=psd_sizes,
        block_coefficient_domains=block_coefficient_domains,
        block_exact_coefficient_domains=block_exact_coefficient_domains,
        block_coefficient_domain_histogram=_value_histogram_pairs(block_coefficient_domains),
        block_exact_coefficient_domain_histogram=_value_histogram_pairs(block_exact_coefficient_domains),
        n_blocks=length(psd_sizes),
        logical_max_block=isempty(logical_sizes) ? 0 : maximum(logical_sizes),
        psd_max_block=isempty(psd_sizes) ? 0 : maximum(psd_sizes),
        logical_total_block_side=sum(logical_sizes; init=0),
        psd_total_block_side=sum(psd_sizes; init=0),
        logical_block_histogram=_histogram_pairs(logical_sizes),
        psd_block_histogram=_histogram_pairs(psd_sizes),
        logical_block_feature_histogram=_pauli_target_block_feature_histogram(
            base_logical_sizes,
            psd_state_logical_sizes,
            rdm_logical_sizes;
            base_decomposition,
            rdm_decomposition,
        ),
        psd_block_feature_histogram=_pauli_target_block_feature_histogram(
            base_psd_sizes,
            psd_state_psd_sizes,
            rdm_psd_sizes;
            base_decomposition,
            rdm_decomposition,
        ),
        psd_dense_entries=psd_dense_entries,
        psd_symmetric_entries=psd_symmetric_entries,
        psd_dense_bytes=psd_dense_entries * bytes,
        psd_symmetric_bytes=psd_symmetric_entries * bytes,
        product_cache_stats...,
        max_contiguous_rdm_k=max_contiguous_rdm_k,
        max_linear_state_opt_width=max_linear_state_opt_width,
        max_psd_state_opt_width=max_psd_state_opt_width,
        axis_orbit_closed=!qmbcertify_base_construct,
        axis_orbit_basis_size=qmbcertify_base_construct ? 0 : axis_summary.axis_orbit_count,
        axis_orbit_size_histogram=qmbcertify_base_construct ?
            Pair{Int,Int}[] :
            axis_summary.axis_orbit_size_histogram,
        max_axis_orbit_size=qmbcertify_base_construct ? 0 : axis_summary.max_axis_orbit_size,
        axis_reduction_ratio=qmbcertify_base_construct ? 1.0 : axis_summary.axis_orbit_count == 0 ?
            Inf :
            orbit_basis_size / axis_summary.axis_orbit_count,
        base_logical_block_sizes=base_logical_sizes,
        base_psd_block_sizes=base_psd_sizes,
        rdm_logical_block_sizes=rdm_logical_sizes,
        rdm_psd_block_sizes=rdm_psd_sizes,
        contiguous_rdm_zero_row_count=rdm_zero_row_count,
        psd_state_opt_logical_block_sizes=psd_state_logical_sizes,
        psd_state_opt_psd_block_sizes=psd_state_psd_sizes,
        contiguous_rdm_k=copy(rdm_ks),
        contiguous_rdm_decomposition=rdm_decomposition,
        contiguous_rdm_support=rdm_support,
        base_su2_extend_rdm=base_su2_extend_rdm,
        su2_base_accounting...,
        psd_state_opt_support_policy=extend_addon_support ? :extend : :closed,
        psd_state_opt_base_support_checked=psd_state_support.checked,
        psd_state_opt_base_support_closed=psd_state_support.base_support_closed,
        psd_state_opt_missing_moment_count=psd_state_support.missing_moment_count,
        psd_state_opt_missing_moment_sample=psd_state_support.missing_moment_sample,
        axis_rotation_equalities=effective_axis_rotation_equalities,
        axis_rotation_equality_row_count=axis_rotation_equality_row_count,
        axis_rotation_raw_moment_key_count=axis_rotation_quotient_stats.raw_moment_key_count,
        axis_rotation_moment_class_count=axis_rotation_quotient_stats.moment_class_count,
        axis_rotation_quotient_moment_key_count=axis_rotation_quotient_stats.quotient_moment_key_count,
        axis_rotation_forced_zero_moment_class_count=axis_rotation_quotient_stats.forced_zero_moment_class_count,
        axis_rotation_moment_quotient_reduction_ratio=axis_rotation_quotient_stats.reduction_ratio,
        moment_eq_h2=moment_eq_h2,
        moment_equality_row_count=moment_equality_row_count,
        add_on_zero_row_count=add_on_zero_row_count,
        known_zero_constraint_feature_histogram=known_zero_histogram,
        qmbcertify_rdm_reference_blocks=qmbcertify_rdm_references,
        qmbcertify_rdm_reference_block_sizes=qmbcertify_rdm_reference_summary.block_sizes,
        qmbcertify_rdm_reference_n_blocks=qmbcertify_rdm_reference_summary.n_blocks,
        qmbcertify_rdm_reference_max_block=qmbcertify_rdm_reference_summary.max_block,
        qmbcertify_rdm_reference_total_block_side=qmbcertify_rdm_reference_summary.total_block_side,
        qmbcertify_rdm_reference_dense_entries=qmbcertify_rdm_reference_summary.dense_entries,
        qmbcertify_rdm_reference_symmetric_entries=qmbcertify_rdm_reference_summary.symmetric_entries,
        qmbcertify_rdm_reference_dense_bytes=qmbcertify_rdm_reference_summary.dense_bytes,
        qmbcertify_rdm_reference_symmetric_bytes=qmbcertify_rdm_reference_summary.symmetric_bytes,
        qmbcertify_rdm_reference_requires_construction=qmbcertify_rdm_reference_summary.requires_construction,
        qmbcertify_base_reference_extra=qmbcertify_base_reference.extra,
        qmbcertify_base_reference_family_count_by_parity=
            qmbcertify_base_reference.family_count_by_parity,
        qmbcertify_base_reference_blocks=qmbcertify_base_reference.block_records,
        qmbcertify_base_reference_block_sizes=qmbcertify_base_reference.block_sizes,
        qmbcertify_base_reference_block_histogram=qmbcertify_base_reference.block_histogram,
        qmbcertify_base_reference_n_blocks=qmbcertify_base_reference.n_blocks,
        qmbcertify_base_reference_max_block=qmbcertify_base_reference.max_block,
        qmbcertify_base_reference_total_block_side=
            qmbcertify_base_reference.total_block_side,
        qmbcertify_base_reference_dense_entries=qmbcertify_base_reference.dense_entries,
        qmbcertify_base_reference_symmetric_entries=
            qmbcertify_base_reference.symmetric_entries,
        qmbcertify_base_reference_dense_bytes=qmbcertify_base_reference.dense_bytes,
        qmbcertify_base_reference_symmetric_bytes=
            qmbcertify_base_reference.symmetric_bytes,
        qmbcertify_base_reference_support_nonzero_row_count=
            qmbcertify_base_reference.support_nonzero_row_count,
        qmbcertify_base_reference_support_zero_row_count=
            qmbcertify_base_reference.support_zero_row_count,
        qmbcertify_base_reference_support_diagonal_nonzero_row_count=
            qmbcertify_base_reference.support_diagonal_nonzero_row_count,
        qmbcertify_base_reference_support_offdiagonal_nonzero_row_count=
            qmbcertify_base_reference.support_offdiagonal_nonzero_row_count,
        qmbcertify_base_reference_support_unique_count=
            qmbcertify_base_reference.support_unique_count,
        qmbcertify_base_reference_support_diagonal_unique_count=
            qmbcertify_base_reference.support_diagonal_unique_count,
        qmbcertify_base_reference_support_offdiagonal_unique_count=
            qmbcertify_base_reference.support_offdiagonal_unique_count,
        qmbcertify_base_reference_support_word_length_histogram=
            qmbcertify_base_reference.support_word_length_histogram,
        qmbcertify_base_reference_active=qmbcertify_base_reference.active,
        qmbcertify_base_reference_requires_construction=
            qmbcertify_base_reference.requires_construction,
        linear_state_opt_width=linear_state_width,
        linear_state_opt_row_count=linear_state_opt_row_count,
        linear_state_opt_candidate_count=qmbcertify_base_construct ?
            qmbcertify_linear_state_target.candidate_count :
            _pauli_state_opt_test_count(
                linear_state_width;
                sign_symmetry,
            ),
        psd_state_opt_width=psd_state_width,
        psd_state_opt_candidate_count=qmbcertify_base_construct ?
            qmbcertify_psd_state_target.candidate_count :
            _pauli_state_opt_test_count(
                psd_state_width;
                sign_symmetry,
            ),
        solve_supported=false,
        solve_blocker=:structural_target_only,
        solve_blocker_reason=_TRANSLATION_SOLVE_STRUCTURAL_TARGET_REASON,
        estimated_model_size_gate_status=:blocked_missing_scalar_equality_estimate,
        estimated_model_size_gate_reason=_TRANSLATION_STRUCTURAL_MODEL_SIZE_GATE_REASON,
        solve_unsupported_block_features=Symbol[],
        solve_unsupported_zero_features=Symbol[],
        requires_construction=false,
        assumptions=(
            periodic_contiguous_basis=!qmbcertify_base_construct,
            qmbcertify_source_basis=qmbcertify_base_construct,
            no_short_orbits=!qmbcertify_base_construct,
            n_sites_gt_2order=true,
            reflection_symmetry=reflection_symmetry,
            hamiltonian_width=h,
            add_on_zero_row_counts=true,
            complete_add_on_zero_row_count=(iszero(linear_state_width) || h == 2) &&
                rdm_decomposition in (:full, :u1),
        ),
    )
end

"""
    TranslationInvariantResult

Result returned by [`pauli_translation_invariant_nctssos`](@ref).
"""
struct TranslationInvariantResult{T,MP}
    objective::T
    model::GenericModel{T}
    moment_problem::MP
    report::TranslationInvariantReport
    solve_time_ns::Int
end

translation_zero_origin_histogram(result::TranslationInvariantResult) =
    translation_zero_origin_histogram(result.moment_problem)

translation_linear_provenance(result::TranslationInvariantResult) =
    translation_linear_provenance(result.moment_problem)

function Base.show(io::IO, result::TranslationInvariantResult)
    println(io, "Translation-Invariant Optimization Result")
    println(io, "Objective: ", result.objective)
    println(io, "Solve time (ns): ", result.solve_time_ns)
    print(io, "Report: ", result.report)
end

translation_block_histogram(
    result::TranslationInvariantResult;
    kind::Symbol=:logical,
) = translation_block_histogram(result.report; kind)

translation_block_feature_histogram(
    result::TranslationInvariantResult;
    kind::Symbol=:logical,
) = translation_block_feature_histogram(result.report; kind)

translation_symmetry_profile(result::TranslationInvariantResult) =
    translation_symmetry_profile(result.report)

translation_solve_support(result::TranslationInvariantResult) =
    translation_solve_support(result.report)

translation_report_metrics(
    result::TranslationInvariantResult;
    scalar_bytes::Integer=sizeof(Float64),
) = translation_report_metrics(result.report; scalar_bytes)

"""
    pauli_translation_invariant_moment_relaxation(pop, (σx, σy, σz), order; sign_symmetry=true, reflection_symmetry=false, axis_rotation_symmetry=false, axis_rotation_equalities=false, momenta=nothing, contiguous_rdm_k=nothing, contiguous_rdm_decomposition=:full, contiguous_rdm_support=:closed, u1_symmetry=false, su2_symmetry=false, linear_state_opt_width=nothing, psd_state_opt_width=nothing)

Construct a periodic-chain Pauli moment relaxation directly in translation
(momentum) sectors, without building the full site-space moment matrix.

This is an intentionally narrow large-spin-chain path:
- `PauliAlgebra` polynomial objectives with translation-compatible scalar
  equality, inequality, and one-sided moment-equality constraints;
- periodic translation by one site on the declared chain `1:N`;
- a contiguous local half-basis from [`pauli_contiguous_chain_basis`](@ref);
- optional contiguous `k`-site RDM PSD blocks via `contiguous_rdm_k`, either
  full, U(1)-decomposed with `contiguous_rdm_decomposition=:u1`, or
  SU(2)-decomposed with `contiguous_rdm_decomposition=:su2`;
  by default RDM moments must be covered by the base moment PSD blocks, while
  `contiguous_rdm_support=:extend` admits new moment keys from the RDM block
  and from PSD-state-opt blocks in the same extending-support add-on profile;
  `contiguous_rdm_decomposition=:qmbcertify` is supported only on the
  direct-linear finite-axis quotient path, where the shared QMBCertify-style
  RDM PSD blocks use the same translation/reflection/global-axis
  identifications as the base moment quotient;
- optional first-order linear state-opt rows `⟨i[H,A]⟩ = 0` via
  `linear_state_opt_width`;
- optional PSD state-opt blocks
  `⟨vHw* - 1/2(Hvw* + vw*H)⟩` via `psd_state_opt_width`;
- optional global Pauli-axis rotation moment equalities via
  `axis_rotation_equalities=true`; these are scalar equality rows, not PSD
  block reductions;
- optional `(ℤ₂)^2` Heisenberg sign-symmetry splitting, enabled by default.
- optional reflection splitting for momentum sectors fixed by reflection
  (`k=0`, and `k=N/2` when `N` is even) via `reflection_symmetry=true`.
- global Pauli-axis rotation diagnostics and guarded PSD splitting.  The
  symbolic and direct-linear base paths support the validated
  `axis_rotation_symmetry=true` case, including reflection-adapted sectors;
  first-order linear state-opt rows, contiguous RDM blocks, and PSD state-opt
  blocks compose as unsplit add-ons.  `sign_symmetry=true` is treated as
  subsumed by the finite axis group for this split.  With `su2_symmetry=true`,
  finite axis rotations are treated as redundant because they are subsumed by
  SU(2).

Set `u1_symmetry=true` only when the objective and constraints are known to be
globally U(1)-invariant.  The U(1)-decomposed RDM path verifies charge
neutrality and enforces off-sector RDM entries as zero before adding the
Hamming-weight PSD blocks.

Set `su2_symmetry=true` only for globally SU(2)-invariant Pauli spin problems.
The SU(2)-decomposed RDM path checks global Pauli-axis rotation invariance,
enforces the local RDM Schur-basis zero/equality rows, and adds PSD constraints
only on multiplicity blocks.  The base translation moment matrix is SU(2)
reduced for the guarded base branch, closed-support full/U(1)/SU(2) RDM
  add-ons, and `psd_state_opt_width` add-ons.  Extending-support RDM add-ons keep ordinary
translation moment blocks because they introduce moments outside the reduced
base basis.  In the same extending-support add-on profile, linear state-opt
rows may register additional free moment keys just like RDM and PSD state-opt
add-ons.

By default this builds the paper-style real primal moment matrix: complex
momentum blocks are realified once, conjugate momenta are not duplicated, and
moment variables are real.  Set `real_moment_matrix=false` only for debugging the
older complex Hermitian block form.

When `sign_symmetry=true`, the objective must be invariant under the global
Heisenberg sign flips.  When `reflection_symmetry=true`, the objective and scalar
constraints must be invariant under chain reflection after translation reduction.
If `momenta` is supplied with `real_moment_matrix=false`, it must include sector
`0` because the normalized identity moment lives there.

For the XXX chain with `N=100, order=4`, the default basis has 12,001 site-space
monomials; the logical complex momentum blocks have side at most 31 and the
solver-facing real PSD blocks have side at most 62.
"""
function _translation_reject_direct_base_extension(name::Symbol, value)
    isnothing(value) && return nothing
    if value isa Bool
        value || return nothing
    end
    throw(ArgumentError(
        "Direct base translation linear relaxation does not yet support `$name`; " *
        "use `pauli_translation_invariant_moment_relaxation` for the symbolic path."
    ))
end

function _pauli_monomial_from_moment_key(
    ::Type{M},
    key::AbstractVector{T},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    return M(copy(key))
end

_translation_registration_key_token(key) = key
_translation_registration_key_token(key::AbstractVector) = key

function _translation_registered_key_tokens(builder::MomentLinearBuilder)
    tokens = Set{Any}()
    for key in keys(builder.key_to_monomial)
        push!(tokens, _translation_registration_key_token(key))
    end
    return tokens
end

function _deduplicate_translation_zero_constraints!(
    builder::MomentLinearBuilder{K,C,M},
) where {K,C,M}
    seen = Set{_LinearMomentFormCacheKey{K,C}}()
    sizehint!(seen, length(builder.zero_constraints))
    deduped = ScalarLinearConstraint{K,C}[]
    sizehint!(deduped, length(builder.zero_constraints))
    for zc in builder.zero_constraints
        if zc.kind != :zero
            push!(deduped, zc)
            continue
        end
        key = _linear_moment_form_cache_key(zc.form)
        # push! already performs the membership probe; avoid a separate lookup.
        seen_count = length(seen)
        push!(seen, key)
        length(seen) == seen_count && continue
        push!(deduped, zc)
    end
    builder.zero_constraints = deduped
    return length(seen)
end

function _translation_builder_from_linear(
    linear::MomentLinearData{K,C,M},
) where {K,C,M<:NormalMonomial}
    builder = MomentLinearBuilder(K, C, M)
    builder.identity = _owned_moment_key(K, linear.identity)

    key_to_monomial = Dict{K,M}()
    sizehint!(key_to_monomial, length(linear.key_to_monomial))
    for (key, mono) in linear.key_to_monomial
        _register_builder_moment_key!(key_to_monomial, K, key, mono)
    end
    builder.key_to_monomial = key_to_monomial

    builder.objective_terms = _owned_builder_pairs(K, C, linear.objective_lin.terms)

    blocks = PSDBlockLin{K,C,M}[]
    sizehint!(blocks, length(linear.psd_blocks_lin))
    for block in linear.psd_blocks_lin
        entries = Matrix{LinearMomentForm{K,C}}(undef, size(block.entries))
        for idx in eachindex(block.entries)
            entries[idx] = _owned_builder_form(K, C, block.entries[idx])
        end
        push!(blocks, PSDBlockLin{K,C,M}(block.size, entries, block.meta))
    end
    builder.psd_blocks_lin = blocks
    builder.psd_block_constraint_idx = copy(linear.psd_block_constraint_idx)

    zero_constraints = ScalarLinearConstraint{K,C}[]
    sizehint!(zero_constraints, length(linear.zero_constraints))
    for zc in linear.zero_constraints
        push!(
            zero_constraints,
            ScalarLinearConstraint(
                _owned_builder_form(K, C, zc.form),
                zc.kind,
                zc.origin;
                trusted_self_adjoint=zc.trusted_self_adjoint,
            ),
        )
    end
    builder.zero_constraints = zero_constraints
    return builder
end

function _translation_append_builder_from_linear(
    linear::MomentLinearData{K,C,M},
) where {K,C,M<:NormalMonomial}
    builder = MomentLinearBuilder(K, C, M)
    builder.identity = _owned_moment_key(K, linear.identity)

    # Append-only callers never mutate existing blocks/forms; reuse them and
    # copy only the containers that receive appended constraints.
    builder.key_to_monomial = copy(linear.key_to_monomial)
    builder.objective_terms = copy(linear.objective_lin.terms)
    builder.psd_blocks_lin = copy(linear.psd_blocks_lin)
    builder.psd_block_constraint_idx = copy(linear.psd_block_constraint_idx)
    builder.zero_constraints = copy(linear.zero_constraints)
    return builder
end

function _register_translation_linear_form_keys!(
    builder::MomentLinearBuilder{K,C,M},
    form::LinearMomentForm{K,FC},
) where {K,C,FC,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    for (key, _) in form
        register_moment!(builder, key, _pauli_monomial_from_moment_key(M, key))
    end
    return nothing
end

function _register_translation_linear_entries!(
    builder::MomentLinearBuilder{K,C,M},
    entries::AbstractMatrix{LinearMomentForm{K,FC}},
) where {K,C,FC,M}
    for form in entries
        _register_translation_linear_form_keys!(builder, form)
    end
    return nothing
end

function _register_translation_linear_form_keys!(
    builder::MomentLinearBuilder{K,C,M},
    form::LinearMomentForm{K,FC},
    registered_key_tokens::Set{Any},
) where {K,C,FC,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    for (key, _) in form
        token = _translation_registration_key_token(key)
        token in registered_key_tokens && continue
        register_moment!(builder, key, _pauli_monomial_from_moment_key(M, key))
        push!(registered_key_tokens, token)
    end
    return nothing
end

function _register_translation_linear_entries!(
    builder::MomentLinearBuilder{K,C,M},
    entries::AbstractMatrix{LinearMomentForm{K,FC}},
    registered_key_tokens::Set{Any},
) where {K,C,FC,M}
    for form in entries
        _register_translation_linear_form_keys!(builder, form, registered_key_tokens)
    end
    return nothing
end

function _quotient_translation_linear_entries!(
    builder::MomentLinearBuilder{K,C,M},
    entries::AbstractMatrix{<:LinearMomentForm},
    quotient,
) where {K,C,M}
    rewritten = Matrix{LinearMomentForm{K,C}}(undef, size(entries))
    for idx in eachindex(entries)
        rewritten[idx] = _axis_quotient_rewrite_form!(
            builder.key_to_monomial,
            K,
            C,
            entries[idx],
            quotient,
        )
    end
    return rewritten
end

function _quotient_translation_linear_entries!(
    builder::MomentLinearBuilder{K,C,M},
    entries::AbstractMatrix{LinearMomentForm{K,FC}},
    quotient,
    cache::Dict{_LinearMomentFormCacheKey{K,FC},LinearMomentForm{K,C}},
) where {K,C,FC,M}
    rewritten = Matrix{LinearMomentForm{K,C}}(undef, size(entries))
    for idx in eachindex(entries)
        form = entries[idx]
        cache_key = _linear_moment_form_cache_key(form)
        rewritten[idx] = get!(cache, cache_key) do
            _axis_quotient_rewrite_form!(
                builder.key_to_monomial,
                K,
                C,
                form,
                quotient,
            )
        end
    end
    return rewritten
end

function _prepare_translation_linear_entries!(
    builder::MomentLinearBuilder{K,C,M},
    entries::AbstractMatrix{<:LinearMomentForm},
    registered_key_tokens::Set{Any},
    quotient,
    quotient_cache=nothing,
) where {K,C,M}
    if quotient === nothing
        _register_translation_linear_entries!(builder, entries, registered_key_tokens)
        return entries
    end
    quotient_cache === nothing ||
        return _quotient_translation_linear_entries!(builder, entries, quotient, quotient_cache)
    return _quotient_translation_linear_entries!(builder, entries, quotient)
end

function _quotient_translation_linear_entry_forms(
    ::Type{M},
    entries::AbstractMatrix{LinearMomentForm{K,C}},
    quotient,
) where {K,C,M}
    quotient === nothing && return entries
    key_to_monomial = Dict{K,M}()
    rewritten = Matrix{LinearMomentForm{K,C}}(undef, size(entries))
    for idx in eachindex(entries)
        rewritten[idx] = _axis_quotient_rewrite_form!(
            key_to_monomial,
            K,
            C,
            entries[idx],
            quotient,
        )
    end
    return rewritten
end

function _translation_stored_adjoint_key!(
    builder::MomentLinearBuilder{K,C,M},
    stored_key::K,
    registered_key_tokens=nothing,
    adjoint_key_cache=nothing,
) where {K,C,A<:AlgebraType,T<:Unsigned,M<:NormalMonomial{A,T}}
    if adjoint_key_cache !== nothing
        cached = get(adjoint_key_cache, stored_key, nothing)
        cached === nothing || return cached
    end

    mono = builder.key_to_monomial[stored_key]
    adjoint_mono = _moment_linear_adjoint_monomial(mono)
    adjoint_key = register_moment!(builder, adjoint_mono)
    adjoint_stored_key = _builder_stored_key(builder.key_to_monomial, adjoint_key)
    adjoint_stored_key === nothing && (adjoint_stored_key = adjoint_key)
    if registered_key_tokens !== nothing
        push!(registered_key_tokens, _translation_registration_key_token(adjoint_stored_key))
    end
    if adjoint_key_cache !== nothing
        adjoint_key_cache[stored_key] = adjoint_stored_key
        reverse_stored_key = _builder_stored_key(builder.key_to_monomial, adjoint_stored_key)
        reverse_stored_key === nothing ||
            (adjoint_key_cache[reverse_stored_key] = stored_key)
    end
    return adjoint_stored_key
end

function _translation_complex_zero_linear_parts!(
    builder::MomentLinearBuilder{K,C,M},
    form::LinearMomentForm{K,FC};
    registered_key_tokens=nothing,
    register_keys::Bool=true,
    adjoint_key_cache=nothing,
) where {K,C,FC,A<:AlgebraType,T<:Unsigned,M<:NormalMonomial{A,T}}
    if register_keys
        if registered_key_tokens === nothing
            _register_translation_linear_form_keys!(builder, form)
        else
            _register_translation_linear_form_keys!(builder, form, registered_key_tokens)
        end
    end

    re_pairs = Pair{K,C}[]
    im_pairs = Pair{K,C}[]
    sizehint!(re_pairs, 2 * length(form))
    sizehint!(im_pairs, 2 * length(form))
    half = convert(C, 0.5)
    neg_half_im = convert(C, -0.5im)
    pos_half_im = convert(C, 0.5im)
    for (key, coef) in form
        stored_key = _builder_stored_key(builder.key_to_monomial, key)
        stored_key === nothing && throw(ArgumentError(
            "Internal error: direct translation zero form was not registered at key $(repr(key))."
        ))
        adjoint_stored_key = _translation_stored_adjoint_key!(
            builder,
            stored_key,
            registered_key_tokens,
            adjoint_key_cache,
        )
        converted_coef = convert(C, coef)
        adjoint_coef = convert(C, conj(converted_coef))
        push!(re_pairs, stored_key => half * converted_coef)
        push!(re_pairs, adjoint_stored_key => half * adjoint_coef)
        push!(im_pairs, stored_key => neg_half_im * converted_coef)
        push!(im_pairs, adjoint_stored_key => pos_half_im * adjoint_coef)
    end
    return (
        _linear_moment_form_from_owned_pairs!(re_pairs),
        _linear_moment_form_from_owned_pairs!(im_pairs),
    )
end

function _translation_builder_psd_moment_basis(
    builder::MomentLinearBuilder{K,C,M},
) where {K,C,M}
    seen = Set{K}()
    basis = M[]
    for block in builder.psd_blocks_lin
        for form in block.entries
            for (key, _) in form
                key in seen && continue
                push!(seen, key)
                push!(basis, builder.key_to_monomial[key])
            end
        end
    end
    return sorted_unique!(basis)
end

function _linear_entries_moment_basis(
    ::Type{M},
    entries::AbstractMatrix{LinearMomentForm{K,C}},
) where {K,C,M<:NormalMonomial}
    seen = Set{Any}()
    basis = M[]
    for form in entries
        for (key, _) in form
            token = _translation_registration_key_token(key)
            token in seen && continue
            push!(seen, token)
            push!(basis, _pauli_monomial_from_moment_key(M, key))
        end
    end
    return basis
end

_translation_zero_linear_origin(label) = TranslationZeroOrigin(label, 1, 1, :scalar)

_pauli_su2_singlet_channel_linear_origin(label) =
    PauliSU2SingletChannelEqualityOrigin(label, 1, 1, :scalar)

function _indexed_real_part_form_coefficients(
    form::LinearMomentForm{Int,FC},
    index_to_key::AbstractVector{K},
    ::Type{R};
    atol,
) where {K,FC<:Number,R<:Real}
    pairs = Pair{K,R}[]
    sizehint!(pairs, length(form))
    tolerance = R(atol)
    for (idx, coef) in form
        value = _clean_real_part(real(coef), R, tolerance)
        iszero(value) || push!(pairs, index_to_key[idx] => value)
    end
    return LinearMomentForm{K,R}(pairs, Val(:trusted))
end

function _indexed_imag_part_form_coefficients(
    form::LinearMomentForm{Int,FC},
    index_to_key::AbstractVector{K},
    ::Type{R};
    atol,
) where {K,FC<:Number,R<:Real}
    pairs = Pair{K,R}[]
    sizehint!(pairs, length(form))
    tolerance = R(atol)
    for (idx, coef) in form
        value = _clean_real_part(imag(coef), R, tolerance)
        iszero(value) || push!(pairs, index_to_key[idx] => value)
    end
    return LinearMomentForm{K,R}(pairs, Val(:trusted))
end

function _add_translation_zero_indexed_linear_form!(
    builder::MomentLinearBuilder{K,C,M},
    form::LinearMomentForm{Int,FC},
    index_to_key::AbstractVector{K};
    phase_atol,
    label,
    registered_key_tokens=nothing,
    origin_factory=_translation_zero_linear_origin,
    register_keys::Bool=true,
    adjoint_key_cache=nothing,
) where {K,C,FC,M}
    if C <: Real
        re = _indexed_real_part_form_coefficients(form, index_to_key, C; atol=phase_atol)
        if !isempty(re)
            if register_keys
                if registered_key_tokens === nothing
                    _register_translation_linear_form_keys!(builder, re)
                else
                    _register_translation_linear_form_keys!(builder, re, registered_key_tokens)
                end
            end
            _add_zero_constraint_trusted!(
                builder,
                re,
                origin_factory(_translation_zero_component_label(label, :real)),
            )
        end
        im = _indexed_imag_part_form_coefficients(form, index_to_key, C; atol=phase_atol)
        if !isempty(im)
            if register_keys
                if registered_key_tokens === nothing
                    _register_translation_linear_form_keys!(builder, im)
                else
                    _register_translation_linear_form_keys!(builder, im, registered_key_tokens)
                end
            end
            _add_zero_constraint_trusted!(
                builder,
                im,
                origin_factory(_translation_zero_component_label(label, :imag)),
            )
        end
    else
        rekeyed = _rekey_indexed_linear_form(form, index_to_key)
        _add_translation_zero_linear_form!(
            builder,
            rekeyed;
            phase_atol,
            label,
            registered_key_tokens,
            origin_factory,
            register_keys,
            adjoint_key_cache,
        )
    end
    return nothing
end

function _add_translation_zero_linear_form!(
    builder::MomentLinearBuilder{K,C,M},
    form::LinearMomentForm{K,FC};
    phase_atol,
    label,
    registered_key_tokens=nothing,
    origin_factory=_translation_zero_linear_origin,
    register_keys::Bool=true,
    adjoint_key_cache=nothing,
) where {K,C,FC,M}
    if C <: Real
        re = _real_part_form_coefficients(form, C; atol=phase_atol)
        if !isempty(re)
            if register_keys
                if registered_key_tokens === nothing
                    _register_translation_linear_form_keys!(builder, re)
                else
                    _register_translation_linear_form_keys!(builder, re, registered_key_tokens)
                end
            end
            _add_zero_constraint_trusted!(
                builder,
                re,
                origin_factory(_translation_zero_component_label(label, :real)),
            )
        end
        im = _imag_part_form_coefficients(form, C; atol=phase_atol)
        if !isempty(im)
            if register_keys
                if registered_key_tokens === nothing
                    _register_translation_linear_form_keys!(builder, im)
                else
                    _register_translation_linear_form_keys!(builder, im, registered_key_tokens)
                end
            end
            _add_zero_constraint_trusted!(
                builder,
                im,
                origin_factory(_translation_zero_component_label(label, :imag)),
            )
        end
    elseif !isempty(form)
        re, im = _translation_complex_zero_linear_parts!(
            builder,
            form;
            registered_key_tokens,
            register_keys,
            adjoint_key_cache,
        )
        if !isempty(re)
            _add_zero_constraint_trusted!(
                builder,
                re,
                origin_factory(_translation_zero_component_label(label, :complex));
                trusted_self_adjoint=true,
            )
        end
        if !isempty(im)
            _add_zero_constraint_trusted!(
                builder,
                im,
                origin_factory(_translation_zero_component_label(label, :complex));
                trusted_self_adjoint=true,
            )
        end
    end
    return nothing
end

function _add_translation_zero_polynomial_linear_form!(
    builder::MomentLinearBuilder{K,C,M},
    poly::Polynomial{PauliAlgebra,T},
    label;
    registered_key_tokens=nothing,
) where {K,C,T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    form = _linearize_moment_polynomial(K, C, poly)
    isempty(form) && return nothing
    if registered_key_tokens === nothing
        _register_translation_linear_form_keys!(builder, form)
    else
        _register_translation_linear_form_keys!(builder, form, registered_key_tokens)
    end
    _add_zero_constraint_trusted!(builder, form, TranslationZeroOrigin(label, 1, 1, :scalar))
    return nothing
end

function _add_translation_scalar_equality_linear_constraints!(
    builder::MomentLinearBuilder{K,LC,M},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
) where {
    K,
    LC,
    T<:Unsigned,
    PC<:Number,
    P<:Polynomial{PauliAlgebra,T,PC},
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T},
}
    for (idx, poly) in pairs(pop.eq_constraints)
        context = "Equality constraint $idx"
        reduced = _translation_reduced_constraint_polynomial(poly, n_sites, MP; context)
        _check_polynomial_moments_covered(reduced, moment_basis, context)
        _add_translation_zero_polynomial_linear_form!(
            builder,
            reduced,
            (feature=:scalar_equality, index=idx),
        )
    end
    return nothing
end

function _add_qmbcertify_scalar_equality_linear_constraints!(
    builder::MomentLinearBuilder{K,LC,M},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
) where {
    K,
    LC,
    T<:Unsigned,
    PC<:Number,
    P<:Polynomial{PauliAlgebra,T,PC},
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T},
}
    for (idx, poly) in pairs(pop.eq_constraints)
        context = "Equality constraint $idx"
        reduced = _qmbcertify_reduce_chain_polynomial(poly, Int(n_sites), MP)
        _check_polynomial_moments_covered(reduced, moment_basis, context)
        _add_translation_zero_polynomial_linear_form!(
            builder,
            reduced,
            (feature=:scalar_equality, index=idx),
        )
    end
    return nothing
end

function _add_translation_scalar_inequality_linear_blocks!(
    builder::MomentLinearBuilder{K,LC,M},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
    constraint_idx::Int,
    stage_times_ns=nothing,
) where {
    K,
    LC,
    T<:Unsigned,
    PC<:Number,
    P<:Polynomial{PauliAlgebra,T,PC},
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T,LC},
}
    cone = LC <: Real ? :PSD : :HPSD
    for (idx, poly) in pairs(pop.ineq_constraints)
        context = "Inequality constraint $idx"
        reduced = _translation_reduced_constraint_polynomial(poly, n_sites, MP; context)
        _check_polynomial_moments_covered(reduced, moment_basis, context)
        form = _linearize_moment_polynomial(K, LC, reduced)
        _register_translation_linear_form_keys!(builder, form)
        label = (feature=:scalar_inequality, index=idx)
        logical_rows = Any[one(M)]
        constraint_idx += 1
        add_psd_block!(
            builder,
            cone,
            reshape([form], 1, 1),
            _translation_block_meta(M, cone, 1, label, logical_rows, nothing);
            constraint_idx,
        )
        push!(logical_block_sizes, 1)
        push!(block_sizes, 1)
        push!(block_labels, label)
        push!(block_logical_row_labels, logical_rows)
        push!(block_transforms, nothing)
    end
    return constraint_idx
end

function _add_qmbcertify_scalar_inequality_linear_blocks!(
    builder::MomentLinearBuilder{K,LC,M},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
    constraint_idx::Int,
) where {
    K,
    LC,
    T<:Unsigned,
    PC<:Number,
    P<:Polynomial{PauliAlgebra,T,PC},
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T,LC},
}
    cone = LC <: Real ? :PSD : :HPSD
    for (idx, poly) in pairs(pop.ineq_constraints)
        context = "Inequality constraint $idx"
        reduced = _qmbcertify_reduce_chain_polynomial(poly, Int(n_sites), MP)
        _check_polynomial_moments_covered(reduced, moment_basis, context)
        form = _linearize_moment_polynomial(K, LC, reduced)
        _register_translation_linear_form_keys!(builder, form)
        label = (feature=:scalar_inequality, index=idx)
        logical_rows = Any[one(M)]
        constraint_idx += 1
        add_psd_block!(
            builder,
            cone,
            reshape([form], 1, 1),
            _translation_block_meta(M, cone, 1, label, logical_rows, nothing);
            constraint_idx,
        )
        push!(logical_block_sizes, 1)
        push!(block_sizes, 1)
        push!(block_labels, label)
        push!(block_logical_row_labels, logical_rows)
        push!(block_transforms, nothing)
    end
    return constraint_idx
end

function _add_translation_moment_eq_linear_constraints!(
    builder::MomentLinearBuilder{K,LC,M},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    orbit_reps::Vector{NormalMonomial{PauliAlgebra,T}},
    sign_symmetry::Bool,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
) where {
    K,
    LC,
    T<:Unsigned,
    PC<:Number,
    P<:Polynomial{PauliAlgebra,T,PC},
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T},
}
    isempty(pop.moment_eq_constraints) && return nothing

    row_bases, row_degrees = _translation_moment_eq_row_bases(
        orbit_reps;
        sign_symmetry,
    )
    one_mono = one(M)
    work_coeff_type = LC <: Real ? Complex{LC} : LC
    buf = T[]

    for (idx, g) in pairs(pop.moment_eq_constraints)
        rows = _truncate_moment_eq_row_bases(row_bases, row_degrees, g)
        isempty(rows) && continue

        for (row_idx, row_mono) in pairs(rows)
            terms = Tuple{work_coeff_type,NormalMonomial{PauliAlgebra,T}}[]
            sizehint!(terms, length(row_mono) * length(g.terms))
            for (c_row, row_word) in row_mono
                conj_row = _conj_coef(PauliAlgebra, c_row)
                for (coef, mono) in g.terms
                    _push_scaled_buffered_terms!(
                        terms,
                        conj_row * coef,
                        PauliAlgebra,
                        simplify!(PauliAlgebra, _neat_dot3!(buf, row_word, mono, one_mono)),
                        T,
                        work_coeff_type,
                    )
                end
            end

            poly = _polynomial_from_owned_terms!(terms)
            iszero(poly) && continue
            context = "Moment equality constraint $idx"
            reduced = _translation_reduced_constraint_polynomial(poly, n_sites, MP; context)
            _check_polynomial_moments_covered(reduced, moment_basis, context)
            _add_translation_zero_polynomial_linear_form!(
                builder,
                reduced,
                (feature=:moment_equality, index=idx, row=row_idx, row_monomial=row_mono),
            )
        end
    end
    return nothing
end

function _translation_axis_rotation_equality_key(
    source::M,
    coefficient::Integer,
    target::M,
) where {M<:NormalMonomial{PauliAlgebra}}
    terms = Dict{M,Int}()
    terms[source] = 1
    terms[target] = get(terms, target, 0) - Int(coefficient)
    pairs = sort!(
        [(mono, coef) for (mono, coef) in terms if !iszero(coef)];
        by=first,
    )
    isempty(pairs) && return ()
    if last(first(pairs)) < 0
        pairs = [(mono, -coef) for (mono, coef) in pairs]
    end
    return Tuple((coef, mono) for (mono, coef) in pairs)
end

function _translation_axis_rotation_equality_rows(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    n_sites::Integer,
    moment_basis::Vector{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    basis_set = Set(moment_basis)
    seen = Set{Any}()
    rows = NamedTuple[]
    for (generator_idx, rotation) in pairs(pauli_global_axis_rotation_generators(ops))
        for source in moment_basis
            coef, image = _act_monomial(rotation, source)
            target = _translation_orbit_representative(image, n_sites)
            target in basis_set || throw(ArgumentError(
                "axis_rotation_equalities=true requires the translation moment basis " *
                "to be closed under global Pauli-axis rotations; generator $generator_idx " *
                "maps $source to $target, which is missing."
            ))
            coef == 1 && target == source && continue
            key = _translation_axis_rotation_equality_key(source, coef, target)
            key in seen && continue
            push!(seen, key)
            push!(rows, (
                generator=generator_idx,
                source=source,
                coefficient=coef,
                target=target,
            ))
        end
    end
    return rows
end

function _translation_axis_rotation_equality_polynomial(
    row,
    ::Type{MP},
) where {T<:Unsigned,CMP<:Number,MP<:Polynomial{PauliAlgebra,T,CMP}}
    M = NormalMonomial{PauliAlgebra,T}
    return MP(Tuple{CMP,M}[
        (one(CMP), row.source),
        (-convert(CMP, row.coefficient), row.target),
    ])
end

function _add_translation_axis_rotation_linear_constraints!(
    builder::MomentLinearBuilder{K,LC,M},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    n_sites::Integer,
    moment_basis::Vector{M},
    ::Type{MP},
) where {
    K,
    LC,
    T<:Unsigned,
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T},
}
    rows = _translation_axis_rotation_equality_rows(ops, n_sites, moment_basis)
    for (row_idx, row) in pairs(rows)
        poly = _translation_axis_rotation_equality_polynomial(row, MP)
        iszero(poly) && continue
        _add_translation_zero_polynomial_linear_form!(
            builder,
            poly,
            (
                feature=:axis_rotation_equality,
                generator=row.generator,
                row=row_idx,
                source=row.source,
                target=row.target,
                coefficient=row.coefficient,
            ),
        )
    end
    return nothing
end

function _add_translation_objective_terms!(
    builder::MomentLinearBuilder{K,C,M},
    objective::Polynomial{PauliAlgebra,T,PC},
) where {K,C,T<:Unsigned,PC<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(objective.terms))
    for (coef, mono) in simplify(objective)
        converted = convert(C, coef)
        iszero(converted) && continue
        key = register_moment!(builder, mono)
        push!(pairs, key => converted)
    end
    add_objective_terms!(builder, pairs)
    return nothing
end

function _add_translation_objective_terms!(
    builder::MomentLinearBuilder{K,C,M},
    objective::Polynomial{PauliAlgebra,T,PC},
    quotient,
) where {K,C,T<:Unsigned,PC<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    quotient === nothing && return _add_translation_objective_terms!(builder, objective)
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(objective.terms))
    for (coef, mono) in simplify(objective)
        converted = convert(C, coef)
        iszero(converted) && continue
        key = _moment_key(K, mono)
        info = _axis_quotient_info(quotient, key)
        info.forced_zero && continue
        stored = _register_builder_moment_key!(
            builder.key_to_monomial,
            K,
            info.key,
            info.mono,
        )
        push!(pairs, stored => converted * convert(C, info.sign))
    end
    add_objective_terms!(builder, pairs)
    return nothing
end

function _add_translation_linear_state_opt_linear_constraints!(
    builder::MomentLinearBuilder{K,LC,M},
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    n_sites::Integer,
    test_width::Int,
    sign_symmetry::Bool,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    allow_free_moments::Bool=false,
    ;
    registered_key_tokens=nothing,
    linear_state_opt_mode::Symbol=:contiguous,
) where {
    K,
    LC,
    T<:Unsigned,
    H<:Number,
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T},
}
    test_width == 0 && return nothing

    mode = _normalize_linear_state_opt_mode(linear_state_opt_mode)
    seen_rows = mode == :qmbcertify ? Set{MP}() : nothing
    moment_set = (!allow_free_moments || mode == :qmbcertify) ?
        Set(moment_basis) :
        nothing
    for test_mono in _linear_state_opt_tests(M, Int(n_sites), test_width; sign_symmetry, mode)
        row = im * (hamiltonian * test_mono - test_mono * hamiltonian)
        iszero(row) && continue
        context = "Linear state-opt width $test_width row"
        reduced = _translation_reduced_constraint_polynomial(row, n_sites, MP; context)
        iszero(reduced) && continue
        if seen_rows !== nothing
            reduced in seen_rows && continue
            push!(seen_rows, reduced)
        end
        if moment_set !== nothing
            covered = all(mono -> mono in moment_set, monomials(reduced))
            if !covered
                mode == :qmbcertify && continue
                _check_polynomial_moments_covered(reduced, moment_basis, context)
            end
        end
        _add_translation_zero_polynomial_linear_form!(
            builder,
            reduced,
            (
                feature=:linear_state_opt,
                width=test_width,
                mode=mode,
                test_monomial=test_mono,
            );
            registered_key_tokens,
        )
    end
    return nothing
end

function _check_linear_entries_moments_covered(
    entries::AbstractMatrix{LinearMomentForm{K,C}},
    moment_basis::Vector{M},
    context::AbstractString,
) where {K,C,M<:NormalMonomial}
    moment_set = Set(symmetric_canon(expval(mono)) for mono in moment_basis)
    for form in entries
        for (key, _) in form
            key in moment_set || throw(ArgumentError(
                "$context references moment key $(repr(key)) outside the current moment basis. " *
                "Increase the half-basis order or use contiguous_rdm_support=:extend."
            ))
        end
    end
    return nothing
end

function _add_translation_su2_extend_rdm_linear_blocks!(
    builder::MomentLinearBuilder{K,LC,M},
    n_sites::Integer,
    rdm_k::Integer,
    ::Type{MP_C},
    real_moment_matrix::Bool,
    phase_atol::R,
    registered_key_tokens::Set{Any},
    quotient,
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
    constraint_idx::Int,
    stage_times_ns=nothing,
) where {
    K,
    LC,
    MP_C<:Number,
    R<:Real,
    T<:Unsigned,
    M<:NormalMonomial{PauliAlgebra,T},
}
    stage_start_ns = time_ns()
    n = Int(n_sites)
    schur_transform, states = _pauli_su2_schur_matrix(rdm_k)
    schur_rows = _pauli_sparse_transform_rows(
        transpose(schur_transform);
        atol=R(phase_atol),
    )
    columns = _pauli_su2_state_columns(states)
    blocks = pauli_su2_rdm_blocks(rdm_k)
    orbit_monos, orbit_indices = _contiguous_rdm_reduced_orbit_data(n, rdm_k, M)
    code_actions = _pauli_rdm_code_actions(MP_C, rdm_k)
    cone = real_moment_matrix ? :PSD : :HPSD
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_rdm_setup] =
            get(stage_times_ns, :su2_extend_rdm_setup, 0) + Int(time_ns() - stage_start_ns)
    end

    for block in blocks
        stage_start_ns = time_ns()
        reference_m2 = -block.spin2
        reference_rows = [
            columns[(block.spin2, reference_m2, mult)]
            for mult in 1:block.multiplicity
        ]
        reference_schur_rows = _select_sparse_transform_rows(
            schur_rows,
            reference_rows,
        )
        if stage_times_ns !== nothing
            stage_times_ns[:su2_extend_rdm_block_setup] =
                get(stage_times_ns, :su2_extend_rdm_block_setup, 0) +
                Int(time_ns() - stage_start_ns)
        end

        stage_start_ns = time_ns()
        complex_block = _translation_contiguous_rdm_su2_reduced_linear_block(
            orbit_monos,
            orbit_indices,
            rdm_k,
            reference_schur_rows,
            code_actions,
            K,
            MP_C;
            atol=R(phase_atol),
            stage_times_ns,
        )
        if stage_times_ns !== nothing
            stage_times_ns[:su2_extend_rdm_block_build] =
                get(stage_times_ns, :su2_extend_rdm_block_build, 0) +
                Int(time_ns() - stage_start_ns)
        end

        stage_start_ns = time_ns()
        entries = real_moment_matrix ?
            _realify_hermitian_linear_block(complex_block, LC; atol=R(phase_atol)) :
            complex_block
        if stage_times_ns !== nothing
            stage_times_ns[:su2_extend_rdm_realify] =
                get(stage_times_ns, :su2_extend_rdm_realify, 0) + Int(time_ns() - stage_start_ns)
        end

        stage_start_ns = time_ns()
        label = (
            feature=:contiguous_rdm,
            k=Int(rdm_k),
            decomposition=:su2,
            spin2=block.spin2,
        )
        logical_rows = _contiguous_rdm_su2_row_labels(rdm_k, block)
        transform = TranslationSU2RDMTransform(rdm_k, block, schur_transform)
        entries = _prepare_translation_linear_entries!(
            builder,
            entries,
            registered_key_tokens,
            quotient,
        )
        if stage_times_ns !== nothing
            stage_times_ns[:su2_extend_rdm_prepare] =
                get(stage_times_ns, :su2_extend_rdm_prepare, 0) + Int(time_ns() - stage_start_ns)
        end

        stage_start_ns = time_ns()
        constraint_idx += 1
        add_psd_block!(
            builder,
            cone,
            entries,
            _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
            constraint_idx,
        )
        push!(logical_block_sizes, size(complex_block, 1))
        push!(block_sizes, size(entries, 1))
        push!(block_labels, label)
        push!(block_logical_row_labels, logical_rows)
        push!(block_transforms, transform)
        if stage_times_ns !== nothing
            stage_times_ns[:su2_extend_rdm_append] =
                get(stage_times_ns, :su2_extend_rdm_append, 0) + Int(time_ns() - stage_start_ns)
        end
    end
    return constraint_idx
end

function _translation_contiguous_rdm_linear_block(
    n_sites::Integer,
    k::Integer,
    ::Type{K},
    ::Type{C},
) where {K,C<:Number}
    kk = Int(k)
    1 <= kk <= Int(n_sites) || throw(DomainError(k, "`contiguous_rdm_k` must satisfy 1 <= k <= n_sites."))

    T = eltype(K)
    M = NormalMonomial{PauliAlgebra,T}
    dim = 1 << kk
    R = typeof(real(one(C)))
    scale = C(inv(R(2)^kk))
    reduced_by_code = _contiguous_rdm_reduced_monomials(n_sites, kk, M)

    entries = Matrix{LinearMomentForm{K,C}}(undef, dim, dim)
    for col_state in 0:(dim - 1), row_state in 0:(dim - 1)
        entries[row_state + 1, col_state + 1] = _translation_terms_to_linear_form(
            K,
            C,
            _contiguous_rdm_entry_terms(reduced_by_code, kk, row_state, col_state, C, scale),
        )
    end
    return entries
end

function _translation_qmbcertify_rdm_linear_block(
    n_sites::Integer,
    k::Integer,
    rows::Vector{Int},
    ::Type{K},
    ::Type{C},
) where {K,C<:Number}
    return only(_translation_qmbcertify_rdm_linear_blocks(n_sites, k, [rows], K, C))
end

function _translation_qmbcertify_rdm_linear_blocks(
    n_sites::Integer,
    k::Integer,
    row_blocks::AbstractVector{<:AbstractVector{<:Integer}},
    ::Type{K},
    ::Type{C},
) where {K,C<:Number}
    kk = Int(k)
    1 <= kk <= Int(n_sites) || throw(DomainError(k, "`contiguous_rdm_k` must satisfy 1 <= k <= n_sites."))

    T = eltype(K)
    M = NormalMonomial{PauliAlgebra,T}
    dim = 1 << kk

    row_to_block = zeros(Int, dim)
    row_to_local = zeros(Int, dim)
    active_rows = Tuple{Int,Int,Int}[]
    sizehint!(active_rows, sum(length, row_blocks; init=0))
    entry_terms = Vector{Matrix{Vector{Pair{K,C}}}}(undef, length(row_blocks))

    for (block_idx, rows) in pairs(row_blocks)
        isempty(rows) && throw(ArgumentError("QMBCertify RDM row block must be non-empty."))
        block_dim = length(rows)
        entry_terms[block_idx] = [Pair{K,C}[] for _ in 1:block_dim, _ in 1:block_dim]
        for (local_idx, row_value) in pairs(rows)
            row = Int(row_value)
            1 <= row <= dim || throw(ArgumentError(
                "QMBCertify RDM row block $block_idx contains row $row outside 1:$dim for k=$kk."
            ))
            if !iszero(row_to_block[row])
                row_to_block[row] == block_idx && throw(ArgumentError(
                    "QMBCertify RDM row block $block_idx contains duplicate row $row."
                ))
                throw(ArgumentError(
                    "QMBCertify RDM row $row appears in both block $(row_to_block[row]) and block $block_idx."
                ))
            end
            row_to_block[row] = block_idx
            row_to_local[row] = local_idx
            push!(active_rows, (block_idx, local_idx, row))
        end
    end

    R = typeof(real(one(C)))
    scale = C(inv(R(2)^kk))

    for code0 in 0:(4^kk - 1)
        word = _qmbcertify_even_axis_word(kk, code0)
        word === nothing && continue
        mono = isempty(word) ? one(M) : M(T.(word))
        reduced_key = _moment_key(K, _translation_orbit_representative(mono, n_sites))
        for (block_idx, row_local, row) in active_rows
            col0, sign = _qmbcertify_rdm_entry_sign(row - 1, code0, kk)
            col = col0 + 1
            row_to_block[col] == block_idx || continue
            push!(entry_terms[block_idx][row_local, row_to_local[col]], reduced_key => scale * convert(C, sign))
        end
    end

    blocks = Vector{Matrix{LinearMomentForm{K,C}}}(undef, length(row_blocks))
    for (block_idx, terms) in pairs(entry_terms)
        entries = Matrix{LinearMomentForm{K,C}}(undef, size(terms))
        for idx in eachindex(terms)
            entries[idx] = _linear_moment_form_from_owned_pairs!(terms[idx])
        end
        blocks[block_idx] = entries
    end
    return blocks
end

function _translation_qmbcertify_rdm_addons(
    n_sites::Integer,
    rdm_ks::AbstractVector{<:Integer},
    ::Type{K},
    ::Type{C},
) where {K,C<:Number}
    T = eltype(K)
    M = NormalMonomial{PauliAlgebra,T}
    row_blocks_by_k = Dict{Int,Vector{Vector{Int}}}()
    linear_blocks_by_k = Dict{Int,Vector{Matrix{LinearMomentForm{K,C}}}}()
    moment_basis = M[]

    for rdm_k in rdm_ks
        kk = Int(rdm_k)
        row_blocks = _qmbcertify_rdm_block_templates(kk; ambient_sites=n_sites)
        linear_blocks = _translation_qmbcertify_rdm_linear_blocks(
            n_sites,
            kk,
            row_blocks,
            K,
            C,
        )
        row_blocks_by_k[kk] = row_blocks
        linear_blocks_by_k[kk] = linear_blocks
        for entries in linear_blocks
            append!(moment_basis, _linear_entries_moment_basis(M, entries))
        end
    end

    return (
        row_blocks_by_k=row_blocks_by_k,
        linear_blocks_by_k=linear_blocks_by_k,
        moment_basis=sorted_unique!(moment_basis),
    )
end

function _qmbcertify_chain_rdm_linear_blocks(
    n_sites::Integer,
    k::Integer,
    row_blocks::AbstractVector{<:AbstractVector{<:Integer}},
    ::Type{K},
    ::Type{C},
) where {K,C<:Number}
    kk = Int(k)
    1 <= kk <= Int(n_sites) || throw(DomainError(k, "`contiguous_rdm_k` must satisfy 1 <= k <= n_sites."))

    T = eltype(K)
    M = NormalMonomial{PauliAlgebra,T}
    dim = 1 << kk

    row_to_block = zeros(Int, dim)
    row_to_local = zeros(Int, dim)
    active_rows = Tuple{Int,Int,Int}[]
    sizehint!(active_rows, sum(length, row_blocks; init=0))
    entry_terms = Vector{Matrix{Vector{Pair{K,C}}}}(undef, length(row_blocks))

    for (block_idx, rows) in pairs(row_blocks)
        isempty(rows) && throw(ArgumentError("QMBCertify RDM row block must be non-empty."))
        block_dim = length(rows)
        entry_terms[block_idx] = [Pair{K,C}[] for _ in 1:block_dim, _ in 1:block_dim]
        for (local_idx, row_value) in pairs(rows)
            row = Int(row_value)
            1 <= row <= dim || throw(ArgumentError(
                "QMBCertify RDM row block $block_idx contains row $row outside 1:$dim for k=$kk."
            ))
            if !iszero(row_to_block[row])
                row_to_block[row] == block_idx && throw(ArgumentError(
                    "QMBCertify RDM row block $block_idx contains duplicate row $row."
                ))
                throw(ArgumentError(
                    "QMBCertify RDM row $row appears in both block $(row_to_block[row]) and block $block_idx."
                ))
            end
            row_to_block[row] = block_idx
            row_to_local[row] = local_idx
            push!(active_rows, (block_idx, local_idx, row))
        end
    end

    R = typeof(real(one(C)))
    scale = C(inv(R(2)^kk))

    for code0 in 0:(4^kk - 1)
        word = _qmbcertify_even_axis_word(kk, code0)
        word === nothing && continue
        reduced = _qmbcertify_chain_support_monomial(M, word, Int(n_sites))
        reduced === nothing && continue
        reduced_key = _moment_key(K, reduced)
        for (block_idx, row_local, row) in active_rows
            col0, sign = _qmbcertify_rdm_entry_sign(row - 1, code0, kk)
            col = col0 + 1
            row_to_block[col] == block_idx || continue
            push!(
                entry_terms[block_idx][row_local, row_to_local[col]],
                reduced_key => scale * convert(C, sign),
            )
        end
    end

    blocks = Vector{Matrix{LinearMomentForm{K,C}}}(undef, length(row_blocks))
    for (block_idx, terms) in pairs(entry_terms)
        entries = Matrix{LinearMomentForm{K,C}}(undef, size(terms))
        for idx in eachindex(terms)
            entries[idx] = _linear_moment_form_from_owned_pairs!(terms[idx])
        end
        blocks[block_idx] = entries
    end
    return blocks
end

function _qmbcertify_chain_psd_state_opt_rows(qmb_basis, family_rows_by_parity, width::Int)
    rows_by_parity = [eltype(first(family_rows_by_parity))[] for _ in 1:2]
    labels_by_parity = [Any[] for _ in 1:2]
    family_count_by_parity = zeros(Int, 2)
    for record in qmb_basis.family_records
        record.word_length <= width || continue
        parity_index = record.parity + 1
        row = family_rows_by_parity[parity_index][record.family]
        family_count_by_parity[parity_index] += 1
        push!(rows_by_parity[parity_index], row)
        push!(
            labels_by_parity[parity_index],
            (
                feature=:qmbcertify_psd_state_opt_source_family,
                parity=record.parity,
                family=record.family,
                monomial=row,
            ),
        )
    end
    return (
        rows_by_parity=rows_by_parity,
        labels_by_parity=labels_by_parity,
        family_count_by_parity=family_count_by_parity,
    )
end

function _qmbcertify_chain_psd_state_opt_block_records(
    family_counts::AbstractVector{<:Integer},
    n_sites::Int,
)
    iseven(n_sites) || throw(ArgumentError(
        "QMBCertify chain PSD state-opt block accounting currently requires even n_sites; got $n_sites."
    ))
    length(family_counts) == 2 || throw(ArgumentError(
        "QMBCertify chain PSD state-opt block accounting expects two parity family counts."
    ))
    half = div(n_sites, 2)
    records = NamedTuple[]
    for parity_index in 1:2
        parity = parity_index - 1
        family_count = Int(family_counts[parity_index])
        iszero(family_count) && continue
        for momentum in 0:half
            realified = !(momentum == 0 || momentum == half)
            push!(
                records,
                (
                    parity=parity,
                    momentum=momentum,
                    block_size=realified ? 2family_count : family_count,
                    family_count=family_count,
                    realified=realified,
                ),
            )
        end
    end
    return records
end

function _qmbcertify_chain_psd_state_opt_linear_entry(
    row::M,
    col::M,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    k::Int,
    n::Int,
    col_translates::Vector{M},
    ::Type{K},
    ::Type{C},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    R = typeof(real(one(C)))
    terms = Tuple{C,M}[]
    sizehint!(terms, length(col_translates) * max(1, length(hamiltonian.terms)))

    for r in 0:(n - 1)
        phase = _momentum_phase(R, k, r, n)
        translated_col = col_translates[r + 1]
        entry = row * hamiltonian * translated_col -
            (1 // 2) * (hamiltonian * row * translated_col + row * translated_col * hamiltonian)
        for (coef, mono) in entry
            iszero(coef) && continue
            reduced = _qmbcertify_chain_support_monomial(M, mono.word, n)
            reduced === nothing && continue
            push!(terms, (convert(C, phase * coef), reduced))
        end
    end
    return _translation_terms_to_linear_form(K, C, terms)
end

function _qmbcertify_chain_psd_state_opt_term_cache(
    block_basis::Vector{M},
    n::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    translated::Dict{M,Vector{M}},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    Coef = promote_type(H, Rational{Int})
    entries = Matrix{Vector{Tuple{Int,Coef,M}}}(undef, length(block_basis), length(block_basis))
    for (col_idx, col) in pairs(block_basis), (row_idx, row) in pairs(block_basis)
        col_translates = translated[col]
        terms = Tuple{Int,Coef,M}[]
        sizehint!(terms, length(col_translates) * max(1, length(hamiltonian.terms)))
        for r in 0:(n - 1)
            translated_col = col_translates[r + 1]
            entry = row * hamiltonian * translated_col -
                (1 // 2) * (hamiltonian * row * translated_col + row * translated_col * hamiltonian)
            for (coef, mono) in entry
                iszero(coef) && continue
                reduced = _qmbcertify_chain_support_monomial(M, mono.word, n)
                reduced === nothing && continue
                push!(terms, (r, convert(Coef, coef), reduced))
            end
        end
        entries[row_idx, col_idx] = terms
    end
    return entries
end

function _qmbcertify_chain_psd_state_opt_linear_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    translated::Dict{M,Vector{M}},
    ::Type{K},
    ::Type{C},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    entries = Matrix{LinearMomentForm{K,C}}(undef, length(block_basis), length(block_basis))
    for (col_idx, col) in pairs(block_basis), (row_idx, row) in pairs(block_basis)
        entries[row_idx, col_idx] = _qmbcertify_chain_psd_state_opt_linear_entry(
            row,
            col,
            hamiltonian,
            k,
            n,
            translated[col],
            K,
            C,
        )
    end
    return entries
end

function _qmbcertify_chain_psd_state_opt_linear_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    term_cache::AbstractMatrix{<:Vector{Tuple{Int,H,M}}},
    ::Type{K},
    ::Type{C},
) where {H<:Number,M<:NormalMonomial{PauliAlgebra},K,C<:Number}
    return _translation_psd_state_opt_linear_block(
        block_basis,
        k,
        n,
        term_cache,
        K,
        C,
    )
end

function _translation_psd_state_opt_linear_entry(
    row::M,
    col::M,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    k::Int,
    n::Int,
    col_translates::Vector{M},
    rep_cache::Dict{M,M},
    ::Type{K},
    ::Type{C},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    R = typeof(real(one(C)))
    terms = Tuple{C,M}[]
    sizehint!(terms, length(col_translates) * max(1, length(hamiltonian.terms)))

    for r in 0:(n - 1)
        phase = _momentum_phase(R, k, r, n)
        translated_col = col_translates[r + 1]
        entry = row * hamiltonian * translated_col -
            (1 // 2) * (hamiltonian * row * translated_col + row * translated_col * hamiltonian)
        for (coef, mono) in entry
            iszero(coef) && continue
            rep = get!(rep_cache, mono) do
                _translation_orbit_representative(mono, n)
            end
            push!(terms, (convert(C, phase * coef), rep))
        end
    end
    return _translation_terms_to_linear_form(K, C, terms)
end

function _translation_psd_state_opt_term_cache(
    block_basis::Vector{M},
    n::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    translated::Dict{M,Vector{M}},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T}}
    Coef = promote_type(H, Rational{Int})
    entries = Matrix{Vector{Tuple{Int,Coef,M}}}(undef, length(block_basis), length(block_basis))
    rep_cache = Dict{M,M}()
    for (col_idx, col) in pairs(block_basis), (row_idx, row) in pairs(block_basis)
        col_translates = translated[col]
        terms = Tuple{Int,Coef,M}[]
        sizehint!(terms, length(col_translates) * max(1, length(hamiltonian.terms)))
        for r in 0:(n - 1)
            translated_col = col_translates[r + 1]
            entry = row * hamiltonian * translated_col -
                (1 // 2) * (hamiltonian * row * translated_col + row * translated_col * hamiltonian)
            for (coef, mono) in entry
                iszero(coef) && continue
                rep = get!(rep_cache, mono) do
                    _translation_orbit_representative(mono, n)
                end
                push!(terms, (r, convert(Coef, coef), rep))
            end
        end
        entries[row_idx, col_idx] = terms
    end
    return entries
end

function _translation_psd_state_opt_required_moments(
    term_cache::AbstractMatrix{<:Vector{Tuple{Int,H,M}}},
) where {H<:Number,M<:NormalMonomial{PauliAlgebra}}
    required = Set{M}()
    for terms in term_cache
        for (_, coef, rep) in terms
            iszero(coef) && continue
            push!(required, rep)
        end
    end
    return required
end

function _translation_psd_state_opt_linear_entry(
    terms::Vector{Tuple{Int,H,M}},
    k::Int,
    n::Int,
    ::Type{K},
    ::Type{C},
) where {H<:Number,M<:NormalMonomial{PauliAlgebra},K,C<:Number}
    R = typeof(real(one(C)))
    phased_terms = Tuple{C,M}[]
    sizehint!(phased_terms, length(terms))
    for (r, coef, rep) in terms
        phase = _momentum_phase(R, k, r, n)
        push!(phased_terms, (convert(C, phase * coef), rep))
    end
    return _translation_terms_to_linear_form(K, C, phased_terms)
end

function _translation_psd_state_opt_linear_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{K},
    ::Type{C},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    entries = Matrix{LinearMomentForm{K,C}}(undef, length(block_basis), length(block_basis))
    for (col_idx, col) in pairs(block_basis), (row_idx, row) in pairs(block_basis)
        entries[row_idx, col_idx] = _translation_psd_state_opt_linear_entry(
            row,
            col,
            hamiltonian,
            k,
            n,
            translated[col],
            rep_cache,
            K,
            C,
        )
    end
    return entries
end

function _translation_psd_state_opt_linear_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    term_cache::AbstractMatrix{<:Vector{Tuple{Int,H,M}}},
    ::Type{K},
    ::Type{C},
) where {H<:Number,M<:NormalMonomial{PauliAlgebra},K,C<:Number}
    entries = Matrix{LinearMomentForm{K,C}}(undef, length(block_basis), length(block_basis))
    for col_idx in eachindex(block_basis), row_idx in eachindex(block_basis)
        entries[row_idx, col_idx] = _translation_psd_state_opt_linear_entry(
            term_cache[row_idx, col_idx],
            k,
            n,
            K,
            C,
        )
    end
    return entries
end

function _subtract_linear_forms(
    left::LinearMomentForm{K,C},
    right::LinearMomentForm{K,C};
    atol,
) where {K,C}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(left) + length(right))
    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(atol) : atol
    left_terms = left.terms
    right_terms = right.terms
    left_idx = firstindex(left_terms)
    right_idx = firstindex(right_terms)
    while left_idx <= lastindex(left_terms) && right_idx <= lastindex(right_terms)
        left_term = left_terms[left_idx]
        right_term = right_terms[right_idx]
        if key_lt(left_term.first, right_term.first)
            abs(left_term.second) <= tolerance || push!(pairs, left_term)
            left_idx += 1
        elseif key_lt(right_term.first, left_term.first)
            value = -right_term.second
            abs(value) <= tolerance || push!(pairs, right_term.first => value)
            right_idx += 1
        else
            left_active = abs(left_term.second) > tolerance
            right_active = abs(right_term.second) > tolerance
            value = (left_active ? left_term.second : zero(C)) -
                (right_active ? right_term.second : zero(C))
            iszero(value) || push!(pairs, left_term.first => value)
            left_idx += 1
            right_idx += 1
        end
    end
    while left_idx <= lastindex(left_terms)
        term = left_terms[left_idx]
        abs(term.second) <= tolerance || push!(pairs, term)
        left_idx += 1
    end
    while right_idx <= lastindex(right_terms)
        term = right_terms[right_idx]
        value = -term.second
        abs(value) <= tolerance || push!(pairs, term.first => value)
        right_idx += 1
    end
    return LinearMomentForm{K,C}(pairs, Val(:trusted))
end

function _translation_block_meta(
    ::Type{M},
    cone::Symbol,
    block_size::Integer,
    label,
    logical_row_labels::Vector{Any},
    transform,
) where {M<:NormalMonomial{PauliAlgebra}}
    return BlockMeta{M}(
        cone,
        TranslationInvariantBlockOrigin(label, logical_row_labels; transform),
        fill(one(M), Int(block_size)),
    )
end

function _axis_isotypic_row_labels(transform, block_idx::Integer)
    return Any[
        (
            feature=:axis_isotypic_row,
            irrep=transform.irrep_labels[block_idx],
            irrep_dimension=transform.irrep_dimensions[block_idx],
            irrep_multiplicity=transform.irrep_multiplicities[block_idx],
            row=row_idx,
        )
        for row_idx in 1:transform.isotypic_block_sizes[block_idx]
    ]
end

function _translation_axis_isotypic_adapted_blocks(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    block_basis::Vector{M},
    order::Integer;
    atol::Real,
) where {M<:NormalMonomial{PauliAlgebra}}
    transform = pauli_axis_rotation_isotypic_transform(
        ops,
        order;
        basis=block_basis,
        atol,
        basis_is_orbit_representatives=true,
    )
    blocks = NamedTuple[]
    for block_idx in eachindex(transform.transform_blocks)
        block = transform.transform_blocks[block_idx]
        block_size = size(block, 2)
        iszero(block_size) && continue
        push!(
            blocks,
            (
                axis_irrep=transform.irrep_labels[block_idx],
                irrep_dimension=transform.irrep_dimensions[block_idx],
                irrep_multiplicity=transform.irrep_multiplicities[block_idx],
                axis_group_order=transform.axis_group_order,
                row_basis=Matrix(transpose(block)),
                row_labels=_axis_isotypic_row_labels(transform, block_idx),
            ),
        )
    end
    return blocks
end

function _axis_restricted_projector(
    projector::AbstractMatrix{<:Real},
    row_basis::AbstractMatrix{<:Number};
    realified::Bool,
)
    if realified
        n = size(projector, 1)
        size(row_basis, 2) == 2n || throw(DimensionMismatch(
            "realified reflection row basis has $(size(row_basis, 2)) columns but axis projector has size $n"
        ))
        real_projector = zeros(Float64, 2n, 2n)
        real_projector[1:n, 1:n] .= projector
        real_projector[(n + 1):(2n), (n + 1):(2n)] .= projector
        return row_basis * real_projector * transpose(row_basis)
    end

    size(row_basis, 2) == size(projector, 1) || throw(DimensionMismatch(
        "reflection row basis has $(size(row_basis, 2)) columns but axis projector has size $(size(projector, 1))"
    ))
    return conj(row_basis) * projector * transpose(row_basis)
end

function _axis_reflection_isotypic_row_labels(
    axis_label::Symbol,
    irrep_dimension::Integer,
    irrep_multiplicity::Integer,
    axis_group_order::Integer,
    reflection::Integer,
    row_count::Integer,
)
    return Any[
        (
            feature=:axis_reflection_isotypic_row,
            axis_irrep=axis_label,
            irrep_dimension=Int(irrep_dimension),
            irrep_multiplicity=Int(irrep_multiplicity),
            axis_group_order=Int(axis_group_order),
            reflection=Int(reflection),
            row=row_idx,
        )
        for row_idx in 1:Int(row_count)
    ]
end

function _restricted_projector_fractional_eigen_count(
    projector::AbstractMatrix{<:Number},
    atol::Real,
)
    mat = Matrix(projector)
    vals = eigvals(Hermitian(0.5 .* (mat .+ mat')))
    tol = max(Float64(atol), 1e-8)
    return count(value -> value > tol && value < 1 - tol, vals)
end

function _reflection_axis_projector_groups(
    projectors,
    row_basis::AbstractMatrix{<:Number};
    realified::Bool,
    atol::Real,
)
    standard_idx = findfirst(==(:standard), projectors.irrep_labels)
    axis_vector_idx = findfirst(==(:axis_vector), projectors.irrep_labels)
    combine_standard_axis = false
    if standard_idx !== nothing && axis_vector_idx !== nothing
        for block_idx in (standard_idx, axis_vector_idx)
            restricted = _axis_restricted_projector(
                projectors.projector_matrices[block_idx],
                row_basis;
                realified,
            )
            if _restricted_projector_fractional_eigen_count(restricted, atol) > 0
                combine_standard_axis = true
                break
            end
        end
    end

    groups = NamedTuple[]
    skip_axis_vector = false
    for block_idx in eachindex(projectors.projector_matrices)
        if combine_standard_axis && block_idx == standard_idx
            push!(
                groups,
                (
                    label=:standard_axis_vector_family,
                    irrep_dimension=1,
                    projector=projectors.projector_matrices[standard_idx] .+
                        projectors.projector_matrices[axis_vector_idx],
                    require_divisible=false,
                ),
            )
            skip_axis_vector = true
            continue
        end
        combine_standard_axis && skip_axis_vector && block_idx == axis_vector_idx &&
            continue
        push!(
            groups,
            (
                label=projectors.irrep_labels[block_idx],
                irrep_dimension=projectors.irrep_dimensions[block_idx],
                projector=projectors.projector_matrices[block_idx],
                require_divisible=true,
            ),
        )
    end

    restricted_groups = NamedTuple[]
    for group in groups
        restricted = _axis_restricted_projector(group.projector, row_basis; realified)
        fractional = _restricted_projector_fractional_eigen_count(restricted, atol)
        if fractional == 0
            push!(
                restricted_groups,
                merge(group, (restricted_projector=restricted, fractional_eigen_count=0)),
            )
            continue
        end

        n = size(projectors.projector_matrices[1], 1)
        identity_projector = Matrix{Float64}(I, n, n)
        identity_restricted = _axis_restricted_projector(
            identity_projector,
            row_basis;
            realified,
        )
        identity_fractional = _restricted_projector_fractional_eigen_count(
            identity_restricted,
            atol,
        )
        return [(
            label=:axis_reflection_mixed_family,
            irrep_dimension=1,
            projector=identity_projector,
            require_divisible=false,
            restricted_projector=identity_restricted,
            fractional_eigen_count=identity_fractional,
        )]
    end
    return restricted_groups
end

function _translation_reflection_axis_isotypic_adapted_blocks(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    block_basis::Vector{M},
    order::Integer,
    reflection::Integer,
    row_basis::AbstractMatrix{<:Number};
    realified::Bool,
    atol::Real,
) where {M<:NormalMonomial{PauliAlgebra}}
    projectors = pauli_axis_rotation_irrep_projectors(
        ops,
        order;
        basis=block_basis,
        atol,
        basis_is_orbit_representatives=true,
    )
    blocks = NamedTuple[]
    for group in _reflection_axis_projector_groups(
        projectors,
        row_basis;
        realified,
        atol,
    )
        group.fractional_eigen_count == 0 || throw(ArgumentError(
            "Reflected Pauli axis family $(group.label) still has " *
            "$(group.fractional_eigen_count) " *
            "fractional restricted-projector eigenvalue(s); the reflected axis " *
            "split is not invariant for this sector."
        ))
        restricted = group.restricted_projector
        range_basis = _pauli_projector_range_basis(restricted, atol)
        block_size = size(range_basis, 2)
        iszero(block_size) && continue
        irrep_dimension = group.irrep_dimension
        if group.require_divisible
            block_size % irrep_dimension == 0 || throw(ArgumentError(
                "Reflected Pauli axis irrep $(group.label) projector rank " *
                "$block_size is not divisible by irrep dimension $irrep_dimension."
            ))
        end
        irrep_multiplicity = group.require_divisible ?
            div(block_size, irrep_dimension) :
            block_size
        row_transform = transpose(range_basis) * row_basis
        push!(
            blocks,
            (
                reflection=Int(reflection),
                axis_irrep=group.label,
                irrep_dimension=irrep_dimension,
                irrep_multiplicity=irrep_multiplicity,
                axis_group_order=projectors.axis_group_order,
                row_basis=row_transform,
                row_labels=_axis_reflection_isotypic_row_labels(
                    group.label,
                    irrep_dimension,
                    irrep_multiplicity,
                    projectors.axis_group_order,
                    reflection,
                    block_size,
                ),
            ),
        )
    end
    return blocks
end

function _qmbcertify_chain_family_representatives(qmb_basis)
    n = qmb_basis.n_sites
    return [
        [
            qmb_basis.basis_by_parity[parity][n * (family - 1) + 1]
            for family in 1:qmb_basis.family_count_by_parity[parity]
        ]
        for parity in 1:2
    ]
end

function _qmbcertify_base_row_labels(record, row_labels::Vector{Any})
    labels = copy(row_labels)
    if record.realified
        real_labels = Any[(part=:real, row=label) for label in labels]
        append!(real_labels, Any[(part=:imag, row=label) for label in labels])
        return real_labels
    end
    return labels
end

_qmbcertify_linear_state_opt_row_sign(coef::Real) = coef < 0 ? -1 : 1

function _qmbcertify_linear_state_opt_row_sign(coef::Complex)
    real_part = real(coef)
    !iszero(real_part) && return real_part < 0 ? -1 : 1
    return imag(coef) < 0 ? -1 : 1
end

function _qmbcertify_linear_state_opt_row_key(poly::Polynomial)
    isempty(poly.terms) && return poly
    coef = first(poly.terms)[1]
    return _qmbcertify_linear_state_opt_row_sign(coef) < 0 ? -poly : poly
end

function _add_qmbcertify_chain_linear_state_opt_constraints!(
    builder::MomentLinearBuilder{K,LC,M},
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    n_sites::Integer,
    test_width::Int,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP};
    registered_key_tokens=nothing,
) where {
    K,
    LC,
    T<:Unsigned,
    H<:Number,
    M<:NormalMonomial{PauliAlgebra,T},
    MP<:Polynomial{PauliAlgebra,T},
}
    test_width == 0 && return nothing

    moment_set = Set(moment_basis)
    seen_rows = Set{MP}()
    for test_mono in _linear_state_opt_tests(
        M,
        Int(n_sites),
        test_width;
        sign_symmetry=true,
        mode=:qmbcertify,
    )
        row = im * (hamiltonian * test_mono - test_mono * hamiltonian)
        iszero(row) && continue
        reduced = _qmbcertify_reduce_chain_polynomial(row, Int(n_sites), MP)
        iszero(reduced) && continue
        all(mono -> mono in moment_set, monomials(reduced)) || continue
        row_key = _qmbcertify_linear_state_opt_row_key(reduced)
        row_key in seen_rows && continue
        push!(seen_rows, row_key)
        _add_translation_zero_polynomial_linear_form!(
            builder,
            reduced,
            (
                feature=:linear_state_opt,
                width=test_width,
                mode=:qmbcertify,
                test_monomial=test_mono,
            );
            registered_key_tokens,
        )
    end
    return nothing
end

function _pauli_qmbcertify_chain_base_linear_relaxation(
    pop::PolyOpt{PauliAlgebra,T,P},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    extra::Integer=0,
    three_type::Tuple{<:Integer,<:Integer}=(1, 1),
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition::Symbol=:qmbcertify,
    contiguous_rdm_support::Symbol=:extend,
    linear_state_opt_width=nothing,
    linear_state_opt_mode=nothing,
    psd_state_opt_width=nothing,
) where {T<:Unsigned,C0<:Number,P<:Polynomial{PauliAlgebra,T,C0}}
    isempty(pop.moment_eq_constraints) || throw(ArgumentError(
        "QMBCertify base construction currently supports base moment PSD blocks, scalar equality/inequality rows, and source-family add-ons; moment equality constraints are not implemented."
    ))

    construction_start_ns = time_ns()
    stage_times_ns = Dict{Symbol,Int}()
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    eltype(σx) == NormalMonomial{PauliAlgebra,T} || throw(ArgumentError(
        "PolyOpt and Pauli chain operators must use the same Pauli integer type; got objective type $T and operator type $(eltype(σx))."
    ))
    _check_pauli_chain_support(pop.objective, n; context="objective")
    _check_translation_scalar_constraints_compatible(
        pop,
        n;
        check_invariance=true,
        sign_symmetry=true,
        reflection_symmetry=false,
        real_moment_matrix,
    )
    setup_done_ns = time_ns()
    stage_times_ns[:setup] = Int(setup_done_ns - construction_start_ns)

    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    rdm_ks = _normalize_contiguous_rdm_k(contiguous_rdm_k)
    rdm_decomposition = _normalize_contiguous_rdm_decomposition(
        contiguous_rdm_decomposition,
        false,
        false,
    )
    rdm_support = _normalize_contiguous_rdm_support(contiguous_rdm_support)
    isempty(rdm_ks) || rdm_decomposition == :qmbcertify || throw(ArgumentError(
        "QMBCertify base construction supports contiguous RDM blocks only with `contiguous_rdm_decomposition=:qmbcertify`."
    ))
    isempty(rdm_ks) || rdm_support == :extend || throw(ArgumentError(
        "QMBCertify base construction supports contiguous RDM blocks only with `contiguous_rdm_support=:extend`."
    ))
    isempty(rdm_ks) || real_moment_matrix || throw(ArgumentError(
        "QMBCertify base construction currently supports contiguous RDM blocks only with `real_moment_matrix=true`."
    ))
    linear_state_width = _normalize_linear_state_opt_width(linear_state_opt_width)
    linear_state_mode = _normalize_linear_state_opt_mode(
        something(linear_state_opt_mode, :qmbcertify),
    )
    psd_state_width = _normalize_psd_state_opt_width(psd_state_opt_width)
    linear_state_width == 0 || linear_state_mode == :qmbcertify || throw(ArgumentError(
        "QMBCertify base construction supports linear state-opt rows only with `linear_state_opt_mode=:qmbcertify`."
    ))
    psd_state_width == 0 || real_moment_matrix || throw(ArgumentError(
        "QMBCertify base construction currently supports PSD state-opt blocks only with `real_moment_matrix=true`."
    ))
    qmb_basis = pauli_qmbcertify_chain_basis(ops, d; extra, three_type)
    family_rows_by_parity = _qmbcertify_chain_family_representatives(qmb_basis)
    one_mono = one(first(σx))
    all_rows = NormalMonomial{PauliAlgebra,T}[one_mono]
    for rows in family_rows_by_parity
        append!(all_rows, rows)
    end
    sorted_unique!(all_rows)
    basis_done_ns = time_ns()
    stage_times_ns[:basis] = Int(basis_done_ns - setup_done_ns)

    MP_R = _pauli_chain_real_coeff_type(C0)
    MP_C = Complex{MP_R}
    LC = real_moment_matrix ? MP_R : MP_C
    MP_P = Polynomial{PauliAlgebra,T,LC}
    objective_mp = _qmbcertify_reduce_chain_polynomial(pop.objective, n, MP_P)
    K = typeof(symmetric_canon(expval(one_mono)))
    M = NormalMonomial{PauliAlgebra,T}
    builder = MomentLinearBuilder(K, LC, M)
    _add_translation_objective_terms!(builder, objective_mp)
    registered_key_tokens = _translation_registered_key_tokens(builder)

    translated = Dict{M,Vector{M}}()
    for row in all_rows
        translated[row] = isone(row) ?
            fill(row, n) :
            [_translate_pauli_monomial(row, r, n) for r in 0:(n - 1)]
    end
    product_cache = _TranslationProductCache(
        Dict{Tuple{M,M},Vector{Tuple{Int,MP_C,M}}}(),
        0,
        0,
    )

    logical_block_sizes = Int[]
    block_sizes = Int[]
    block_labels = Any[]
    block_logical_row_labels = Vector{Any}[]
    block_transforms = Any[]
    constraint_idx = 0
    cone = real_moment_matrix ? :PSD : :HPSD
    atol = MP_R(phase_atol)
    precompute_done_ns = time_ns()
    stage_times_ns[:precompute] = Int(precompute_done_ns - basis_done_ns)

    for record in qmb_basis.base_block_records
        family_rows = family_rows_by_parity[record.parity + 1]
        block_basis = copy(family_rows)
        source_labels = Any[
            (
                feature=:qmbcertify_base_source_family,
                parity=record.parity,
                family=family,
                monomial=family_rows[family],
            )
            for family in eachindex(family_rows)
        ]
        if record.identity_row
            pushfirst!(block_basis, one_mono)
            pushfirst!(source_labels, (feature=:identity, monomial=one_mono))
        end
        complex_entries = _qmbcertify_chain_momentum_block_complex_linear_entries(
            block_basis,
            record.momentum,
            n,
            translated,
            K,
            MP_C,
            product_cache,
        )
        entries = if real_moment_matrix
            record.realified ?
                _realify_hermitian_linear_block_full(complex_entries, LC; atol) :
                _qmbcertify_real_fixed_linear_block(complex_entries, LC; atol)
        else
            complex_entries
        end
        size(entries, 1) == record.block_size || throw(ArgumentError(
            "QMBCertify base block parity=$(record.parity), momentum=$(record.momentum) " *
            "constructed size $(size(entries, 1)) but expected $(record.block_size)."
        ))
        entries = _prepare_translation_linear_entries!(
            builder,
            entries,
            registered_key_tokens,
            nothing,
        )
        label = (
            feature=:qmbcertify_base,
            parity=record.parity,
            momentum=record.momentum,
            identity_row=record.identity_row,
            realified=record.realified,
        )
        logical_rows = _qmbcertify_base_row_labels(record, source_labels)
        transform = TranslationDFTTransform(n, record.momentum, real_moment_matrix)
        constraint_idx += 1
        add_psd_block!(
            builder,
            cone,
            entries,
            _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
            constraint_idx,
        )
        push!(logical_block_sizes, length(logical_rows))
        push!(block_sizes, size(entries, 1))
        push!(block_labels, label)
        push!(block_logical_row_labels, logical_rows)
        push!(block_transforms, transform)
    end
    block_done_ns = time_ns()
    stage_times_ns[:block_assembly] = Int(block_done_ns - precompute_done_ns)

    base_moment_basis = _translation_builder_psd_moment_basis(builder)
    _check_objective_moments_covered(objective_mp, base_moment_basis)
    constraint_start_ns = time_ns()
    _add_qmbcertify_scalar_equality_linear_constraints!(
        builder,
        pop,
        n,
        base_moment_basis,
        MP_P,
    )
    constraint_idx = _add_qmbcertify_scalar_inequality_linear_blocks!(
        builder,
        pop,
        n,
        base_moment_basis,
        MP_P,
        logical_block_sizes,
        block_sizes,
        block_labels,
        block_logical_row_labels,
        block_transforms,
        constraint_idx,
    )
    constraint_done_ns = time_ns()
    stage_times_ns[:constraint_append] = Int(constraint_done_ns - constraint_start_ns)
    linear_state_moment_basis = base_moment_basis
    if !isempty(rdm_ks)
        rdm_start_ns = time_ns()
        for rdm_k in rdm_ks
            row_blocks = _qmbcertify_rdm_block_templates(rdm_k; ambient_sites=n)
            complex_blocks = _qmbcertify_chain_rdm_linear_blocks(
                n,
                rdm_k,
                row_blocks,
                K,
                MP_C,
            )
            for (block_idx, rows) in pairs(row_blocks)
                complex_block = complex_blocks[block_idx]
                entries = _realify_hermitian_linear_block(
                    complex_block,
                    LC;
                    atol=MP_R(phase_atol),
                )
                label = (
                    feature=:contiguous_rdm,
                    k=rdm_k,
                    decomposition=:qmbcertify,
                    block=block_idx,
                )
                logical_rows = _contiguous_rdm_state_labels(rdm_k, rows)
                entries = _prepare_translation_linear_entries!(
                    builder,
                    entries,
                    registered_key_tokens,
                    nothing,
                )
                constraint_idx += 1
                add_psd_block!(
                    builder,
                    cone,
                    entries,
                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, nothing);
                    constraint_idx,
                )
                push!(logical_block_sizes, size(complex_block, 1))
                push!(block_sizes, size(entries, 1))
                push!(block_labels, label)
                push!(block_logical_row_labels, logical_rows)
                push!(block_transforms, nothing)
            end
        end
        stage_times_ns[:contiguous_rdm] = Int(time_ns() - rdm_start_ns)
        linear_state_moment_basis = _translation_builder_psd_moment_basis(builder)
    end
    if psd_state_width > 0
        psd_state_start_ns = time_ns()
        psd_state_rows = _qmbcertify_chain_psd_state_opt_rows(
            qmb_basis,
            family_rows_by_parity,
            psd_state_width,
        )
        psd_state_records = _qmbcertify_chain_psd_state_opt_block_records(
            psd_state_rows.family_count_by_parity,
            n,
        )
        psd_state_term_caches = Vector{Any}(undef, 2)
        for parity_index in 1:2
            block_basis = psd_state_rows.rows_by_parity[parity_index]
            psd_state_term_caches[parity_index] = isempty(block_basis) ?
                nothing :
                _qmbcertify_chain_psd_state_opt_term_cache(
                    block_basis,
                    n,
                    pop.objective,
                    translated,
                )
        end
        for record in psd_state_records
            block_basis = psd_state_rows.rows_by_parity[record.parity + 1]
            term_cache = something(psd_state_term_caches[record.parity + 1])
            complex_entries = _hermitianize_pauli_linear_block(
                _qmbcertify_chain_psd_state_opt_linear_block(
                    block_basis,
                    record.momentum,
                    n,
                    term_cache,
                    K,
                    MP_C,
                ),
            )
            entries = if record.realified
                _realify_hermitian_linear_block_full(complex_entries, LC; atol)
            else
                _qmbcertify_real_fixed_linear_block(complex_entries, LC; atol)
            end
            size(entries, 1) == record.block_size || throw(ArgumentError(
                "QMBCertify PSD state-opt block parity=$(record.parity), momentum=$(record.momentum) " *
                "constructed size $(size(entries, 1)) but expected $(record.block_size)."
            ))
            entries = _prepare_translation_linear_entries!(
                builder,
                entries,
                registered_key_tokens,
                nothing,
            )
            label = (
                feature=:psd_state_opt,
                width=psd_state_width,
                mode=:qmbcertify,
                parity=record.parity,
                momentum=record.momentum,
            )
            logical_rows = psd_state_rows.labels_by_parity[record.parity + 1]
            transform = TranslationDFTTransform(n, record.momentum, real_moment_matrix)
            constraint_idx += 1
            add_psd_block!(
                builder,
                cone,
                entries,
                _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                constraint_idx,
            )
            push!(logical_block_sizes, length(logical_rows))
            push!(block_sizes, size(entries, 1))
            push!(block_labels, label)
            push!(block_logical_row_labels, logical_rows)
            push!(block_transforms, transform)
        end
        stage_times_ns[:psd_state_opt] = Int(time_ns() - psd_state_start_ns)
        linear_state_moment_basis = _translation_builder_psd_moment_basis(builder)
    end
    if linear_state_width > 0
        linear_state_start_ns = time_ns()
        _add_qmbcertify_chain_linear_state_opt_constraints!(
            builder,
            pop.objective,
            n,
            linear_state_width,
            linear_state_moment_basis,
            MP_P;
            registered_key_tokens,
        )
        stage_times_ns[:linear_state_opt] =
            Int(time_ns() - linear_state_start_ns)
    end
    linearization_start_ns = time_ns()
    linear = finalize!(builder; stage_times_ns, stage_prefix=:qmbcertify_base_finalize)
    linearization_done_ns = time_ns()
    stage_times_ns[:linearization] = Int(linearization_done_ns - linearization_start_ns)

    report = TranslationInvariantReport(
        n,
        d,
        qmb_basis.source_row_count + 1,
        sum(qmb_basis.family_count_by_parity) + 1,
        false,
        0,
        Pair{Int,Int}[],
        1.0,
        collect(0:div(n, 2)),
        true,
        logical_block_sizes,
        block_sizes,
        block_labels,
        _translation_block_coefficient_domains(block_transforms),
        _translation_block_exact_coefficient_domains(block_transforms),
        length(base_moment_basis),
        length(linear.moments),
        length(linear.zero_constraints),
        Int(linearization_done_ns - construction_start_ns),
        stage_times_ns,
        real_moment_matrix,
        product_cache.hits,
        product_cache.misses,
        length(product_cache.terms),
        _translation_product_cache_hit_rate(product_cache),
        _translation_zero_report_histogram(linear.zero_constraints),
        false,
        0,
        0,
        0,
        0,
        1.0,
    )
    return linear, report
end

function _pauli_translation_base_linear_relaxation(
    pop::PolyOpt{PauliAlgebra,T,P},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis::Union{Nothing,Vector{NormalMonomial{PauliAlgebra,T}}}=nothing,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    sign_symmetry::Bool=true,
    reflection_symmetry::Bool=false,
    axis_rotation_symmetry::Bool=false,
    check_invariance::Bool=true,
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition::Symbol=:full,
    contiguous_rdm_support::Symbol=:closed,
    u1_symmetry::Bool=false,
    su2_symmetry::Bool=false,
    base_su2_extend_rdm::Bool=false,
    su2_moment_quotient::Bool=false,
    su2_moment_quotient_atol::Real=1e-11,
    su2_moment_quotient_condition_limit::Real=1e10,
    axis_rotation_equalities::Bool=false,
    axis_rotation_quotient::Bool=false,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
    linear_state_opt_width=nothing,
    linear_state_opt_mode::Symbol=:contiguous,
    psd_state_opt_width=nothing,
) where {T<:Unsigned,C0<:Number,P<:Polynomial{PauliAlgebra,T,C0}}
    rdm_ks = _normalize_contiguous_rdm_k(contiguous_rdm_k)
    rdm_decomposition = _normalize_contiguous_rdm_decomposition(
        contiguous_rdm_decomposition,
        u1_symmetry,
        su2_symmetry,
    )
    rdm_support = _normalize_contiguous_rdm_support(contiguous_rdm_support)
    rdm_decomposition in (:full, :u1, :su2, :qmbcertify) || throw(ArgumentError(
        "Direct base translation linear relaxation currently supports contiguous_rdm_decomposition=:full, :u1, :su2, or :qmbcertify."
    ))
    linear_state_width = _normalize_linear_state_opt_width(linear_state_opt_width)
    linear_state_mode = _normalize_linear_state_opt_mode(linear_state_opt_mode)
    psd_state_width = _normalize_psd_state_opt_width(psd_state_opt_width)
    su2_quotient_atol = Float64(su2_moment_quotient_atol)
    su2_quotient_atol >= 0 || throw(DomainError(
        su2_moment_quotient_atol,
        "`su2_moment_quotient_atol` must be non-negative.",
    ))
    su2_quotient_condition_limit = Float64(su2_moment_quotient_condition_limit)
    isfinite(su2_quotient_condition_limit) && su2_quotient_condition_limit > 0 ||
        throw(DomainError(
            su2_moment_quotient_condition_limit,
            "`su2_moment_quotient_condition_limit` must be positive and finite.",
        ))
    su2_moment_quotient && !su2_symmetry && throw(ArgumentError(
        "`su2_moment_quotient=true` requires `su2_symmetry=true`.",
    ))
    su2_moment_quotient && !base_su2_extend_rdm && throw(ArgumentError(
        "`su2_moment_quotient=true` requires `base_su2_extend_rdm=true`.",
    ))
    standalone_axis_reduction = axis_rotation_symmetry && !su2_symmetry
    axis_constraints_requested = axis_rotation_equalities || standalone_axis_reduction
    axis_rotation_quotient && !axis_constraints_requested && throw(ArgumentError(
        "`axis_rotation_quotient=true` requires `axis_rotation_equalities=true` or `axis_rotation_symmetry=true`."
    ))
    axis_rotation_quotient && su2_symmetry && throw(ArgumentError(
        "`axis_rotation_quotient=true` is for the finite-axis direct-linear path; `su2_symmetry=true` already subsumes finite axis rotations."
    ))
    axis_rotation_quotient && !real_moment_matrix && throw(ArgumentError(
        "`axis_rotation_quotient=true` currently requires `real_moment_matrix=true`."
    ))
    if rdm_decomposition == :qmbcertify
        axis_rotation_quotient || throw(ArgumentError(
            "`contiguous_rdm_decomposition=:qmbcertify` is structural-target-only unless " *
            "`axis_rotation_quotient=true` is active on the direct-linear finite-axis path. " *
            "QMBCertify's shared RDM PSD variables require the same translation, " *
            "reflection, and global-axis identifications used by that quotient."
        ))
        reflection_symmetry || throw(ArgumentError(
            "`contiguous_rdm_decomposition=:qmbcertify` requires `reflection_symmetry=true` " *
            "because QMBCertify's RDM row blocks are built from reduce4 translation/reflection classes."
        ))
    end
    axis_constraints_requested && check_invariance &&
        _check_pauli_su2_symmetry_compatible(pop, ops)
    u1_symmetry && _check_pauli_charge_neutral(pop)
    su2_symmetry && check_invariance && _check_pauli_su2_symmetry_compatible(pop, ops)
    linear_state_mode == :qmbcertify && su2_symmetry && throw(ArgumentError(
        "`linear_state_opt_mode=:qmbcertify` is currently supported only by the finite-axis direct-linear path without `su2_symmetry=true`."
    ))

    construction_start_ns = time_ns()
    construction_stage_times_ns = Dict{Symbol,Int}()
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    eltype(σx) == NormalMonomial{PauliAlgebra,T} || throw(ArgumentError(
        "PolyOpt and Pauli chain operators must use the same Pauli integer type; got objective type $T and operator type $(eltype(σx))."
    ))

    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    one_mono = one(first(σx))

    _check_pauli_chain_support(pop.objective, n; context="objective")
    _check_translation_complex_reflection_sectors(
        n,
        momenta;
        reflection_symmetry,
        real_moment_matrix,
    )
    _check_translation_scalar_constraints_compatible(
        pop,
        n;
        check_invariance,
        sign_symmetry,
        reflection_symmetry,
        real_moment_matrix,
    )
    setup_done_ns = time_ns()
    construction_stage_times_ns[:setup] = Int(setup_done_ns - construction_start_ns)

    local_basis_size = 0
    orbit_reps = NormalMonomial{PauliAlgebra,T}[]
    closure_basis = nothing
    if isnothing(basis) && d < n
        local_basis_size = _pauli_contiguous_chain_basis_size_hint(n, d; periodic=true)
        orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    else
        local_basis = isnothing(basis) ? pauli_contiguous_chain_basis(ops, d; periodic=true) : basis
        one_mono in local_basis || throw(ArgumentError("Translation-invariant Pauli basis must include the identity."))
        _check_pauli_chain_support(local_basis, n; context="basis")
        _check_translation_basis_closure(local_basis, n)
        local_basis_size = length(local_basis)
        orbit_reps = _translation_orbit_representatives(local_basis, n)
        closure_basis = isnothing(basis) ? nothing : local_basis
    end
    _check_pauli_translation_constraint_closure(
        closure_basis,
        d,
        n,
        pop.objective,
        rdm_ks,
        rdm_support,
        linear_state_width,
        psd_state_width;
        sign_symmetry,
    )
    basis_done_ns = time_ns()
    construction_stage_times_ns[:basis] = Int(basis_done_ns - setup_done_ns)

    axis_orbit_summary = _pauli_axis_orbit_summary(ops, orbit_reps, d; strict=false)
    axis_done_ns = time_ns()
    construction_stage_times_ns[:axis_diagnostics] = Int(axis_done_ns - basis_done_ns)

    nontrivial_reps = [m for m in orbit_reps if !isone(m)]
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    real_moment_matrix || 0 in sectors || throw(ArgumentError(
        "Momentum sector 0 is required because it carries the normalized identity moment."
    ))

    MP_R = _pauli_chain_real_coeff_type(C0)
    base_su2_extend_rdm && !(su2_symmetry && rdm_support == :extend && rdm_decomposition == :su2) &&
        throw(ArgumentError(
            "`base_su2_extend_rdm=true` requires `su2_symmetry=true`, " *
            "`contiguous_rdm_decomposition=:su2`, and `contiguous_rdm_support=:extend`."
        ))
    base_su2_extend_rdm && !real_moment_matrix && reflection_symmetry && throw(ArgumentError(
        "`base_su2_extend_rdm=true` with `real_moment_matrix=false` does not " *
        "yet support `reflection_symmetry=true`."
    ))
    direct_su2_extend_rdm =
        base_su2_extend_rdm && su2_symmetry && rdm_support == :extend && rdm_decomposition == :su2
    use_su2_base_reduction =
        su2_symmetry && (isempty(rdm_ks) || rdm_support == :closed || direct_su2_extend_rdm)
    if use_su2_base_reduction
        use_direct_su2_base_linear =
            direct_su2_extend_rdm &&
            isempty(pop.eq_constraints) &&
            isempty(pop.ineq_constraints) &&
            isempty(pop.moment_eq_constraints) &&
            !axis_rotation_equalities &&
            (real_moment_matrix || !reflection_symmetry)
        if use_direct_su2_base_linear
            direct_base = _pauli_su2_translation_base_linear_relaxation(
                pop.objective,
                orbit_reps,
                n,
                sectors;
                real_moment_matrix,
                reflection_symmetry,
                phase_atol=MP_R(phase_atol),
                singlet_channel_equalities=
                    singlet_channel_equalities && !su2_moment_quotient,
                singlet_channel_atol,
                emit_invariance_rows=!su2_moment_quotient,
            )
            su2_done_ns = time_ns()
            construction_stage_times_ns[:su2_linear_base] = Int(su2_done_ns - axis_done_ns)
            for (stage, elapsed_ns) in direct_base.stage_times_ns
                construction_stage_times_ns[stage] = elapsed_ns
            end
            base_report = _translation_su2_base_report(
                direct_base.linear,
                n,
                d,
                local_basis_size,
                orbit_reps,
                axis_orbit_summary,
                sectors,
                sign_symmetry,
                Int(su2_done_ns - construction_start_ns),
                construction_stage_times_ns,
                real_moment_matrix,
                direct_base.product_cache_hits,
                direct_base.product_cache_misses,
                direct_base.product_cache_entries,
                direct_base.product_cache_hit_rate,
            )
            extended = _append_translation_direct_extend_addons_to_su2_linear(
                pop,
                direct_base.linear,
                base_report,
                n;
                sign_symmetry,
                real_moment_matrix,
                linear_state_width,
                linear_state_opt_mode=linear_state_mode,
                psd_state_width,
                rdm_ks,
                sectors,
                phase_atol=MP_R(phase_atol),
                stage_times_ns=construction_stage_times_ns,
            )
            extended_done_ns = time_ns()
            construction_stage_times_ns[:su2_extend_addons] = Int(extended_done_ns - su2_done_ns)
            report = _translation_report_with_linear(
                base_report,
                extended.linear,
                extended.logical_block_sizes,
                extended.block_sizes,
                extended.block_labels,
                extended.block_transforms,
                Int(extended_done_ns - construction_start_ns),
                construction_stage_times_ns,
            )
            quotiented = _apply_translation_su2_moment_quotient(
                extended.linear,
                report,
                n;
                enabled=su2_moment_quotient,
                allow_registered_projection=isnothing(basis),
                atol=su2_quotient_atol,
                condition_limit=su2_quotient_condition_limit,
                construction_start_ns,
                stage_times_ns=construction_stage_times_ns,
            )
            return quotiented.linear, quotiented.report
        end

        use_wigner_su2_base = !reflection_symmetry && 0 in sectors
        su2_base_mp = if use_wigner_su2_base
            pauli_su2_translation_orbit_wigner_eckart_moment_problem(
                pop.objective,
                orbit_reps,
                n;
                assume_su2_invariant=!check_invariance,
                ops,
                momenta=sectors,
                phase_atol,
                singlet_channel_equalities=
                    singlet_channel_equalities && !su2_moment_quotient,
                singlet_channel_atol,
            )
        else
            pauli_su2_translation_orbit_moment_problem(
                pop.objective,
                orbit_reps,
                n;
                assume_su2_invariant=!check_invariance,
                ops,
                real_moment_matrix,
                momenta,
                reflection_symmetry,
                phase_atol,
                singlet_channel_equalities=
                    singlet_channel_equalities && !su2_moment_quotient,
                singlet_channel_atol,
            )
        end
        su2_mp = if real_moment_matrix && any(cone == :HPSD for (cone, _) in su2_base_mp.constraints)
            MP_R = _pauli_chain_real_coeff_type(C0)
            MP_P = Polynomial{PauliAlgebra,T,MP_R}
            _realify_pauli_su2_translation_orbit_moment_problem(
                su2_base_mp,
                MP_P;
                phase_atol=MP_R(phase_atol),
            )
        else
            su2_base_mp
        end
        append_rdm_ks = direct_su2_extend_rdm ? Int[] : rdm_ks
        append_linear_state_width = direct_su2_extend_rdm ? 0 : linear_state_width
        append_psd_state_width = direct_su2_extend_rdm ? 0 : psd_state_width
        su2_mp = _append_translation_extra_constraints_to_su2_base_moment_problem(
            pop,
            su2_mp,
            orbit_reps,
            n;
            real_moment_matrix,
            sign_symmetry,
            linear_state_width=append_linear_state_width,
            axis_rotation_equalities,
            rdm_ks=append_rdm_ks,
            rdm_decomposition,
            rdm_support,
            psd_state_width=append_psd_state_width,
            sectors,
            phase_atol,
            ops,
        )
        su2_done_ns = time_ns()
        construction_stage_times_ns[:su2_moment_problem] = Int(su2_done_ns - axis_done_ns)
        base_report = _translation_su2_base_report(
            su2_mp,
            n,
            d,
            local_basis_size,
            orbit_reps,
            axis_orbit_summary,
            sectors,
            sign_symmetry,
            Int(su2_done_ns - construction_start_ns),
            construction_stage_times_ns,
            real_moment_matrix,
        )
        if direct_su2_extend_rdm
            extended = _append_translation_direct_extend_addons_to_su2_linear(
                pop,
                su2_mp.linear,
                base_report,
                n;
                sign_symmetry,
                real_moment_matrix,
                linear_state_width,
                linear_state_opt_mode=linear_state_mode,
                psd_state_width,
                rdm_ks,
                sectors,
                phase_atol=MP_R(phase_atol),
                stage_times_ns=construction_stage_times_ns,
            )
            extended_done_ns = time_ns()
            construction_stage_times_ns[:su2_extend_addons] = Int(extended_done_ns - su2_done_ns)
            report = _translation_report_with_linear(
                base_report,
                extended.linear,
                extended.logical_block_sizes,
                extended.block_sizes,
                extended.block_labels,
                extended.block_transforms,
                Int(extended_done_ns - construction_start_ns),
                construction_stage_times_ns,
            )
            quotiented = _apply_translation_su2_moment_quotient(
                extended.linear,
                report,
                n;
                enabled=su2_moment_quotient,
                allow_registered_projection=isnothing(basis),
                atol=su2_quotient_atol,
                condition_limit=su2_quotient_condition_limit,
                construction_start_ns,
                stage_times_ns=construction_stage_times_ns,
            )
            return quotiented.linear, quotiented.report
        end
        report = base_report
        return su2_mp.linear, report
    end

    MP_R = _pauli_chain_real_coeff_type(C0)
    MP_C = Complex{MP_R}
    LC = real_moment_matrix ? MP_R : MP_C
    MP_P = Polynomial{PauliAlgebra,T,LC}
    objective_mp = convert(MP_P, _translation_orbit_reduce_polynomial(pop.objective, n))

    translated = Dict{NormalMonomial{PauliAlgebra,T},Vector{NormalMonomial{PauliAlgebra,T}}}()
    for rep in nontrivial_reps
        translated[rep] = [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
    end
    translated[one_mono] = fill(one_mono, n)

    rep_cache = Dict{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}}()
    product_cache = _TranslationProductCache(
        Dict{
            Tuple{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}},
            Vector{Tuple{Int,MP_C,NormalMonomial{PauliAlgebra,T}}},
        }(),
        0,
        0,
    )

    M = NormalMonomial{PauliAlgebra,T}
    axis_equality_moment_basis = M[]
    identity = symmetric_canon(expval(one(M)))
    K = typeof(identity)
    qmbcertify_rdm_addons = rdm_decomposition == :qmbcertify && !isempty(rdm_ks) ?
        _translation_qmbcertify_rdm_addons(n, rdm_ks, K, MP_C) :
        nothing
    psd_state_opt_basis = psd_state_width > 0 ?
        _contiguous_state_opt_tests(M, n, psd_state_width; sign_symmetry) :
        M[]
    psd_state_opt_translated = Dict{M,Vector{M}}()
    for rep in psd_state_opt_basis
        psd_state_opt_translated[rep] = [
            _translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)
        ]
    end
    psd_state_opt_term_cache = isempty(psd_state_opt_basis) ?
        nothing :
        _translation_psd_state_opt_term_cache(
            psd_state_opt_basis,
            n,
            pop.objective,
            psd_state_opt_translated,
        )
    psd_state_opt_moment_basis = if axis_rotation_quotient && psd_state_opt_term_cache !== nothing
        sorted_unique!(collect(_translation_psd_state_opt_required_moments(
            something(psd_state_opt_term_cache),
        )))
    else
        M[]
    end
    extend_addon_support = rdm_support == :extend && !isempty(rdm_ks)
    can_pre_axis_quotient = axis_rotation_quotient &&
        isnothing(basis) &&
        d < n &&
        isempty(pop.eq_constraints) &&
        isempty(pop.ineq_constraints) &&
        isempty(pop.moment_eq_constraints) &&
        iszero(linear_state_width)
    use_pre_axis_quotient = can_pre_axis_quotient &&
        (isempty(rdm_ks) || qmbcertify_rdm_addons !== nothing)
    pre_axis_quotient_base_basis = use_pre_axis_quotient ?
        _pauli_translation_target_base_moment_basis(
            orbit_reps,
            n;
            sign_symmetry=sign_symmetry && !standalone_axis_reduction,
            real_moment_matrix,
            momenta,
        ) :
        nothing
    pre_axis_quotient_basis = if use_pre_axis_quotient
        quotient_basis = copy(something(pre_axis_quotient_base_basis))
        qmbcertify_rdm_addons === nothing ||
            append!(quotient_basis, qmbcertify_rdm_addons.moment_basis)
        append!(quotient_basis, psd_state_opt_moment_basis)
        sorted_unique!(quotient_basis)
    else
        nothing
    end
    pre_axis_quotient = use_pre_axis_quotient ?
        _translation_axis_rotation_quotient_map(
            K,
            M,
            ops,
            n,
            pre_axis_quotient_basis,
        ) :
        nothing
    entry_axis_quotient = use_pre_axis_quotient ? pre_axis_quotient : nothing
    prepare_axis_quotient = use_pre_axis_quotient ? nothing : pre_axis_quotient
    addon_axis_quotient = use_pre_axis_quotient ? pre_axis_quotient : prepare_axis_quotient
    rdm_axis_quotient = addon_axis_quotient
    rdm_quotient_form_cache = rdm_axis_quotient === nothing ?
        nothing :
        Dict{_LinearMomentFormCacheKey{K,LC},LinearMomentForm{K,LC}}()
    builder = MomentLinearBuilder(K, LC, M)
    _add_translation_objective_terms!(builder, objective_mp, pre_axis_quotient)
    registered_key_tokens = _translation_registered_key_tokens(builder)

    logical_block_sizes = Int[]
    block_sizes = Int[]
    block_labels = Any[]
    block_logical_row_labels = Vector{Any}[]
    block_transforms = Any[]
    precompute_done_ns = time_ns()
    construction_stage_times_ns[:precompute] = Int(precompute_done_ns - axis_done_ns)

    constraint_idx = 0
    cone = real_moment_matrix ? :PSD : :HPSD
    for k in sectors
        sector_basis = k == 0 ? orbit_reps : nontrivial_reps
        blocks = sign_symmetry && !standalone_axis_reduction ?
            _pauli_signature_blocks(sector_basis) :
            [(:all, sector_basis)]
        for (signature, block_basis) in blocks
            isempty(block_basis) && continue
            if reflection_symmetry
                if _translation_reflection_fixed_momentum(k, n)
                    complex_entries = _translation_momentum_block_linear_entries(
                        block_basis,
                        k,
                        n,
                        translated,
                        rep_cache,
                        K,
                        MP_C,
                        product_cache;
                        real_moment_matrix=false,
                        atol=MP_R(phase_atol),
                        quotient=entry_axis_quotient,
                    )
                    standalone_axis_reduction && !use_pre_axis_quotient &&
                        append!(axis_equality_moment_basis, _linear_entries_moment_basis(M, complex_entries))
                    for adapted in _translation_reflection_adapted_blocks(block_basis, k, n; atol=Float64(phase_atol))
                        if standalone_axis_reduction
                            for axis_adapted in _translation_reflection_axis_isotypic_adapted_blocks(
                                ops,
                                block_basis,
                                d,
                                adapted.reflection,
                                adapted.row_basis;
                                realified=false,
                                atol=Float64(phase_atol),
                            )
                                adapted_complex_entries = _hermitianize_pauli_linear_block(
                                    _transform_linear_block(
                                        complex_entries,
                                        axis_adapted.row_basis,
                                        axis_adapted.row_basis;
                                        atol=Float64(phase_atol),
                                    ),
                                )
                                entries = real_moment_matrix ?
                                    _realify_hermitian_linear_block(adapted_complex_entries, LC; atol=MP_R(phase_atol)) :
                                    adapted_complex_entries
                                label = (
                                    momentum=k,
                                    signature=signature,
                                    reflection=axis_adapted.reflection,
                                    axis_irrep=axis_adapted.axis_irrep,
                                    axis_irrep_dimension=axis_adapted.irrep_dimension,
                                    axis_irrep_multiplicity=axis_adapted.irrep_multiplicity,
                                    axis_group_order=axis_adapted.axis_group_order,
                                )
                                logical_rows = axis_adapted.row_labels
                                transform = TranslationAxisIrrepTransform(
                                    TranslationReflectionTransform(
                                        n,
                                        k,
                                        axis_adapted.reflection,
                                        real_moment_matrix,
                                        false,
                                    ),
                                    axis_adapted.axis_irrep,
                                    axis_adapted.irrep_dimension,
                                    axis_adapted.irrep_multiplicity,
                                    axis_adapted.axis_group_order,
                                )
                                entries = _prepare_translation_linear_entries!(
                                    builder,
                                    entries,
                                    registered_key_tokens,
                                    prepare_axis_quotient,
                                )
                                constraint_idx += 1
                                add_psd_block!(
                                    builder,
                                    cone,
                                    entries,
                                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                                    constraint_idx,
                                )
                                push!(logical_block_sizes, size(adapted_complex_entries, 1))
                                push!(block_sizes, size(entries, 1))
                                push!(block_labels, label)
                                push!(block_logical_row_labels, logical_rows)
                                push!(block_transforms, transform)
                            end
                            continue
                        end
                        adapted_complex_entries = _hermitianize_pauli_linear_block(
                            _transform_linear_block(
                                complex_entries,
                                adapted.row_basis_sparse,
                                adapted.row_basis_sparse;
                                atol=Float64(phase_atol),
                            ),
                        )
                        entries = real_moment_matrix ?
                            _realify_hermitian_linear_block(adapted_complex_entries, LC; atol=MP_R(phase_atol)) :
                            adapted_complex_entries
                        label = (momentum=k, signature=signature, reflection=adapted.reflection)
                        logical_rows = adapted.row_labels
                        transform = _translation_transform_descriptor(n, k, real_moment_matrix, label)
                        entries = _prepare_translation_linear_entries!(
                            builder,
                            entries,
                            registered_key_tokens,
                            prepare_axis_quotient,
                        )
                        constraint_idx += 1
                        add_psd_block!(
                            builder,
                            cone,
                            entries,
                            _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                            constraint_idx,
                        )
                        push!(logical_block_sizes, length(logical_rows))
                        push!(block_sizes, size(entries, 1))
                        push!(block_labels, label)
                        push!(block_logical_row_labels, logical_rows)
                        push!(block_transforms, transform)
                    end
                else
                    real_moment_matrix || throw(ArgumentError(
                        "Direct base translation linear relaxation supports non-fixed reflection sectors only for real_moment_matrix=true; got k=$k, n=$n."
                    ))
                    real_entries = _translation_momentum_block_linear_entries(
                        block_basis,
                        k,
                        n,
                        translated,
                        rep_cache,
                        K,
                        LC,
                        product_cache;
                        real_moment_matrix=true,
                        atol=MP_R(phase_atol),
                        quotient=entry_axis_quotient,
                        full_realification=true,
                    )
                    standalone_axis_reduction && !use_pre_axis_quotient &&
                        append!(axis_equality_moment_basis, _linear_entries_moment_basis(M, real_entries))
                    for adapted in _translation_real_reflection_adapted_blocks(block_basis, k, n; atol=Float64(phase_atol))
                        if standalone_axis_reduction
                            for axis_adapted in _translation_reflection_axis_isotypic_adapted_blocks(
                                ops,
                                block_basis,
                                d,
                                adapted.reflection,
                                adapted.row_basis;
                                realified=true,
                                atol=Float64(phase_atol),
                            )
                                entries = _symmetrize_real_linear_block(
                                    _transform_linear_block(
                                        real_entries,
                                        axis_adapted.row_basis,
                                        axis_adapted.row_basis;
                                        atol=Float64(phase_atol),
                                    ),
                                )
                                label = (
                                    momentum=k,
                                    signature=signature,
                                    reflection=axis_adapted.reflection,
                                    axis_irrep=axis_adapted.axis_irrep,
                                    axis_irrep_dimension=axis_adapted.irrep_dimension,
                                    axis_irrep_multiplicity=axis_adapted.irrep_multiplicity,
                                    axis_group_order=axis_adapted.axis_group_order,
                                )
                                logical_rows = axis_adapted.row_labels
                                transform = TranslationAxisIrrepTransform(
                                    TranslationReflectionTransform(
                                        n,
                                        k,
                                        axis_adapted.reflection,
                                        real_moment_matrix,
                                        true,
                                    ),
                                    axis_adapted.axis_irrep,
                                    axis_adapted.irrep_dimension,
                                    axis_adapted.irrep_multiplicity,
                                    axis_adapted.axis_group_order,
                                )
                                entries = _prepare_translation_linear_entries!(
                                    builder,
                                    entries,
                                    registered_key_tokens,
                                    prepare_axis_quotient,
                                )
                                constraint_idx += 1
                                add_psd_block!(
                                    builder,
                                    cone,
                                    entries,
                                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                                    constraint_idx,
                                )
                                push!(logical_block_sizes, size(entries, 1))
                                push!(block_sizes, size(entries, 1))
                                push!(block_labels, label)
                                push!(block_logical_row_labels, logical_rows)
                                push!(block_transforms, transform)
                            end
                            continue
                        end
                        entries = _symmetrize_real_linear_block(
                            _transform_linear_block(
                                real_entries,
                                adapted.row_basis_sparse,
                                adapted.row_basis_sparse;
                                atol=Float64(phase_atol),
                            ),
                        )
                        label = (momentum=k, signature=signature, reflection=adapted.reflection)
                        logical_rows = adapted.row_labels
                        transform = _translation_transform_descriptor(n, k, real_moment_matrix, label)
                        entries = _prepare_translation_linear_entries!(
                            builder,
                            entries,
                            registered_key_tokens,
                            prepare_axis_quotient,
                        )
                        constraint_idx += 1
                        add_psd_block!(
                            builder,
                            cone,
                            entries,
                            _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                            constraint_idx,
                        )
                        push!(logical_block_sizes, size(entries, 1))
                        push!(block_sizes, size(entries, 1))
                        push!(block_labels, label)
                        push!(block_logical_row_labels, logical_rows)
                        push!(block_transforms, transform)
                    end
                end
            elseif standalone_axis_reduction
                complex_entries = _translation_momentum_block_linear_entries(
                    block_basis,
                    k,
                    n,
                    translated,
                    rep_cache,
                    K,
                    MP_C,
                    product_cache;
                    real_moment_matrix=false,
                    atol=MP_R(phase_atol),
                    quotient=entry_axis_quotient,
                )
                !use_pre_axis_quotient &&
                    append!(axis_equality_moment_basis, _linear_entries_moment_basis(M, complex_entries))
                for adapted in _translation_axis_isotypic_adapted_blocks(
                    ops,
                    block_basis,
                    d;
                    atol=Float64(phase_atol),
                )
                    adapted_complex_entries = _hermitianize_pauli_linear_block(
                        _transform_linear_block(
                            complex_entries,
                            adapted.row_basis,
                            adapted.row_basis;
                            atol=Float64(phase_atol),
                        ),
                    )
                    entries = real_moment_matrix ?
                        _realify_hermitian_linear_block(adapted_complex_entries, LC; atol=MP_R(phase_atol)) :
                        adapted_complex_entries
                    label = (
                        momentum=k,
                        signature=signature,
                        axis_irrep=adapted.axis_irrep,
                        axis_irrep_dimension=adapted.irrep_dimension,
                        axis_irrep_multiplicity=adapted.irrep_multiplicity,
                        axis_group_order=adapted.axis_group_order,
                    )
                    logical_rows = adapted.row_labels
                    transform = TranslationAxisIrrepTransform(
                        TranslationDFTTransform(n, k, real_moment_matrix),
                        adapted.axis_irrep,
                        adapted.irrep_dimension,
                        adapted.irrep_multiplicity,
                        adapted.axis_group_order,
                    )
                    entries = _prepare_translation_linear_entries!(
                        builder,
                        entries,
                        registered_key_tokens,
                        prepare_axis_quotient,
                    )
                    constraint_idx += 1
                    add_psd_block!(
                        builder,
                        cone,
                        entries,
                        _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                        constraint_idx,
                    )
                    push!(logical_block_sizes, size(adapted_complex_entries, 1))
                    push!(block_sizes, size(entries, 1))
                    push!(block_labels, label)
                    push!(block_logical_row_labels, logical_rows)
                    push!(block_transforms, transform)
                end
            else
                entries = _translation_momentum_block_linear_entries(
                    block_basis,
                    k,
                    n,
                    translated,
                    rep_cache,
                    K,
                    LC,
                    product_cache;
                    real_moment_matrix,
                    atol=MP_R(phase_atol),
                    quotient=entry_axis_quotient,
                )
                label = (momentum=k, signature=signature)
                logical_rows = Any[mono for mono in block_basis]
                transform = _translation_transform_descriptor(n, k, real_moment_matrix, label)
                entries = _prepare_translation_linear_entries!(
                    builder,
                    entries,
                    registered_key_tokens,
                    prepare_axis_quotient,
                )
                constraint_idx += 1
                add_psd_block!(
                    builder,
                    cone,
                    entries,
                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                    constraint_idx,
                )
                push!(logical_block_sizes, length(block_basis))
                push!(block_sizes, size(entries, 1))
                push!(block_labels, label)
                push!(block_logical_row_labels, logical_rows)
                push!(block_transforms, transform)
            end
        end
    end

    base_moment_basis = use_pre_axis_quotient ?
        something(pre_axis_quotient_base_basis) :
        _translation_builder_psd_moment_basis(builder)
    base_moment_count = use_pre_axis_quotient ?
        length(base_moment_basis) :
        length(base_moment_basis)
    constraint_moment_basis = standalone_axis_reduction ?
        sorted_unique!(vcat(copy(base_moment_basis), axis_equality_moment_basis)) :
        base_moment_basis
    _check_objective_moments_covered(objective_mp, base_moment_basis)
    axis_rotation_quotient_stats = if use_pre_axis_quotient
        pre_axis_quotient
    elseif axis_constraints_requested
        _translation_axis_rotation_quotient_stats(ops, n, constraint_moment_basis)
    else
        _EMPTY_AXIS_ROTATION_QUOTIENT_STATS
    end
    if axis_constraints_requested && !axis_rotation_quotient
        _add_translation_axis_rotation_linear_constraints!(
            builder,
            ops,
            n,
            constraint_moment_basis,
            MP_P,
        )
    end
    _add_translation_scalar_equality_linear_constraints!(
        builder,
        pop,
        n,
        constraint_moment_basis,
        MP_P,
    )
    constraint_idx = _add_translation_scalar_inequality_linear_blocks!(
        builder,
        pop,
        n,
        constraint_moment_basis,
        MP_P,
        logical_block_sizes,
        block_sizes,
        block_labels,
        block_logical_row_labels,
        block_transforms,
        constraint_idx,
    )
    _add_translation_moment_eq_linear_constraints!(
        builder,
        pop,
        n,
        orbit_reps,
        sign_symmetry,
        constraint_moment_basis,
        MP_P,
    )
    linear_state_moment_basis = if linear_state_mode == :qmbcertify && extend_addon_support
        addon_basis = copy(constraint_moment_basis)
        qmbcertify_rdm_addons === nothing ||
            append!(addon_basis, qmbcertify_rdm_addons.moment_basis)
        append!(addon_basis, psd_state_opt_moment_basis)
        sorted_unique!(addon_basis)
    else
        constraint_moment_basis
    end
    linear_state_allow_free_moments =
        extend_addon_support && linear_state_mode != :qmbcertify
    _add_translation_linear_state_opt_linear_constraints!(
        builder,
        pop.objective,
        n,
        linear_state_width,
        sign_symmetry,
        linear_state_moment_basis,
        MP_P,
        linear_state_allow_free_moments;
        registered_key_tokens,
        linear_state_opt_mode=linear_state_mode,
    )
    if psd_state_width > 0
        state_opt_basis = psd_state_opt_basis
        if !isempty(state_opt_basis)
            state_opt_term_cache = something(psd_state_opt_term_cache)
            for k in sectors
                complex_entries = _hermitianize_pauli_linear_block(
                    _translation_psd_state_opt_linear_block(
                        state_opt_basis,
                        k,
                        n,
                        state_opt_term_cache,
                        K,
                        MP_C,
                    ),
                )
                if !extend_addon_support
                    _check_linear_entries_moments_covered(
                        complex_entries,
                        constraint_moment_basis,
                        "PSD state-opt width=$psd_state_width momentum=$k block",
                    )
                end
                entries = real_moment_matrix ?
                    _realify_hermitian_linear_block(complex_entries, LC; atol=MP_R(phase_atol)) :
                    complex_entries
                label = (feature=:psd_state_opt, width=psd_state_width, momentum=k)
                logical_rows = Any[mono for mono in state_opt_basis]
                transform = TranslationDFTTransform(n, k, real_moment_matrix)
                entries = _prepare_translation_linear_entries!(
                    builder,
                    entries,
                    registered_key_tokens,
                    addon_axis_quotient,
                )
                constraint_idx += 1
                add_psd_block!(
                    builder,
                    cone,
                    entries,
                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                    constraint_idx,
                )
                push!(logical_block_sizes, size(complex_entries, 1))
                push!(block_sizes, size(entries, 1))
                push!(block_labels, label)
                push!(block_logical_row_labels, logical_rows)
                push!(block_transforms, transform)
            end
        end
    end

    if !isempty(rdm_ks)
        for rdm_k in rdm_ks
            if rdm_decomposition == :qmbcertify
                qmbcertify_rdm_blocks = qmbcertify_rdm_addons.row_blocks_by_k[Int(rdm_k)]
                qmbcertify_linear_blocks =
                    qmbcertify_rdm_addons.linear_blocks_by_k[Int(rdm_k)]
                for (block_idx, rows) in pairs(qmbcertify_rdm_blocks)
                    complex_block = qmbcertify_linear_blocks[block_idx]
                    if rdm_support == :closed
                        _check_linear_entries_moments_covered(
                            complex_block,
                            constraint_moment_basis,
                            "QMBCertify contiguous RDM k=$rdm_k block $block_idx",
                        )
                    end
                    entries = _realify_hermitian_linear_block(
                        complex_block,
                        LC;
                        atol=MP_R(phase_atol),
                    )
                    label = (
                        feature=:contiguous_rdm,
                        k=rdm_k,
                        decomposition=:qmbcertify,
                        block=block_idx,
                    )
                    logical_rows = _contiguous_rdm_state_labels(rdm_k, rows)
                    entries = _prepare_translation_linear_entries!(
                        builder,
                        entries,
                        registered_key_tokens,
                        rdm_axis_quotient,
                        rdm_quotient_form_cache,
                    )
                    constraint_idx += 1
                    add_psd_block!(
                        builder,
                        cone,
                        entries,
                        _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, nothing);
                        constraint_idx,
                    )
                    push!(logical_block_sizes, size(complex_block, 1))
                    push!(block_sizes, size(entries, 1))
                    push!(block_labels, label)
                    push!(block_logical_row_labels, logical_rows)
                    push!(block_transforms, nothing)
                end
                continue
            end

            if rdm_decomposition == :su2 && rdm_support == :extend
                constraint_idx = _add_translation_su2_extend_rdm_linear_blocks!(
                    builder,
                    n,
                    rdm_k,
                    MP_C,
                    real_moment_matrix,
                    MP_R(phase_atol),
                    registered_key_tokens,
                    prepare_axis_quotient,
                    logical_block_sizes,
                    block_sizes,
                    block_labels,
                    block_logical_row_labels,
                    block_transforms,
                    constraint_idx,
                )
                continue
            end

            complex_entries = _translation_contiguous_rdm_linear_block(n, rdm_k, K, MP_C)
            if rdm_support == :closed
                _check_linear_entries_moments_covered(
                    complex_entries,
                    constraint_moment_basis,
                    "Contiguous RDM k=$rdm_k block",
                )
            end
            if rdm_decomposition == :full
                entries = real_moment_matrix ?
                    _realify_hermitian_linear_block(complex_entries, LC; atol=MP_R(phase_atol)) :
                    complex_entries
                label = (feature=:contiguous_rdm, k=rdm_k)
                logical_rows = _contiguous_rdm_state_labels(rdm_k, collect(1:size(complex_entries, 1)))
                entries = _prepare_translation_linear_entries!(
                    builder,
                    entries,
                    registered_key_tokens,
                    prepare_axis_quotient,
                )
                constraint_idx += 1
                add_psd_block!(
                    builder,
                    cone,
                    entries,
                    _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, nothing);
                    constraint_idx,
                )
                push!(logical_block_sizes, size(complex_entries, 1))
                push!(block_sizes, size(entries, 1))
                push!(block_labels, label)
                push!(block_logical_row_labels, logical_rows)
                push!(block_transforms, nothing)
            elseif rdm_decomposition == :u1
                state_count = size(complex_entries, 1)
                magnetizations = [
                    _contiguous_rdm_magnetization(state) for state in 0:(state_count - 1)
                ]
                for col in 1:state_count, row in 1:state_count
                    magnetizations[row] == magnetizations[col] && continue
                    _add_translation_zero_linear_form!(
                        builder,
                        complex_entries[row, col];
                        phase_atol=MP_R(phase_atol),
                        label=(
                            feature=:contiguous_rdm_zero,
                            k=rdm_k,
                            decomposition=:u1,
                            reason=:magnetization_offblock,
                            row_state=row - 1,
                            row_magnetization=magnetizations[row],
                            col_state=col - 1,
                            col_magnetization=magnetizations[col],
                        ),
                    )
                end
                for (weight0, sector) in enumerate(_contiguous_rdm_u1_sectors(rdm_k))
                    complex_block = complex_entries[sector, sector]
                    entries = real_moment_matrix ?
                        _realify_hermitian_linear_block(
                            complex_block,
                            LC;
                            atol=MP_R(phase_atol),
                        ) :
                        complex_block
                    label = (
                        feature=:contiguous_rdm,
                        k=rdm_k,
                        decomposition=:u1,
                        magnetization=weight0 - 1,
                    )
                    logical_rows = _contiguous_rdm_state_labels(rdm_k, sector)
                    entries = _prepare_translation_linear_entries!(
                        builder,
                        entries,
                        registered_key_tokens,
                        prepare_axis_quotient,
                    )
                    constraint_idx += 1
                    add_psd_block!(
                        builder,
                        cone,
                        entries,
                        _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, nothing);
                        constraint_idx,
                    )
                    push!(logical_block_sizes, size(complex_block, 1))
                    push!(block_sizes, size(entries, 1))
                    push!(block_labels, label)
                    push!(block_logical_row_labels, logical_rows)
                    push!(block_transforms, nothing)
                end
            else
                schur_transform, states = _pauli_su2_schur_matrix(rdm_k)
                schur_rows = _pauli_sparse_transform_rows(
                    transpose(schur_transform);
                    atol=MP_R(phase_atol),
                )
                schur_entries = _transform_linear_block(
                    complex_entries,
                    schur_rows,
                    schur_rows;
                    atol=MP_R(phase_atol),
                )
                columns = _pauli_su2_state_columns(states)
                blocks = pauli_su2_rdm_blocks(rdm_k)

                for col_state in states, row_state in states
                    same_spin_m = row_state.spin2 == col_state.spin2 &&
                        row_state.m2 == col_state.m2
                    same_spin_m && continue
                    _add_translation_zero_linear_form!(
                        builder,
                        schur_entries[
                            columns[(row_state.spin2, row_state.m2, row_state.multiplicity)],
                            columns[(col_state.spin2, col_state.m2, col_state.multiplicity)],
                        ];
                        phase_atol=MP_R(phase_atol),
                        label=(
                            feature=:contiguous_rdm_zero,
                            k=rdm_k,
                            decomposition=:su2,
                            reason=:spin_magnetic_offblock,
                            row_spin2=row_state.spin2,
                            row_m2=row_state.m2,
                            row_multiplicity=row_state.multiplicity,
                            col_spin2=col_state.spin2,
                            col_m2=col_state.m2,
                            col_multiplicity=col_state.multiplicity,
                        ),
                    )
                end

                reference_rows_by_block = Vector{Vector{Int}}(undef, length(blocks))
                for (block_idx, block) in pairs(blocks)
                    reference_m2 = -block.spin2
                    reference_rows = [
                        columns[(block.spin2, reference_m2, mult)]
                        for mult in 1:block.multiplicity
                    ]
                    reference_rows_by_block[block_idx] = reference_rows
                    for m2 in (reference_m2 + 2):2:block.spin2
                        rows = [columns[(block.spin2, m2, mult)] for mult in 1:block.multiplicity]
                        for col in 1:block.multiplicity, row in 1:block.multiplicity
                            _add_translation_zero_linear_form!(
                                builder,
                                _subtract_linear_forms(
                                    schur_entries[rows[row], rows[col]],
                                    schur_entries[reference_rows[row], reference_rows[col]];
                                    atol=MP_R(phase_atol),
                                );
                                phase_atol=MP_R(phase_atol),
                                label=(
                                    feature=:contiguous_rdm_zero,
                                    k=rdm_k,
                                    decomposition=:su2,
                                    reason=:magnetic_copy,
                                    spin2=block.spin2,
                                    m2=m2,
                                    reference_m2=reference_m2,
                                    row_multiplicity=row,
                                    col_multiplicity=col,
                                ),
                            )
                        end
                    end
                end

                for (block, reference_rows) in zip(blocks, reference_rows_by_block)
                    complex_block = schur_entries[reference_rows, reference_rows]
                    entries = real_moment_matrix ?
                        _realify_hermitian_linear_block(complex_block, LC; atol=MP_R(phase_atol)) :
                        complex_block
                    label = (feature=:contiguous_rdm, k=rdm_k, decomposition=:su2, spin2=block.spin2)
                    logical_rows = _contiguous_rdm_su2_row_labels(rdm_k, block)
                    transform = TranslationSU2RDMTransform(rdm_k, block, schur_transform)
                    entries = _prepare_translation_linear_entries!(
                        builder,
                        entries,
                        registered_key_tokens,
                        prepare_axis_quotient,
                    )
                    constraint_idx += 1
                    add_psd_block!(
                        builder,
                        cone,
                        entries,
                        _translation_block_meta(M, cone, size(entries, 1), label, logical_rows, transform);
                        constraint_idx,
                    )
                    push!(logical_block_sizes, size(complex_block, 1))
                    push!(block_sizes, size(entries, 1))
                    push!(block_labels, label)
                    push!(block_logical_row_labels, logical_rows)
                    push!(block_transforms, transform)
                end
            end
        end
    end

    if axis_rotation_quotient && !use_pre_axis_quotient
        quotient_basis = sorted_unique!(collect(values(builder.key_to_monomial)))
        axis_rotation_quotient_stats = _translation_axis_rotation_quotient_map(
            K,
            M,
            ops,
            n,
            quotient_basis,
        )
        _apply_translation_axis_rotation_quotient!(builder, axis_rotation_quotient_stats)
    end

    block_assembly_done_ns = time_ns()
    construction_stage_times_ns[:block_assembly] = Int(block_assembly_done_ns - precompute_done_ns)
    construction_stage_times_ns[:constraint_append] = 0

    raw_linear = finalize!(builder)
    quotient_result = nothing
    linear = raw_linear
    if su2_moment_quotient
        quotient_start_ns = time_ns()
        quotient_result = _pauli_su2_quotient_linear_data(
            raw_linear,
            n;
            sign_symmetry,
            allow_registered_projection=isnothing(basis),
            atol=su2_quotient_atol,
            condition_limit=su2_quotient_condition_limit,
            stage_times_ns=construction_stage_times_ns,
        )
        linear = quotient_result.linear
        construction_stage_times_ns[:su2_moment_quotient] =
            Int(time_ns() - quotient_start_ns)
    end
    linearization_done_ns = time_ns()
    construction_stage_times_ns[:linearization] = Int(linearization_done_ns - block_assembly_done_ns)

    report = TranslationInvariantReport(
        n,
        d,
        local_basis_size,
        length(orbit_reps),
        axis_orbit_summary.axis_orbit_closed,
        axis_orbit_summary.axis_orbit_count,
        axis_orbit_summary.axis_orbit_size_histogram,
        axis_orbit_summary.reduction_ratio,
        sectors,
        sign_symmetry,
        logical_block_sizes,
        block_sizes,
        block_labels,
        _translation_block_coefficient_domains(block_transforms),
        _translation_block_exact_coefficient_domains(block_transforms),
        base_moment_count,
        length(linear.moments),
        length(linear.zero_constraints),
        Int(linearization_done_ns - construction_start_ns),
        construction_stage_times_ns,
        real_moment_matrix,
        product_cache.hits,
        product_cache.misses,
        length(product_cache.terms),
        _translation_product_cache_hit_rate(product_cache),
        _translation_zero_report_histogram(linear.zero_constraints),
        axis_rotation_quotient,
        axis_rotation_quotient_stats.raw_moment_key_count,
        axis_rotation_quotient_stats.moment_class_count,
        axis_rotation_quotient_stats.quotient_moment_key_count,
        axis_rotation_quotient_stats.forced_zero_moment_class_count,
        axis_rotation_quotient_stats.reduction_ratio,
    )
    if quotient_result !== nothing
        report = _translation_report_with_su2_moment_quotient(
            report,
            quotient_result,
        )
    end
    return linear, report
end

function pauli_translation_invariant_moment_relaxation(
    pop::PolyOpt{PauliAlgebra,T,P},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis::Union{Nothing,Vector{NormalMonomial{PauliAlgebra,T}}}=nothing,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    sign_symmetry::Bool=true,
    reflection_symmetry::Bool=false,
    axis_rotation_symmetry::Bool=false,
    check_invariance::Bool=true,
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition::Symbol=:full,
    contiguous_rdm_support::Symbol=:closed,
    u1_symmetry::Bool=false,
    su2_symmetry::Bool=false,
    axis_rotation_equalities::Bool=false,
    axis_rotation_quotient::Bool=false,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
    linear_state_opt_width=nothing,
    psd_state_opt_width=nothing,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    axis_rotation_quotient && throw(ArgumentError(
        "`axis_rotation_quotient=true` is currently supported only by the direct-linear Pauli translation constructor."
    ))
    construction_start_ns = time_ns()
    construction_stage_times_ns = Dict{Symbol,Int}()
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    eltype(σx) == NormalMonomial{PauliAlgebra,T} || throw(ArgumentError(
        "PolyOpt and Pauli chain operators must use the same Pauli integer type; got objective type $T and operator type $(eltype(σx))."
    ))

    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    one_mono = one(first(σx))

    _check_pauli_chain_support(pop.objective, n; context="objective")
    _check_translation_complex_reflection_sectors(
        n,
        momenta;
        reflection_symmetry,
        real_moment_matrix,
    )
    _check_translation_scalar_constraints_compatible(
        pop,
        n;
        check_invariance,
        sign_symmetry,
        reflection_symmetry,
        real_moment_matrix,
    )
    rdm_decomposition = _normalize_contiguous_rdm_decomposition(
        contiguous_rdm_decomposition,
        u1_symmetry,
        su2_symmetry,
    )
    rdm_support = _normalize_contiguous_rdm_support(contiguous_rdm_support)
    rdm_decomposition == :qmbcertify && throw(ArgumentError(
        "`contiguous_rdm_decomposition=:qmbcertify` is currently structural-target-only. " *
        "The constructor does not yet emit QMBCertify's shared RDM PSD variables."
    ))
    rdm_ks = _normalize_contiguous_rdm_k(contiguous_rdm_k)
    linear_state_width = _normalize_linear_state_opt_width(linear_state_opt_width)
    psd_state_width = _normalize_psd_state_opt_width(psd_state_opt_width)
    standalone_axis_reduction = axis_rotation_symmetry && !su2_symmetry
    (axis_rotation_equalities || standalone_axis_reduction) && check_invariance &&
        _check_pauli_su2_symmetry_compatible(pop, ops)
    u1_symmetry && _check_pauli_charge_neutral(pop)
    su2_symmetry && check_invariance && _check_pauli_su2_symmetry_compatible(pop, ops)
    setup_done_ns = time_ns()
    construction_stage_times_ns[:setup] = Int(setup_done_ns - construction_start_ns)

    local_basis_size = 0
    orbit_reps = NormalMonomial{PauliAlgebra,T}[]
    closure_basis = nothing
    if isnothing(basis) && d < n
        local_basis_size = _pauli_contiguous_chain_basis_size_hint(n, d; periodic=true)
        orbit_reps = _pauli_contiguous_chain_orbit_representatives(ops, d; periodic=true)
    else
        local_basis = isnothing(basis) ? pauli_contiguous_chain_basis(ops, d; periodic=true) : basis
        one_mono in local_basis || throw(ArgumentError("Translation-invariant Pauli basis must include the identity."))
        _check_pauli_chain_support(local_basis, n; context="basis")
        _check_translation_basis_closure(local_basis, n)
        local_basis_size = length(local_basis)
        orbit_reps = _translation_orbit_representatives(local_basis, n)
        closure_basis = isnothing(basis) ? nothing : local_basis
    end
    _check_pauli_translation_constraint_closure(
        closure_basis,
        d,
        n,
        pop.objective,
        rdm_ks,
        rdm_support,
        linear_state_width,
        psd_state_width;
        sign_symmetry,
    )
    basis_done_ns = time_ns()
    construction_stage_times_ns[:basis] = Int(basis_done_ns - setup_done_ns)
    axis_orbit_summary = _pauli_axis_orbit_summary(ops, orbit_reps, d; strict=false)
    axis_done_ns = time_ns()
    construction_stage_times_ns[:axis_diagnostics] = Int(axis_done_ns - basis_done_ns)
    nontrivial_reps = [m for m in orbit_reps if !isone(m)]
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    real_moment_matrix || 0 in sectors || throw(ArgumentError("Momentum sector 0 is required because it carries the normalized identity moment."))

    use_su2_base_reduction = su2_symmetry && (isempty(rdm_ks) || rdm_support == :closed)
    if use_su2_base_reduction
        use_wigner_su2_base = !reflection_symmetry && 0 in sectors
        su2_base_mp = if use_wigner_su2_base
            pauli_su2_translation_orbit_wigner_eckart_moment_problem(
                pop.objective,
                orbit_reps,
                n;
                assume_su2_invariant=!check_invariance,
                ops,
                momenta=sectors,
                phase_atol,
                singlet_channel_equalities,
                singlet_channel_atol,
            )
        else
            pauli_su2_translation_orbit_moment_problem(
                pop.objective,
                orbit_reps,
                n;
                assume_su2_invariant=!check_invariance,
                ops,
                real_moment_matrix,
                momenta,
                reflection_symmetry,
                phase_atol,
                singlet_channel_equalities,
                singlet_channel_atol,
            )
        end
        su2_mp = if real_moment_matrix && any(cone == :HPSD for (cone, _) in su2_base_mp.constraints)
            MP_R = _pauli_chain_real_coeff_type(C)
            MP_P = Polynomial{PauliAlgebra,T,MP_R}
            _realify_pauli_su2_translation_orbit_moment_problem(
                su2_base_mp,
                MP_P;
                phase_atol=MP_R(phase_atol),
            )
        else
            su2_base_mp
        end
        su2_mp = _append_translation_extra_constraints_to_su2_base_moment_problem(
            pop,
            su2_mp,
            orbit_reps,
            n;
            real_moment_matrix,
            sign_symmetry,
            linear_state_width,
            axis_rotation_equalities,
            rdm_ks,
            rdm_decomposition,
            rdm_support,
            psd_state_width,
            sectors,
            phase_atol,
            ops,
        )
        su2_done_ns = time_ns()
        construction_stage_times_ns[:su2_moment_problem] = Int(su2_done_ns - axis_done_ns)
        report = _translation_su2_base_report(
            su2_mp,
            n,
            d,
            local_basis_size,
            orbit_reps,
            axis_orbit_summary,
            sectors,
            sign_symmetry,
            Int(su2_done_ns - construction_start_ns),
            construction_stage_times_ns,
            real_moment_matrix,
        )
        return su2_mp, report
    end
    singlet_channel_equalities && throw(ArgumentError(
        "`singlet_channel_equalities=true` requires `su2_symmetry=true` and the " *
        "translation/SU(2) base reducer."
    ))

    MP_R = _pauli_chain_real_coeff_type(C)
    MP_C = Complex{MP_R}
    MP_P = Polynomial{PauliAlgebra,T,real_moment_matrix ? MP_R : MP_C}
    BLOCK_P = Polynomial{PauliAlgebra,T,MP_C}
    objective_mp = convert(MP_P, _translation_orbit_reduce_polynomial(pop.objective, n))

    translated = Dict{NormalMonomial{PauliAlgebra,T},Vector{NormalMonomial{PauliAlgebra,T}}}()
    for rep in nontrivial_reps
        translated[rep] = [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
    end
    translated[one_mono] = fill(one_mono, n)

    rep_cache = Dict{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}}()
    product_cache = _TranslationProductCache(
        Dict{
            Tuple{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}},
            Vector{Tuple{Int,MP_C,NormalMonomial{PauliAlgebra,T}}},
        }(),
        0,
        0,
    )
    constraints = Tuple{Symbol,Matrix{MP_P}}[]
    zero_origin_by_constraint = Dict{Int,Any}()
    moment_terms = NormalMonomial{PauliAlgebra,T}[]
    axis_equality_moment_terms = NormalMonomial{PauliAlgebra,T}[]
    logical_block_sizes = Int[]
    block_sizes = Int[]
    block_labels = Any[]
    block_logical_row_labels = Vector{Any}[]
    block_transforms = Any[]
    precompute_done_ns = time_ns()
    construction_stage_times_ns[:precompute] = Int(precompute_done_ns - axis_done_ns)

    for k in sectors
        sector_basis = k == 0 ? orbit_reps : nontrivial_reps
        blocks = sign_symmetry && !standalone_axis_reduction ?
            _pauli_signature_blocks(sector_basis) :
            [(:all, sector_basis)]
        for (signature, block_basis) in blocks
            isempty(block_basis) && continue
            complex_mat = _translation_momentum_block(
                block_basis,
                k,
                n,
                translated,
                rep_cache,
                BLOCK_P,
                product_cache,
            )
            if standalone_axis_reduction && reflection_symmetry
                for entry in complex_mat
                    append!(axis_equality_moment_terms, monomials(entry))
                end
                if _translation_reflection_fixed_momentum(k, n)
                    for reflected in _translation_reflection_adapted_blocks(block_basis, k, n; atol=Float64(phase_atol))
                        for axis_adapted in _translation_reflection_axis_isotypic_adapted_blocks(
                            ops,
                            block_basis,
                            d,
                            reflected.reflection,
                            reflected.row_basis,
                            ;
                            realified=false,
                            atol=Float64(phase_atol),
                        )
                            adapted_complex_mat = _hermitianize_pauli_polynomial_block(
                                _transform_polynomial_block(
                                    complex_mat,
                                    axis_adapted.row_basis,
                                    axis_adapted.row_basis;
                                    atol=Float64(phase_atol),
                                ),
                                BLOCK_P,
                            )
                            cone, final_mat = real_moment_matrix ?
                                (:PSD, _realify_hermitian_block(adapted_complex_mat, MP_P; atol=MP_R(phase_atol))) :
                                (:HPSD, map(p -> convert(MP_P, p), adapted_complex_mat))
                            final_label = (
                                momentum=k,
                                signature=signature,
                                reflection=axis_adapted.reflection,
                                axis_irrep=axis_adapted.axis_irrep,
                                axis_irrep_dimension=axis_adapted.irrep_dimension,
                                axis_irrep_multiplicity=axis_adapted.irrep_multiplicity,
                                axis_group_order=axis_adapted.axis_group_order,
                            )
                            transform = TranslationAxisIrrepTransform(
                                TranslationReflectionTransform(
                                    n,
                                    k,
                                    axis_adapted.reflection,
                                    real_moment_matrix,
                                    false,
                                ),
                                axis_adapted.axis_irrep,
                                axis_adapted.irrep_dimension,
                                axis_adapted.irrep_multiplicity,
                                axis_adapted.axis_group_order,
                            )
                            push!(constraints, (cone, final_mat))
                            push!(logical_block_sizes, size(adapted_complex_mat, 1))
                            push!(block_sizes, size(final_mat, 1))
                            push!(block_labels, final_label)
                            push!(block_logical_row_labels, axis_adapted.row_labels)
                            push!(block_transforms, transform)
                            for entry in final_mat
                                append!(moment_terms, monomials(entry))
                            end
                        end
                    end
                else
                    real_moment_matrix || throw(ArgumentError(
                        "Symbolic axis/reflection reduction supports non-fixed reflected sectors only for real_moment_matrix=true; got k=$k, n=$n."
                    ))
                    real_mat = _realify_hermitian_block(complex_mat, MP_P; atol=MP_R(phase_atol))
                    for reflected in _translation_real_reflection_adapted_blocks(block_basis, k, n; atol=Float64(phase_atol))
                        for axis_adapted in _translation_reflection_axis_isotypic_adapted_blocks(
                            ops,
                            block_basis,
                            d,
                            reflected.reflection,
                            reflected.row_basis;
                            realified=true,
                            atol=Float64(phase_atol),
                        )
                            final_mat = _symmetrize_real_polynomial_block(
                                _transform_polynomial_block(
                                    real_mat,
                                    axis_adapted.row_basis,
                                    axis_adapted.row_basis;
                                    atol=Float64(phase_atol),
                                ),
                                MP_P,
                            )
                            final_label = (
                                momentum=k,
                                signature=signature,
                                reflection=axis_adapted.reflection,
                                axis_irrep=axis_adapted.axis_irrep,
                                axis_irrep_dimension=axis_adapted.irrep_dimension,
                                axis_irrep_multiplicity=axis_adapted.irrep_multiplicity,
                                axis_group_order=axis_adapted.axis_group_order,
                            )
                            transform = TranslationAxisIrrepTransform(
                                TranslationReflectionTransform(
                                    n,
                                    k,
                                    axis_adapted.reflection,
                                    real_moment_matrix,
                                    true,
                                ),
                                axis_adapted.axis_irrep,
                                axis_adapted.irrep_dimension,
                                axis_adapted.irrep_multiplicity,
                                axis_adapted.axis_group_order,
                            )
                            push!(constraints, (:PSD, final_mat))
                            push!(logical_block_sizes, size(final_mat, 1))
                            push!(block_sizes, size(final_mat, 1))
                            push!(block_labels, final_label)
                            push!(block_logical_row_labels, axis_adapted.row_labels)
                            push!(block_transforms, transform)
                            for entry in final_mat
                                append!(moment_terms, monomials(entry))
                            end
                        end
                    end
                end
                continue
            end
            if standalone_axis_reduction
                for entry in complex_mat
                    append!(axis_equality_moment_terms, monomials(entry))
                end
                for adapted in _translation_axis_isotypic_adapted_blocks(
                    ops,
                    block_basis,
                    d;
                    atol=Float64(phase_atol),
                )
                    adapted_complex_mat = _hermitianize_pauli_polynomial_block(
                        _transform_polynomial_block(
                            complex_mat,
                            adapted.row_basis,
                            adapted.row_basis;
                            atol=Float64(phase_atol),
                        ),
                        BLOCK_P,
                    )
                    cone, final_mat = real_moment_matrix ?
                        (:PSD, _realify_hermitian_block(adapted_complex_mat, MP_P; atol=MP_R(phase_atol))) :
                        (:HPSD, map(p -> convert(MP_P, p), adapted_complex_mat))
                    final_label = (
                        momentum=k,
                        signature=signature,
                        axis_irrep=adapted.axis_irrep,
                        axis_irrep_dimension=adapted.irrep_dimension,
                        axis_irrep_multiplicity=adapted.irrep_multiplicity,
                        axis_group_order=adapted.axis_group_order,
                    )
                    transform = TranslationAxisIrrepTransform(
                        TranslationDFTTransform(n, k, real_moment_matrix),
                        adapted.axis_irrep,
                        adapted.irrep_dimension,
                        adapted.irrep_multiplicity,
                        adapted.axis_group_order,
                    )
                    push!(constraints, (cone, final_mat))
                    push!(logical_block_sizes, size(adapted_complex_mat, 1))
                    push!(block_sizes, size(final_mat, 1))
                    push!(block_labels, final_label)
                    push!(block_logical_row_labels, adapted.row_labels)
                    push!(block_transforms, transform)
                    for entry in final_mat
                        append!(moment_terms, monomials(entry))
                    end
                end
                continue
            end
            adapted_blocks = if reflection_symmetry && _translation_reflection_fixed_momentum(k, n)
                _translation_reflection_adapted_blocks(block_basis, k, n; atol=Float64(phase_atol))
            else
                [(reflection=nothing, row_basis=nothing, row_labels=Any[mono for mono in block_basis])]
            end
            for adapted in adapted_blocks
                adapted_complex_mat = if adapted.row_basis === nothing
                    complex_mat
                else
                    _hermitianize_pauli_polynomial_block(
                        _transform_polynomial_block(
                            complex_mat,
                            adapted.row_basis,
                            adapted.row_basis;
                            atol=Float64(phase_atol),
                        ),
                        BLOCK_P,
                    )
                end
                cone, mat = real_moment_matrix ?
                    (:PSD, _realify_hermitian_block(adapted_complex_mat, MP_P; atol=MP_R(phase_atol))) :
                    (:HPSD, map(p -> convert(MP_P, p), adapted_complex_mat))
                label = adapted.reflection === nothing ?
                    (momentum=k, signature=signature) :
                    (momentum=k, signature=signature, reflection=adapted.reflection)
                real_reflection_blocks = if reflection_symmetry &&
                        real_moment_matrix &&
                        adapted.reflection === nothing &&
                        !_translation_reflection_fixed_momentum(k, n) &&
                        size(mat, 1) == 2 * size(adapted_complex_mat, 1)
                    _translation_real_reflection_adapted_blocks(block_basis, k, n; atol=Float64(phase_atol))
                else
                    [(reflection=nothing, row_basis=nothing, row_labels=adapted.row_labels)]
                end

                for real_adapted in real_reflection_blocks
                    final_mat = real_adapted.row_basis === nothing ?
                        mat :
                        _symmetrize_real_polynomial_block(
                            _transform_polynomial_block(
                                mat,
                                real_adapted.row_basis,
                                real_adapted.row_basis;
                                atol=Float64(phase_atol),
                            ),
                            MP_P,
                        )
                    final_label = real_adapted.reflection === nothing ?
                        label :
                        merge(label, (reflection=real_adapted.reflection,))
                    logical_size = real_adapted.row_basis === nothing ?
                        size(adapted_complex_mat, 1) :
                        size(final_mat, 1)
                    push!(constraints, (cone, final_mat))
                    push!(logical_block_sizes, logical_size)
                    push!(block_sizes, size(final_mat, 1))
                    push!(block_labels, final_label)
                    push!(block_logical_row_labels, real_adapted.row_labels)
                    push!(
                        block_transforms,
                        _translation_transform_descriptor(n, k, real_moment_matrix, final_label),
                    )
                    for entry in final_mat
                        append!(moment_terms, monomials(entry))
                    end
                end
            end
        end
    end
    block_assembly_done_ns = time_ns()
    construction_stage_times_ns[:block_assembly] = Int(block_assembly_done_ns - precompute_done_ns)

    moment_basis = sorted_unique!(moment_terms)
    constraint_moment_basis = standalone_axis_reduction ?
        sorted_unique!(vcat(copy(moment_basis), axis_equality_moment_terms)) :
        moment_basis
    _check_objective_moments_covered(objective_mp, moment_basis)
    axis_rotation_quotient_stats = axis_rotation_equalities || standalone_axis_reduction ?
        _translation_axis_rotation_quotient_stats(ops, n, constraint_moment_basis) :
        _EMPTY_AXIS_ROTATION_QUOTIENT_STATS
    if axis_rotation_equalities || standalone_axis_reduction
        _append_translation_axis_rotation_constraints!(
            constraints,
            ops,
            n,
            constraint_moment_basis,
            MP_P,
            zero_origin_by_constraint,
        )
    end
    _append_translation_scalar_constraints!(
        constraints,
        pop,
        n,
        constraint_moment_basis,
        MP_P,
        logical_block_sizes,
        block_sizes,
        block_labels,
        block_logical_row_labels,
        block_transforms,
        zero_origin_by_constraint,
    )
    _append_translation_moment_eq_constraints!(
        constraints,
        pop,
        n,
        orbit_reps,
        sign_symmetry,
        constraint_moment_basis,
        MP_P,
        zero_origin_by_constraint,
    )
    _append_translation_linear_state_opt_constraints!(
        constraints,
        pop.objective,
        n,
        linear_state_width,
        sign_symmetry,
        constraint_moment_basis,
        MP_P,
        zero_origin_by_constraint,
    )
    _append_translation_psd_state_opt_constraints!(
        constraints,
        pop.objective,
        n,
        psd_state_width,
        sectors,
        sign_symmetry,
        constraint_moment_basis,
        BLOCK_P,
        MP_P,
        real_moment_matrix,
        MP_R(phase_atol),
        logical_block_sizes,
        block_sizes,
        block_labels,
        block_logical_row_labels,
        block_transforms,
    )
    _append_translation_contiguous_rdm_constraints!(
        constraints,
        n,
        rdm_ks,
        rdm_decomposition,
        rdm_support,
        constraint_moment_basis,
        BLOCK_P,
        MP_P,
        real_moment_matrix,
        MP_R(phase_atol),
        zero_origin_by_constraint,
        logical_block_sizes,
        block_sizes,
        block_labels,
        block_logical_row_labels,
        block_transforms,
    )
    constraint_append_done_ns = time_ns()
    construction_stage_times_ns[:constraint_append] = Int(
        constraint_append_done_ns - block_assembly_done_ns
    )
    total_basis = _collect_translation_moment_basis(objective_mp, constraints)
    n_unique = length(moment_basis)
    M = NormalMonomial{PauliAlgebra,T}
    block_meta_by_constraint = _translation_block_meta_by_constraint(
        M,
        constraints,
        block_labels,
        block_logical_row_labels,
        block_transforms,
    )

    mp = MomentProblem{PauliAlgebra,T,NormalMonomial{PauliAlgebra,T},MP_P}(
        objective_mp,
        constraints,
        total_basis,
        n_unique;
        block_meta_by_constraint=block_meta_by_constraint,
        zero_origin_by_constraint=zero_origin_by_constraint,
        real_moments=real_moment_matrix,
    )
    linearization_done_ns = time_ns()
    construction_stage_times_ns[:linearization] = Int(
        linearization_done_ns - constraint_append_done_ns
    )

    report = TranslationInvariantReport(
        n,
        d,
        local_basis_size,
        length(orbit_reps),
        axis_orbit_summary.axis_orbit_closed,
        axis_orbit_summary.axis_orbit_count,
        axis_orbit_summary.axis_orbit_size_histogram,
        axis_orbit_summary.reduction_ratio,
        sectors,
        sign_symmetry,
        logical_block_sizes,
        block_sizes,
        block_labels,
        _translation_block_coefficient_domains(block_transforms),
        _translation_block_exact_coefficient_domains(block_transforms),
        n_unique,
        length(mp.linear.moments),
        length(mp.linear.zero_constraints),
        Int(linearization_done_ns - construction_start_ns),
        construction_stage_times_ns,
        real_moment_matrix,
        product_cache.hits,
        product_cache.misses,
        length(product_cache.terms),
        _translation_product_cache_hit_rate(product_cache),
        _translation_zero_report_histogram(mp.linear.zero_constraints),
        false,
        axis_rotation_quotient_stats.raw_moment_key_count,
        axis_rotation_quotient_stats.moment_class_count,
        axis_rotation_quotient_stats.quotient_moment_key_count,
        axis_rotation_quotient_stats.forced_zero_moment_class_count,
        axis_rotation_quotient_stats.reduction_ratio,
    )
    return mp, report
end

"""
    pauli_translation_invariant_nctssos(pop, (σx, σy, σz), order, optimizer; kwargs...)

Build and solve the translation-invariant Pauli-chain relaxation from
[`pauli_translation_invariant_moment_relaxation`](@ref).

Set `qmbcertify_base_construct=true` to solve the QMBCertify-style
source-family base construction directly as linear SDP data.  That route also
supports QMBCertify-mode `contiguous_rdm_k` blocks, `linear_state_opt_width`
rows, and `psd_state_opt_width` blocks.
"""
function pauli_translation_invariant_nctssos(
    pop::PolyOpt{PauliAlgebra,T,P},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer,
    optimizer;
    dualize::Bool=false,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
    sos_hermitian_representation::Symbol=:real_lift,
    direct_linear::Bool=false,
    base_su2_extend_rdm::Bool=false,
    su2_moment_quotient::Bool=false,
    su2_moment_quotient_atol::Real=1e-11,
    su2_moment_quotient_condition_limit::Real=1e10,
    qmbcertify_base_construct::Bool=false,
    qmbcertify_base_extra=nothing,
    qmbcertify_base_three_type::Tuple{<:Integer,<:Integer}=(1, 1),
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition=nothing,
    contiguous_rdm_support=nothing,
    linear_state_opt_width=nothing,
    linear_state_opt_mode=nothing,
    psd_state_opt_width=nothing,
    kwargs...,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    !su2_moment_quotient &&
        (su2_moment_quotient_atol != 1e-11 ||
            su2_moment_quotient_condition_limit != 1e10) && throw(ArgumentError(
                "SU(2) moment quotient numerical options require " *
                "`su2_moment_quotient=true`."
            ))
    moment_problem, report = if qmbcertify_base_construct
        direct_linear && throw(ArgumentError(
            "`qmbcertify_base_construct=true` already selects a direct linear-data constructor; do not also set `direct_linear=true`."
        ))
        base_su2_extend_rdm && throw(ArgumentError(
            "`qmbcertify_base_construct=true` cannot be combined with `base_su2_extend_rdm=true`."
        ))
        su2_moment_quotient && throw(ArgumentError(
            "`su2_moment_quotient=true` is not supported by the QMBCertify-base route."
        ))
        (su2_moment_quotient_atol != 1e-11 ||
            su2_moment_quotient_condition_limit != 1e10) && throw(ArgumentError(
                "SU(2) moment quotient numerical options require " *
                "`direct_linear=true` and `su2_moment_quotient=true`."
            ))
        if !isempty(kwargs)
            unsupported = join(sort!(String[string(k) for k in keys(kwargs)]), ", ")
            throw(ArgumentError(
                "`qmbcertify_base_construct=true` currently supports the source-base route plus `qmbcertify_base_extra`, `qmbcertify_base_three_type`, `real_moment_matrix`, `phase_atol`, `contiguous_rdm_k`, `contiguous_rdm_decomposition`, `contiguous_rdm_support`, `linear_state_opt_width`, `linear_state_opt_mode`, and `psd_state_opt_width`; unsupported keyword(s): $unsupported."
            ))
        end
        qmb_rdm_decomposition = something(contiguous_rdm_decomposition, :qmbcertify)
        qmb_rdm_support = something(contiguous_rdm_support, :extend)
        extra = if isnothing(qmbcertify_base_extra)
            0
        elseif qmbcertify_base_extra isa Integer
            Int(qmbcertify_base_extra)
        else
            throw(ArgumentError(
                "`qmbcertify_base_extra` must be an integer or `nothing`; got $(typeof(qmbcertify_base_extra))."
            ))
        end
        _pauli_qmbcertify_chain_base_linear_relaxation(
            pop,
            ops,
            order;
            extra,
            three_type=(
                Int(qmbcertify_base_three_type[1]),
                Int(qmbcertify_base_three_type[2]),
            ),
            real_moment_matrix,
            phase_atol,
            contiguous_rdm_k,
            contiguous_rdm_decomposition=qmb_rdm_decomposition,
            contiguous_rdm_support=qmb_rdm_support,
            linear_state_opt_width,
            linear_state_opt_mode=something(linear_state_opt_mode, :qmbcertify),
            psd_state_opt_width,
        )
    elseif direct_linear
        direct_rdm_decomposition = something(contiguous_rdm_decomposition, :full)
        direct_rdm_support = something(contiguous_rdm_support, :closed)
        direct_linear_state_mode =
            _normalize_linear_state_opt_mode(something(linear_state_opt_mode, :contiguous))
        _pauli_translation_base_linear_relaxation(
            pop,
            ops,
            order;
            base_su2_extend_rdm,
            su2_moment_quotient,
            su2_moment_quotient_atol,
            su2_moment_quotient_condition_limit,
            real_moment_matrix,
            phase_atol,
            contiguous_rdm_k,
            contiguous_rdm_decomposition=direct_rdm_decomposition,
            contiguous_rdm_support=direct_rdm_support,
            linear_state_opt_width,
            linear_state_opt_mode=direct_linear_state_mode,
            psd_state_opt_width,
            kwargs...,
        )
    else
        base_su2_extend_rdm && throw(ArgumentError(
            "`base_su2_extend_rdm=true` requires `direct_linear=true`."
        ))
        su2_moment_quotient && throw(ArgumentError(
            "`su2_moment_quotient=true` requires `direct_linear=true`."
        ))
        (su2_moment_quotient_atol != 1e-11 ||
            su2_moment_quotient_condition_limit != 1e10) && throw(ArgumentError(
                "SU(2) moment quotient numerical options require " *
                "`direct_linear=true` and `su2_moment_quotient=true`."
            ))
        symbolic_rdm_decomposition = something(contiguous_rdm_decomposition, :full)
        symbolic_rdm_support = something(contiguous_rdm_support, :closed)
        symbolic_linear_state_mode =
            _normalize_linear_state_opt_mode(something(linear_state_opt_mode, :contiguous))
        symbolic_linear_state_mode == :contiguous || throw(ArgumentError(
            "`linear_state_opt_mode=:qmbcertify` requires `direct_linear=true` or `qmbcertify_base_construct=true`."
        ))
        pauli_translation_invariant_moment_relaxation(
            pop,
            ops,
            order;
            real_moment_matrix,
            phase_atol,
            contiguous_rdm_k,
            contiguous_rdm_decomposition=symbolic_rdm_decomposition,
            contiguous_rdm_support=symbolic_rdm_support,
            linear_state_opt_width,
            psd_state_opt_width,
            kwargs...,
        )
    end
    support = translation_solve_support(report)
    if !support.supported
        throw(ArgumentError(support.reason))
    end
    linear = (direct_linear || qmbcertify_base_construct) ? moment_problem : moment_problem.linear
    solve_start_ns = time_ns()
    result = solve_sdp(
        linear,
        optimizer;
        dualize,
        formulation,
        representation,
        orphan_policy,
        sos_hermitian_representation,
    )
    solve_time_ns = Int(time_ns() - solve_start_ns)
    return TranslationInvariantResult(result.objective, result.model, moment_problem, report, solve_time_ns)
end

# -----------------------------------------------------------------------------
# Internals
# -----------------------------------------------------------------------------

function _validate_pauli_chain_ops(ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector})
    σx, σy, σz = ops
    n = length(σx)
    n > 0 || throw(ArgumentError("Pauli chain needs at least one site."))
    length(σy) == n && length(σz) == n || throw(ArgumentError(
        "Pauli chain operator vectors must have equal length; got $(length(σx)), $(length(σy)), $(length(σz))."
    ))
    eltype(σx) == eltype(σy) == eltype(σz) || throw(ArgumentError("Pauli chain operator vectors must have the same monomial type."))
    M = eltype(σx)
    M <: NormalMonomial{PauliAlgebra} || throw(ArgumentError("Pauli chain helpers require `PauliAlgebra` monomials."))
    all(m -> degree(m) == 1, Iterators.flatten(ops)) || throw(ArgumentError("Pauli chain operator vectors must contain single Pauli-letter monomials."))

    seen = Set{eltype(first(σx).word)}()
    for (label, vec, ptype) in (("σx", σx, _PAULI_X_TYPE), ("σy", σy, _PAULI_Y_TYPE), ("σz", σz, _PAULI_Z_TYPE))
        for i in 1:n
            idx = only(vec[i].word)
            _pauli_site(idx) == i || throw(ArgumentError(
                "Pauli chain operator $label[$i] must use canonical encoded site $i, got site $(_pauli_site(idx))."
            ))
            _pauli_type(idx) == ptype || throw(ArgumentError(
                "Pauli chain operator $label[$i] has Pauli type $(_pauli_type(idx)); expected $ptype."
            ))
            idx in seen && throw(ArgumentError("Pauli chain operator vectors contain duplicate Pauli index $idx."))
            push!(seen, idx)
        end
    end
    return σx, σy, σz, n
end

_pauli_chain_real_coeff_type(::Type{C}) where {C<:Number} = typeof(float(real(one(C))))
_pauli_chain_complex_coeff_type(::Type{C}) where {C<:Number} = Complex{_pauli_chain_real_coeff_type(C)}

function _pauli_chain_momentum_sectors(
    n::Integer,
    momenta::Union{Nothing,AbstractVector{<:Integer}};
    real_moment_matrix::Bool,
)
    if isnothing(momenta)
        return real_moment_matrix ? collect(0:fld(Int(n), 2)) : collect(0:(Int(n) - 1))
    end

    sectors = if real_moment_matrix
        [min(mod(Int(k), Int(n)), mod(-Int(k), Int(n))) for k in momenta]
    else
        [mod(Int(k), Int(n)) for k in momenta]
    end
    isempty(sectors) && throw(ArgumentError("At least one momentum sector is required."))
    unique!(sectors)
    sort!(sectors)
    return sectors
end

function _check_translation_complex_reflection_sectors(
    n_sites::Integer,
    momenta;
    reflection_symmetry::Bool,
    real_moment_matrix::Bool,
)
    reflection_symmetry || return nothing
    real_moment_matrix && return nothing

    n = Int(n_sites)
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    nonfixed = Int[k for k in sectors if !_translation_reflection_fixed_momentum(k, n)]
    isempty(nonfixed) && return nothing

    throw(ArgumentError(
        "`real_moment_matrix=false` with reflected Pauli-chain symmetry is supported " *
        "only for reflection-fixed momentum sectors; got non-fixed sector(s) $nonfixed. " *
        "Pass `reflection_symmetry=false`, keep `real_moment_matrix=true`, or restrict " *
        "`momenta` to reflection-fixed sectors."
    ))
end

function _check_real_pauli_chain_polynomial(
    poly::Polynomial{PauliAlgebra,T,C};
    context::AbstractString,
) where {T<:Unsigned,C<:Number}
    for (coef, mono) in poly.terms
        iszero(imag(coef)) || throw(ArgumentError(
            "real_moment_matrix=true requires real $context coefficients; term $mono has coefficient $coef."
        ))
    end
    return nothing
end

_check_real_pauli_chain_objective(poly::Polynomial{PauliAlgebra}) =
    _check_real_pauli_chain_polynomial(poly; context="objective")

@inline _pauli_index_from_site_type(::Type{T}, site::Integer, pauli_type::Integer) where {T<:Integer} =
    convert(T, 3 * (site - 1) + pauli_type + 1)

function _translate_pauli_monomial(
    mono::NormalMonomial{PauliAlgebra,T},
    shift::Integer,
    n_sites::Integer,
) where {T<:Unsigned}
    isempty(mono.word) && return mono
    n = Int(n_sites)
    s = mod(Int(shift), n)
    word = Vector{T}(undef, length(mono.word))
    for (i, idx) in pairs(mono.word)
        site = _pauli_site(idx)
        ptype = _pauli_type(idx)
        word[i] = _pauli_index_from_site_type(T, mod1(site + s, n), ptype)
    end
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError("Internal error: Pauli translation produced phase $phase."))
    return NormalMonomial{PauliAlgebra,T}(copy(simplified))
end

function _translation_orbit_representative(
    mono::NormalMonomial{PauliAlgebra,T},
    n_sites::Integer,
) where {T<:Unsigned}
    rep = mono
    for shift in 1:(Int(n_sites) - 1)
        image = _translate_pauli_monomial(mono, shift, n_sites)
        image < rep && (rep = image)
    end
    return rep
end

function _translation_orbit_length(mono::NormalMonomial{PauliAlgebra}, n_sites::Integer)
    isone(mono) && return 1
    for shift in 1:(Int(n_sites) - 1)
        _translate_pauli_monomial(mono, shift, n_sites) == mono && return shift
    end
    return Int(n_sites)
end

function _translation_orbit_representatives(
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
) where {T<:Unsigned}
    reps = Dict{NormalMonomial{PauliAlgebra,T},Nothing}()
    for mono in basis
        rep = _translation_orbit_representative(mono, n_sites)
        len = _translation_orbit_length(rep, n_sites)
        (isone(rep) || len == Int(n_sites)) || throw(ArgumentError(
            "Translation-invariant Pauli path currently supports only identity and full-length translation orbits; orbit representative $rep has length $len."
        ))
        reps[rep] = nothing
    end
    return sort!(collect(keys(reps)))
end

function _check_pauli_chain_support(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer;
    context::AbstractString,
) where {T<:Unsigned,C<:Number}
    n = Int(n_sites)
    for (_, mono) in poly.terms
        _check_pauli_chain_support(mono, n; context)
    end
    return nothing
end

function _check_pauli_chain_support(
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer;
    context::AbstractString,
) where {T<:Unsigned}
    n = Int(n_sites)
    for mono in basis
        _check_pauli_chain_support(mono, n; context)
    end
    return nothing
end

function _check_pauli_chain_support(
    mono::NormalMonomial{PauliAlgebra},
    n_sites::Integer;
    context::AbstractString,
)
    n = Int(n_sites)
    for idx in mono.word
        site = _pauli_site(idx)
        1 <= site <= n || throw(ArgumentError(
            "Translation-invariant Pauli chain $context contains site $site outside the declared chain 1:$n."
        ))
    end
    return nothing
end

function _check_translation_basis_closure(
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
) where {T<:Unsigned}
    basis_set = Set(basis)
    for mono in basis
        image = _translate_pauli_monomial(mono, 1, n_sites)
        image in basis_set || throw(ArgumentError(
            "Translation-invariant Pauli basis is not closed under one-site translation; $mono maps to $image."
        ))
    end
    return nothing
end

function _check_translation_invariance(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer;
    context::AbstractString="objective",
) where {T<:Unsigned,C<:Number}
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:Int(n_sites), ptype in 0:2
        src = NormalMonomial{PauliAlgebra,T}(T[_pauli_index_from_site_type(T, site, ptype)])
        dst = NormalMonomial{PauliAlgebra,T}(T[_pauli_index_from_site_type(T, mod1(site + 1, Int(n_sites)), ptype)])
        images[src] = (1, dst)
    end
    translation = CliffordSymmetry(images; nqubits=Int(n_sites))
    _act_polynomial(translation, poly) == poly || throw(ArgumentError(
        "Translation-invariant Pauli relaxation requires a one-site translation-invariant $context."
    ))
    return nothing
end

function _translation_orbit_reduce_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer,
) where {T<:Unsigned,C<:Number}
    terms = Tuple{C,NormalMonomial{PauliAlgebra,T}}[]
    sizehint!(terms, length(poly.terms))
    cache = Dict{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}}()
    for (coef, mono) in poly.terms
        rep = get!(cache, mono) do
            _translation_orbit_representative(mono, n_sites)
        end
        push!(terms, (coef, rep))
    end
    return Polynomial(terms)
end

function _translation_reflection_reduce_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer,
) where {T<:Unsigned,C<:Number}
    terms = Tuple{C,NormalMonomial{PauliAlgebra,T}}[]
    sizehint!(terms, length(poly.terms))
    cache = Dict{NormalMonomial{PauliAlgebra,T},NormalMonomial{PauliAlgebra,T}}()
    for (coef, mono) in poly.terms
        rep = get!(cache, mono) do
            reflected = _reflect_pauli_monomial(mono, n_sites)
            _translation_orbit_representative(reflected, n_sites)
        end
        push!(terms, (coef, rep))
    end
    return Polynomial(terms)
end

function _check_translation_reflection_invariance(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer;
    context::AbstractString="objective",
) where {T<:Unsigned,C<:Number}
    reduced = _translation_orbit_reduce_polynomial(poly, n_sites)
    reflected = _translation_reflection_reduce_polynomial(poly, n_sites)
    reflected == reduced || throw(ArgumentError(
        "reflection_symmetry=true requires a reflection-invariant $context."
    ))
    return nothing
end

function _check_pauli_sign_invariance(
    poly::Polynomial{PauliAlgebra,T,C};
    context::AbstractString="objective",
) where {T<:Unsigned,C<:Number}
    for (coef, mono) in poly.terms
        iszero(coef) && continue
        _pauli_sign_signature(mono) == 0x00 || throw(ArgumentError(
            "sign_symmetry=true requires the $context to be invariant under global Heisenberg sign flips; term $mono has nontrivial signature. Pass sign_symmetry=false for non-invariant objectives."
        ))
    end
    return nothing
end

function _check_translation_compatible_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer;
    context::AbstractString,
    check_invariance::Bool,
    sign_symmetry::Bool,
    reflection_symmetry::Bool,
    real_moment_matrix::Bool,
) where {T<:Unsigned,C<:Number}
    _check_pauli_chain_support(poly, n_sites; context)
    check_invariance && _check_translation_invariance(poly, n_sites; context)
    check_invariance && reflection_symmetry && _check_translation_reflection_invariance(
        poly,
        n_sites;
        context,
    )
    sign_symmetry && _check_pauli_sign_invariance(poly; context)
    real_moment_matrix && _check_real_pauli_chain_polynomial(poly; context)
    return nothing
end

function _check_translation_scalar_constraints_compatible(
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer;
    check_invariance::Bool,
    sign_symmetry::Bool,
    reflection_symmetry::Bool,
    real_moment_matrix::Bool,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    _check_translation_compatible_polynomial(
        pop.objective,
        n_sites;
        context="objective",
        check_invariance,
        sign_symmetry,
        reflection_symmetry,
        real_moment_matrix,
    )
    for (idx, poly) in pairs(pop.eq_constraints)
        _check_translation_compatible_polynomial(
            poly,
            n_sites;
            context="equality constraint $idx",
            check_invariance,
            sign_symmetry,
            reflection_symmetry,
            real_moment_matrix,
        )
    end
    for (idx, poly) in pairs(pop.ineq_constraints)
        _check_translation_compatible_polynomial(
            poly,
            n_sites;
            context="inequality constraint $idx",
            check_invariance,
            sign_symmetry,
            reflection_symmetry,
            real_moment_matrix,
        )
    end
    for (idx, poly) in pairs(pop.moment_eq_constraints)
        _check_translation_compatible_polynomial(
            poly,
            n_sites;
            context="moment equality constraint $idx",
            check_invariance,
            sign_symmetry,
            reflection_symmetry,
            real_moment_matrix,
        )
    end
    return nothing
end

function _check_pauli_axis_rotation_invariance(
    poly::Polynomial{PauliAlgebra,T,C},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector};
    context::AbstractString,
) where {T<:Unsigned,C<:Number}
    for rotation in pauli_global_axis_rotation_generators(ops)
        _act_polynomial(rotation, poly) == poly || throw(ArgumentError(
            "su2_symmetry=true requires a global Pauli-axis-rotation-invariant $context."
        ))
    end
    return nothing
end

function _check_pauli_su2_symmetry_compatible(
    pop::PolyOpt{PauliAlgebra,T,P},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    _check_pauli_axis_rotation_invariance(pop.objective, ops; context="objective")
    for (idx, poly) in pairs(pop.eq_constraints)
        _check_pauli_axis_rotation_invariance(poly, ops; context="equality constraint $idx")
    end
    for (idx, poly) in pairs(pop.ineq_constraints)
        _check_pauli_axis_rotation_invariance(poly, ops; context="inequality constraint $idx")
    end
    for (idx, poly) in pairs(pop.moment_eq_constraints)
        _check_pauli_axis_rotation_invariance(poly, ops; context="moment equality constraint $idx")
    end
    return nothing
end

function _check_polynomial_moments_covered(
    poly::Polynomial{PauliAlgebra,T,C},
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    context::AbstractString,
) where {T<:Integer,C<:Number}
    moment_set = Set(moment_basis)
    missing = NormalMonomial{PauliAlgebra,T}[]
    for mono in monomials(poly)
        mono in moment_set || push!(missing, mono)
    end
    isempty(missing) && return nothing
    shown = join((sprint(show, mono) for mono in Iterators.take(missing, 5)), ", ")
    length(missing) > 5 && (shown *= ", ...")
    throw(ArgumentError(
        "$context contains $(length(missing)) moment(s) not generated by the translation-invariant PSD blocks: [$shown]. Increase `order`, adjust `basis`, or disable incompatible symmetry reductions."
    ))
end

function _check_objective_moments_covered(
    objective::Polynomial{PauliAlgebra,T,C},
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
) where {T<:Integer,C<:Number}
    return _check_polynomial_moments_covered(objective, moment_basis, "Objective")
end

function _translation_reduced_constraint_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer,
    ::Type{MP};
    context::AbstractString,
) where {T<:Unsigned,C<:Number,CMP<:Number,MP<:Polynomial{PauliAlgebra,T,CMP}}
    reduced = _translation_orbit_reduce_polynomial(poly, n_sites)
    CMP <: Real && _check_real_pauli_chain_polynomial(reduced; context)
    return convert(MP, reduced)
end

function _append_translation_scalar_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
    zero_origin_by_constraint::Dict{Int,Any},
) where {
    T<:Unsigned,
    CMP<:Number,
    CPP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
    P<:Polynomial{PauliAlgebra,T,CPP},
}
    for (idx, poly) in pairs(pop.eq_constraints)
        reduced = _translation_reduced_constraint_polynomial(
            poly,
            n_sites,
            MP;
            context="equality constraint $idx",
        )
        _check_polynomial_moments_covered(reduced, moment_basis, "Equality constraint $idx")
        before = length(constraints)
        _append_constraint!(constraints, :Zero, reshape([reduced], 1, 1), MP)
        zero_origin_by_constraint[before + 1] = TranslationZeroOriginSeed(
            (feature=:scalar_equality, index=idx),
        )
    end

    psd_cone = CMP <: Real ? :PSD : :HPSD
    for (idx, poly) in pairs(pop.ineq_constraints)
        reduced = _translation_reduced_constraint_polynomial(
            poly,
            n_sites,
            MP;
            context="inequality constraint $idx",
        )
        _check_polynomial_moments_covered(reduced, moment_basis, "Inequality constraint $idx")
        _append_constraint!(constraints, psd_cone, reshape([reduced], 1, 1), MP)
        push!(logical_block_sizes, 1)
        push!(block_sizes, 1)
        push!(block_labels, (feature=:scalar_inequality, index=idx))
        push!(block_logical_row_labels, Any[one(NormalMonomial{PauliAlgebra,T})])
        push!(block_transforms, nothing)
    end
    return nothing
end

function _translation_moment_eq_row_bases(
    orbit_reps::Vector{NormalMonomial{PauliAlgebra,T}},
    ;
    sign_symmetry::Bool,
) where {T<:Unsigned}
    row_bases = if sign_symmetry
        NormalMonomial{PauliAlgebra,T}[
            mono for mono in orbit_reps if _pauli_sign_signature(mono) == 0x00
        ]
    else
        copy(orbit_reps)
    end
    sort!(row_bases; by=mono -> (degree(mono), mono))
    return row_bases, degree.(row_bases)
end

function _append_translation_moment_eq_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    pop::PolyOpt{PauliAlgebra,T,P},
    n_sites::Integer,
    orbit_reps::Vector{NormalMonomial{PauliAlgebra,T}},
    sign_symmetry::Bool,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    zero_origin_by_constraint::Dict{Int,Any},
) where {
    T<:Unsigned,
    CMP<:Number,
    CPP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
    P<:Polynomial{PauliAlgebra,T,CPP},
}
    isempty(pop.moment_eq_constraints) && return nothing

    row_bases, row_degrees = _translation_moment_eq_row_bases(
        orbit_reps;
        sign_symmetry,
    )
    one_mono = one(NormalMonomial{PauliAlgebra,T})
    work_coeff_type = CMP <: Real ? Complex{CMP} : CMP
    buf = T[]

    for (idx, g) in pairs(pop.moment_eq_constraints)
        rows = _truncate_moment_eq_row_bases(row_bases, row_degrees, g)
        isempty(rows) && continue

        for (row_idx, row_mono) in pairs(rows)
            terms = Tuple{work_coeff_type,NormalMonomial{PauliAlgebra,T}}[]
            sizehint!(terms, length(row_mono) * length(g.terms))
            for (c_row, row_word) in row_mono
                conj_row = _conj_coef(PauliAlgebra, c_row)
                for (coef, mono) in g.terms
                    _push_scaled_buffered_terms!(
                        terms,
                        conj_row * coef,
                        PauliAlgebra,
                        simplify!(PauliAlgebra, _neat_dot3!(buf, row_word, mono, one_mono)),
                        T,
                        work_coeff_type,
                    )
                end
            end

            poly = _polynomial_from_owned_terms!(terms)
            iszero(poly) && continue
            context = "Moment equality constraint $idx"
            reduced = _translation_reduced_constraint_polynomial(poly, n_sites, MP; context)
            _check_polynomial_moments_covered(reduced, moment_basis, context)
            before = length(constraints)
            _append_constraint!(constraints, :Zero, reshape([reduced], 1, 1), MP)
            zero_origin_by_constraint[before + 1] = TranslationZeroOriginSeed(
                (feature=:moment_equality, index=idx, row=row_idx, row_monomial=row_mono),
            )
        end
    end
    return nothing
end

function _append_translation_axis_rotation_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    n_sites::Integer,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    zero_origin_by_constraint::Dict{Int,Any},
) where {
    T<:Unsigned,
    CMP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
}
    rows = _translation_axis_rotation_equality_rows(ops, n_sites, moment_basis)
    for (row_idx, row) in pairs(rows)
        poly = _translation_axis_rotation_equality_polynomial(row, MP)
        iszero(poly) && continue
        before = length(constraints)
        _append_constraint!(constraints, :Zero, reshape([poly], 1, 1), MP)
        zero_origin_by_constraint[before + 1] = TranslationZeroOriginSeed(
            (
                feature=:axis_rotation_equality,
                generator=row.generator,
                row=row_idx,
                source=row.source,
                target=row.target,
                coefficient=row.coefficient,
            ),
        )
    end
    return nothing
end

_normalize_linear_state_opt_width(::Nothing) = 0

function _normalize_linear_state_opt_width(width::Integer)
    w = Int(width)
    w >= 0 || throw(DomainError(width, "`linear_state_opt_width` must be non-negative."))
    return w
end

function _normalize_linear_state_opt_mode(mode::Symbol)
    mode in (:contiguous, :qmbcertify) || throw(ArgumentError(
        "`linear_state_opt_mode` must be `:contiguous` or `:qmbcertify`; got $mode."
    ))
    return mode
end

function _append_translation_linear_state_opt_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    hamiltonian::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer,
    test_width::Int,
    sign_symmetry::Bool,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{MP},
    zero_origin_by_constraint::Dict{Int,Any},
    ;
    linear_state_opt_mode::Symbol=:contiguous,
) where {
    T<:Unsigned,
    C<:Number,
    CMP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
}
    test_width == 0 && return nothing

    M = NormalMonomial{PauliAlgebra,T}
    mode = _normalize_linear_state_opt_mode(linear_state_opt_mode)
    seen_rows = mode == :qmbcertify ? Set{MP}() : nothing
    moment_set = Set(moment_basis)
    for test_mono in _linear_state_opt_tests(M, Int(n_sites), test_width; sign_symmetry, mode)
        row = im * (hamiltonian * test_mono - test_mono * hamiltonian)
        iszero(row) && continue
        context = "Linear state-opt width $test_width row"
        reduced = _translation_reduced_constraint_polynomial(row, n_sites, MP; context)
        iszero(reduced) && continue
        covered = all(mono -> mono in moment_set, monomials(reduced))
        if !covered
            mode == :qmbcertify && continue
            _check_polynomial_moments_covered(reduced, moment_basis, context)
        end
        if seen_rows !== nothing
            row_key = _qmbcertify_linear_state_opt_row_key(reduced)
            row_key in seen_rows && continue
            push!(seen_rows, row_key)
        end
        before = length(constraints)
        _append_constraint!(constraints, :Zero, reshape([reduced], 1, 1), MP)
        zero_origin_by_constraint[before + 1] = TranslationZeroOriginSeed(
            (
                feature=:linear_state_opt,
                width=test_width,
                mode=mode,
                test_monomial=test_mono,
            ),
        )
    end
    return nothing
end

_normalize_psd_state_opt_width(::Nothing) = 0

function _normalize_psd_state_opt_width(width::Integer)
    w = Int(width)
    w >= 0 || throw(DomainError(width, "`psd_state_opt_width` must be non-negative."))
    return w
end

function _translation_psd_state_opt_entry(
    row::M,
    col::M,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    k::Int,
    n::Int,
    col_translates::Vector{M},
    rep_cache::Dict{M,M},
    ::Type{P},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T},P<:Polynomial{PauliAlgebra,T}}
    C = _coefficient_type(P)
    R = typeof(real(one(C)))
    terms = Tuple{C,M}[]
    sizehint!(terms, length(col_translates) * max(1, length(hamiltonian.terms)))

    for r in 0:(n - 1)
        phase = _momentum_phase(R, k, r, n)
        translated_col = col_translates[r + 1]
        entry = row * hamiltonian * translated_col -
            (1 // 2) * (hamiltonian * row * translated_col + row * translated_col * hamiltonian)
        for (coef, mono) in entry
            iszero(coef) && continue
            rep = get!(rep_cache, mono) do
                _translation_orbit_representative(mono, n)
            end
            push!(terms, (convert(C, phase * coef), rep))
        end
    end
    return _polynomial_from_owned_terms!(terms)
end

function _translation_psd_state_opt_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{P},
) where {T<:Unsigned,H<:Number,M<:NormalMonomial{PauliAlgebra,T},P<:Polynomial{PauliAlgebra,T}}
    mat = fill(zero(P), length(block_basis), length(block_basis))
    for (col_idx, col) in pairs(block_basis), (row_idx, row) in pairs(block_basis)
        mat[row_idx, col_idx] = _translation_psd_state_opt_entry(
            row,
            col,
            hamiltonian,
            k,
            n,
            translated[col],
            rep_cache,
            P,
        )
    end
    return mat
end

function _append_translation_psd_state_opt_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    hamiltonian::Polynomial{PauliAlgebra,T,H},
    n_sites::Integer,
    test_width::Int,
    sectors::Vector{Int},
    sign_symmetry::Bool,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{BLOCK_P},
    ::Type{MP},
    real_moment_matrix::Bool,
    phase_atol::R,
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
) where {
    T<:Unsigned,
    H<:Number,
    R<:Real,
    CMP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
    BLOCK_P<:Polynomial{PauliAlgebra,T},
}
    test_width == 0 && return nothing

    n = Int(n_sites)
    M = NormalMonomial{PauliAlgebra,T}
    block_basis = _contiguous_state_opt_tests(M, n, test_width; sign_symmetry)
    isempty(block_basis) && return nothing

    translated = Dict{M,Vector{M}}()
    for rep in block_basis
        translated[rep] = [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
    end
    rep_cache = Dict{M,M}()

    for k in sectors
        complex_mat = _translation_psd_state_opt_block(
            block_basis,
            k,
            n,
            hamiltonian,
            translated,
            rep_cache,
            BLOCK_P,
        )
        complex_mat = _hermitianize_pauli_polynomial_block(complex_mat, BLOCK_P)
        context = "PSD state-opt width=$test_width momentum=$k block"
        _check_polynomial_matrix_moments_covered(complex_mat, moment_basis, context)
        cone, mat = real_moment_matrix ?
            (:PSD, _realify_hermitian_block(complex_mat, MP; atol=phase_atol)) :
            (:HPSD, map(p -> convert(MP, p), complex_mat))
        _append_constraint!(constraints, cone, mat, MP)
        push!(logical_block_sizes, size(complex_mat, 1))
        push!(block_sizes, size(mat, 1))
        push!(block_labels, (feature=:psd_state_opt, width=test_width, momentum=k))
        push!(block_logical_row_labels, Any[mono for mono in block_basis])
        push!(block_transforms, TranslationDFTTransform(n, k, real_moment_matrix))
    end
    return nothing
end

function _translation_zero_component_label(label, component::Symbol)
    label === nothing && return nothing
    label isa NamedTuple && return merge(label, (component=component,))
    return (feature=:translation_zero, label=label, component=component)
end

function _record_translation_zero_origin!(
    zero_origin_by_constraint,
    constraint_idx::Integer,
    label,
    seed_factory=TranslationZeroOriginSeed,
)
    (zero_origin_by_constraint === nothing || label === nothing) && return nothing
    zero_origin_by_constraint[Int(constraint_idx)] = seed_factory(label)
    return nothing
end

_normalize_contiguous_rdm_k(::Nothing) = Int[]

function _normalize_contiguous_rdm_k(k::Integer)
    return Int[Int(k)]
end

function _normalize_contiguous_rdm_k(ks)
    return Int[Int(k) for k in ks]
end

function _normalize_contiguous_rdm_decomposition(
    decomposition::Symbol,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
)
    decomposition in (:full, :u1, :su2, :qmbcertify) || throw(ArgumentError(
        "`contiguous_rdm_decomposition` must be `:full`, `:u1`, `:su2`, or `:qmbcertify`; got $decomposition."
    ))
    decomposition == :u1 && !u1_symmetry && throw(ArgumentError(
        "`contiguous_rdm_decomposition=:u1` requires `u1_symmetry=true` so the off-sector RDM entries are valid zero constraints."
    ))
    decomposition == :su2 && !su2_symmetry && throw(ArgumentError(
        "`contiguous_rdm_decomposition=:su2` requires `su2_symmetry=true` so the Schur-basis RDM reduction is valid."
    ))
    return decomposition
end

function _normalize_contiguous_rdm_support(support::Symbol)
    support in (:closed, :extend) || throw(ArgumentError(
        "`contiguous_rdm_support` must be `:closed` or `:extend`; got $support."
    ))
    return support
end

@inline _state_bit(state::Int, offset::Int) = (state >> offset) & 1

function _contiguous_rdm_monomial(
    ::Type{M},
    pauli_codes::Vector{Int},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    word = T[]
    sizehint!(word, length(pauli_codes))
    for (offset, pauli_code) in pairs(pauli_codes)
        pauli_code == 0 && continue
        push!(word, _pauli_index_from_site_type(T, offset, pauli_code - 1))
    end
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError(
        "Internal error: contiguous RDM monomial enumeration produced non-real phase $phase."
    ))
    return M(copy(simplified))
end

function _contiguous_rdm_reduced_monomials(
    n_sites::Integer,
    k::Int,
    ::Type{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    pauli_codes = fill(0, k)
    reduced = Vector{M}(undef, 4^k)
    for code0 in 0:(4^k - 1)
        code = code0
        for idx in 1:k
            pauli_codes[idx] = mod(code, 4)
            code = div(code, 4)
        end
        mono = _contiguous_rdm_monomial(M, pauli_codes)
        reduced[code0 + 1] = _translation_orbit_representative(mono, n_sites)
    end
    return reduced
end

function _contiguous_rdm_reduced_orbit_data(
    n_sites::Integer,
    k::Int,
    ::Type{M},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T}}
    reduced_by_code = _contiguous_rdm_reduced_monomials(n_sites, k, M)
    orbit_monos = M[]
    orbit_lookup = Dict{M,Int}()
    orbit_indices = Vector{Int}(undef, length(reduced_by_code))
    for (code_idx, mono) in pairs(reduced_by_code)
        orbit_indices[code_idx] = get!(orbit_lookup, mono) do
            push!(orbit_monos, mono)
            length(orbit_monos)
        end
    end
    return orbit_monos, orbit_indices
end

function _contiguous_rdm_entry_terms(
    reduced_by_code::Vector{M},
    k::Int,
    row_state::Int,
    col_state::Int,
    ::Type{C},
    scale::C,
) where {C<:Number,M<:NormalMonomial}
    terms = Tuple{C,M}[]
    sizehint!(terms, 1 << k)
    flip_mask = xor(row_state, col_state)
    imag_unit = im * one(C)
    for choice0 in 0:((1 << k) - 1)
        choice = choice0
        code0 = 0
        place = 1
        coeff = scale
        for offset in 0:(k - 1)
            choose_second = !iszero(choice & 1)
            choice >>= 1
            row_bit = _state_bit(row_state, offset)
            if iszero((flip_mask >> offset) & 1)
                if choose_second
                    code = 3
                    coeff *= row_bit == 0 ? one(C) : -one(C)
                else
                    code = 0
                end
            else
                if choose_second
                    code = 2
                    coeff *= row_bit == 0 ? -imag_unit : imag_unit
                else
                    code = 1
                end
            end
            code0 += code * place
            place *= 4
        end
        push!(terms, (coeff, reduced_by_code[code0 + 1]))
    end
    return terms
end

function _pauli_rdm_code_action(
    ::Type{C},
    code0::Int,
    k::Int,
    col_state::Int,
) where {C<:Number}
    row_state = col_state
    coeff = one(C)
    code = code0
    imag_unit = im * one(C)
    for offset in 0:(k - 1)
        pauli_code = code & 0x03
        code >>= 2
        bit = _state_bit(col_state, offset)
        if pauli_code == 0
            continue
        elseif pauli_code == 1
            row_state = xor(row_state, 1 << offset)
        elseif pauli_code == 2
            row_state = xor(row_state, 1 << offset)
            coeff *= bit == 0 ? imag_unit : -imag_unit
        elseif pauli_code == 3
            coeff *= bit == 0 ? one(C) : -one(C)
        else
            throw(ArgumentError("Invalid Pauli RDM code $pauli_code."))
        end
    end
    return row_state, coeff
end

function _pauli_rdm_code_actions(::Type{C}, k::Int) where {C<:Number}
    count = 4^k
    flip_masks = Vector{Int}(undef, count)
    sign_masks = Vector{Int}(undef, count)
    phases = Vector{C}(undef, count)
    imag_unit = im * one(C)

    for code0 in 0:(count - 1)
        code = code0
        flip_mask = 0
        sign_mask = 0
        phase = one(C)
        for offset in 0:(k - 1)
            pauli_code = code & 0x03
            code >>= 2
            if pauli_code == 1
                flip_mask |= 1 << offset
            elseif pauli_code == 2
                flip_mask |= 1 << offset
                sign_mask |= 1 << offset
                phase *= imag_unit
            elseif pauli_code == 3
                sign_mask |= 1 << offset
            end
        end
        flip_masks[code0 + 1] = flip_mask
        sign_masks[code0 + 1] = sign_mask
        phases[code0 + 1] = phase
    end
    return (flip_masks=flip_masks, sign_masks=sign_masks, phases=phases)
end

function _pauli_rdm_state_sign_parities(k::Int)
    dim = 1 << k
    parities = Matrix{Bool}(undef, dim, dim)
    @inbounds for sign_mask in 0:(dim - 1), state0 in 0:(dim - 1)
        parities[state0 + 1, sign_mask + 1] = isodd(count_ones(sign_mask & state0))
    end
    return parities
end

function _sparse_transform_columns_by_state(
    transform::_SparseTransformRows{TC},
) where {TC<:Number}
    rows_by_state = [Tuple{Int,TC}[] for _ in 1:transform.ncols]
    for (local_idx, row) in pairs(transform.rows)
        for (state_idx, coeff) in row
            push!(rows_by_state[state_idx], (local_idx, coeff))
        end
    end
    return rows_by_state
end

const _PAULI_SU2_RDM_DENSE_ACCUMULATOR_SLOT_LIMIT = 8_000_000

function _translation_contiguous_rdm_su2_reduced_linear_block_dense(
    orbit_monos::Vector{M},
    orbit_indices::Vector{Int},
    k::Int,
    transform::_SparseTransformRows,
    rows_by_state,
    code_actions,
    ::Type{K},
    ::Type{C};
    atol::Real,
    stage_times_ns=nothing,
    state_sign_parities=_pauli_rdm_state_sign_parities(k),
) where {K,C<:Number,M<:NormalMonomial}
    block_dim = length(transform.rows)
    orbit_count = length(orbit_monos)
    entry_count = block_dim * block_dim
    values = zeros(C, entry_count, orbit_count)
    seen = falses(entry_count, orbit_count)
    touched = [Int[] for _ in 1:entry_count]
    R = typeof(real(one(C)))
    scale = C(inv(R(2)^k))
    tolerance = R(atol)

    stage_start_ns = time_ns()
    for code0 in 0:(length(orbit_indices) - 1)
        action_idx = code0 + 1
        @inbounds orbit_idx = orbit_indices[action_idx]
        @inbounds flip_mask = code_actions.flip_masks[action_idx]
        @inbounds sign_mask = code_actions.sign_masks[action_idx]
        @inbounds phase = code_actions.phases[action_idx]
        for (col_local, right_row) in pairs(transform.rows)
            entry_col_offset = (col_local - 1) * block_dim
            for (col_state, right_coeff) in right_row
                col_state0 = col_state - 1
                row_state0 = xor(col_state0, flip_mask)
                @inbounds is_negative = state_sign_parities[col_state0 + 1, sign_mask + 1]
                matrix_entry = is_negative ? -phase : phase
                @inbounds row_entries = rows_by_state[row_state0 + 1]
                isempty(row_entries) && continue
                for (row_local, left_coeff) in row_entries
                    value = scale * convert(C, left_coeff * conj(right_coeff)) * matrix_entry
                    abs(value) <= tolerance && continue
                    entry_idx = entry_col_offset + row_local
                    if !seen[entry_idx, orbit_idx]
                        seen[entry_idx, orbit_idx] = true
                        push!(touched[entry_idx], orbit_idx)
                    end
                    @inbounds values[entry_idx, orbit_idx] += value
                end
            end
        end
    end
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_rdm_block_accumulate] =
            get(stage_times_ns, :su2_extend_rdm_block_accumulate, 0) +
            Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    entries = Matrix{LinearMomentForm{K,C}}(undef, block_dim, block_dim)
    for entry_idx in 1:entry_count
        pairs = Pair{K,C}[]
        sizehint!(pairs, length(touched[entry_idx]))
        for orbit_idx in touched[entry_idx]
            @inbounds coef = values[entry_idx, orbit_idx]
            abs(coef) <= tolerance && continue
            @inbounds push!(pairs, _moment_key(K, orbit_monos[orbit_idx]) => coef)
        end
        entries[entry_idx] = _linear_moment_form_from_owned_pairs!(pairs)
    end
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_rdm_block_finalize] =
            get(stage_times_ns, :su2_extend_rdm_block_finalize, 0) +
            Int(time_ns() - stage_start_ns)
    end
    return entries
end

function _translation_contiguous_rdm_su2_reduced_linear_block_dict(
    orbit_monos::Vector{M},
    orbit_indices::Vector{Int},
    k::Int,
    transform::_SparseTransformRows,
    rows_by_state,
    code_actions,
    ::Type{K},
    ::Type{C};
    atol::Real,
    stage_times_ns=nothing,
    state_sign_parities=_pauli_rdm_state_sign_parities(k),
) where {K,C<:Number,M<:NormalMonomial}
    block_dim = length(transform.rows)
    accumulators = [Dict{Int,C}() for _ in 1:block_dim, _ in 1:block_dim]
    R = typeof(real(one(C)))
    scale = C(inv(R(2)^k))
    tolerance = R(atol)

    stage_start_ns = time_ns()
    for code0 in 0:(length(orbit_indices) - 1)
        orbit_idx = orbit_indices[code0 + 1]
        flip_mask = code_actions.flip_masks[code0 + 1]
        sign_mask = code_actions.sign_masks[code0 + 1]
        phase = code_actions.phases[code0 + 1]
        for (col_local, right_row) in pairs(transform.rows)
            for (col_state, right_coeff) in right_row
                col_state0 = col_state - 1
                row_state0 = xor(col_state0, flip_mask)
                @inbounds is_negative = state_sign_parities[col_state0 + 1, sign_mask + 1]
                matrix_entry = is_negative ? -phase : phase
                row_entries = rows_by_state[row_state0 + 1]
                isempty(row_entries) && continue
                for (row_local, left_coeff) in row_entries
                    value = scale * convert(C, left_coeff * conj(right_coeff)) * matrix_entry
                    abs(value) <= tolerance && continue
                    acc = accumulators[row_local, col_local]
                    acc[orbit_idx] = get(acc, orbit_idx, zero(C)) + value
                end
            end
        end
    end
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_rdm_block_accumulate] =
            get(stage_times_ns, :su2_extend_rdm_block_accumulate, 0) +
            Int(time_ns() - stage_start_ns)
    end

    stage_start_ns = time_ns()
    entries = Matrix{LinearMomentForm{K,C}}(undef, block_dim, block_dim)
    for idx in eachindex(accumulators)
        acc = accumulators[idx]
        pairs = Pair{K,C}[]
        sizehint!(pairs, length(acc))
        for (orbit_idx, coef) in acc
            abs(coef) <= tolerance && continue
            push!(pairs, _moment_key(K, orbit_monos[orbit_idx]) => coef)
        end
        entries[idx] = _linear_moment_form_from_owned_pairs!(pairs)
    end
    if stage_times_ns !== nothing
        stage_times_ns[:su2_extend_rdm_block_finalize] =
            get(stage_times_ns, :su2_extend_rdm_block_finalize, 0) +
            Int(time_ns() - stage_start_ns)
    end
    return entries
end

function _translation_contiguous_rdm_su2_reduced_linear_block(
    orbit_monos::Vector{M},
    orbit_indices::Vector{Int},
    k::Int,
    transform::_SparseTransformRows,
    ::Type{K},
    ::Type{C};
    atol::Real,
    stage_times_ns=nothing,
) where {K,C<:Number,M<:NormalMonomial}
    return _translation_contiguous_rdm_su2_reduced_linear_block(
        orbit_monos,
        orbit_indices,
        k,
        transform,
        _pauli_rdm_code_actions(C, k),
        K,
        C;
        atol,
        stage_times_ns,
    )
end

function _translation_contiguous_rdm_su2_reduced_linear_block(
    orbit_monos::Vector{M},
    orbit_indices::Vector{Int},
    k::Int,
    transform::_SparseTransformRows,
    code_actions,
    ::Type{K},
    ::Type{C};
    atol::Real,
    stage_times_ns=nothing,
) where {K,C<:Number,M<:NormalMonomial}
    block_dim = length(transform.rows)
    rows_by_state = _sparse_transform_columns_by_state(transform)
    state_sign_parities = _pauli_rdm_state_sign_parities(k)
    dense_slots = length(orbit_monos) * block_dim * block_dim
    if dense_slots <= _PAULI_SU2_RDM_DENSE_ACCUMULATOR_SLOT_LIMIT
        return _translation_contiguous_rdm_su2_reduced_linear_block_dense(
            orbit_monos,
            orbit_indices,
            k,
            transform,
            rows_by_state,
            code_actions,
            K,
            C;
            atol,
            stage_times_ns,
            state_sign_parities,
        )
    end

    return _translation_contiguous_rdm_su2_reduced_linear_block_dict(
        orbit_monos,
        orbit_indices,
        k,
        transform,
        rows_by_state,
        code_actions,
        K,
        C;
        atol,
        stage_times_ns,
        state_sign_parities,
    )
end

function _translation_contiguous_rdm_block(
    n_sites::Integer,
    k::Integer,
    ::Type{P},
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    kk = Int(k)
    1 <= kk <= Int(n_sites) || throw(DomainError(k, "`contiguous_rdm_k` must satisfy 1 <= k <= n_sites."))

    M = NormalMonomial{PauliAlgebra,T}
    dim = 1 << kk
    mat = fill(zero(P), dim, dim)
    R = typeof(real(one(C)))
    scale = C(inv(R(2)^kk))
    reduced_by_code = _contiguous_rdm_reduced_monomials(n_sites, kk, M)

    for col_state in 0:(dim - 1), row_state in 0:(dim - 1)
        mat[row_state + 1, col_state + 1] = P(
            _contiguous_rdm_entry_terms(reduced_by_code, kk, row_state, col_state, C, scale),
        )
    end
    return mat
end

function _check_polynomial_matrix_moments_covered(
    mat::AbstractMatrix{<:Polynomial{PauliAlgebra,T}},
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    context::AbstractString,
) where {T<:Unsigned}
    for entry in mat
        _check_polynomial_moments_covered(entry, moment_basis, context)
    end
    return nothing
end

function _append_translation_zero_scalar_constraint!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    scalar_poly::Polynomial{PauliAlgebra,T},
    ::Type{MP};
    zero_origin_by_constraint,
    origin_label,
    origin_seed_factory=TranslationZeroOriginSeed,
) where {T<:Unsigned,MP<:Polynomial{PauliAlgebra,T}}
    row = convert(MP, scalar_poly)
    iszero(row) && return nothing

    row_adjoint = adjoint(row)
    if iszero(row - row_adjoint)
        before = length(constraints)
        push!(constraints, (:Zero, reshape([row], 1, 1)))
        _record_translation_zero_origin!(
            zero_origin_by_constraint,
            before + 1,
            origin_label,
            origin_seed_factory,
        )
        return nothing
    end

    hermitian_part = convert(MP, (row + row_adjoint) * (1 // 2))
    if !iszero(hermitian_part)
        before = length(constraints)
        push!(constraints, (:Zero, reshape([hermitian_part], 1, 1)))
        _record_translation_zero_origin!(
            zero_origin_by_constraint,
            before + 1,
            origin_label,
            origin_seed_factory,
        )
    end

    skewhermitian_part = convert(MP, (row - row_adjoint) / (2im))
    if !iszero(skewhermitian_part)
        before = length(constraints)
        push!(constraints, (:Zero, reshape([skewhermitian_part], 1, 1)))
        _record_translation_zero_origin!(
            zero_origin_by_constraint,
            before + 1,
            origin_label,
            origin_seed_factory,
        )
    end
    return nothing
end

function _append_translation_zero_polynomial!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    poly::Polynomial{PauliAlgebra,T},
    ::Type{MP};
    phase_atol,
    zero_origin_by_constraint=nothing,
    label=nothing,
    origin_seed_factory=TranslationZeroOriginSeed,
) where {T<:Unsigned,CMP<:Number,MP<:Polynomial{PauliAlgebra,T,CMP}}
    if CMP <: Real
        re = _real_part_polynomial(poly, MP; atol=phase_atol)
        _append_translation_zero_scalar_constraint!(
            constraints,
            re,
            MP;
            zero_origin_by_constraint,
            origin_label=_translation_zero_component_label(label, :real),
            origin_seed_factory,
        )
        im = _imag_part_polynomial(poly, MP; atol=phase_atol)
        _append_translation_zero_scalar_constraint!(
            constraints,
            im,
            MP;
            zero_origin_by_constraint,
            origin_label=_translation_zero_component_label(label, :imag),
            origin_seed_factory,
        )
    else
        _append_translation_zero_scalar_constraint!(
            constraints,
            poly,
            MP;
            zero_origin_by_constraint,
            origin_label=_translation_zero_component_label(label, :complex),
            origin_seed_factory,
        )
    end
    return nothing
end

@inline _contiguous_rdm_magnetization(state::Integer) = count_ones(Int(state))

function _contiguous_rdm_u1_sectors(k::Integer)
    kk = Int(k)
    sectors = [Int[] for _ in 0:kk]
    for state in 0:((1 << kk) - 1)
        push!(sectors[_contiguous_rdm_magnetization(state) + 1], state + 1)
    end
    return sectors
end

function _contiguous_rdm_state_labels(k::Integer, rows::AbstractVector{<:Integer})
    kk = Int(k)
    return Any[(feature=:contiguous_rdm_state, k=kk, state=Int(row) - 1) for row in rows]
end

struct _PauliSU2SchurState
    spin2::Int
    m2::Int
    multiplicity::Int
    vector::Vector{Float64}
end

function _pauli_su2_lift_vector(
    vector::Vector{Float64},
    previous_sites::Int,
    new_bit::Int,
)
    out = zeros(Float64, 1 << (previous_sites + 1))
    offset = new_bit << previous_sites
    for state in 0:((1 << previous_sites) - 1)
        out[offset + state + 1] = vector[state + 1]
    end
    return out
end

function _pauli_su2_add_scaled!(dest::Vector{Float64}, scale::Float64, src::Vector{Float64})
    iszero(scale) && return dest
    @inbounds for idx in eachindex(dest, src)
        dest[idx] += scale * src[idx]
    end
    return dest
end

function _pauli_su2_spin_half_coupling_coefficients(
    parent_spin2::Integer,
    child_spin2::Integer,
    child_m2::Integer,
)
    parent = Int(parent_spin2)
    child = Int(child_spin2)
    m2 = Int(child_m2)
    parent >= 0 || throw(DomainError(parent_spin2, "`parent_spin2` must be non-negative."))
    child in (parent - 1, parent + 1) || throw(ArgumentError(
        "Spin-1/2 coupling requires child_spin2 == parent_spin2 ± 1; got parent_spin2=$parent, child_spin2=$child."
    ))
    child >= 0 || throw(DomainError(child_spin2, "`child_spin2` must be non-negative."))
    -child <= m2 <= child || throw(DomainError(
        child_m2,
        "`child_m2` must satisfy -child_spin2 <= child_m2 <= child_spin2.",
    ))
    iseven(child - m2) || throw(DomainError(
        child_m2,
        "`child_m2` must have the same parity as child_spin2.",
    ))

    denominator = 2 * (parent + 1)
    up_numerator = child == parent + 1 ? parent + m2 + 1 : parent - m2 + 1
    down_numerator = child == parent + 1 ? parent - m2 + 1 : parent + m2 + 1
    return (
        denominator=denominator,
        up_numerator=up_numerator,
        up_sign=child == parent - 1 ? -1 : 1,
        down_numerator=down_numerator,
        down_sign=1,
    )
end

function _pauli_su2_float_coupling_scale(sign::Integer, numerator::Integer, denominator::Integer)
    return Int(sign) * sqrt(max(0.0, Int(numerator) / Int(denominator)))
end

function _pauli_su2_coupled_vector(
    states_by_m::Dict{Int,_PauliSU2SchurState},
    previous_sites::Int,
    parent_spin2::Int,
    child_spin2::Int,
    child_m2::Int,
)
    coeffs = _pauli_su2_spin_half_coupling_coefficients(
        parent_spin2,
        child_spin2,
        child_m2,
    )
    vector = zeros(Float64, 1 << (previous_sites + 1))

    up_parent = get(states_by_m, child_m2 - 1, nothing)
    if up_parent !== nothing
        scale = _pauli_su2_float_coupling_scale(
            coeffs.up_sign,
            coeffs.up_numerator,
            coeffs.denominator,
        )
        _pauli_su2_add_scaled!(
            vector,
            scale,
            _pauli_su2_lift_vector(up_parent.vector, previous_sites, 0),
        )
    end

    down_parent = get(states_by_m, child_m2 + 1, nothing)
    if down_parent !== nothing
        scale = _pauli_su2_float_coupling_scale(
            coeffs.down_sign,
            coeffs.down_numerator,
            coeffs.denominator,
        )
        _pauli_su2_add_scaled!(
            vector,
            scale,
            _pauli_su2_lift_vector(down_parent.vector, previous_sites, 1),
        )
    end

    nrm = norm(vector)
    nrm > 0 || throw(ArgumentError("Internal error: SU(2) coupling produced a zero Schur vector."))
    abs(nrm - 1.0) <= 1e-12 || (vector ./= nrm)
    return vector
end

function _pauli_su2_couple_next(
    states::Vector{_PauliSU2SchurState},
    previous_sites::Int,
)
    multiplets = Dict{Tuple{Int,Int},Dict{Int,_PauliSU2SchurState}}()
    for state in states
        by_m = get!(multiplets, (state.spin2, state.multiplicity)) do
            Dict{Int,_PauliSU2SchurState}()
        end
        by_m[state.m2] = state
    end

    child_counts = Dict{Int,Int}()
    children = _PauliSU2SchurState[]
    for (parent_spin2, parent_mult) in sort!(collect(keys(multiplets)))
        states_by_m = multiplets[(parent_spin2, parent_mult)]
        for child_spin2 in (parent_spin2 + 1, parent_spin2 - 1)
            child_spin2 < 0 && continue
            child_mult = get(child_counts, child_spin2, 0) + 1
            child_counts[child_spin2] = child_mult
            for child_m2 in -child_spin2:2:child_spin2
                vector = _pauli_su2_coupled_vector(
                    states_by_m,
                    previous_sites,
                    parent_spin2,
                    child_spin2,
                    child_m2,
                )
                push!(
                    children,
                    _PauliSU2SchurState(child_spin2, child_m2, child_mult, vector),
                )
            end
        end
    end

    return sort!(children; by=state -> (state.spin2, state.m2, state.multiplicity))
end

function _pauli_su2_schur_states(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))
    kk == 0 && return _PauliSU2SchurState[_PauliSU2SchurState(0, 0, 1, [1.0])]

    states = _PauliSU2SchurState[
        _PauliSU2SchurState(1, -1, 1, [0.0, 1.0]),
        _PauliSU2SchurState(1, 1, 1, [1.0, 0.0]),
    ]
    for previous_sites in 1:(kk - 1)
        states = _pauli_su2_couple_next(states, previous_sites)
    end
    return sort!(states; by=state -> (state.spin2, state.m2, state.multiplicity))
end

function _pauli_su2_schur_matrix(k::Integer)
    states = _pauli_su2_schur_states(k)
    dim = 1 << Int(k)
    mat = Matrix{Float64}(undef, dim, length(states))
    for (col, state) in pairs(states)
        mat[:, col] = state.vector
    end
    return mat, states
end

"""
    pauli_su2_schur_diagnostics(k)

Return numerical diagnostics for the generated `k`-qubit SU(2) Schur basis:
dimension, state count, spin multiplicities, and the maximum entrywise
orthonormality residual of the transform.  The `sz_residual` and
`casimir_residual` fields check that the same basis diagonalizes total `J_z`
and `J^2`.
"""
function pauli_su2_schur_diagnostics(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    transform, states = _pauli_su2_schur_matrix(kk)
    gram_residual = transpose(transform) * transform - I
    sz_residual, casimir_residual = _pauli_su2_spin_residuals(states, kk)
    spin_counts = Dict{Int,Set{Int}}()
    for state in states
        multiplicities = get!(spin_counts, state.spin2) do
            Set{Int}()
        end
        push!(multiplicities, state.multiplicity)
    end
    spin_multiplicities = [
        spin2 => length(spin_counts[spin2])
        for spin2 in sort!(collect(keys(spin_counts)))
    ]

    return (
        k=kk,
        dimension=size(transform, 1),
        state_count=length(states),
        spin_multiplicities=spin_multiplicities,
        coefficient_domain=:algebraic_float64,
        exact_coefficient_domain=:sqrt_rational,
        unitarity_residual=maximum(abs, gram_residual),
        sz_residual=sz_residual,
        casimir_residual=casimir_residual,
        max_residual=max(maximum(abs, gram_residual), sz_residual, casimir_residual),
    )
end

function _pauli_su2_basis_m2(state::Int, k::Int)
    m2 = 0
    for bit in 0:(k - 1)
        m2 += iszero((state >> bit) & 1) ? 1 : -1
    end
    return m2
end

function _pauli_su2_apply_jplus!(
    out::Vector{Float64},
    vector::Vector{Float64},
    k::Int,
)
    Base.fill!(out, 0.0)
    for state in 0:(length(vector) - 1)
        amp = vector[state + 1]
        iszero(amp) && continue
        for bit in 0:(k - 1)
            ((state >> bit) & 1) == 1 || continue
            out[(state & ~(1 << bit)) + 1] += amp
        end
    end
    return out
end

function _pauli_su2_apply_jminus!(
    out::Vector{Float64},
    vector::Vector{Float64},
    k::Int,
)
    Base.fill!(out, 0.0)
    for state in 0:(length(vector) - 1)
        amp = vector[state + 1]
        iszero(amp) && continue
        for bit in 0:(k - 1)
            ((state >> bit) & 1) == 0 || continue
            out[(state | (1 << bit)) + 1] += amp
        end
    end
    return out
end

function _pauli_su2_spin_residuals(
    states::Vector{_PauliSU2SchurState},
    k::Int,
)
    dim = 1 << k
    jplus = zeros(Float64, dim)
    jminus_jplus = zeros(Float64, dim)
    casimir_image = zeros(Float64, dim)
    max_sz_residual = 0.0
    max_casimir_residual = 0.0

    for state in states
        sz_residual_sq = 0.0
        for basis_state in 0:(dim - 1)
            delta = 0.5 * (_pauli_su2_basis_m2(basis_state, k) - state.m2)
            sz_residual_sq += abs2(delta * state.vector[basis_state + 1])
        end
        max_sz_residual = max(max_sz_residual, sqrt(sz_residual_sq))

        _pauli_su2_apply_jplus!(jplus, state.vector, k)
        _pauli_su2_apply_jminus!(jminus_jplus, jplus, k)
        for basis_state in 0:(dim - 1)
            m = 0.5 * _pauli_su2_basis_m2(basis_state, k)
            casimir_image[basis_state + 1] =
                jminus_jplus[basis_state + 1] + (m^2 + m) * state.vector[basis_state + 1]
        end
        eigenvalue = 0.25 * state.spin2 * (state.spin2 + 2)
        max_casimir_residual = max(
            max_casimir_residual,
            norm(casimir_image .- eigenvalue .* state.vector),
        )
    end
    return max_sz_residual, max_casimir_residual
end

"""
    pauli_su2_spin_diagnostics(k)

Return residuals showing that the generated `k`-qubit Schur basis diagonalizes
the total `J_z` and `J²` operators.
"""
function pauli_su2_spin_diagnostics(k::Integer)
    kk = Int(k)
    kk >= 0 || throw(DomainError(k, "`k` must be non-negative."))

    states = _pauli_su2_schur_states(kk)
    sz_residual, casimir_residual = _pauli_su2_spin_residuals(states, kk)
    return (
        k=kk,
        dimension=1 << kk,
        state_count=length(states),
        sz_residual=sz_residual,
        casimir_residual=casimir_residual,
    )
end

function _pauli_su2_transform_block(
    mat::Matrix{P},
    k::Int,
    ::Type{P},
    transform::Matrix{Float64},
    atol::Real,
) where {P<:Polynomial{PauliAlgebra}}
    dim = size(mat, 1)
    size(mat, 2) == dim || throw(DimensionMismatch("SU(2) RDM block must be square, got $(size(mat))."))
    size(transform, 1) == dim || throw(DimensionMismatch(
        "SU(2) Schur transform for k=$k has $(size(transform, 1)) rows but RDM block has size $dim."
    ))

    out = fill(zero(P), dim, dim)
    for col in 1:dim, row in 1:dim
        acc = zero(P)
        for source_col in 1:dim, source_row in 1:dim
            coef = transform[source_row, row] * transform[source_col, col]
            abs(coef) <= atol && continue
            acc += coef * mat[source_row, source_col]
        end
        out[row, col] = convert(P, acc)
    end
    return out
end

function _pauli_su2_state_columns(states::Vector{_PauliSU2SchurState})
    columns = Dict{Tuple{Int,Int,Int},Int}()
    for (idx, state) in pairs(states)
        columns[(state.spin2, state.m2, state.multiplicity)] = idx
    end
    return columns
end

function _contiguous_rdm_su2_row_labels(k::Integer, block::PauliSU2RDMBlock)
    kk = Int(k)
    return Any[
        (feature=:contiguous_rdm_su2_state, k=kk, spin2=block.spin2, multiplicity=mult)
        for mult in 1:block.multiplicity
    ]
end

function _append_translation_contiguous_rdm_su2_blocks!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    complex_mat::Matrix{BLOCK_P},
    k::Int,
    ::Type{MP},
    real_moment_matrix::Bool,
    phase_atol::R,
    zero_origin_by_constraint::Dict{Int,Any},
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
) where {
    R<:Real,
    MP<:Polynomial{PauliAlgebra},
    BLOCK_P<:Polynomial{PauliAlgebra},
}
    schur_transform, states = _pauli_su2_schur_matrix(k)
    schur_mat = _pauli_su2_transform_block(complex_mat, k, BLOCK_P, schur_transform, phase_atol)
    columns = _pauli_su2_state_columns(states)
    blocks = pauli_su2_rdm_blocks(k)

    for col_state in states, row_state in states
        same_spin_m = row_state.spin2 == col_state.spin2 && row_state.m2 == col_state.m2
        same_spin_m && continue
        _append_translation_zero_polynomial!(
            constraints,
            schur_mat[
                columns[(row_state.spin2, row_state.m2, row_state.multiplicity)],
                columns[(col_state.spin2, col_state.m2, col_state.multiplicity)],
            ],
            MP;
            phase_atol,
            zero_origin_by_constraint,
            label=(
                feature=:contiguous_rdm_zero,
                k=k,
                decomposition=:su2,
                reason=:spin_magnetic_offblock,
                row_spin2=row_state.spin2,
                row_m2=row_state.m2,
                row_multiplicity=row_state.multiplicity,
                col_spin2=col_state.spin2,
                col_m2=col_state.m2,
                col_multiplicity=col_state.multiplicity,
            ),
        )
    end

    reference_rows_by_block = Vector{Vector{Int}}(undef, length(blocks))
    for (block_idx, block) in pairs(blocks)
        reference_m2 = -block.spin2
        reference_rows = [
            columns[(block.spin2, reference_m2, mult)]
            for mult in 1:block.multiplicity
        ]
        reference_rows_by_block[block_idx] = reference_rows
        for m2 in (reference_m2 + 2):2:block.spin2
            rows = [columns[(block.spin2, m2, mult)] for mult in 1:block.multiplicity]
            for col in 1:block.multiplicity, row in 1:block.multiplicity
                _append_translation_zero_polynomial!(
                    constraints,
                    schur_mat[rows[row], rows[col]] -
                    schur_mat[reference_rows[row], reference_rows[col]],
                    MP;
                    phase_atol,
                    zero_origin_by_constraint,
                    label=(
                        feature=:contiguous_rdm_zero,
                        k=k,
                        decomposition=:su2,
                        reason=:magnetic_copy,
                        spin2=block.spin2,
                        m2=m2,
                        reference_m2=reference_m2,
                        row_multiplicity=row,
                        col_multiplicity=col,
                    ),
                )
            end
        end
    end

    for (block, reference_rows) in zip(blocks, reference_rows_by_block)
        complex_block = schur_mat[reference_rows, reference_rows]
        cone, mat = real_moment_matrix ?
            (:PSD, _realify_hermitian_block(complex_block, MP; atol=phase_atol)) :
            (:HPSD, map(p -> convert(MP, p), complex_block))
        _append_constraint!(constraints, cone, mat, MP)
        push!(logical_block_sizes, size(complex_block, 1))
        push!(block_sizes, size(mat, 1))
        push!(
            block_labels,
            (feature=:contiguous_rdm, k=k, decomposition=:su2, spin2=block.spin2),
        )
        push!(block_logical_row_labels, _contiguous_rdm_su2_row_labels(k, block))
        push!(block_transforms, TranslationSU2RDMTransform(k, block, schur_transform))
    end
    return nothing
end

function _append_translation_contiguous_rdm_u1_blocks!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    complex_mat::Matrix{BLOCK_P},
    k::Int,
    ::Type{MP},
    real_moment_matrix::Bool,
    phase_atol::R,
    zero_origin_by_constraint::Dict{Int,Any},
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
) where {
    R<:Real,
    MP<:Polynomial{PauliAlgebra},
    BLOCK_P<:Polynomial{PauliAlgebra},
}
    state_count = size(complex_mat, 1)
    magnetizations = [_contiguous_rdm_magnetization(state) for state in 0:(state_count - 1)]

    for col in 1:state_count, row in 1:state_count
        magnetizations[row] == magnetizations[col] && continue
        _append_translation_zero_polynomial!(
            constraints,
            complex_mat[row, col],
            MP;
            phase_atol,
            zero_origin_by_constraint,
            label=(
                feature=:contiguous_rdm_zero,
                k=k,
                decomposition=:u1,
                reason=:magnetization_offblock,
                row_state=row - 1,
                row_magnetization=magnetizations[row],
                col_state=col - 1,
                col_magnetization=magnetizations[col],
            ),
        )
    end

    for (weight0, sector) in enumerate(_contiguous_rdm_u1_sectors(k))
        complex_block = complex_mat[sector, sector]
        cone, mat = real_moment_matrix ?
            (:PSD, _realify_hermitian_block(complex_block, MP; atol=phase_atol)) :
            (:HPSD, map(p -> convert(MP, p), complex_block))
        _append_constraint!(constraints, cone, mat, MP)
        push!(logical_block_sizes, size(complex_block, 1))
        push!(block_sizes, size(mat, 1))
        push!(
            block_labels,
            (feature=:contiguous_rdm, k=k, decomposition=:u1, magnetization=weight0 - 1),
        )
        push!(block_logical_row_labels, _contiguous_rdm_state_labels(k, sector))
        push!(block_transforms, nothing)
    end
    return nothing
end

function _append_translation_contiguous_rdm_constraints!(
    constraints::Vector{Tuple{Symbol, Matrix{MP}}},
    n_sites::Integer,
    rdm_ks::Vector{Int},
    rdm_decomposition::Symbol,
    rdm_support::Symbol,
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
    ::Type{BLOCK_P},
    ::Type{MP},
    real_moment_matrix::Bool,
    phase_atol::R,
    zero_origin_by_constraint::Dict{Int,Any},
    logical_block_sizes::Vector{Int},
    block_sizes::Vector{Int},
    block_labels::Vector{Any},
    block_logical_row_labels::Vector{Vector{Any}},
    block_transforms::Vector{Any},
) where {
    T<:Unsigned,
    R<:Real,
    CMP<:Number,
    MP<:Polynomial{PauliAlgebra,T,CMP},
    BLOCK_P<:Polynomial{PauliAlgebra,T},
}
    for k in rdm_ks
        complex_mat = _translation_contiguous_rdm_block(n_sites, k, BLOCK_P)
        if rdm_support == :closed
            _check_polynomial_matrix_moments_covered(
                complex_mat,
                moment_basis,
                "Contiguous RDM k=$k block",
            )
        end
        if rdm_decomposition == :full
            cone, mat = real_moment_matrix ?
                (:PSD, _realify_hermitian_block(complex_mat, MP; atol=phase_atol)) :
                (:HPSD, map(p -> convert(MP, p), complex_mat))
            _append_constraint!(constraints, cone, mat, MP)
            push!(logical_block_sizes, size(complex_mat, 1))
            push!(block_sizes, size(mat, 1))
            push!(block_labels, (feature=:contiguous_rdm, k=k))
            push!(
                block_logical_row_labels,
                _contiguous_rdm_state_labels(k, collect(1:size(complex_mat, 1))),
            )
            push!(block_transforms, nothing)
        elseif rdm_decomposition == :u1
            _append_translation_contiguous_rdm_u1_blocks!(
                constraints,
                complex_mat,
                k,
                MP,
                real_moment_matrix,
                phase_atol,
                zero_origin_by_constraint,
                logical_block_sizes,
                block_sizes,
                block_labels,
                block_logical_row_labels,
                block_transforms,
            )
        else
            _append_translation_contiguous_rdm_su2_blocks!(
                constraints,
                complex_mat,
                k,
                MP,
                real_moment_matrix,
                phase_atol,
                zero_origin_by_constraint,
                logical_block_sizes,
                block_sizes,
                block_labels,
                block_logical_row_labels,
                block_transforms,
            )
        end
    end
    return nothing
end

function _collect_translation_moment_basis(
    objective::Polynomial{PauliAlgebra,T,C},
    constraints::Vector{Tuple{Symbol, Matrix{Polynomial{PauliAlgebra,T,C}}}},
) where {T<:Unsigned,C<:Number}
    basis = NormalMonomial{PauliAlgebra,T}[]
    append!(basis, monomials(objective))
    for (_, mat) in constraints, poly in mat
        append!(basis, monomials(poly))
    end
    return sorted_unique!(basis)
end

function _pauli_sign_signature(mono::NormalMonomial{PauliAlgebra})
    px = false
    py = false
    pz = false
    for idx in mono.word
        ptype = _pauli_type(idx)
        if ptype == _PAULI_X_TYPE
            px = !px
        elseif ptype == _PAULI_Y_TYPE
            py = !py
        else
            pz = !pz
        end
    end
    return UInt8((xor(px, py) ? 0x01 : 0x00) | (xor(py, pz) ? 0x02 : 0x00))
end

function _pauli_signature_blocks(basis::Vector{M}) where {M<:NormalMonomial{PauliAlgebra}}
    buckets = Dict{UInt8,Vector{M}}()
    order = UInt8[]
    for mono in basis
        sig = _pauli_sign_signature(mono)
        if !haskey(buckets, sig)
            buckets[sig] = M[]
            push!(order, sig)
        end
        push!(buckets[sig], mono)
    end
    sort!(order)
    return [(sig, buckets[sig]) for sig in order]
end

@inline _translation_reflection_fixed_momentum(k::Int, n::Int) = iszero(k) || 2k == n

function _reflect_pauli_monomial(
    mono::NormalMonomial{PauliAlgebra,T},
    n_sites::Integer,
) where {T<:Unsigned}
    n = Int(n_sites)
    word = Vector{T}(undef, length(mono.word))
    for (idx, pauli_idx) in pairs(mono.word)
        reflected_site = n + 1 - _pauli_site(pauli_idx)
        word[idx] = _pauli_index_from_site_type(T, reflected_site, _pauli_type(pauli_idx))
    end
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError("Internal error: Pauli reflection produced phase $phase."))
    return NormalMonomial{PauliAlgebra,T}(copy(simplified))
end

function _translation_orbit_representative_shift(
    mono::NormalMonomial{PauliAlgebra,T},
    n_sites::Integer,
) where {T<:Unsigned}
    n = Int(n_sites)
    rep = _translation_orbit_representative(mono, n)
    for shift in 0:(n - 1)
        _translate_pauli_monomial(mono, shift, n) == rep && return rep, shift
    end
    throw(ArgumentError("Internal error: failed to locate translation shift for orbit representative $rep."))
end

function _translation_reflection_image(
    mono::NormalMonomial{PauliAlgebra,T},
    k::Int,
    n::Int,
) where {T<:Unsigned}
    reflected = _reflect_pauli_monomial(mono, n)
    rep, shift = _translation_orbit_representative_shift(reflected, n)
    return ComplexF64(_momentum_phase(Float64, k, shift, n)), rep
end

function _reflection_row_label(terms::Vector{Tuple{ComplexF64,M}}) where {M<:NormalMonomial{PauliAlgebra}}
    return (
        feature=:reflection_adapted_row,
        terms=Any[(coefficient=coef, monomial=mono) for (coef, mono) in terms],
    )
end

function _rows_to_matrix(rows::Vector{Vector{C}}, ncols::Int) where {C<:Number}
    mat = Matrix{C}(undef, length(rows), ncols)
    for (row_idx, row) in pairs(rows)
        mat[row_idx, :] = row
    end
    return mat
end

function _translation_reflection_adapted_blocks(
    block_basis::Vector{M},
    k::Int,
    n::Int;
    atol::Float64,
) where {M<:NormalMonomial{PauliAlgebra}}
    _translation_reflection_fixed_momentum(k, n) || throw(ArgumentError(
        "Reflection adaptation is only supported for momentum sectors fixed by reflection; got k=$k, n=$n."
    ))

    basis_index = Dict{M,Int}(mono => idx for (idx, mono) in pairs(block_basis))
    visited = falses(length(block_basis))
    rows_by_parity = Dict(1 => Vector{ComplexF64}[], -1 => Vector{ComplexF64}[])
    sparse_rows_by_parity = Dict(
        1 => Vector{Vector{Tuple{Int,ComplexF64}}}(),
        -1 => Vector{Vector{Tuple{Int,ComplexF64}}}(),
    )
    labels_by_parity = Dict(1 => Any[], -1 => Any[])
    invsqrt2 = inv(sqrt(2.0))

    for i in eachindex(block_basis)
        visited[i] && continue
        mono = block_basis[i]
        phase, image = _translation_reflection_image(mono, k, n)
        j = get(basis_index, image, 0)
        j == 0 && throw(ArgumentError(
            "Reflection image $image of translation basis row $mono is missing from the momentum block basis."
        ))

        if i == j
            abs(imag(phase)) <= atol && abs(abs(real(phase)) - 1) <= atol || throw(ArgumentError(
                "Reflection fixed row $mono in momentum $k has non-real parity phase $phase."
            ))
            parity = real(phase) >= 0 ? 1 : -1
            row = zeros(ComplexF64, length(block_basis))
            row[i] = 1.0 + 0.0im
            sparse_row = Tuple{Int,ComplexF64}[(i, 1.0 + 0.0im)]
            push!(rows_by_parity[parity], row)
            push!(sparse_rows_by_parity[parity], sparse_row)
            push!(labels_by_parity[parity], _reflection_row_label([(1.0 + 0.0im, mono)]))
            visited[i] = true
        else
            phase_back, back = _translation_reflection_image(image, k, n)
            back == mono || throw(ArgumentError(
                "Reflection action is not involutive on momentum basis rows $mono and $image."
            ))
            abs(phase * phase_back - (1.0 + 0.0im)) <= atol || throw(ArgumentError(
                "Reflection phases for $mono and $image are inconsistent: $phase, $phase_back."
            ))
            visited[i] = true
            visited[j] = true
            for parity in (1, -1)
                row = zeros(ComplexF64, length(block_basis))
                row[i] = invsqrt2
                row[j] = parity * phase * invsqrt2
                sparse_row = Tuple{Int,ComplexF64}[
                    (i, ComplexF64(invsqrt2)),
                    (j, ComplexF64(parity * phase * invsqrt2)),
                ]
                push!(rows_by_parity[parity], row)
                push!(sparse_rows_by_parity[parity], sparse_row)
                push!(
                    labels_by_parity[parity],
                    _reflection_row_label([
                        (ComplexF64(invsqrt2), mono),
                        (ComplexF64(parity * phase * invsqrt2), image),
                    ]),
                )
            end
        end
    end

    blocks = NamedTuple[]
    for parity in (1, -1)
        rows = rows_by_parity[parity]
        isempty(rows) && continue
        push!(
            blocks,
            (
                reflection=parity,
                row_basis=_rows_to_matrix(rows, length(block_basis)),
                row_basis_sparse=_sparse_transform_rows(
                    sparse_rows_by_parity[parity],
                    length(block_basis),
                ),
                row_labels=labels_by_parity[parity],
            ),
        )
    end
    return blocks
end

function _real_reflection_row_label(
    parity::Int,
    terms::Vector{Tuple{ComplexF64,M}},
) where {M<:NormalMonomial{PauliAlgebra}}
    return (
        feature=:real_reflection_adapted_row,
        reflection=parity,
        terms=Any[(coefficient=coef, monomial=mono) for (coef, mono) in terms],
    )
end

function _realified_coordinate_row(
    terms::Vector{Tuple{ComplexF64,Int}},
    n_complex::Int,
)
    row = zeros(Float64, 2n_complex)
    for (coef, idx) in terms
        row[idx] += real(coef)
        row[n_complex + idx] += imag(coef)
    end
    return row
end

function _realified_coordinate_sparse_row(
    terms::Vector{Tuple{ComplexF64,Int}},
    n_complex::Int,
)
    sparse = Tuple{Int,Float64}[]
    for (coef, idx) in terms
        re = real(coef)
        im = imag(coef)
        iszero(re) || push!(sparse, (idx, re))
        iszero(im) || push!(sparse, (n_complex + idx, im))
    end
    return sparse
end

function _translation_real_reflection_adapted_blocks(
    block_basis::Vector{M},
    k::Int,
    n::Int;
    atol::Float64,
) where {M<:NormalMonomial{PauliAlgebra}}
    _translation_reflection_fixed_momentum(k, n) && throw(ArgumentError(
        "Realified reflection adaptation is for conjugate momentum pairs only; got fixed momentum k=$k, n=$n."
    ))

    basis_index = Dict{M,Int}(mono => idx for (idx, mono) in pairs(block_basis))
    visited = falses(length(block_basis))
    rows_by_parity = Dict(1 => Vector{Float64}[], -1 => Vector{Float64}[])
    sparse_rows_by_parity = Dict(
        1 => Vector{Vector{Tuple{Int,Float64}}}(),
        -1 => Vector{Vector{Tuple{Int,Float64}}}(),
    )
    labels_by_parity = Dict(1 => Any[], -1 => Any[])
    invsqrt2 = inv(sqrt(2.0))

    for i in eachindex(block_basis)
        visited[i] && continue
        mono = block_basis[i]
        phase, image = _translation_reflection_image(mono, k, n)
        j = get(basis_index, image, 0)
        j == 0 && throw(ArgumentError(
            "Reflection image $image of translation basis row $mono is missing from the momentum block basis."
        ))

        if i == j
            phase_back, back = _translation_reflection_image(image, k, n)
            back == mono || throw(ArgumentError(
                "Reflection action is not involutive on momentum basis row $mono."
            ))
            abs(phase_back * conj(phase) - (1.0 + 0.0im)) <= atol || throw(ArgumentError(
                "Antiunitary reflection phase for $mono is inconsistent: $phase, $phase_back."
            ))
            root_phase = sqrt(phase)
            for parity in (1, -1)
                coef = parity == 1 ? root_phase : im * root_phase
                push!(
                    rows_by_parity[parity],
                    _realified_coordinate_row([(coef, i)], length(block_basis)),
                )
                push!(
                    sparse_rows_by_parity[parity],
                    _realified_coordinate_sparse_row([(coef, i)], length(block_basis)),
                )
                push!(
                    labels_by_parity[parity],
                    _real_reflection_row_label(parity, [(coef, mono)]),
                )
            end
            visited[i] = true
        else
            phase_back, back = _translation_reflection_image(image, k, n)
            back == mono || throw(ArgumentError(
                "Reflection action is not involutive on momentum basis rows $mono and $image."
            ))
            abs(phase_back * conj(phase) - (1.0 + 0.0im)) <= atol || throw(ArgumentError(
                "Antiunitary reflection phases for $mono and $image are inconsistent: $phase, $phase_back."
            ))
            visited[i] = true
            visited[j] = true
            for parity in (1, -1)
                a_terms = [
                    (ComplexF64(invsqrt2), i),
                    (ComplexF64(parity * phase * invsqrt2), j),
                ]
                b_terms = [
                    (ComplexF64(im * invsqrt2), i),
                    (ComplexF64(-parity * im * phase * invsqrt2), j),
                ]
                push!(rows_by_parity[parity], _realified_coordinate_row(a_terms, length(block_basis)))
                push!(rows_by_parity[parity], _realified_coordinate_row(b_terms, length(block_basis)))
                push!(
                    sparse_rows_by_parity[parity],
                    _realified_coordinate_sparse_row(a_terms, length(block_basis)),
                )
                push!(
                    sparse_rows_by_parity[parity],
                    _realified_coordinate_sparse_row(b_terms, length(block_basis)),
                )
                push!(
                    labels_by_parity[parity],
                    _real_reflection_row_label(parity, [
                        (a_terms[1][1], mono),
                        (a_terms[2][1], image),
                    ]),
                )
                push!(
                    labels_by_parity[parity],
                    _real_reflection_row_label(parity, [
                        (b_terms[1][1], mono),
                        (b_terms[2][1], image),
                    ]),
                )
            end
        end
    end

    blocks = NamedTuple[]
    for parity in (1, -1)
        rows = rows_by_parity[parity]
        isempty(rows) && continue
        push!(
            blocks,
            (
                reflection=parity,
                row_basis=_rows_to_matrix(rows, 2 * length(block_basis)),
                row_basis_sparse=_sparse_transform_rows(
                    sparse_rows_by_parity[parity],
                    2 * length(block_basis),
                ),
                row_labels=labels_by_parity[parity],
            ),
        )
    end
    return blocks
end

@inline function _momentum_phase(::Type{R}, k::Int, r::Int, n::Int) where {R<:Real}
    (iszero(k) || iszero(r)) && return complex(one(R), zero(R))
    θ = -R(2) * R(pi) * R(k) * R(r) / R(n)
    return cis(θ)
end

mutable struct _TranslationProductCache{M,C,T<:Unsigned}
    terms::Dict{Tuple{M,M},Vector{Tuple{Int,C,M}}}
    hits::Int
    misses::Int
    buf::Vector{T}
end

function _TranslationProductCache(
    terms::Dict{Tuple{M,M},Vector{Tuple{Int,C,M}}},
    hits::Integer,
    misses::Integer,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    return _TranslationProductCache{M,C,T}(terms, Int(hits), Int(misses), T[])
end

function _translation_product_cache_hit_rate(cache::_TranslationProductCache)
    total = cache.hits + cache.misses
    return iszero(total) ? 0.0 : cache.hits / total
end

function _translation_momentum_entry(
    row::M,
    col::M,
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{P},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},P<:Polynomial{PauliAlgebra,T}}
    C = _coefficient_type(P)
    product_cache = Dict{Tuple{M,M},Vector{Tuple{Int,C,M}}}()
    return _translation_momentum_entry(
        row,
        col,
        k,
        n,
        translated,
        rep_cache,
        P,
        product_cache,
    )
end

function _translation_momentum_product_terms(
    row::M,
    col::M,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{C},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    return _translation_momentum_product_terms(
        row,
        col,
        n,
        translated,
        rep_cache,
        C,
        T[],
    )
end

function _translation_momentum_product_terms(
    row::M,
    col::M,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{C},
    buf::Vector{T},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    if isone(row) && isone(col)
        return Tuple{Int,C,M}[(0, one(C), row)]
    end

    R = typeof(real(one(C)))
    cross_scale = xor(isone(row), isone(col)) ? C(inv(sqrt(R(n)))) : one(C)
    one_mono = one(row)
    terms = Tuple{Int,C,M}[]
    sizehint!(terms, n)
    col_translates = translated[col]
    for r in 0:(n - 1)
        word, phase_k = simplify!(
            PauliAlgebra,
            _neat_dot3!(buf, one_mono, row, col_translates[r + 1]),
        )
        phase_k == 0x04 && continue
        coef = cross_scale * C(_coeff_to_number(PauliAlgebra, phase_k))
        if !iszero(coef)
            mono = _unchecked_monomial(PauliAlgebra, copy(word))
            rep = get!(rep_cache, mono) do
                _translation_orbit_representative(mono, n)
            end
            push!(terms, (r, coef, rep))
        end
    end
    return terms
end

function _translation_momentum_product_terms!(
    cache::Dict{Tuple{M,M},Vector{Tuple{Int,C,M}}},
    row::M,
    col::M,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{C},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    return get!(cache, (row, col)) do
        _translation_momentum_product_terms(row, col, n, translated, rep_cache, C)
    end
end

function _translation_momentum_product_terms!(
    cache::_TranslationProductCache{M,C,T},
    row::M,
    col::M,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{C},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    key = (row, col)
    cached = get(cache.terms, key, nothing)
    if cached !== nothing
        cache.hits += 1
        return cached
    end

    cache.misses += 1
    terms = _translation_momentum_product_terms(
        row,
        col,
        n,
        translated,
        rep_cache,
        C,
        cache.buf,
    )
    cache.terms[key] = terms
    return terms
end

function _qmbcertify_chain_momentum_product_terms(
    row::M,
    col::M,
    n::Int,
    translated::Dict{M,Vector{M}},
    ::Type{C},
    buf::Vector{T},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    if isone(row) && isone(col)
        return Tuple{Int,C,M}[(0, one(C), row)]
    end

    R = typeof(real(one(C)))
    cross_scale = xor(isone(row), isone(col)) ? C(inv(sqrt(R(n)))) : one(C)
    one_mono = one(row)
    terms = Tuple{Int,C,M}[]
    sizehint!(terms, n)
    col_translates = translated[col]
    for r in 0:(n - 1)
        word, phase_k = simplify!(
            PauliAlgebra,
            _neat_dot3!(buf, one_mono, row, col_translates[r + 1]),
        )
        mono, support_coef = _qmbcertify_chain_support_term(M, word, phase_k, n; realify=true)
        mono === nothing && continue
        coef = cross_scale * C(support_coef)
        iszero(coef) && continue
        push!(terms, (r, coef, mono))
    end
    return terms
end

function _qmbcertify_chain_momentum_product_terms!(
    cache::_TranslationProductCache{M,C,T},
    row::M,
    col::M,
    n::Int,
    translated::Dict{M,Vector{M}},
    ::Type{C},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    key = (row, col)
    cached = get(cache.terms, key, nothing)
    if cached !== nothing
        cache.hits += 1
        return cached
    end

    cache.misses += 1
    terms = _qmbcertify_chain_momentum_product_terms(
        row,
        col,
        n,
        translated,
        C,
        cache.buf,
    )
    cache.terms[key] = terms
    return terms
end

function _qmbcertify_realified_phase_number(phase_k::UInt8)
    phase_k == 0x04 && return 0.0
    coef = _coeff_to_number(PauliAlgebra, phase_k)
    return isreal(coef) ? Float64(real(coef)) : Float64(imag(coef))
end

function _qmbcertify_chain_momentum_block_complex_linear_entries(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    ::Type{K},
    ::Type{C},
    product_cache::_TranslationProductCache{M,C,T},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    m = length(block_basis)
    entries = Matrix{LinearMomentForm{K,C}}(undef, m, m)
    for j in 1:m, i in 1:j
        base_terms = _qmbcertify_chain_momentum_product_terms!(
            product_cache,
            block_basis[i],
            block_basis[j],
            n,
            translated,
            C,
        )
        if i == j
            entries[i, i] = _translation_momentum_diagonal_terms_to_linear_form(
                K,
                C,
                base_terms,
                k,
                n,
            )
        else
            entries[i, j] = _translation_momentum_terms_to_linear_form(
                K,
                C,
                base_terms,
                k,
                n,
            )
            entries[j, i] = _translation_momentum_adjoint_terms_to_linear_form(
                K,
                C,
                base_terms,
                k,
                n,
            )
        end
    end
    return entries
end

function _qmbcertify_chain_momentum_block_linear_entries(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    ::Type{K},
    ::Type{C},
    product_cache;
    real_moment_matrix::Bool=true,
    atol=1e-12,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    entry_coeff_type = C <: Real ? Complex{C} : C
    complex_entries = _qmbcertify_chain_momentum_block_complex_linear_entries(
        block_basis,
        k,
        n,
        translated,
        K,
        entry_coeff_type,
        product_cache,
    )
    real_moment_matrix || return complex_entries
    C <: Real || throw(ArgumentError("real_moment_matrix=true requires real linear coefficient type, got $C."))
    return _realify_hermitian_linear_block(complex_entries, C; atol)
end

function _qmbcertify_real_fixed_linear_block(
    entries::Matrix{LinearMomentForm{K,C}},
    ::Type{R};
    atol,
) where {K,C<:Number,R<:Real}
    real_entries = Matrix{LinearMomentForm{K,R}}(undef, size(entries))
    for idx in eachindex(entries)
        real_entries[idx] = _real_part_form_coefficients(entries[idx], R; atol)
    end
    return _symmetrize_real_linear_block(real_entries)
end

function _translation_momentum_entry(
    row::M,
    col::M,
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{P},
    product_cache,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    R = typeof(real(one(C)))
    base_terms = _translation_momentum_product_terms!(
        product_cache,
        row,
        col,
        n,
        translated,
        rep_cache,
        C,
    )

    terms = Tuple{C,M}[]
    sizehint!(terms, length(base_terms))
    for (r, coef, rep) in base_terms
        phase = C(_momentum_phase(R, k, r, n))
        push!(terms, (phase * coef, rep))
    end
    return convert(P, Polynomial(terms))
end

function _translation_momentum_entry_terms(
    row::M,
    col::M,
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{C},
    product_cache,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number}
    R = typeof(real(one(C)))
    base_terms = _translation_momentum_product_terms!(
        product_cache,
        row,
        col,
        n,
        translated,
        rep_cache,
        C,
    )

    terms = Tuple{C,M}[]
    sizehint!(terms, length(base_terms))
    for (r, coef, rep) in base_terms
        phase = C(_momentum_phase(R, k, r, n))
        value = phase * coef
        iszero(value) || push!(terms, (value, rep))
    end
    return terms
end

function _translation_terms_to_linear_form(
    ::Type{K},
    ::Type{C},
    terms::AbstractVector{<:Tuple{<:Number,M}},
) where {K,C,M<:NormalMonomial}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(terms))
    for (coef, mono) in terms
        converted = convert(C, coef)
        iszero(converted) && continue
        push!(pairs, _moment_key(K, mono) => converted)
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _translation_terms_to_linear_form(
    ::Type{K},
    ::Type{C},
    terms::AbstractVector{<:Tuple{<:Number,M}},
    quotient,
) where {K,C,M<:NormalMonomial}
    quotient === nothing && return _translation_terms_to_linear_form(K, C, terms)
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(terms))
    for (coef, mono) in terms
        converted = convert(C, coef)
        iszero(converted) && continue
        info = _axis_quotient_info(quotient, _moment_key(K, mono))
        info.forced_zero && continue
        value = converted * convert(C, info.sign)
        iszero(value) && continue
        push!(pairs, _owned_moment_key(K, info.key) => value)
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

@inline function _push_translation_linear_pair!(
    pairs::Vector{Pair{K,C}},
    ::Type{K},
    ::Type{C},
    mono::M,
    coef,
) where {K,C,M<:NormalMonomial}
    value = convert(C, coef)
    iszero(value) || push!(pairs, _moment_key(K, mono) => value)
    return nothing
end

@inline function _push_translation_quotient_linear_pair!(
    pairs::Vector{Pair{K,C}},
    ::Type{K},
    ::Type{C},
    mono::M,
    coef,
    quotient,
) where {K,C,M<:NormalMonomial}
    converted = convert(C, coef)
    iszero(converted) && return nothing
    info = _axis_quotient_info(quotient, _moment_key(K, mono))
    info.forced_zero && return nothing
    value = converted * convert(C, info.sign)
    iszero(value) || push!(pairs, _owned_moment_key(K, info.key) => value)
    return nothing
end

function _translation_momentum_terms_to_linear_form(
    ::Type{K},
    ::Type{C},
    base_terms::AbstractVector{<:Tuple{Int,<:Number,M}},
    k::Int,
    n::Int;
    quotient=nothing,
) where {K,C<:Number,M<:NormalMonomial}
    R = typeof(real(one(C)))
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(base_terms))
    if quotient === nothing
        for (r, coef, mono) in base_terms
            phase = C(_momentum_phase(R, k, r, n))
            _push_translation_linear_pair!(pairs, K, C, mono, phase * coef)
        end
    else
        for (r, coef, mono) in base_terms
            phase = C(_momentum_phase(R, k, r, n))
            _push_translation_quotient_linear_pair!(pairs, K, C, mono, phase * coef, quotient)
        end
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _translation_momentum_adjoint_terms_to_linear_form(
    ::Type{K},
    ::Type{C},
    base_terms::AbstractVector{<:Tuple{Int,<:Number,M}},
    k::Int,
    n::Int;
    quotient=nothing,
) where {K,C<:Number,M<:NormalMonomial}
    R = typeof(real(one(C)))
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(base_terms))
    if quotient === nothing
        for (r, coef, mono) in base_terms
            phase = C(_momentum_phase(R, k, r, n))
            value = conj(phase * coef)
            _push_translation_linear_pair!(
                pairs,
                K,
                C,
                _moment_linear_adjoint_monomial(mono),
                value,
            )
        end
    else
        for (r, coef, mono) in base_terms
            phase = C(_momentum_phase(R, k, r, n))
            value = conj(phase * coef)
            _push_translation_quotient_linear_pair!(
                pairs,
                K,
                C,
                _moment_linear_adjoint_monomial(mono),
                value,
                quotient,
            )
        end
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _translation_momentum_diagonal_terms_to_linear_form(
    ::Type{K},
    ::Type{C},
    base_terms::AbstractVector{<:Tuple{Int,<:Number,M}},
    k::Int,
    n::Int;
    quotient=nothing,
) where {K,C<:Number,M<:NormalMonomial}
    R = typeof(real(one(C)))
    half = convert(C, 0.5)
    pairs = Pair{K,C}[]
    sizehint!(pairs, 2 * length(base_terms))
    if quotient === nothing
        for (r, coef, mono) in base_terms
            phase = C(_momentum_phase(R, k, r, n))
            value = phase * coef
            _push_translation_linear_pair!(pairs, K, C, mono, half * value)
            _push_translation_linear_pair!(
                pairs,
                K,
                C,
                _moment_linear_adjoint_monomial(mono),
                half * conj(value),
            )
        end
    else
        for (r, coef, mono) in base_terms
            phase = C(_momentum_phase(R, k, r, n))
            value = phase * coef
            _push_translation_quotient_linear_pair!(pairs, K, C, mono, half * value, quotient)
            _push_translation_quotient_linear_pair!(
                pairs,
                K,
                C,
                _moment_linear_adjoint_monomial(mono),
                half * conj(value),
                quotient,
            )
        end
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _translation_adjoint_entry_terms(
    ::Type{C},
    terms::AbstractVector{<:Tuple{<:Number,M}},
) where {C,M<:NormalMonomial}
    adjoint_terms = Tuple{C,M}[]
    sizehint!(adjoint_terms, length(terms))
    for (coef, mono) in terms
        converted = convert(C, conj(coef))
        iszero(converted) && continue
        push!(adjoint_terms, (converted, _moment_linear_adjoint_monomial(mono)))
    end
    return adjoint_terms
end

function _translation_hermitian_diagonal_terms(
    ::Type{C},
    terms::AbstractVector{<:Tuple{<:Number,M}},
) where {C,M<:NormalMonomial}
    half = convert(C, 0.5)
    diag_terms = Tuple{C,M}[]
    sizehint!(diag_terms, 2 * length(terms))
    for (coef, mono) in terms
        converted = convert(C, coef)
        iszero(converted) || push!(diag_terms, (half * converted, mono))
        adjoint_converted = convert(C, conj(coef))
        iszero(adjoint_converted) ||
            push!(diag_terms, (half * adjoint_converted, _moment_linear_adjoint_monomial(mono)))
    end
    return diag_terms
end

function _translation_momentum_block_complex_linear_entries(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{K},
    ::Type{C},
    product_cache;
    quotient=nothing,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    m = length(block_basis)
    entries = Matrix{LinearMomentForm{K,C}}(undef, m, m)
    for j in 1:m, i in 1:j
        base_terms = _translation_momentum_product_terms!(
            product_cache,
            block_basis[i],
            block_basis[j],
            n,
            translated,
            rep_cache,
            C,
        )
        if i == j
            entries[i, i] = _translation_momentum_diagonal_terms_to_linear_form(
                K,
                C,
                base_terms,
                k,
                n;
                quotient,
            )
        else
            entries[i, j] = _translation_momentum_terms_to_linear_form(
                K,
                C,
                base_terms,
                k,
                n;
                quotient,
            )
            entries[j, i] = _translation_momentum_adjoint_terms_to_linear_form(
                K,
                C,
                base_terms,
                k,
                n;
                quotient,
            )
        end
    end
    return entries
end

function _real_part_form_coefficients(
    form::LinearMomentForm{K,C},
    ::Type{R};
    atol,
) where {K,C<:Number,R<:Real}
    pairs = Pair{K,R}[]
    sizehint!(pairs, length(form))
    tolerance = R(atol)
    for (key, coef) in form
        value = _clean_real_part(real(coef), R, tolerance)
        iszero(value) || push!(pairs, key => value)
    end
    return LinearMomentForm{K,R}(pairs, Val(:trusted))
end

function _imag_part_form_coefficients(
    form::LinearMomentForm{K,C},
    ::Type{R};
    atol,
) where {K,C<:Number,R<:Real}
    pairs = Pair{K,R}[]
    sizehint!(pairs, length(form))
    tolerance = R(atol)
    for (key, coef) in form
        value = _clean_real_part(imag(coef), R, tolerance)
        iszero(value) || push!(pairs, key => value)
    end
    return LinearMomentForm{K,R}(pairs, Val(:trusted))
end

function _scale_linear_form(
    form::LinearMomentForm{K,C},
    scale,
) where {K,C}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(form))
    converted_scale = convert(C, scale)
    for (key, coef) in form
        value = converted_scale * coef
        iszero(value) || push!(pairs, key => value)
    end
    return LinearMomentForm{K,C}(pairs, Val(:trusted))
end

@inline _within_abs_atol(value::Real, atol) = abs(value) <= atol
@inline _within_abs_atol(value::Complex, atol) = abs2(value) <= atol * atol

function _add_scaled_linear_form_terms!(
    pairs::Vector{Pair{K,C}},
    factor,
    form::LinearMomentForm{K},
    ::Type{C};
    atol,
) where {K,C}
    _within_abs_atol(factor, atol) && return nothing
    converted_factor = convert(C, factor)
    for (key, coef) in form
        value = converted_factor * coef
        _within_abs_atol(value, atol) && continue
        push!(pairs, key => value)
    end
    return nothing
end

function _add_scaled_linear_form_terms!(
    acc::Dict{K,C},
    factor,
    form::LinearMomentForm{K},
    ::Type{C};
    atol,
) where {K,C}
    _within_abs_atol(factor, atol) && return nothing
    converted_factor = convert(C, factor)
    for (key, coef) in form
        value = converted_factor * coef
        _within_abs_atol(value, atol) && continue
        updated = get(acc, key, zero(C)) + value
        if _within_abs_atol(updated, atol)
            delete!(acc, key)
        else
            acc[key] = updated
        end
    end
    return nothing
end

function _linear_moment_form_from_unique_owned_pairs!(terms::Vector{Pair{K,C}}) where {K,C}
    if isempty(terms)
        return LinearMomentForm{K,C}(terms, Val(:trusted))
    end

    write_idx = 0
    for term in terms
        iszero(term.second) && continue
        write_idx += 1
        terms[write_idx] = term
    end
    resize!(terms, write_idx)
    write_idx == 0 && return LinearMomentForm{K,C}(terms, Val(:trusted))
    write_idx == 1 && return LinearMomentForm{K,C}(terms, Val(:trusted))

    sort!(terms; by=x -> x.first, lt=key_lt)
    @inbounds for idx in 2:write_idx
        if key_isequal(terms[idx - 1].first, terms[idx].first)
            return _linear_moment_form_from_owned_pairs!(terms)
        end
    end
    return LinearMomentForm{K,C}(terms, Val(:trusted))
end

function _linear_moment_form_from_distinct_owned_pairs!(terms::Vector{Pair{K,C}}) where {K,C}
    if isempty(terms)
        return LinearMomentForm{K,C}(terms, Val(:trusted))
    end

    write_idx = 0
    for term in terms
        iszero(term.second) && continue
        write_idx += 1
        terms[write_idx] = term
    end
    resize!(terms, write_idx)
    write_idx <= 1 && return LinearMomentForm{K,C}(terms, Val(:trusted))

    sort!(terms; by=x -> x.first, lt=key_lt)
    return LinearMomentForm{K,C}(terms, Val(:trusted))
end

function _linear_moment_form_from_accumulator!(
    acc::Dict{K,C};
    atol,
) where {K,C}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(acc))
    for (key, coef) in acc
        abs(coef) <= atol && continue
        push!(pairs, key => coef)
    end
    return _linear_moment_form_from_distinct_owned_pairs!(pairs)
end

function _transform_linear_block(
    entries::Matrix{LinearMomentForm{K,C}},
    U_left::AbstractMatrix{<:Number},
    U_right::AbstractMatrix{<:Number};
    atol::Real=1e-12,
) where {K,C}
    source_size = size(entries)
    size(U_left, 2) == source_size[1] || throw(DimensionMismatch(
        "left transform has $(size(U_left, 2)) columns but linear block has $(source_size[1]) rows"
    ))
    size(U_right, 2) == source_size[2] || throw(DimensionMismatch(
        "right transform has $(size(U_right, 2)) columns but linear block has $(source_size[2]) columns"
    ))

    transformed = Matrix{LinearMomentForm{K,C}}(undef, size(U_left, 1), size(U_right, 1))
    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(atol) : atol
    entry_lengths = map(length, entries)
    for j in axes(transformed, 2), i in axes(transformed, 1)
        pair_count_hint = 0
        for b in axes(entries, 2)
            right_coeff = conj(U_right[j, b])
            abs(right_coeff) <= tolerance && continue
            for a in axes(entries, 1)
                left_coeff = U_left[i, a]
                abs(left_coeff) <= tolerance && continue
                pair_count_hint += entry_lengths[a, b]
            end
        end
        transformed[i, j] = if pair_count_hint >= 128
            acc = Dict{K,C}()
            sizehint!(acc, min(pair_count_hint >>> 1, 8192))
            for b in axes(entries, 2)
                right_coeff = conj(U_right[j, b])
                abs(right_coeff) <= tolerance && continue
                for a in axes(entries, 1)
                    left_coeff = U_left[i, a]
                    abs(left_coeff) <= tolerance && continue
                    entry_form = entries[a, b]
                    isempty(entry_form) && continue
                    coeff = left_coeff * right_coeff
                    _add_scaled_linear_form_terms!(
                        acc,
                        coeff,
                        entry_form,
                        C;
                        atol=tolerance,
                    )
                end
            end
            _linear_moment_form_from_accumulator!(acc; atol=tolerance)
        else
            pairs = Pair{K,C}[]
            sizehint!(pairs, pair_count_hint)
            for b in axes(entries, 2)
                right_coeff = conj(U_right[j, b])
                abs(right_coeff) <= tolerance && continue
                for a in axes(entries, 1)
                    left_coeff = U_left[i, a]
                    abs(left_coeff) <= tolerance && continue
                    entry_form = entries[a, b]
                    isempty(entry_form) && continue
                    coeff = left_coeff * right_coeff
                    _add_scaled_linear_form_terms!(
                        pairs,
                        coeff,
                        entry_form,
                        C;
                        atol=tolerance,
                    )
                end
            end
            _linear_moment_form_from_owned_pairs!(pairs)
        end
    end
    return transformed
end

function _transform_linear_block(
    entries::Matrix{LinearMomentForm{K,C}},
    U_left::_SparseTransformRows,
    U_right::_SparseTransformRows;
    atol::Real=1e-12,
) where {K,C}
    source_size = size(entries)
    size(U_left, 2) == source_size[1] || throw(DimensionMismatch(
        "left transform has $(size(U_left, 2)) columns but linear block has $(source_size[1]) rows"
    ))
    size(U_right, 2) == source_size[2] || throw(DimensionMismatch(
        "right transform has $(size(U_right, 2)) columns but linear block has $(source_size[2]) columns"
    ))

    transformed = Matrix{LinearMomentForm{K,C}}(undef, size(U_left, 1), size(U_right, 1))
    tolerance = real(one(C)) isa Real ? typeof(real(one(C)))(atol) : atol
    entry_lengths = map(length, entries)
    for j in axes(transformed, 2), i in axes(transformed, 1)
        pair_count_hint = 0
        for (a, _) in U_left.rows[i]
            for (b, _) in U_right.rows[j]
                pair_count_hint += entry_lengths[a, b]
            end
        end
        pairs = Pair{K,C}[]
        sizehint!(pairs, pair_count_hint)
        for (a, left_coeff) in U_left.rows[i]
            for (b, right_coeff) in U_right.rows[j]
                entry_form = entries[a, b]
                isempty(entry_form) && continue
                coeff = left_coeff * conj(right_coeff)
                _add_scaled_linear_form_terms!(
                    pairs,
                    coeff,
                    entry_form,
                    C;
                    atol=tolerance,
                )
            end
        end
        transformed[i, j] = _linear_moment_form_from_owned_pairs!(pairs)
    end
    return transformed
end

function _transform_real_symmetric_linear_block(
    entries::Matrix{LinearMomentForm{K,R}},
    U::_SparseTransformRows{TC};
    atol::Real=1e-12,
) where {K,R<:Real,TC<:Real}
    source_size = size(entries)
    source_size[1] == source_size[2] || throw(DimensionMismatch(
        "symmetric linear block must be square, got $source_size"
    ))
    size(U, 2) == source_size[1] || throw(DimensionMismatch(
        "transform has $(size(U, 2)) columns but linear block has $(source_size[1]) rows"
    ))

    transformed = Matrix{LinearMomentForm{K,R}}(undef, size(U, 1), size(U, 1))
    tolerance = R(atol)
    entry_lengths = map(length, entries)
    for j in axes(transformed, 2), i in 1:j
        pair_count_hint = 0
        for (a, _) in U.rows[i]
            for (b, _) in U.rows[j]
                pair_count_hint += entry_lengths[a, b]
            end
        end
        form = if pair_count_hint >= 128
            acc = Dict{K,R}()
            sizehint!(acc, min(pair_count_hint >>> 1, 8192))
            for (a, left_coeff) in U.rows[i]
                for (b, right_coeff) in U.rows[j]
                    entry_form = entries[a, b]
                    isempty(entry_form) && continue
                    coeff = left_coeff * right_coeff
                    _add_scaled_linear_form_terms!(
                        acc,
                        coeff,
                        entry_form,
                        R;
                        atol=tolerance,
                    )
                end
            end
            _linear_moment_form_from_accumulator!(acc; atol=tolerance)
        else
            pairs = Pair{K,R}[]
            sizehint!(pairs, pair_count_hint)
            for (a, left_coeff) in U.rows[i]
                for (b, right_coeff) in U.rows[j]
                    entry_form = entries[a, b]
                    isempty(entry_form) && continue
                    coeff = left_coeff * right_coeff
                    _add_scaled_linear_form_terms!(
                        pairs,
                        coeff,
                        entry_form,
                        R;
                        atol=tolerance,
                    )
                end
            end
            _linear_moment_form_from_owned_pairs!(pairs)
        end
        transformed[i, j] = form
        transformed[j, i] = form
    end
    return transformed
end

function _conj_pauli_linear_form(form::LinearMomentForm{K,C}) where {K,C}
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(form))
    for (key, coef) in form
        converted = convert(C, conj(coef))
        iszero(converted) || push!(pairs, key => converted)
    end
    return LinearMomentForm{K,C}(pairs, Val(:trusted))
end

function _average_pauli_linear_forms(
    left::LinearMomentForm{K,C},
    right::LinearMomentForm{K,C},
) where {K,C}
    half = convert(C, 0.5)
    pairs = Pair{K,C}[]
    sizehint!(pairs, length(left) + length(right))
    for (key, coef) in left
        push!(pairs, key => half * coef)
    end
    for (key, coef) in right
        push!(pairs, key => half * coef)
    end
    return _linear_moment_form_from_owned_pairs!(pairs)
end

function _symmetrize_real_linear_block(
    entries::Matrix{LinearMomentForm{K,R}},
) where {K,R<:Real}
    n = size(entries, 1)
    size(entries, 2) == n || throw(DimensionMismatch("Real linear PSD block must be square, got $(size(entries))."))
    out = Matrix{LinearMomentForm{K,R}}(undef, n, n)
    for i in 1:n
        out[i, i] = entries[i, i]
        for j in (i + 1):n
            avg = _average_pauli_linear_forms(entries[i, j], entries[j, i])
            out[i, j] = avg
            out[j, i] = avg
        end
    end
    return out
end

function _hermitianize_pauli_linear_block(
    entries::Matrix{LinearMomentForm{K,C}},
) where {K,C}
    n = size(entries, 1)
    size(entries, 2) == n || throw(DimensionMismatch("Hermitian linear block must be square, got $(size(entries))."))
    out = Matrix{LinearMomentForm{K,C}}(undef, n, n)
    for i in 1:n
        out[i, i] = _average_pauli_linear_forms(
            entries[i, i],
            _conj_pauli_linear_form(entries[i, i]),
        )
        for j in (i + 1):n
            avg = _average_pauli_linear_forms(
                entries[i, j],
                _conj_pauli_linear_form(entries[j, i]),
            )
            out[i, j] = avg
            out[j, i] = _conj_pauli_linear_form(avg)
        end
    end
    return out
end

function _realify_hermitian_linear_block(
    entries::Matrix{LinearMomentForm{K,C}},
    ::Type{R};
    atol,
) where {K,C<:Number,R<:Real}
    n = size(entries, 1)
    size(entries, 2) == n || throw(DimensionMismatch("Hermitian linear block must be square, got $(size(entries))."))

    re = Matrix{LinearMomentForm{K,R}}(undef, n, n)
    im = Matrix{LinearMomentForm{K,R}}(undef, n, n)
    any_imag = false
    for j in 1:n, i in 1:n
        im_entry = _imag_part_form_coefficients(entries[i, j], R; atol)
        re[i, j] = _real_part_form_coefficients(entries[i, j], R; atol)
        im[i, j] = im_entry
        any_imag |= !isempty(im_entry)
    end
    !any_imag && return re

    realified = Matrix{LinearMomentForm{K,R}}(undef, 2n, 2n)
    for j in 1:n, i in 1:n
        realified[i, j] = re[i, j]
        realified[i, n + j] = _scale_linear_form(im[i, j], -one(R))
        realified[n + i, j] = im[i, j]
        realified[n + i, n + j] = re[i, j]
    end
    return realified
end

function _realify_hermitian_linear_block_full(
    entries::Matrix{LinearMomentForm{K,C}},
    ::Type{R};
    atol,
) where {K,C<:Number,R<:Real}
    n = size(entries, 1)
    size(entries, 2) == n || throw(DimensionMismatch("Hermitian linear block must be square, got $(size(entries))."))

    re = Matrix{LinearMomentForm{K,R}}(undef, n, n)
    im = Matrix{LinearMomentForm{K,R}}(undef, n, n)
    for j in 1:n, i in 1:n
        re[i, j] = _real_part_form_coefficients(entries[i, j], R; atol)
        im[i, j] = _imag_part_form_coefficients(entries[i, j], R; atol)
    end

    realified = Matrix{LinearMomentForm{K,R}}(undef, 2n, 2n)
    for j in 1:n, i in 1:n
        realified[i, j] = re[i, j]
        realified[i, n + j] = _scale_linear_form(im[i, j], -one(R))
        realified[n + i, j] = im[i, j]
        realified[n + i, n + j] = re[i, j]
    end
    return realified
end

function _translation_momentum_block_linear_entries(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{K},
    ::Type{C},
    product_cache;
    real_moment_matrix::Bool=true,
    atol=1e-12,
    quotient=nothing,
    full_realification::Bool=false,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},K,C<:Number}
    entry_coeff_type = C <: Real ? Complex{C} : C
    complex_entries = _translation_momentum_block_complex_linear_entries(
        block_basis,
        k,
        n,
        translated,
        rep_cache,
        K,
        entry_coeff_type,
        product_cache,
        quotient=quotient,
    )
    real_moment_matrix || return complex_entries
    C <: Real || throw(ArgumentError("real_moment_matrix=true requires real linear coefficient type, got $C."))
    full_realification && return _realify_hermitian_linear_block_full(complex_entries, C; atol)
    return _realify_hermitian_linear_block(complex_entries, C; atol)
end

@inline function _clean_real_part(x, ::Type{R}, atol::R) where {R<:Real}
    y = R(x)
    return abs(y) <= atol ? zero(R) : y
end

function _real_part_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    ::Type{P};
    atol,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T}}
    R = _coefficient_type(P)
    terms = Tuple{R,NormalMonomial{PauliAlgebra,T}}[]
    for (coef, mono) in poly.terms
        value = _clean_real_part(real(coef), R, atol)
        iszero(value) || push!(terms, (value, mono))
    end
    return P(terms)
end

function _imag_part_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    ::Type{P};
    atol,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T}}
    R = _coefficient_type(P)
    terms = Tuple{R,NormalMonomial{PauliAlgebra,T}}[]
    for (coef, mono) in poly.terms
        value = _clean_real_part(imag(coef), R, atol)
        iszero(value) || push!(terms, (value, mono))
    end
    return P(terms)
end

function _conj_pauli_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    ::Type{P},
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T}}
    CP = _coefficient_type(P)
    terms = Tuple{CP,NormalMonomial{PauliAlgebra,T}}[]
    sizehint!(terms, length(poly.terms))
    for (coef, mono) in poly.terms
        push!(terms, (convert(CP, conj(coef)), mono))
    end
    return P(terms)
end

function _hermitianize_pauli_polynomial_block(
    mat::Matrix{P},
    ::Type{P},
) where {T<:Unsigned,P<:Polynomial{PauliAlgebra,T}}
    n = size(mat, 1)
    size(mat, 2) == n || throw(DimensionMismatch("Hermitian block must be square, got $(size(mat))."))
    out = Matrix{P}(undef, n, n)
    for i in 1:n
        out[i, i] = convert(P, (mat[i, i] + _conj_pauli_polynomial(mat[i, i], P)) * (1 // 2))
        for j in (i + 1):n
            avg = convert(P, (mat[i, j] + _conj_pauli_polynomial(mat[j, i], P)) * (1 // 2))
            out[i, j] = avg
            out[j, i] = _conj_pauli_polynomial(avg, P)
        end
    end
    return out
end

function _symmetrize_real_polynomial_block(
    mat::Matrix{P},
    ::Type{P},
    ;
    atol=nothing,
) where {T<:Unsigned,P<:Polynomial{PauliAlgebra,T}}
    n = size(mat, 1)
    size(mat, 2) == n || throw(DimensionMismatch("Real polynomial PSD block must be square, got $(size(mat))."))
    if atol !== nothing
        threshold = 100 * Float64(atol)
        for j in 1:n, i in 1:(j - 1)
            residual = _max_abs_polynomial_coefficient(mat[i, j] - mat[j, i])
            residual <= threshold || throw(ArgumentError(
                "Real polynomial PSD block is not symmetric within tolerance; " *
                "entries ($i, $j) and ($j, $i) differ by coefficient residual $residual."
            ))
        end
    end
    out = Matrix{P}(undef, n, n)
    for i in 1:n
        out[i, i] = mat[i, i]
        for j in (i + 1):n
            avg = convert(P, (mat[i, j] + mat[j, i]) * (1 // 2))
            out[i, j] = avg
            out[j, i] = avg
        end
    end
    return out
end

function _max_abs_polynomial_coefficient(poly::Polynomial)
    isempty(poly.terms) && return 0.0
    return Float64(maximum(abs(coef) for (coef, _) in poly.terms))
end

function _realify_hermitian_block(
    mat::Matrix{PComplex},
    ::Type{PReal};
    atol,
) where {T<:Unsigned,PComplex<:Polynomial{PauliAlgebra,T},PReal<:Polynomial{PauliAlgebra,T}}
    n = size(mat, 1)
    size(mat, 2) == n || throw(DimensionMismatch("Hermitian block must be square, got $(size(mat))."))

    re = Matrix{PReal}(undef, n, n)
    im = Matrix{PReal}(undef, n, n)
    any_imag = false
    for j in 1:n, i in 1:n
        im_entry = _imag_part_polynomial(mat[i, j], PReal; atol)
        re[i, j] = _real_part_polynomial(mat[i, j], PReal; atol)
        im[i, j] = im_entry
        any_imag |= !iszero(im_entry)
    end
    !any_imag && return re

    realified = Matrix{PReal}(undef, 2n, 2n)
    for j in 1:n, i in 1:n
        realified[i, j] = re[i, j]
        realified[i, n + j] = -im[i, j]
        realified[n + i, j] = im[i, j]
        realified[n + i, n + j] = re[i, j]
    end
    return realified
end

function _realify_hermitian_block_full(
    mat::Matrix{PComplex},
    ::Type{PReal};
    atol,
) where {T<:Unsigned,PComplex<:Polynomial{PauliAlgebra,T},PReal<:Polynomial{PauliAlgebra,T}}
    n = size(mat, 1)
    size(mat, 2) == n || throw(DimensionMismatch("Hermitian block must be square, got $(size(mat))."))

    re = Matrix{PReal}(undef, n, n)
    im = Matrix{PReal}(undef, n, n)
    for j in 1:n, i in 1:n
        re[i, j] = _real_part_polynomial(mat[i, j], PReal; atol)
        im[i, j] = _imag_part_polynomial(mat[i, j], PReal; atol)
    end

    realified = Matrix{PReal}(undef, 2n, 2n)
    for j in 1:n, i in 1:n
        realified[i, j] = re[i, j]
        realified[i, n + j] = -im[i, j]
        realified[n + i, j] = im[i, j]
        realified[n + i, n + j] = re[i, j]
    end
    return realified
end

function _translation_momentum_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{P},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},P<:Polynomial{PauliAlgebra,T}}
    C = _coefficient_type(P)
    product_cache = Dict{Tuple{M,M},Vector{Tuple{Int,C,M}}}()
    return _translation_momentum_block(
        block_basis,
        k,
        n,
        translated,
        rep_cache,
        P,
        product_cache,
    )
end

function _translation_momentum_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    rep_cache::Dict{M,M},
    ::Type{P},
    product_cache,
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    m = length(block_basis)
    mat = Matrix{P}(undef, m, m)
    for j in 1:m, i in 1:j
        entry = _translation_momentum_entry(
            block_basis[i],
            block_basis[j],
            k,
            n,
            translated,
            rep_cache,
            P,
            product_cache,
        )
        if i == j
            mat[i, i] = (entry + adjoint(entry)) / 2
        else
            mat[i, j] = entry
            mat[j, i] = adjoint(entry)
        end
    end
    return mat
end

_coefficient_type(::Type{<:Polynomial{A,T,C}}) where {A,T,C} = C
