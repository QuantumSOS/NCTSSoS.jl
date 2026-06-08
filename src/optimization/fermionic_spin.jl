"""
    FermionicSpinBlockLabel

Label attached to one Casimir-adapted SU(2) block produced by
[`FermionicSpinAdaptationSpec`](@ref). Carries:

- `sector::FermionicSectorLabel` — the underlying fermionic-sector label
  the SU(2) split was layered on top of.
- `total_spin2::Int` — twice the total spin quantum number `2 S` (so
  `0` = singlet `S = 0`, `1` = doublet `S = 1/2`, `2` = triplet
  `S = 1`, ...). The associated Casimir eigenvalue is
  `S(S+1) = total_spin2 * (total_spin2 + 2) / 4`.
- `multiplicity_tag::Int` — opaque tag distinguishing equivalent SU(2)
  blocks of the same total spin within the same sector.
"""
struct FermionicSpinBlockLabel
    sector::FermionicSectorLabel
    total_spin2::Int
    multiplicity_tag::Int
end

function Base.show(io::IO, label::FermionicSpinBlockLabel)
    print(
        io,
        "FermionicSpinBlockLabel(sector=$(label.sector), total_spin2=$(label.total_spin2), multiplicity_tag=$(label.multiplicity_tag))",
    )
end

function _validate_matching_mode_layouts(
    sector_layout::FermionicModeLayout,
    spin_layout::FermionicModeLayout,
    active_modes::AbstractVector{<:Integer},
)
    for mode_raw in active_modes
        mode = Int(mode_raw)
        _fermionic_mode_orbital(sector_layout, mode; context="Fermionic spin adaptation") ==
            _fermionic_mode_orbital(spin_layout, mode; context="Fermionic spin adaptation") || throw(ArgumentError(
                "Fermionic spin adaptation needs `mode_layout.orbital_of[$mode]` to match the fermionic sector metadata."
            ))

        sector_spin = get(sector_layout.spin2_of, mode, nothing)
        spin_spin = get(spin_layout.spin2_of, mode, nothing)
        sector_spin == spin_spin || throw(ArgumentError(
            "Fermionic spin adaptation needs `mode_layout.spin2_of[$mode]` to match the fermionic sector metadata."
        ))
    end

    return nothing
end

function _validate_fermionic_spin_generators(
    generators::AbstractVector{<:FermionicModePermutation},
    mode_layout::FermionicModeLayout,
    active_modes::AbstractVector{<:Integer},
)
    for g in generators, mode_raw in active_modes
        mode = Int(mode_raw)
        dst = get(g.images, mode, mode)
        get(mode_layout.spin2_of, dst, nothing) == get(mode_layout.spin2_of, mode, nothing) || throw(ArgumentError(
            "Fermionic spin adaptation currently requires every `FermionicModePermutation` to preserve spin labels on active modes."
        ))
    end
    return nothing
end

function _fermionic_spin_orbital_pairs(
    mode_layout::FermionicModeLayout,
    active_modes::AbstractVector{<:Integer},
)
    by_orbital = Dict{Int,Dict{Int,Int}}()

    for mode_raw in active_modes
        mode = Int(mode_raw)
        orbital = _fermionic_mode_orbital(mode_layout, mode; context="Fermionic spin adaptation")
        spin2 = get(mode_layout.spin2_of, mode, nothing)
        isnothing(spin2) && throw(ArgumentError(
            "Fermionic spin adaptation needs `spin2_of[$mode]` for every active mode."
        ))
        spin2 in (-1, 1) || throw(ArgumentError(
            "Fermionic spin adaptation currently supports only spin-1/2 doublets encoded by `spin2_of[mode] ∈ {-1, +1}`. Got $spin2 for mode $mode."
        ))

        bucket = get!(by_orbital, orbital, Dict{Int,Int}())
        haskey(bucket, spin2) && bucket[spin2] != mode && throw(ArgumentError(
            "Fermionic spin adaptation needs exactly one active mode per `(orbital, spin)` pair. Orbital $orbital has duplicate spin label $spin2."
        ))
        bucket[spin2] = mode
    end

    pairs = Pair{Int,Tuple{Int,Int}}[]
    for orbital in sort!(collect(keys(by_orbital)))
        bucket = by_orbital[orbital]
        haskey(bucket, 1) && haskey(bucket, -1) || throw(ArgumentError(
            "Fermionic spin adaptation needs a complete `(↑, ↓)` doublet for every active orbital. Orbital $orbital is incomplete."
        ))
        push!(pairs, orbital => (bucket[1], bucket[-1]))
    end

    return pairs
end

function _validate_fermionic_spin_adaptation_spec(
    symmetry::SymmetrySpec,
    active_modes::AbstractVector{<:Integer},
)
    spin_adaptation = symmetry.spin_adaptation
    isnothing(spin_adaptation) && return nothing

    sector = symmetry.sector
    isnothing(sector) && throw(ArgumentError(
        "Fermionic spin adaptation needs a fermionic sector split; supply `sector=FermionicSectorSpec(..., split_spin=true)` as well."
    ))
    sector.split_spin || throw(ArgumentError(
        "Fermionic spin adaptation needs `split_spin=true` in the supplied `FermionicSectorSpec`."
    ))
    spin_adaptation.casimir_tol > 0 || throw(ArgumentError(
        "`FermionicSpinAdaptationSpec.casimir_tol` must be positive."
    ))
    spin_adaptation.compress_multiplets && throw(ArgumentError(
        "`compress_multiplets=true` is not implemented yet; keep multiplicity explicit for now."
    ))

    _validate_matching_mode_layouts(sector.mode_layout, spin_adaptation.mode_layout, active_modes)
    _validate_fermionic_spin_generators(
        symmetry.fermionic_generators,
        spin_adaptation.mode_layout,
        active_modes,
    )
    _fermionic_spin_orbital_pairs(spin_adaptation.mode_layout, active_modes)

    return nothing
end

@inline _fermionic_mode_operator(::Type{T}, idx::Integer) where {T<:Integer} =
    NormalMonomial{FermionicAlgebra,T}(T[idx])

function _fermionic_total_spin_operators(
    ::Type{T},
    spin_adaptation::FermionicSpinAdaptationSpec,
    active_modes::AbstractVector{<:Integer},
) where {T<:Integer}
    pairs = _fermionic_spin_orbital_pairs(spin_adaptation.mode_layout, active_modes)
    isempty(pairs) && throw(ArgumentError("Fermionic spin adaptation needs at least one active spin doublet."))

    up_probe, dn_probe = last(first(pairs))
    probe_poly = _fermionic_mode_operator(T, -up_probe) * _fermionic_mode_operator(T, up_probe)
    sz = zero(probe_poly)
    splus = zero(probe_poly)
    sminus = zero(probe_poly)

    for (_, (up_mode, dn_mode)) in pairs
        up = _fermionic_mode_operator(T, up_mode)
        dn = _fermionic_mode_operator(T, dn_mode)
        up_dag = _fermionic_mode_operator(T, -up_mode)
        dn_dag = _fermionic_mode_operator(T, -dn_mode)

        sz += 0.5 * (up_dag * up - dn_dag * dn)
        splus += up_dag * dn
        sminus += dn_dag * up
    end

    return sz, splus, sminus
end

@inline _adjoint_action(op::Polynomial{FermionicAlgebra,T,C}, mono::NormalMonomial{FermionicAlgebra,T}) where {T<:Integer,C<:Number} =
    op * mono - mono * op
@inline _adjoint_action(op::Polynomial{FermionicAlgebra,T,C1}, poly::Polynomial{FermionicAlgebra,T,C2}) where {T<:Integer,C1<:Number,C2<:Number} =
    op * poly - poly * op

function _fermionic_polynomial_coordinate_vector(
    poly::Polynomial{FermionicAlgebra,T,C},
    basis_index::Dict{NormalMonomial{FermionicAlgebra,T},Int},
    context::AbstractString;
    atol::Float64=_SYMMETRY_ATOL,
) where {T<:Integer,C<:Number}
    coords = zeros(ComplexF64, length(basis_index))
    for (coef, mono) in poly.terms
        idx = get(basis_index, mono, 0)
        if idx == 0
            abs(coef) <= atol || throw(ArgumentError(
                "$context leaves the supplied basis. Missing basis element $mono with coefficient $coef."
            ))
        else
            coords[idx] += coef
        end
    end
    return coords
end

function _fermionic_spin_casimir_matrix(
    sector_basis::AbstractVector{<:NormalMonomial{FermionicAlgebra,T}},
    spin_adaptation::FermionicSpinAdaptationSpec,
    active_modes::AbstractVector{<:Integer};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="fermionic spin block",
) where {T<:Integer}
    sz, splus, sminus = _fermionic_total_spin_operators(T, spin_adaptation, active_modes)
    basis_index = Dict{NormalMonomial{FermionicAlgebra,T},Int}(mono => idx for (idx, mono) in enumerate(sector_basis))
    casimir = zeros(ComplexF64, length(sector_basis), length(sector_basis))

    for (col, mono) in enumerate(sector_basis)
        acted = _adjoint_action(sz, _adjoint_action(sz, mono)) +
            0.5 * (
                _adjoint_action(splus, _adjoint_action(sminus, mono)) +
                _adjoint_action(sminus, _adjoint_action(splus, mono))
            )
        casimir[:, col] .= _fermionic_polynomial_coordinate_vector(
            acted,
            basis_index,
            "$context Casimir action on basis element $col";
            atol,
        )
    end

    return 0.5 * (casimir + casimir')
end

function _group_eigenvalue_indices(values::AbstractVector{<:Real}, tol::Real)
    isempty(values) && return Vector{UnitRange{Int}}()

    groups = UnitRange{Int}[]
    start = 1
    for idx in 2:length(values)
        if abs(values[idx] - values[idx - 1]) > tol
            push!(groups, start:(idx - 1))
            start = idx
        end
    end
    push!(groups, start:length(values))
    return groups
end

function _normalize_phase!(vec::AbstractVector{ComplexF64}; atol::Float64=_SYMMETRY_ATOL)
    for entry in vec
        abs(entry) > atol || continue
        vec ./= entry / abs(entry)
        return vec
    end
    return vec
end

function _canonical_eigenspace_row_basis(
    eigvecs::AbstractMatrix{<:Number};
    atol::Float64=_SYMMETRY_ATOL,
)
    projector = ComplexF64.(eigvecs) * ComplexF64.(eigvecs)'
    qcols = Matrix{ComplexF64}(undef, size(projector, 1), 0)

    for col in axes(projector, 2)
        candidate = projector[:, col]
        if size(qcols, 2) > 0
            candidate -= qcols * (qcols' * candidate)
            candidate -= qcols * (qcols' * candidate)
        end

        norm_candidate = norm(candidate)
        norm_candidate > atol || continue
        candidate ./= norm_candidate
        _normalize_phase!(candidate; atol)
        qcols = hcat(qcols, candidate)
        size(qcols, 2) == size(eigvecs, 2) && break
    end

    size(qcols, 2) == size(eigvecs, 2) || throw(ArgumentError(
        "Failed to build a deterministic basis for a degenerate SU(2) eigenspace of dimension $(size(eigvecs, 2))."
    ))

    return qcols'
end

function _casimir_eigenvalue_to_spin2(value::Real, tol::Real)
    value < -tol && throw(ArgumentError(
        "Encountered a negative SU(2) Casimir eigenvalue $value beyond tolerance $tol."
    ))
    twoS = round(Int, sqrt(max(0.0, 1 + 4 * value)) - 1)
    expected = twoS * (twoS + 2) / 4
    abs(value - expected) <= tol || throw(ArgumentError(
        "Casimir eigenvalue $value does not match any exact `S(S+1)` value within tolerance $tol."
    ))
    return twoS
end

@inline _spin_transform_provenance(prev::Symbol) =
    prev == :sector_wedderburn ? :sector_spin_wedderburn : :sector_spin

function _fermionic_spin_transform_blocks(
    sector_basis::AbstractVector{<:NormalMonomial{FermionicAlgebra,T}},
    transform_blocks::Vector{<:_BasisTransformBlock},
    sector_label::FermionicSectorLabel,
    spin_adaptation::FermionicSpinAdaptationSpec,
    active_modes::AbstractVector{<:Integer};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="fermionic spin block",
) where {T<:Integer}
    casimir = _fermionic_spin_casimir_matrix(
        sector_basis,
        spin_adaptation,
        active_modes;
        atol,
        context,
    )
    counts = Dict{Int,Int}()
    adapted = _BasisTransformBlock[]

    for block in transform_blocks
        row_basis = ComplexF64.(block.row_basis)
        restricted = 0.5 * (row_basis * casimir * row_basis' + (row_basis * casimir * row_basis')')
        eig = eigen(Hermitian(restricted))
        groups = _group_eigenvalue_indices(real.(eig.values), spin_adaptation.casimir_tol)

        for group in groups
            eigvalue = sum(real.(eig.values[group])) / length(group)
            total_spin2 = _casimir_eigenvalue_to_spin2(eigvalue, spin_adaptation.casimir_tol)
            local_row_basis = _canonical_eigenspace_row_basis(
                eig.vectors[:, group];
                atol=max(atol, spin_adaptation.casimir_tol),
            )
            multiplicity_tag = get(counts, total_spin2, 0) + 1
            counts[total_spin2] = multiplicity_tag
            push!(adapted, _BasisTransformBlock(
                local_row_basis * row_basis,
                FermionicSpinBlockLabel(sector_label, total_spin2, multiplicity_tag),
                _spin_transform_provenance(block.provenance),
            ))
        end
    end

    return adapted
end

function _fermionic_sector_transform_blocks(
    sector_basis::AbstractVector{<:NormalMonomial{FermionicAlgebra,T}},
    sector_label::FermionicSectorLabel,
    sw_group,
    spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec},
    active_modes::AbstractVector{<:Integer};
    atol::Float64=_SYMMETRY_ATOL,
    context::AbstractString="fermionic sector block",
) where {T<:Integer}
    n = length(sector_basis)
    transform_blocks = if isnothing(sw_group) || n == 1
        _BasisTransformBlock[_BasisTransformBlock(Matrix{Float64}(I, n, n), sector_label, :sector_split)]
    else
        _BasisTransformBlock[
            _BasisTransformBlock(Matrix(block), sector_label, :sector_wedderburn)
            for block in _sw_decompose_half_basis(collect(sector_basis), sw_group)
        ]
    end

    isnothing(spin_adaptation) && return transform_blocks
    return _fermionic_spin_transform_blocks(
        sector_basis,
        transform_blocks,
        sector_label,
        spin_adaptation,
        active_modes;
        atol,
        context,
    )
end
