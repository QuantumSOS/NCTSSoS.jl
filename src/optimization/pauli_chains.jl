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
    heisenberg_chain_symmetry_spec((σx, σy, σz); translation=true, reflection=true, axis_rotations=true, ...)

Build a `SymmetrySpec` for a periodic Heisenberg XXX chain from the common
structural generators: lattice translation, reflection, and global Pauli-axis
rotations.  This is a convenience wrapper around `CliffordSymmetry`.
"""
function heisenberg_chain_symmetry_spec(
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector};
    translation::Bool=true,
    reflection::Bool=true,
    axis_rotations::Bool=true,
    check_invariance::Bool=true,
    offblock_check::Symbol=:randomized,
)
    generators = CliffordSymmetry[]
    translation && push!(generators, pauli_chain_translation(ops))
    reflection && push!(generators, pauli_chain_reflection(ops))
    axis_rotations && append!(generators, pauli_global_axis_rotation_generators(ops))
    isempty(generators) && throw(ArgumentError("At least one Heisenberg chain symmetry generator must be enabled."))
    return SymmetrySpec(generators; check_invariance, offblock_check)
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
    momentum_sectors::Vector{Int}
    sign_symmetry::Bool
    reflection::Bool
    conjugate_symmetry::Bool
    psd_block_sizes::Vector{Int}
    block_labels::Vector{Any}
    n_unique_moment_matrix_elements::Int
    real_moment_matrix::Bool
end

function Base.show(io::IO, report::TranslationInvariantReport)
    print(
        io,
        "TranslationInvariantReport(n_sites=$(report.n_sites), order=$(report.order), " *
        "basis_size=$(report.basis_size), orbit_basis_size=$(report.orbit_basis_size), " *
        "momentum_sectors=$(report.momentum_sectors), sign_symmetry=$(report.sign_symmetry), " *
        "reflection=$(report.reflection), conjugate_symmetry=$(report.conjugate_symmetry), " *
        "psd_block_sizes=$(report.psd_block_sizes), " *
        "n_unique_moment_matrix_elements=$(report.n_unique_moment_matrix_elements), " *
        "real_moment_matrix=$(report.real_moment_matrix))"
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
end

function Base.show(io::IO, result::TranslationInvariantResult)
    println(io, "Translation-Invariant Optimization Result")
    println(io, "Objective: ", result.objective)
    print(io, "Report: ", result.report)
end

"""
    pauli_translation_invariant_moment_relaxation(pop, (σx, σy, σz), order; sign_symmetry=true, momenta=nothing, rdm_levels=Int[])

Construct a periodic-chain Pauli moment relaxation directly in translation
(momentum) sectors, without building the full site-space moment matrix.

This is an intentionally narrow large-spin-chain path:
- ordinary unconstrained `PauliAlgebra` polynomial objectives only;
- periodic translation by one site on the declared chain `1:N`;
- a contiguous local half-basis from [`pauli_contiguous_chain_basis`](@ref);
- optional `(ℤ₂)^2` Heisenberg sign-symmetry splitting, enabled by default;
- optional mirror (`reflection`) and conjugation (`conjugate_symmetry`) moment
  replacement rules, both enabled by default.

By default this builds the paper-style real primal moment matrix: conjugate
momenta are not duplicated and moment variables are real.  Set
`real_moment_matrix=false` only for debugging the older complex Hermitian
block form.

When `sign_symmetry=true`, the objective must be invariant under the global
Heisenberg sign flips.  When `reflection=true`, the objective must be invariant
under the mirror `site i ↦ N + 1 - i` and the basis must be mirror-closed;
moments are then identified over dihedral (translation + mirror) orbits.  When
`conjugate_symmetry=true`, every objective term must have a real coefficient and
an even number of σʸ letters; moments of monomials with an odd σʸ count are then
identically zero and are dropped.  With both rules and `real_moment_matrix=true`,
each complex momentum block is rotated by an explicit phase-adapted mirror
transform into an exactly real PSD block of the same side (instead of the
doubled realified side), and the `k = 0` sector additionally splits into
mirror-even/odd blocks.  If `momenta` is supplied with
`real_moment_matrix=false`, it must include sector `0` because the normalized
identity moment lives there.

`check_invariance=false` skips the translation *and* reflection invariance
checks on the objective while the corresponding moment identifications are
still applied, so it asserts that the caller has independently verified every
enabled symmetry; passing a non-invariant objective with it produces an
invalid relaxation.

Each `k` in `rdm_levels` adds positivity of the contiguous `k`-site reduced
density matrix on sites `1:k`, split into fixed-magnetization principal blocks.
Translation invariance makes one window sufficient.  The blocks use the full
complex Pauli representation and are realified exactly when
`real_moment_matrix=true`; the positive factor `2^-k` is omitted because it
does not change the PSD constraint.

For the XXX chain with `N=100, order=4`, the default basis has 12,001 site-space
monomials and the solver-facing real PSD blocks have side at most 31, matching
the specialized construction of arXiv:2604.01555.
"""
function pauli_translation_invariant_moment_relaxation(
    pop::PolyOpt{PauliAlgebra,T,P},
    ops::Tuple{<:AbstractVector,<:AbstractVector,<:AbstractVector},
    order::Integer;
    basis::Union{Nothing,Vector{NormalMonomial{PauliAlgebra,T}}}=nothing,
    momenta::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    sign_symmetry::Bool=true,
    reflection::Bool=true,
    conjugate_symmetry::Bool=true,
    check_invariance::Bool=true,
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    rdm_levels::AbstractVector{<:Integer}=Int[],
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    σx, _, _, n = _validate_pauli_chain_ops(ops)
    eltype(σx) == NormalMonomial{PauliAlgebra,T} || throw(ArgumentError(
        "PolyOpt and Pauli chain operators must use the same Pauli integer type; got objective type $T and operator type $(eltype(σx))."
    ))
    isempty(pop.eq_constraints) || throw(ArgumentError("Translation-invariant Pauli path currently supports unconstrained objectives only; equality constraints are not yet supported."))
    isempty(pop.ineq_constraints) || throw(ArgumentError("Translation-invariant Pauli path currently supports unconstrained objectives only; inequality constraints are not yet supported."))
    isempty(pop.moment_eq_constraints) || throw(ArgumentError("Translation-invariant Pauli path currently supports unconstrained objectives only; moment equality constraints are not yet supported."))

    d = Int(order)
    d >= 0 || throw(DomainError(order, "`order` must be non-negative."))
    rdm_ks = _validate_pauli_rdm_levels(rdm_levels, n)
    local_basis = isnothing(basis) ? pauli_contiguous_chain_basis(ops, d; periodic=true) : basis
    one_mono = one(first(σx))
    one_mono in local_basis || throw(ArgumentError("Translation-invariant Pauli basis must include the identity."))

    _check_pauli_chain_support(pop.objective, n; context="objective")
    _check_pauli_chain_support(local_basis, n; context="basis")
    check_invariance && _check_translation_invariance(pop.objective, n)
    sign_symmetry && _check_pauli_sign_invariance(pop.objective)
    real_moment_matrix && _check_real_pauli_chain_objective(pop.objective)
    reflection && check_invariance && _check_reflection_invariance(pop.objective, n)
    conjugate_symmetry && _check_conjugate_invariant_objective(pop.objective)
    _check_translation_basis_closure(local_basis, n)
    reflection && _check_reflection_basis_closure(local_basis, n)

    orbit_reps = _translation_orbit_representatives(local_basis, n)
    nontrivial_reps = [m for m in orbit_reps if !isone(m)]
    sectors = _pauli_chain_momentum_sectors(n, momenta; real_moment_matrix)
    real_moment_matrix || 0 in sectors || throw(ArgumentError("Momentum sector 0 is required because it carries the normalized identity moment."))

    MP_R = _pauli_chain_real_coeff_type(C)
    MP_C = Complex{MP_R}
    MP_P = Polynomial{PauliAlgebra,T,real_moment_matrix ? MP_R : MP_C}
    BLOCK_P = Polynomial{PauliAlgebra,T,MP_C}
    reducer = _PauliMomentReducer{NormalMonomial{PauliAlgebra,T}}(n, reflection, conjugate_symmetry)
    objective_mp = convert(MP_P, _reduce_moment_polynomial(pop.objective, reducer))

    mirror_transform = reflection && conjugate_symmetry && real_moment_matrix
    pairing = reflection ? _mirror_pairing(orbit_reps, n) : nothing

    translated = Dict{NormalMonomial{PauliAlgebra,T},Vector{NormalMonomial{PauliAlgebra,T}}}()
    for rep in nontrivial_reps
        translated[rep] = [_translate_pauli_monomial(rep, r, n) for r in 0:(n - 1)]
    end
    translated[one_mono] = fill(one_mono, n)

    constraints = Tuple{Symbol,Matrix{MP_P}}[]
    moment_terms = NormalMonomial{PauliAlgebra,T}[]
    block_sizes = Int[]
    block_labels = Any[]

    for k in sectors
        sector_basis = k == 0 ? orbit_reps : nontrivial_reps
        blocks = sign_symmetry ? _pauli_signature_blocks(sector_basis) : [(:all, sector_basis)]
        for (signature, block_basis) in blocks
            isempty(block_basis) && continue
            complex_mat = _translation_momentum_block(block_basis, k, n, translated, reducer, BLOCK_P)
            if mirror_transform
                for (parity, mat) in _mirror_real_blocks(complex_mat, block_basis, pairing, k, n, MP_P; atol=MP_R(phase_atol))
                    push!(constraints, (:PSD, mat))
                    push!(block_sizes, size(mat, 1))
                    label = parity === :none ? (momentum=k, signature=signature) :
                        (momentum=k, signature=signature, parity=parity)
                    push!(block_labels, label)
                    for entry in mat
                        append!(moment_terms, monomials(entry))
                    end
                end
            else
                cone, mat = real_moment_matrix ?
                    (:PSD, _realify_hermitian_block(complex_mat, MP_P; atol=MP_R(phase_atol))) :
                    (:HPSD, map(p -> convert(MP_P, p), complex_mat))
                push!(constraints, (cone, mat))
                push!(block_sizes, size(mat, 1))
                push!(block_labels, (momentum=k, signature=signature))
                for entry in mat
                    append!(moment_terms, monomials(entry))
                end
            end
        end
    end

    for k in rdm_ks
        for (down_spins, complex_mat) in _pauli_rdm_blocks(
            k,
            reducer,
            BLOCK_P;
            sign_symmetry,
        )
            cone, mat = real_moment_matrix ?
                (:PSD, _realify_hermitian_block(complex_mat, MP_P; atol=MP_R(phase_atol))) :
                (:HPSD, map(p -> convert(MP_P, p), complex_mat))
            push!(constraints, (cone, mat))
            push!(block_sizes, size(mat, 1))
            push!(block_labels, (rdm=k, down_spins=down_spins))
            for entry in mat
                append!(moment_terms, monomials(entry))
            end
        end
    end

    moment_basis = sorted_unique!(moment_terms)
    _check_objective_moments_covered(objective_mp, moment_basis)
    total_basis = sorted_unique!(vcat(monomials(objective_mp), moment_basis))
    n_unique = length(moment_basis)

    mp = MomentProblem{PauliAlgebra,T,NormalMonomial{PauliAlgebra,T},MP_P}(
        objective_mp,
        constraints,
        total_basis,
        n_unique;
        real_moments=real_moment_matrix,
    )

    report = TranslationInvariantReport(
        n,
        d,
        length(local_basis),
        length(orbit_reps),
        sectors,
        sign_symmetry,
        reflection,
        conjugate_symmetry,
        block_sizes,
        block_labels,
        n_unique,
        real_moment_matrix,
    )
    return mp, report
end

"""
    pauli_translation_invariant_nctssos(pop, (σx, σy, σz), order, optimizer; kwargs...)

Build and solve the translation-invariant Pauli-chain relaxation from
[`pauli_translation_invariant_moment_relaxation`](@ref).

The `dualize` keyword selects the dual SOS or primal moment formulation; at
`N=100`, the dual form peaked near 270 GiB RSS versus about 9 GiB for the
default moment form.
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
    kwargs...,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    mp, report = pauli_translation_invariant_moment_relaxation(pop, ops, order; kwargs...)
    result = solve_sdp(
        mp,
        optimizer;
        dualize,
        formulation,
        representation,
        orphan_policy,
    )
    return TranslationInvariantResult(result.objective, result.model, mp, report)
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

function _validate_pauli_rdm_levels(
    rdm_levels::AbstractVector{<:Integer},
    n_sites::Integer,
)
    levels = sort!(unique!(Int[k for k in rdm_levels]))
    for k in levels
        1 <= k <= Int(n_sites) || throw(ArgumentError(
            "Each `rdm_levels` entry must lie in 1:n_sites; got k=$k for n_sites=$(Int(n_sites))."
        ))
        k < 8sizeof(Int) - 1 || throw(ArgumentError(
            "RDM level $k is too large for computational-basis indexing on this platform."
        ))
    end
    return levels
end

function _pauli_magnetization_states(k::Int, down_spins::Int)
    return Int[state for state in 0:((1 << k) - 1) if count_ones(UInt(state)) == down_spins]
end

function _pauli_rdm_entry(
    row_state::Int,
    col_state::Int,
    k::Int,
    reducer,
    ::Type{P};
    sign_symmetry::Bool,
) where {
    T<:Unsigned,
    P<:Polynomial{PauliAlgebra,T},
}
    C = _coefficient_type(P)
    M = NormalMonomial{PauliAlgebra,T}
    terms = Tuple{C,M}[]
    sizehint!(terms, 1 << max(k - (sign_symmetry ? 1 : 0) - (reducer.drop_odd_y ? 1 : 0), 0))
    word = T[]

    for choices in 0:((1 << k) - 1)
        y_odd = false
        z_odd = false
        differing = 0
        for site in 1:k
            bit_shift = k - site
            row_bit = !iszero((row_state >> bit_shift) & 1)
            col_bit = !iszero((col_state >> bit_shift) & 1)
            choose_second = !iszero((choices >> bit_shift) & 1)
            if row_bit == col_bit
                choose_second && (z_odd = !z_odd)
            else
                differing += 1
                choose_second && (y_odd = !y_odd)
            end
        end

        x_odd = xor(isodd(differing), y_odd)
        sign_symmetry && !(x_odd == y_odd == z_odd) && continue
        reducer.drop_odd_y && y_odd && continue

        empty!(word)
        coef = one(C)
        for site in 1:k
            bit_shift = k - site
            row_bit = !iszero((row_state >> bit_shift) & 1)
            col_bit = !iszero((col_state >> bit_shift) & 1)
            choose_second = !iszero((choices >> bit_shift) & 1)
            if row_bit == col_bit
                if choose_second
                    push!(word, _pauli_index_from_site_type(T, site, _PAULI_Z_TYPE))
                    row_bit && (coef = -coef)
                end
            elseif choose_second
                push!(word, _pauli_index_from_site_type(T, site, _PAULI_Y_TYPE))
                coef *= row_bit ? C(im) : C(-im)
            else
                push!(word, _pauli_index_from_site_type(T, site, _PAULI_X_TYPE))
            end
        end

        mono = M(copy(word))
        rep = _reduce_moment(reducer, mono)
        rep === nothing || push!(terms, (coef, rep))
    end
    return P(terms)
end

function _pauli_rdm_blocks(
    k::Int,
    reducer,
    ::Type{P};
    sign_symmetry::Bool,
) where {P<:Polynomial}
    blocks = Tuple{Int,Matrix{P}}[]
    for down_spins in 0:fld(k, 2)
        states = _pauli_magnetization_states(k, down_spins)
        dim = length(states)
        mat = Matrix{P}(undef, dim, dim)
        for j in 1:dim, i in 1:j
            entry = _pauli_rdm_entry(
                states[i],
                states[j],
                k,
                reducer,
                P;
                sign_symmetry,
            )
            if i == j
                mat[i, i] = (entry + adjoint(entry)) / 2
            else
                mat[i, j] = entry
                mat[j, i] = adjoint(entry)
            end
        end
        push!(blocks, (down_spins, mat))
    end
    return blocks
end

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

function _check_real_pauli_chain_objective(poly::Polynomial{PauliAlgebra,T,C}) where {T<:Unsigned,C<:Number}
    for (coef, mono) in poly.terms
        iszero(imag(coef)) || throw(ArgumentError(
            "real_moment_matrix=true requires real objective coefficients; term $mono has coefficient $coef."
        ))
    end
    return nothing
end

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

function _check_translation_invariance(poly::Polynomial{PauliAlgebra,T,C}, n_sites::Integer) where {T<:Unsigned,C<:Number}
    images = Dict{NormalMonomial{PauliAlgebra,T},Tuple{Int,NormalMonomial{PauliAlgebra,T}}}()
    for site in 1:Int(n_sites), ptype in 0:2
        src = NormalMonomial{PauliAlgebra,T}(T[_pauli_index_from_site_type(T, site, ptype)])
        dst = NormalMonomial{PauliAlgebra,T}(T[_pauli_index_from_site_type(T, mod1(site + 1, Int(n_sites)), ptype)])
        images[src] = (1, dst)
    end
    translation = CliffordSymmetry(images; nqubits=Int(n_sites))
    _act_polynomial(translation, poly) == poly || throw(ArgumentError(
        "Translation-invariant Pauli relaxation requires a one-site translation-invariant objective."
    ))
    return nothing
end

"""
    _reflect_pauli_monomial(mono, n_sites)

Apply the chain mirror `site i ↦ n_sites + 1 - i` to a normal-form Pauli word.
Site permutations of normal words never produce a phase because all letters
live on distinct sites.
"""
function _reflect_pauli_monomial(
    mono::NormalMonomial{PauliAlgebra,T},
    n_sites::Integer,
) where {T<:Unsigned}
    isempty(mono.word) && return mono
    n = Int(n_sites)
    word = Vector{T}(undef, length(mono.word))
    for (i, idx) in pairs(mono.word)
        site = _pauli_site(idx)
        ptype = _pauli_type(idx)
        word[i] = _pauli_index_from_site_type(T, n + 1 - site, ptype)
    end
    simplified, phase = simplify(PauliAlgebra, word)
    phase == 0x00 || throw(ArgumentError("Internal error: Pauli reflection produced phase $phase."))
    return NormalMonomial{PauliAlgebra,T}(copy(simplified))
end

_pauli_y_count(mono::NormalMonomial{PauliAlgebra}) =
    count(idx -> _pauli_type(idx) == _PAULI_Y_TYPE, mono.word)

"""
    _PauliMomentReducer

Reduces monomials appearing in moment-matrix entries to canonical moment
representatives: translation-orbit representatives, optionally identified over
mirror images (`reflection`), with moments of odd σʸ count dropped as
identically zero (`drop_odd_y`, valid under conjugation symmetry).
"""
struct _PauliMomentReducer{M<:NormalMonomial{PauliAlgebra}}
    n::Int
    reflection::Bool
    drop_odd_y::Bool
    cache::Dict{M,Union{M,Nothing}}
end

function _PauliMomentReducer{M}(n::Integer, reflection::Bool, drop_odd_y::Bool) where {M}
    return _PauliMomentReducer{M}(Int(n), reflection, drop_odd_y, Dict{M,Union{M,Nothing}}())
end

function _reduce_moment(reducer::_PauliMomentReducer{M}, mono::M) where {M}
    return get!(reducer.cache, mono) do
        reducer.drop_odd_y && isodd(_pauli_y_count(mono)) && return nothing
        rep = _translation_orbit_representative(mono, reducer.n)
        if reducer.reflection
            alt = _translation_orbit_representative(
                _reflect_pauli_monomial(mono, reducer.n), reducer.n
            )
            alt < rep && (rep = alt)
        end
        return rep
    end
end

function _reduce_moment_polynomial(
    poly::Polynomial{PauliAlgebra,T,C},
    reducer::_PauliMomentReducer{NormalMonomial{PauliAlgebra,T}},
) where {T<:Unsigned,C<:Number}
    terms = Tuple{C,NormalMonomial{PauliAlgebra,T}}[]
    sizehint!(terms, length(poly.terms))
    for (coef, mono) in poly.terms
        rep = _reduce_moment(reducer, mono)
        rep === nothing && continue
        push!(terms, (coef, rep))
    end
    return Polynomial(terms)
end

"""
    _mirror_pairing(reps, n_sites)

For each translation-orbit representative `a`, find the representative `â` of
its mirror image and the shift `s` with `ω(a) = T^s(â)`.  The map `a ↦ â` is an
involution with matching shifts (`s_â = s_a`).
"""
function _mirror_pairing(
    reps::Vector{M},
    n_sites::Integer,
) where {M<:NormalMonomial{PauliAlgebra}}
    n = Int(n_sites)
    partner = Dict{M,M}()
    shift = Dict{M,Int}()
    for a in reps
        image = _reflect_pauli_monomial(a, n)
        ahat = _translation_orbit_representative(image, n)
        s = -1
        for r in 0:(n - 1)
            if _translate_pauli_monomial(ahat, r, n) == image
                s = r
                break
            end
        end
        s >= 0 || throw(ArgumentError("Internal error: mirror image of $a is not a translate of its orbit representative."))
        partner[a] = ahat
        shift[a] = s
    end
    return (partner=partner, shift=shift)
end

function _check_reflection_invariance(
    poly::Polynomial{PauliAlgebra,T,C},
    n_sites::Integer,
) where {T<:Unsigned,C<:Number}
    n = Int(n_sites)
    reflected = Polynomial([(coef, _reflect_pauli_monomial(mono, n)) for (coef, mono) in poly.terms])
    reflected == poly || throw(ArgumentError(
        "reflection=true requires the objective to be invariant under the chain mirror site i ↦ $(n) + 1 - i. Pass reflection=false for non-invariant objectives."
    ))
    return nothing
end

function _check_reflection_basis_closure(
    basis::Vector{NormalMonomial{PauliAlgebra,T}},
    n_sites::Integer,
) where {T<:Unsigned}
    basis_set = Set(basis)
    for mono in basis
        image = _reflect_pauli_monomial(mono, n_sites)
        image in basis_set || throw(ArgumentError(
            "reflection=true requires a mirror-closed basis; $mono maps to $image outside the basis. Pass reflection=false or supply a mirror-closed basis."
        ))
    end
    return nothing
end

function _check_conjugate_invariant_objective(
    poly::Polynomial{PauliAlgebra,T,C},
) where {T<:Unsigned,C<:Number}
    for (coef, mono) in poly.terms
        iszero(coef) && continue
        (iszero(imag(coef)) && iseven(_pauli_y_count(mono))) || throw(ArgumentError(
            "conjugate_symmetry=true requires a conjugation-invariant objective (real coefficients and an even number of σʸ letters per term); term $mono with coefficient $coef violates this. Pass conjugate_symmetry=false."
        ))
    end
    return nothing
end

function _check_pauli_sign_invariance(poly::Polynomial{PauliAlgebra,T,C}) where {T<:Unsigned,C<:Number}
    for (coef, mono) in poly.terms
        iszero(coef) && continue
        _pauli_sign_signature(mono) == 0x00 || throw(ArgumentError(
            "sign_symmetry=true requires the objective to be invariant under global Heisenberg sign flips; term $mono has nontrivial signature. Pass sign_symmetry=false for non-invariant objectives."
        ))
    end
    return nothing
end

function _check_objective_moments_covered(
    objective::Polynomial{PauliAlgebra,T,C},
    moment_basis::Vector{NormalMonomial{PauliAlgebra,T}},
) where {T<:Unsigned,C<:Number}
    moment_set = Set(moment_basis)
    missing = NormalMonomial{PauliAlgebra,T}[]
    for mono in monomials(objective)
        mono in moment_set || push!(missing, mono)
    end
    isempty(missing) && return nothing
    shown = join((sprint(show, mono) for mono in Iterators.take(missing, 5)), ", ")
    length(missing) > 5 && (shown *= ", ...")
    throw(ArgumentError(
        "Objective contains $(length(missing)) moment(s) not generated by the translation-invariant PSD blocks: [$shown]. Increase `order`, adjust `basis`, or disable incompatible symmetry reductions."
    ))
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

@inline function _momentum_phase(::Type{R}, k::Int, r::Int, n::Int) where {R<:Real}
    (iszero(k) || iszero(r)) && return complex(one(R), zero(R))
    θ = -R(2) * R(pi) * R(k) * R(r) / R(n)
    return cis(θ)
end

function _translation_momentum_entry(
    row::M,
    col::M,
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    reducer::_PauliMomentReducer{M},
    ::Type{P},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},P<:Polynomial{PauliAlgebra,T}}
    if isone(row) && isone(col)
        return P([(one(_coefficient_type(P)), row)])
    end

    C = _coefficient_type(P)
    R = typeof(real(one(C)))
    cross_scale = xor(isone(row), isone(col)) ? C(inv(sqrt(R(n)))) : one(C)
    terms = Tuple{C,M}[]
    sizehint!(terms, n)
    col_translates = translated[col]
    for r in 0:(n - 1)
        phase = C(_momentum_phase(R, k, r, n))
        prod = row * col_translates[r + 1]
        for (coef, mono) in prod.terms
            rep = _reduce_moment(reducer, mono)
            rep === nothing && continue
            push!(terms, (cross_scale * phase * C(coef), rep))
        end
    end
    return convert(P, Polynomial(terms))
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

"""
    _mirror_columns(block_basis, pairing, k, n_sites, CT)

Build the columns of the phase-adapted mirror transform `W` for momentum
sector `k`.  Each orbit representative `a` receives the phase
`exp(i·(π k s_a / n + π y_a / 2))` where `s_a` is its mirror shift and `y_a`
its σʸ count; mirror pairs additionally combine into even (`(u_a + u_â)/√2`)
and odd (`i(u_a - u_â)/√2`) columns.  `W` is unitary and `W' M_k W` is exactly
real under dihedral moment identification with odd-σʸ moments dropped.
"""
function _mirror_columns(
    block_basis::Vector{M},
    pairing,
    k::Int,
    n_sites::Int,
    ::Type{CT},
) where {M<:NormalMonomial{PauliAlgebra},CT<:Complex}
    R = typeof(real(one(CT)))
    index = Dict{M,Int}(mono => i for (i, mono) in enumerate(block_basis))
    cols_even = Vector{Vector{Tuple{Int,CT}}}()
    cols_odd = Vector{Vector{Tuple{Int,CT}}}()
    paired = falses(length(block_basis))
    for (i, a) in enumerate(block_basis)
        paired[i] && continue
        paired[i] = true
        ahat = pairing.partner[a]
        s = pairing.shift[a]
        α = R(pi) * R(k) * R(s) / R(n_sites) + R(pi) / 2 * R(_pauli_y_count(a))
        c = CT(cis(α))
        if ahat == a
            push!(cols_even, [(i, c)])
        else
            j = get(index, ahat, 0)
            j > 0 || throw(ArgumentError("Internal error: mirror partner $ahat of $a is missing from its momentum block."))
            paired[j] = true
            scale = CT(inv(sqrt(R(2))))
            push!(cols_even, [(i, c * scale), (j, c * scale)])
            push!(cols_odd, [(i, im * c * scale), (j, -im * c * scale)])
        end
    end
    return cols_even, cols_odd
end

function _mirror_transformed_entry(
    mat::Matrix{P},
    ci::Vector{Tuple{Int,CT}},
    cj::Vector{Tuple{Int,CT}},
) where {P<:Polynomial,CT<:Complex}
    acc = nothing
    for (p, cp) in ci, (q, cq) in cj
        contrib = (conj(cp) * cq) * mat[p, q]
        acc = acc === nothing ? contrib : acc + contrib
    end
    return acc::P
end

function _mirror_transformed_block(
    mat::Matrix{PC},
    cols::Vector{Vector{Tuple{Int,CT}}},
    ::Type{PReal};
    atol,
) where {PC<:Polynomial,CT<:Complex,PReal<:Polynomial}
    m = length(cols)
    out = Matrix{PReal}(undef, m, m)
    for j in 1:m, i in 1:j
        entry = _mirror_transformed_entry(mat, cols[i], cols[j])
        imag_entry = _imag_part_polynomial(entry, PReal; atol)
        iszero(imag_entry) || throw(ArgumentError(
            "Internal error: mirror-adapted momentum block has a non-real entry; residual imaginary part $(imag_entry)."
        ))
        re = _real_part_polynomial(entry, PReal; atol)
        out[i, j] = re
        i == j || (out[j, i] = re)
    end
    return out
end

function _assert_mirror_zero_cross(
    mat::Matrix{PC},
    cols_even::Vector{Vector{Tuple{Int,CT}}},
    cols_odd::Vector{Vector{Tuple{Int,CT}}},
    ::Type{PReal};
    atol,
) where {PC<:Polynomial,CT<:Complex,PReal<:Polynomial}
    for ce in cols_even, co in cols_odd
        entry = _mirror_transformed_entry(mat, ce, co)
        re = _real_part_polynomial(entry, PReal; atol)
        imag_entry = _imag_part_polynomial(entry, PReal; atol)
        (iszero(re) && iszero(imag_entry)) || throw(ArgumentError(
            "Internal error: mirror-adapted k = 0 block has a nonzero even/odd cross entry $(re + imag_entry)."
        ))
    end
    return nothing
end

"""
    _mirror_real_blocks(complex_mat, block_basis, pairing, k, n_sites, PReal; atol)

Rotate a complex Hermitian momentum block into exactly real PSD blocks using
the phase-adapted mirror transform.  For `k = 0` the block splits into
mirror-even and mirror-odd sub-blocks (their cross entries vanish and are
asserted to); all other sectors yield a single real block of the same side.
"""
function _mirror_real_blocks(
    complex_mat::Matrix{PC},
    block_basis::Vector{M},
    pairing,
    k::Int,
    n_sites::Int,
    ::Type{PReal};
    atol,
) where {PC<:Polynomial,M<:NormalMonomial{PauliAlgebra},PReal<:Polynomial}
    CT = _coefficient_type(PC)
    cols_even, cols_odd = _mirror_columns(block_basis, pairing, k, n_sites, CT)
    if k == 0
        blocks = Tuple{Symbol,Matrix{PReal}}[]
        isempty(cols_even) || push!(blocks, (:even, _mirror_transformed_block(complex_mat, cols_even, PReal; atol)))
        if !isempty(cols_odd)
            push!(blocks, (:odd, _mirror_transformed_block(complex_mat, cols_odd, PReal; atol)))
            _assert_mirror_zero_cross(complex_mat, cols_even, cols_odd, PReal; atol)
        end
        return blocks
    end
    cols = vcat(cols_even, cols_odd)
    return Tuple{Symbol,Matrix{PReal}}[(:none, _mirror_transformed_block(complex_mat, cols, PReal; atol))]
end

function _translation_momentum_block(
    block_basis::Vector{M},
    k::Int,
    n::Int,
    translated::Dict{M,Vector{M}},
    reducer::_PauliMomentReducer{M},
    ::Type{P},
) where {T<:Unsigned,M<:NormalMonomial{PauliAlgebra,T},P<:Polynomial{PauliAlgebra,T}}
    m = length(block_basis)
    mat = Matrix{P}(undef, m, m)
    for j in 1:m, i in 1:j
        entry = _translation_momentum_entry(block_basis[i], block_basis[j], k, n, translated, reducer, P)
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
