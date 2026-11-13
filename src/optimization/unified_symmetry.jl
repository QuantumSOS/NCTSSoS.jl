# Unified symmetry front.
#
# This file adds no reduction kernel. `Symmetry` records what the user supplies,
# `check_symmetry` and `cs_nctssos(pop, config, symmetry)` select one existing
# backend, call it with unchanged arguments, and report which checks that
# backend ran. Compositions no backend implements are rejected here.

const _FINITE_SYMMETRY_GENERATOR = Union{SignedPermutation,CliffordSymmetry,FermionicModePermutation}

"""
    Symmetry(; finite=[], continuous=[], compact_metric=nothing, kwargs...)
    Symmetry(spec::SymmetrySpec)

One description of the symmetry of an optimization problem, used by
[`check_symmetry`](@ref) and `cs_nctssos(pop, config, symmetry)`.

A symmetry is supplied as generators of two kinds:

- `finite`: invertible maps of the operators, all of one type:
  [`SignedPermutation`](@ref) (monoid algebras), [`CliffordSymmetry`](@ref)
  (Pauli algebra; a Clifford map is a finite map, not a separate kind), or
  [`FermionicModePermutation`](@ref) (fermionic algebra).
- `continuous`: `SuppliedDerivation`s, the rates of change of a continuous
  family of maps (for example the three global spin rotations). They need
  `compact_metric`, a positive definite matrix on the registered letters that
  every rate preserves; it proves that averaging over the family is defined.

`Symmetry` only selects and calls an existing backend:

| supplied                              | route                | backend                         |
|:--------------------------------------|:---------------------|:--------------------------------|
| `finite` and/or sector options        | `:finite`            | [`SymmetrySpec`](@ref) reduction |
| `continuous`                          | `:continuous`        | `continuous_moment_relaxation`  |
| `CliffordSymmetry` + `continuous`     | `:finite_continuous` | same, with `clifford_generators` |

# Keyword arguments
- `sector`, `spin_adaptation`, `pauli_charge`, `pauli_singlet`: the sector
  options of [`SymmetrySpec`](@ref), forwarded unchanged (`:finite` route only).
- `check_invariance=true`, `offblock_check=:randomized`: the verification
  controls of [`SymmetrySpec`](@ref), forwarded unchanged (`:finite` route only;
  the continuous backend always validates exactly and has no switch).
- `compress=true`: forwarded to `continuous_moment_relaxation` (continuous
  routes only).

# Not supported here
The constructor throws an `ArgumentError` instead of approximating:
- antiunitary maps such as `PauliAntiunitary` (time reversal); use
  `antiunitary_moment_relaxation`;
- finite generators of different types in one `Symmetry`;
- continuous generators combined with non-Clifford finite generators or with
  sector options.

The chain-specific routes (`pauli_translation_invariant_nctssos`) keep their own
entry points.

# Example
```julia
sym = Symmetry(finite=[pauli_site_permutation([2, 3, 1])],
               continuous=spin_rotations, compact_metric=metric)
report = check_symmetry(pop, config, sym)     # validate and reduce, no solver
solved = cs_nctssos(pop, config, sym)         # solved.objective, solved.report
```
"""
struct Symmetry
    finite::Vector{Any}
    continuous::Vector{SuppliedDerivation}
    compact_metric::Union{Nothing,AbstractMatrix}
    sector::Union{Nothing,FermionicSectorSpec}
    spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}
    pauli_charge::Union{Nothing,PauliChargeSectorSpec}
    pauli_singlet::Union{Nothing,PauliSingletConstraintSpec}
    check_invariance::Bool
    offblock_check::Symbol
    compress::Bool

    function Symmetry(;
        finite::AbstractVector=Any[],
        continuous::AbstractVector=SuppliedDerivation[],
        compact_metric::Union{Nothing,AbstractMatrix}=nothing,
        sector::Union{Nothing,FermionicSectorSpec}=nothing,
        spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}=nothing,
        pauli_charge::Union{Nothing,PauliChargeSectorSpec}=nothing,
        pauli_singlet::Union{Nothing,PauliSingletConstraintSpec}=nothing,
        check_invariance::Bool=true,
        offblock_check::Symbol=:randomized,
        compress::Bool=true,
    )
        for g in finite
            g isa PauliAntiunitary && throw(ArgumentError(
                "`Symmetry` covers unitary symmetries only. An antiunitary map such as `PauliAntiunitary` " *
                "is not a finite generator here; use `antiunitary_moment_relaxation`."
            ))
            g isa _FINITE_SYMMETRY_GENERATOR || throw(ArgumentError(
                "Unsupported finite generator of type `$(typeof(g))`. Supply `SignedPermutation`, " *
                "`CliffordSymmetry`, or `FermionicModePermutation` values."
            ))
        end
        kinds = unique!(Symbol[_finite_generator_kind(g) for g in finite])
        length(kinds) <= 1 || throw(ArgumentError(
            "Finite generators must all have one type; got $(join(kinds, ", ")). " *
            "No backend composes different finite action types."
        ))
        all(a -> a isa SuppliedDerivation, continuous) || throw(ArgumentError(
            "Continuous generators must be `SuppliedDerivation` values."
        ))
        offblock_check in (:full, :randomized, :off) || throw(ArgumentError(
            "`offblock_check` must be `:full`, `:randomized`, or `:off`; got `$(repr(offblock_check))`."
        ))
        options = _symmetry_options(sector, spin_adaptation, pauli_charge, pauli_singlet)
        if isempty(continuous)
            isempty(finite) && isempty(options) && throw(ArgumentError(
                "Supply at least one finite generator, continuous generator, or sector option."
            ))
            isnothing(compact_metric) || throw(ArgumentError(
                "`compact_metric` applies to continuous generators only."
            ))
            compress || throw(ArgumentError(
                "`compress=false` applies to continuous generators only."
            ))
            # Run the existing `SymmetrySpec` combination checks now, not at solve time.
            _finite_symmetry_spec(finite, sector, spin_adaptation, pauli_charge, pauli_singlet,
                check_invariance, offblock_check)
        else
            isempty(options) || throw(ArgumentError(
                "Sector options ($(join(options, ", "))) cannot be combined with continuous generators: " *
                "no backend composes them. Use them on a finite-only `Symmetry`."
            ))
            isempty(kinds) || only(kinds) === :clifford || throw(ArgumentError(
                "Finite and continuous generators compose only for `CliffordSymmetry` on Pauli problems; " *
                "got $(only(kinds)) generators."
            ))
            check_invariance && offblock_check === :randomized || throw(ArgumentError(
                "`check_invariance` and `offblock_check` apply to finite-only symmetries. " *
                "The continuous backend always validates exactly and cannot skip its checks."
            ))
        end
        return new(collect(Any, finite), SuppliedDerivation[a for a in continuous], compact_metric,
            sector, spin_adaptation, pauli_charge, pauli_singlet, check_invariance, offblock_check, compress)
    end
end

function Symmetry(spec::SymmetrySpec)
    return Symmetry(;
        finite=Any[spec.generators; spec.fermionic_generators; spec.clifford_generators],
        sector=spec.sector,
        spin_adaptation=spec.spin_adaptation,
        pauli_charge=spec.pauli_charge,
        pauli_singlet=spec.pauli_singlet,
        check_invariance=spec.check_invariance,
        offblock_check=spec.offblock_check,
    )
end

_finite_generator_kind(::SignedPermutation) = :signed_permutation
_finite_generator_kind(::CliffordSymmetry) = :clifford
_finite_generator_kind(::FermionicModePermutation) = :fermionic_mode_permutation

function _symmetry_options(sector, spin_adaptation, pauli_charge, pauli_singlet)
    options = Symbol[]
    isnothing(sector) || push!(options, :sector)
    isnothing(spin_adaptation) || push!(options, :spin_adaptation)
    isnothing(pauli_charge) || push!(options, :pauli_charge)
    isnothing(pauli_singlet) || push!(options, :pauli_singlet)
    return options
end

_symmetry_options(symmetry::Symmetry) = _symmetry_options(
    symmetry.sector, symmetry.spin_adaptation, symmetry.pauli_charge, symmetry.pauli_singlet)

function _symmetry_route(symmetry::Symmetry)
    isempty(symmetry.continuous) && return :finite
    return isempty(symmetry.finite) ? :continuous : :finite_continuous
end

function Base.show(io::IO, symmetry::Symmetry)
    print(io, "Symmetry(route=$(_symmetry_route(symmetry)), finite=$(length(symmetry.finite)), " *
        "continuous=$(length(symmetry.continuous)), options=$(_symmetry_options(symmetry)))")
end

"""
    SymmetrySpec(symmetry::Symmetry)

The [`SymmetrySpec`](@ref) that the `:finite` route of [`Symmetry`](@ref) hands
to the existing finite backend. Throws for symmetries with continuous generators.
"""
function SymmetrySpec(symmetry::Symmetry)
    isempty(symmetry.continuous) || throw(ArgumentError(
        "A `Symmetry` with continuous generators has no `SymmetrySpec` equivalent."
    ))
    return _finite_symmetry_spec(symmetry.finite, symmetry.sector, symmetry.spin_adaptation,
        symmetry.pauli_charge, symmetry.pauli_singlet, symmetry.check_invariance, symmetry.offblock_check)
end

function _finite_symmetry_spec(finite, sector, spin_adaptation, pauli_charge, pauli_singlet,
    check_invariance, offblock_check)
    of_type(G) = G[g for g in finite if g isa G]
    return SymmetrySpec(;
        generators=of_type(SignedPermutation),
        fermionic_generators=of_type(FermionicModePermutation),
        clifford_generators=of_type(CliffordSymmetry),
        sector,
        spin_adaptation,
        pauli_charge,
        pauli_singlet,
        check_invariance,
        offblock_check,
    )
end

"""
    SymmetryCheckReport

What [`check_symmetry`](@ref) or `cs_nctssos(pop, config, symmetry)` did, in the
same terms for every route. Each label states what the selected backend checked,
not what is true of the problem.

# Fields
- `route`: `:finite`, `:continuous`, or `:finite_continuous`.
- `backend`: the existing function that built the reduction.
- `finite_generators`, `continuous_generators`: supplied generator counts.
- `options`: forwarded [`SymmetrySpec`](@ref) sector options.
- `invariance`: how the objective and constraints were shown to be symmetric.
  `:exact`: compared with their images under every group element, coefficient by
  coefficient, without tolerance; `:numeric`: a sector option is present, whose
  neutrality check uses a floating-point tolerance; `:unverified`:
  `check_invariance=false`.
- `compactness`: `:verified` for continuous generators (exact metric check),
  `:not_applicable` otherwise.
- `block_arithmetic`: number type of the block-splitting change of basis,
  `:float64` (finite backend) or `:complex_rationals` (continuous backend).
- `block_verification`: how the dropped couplings between blocks were shown to
  vanish. `:exact`: exact reconstruction of every matrix from its blocks;
  `:numeric`: checked up to floating-point tolerance (`offblock_check=:full` or
  `:randomized`); `:unverified`: `offblock_check=:off`.
- `joint_blocks`: `:not_applicable`, or `:continuous_only` on the
  `:finite_continuous` route, where finite generators reduce the moments but the
  blocks come from the continuous generators alone, so the blocks can be coarser
  than the joint symmetry allows.
- `relation_to_input`: summary of the above. `:equivalent`: same optimum, proved
  by exact checks; `:equivalent_numeric`: same optimum up to the floating-point
  checks; `:unverified`: a check was switched off or an option asserts a property
  that the backend does not check.
- `invariant_moment_count`: free moment coordinates reported by the backend.
- `psd_block_sizes`: sizes of the positive semidefinite blocks given to the solver.
- `notes`: plain statements of skipped checks, asserted properties, and limits.
- `backend_report`: the untouched backend output ([`SymmetryReport`](@ref), or
  the named tuple of `continuous_moment_relaxation`).
"""
struct SymmetryCheckReport
    route::Symbol
    backend::Symbol
    finite_generators::Int
    continuous_generators::Int
    options::Vector{Symbol}
    invariance::Symbol
    compactness::Symbol
    block_arithmetic::Symbol
    block_verification::Symbol
    joint_blocks::Symbol
    relation_to_input::Symbol
    invariant_moment_count::Int
    psd_block_sizes::Vector{Int}
    notes::Vector{String}
    backend_report::Any
end

function Base.show(io::IO, report::SymmetryCheckReport)
    println(io, "SymmetryCheckReport")
    println(io, "  route:              $(report.route) (backend `$(report.backend)`)")
    println(io, "  generators:         $(report.finite_generators) finite, $(report.continuous_generators) continuous" *
        (isempty(report.options) ? "" : "; options $(join(report.options, ", "))"))
    println(io, "  invariance:         $(report.invariance)")
    println(io, "  compactness:        $(report.compactness)")
    println(io, "  blocks:             $(report.block_arithmetic), verification $(report.block_verification)" *
        (report.joint_blocks === :not_applicable ? "" : ", joint blocks $(report.joint_blocks)"))
    println(io, "  relation to input:  $(report.relation_to_input)")
    println(io, "  invariant moments:  $(report.invariant_moment_count)")
    print(io, "  PSD block sizes:    $(report.psd_block_sizes)")
    for note in report.notes
        print(io, "\n  note: $note")
    end
end

function _finite_check_report(symmetry::Symmetry, backend_report::SymmetryReport)
    options = _symmetry_options(symmetry)
    notes = String[]
    !isempty(symmetry.finite) && !symmetry.check_invariance && push!(notes,
        "`check_invariance=false`: invariance of the objective and constraints was not checked.")
    symmetry.offblock_check === :off && push!(notes,
        "`offblock_check=:off`: the couplings dropped between blocks were not checked to vanish.")
    asserted = intersect(options, (:pauli_singlet, :spin_adaptation))
    :pauli_singlet in options && push!(notes,
        "`pauli_singlet` adds moment equations that hold for spin-rotation-invariant states; " *
        "this route does not check that the problem is spin-rotation invariant.")
    :spin_adaptation in options && push!(notes,
        "`spin_adaptation` assumes spin-rotation invariance; this route does not check it for the problem.")
    isempty(setdiff(options, asserted)) || push!(notes,
        "Sector options $(join(setdiff(options, asserted), ", ")) are forwarded unchanged; " *
        "their own neutrality checks are documented on `SymmetrySpec`.")
    generators_unchecked = !isempty(symmetry.finite) && !symmetry.check_invariance
    verified = !generators_unchecked && symmetry.offblock_check !== :off && isempty(asserted)
    return SymmetryCheckReport(
        :finite,
        :moment_relax_symmetric,
        length(symmetry.finite),
        0,
        options,
        generators_unchecked ? :unverified : isempty(options) ? :exact : :numeric,
        :not_applicable,
        :float64,
        symmetry.offblock_check === :off ? :unverified : :numeric,
        :not_applicable,
        verified ? :equivalent_numeric : :unverified,
        backend_report.invariant_moment_count,
        copy(backend_report.psd_block_sizes),
        notes,
        backend_report,
    )
end

function _continuous_check_report(symmetry::Symmetry, reduction)
    mixed = !isempty(symmetry.finite)
    notes = String[]
    mixed && push!(notes,
        "Finite generators reduce the moment coordinates; the blocks use the continuous generators only " *
        "and may be coarser than the joint symmetry allows.")
    symmetry.compress || push!(notes,
        "`compress=false`: invariant moments only; the positive semidefinite blocks are not split.")
    return SymmetryCheckReport(
        mixed ? :finite_continuous : :continuous,
        :continuous_moment_relaxation,
        length(symmetry.finite),
        length(symmetry.continuous),
        Symbol[],
        :exact,
        reduction.compactness,
        reduction.arithmetic,
        :exact,
        mixed ? :continuous_only : :not_applicable,
        reduction.relation_to_input,
        size(reduction.section, 2),
        Int[block.size for block in reduction.moment_problem.linear.psd_blocks_lin],
        notes,
        reduction,
    )
end

function _require_unset_config_symmetry(config::SolverConfig)
    isnothing(config.symmetry) || throw(ArgumentError(
        "Supply the symmetry once: either `SolverConfig(symmetry=SymmetrySpec(...))` or the `Symmetry` argument, not both."
    ))
    return nothing
end

function _config_with_symmetry_spec(config::SolverConfig, spec::SymmetrySpec)
    return SolverConfig(;
        optimizer=config.optimizer,
        order=config.order,
        moment_basis=config.moment_basis,
        cs_algo=config.cs_algo,
        ts_algo=config.ts_algo,
        symmetry=spec,
    )
end

# Reject before any backend work, so an unsupported composition never reaches a
# code path that would handle part of it.
function _check_continuous_route_support(pop::OptimizationProblem{A,P}, symmetry::Symmetry) where {A<:AlgebraType,P}
    isnothing(symmetry.compact_metric) && throw(ArgumentError(
        "Continuous generators need `compact_metric`: a positive definite matrix on the registered letters " *
        "that every generator preserves. Without it the reduction is not known to be exact."
    ))
    isempty(symmetry.finite) && return nothing
    A === PauliAlgebra || throw(ArgumentError(
        "Finite and continuous generators compose only on Pauli problems; got `$(nameof(A))`. " *
        "Use a finite-only or a continuous-only `Symmetry`."
    ))
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) && isempty(pop.moment_eq_constraints) ||
        throw(ArgumentError(
            "Finite and continuous generators compose only for unconstrained problems. " *
            "Use a finite-only or a continuous-only `Symmetry` for this constrained problem."
        ))
    return nothing
end

_clifford_generators(symmetry::Symmetry) = CliffordSymmetry[g for g in symmetry.finite]

# Keywords that `Symmetry` already carries must not be overridden by a second
# value passed at solve time.
function _reject_symmetry_keywords(kwargs)
    for key in (:compact_metric, :compress, :clifford_generators)
        haskey(kwargs, key) && throw(ArgumentError(
            "`$key` is part of the `Symmetry` argument; do not pass it again as a keyword."
        ))
    end
    return nothing
end

"""
    check_symmetry(pop, config::SolverConfig, symmetry::Symmetry; physics=nothing) -> SymmetryCheckReport

Validate `symmetry` against the whole relaxation of `pop` selected by `config`
(objective, constraints, and row bases) and build the reduced problem, without
calling a solver. Runs exactly the validation that
`cs_nctssos(pop, config, symmetry)` runs, through the same backend.

Throws an `ArgumentError` when the symmetry does not hold or when no backend
supports the composition. `config.symmetry` must be `nothing`. `physics` is
forwarded to the backend, which rejects it where it is unsupported.
"""
function check_symmetry(
    pop::OptimizationProblem,
    config::SolverConfig,
    symmetry::Symmetry;
    physics=nothing,
)
    _require_unset_config_symmetry(config)
    if _symmetry_route(symmetry) === :finite
        finite_config = _config_with_symmetry_spec(config, SymmetrySpec(symmetry))
        sparsity = isnothing(physics) ? compute_sparsity(pop, finite_config) :
            compute_sparsity(pop, finite_config, physics)
        _check_symmetry_mvp_support(pop, finite_config, sparsity)
        _, backend_report = moment_relax_symmetric(
            pop,
            sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities,
            finite_config.symmetry,
        )
        return _finite_check_report(symmetry, backend_report)
    end
    _check_continuous_route_support(pop, symmetry)
    reduction = continuous_moment_relaxation(pop, config, symmetry.continuous;
        compact_metric=symmetry.compact_metric, compress=symmetry.compress, physics,
        clifford_generators=_clifford_generators(symmetry))
    return _continuous_check_report(symmetry, reduction)
end

"""
    cs_nctssos(pop, config::SolverConfig, symmetry::Symmetry; kwargs...)

Solve the relaxation of `pop` reduced by `symmetry`, through the backend that
[`Symmetry`](@ref) selects. Returns `(objective, report, result)`:

- `objective`: the solver's optimal value;
- `report`: a [`SymmetryCheckReport`](@ref);
- `result`: the untouched backend result, a `PolyOptResult` on the `:finite`
  route and the named tuple of `continuous_nctssos` on the continuous routes.

Keyword arguments are forwarded to the backend and keep its defaults: the
`:finite` route accepts those of `cs_nctssos(pop, config)` (`dualize=true`, ...),
the continuous routes those of `continuous_nctssos` (`dualize=false`,
`physics`). With `dualize=true` the continuous routes also return
`result.certificate`, the exact-arithmetic check of the solver's dual solution;
it reports the residual of floating-point output instead of chopping it.

`config.symmetry` must be `nothing`.
"""
function cs_nctssos(
    pop::OptimizationProblem,
    config::SolverConfig,
    symmetry::Symmetry;
    kwargs...,
)
    _require_unset_config_symmetry(config)
    if _symmetry_route(symmetry) === :finite
        result = cs_nctssos(pop, _config_with_symmetry_spec(config, SymmetrySpec(symmetry)); kwargs...)
        return (; objective=result.objective, report=_finite_check_report(symmetry, result.symmetry), result)
    end
    _check_continuous_route_support(pop, symmetry)
    _reject_symmetry_keywords(kwargs)
    result = continuous_nctssos(pop, config, symmetry.continuous;
        compact_metric=symmetry.compact_metric, compress=symmetry.compress,
        clifford_generators=_clifford_generators(symmetry), kwargs...)
    return (; objective=result.objective, report=_continuous_check_report(symmetry, result.reduction), result)
end
