# =============================================================================
# Sparsity Results (pre-solve, for debugging/inspection)
# =============================================================================

"""
    SparsityResult{A<:AlgebraType, TI<:Integer, P, M, ST}

Result of sparsity computation for polynomial optimization, before solving.

Contains all information needed to construct and solve the moment relaxation,
useful for debugging and inspecting the problem structure without running the solver.

# Type Parameters
- `A<:AlgebraType`: The algebra type
- `TI<:Integer`: The index type
- `P`: The polynomial type (Polynomial{A,TI} or NCStatePolynomial)
- `M`: The monomial/basis type (Monomial{A,TI} or NCStateWord)
- `ST`: State type - `Nothing` for regular polynomials, `StateType` subtype for state polynomials

# Fields
- `corr_sparsity::CorrelativeSparsity{A,TI,P,M,ST}`: Correlative sparsity structure (cliques, constraints)
- `initial_activated_supps::Vector{Vector{M}}`: Initial activated supports per clique (user-facing, before term sparsity iteration)
- `cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}`: Term sparsity blocks per clique
"""
struct SparsityResult{A<:AlgebraType, TI<:Integer, P, M, ST}
    corr_sparsity::CorrelativeSparsity{A,TI,P,M,ST}
    initial_activated_supps::Vector{Vector{M}}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
end

# =============================================================================
# Polynomial Optimization Result
# =============================================================================

"""
    PolyOptResult{T, A<:AlgebraType, TI<:Integer, P, M, ST}

Result of a polynomial optimization problem solution.

# Type Parameters
- `T`: Coefficient type for the objective value (supports BigFloat for high precision)
- `A`: Algebra type
- `TI`: Index type
- `P`: Polynomial type (Polynomial{A,TI} or NCStatePolynomial)
- `M`: Monomial/basis type (Monomial{A,TI} or NCStateWord)
- `ST`: State type - `Nothing` for regular polynomials, `StateType` subtype for state polynomials

# Fields
- `objective::T`: Optimal objective value
- `sparsity::SparsityResult{A,TI,P,M,ST}`: Sparsity structure (correlative + term sparsity + initial activated supports)
- `model::GenericModel{T}`: JuMP model used for solving
- `moment_matrix_sizes::Vector{Vector{Int}}`: Per-clique vector of term sparsity block sizes for the moment matrix before any symmetry reduction; when `symmetry !== nothing`, the final solver PSD block sizes live in `symmetry.psd_block_sizes`
- `n_unique_moment_matrix_elements::Int`: Number of unique moment variables in all moment matrices (after canonicalization)
- `symmetry::Union{Nothing,SymmetryReport}`: Summary of an applied symmetry reduction, or `nothing` for the ordinary path
"""
struct PolyOptResult{T, A<:AlgebraType, TI<:Integer, P, M, ST}
    objective::T
    sparsity::SparsityResult{A,TI,P,M,ST}
    model::GenericModel{T}
    moment_matrix_sizes::Vector{Vector{Int}}
    n_unique_moment_matrix_elements::Int
    symmetry::Union{Nothing,SymmetryReport}
end

"""
    _compute_moment_matrix_sizes(cliques_term_sparsities)

Compute the block sizes of the moment matrix for each clique.

Returns a `Vector{Vector{Int}}` where each element is the vector of block sizes
for a clique's moment matrix (from `term_sparsity.block_bases` for the first
term sparsity, which corresponds to the moment matrix).
"""
function _compute_moment_matrix_sizes(cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {M}
    map(cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end
end

"""
    PolyOptResult(objective, sparsity, model, n_unique_elements)

Construct a PolyOptResult from a SparsityResult.
"""
function PolyOptResult(
    objective::T,
    sparsity::SparsityResult{A,TI,P,M,ST},
    model::GenericModel{T},
    n_unique_elements::Int;
    symmetry::Union{Nothing,SymmetryReport}=nothing,
) where {T, A<:AlgebraType, TI<:Integer, P, M, ST}
    moment_matrix_sizes = _compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)

    return PolyOptResult{T,A,TI,P,M,ST}(
        objective,
        sparsity,
        model,
        moment_matrix_sizes,
        n_unique_elements,
        symmetry,
    )
end

# Show for regular polynomials (ST = Nothing)
function Base.show(io::IO, result::PolyOptResult{T,A,TI,P,M,Nothing}) where {T,A,TI,P,M}
    _show_poly_opt_result(io, result, "")
end

# Show for state polynomials (ST <: StateType)
function Base.show(io::IO, result::PolyOptResult{T,A,TI,P,M,ST}) where {T,A,TI,P,M,ST<:StateType}
    _show_poly_opt_result(io, result, "State ")
end

function _show_poly_opt_result(io::IO, result::PolyOptResult, prefix::String)
    !isempty(prefix) && println(io, "$(prefix)Optimization Result")
    println(io, "Objective: ", result.objective)
    show(io, result.sparsity.corr_sparsity)
    println(io, "Term Sparsity:")
    for (i, sparsities) in enumerate(result.sparsity.cliques_term_sparsities)
        println(io, "Clique $i:")
        if isnothing(result.symmetry)
            println(io, "   Moment Matrix Block Sizes: ", result.moment_matrix_sizes[i])
        else
            println(io, "   Moment Matrix Block Sizes (pre-symmetry): ", result.moment_matrix_sizes[i])
            i == 1 && println(io, "   Symmetry-Reduced PSD Block Sizes: ", result.symmetry.psd_block_sizes)
        end
        println(io, "   Moment Matrix:")
        println(io, sparsities[1])
        println(io, "   Localizing Matrix:")
        for sparsity in sparsities[2:end]
            show(io, sparsity)
        end
    end
    println(io, "Unique Moment Matrix Elements: ", result.n_unique_moment_matrix_elements)
    if !isnothing(result.symmetry)
        println(io, "Symmetry: ", result.symmetry)
    end
end

"""
    compute_relaxation_order(pop::OptimizationProblem, user_order::Int) -> Int

Compute the relaxation order for the moment-SOS hierarchy.

If `user_order > 0`, returns it directly. Otherwise, computes the minimum order
needed to capture all polynomial degrees: `ceil(max_degree / 2)`.

Returns 1 for trivial problems (empty or zero-degree polynomials).
"""
function compute_relaxation_order(pop::OptimizationProblem, user_order::Int)::Int
    iszero(user_order) || return user_order
    all_polys = [pop.objective; pop.eq_constraints; pop.ineq_constraints; pop.moment_eq_constraints]
    max_deg = maximum(degree(poly) for poly in all_polys)
    isfinite(max_deg) ? max(1, ceil(Int, max_deg / 2)) : 1
end

# Acceptable solver termination statuses
const _ACCEPTABLE_STATUSES = Set([
    MOI.OPTIMAL,
    MOI.ALMOST_OPTIMAL,
    MOI.LOCALLY_SOLVED,
    MOI.ALMOST_LOCALLY_SOLVED,
])

"""
    SolverStatusError

Exception thrown when the solver termination/result statuses are unacceptable.
Stores structured status codes so tests and callers can inspect them without
depending on message formatting.
"""
struct SolverStatusError <: Exception
    termination
    primal
    dual
end

function Base.showerror(io::IO, err::SolverStatusError)
    print(io, "Solver failed: termination=$(err.termination), primal=$(err.primal), dual=$(err.dual)")
end

"""
    _check_solver_status(model)

Check that the solver terminated successfully. Throws an error if the solver
failed with an unacceptable status (infeasible, unbounded, numerical error, etc.).
"""
function _check_solver_status(model)
    status = termination_status(model)
    status ∈ _ACCEPTABLE_STATUSES && return status

    primal = primal_status(model)
    dual = dual_status(model)

    # Some solvers can terminate early (e.g. iteration limit, slow progress) while still
    # returning a usable (nearly) feasible point. Treat those as acceptable so that
    # downstream callers can validate objective values against known tolerances.
    if status in (MOI.ITERATION_LIMIT, MOI.SLOW_PROGRESS) &&
       (primal ∈ (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT) ||
        dual ∈ (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT))
        return status
    end

    throw(SolverStatusError(status, primal, dual))
end

"""
    project_to_clique(poly::Polynomial, clique_indices) -> Polynomial

Project a polynomial to a clique by keeping only terms whose variables are
all contained in the clique.

Returns a polynomial containing only the terms of `poly` where all variable
indices are in `clique_indices`.
"""
function project_to_clique(poly::Polynomial{A,T,C}, clique_indices) where {A,T,C}
    clique_set = Set(clique_indices)
    result_terms = Tuple{C,NormalMonomial{A,T}}[]
    for (coef, mono) in terms(poly)
        if issubset(variable_indices(mono), clique_set)
            push!(result_terms, (coef, mono))
        end
    end
    isempty(result_terms) ? zero(poly) : Polynomial(result_terms)
end

"""
    project_to_clique(poly::NCStatePolynomial, clique_indices) -> NCStatePolynomial

Project a state polynomial to a clique by keeping only terms whose variables are
all contained in the clique.
"""
function project_to_clique(poly::NCStatePolynomial{C,ST,A,T}, clique_indices) where {C,ST,A,T}
    clique_set = Set(clique_indices)
    result_coeffs = C[]
    result_words = NCStateWord{ST,A,T}[]
    for (coef, word) in zip(coefficients(poly), monomials(poly))
        if issubset(variable_indices(word), clique_set)
            push!(result_coeffs, coef)
            push!(result_words, word)
        end
    end
    isempty(result_coeffs) ? zero(poly) : NCStatePolynomial(result_coeffs, result_words)
end

"""
    solve_sdp(moment_problem, optimizer; dualize::Bool=true, sos_hermitian_representation=:real_lift)

Solve the SDP relaxation, either via SOS dualization or directly as moment problem.

For ordinary complex-algebra polynomial problems (`PauliAlgebra`, `FermionicAlgebra`,
`BosonicAlgebra`):
- `dualize=false` solves the primal moment SDP directly, realifying each Hermitian
  cone by the standard block map `[Re(H) -Im(H); Im(H) Re(H)]` (the note's
  block-doubling / Algorithm 1 route).
- `dualize=true` solves the dual SOS SDP for the Hermitian setting covered by the
  note. Hermitian PSD cones use the lifted real PSD dual formulation by default;
  pass `sos_hermitian_representation=:native` to keep native Hermitian cones
  when the solver backend supports them. Complex zero constraints are first
  split into Hermitian real/imaginary components when needed, so non-Hermitian
  equalities still reach the dual path through an explicit Hermitian
  decomposition instead of being modeled implicitly.

Returns a named tuple `(objective, model, n_unique_elements, status)`.

Throws an error if the solver fails (infeasible, unbounded, numerical error).
"""
_solve_sdp_n_unique_elements(moment_problem) =
    moment_problem.n_unique_moment_matrix_elements
_solve_sdp_n_unique_elements(L::MomentLinearData) = length(L.moments)

function solve_sdp(
    moment_problem,
    optimizer;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
    sos_hermitian_representation::Symbol=:real_lift,
)
    if dualize
        if formulation != :moment_variables || representation != :real || orphan_policy != :error
            throw(ArgumentError("Moment lowering options apply only with dualize=false; SOS relowering is deferred."))
        end

        sos_problem = if moment_problem isa Union{MomentProblem,MomentLinearData}
            sos_dualize(
                moment_problem;
                hermitian_representation=sos_hermitian_representation,
            )
        else
            sos_hermitian_representation == :real_lift || throw(ArgumentError(
                "`sos_hermitian_representation` applies only to ordinary Hermitian " *
                "MomentProblem or MomentLinearData SOS dualization."
            ))
            sos_dualize(moment_problem)
        end
        set_optimizer(sos_problem.model, optimizer)
        optimize!(sos_problem.model)
        status = _check_solver_status(sos_problem.model)
        return (
            objective = objective_value(sos_problem.model),
            model = sos_problem.model,
            n_unique_elements = _solve_sdp_n_unique_elements(moment_problem),
            status = status
        )
    else
        result = if moment_problem isa Union{MomentProblem,MomentLinearData}
            solve_moment_problem(
                moment_problem,
                optimizer;
                formulation=formulation,
                representation=representation,
                orphan_policy=orphan_policy,
            )
        else
            if formulation != :moment_variables || representation != :real || orphan_policy != :error
                throw(ArgumentError("State moment lowering does not support formulation/representation options."))
            end
            solve_moment_problem(moment_problem, optimizer)
        end
        status = _check_solver_status(result.model)
        return (
            objective = result.objective,
            model = result.model,
            n_unique_elements = result.n_unique_elements,
            status = status
        )
    end
end

"""
    SolverConfig(; optimizer, order=0, moment_basis=nothing, cs_algo=NoElimination(), ts_algo=NoElimination(), symmetry=nothing)

Configuration for solving polynomial optimization problems.

# Keyword Arguments
- `optimizer` (required): The optimizer to use for solving the SDP problem (e.g. Clarabel.Optimizer)
- `order::Int`: The order of the moment relaxation (default: 0)
- `moment_basis`: Optional custom basis, replacing automatic order-based basis
  generation. For ordinary polynomial problems, pass monomials (or single-term
  unit-coefficient polynomials). For state/trace polynomial problems, pass
  `NCStateWord`s, `StateWord`s, monomials (interpreted as identity-state
  operator words), or single-term unit-coefficient state polynomials. The basis
  must include the identity element. The supplied basis must also generate every
  objective/constraint moment required by the relaxation; underspecified bases
  raise an error instead of silently dropping terms. Default: `nothing`
- `cs_algo::EliminationAlgorithm`: Algorithm for correlative sparsity exploitation (default: NoElimination())
- `ts_algo::EliminationAlgorithm`: Algorithm for term sparsity exploitation (default: NoElimination())
- `symmetry::Union{Nothing,SymmetrySpec}`: Optional symmetry reduction spec. The
  current implementation is still intentionally narrow: dense ordinary polynomial
  relaxations only (`cs_algo=NoElimination()`, `ts_algo=NoElimination()`), with
  supported combinations limited to either
  - `MonoidAlgebra` + `SignedPermutation`,
  - `PauliAlgebra` + `CliffordSymmetry`, optionally with `PauliChargeSectorSpec` and `PauliSingletConstraintSpec`, or
  - `FermionicAlgebra` + `FermionicModePermutation`, `FermionicSectorSpec`, and optional `FermionicSpinAdaptationSpec` layered on sector blocks.
  Unsupported combinations error instead of silently doing the wrong thing.

# Examples
```jldoctest; setup=:(using NCTSSoS, COSMO)
julia> solver_config = SolverConfig(optimizer=COSMO.Optimizer, order=2) # default elimination algorithms
SolverConfig(COSMO.Optimizer, 2, nothing, NoElimination(), NoElimination(), nothing)
```

"""
@kwdef struct SolverConfig
    optimizer
    order::Int = 0
    moment_basis::Union{Nothing,Vector} = nothing
    cs_algo::EliminationAlgorithm = NoElimination()
    ts_algo::EliminationAlgorithm = NoElimination()
    symmetry::Union{Nothing,SymmetrySpec} = nothing
end

function _resolve_relaxation_spec(pop::OptimizationProblem, solver_config::SolverConfig)
    if isnothing(solver_config.moment_basis)
        return compute_relaxation_order(pop, solver_config.order)
    end

    iszero(solver_config.order) ||
        throw(ArgumentError("Specify either `order` or `moment_basis`, not both."))

    return solver_config.moment_basis
end

@inline _has_active_sparsity(solver_config::SolverConfig) =
    !(solver_config.cs_algo isa NoElimination && solver_config.ts_algo isa NoElimination)

function _check_no_pauli_solverconfig_fast_path_overrides(kwargs)
    for key in (:sign_symmetry, :reflection_symmetry, :axis_rotation_symmetry, :check_invariance)
        haskey(kwargs, key) || continue
        throw(ArgumentError(
            "`pauli_translation_invariant_nctssos(pop, solver_config)` derives `$key` from `solver_config.symmetry`; do not pass it manually."
        ))
    end
    return nothing
end

function _pauli_solverconfig_fast_path_supported(profile; su2_symmetry::Bool=false)
    if su2_symmetry && profile.axis_rotation_symmetry
        unsupported = setdiff(profile.unsupported_features, [:axis_rotation])
        return isempty(unsupported) && isempty(profile.missing_required_features)
    end
    return profile.supported_by_translation_fast_path
end

function _check_pauli_solverconfig_fast_path_supported(profile; su2_symmetry::Bool=false)
    _pauli_solverconfig_fast_path_supported(profile; su2_symmetry) && return nothing
    throw(ArgumentError(
        "`solver_config.symmetry` is not supported by the Pauli translation fast path; " *
        "unsupported_features=$(profile.unsupported_features), " *
        "missing_required_features=$(profile.missing_required_features). " *
        "Axis-rotation generators require the full global H/S generator set, " *
        "or `su2_symmetry=true`, where they are subsumed by the SU(2) reducer."
    ))
end

function _check_pauli_complex_reflection_sectors(
    profile;
    real_moment_matrix::Bool=true,
    momenta=nothing,
    kwargs...,
)
    return _check_translation_complex_reflection_sectors(
        profile.n_sites,
        momenta;
        reflection_symmetry=profile.reflection_symmetry,
        real_moment_matrix,
    )
end

function _pauli_translation_fast_path_options(;
    direct_linear::Bool=false,
    momenta=nothing,
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition=nothing,
    contiguous_rdm_support=nothing,
    u1_symmetry::Bool=false,
    su2_symmetry::Bool=false,
    base_su2_extend_rdm::Bool=false,
    su2_moment_quotient::Bool=false,
    su2_moment_quotient_atol::Real=1e-11,
    su2_moment_quotient_condition_limit::Real=1e10,
    qmbcertify_base_construct::Bool=false,
    qmbcertify_base_extra=nothing,
    qmbcertify_base_three_type::Tuple{<:Integer,<:Integer}=(1, 1),
    axis_rotation_equalities::Bool=false,
    axis_rotation_quotient::Bool=false,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
    linear_state_opt_width=nothing,
    linear_state_opt_mode=nothing,
    psd_state_opt_width=nothing,
)
    return (
        direct_linear=direct_linear,
        momenta=momenta,
        real_moment_matrix=real_moment_matrix,
        phase_atol=phase_atol,
        contiguous_rdm_k=contiguous_rdm_k,
        contiguous_rdm_decomposition=contiguous_rdm_decomposition,
        contiguous_rdm_support=contiguous_rdm_support,
        u1_symmetry=u1_symmetry,
        su2_symmetry=su2_symmetry,
        base_su2_extend_rdm=base_su2_extend_rdm,
        su2_moment_quotient=su2_moment_quotient,
        su2_moment_quotient_atol=su2_moment_quotient_atol,
        su2_moment_quotient_condition_limit=su2_moment_quotient_condition_limit,
        qmbcertify_base_construct=qmbcertify_base_construct,
        qmbcertify_base_extra=qmbcertify_base_extra,
        qmbcertify_base_three_type=qmbcertify_base_three_type,
        axis_rotation_equalities=axis_rotation_equalities,
        axis_rotation_quotient=axis_rotation_quotient,
        singlet_channel_equalities=singlet_channel_equalities,
        singlet_channel_atol=singlet_channel_atol,
        linear_state_opt_width=linear_state_opt_width,
        linear_state_opt_mode=linear_state_opt_mode,
        psd_state_opt_width=psd_state_opt_width,
    )
end

function _nondefault_pauli_translation_fast_path_options(options)
    names = Symbol[]
    isnothing(options.momenta) || push!(names, :momenta)
    options.real_moment_matrix || push!(names, :real_moment_matrix)
    options.phase_atol == 1e-12 || push!(names, :phase_atol)
    isnothing(options.contiguous_rdm_k) || push!(names, :contiguous_rdm_k)
    isnothing(options.contiguous_rdm_decomposition) || push!(names, :contiguous_rdm_decomposition)
    isnothing(options.contiguous_rdm_support) || push!(names, :contiguous_rdm_support)
    options.u1_symmetry && push!(names, :u1_symmetry)
    options.su2_symmetry && push!(names, :su2_symmetry)
    options.base_su2_extend_rdm && push!(names, :base_su2_extend_rdm)
    options.su2_moment_quotient && push!(names, :su2_moment_quotient)
    options.su2_moment_quotient_atol == 1e-11 ||
        push!(names, :su2_moment_quotient_atol)
    options.su2_moment_quotient_condition_limit == 1e10 ||
        push!(names, :su2_moment_quotient_condition_limit)
    options.qmbcertify_base_construct && push!(names, :qmbcertify_base_construct)
    isnothing(options.qmbcertify_base_extra) || push!(names, :qmbcertify_base_extra)
    options.qmbcertify_base_three_type == (1, 1) || push!(names, :qmbcertify_base_three_type)
    options.axis_rotation_equalities && push!(names, :axis_rotation_equalities)
    options.axis_rotation_quotient && push!(names, :axis_rotation_quotient)
    options.singlet_channel_equalities && push!(names, :singlet_channel_equalities)
    options.singlet_channel_atol == 1e-12 || push!(names, :singlet_channel_atol)
    isnothing(options.linear_state_opt_width) || push!(names, :linear_state_opt_width)
    isnothing(options.linear_state_opt_mode) || push!(names, :linear_state_opt_mode)
    isnothing(options.psd_state_opt_width) || push!(names, :psd_state_opt_width)
    return names
end

function _check_qmbcertify_base_bridge_options(
    ;
    momenta,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    su2_moment_quotient::Bool,
    su2_moment_quotient_atol::Real,
    su2_moment_quotient_condition_limit::Real,
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    singlet_channel_equalities::Bool,
    singlet_channel_atol::Real,
)
    unsupported = Symbol[]
    isnothing(momenta) || push!(unsupported, :momenta)
    u1_symmetry && push!(unsupported, :u1_symmetry)
    su2_symmetry && push!(unsupported, :su2_symmetry)
    su2_moment_quotient && push!(unsupported, :su2_moment_quotient)
    su2_moment_quotient_atol == 1e-11 ||
        push!(unsupported, :su2_moment_quotient_atol)
    su2_moment_quotient_condition_limit == 1e10 ||
        push!(unsupported, :su2_moment_quotient_condition_limit)
    axis_rotation_equalities && push!(unsupported, :axis_rotation_equalities)
    axis_rotation_quotient && push!(unsupported, :axis_rotation_quotient)
    singlet_channel_equalities && push!(unsupported, :singlet_channel_equalities)
    singlet_channel_atol == 1e-12 || push!(unsupported, :singlet_channel_atol)
    isempty(unsupported) && return nothing
    keys = join(("`$name`" for name in unsupported), ", ")
    throw(ArgumentError(
        "`qmbcertify_base_construct=true` does not support Pauli translation fast-path keyword(s) $keys."
    ))
end

function _check_no_pauli_translation_fast_path_options(options; context::AbstractString)
    names = _nondefault_pauli_translation_fast_path_options(options)
    isempty(names) && return nothing
    keys = join(("`$name`" for name in names), ", ")
    throw(ArgumentError(
        "Pauli translation fast-path keyword(s) $keys require an ordinary Pauli-chain problem routed through the translation fast path; $context."
    ))
end

function _supports_no_symmetry_direct_linear(::OptimizationProblem{A,P}) where {A<:AlgebraType,P}
    return P <: Polynomial && (A <: MonoidAlgebra || A === PauliAlgebra || A <: PBWAlgebra)
end

function _check_no_symmetry_direct_linear_support(pop::OptimizationProblem{A,P}) where {A<:AlgebraType,P}
    P <: Polynomial || throw(ArgumentError(
        "`direct_linear=true` without a supported Pauli translation symmetry currently supports ordinary polynomial optimization only; state/trace polynomial optimization is not yet supported."
    ))
    _supports_no_symmetry_direct_linear(pop) || throw(ArgumentError(
        "`direct_linear=true` without a supported Pauli translation symmetry currently supports ordinary no-symmetry polynomial problems over `MonoidAlgebra`, `PauliAlgebra`, or `PBWAlgebra`; got `$(nameof(A))`."
    ))
    return nothing
end

function _is_trivial_finite_symmetry(symmetry::SymmetrySpec)
    isnothing(symmetry.sector) || return false
    isnothing(symmetry.spin_adaptation) || return false
    isnothing(symmetry.pauli_charge) || return false
    isnothing(symmetry.pauli_singlet) || return false
    all(generator -> isempty(generator.images), symmetry.generators) || return false
    all(generator -> isempty(generator.images), symmetry.fermionic_generators) || return false
    all(generator -> isempty(generator.images), symmetry.clifford_generators) || return false
    return true
end

function _use_no_symmetry_direct_linear(
    pop::OptimizationProblem,
    solver_config::SolverConfig;
    direct_linear::Bool,
    trivial_finite_symmetry::Bool,
)
    no_active_symmetry = isnothing(solver_config.symmetry) || trivial_finite_symmetry
    if direct_linear
        no_active_symmetry || throw(ArgumentError(
            "`direct_linear=true` with `solver_config.symmetry` is supported only by the Pauli translation fast path."
        ))
        _check_no_symmetry_direct_linear_support(pop)
        return true
    end
    no_active_symmetry || return false
    return _supports_no_symmetry_direct_linear(pop)
end

function pauli_translation_invariant_nctssos(
    pop::PolyOpt{PauliAlgebra,T,P},
    solver_config::SolverConfig;
    dualize::Bool=false,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
    sos_hermitian_representation::Symbol=:real_lift,
    kwargs...,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    isnothing(solver_config.symmetry) && throw(ArgumentError(
        "`pauli_translation_invariant_nctssos(pop, solver_config)` requires `solver_config.symmetry`."
    ))
    _has_active_sparsity(solver_config) && throw(ArgumentError(
        "`pauli_translation_invariant_nctssos(pop, solver_config)` does not use generic CS/TS sparsity; pass `NoElimination()` algorithms."
    ))
    isnothing(solver_config.moment_basis) || throw(ArgumentError(
        "`pauli_translation_invariant_nctssos(pop, solver_config)` currently requires an integer relaxation `order`, not `moment_basis`."
    ))
    _check_no_pauli_solverconfig_fast_path_overrides(kwargs)

    order = _resolve_relaxation_spec(pop, solver_config)
    order isa Integer || throw(ArgumentError(
        "`pauli_translation_invariant_nctssos(pop, solver_config)` requires an integer relaxation order."
    ))
    ops = _pauli_chain_ops_from_registry(pop.registry)
    profile = pauli_chain_fast_path_profile(ops, solver_config.symmetry)
    _check_pauli_solverconfig_fast_path_supported(
        profile;
        su2_symmetry=get(kwargs, :su2_symmetry, false),
    )
    _check_pauli_complex_reflection_sectors(profile; kwargs...)

    if get(kwargs, :qmbcertify_base_construct, false)
        return pauli_translation_invariant_nctssos(
            pop,
            ops,
            Int(order),
            solver_config.optimizer;
            dualize,
            formulation,
            representation,
            orphan_policy,
            sos_hermitian_representation,
            kwargs...,
        )
    end

    return pauli_translation_invariant_nctssos(
        pop,
        ops,
        Int(order),
        solver_config.optimizer;
        dualize,
        formulation,
        representation,
        orphan_policy,
        sos_hermitian_representation,
        sign_symmetry=profile.sign_symmetry,
        reflection_symmetry=profile.reflection_symmetry,
        axis_rotation_symmetry=profile.axis_rotation_symmetry,
        check_invariance=solver_config.symmetry.check_invariance,
        kwargs...,
    )
end

function _maybe_pauli_translation_fast_path(
    pop,
    solver_config::SolverConfig;
    dualize::Bool,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
    sos_hermitian_representation::Symbol,
    kwargs...,
)
    return nothing
end

function _maybe_pauli_translation_fast_path(
    pop::PolyOpt{PauliAlgebra,T,P},
    solver_config::SolverConfig;
    dualize::Bool,
    formulation::Symbol,
    representation::Symbol,
    orphan_policy::Symbol,
    sos_hermitian_representation::Symbol,
    direct_linear::Bool,
    momenta,
    real_moment_matrix::Bool,
    phase_atol::Real,
    contiguous_rdm_k,
    contiguous_rdm_decomposition,
    contiguous_rdm_support,
    u1_symmetry::Bool,
    su2_symmetry::Bool,
    base_su2_extend_rdm::Bool,
    su2_moment_quotient::Bool,
    su2_moment_quotient_atol::Real,
    su2_moment_quotient_condition_limit::Real,
    qmbcertify_base_construct::Bool,
    qmbcertify_base_extra,
    qmbcertify_base_three_type::Tuple{<:Integer,<:Integer},
    axis_rotation_equalities::Bool,
    axis_rotation_quotient::Bool,
    singlet_channel_equalities::Bool,
    singlet_channel_atol::Real,
    linear_state_opt_width,
    linear_state_opt_mode,
    psd_state_opt_width,
) where {T<:Unsigned,C<:Number,P<:Polynomial{PauliAlgebra,T,C}}
    isnothing(solver_config.symmetry) && return nothing
    _has_active_sparsity(solver_config) && return nothing
    isnothing(solver_config.moment_basis) || return nothing

    order = _resolve_relaxation_spec(pop, solver_config)
    order isa Integer || return nothing
    ops = try
        _pauli_chain_ops_from_registry(pop.registry)
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    profile = pauli_chain_fast_path_profile(ops, solver_config.symmetry)
    _pauli_solverconfig_fast_path_supported(profile; su2_symmetry) || return nothing
    _check_pauli_complex_reflection_sectors(profile; real_moment_matrix, momenta)

    if qmbcertify_base_construct
        _check_qmbcertify_base_bridge_options(
            ;
            momenta,
            u1_symmetry,
            su2_symmetry,
            su2_moment_quotient,
            su2_moment_quotient_atol,
            su2_moment_quotient_condition_limit,
            axis_rotation_equalities,
            axis_rotation_quotient,
            singlet_channel_equalities,
            singlet_channel_atol,
        )
        return pauli_translation_invariant_nctssos(
            pop,
            ops,
            Int(order),
            solver_config.optimizer;
            dualize,
            formulation,
            representation,
            orphan_policy,
            sos_hermitian_representation,
            direct_linear,
            base_su2_extend_rdm,
            su2_moment_quotient,
            su2_moment_quotient_atol,
            su2_moment_quotient_condition_limit,
            qmbcertify_base_construct,
            qmbcertify_base_extra,
            qmbcertify_base_three_type,
            real_moment_matrix,
            phase_atol,
            contiguous_rdm_k,
            contiguous_rdm_decomposition,
            contiguous_rdm_support,
            linear_state_opt_width,
            linear_state_opt_mode,
            psd_state_opt_width,
        )
    end

    return pauli_translation_invariant_nctssos(
        pop,
        ops,
        Int(order),
        solver_config.optimizer;
        dualize,
        formulation,
        representation,
        orphan_policy,
        sos_hermitian_representation,
        sign_symmetry=profile.sign_symmetry,
        reflection_symmetry=profile.reflection_symmetry,
        axis_rotation_symmetry=profile.axis_rotation_symmetry,
        check_invariance=solver_config.symmetry.check_invariance,
        direct_linear,
        momenta,
        real_moment_matrix,
        phase_atol,
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        contiguous_rdm_support,
        u1_symmetry,
        su2_symmetry,
        base_su2_extend_rdm,
        su2_moment_quotient,
        su2_moment_quotient_atol,
        su2_moment_quotient_condition_limit,
        qmbcertify_base_construct,
        qmbcertify_base_extra,
        qmbcertify_base_three_type,
        axis_rotation_equalities,
        axis_rotation_quotient,
        singlet_channel_equalities,
        singlet_channel_atol,
        linear_state_opt_width,
        linear_state_opt_mode,
        psd_state_opt_width,
    )
end

function _check_symmetry_mvp_support(
    pop::OptimizationProblem{A,P},
    solver_config::SolverConfig,
    sparsity::SparsityResult,
) where {A<:AlgebraType,P}
    isnothing(solver_config.symmetry) && return nothing

    symmetry = solver_config.symmetry

    P <: Polynomial || throw(ArgumentError(
        "Symmetry reduction currently supports ordinary polynomial problems only; state/trace problems are not yet supported."
    ))
    solver_config.cs_algo isa NoElimination || throw(ArgumentError(
        "Symmetry reduction currently requires `cs_algo=NoElimination()`."
    ))
    solver_config.ts_algo isa NoElimination || throw(ArgumentError(
        "Symmetry reduction currently requires `ts_algo=NoElimination()`."
    ))
    length(sparsity.corr_sparsity.cliques) == 1 || throw(ArgumentError(
        "Symmetry reduction currently supports a single dense clique only."
    ))

    for term_sparsities in sparsity.cliques_term_sparsities, term_sparsity in term_sparsities
        length(term_sparsity.block_bases) == 1 || throw(ArgumentError(
            "Symmetry reduction does not yet compose with term-sparsity block splitting."
        ))
    end

    if !isempty(symmetry.generators)
        A <: MonoidAlgebra || throw(ArgumentError(
            "`SignedPermutation` symmetry is currently supported only for ordinary polynomial problems over `MonoidAlgebra`. Got `$(nameof(A))`."
        ))
        (isnothing(symmetry.sector) && isnothing(symmetry.spin_adaptation)) || throw(ArgumentError(
            "Fermionic sector splitting / spin adaptation cannot be combined with raw `SignedPermutation` symmetry. Use `FermionicModePermutation` for fermionic problems."
        ))
    end

    if !isempty(symmetry.clifford_generators) || !isnothing(symmetry.pauli_charge) || !isnothing(symmetry.pauli_singlet)
        A === PauliAlgebra || throw(ArgumentError(
            "`CliffordSymmetry` / Pauli charge and singlet reductions are currently supported only for ordinary polynomial problems over `PauliAlgebra`. Got `$(nameof(A))`."
        ))
        (isnothing(symmetry.sector) && isnothing(symmetry.spin_adaptation)) || throw(ArgumentError(
            "Fermionic sector splitting / spin adaptation cannot be combined with `CliffordSymmetry` or Pauli charge/singlet symmetry."
        ))
    end

    if !isempty(symmetry.fermionic_generators) || !isnothing(symmetry.sector) || !isnothing(symmetry.spin_adaptation)
        A === FermionicAlgebra || throw(ArgumentError(
            "Fermionic mode permutations / sector splitting / spin adaptation are currently supported only for `FermionicAlgebra`. Got `$(nameof(A))`."
        ))
    end

    if !isnothing(symmetry.spin_adaptation)
        isnothing(symmetry.sector) && throw(ArgumentError(
            "Fermionic spin adaptation currently requires `sector=FermionicSectorSpec(..., split_spin=true)` in the same `SymmetrySpec`."
        ))
    end

    return nothing
end

function _check_symmetry_mvp_support(
    pop::PolyOpt{A,T,P},
    solver_config::SolverConfig,
    _sparsity::SparsityResult,
) where {A<:MonoidAlgebra,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    isnothing(solver_config.symmetry) && return nothing
    throw(ArgumentError(
        "Symmetry reduction MVP does not yet support state/trace polynomial optimization."
    ))
end

"""
    compute_sparsity(pop::PolyOpt, solver_config::SolverConfig) -> SparsityResult

Compute correlative and term sparsity for a polynomial optimization problem.

This function performs all sparsity computations without solving the SDP,
useful for debugging and inspecting the problem structure.

# Returns
- `SparsityResult`: Contains correlative sparsity, initial activated supports, and term sparsities
"""
function compute_sparsity(pop::OP, solver_config::SolverConfig) where {A<:AlgebraType, P, OP<:OptimizationProblem{A,P}}
    relaxation_spec = _resolve_relaxation_spec(pop, solver_config)
    corr_sparsity = correlative_sparsity(pop, relaxation_spec, solver_config.cs_algo)

    # Compute partial objectives for each clique
    cliques_objective = map(c -> project_to_clique(pop.objective, c), corr_sparsity.cliques)

    initial_activated_supps_nm = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps_nm, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    if !isnothing(solver_config.moment_basis)
        total_basis, _, _ = _polynomial_total_basis(pop, corr_sparsity, cliques_term_sparsities)
        _validate_polynomial_relaxation_support(pop, total_basis; source="The supplied `moment_basis`")
    end

    return SparsityResult(corr_sparsity, initial_activated_supps_nm, cliques_term_sparsities)
end

"""
    compute_sparsity(pop::PolyOpt, solver_config::SolverConfig) -> SparsityResult

Compute correlative and term sparsity for a state polynomial optimization problem.
"""
function compute_sparsity(pop::PolyOpt{A,T,P}, solver_config::SolverConfig) where {A<:AlgebraType,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    relaxation_spec = _resolve_relaxation_spec(pop, solver_config)

    # Sparse pure-trace benchmarks in NCTSSOS use the ordinary nc-word basis,
    # not the general state-word basis. Reuse that basis here when sparsity is
    # active so CS/TS block counts match the reference trace hierarchy.
    if relaxation_spec isa Int &&
       ST == MaxEntangled &&
       _has_active_sparsity(solver_config) &&
       _uses_scalar_trace_word_basis(pop)
        relaxation_spec = _embedded_trace_moment_basis(pop.registry, relaxation_spec)
    end

    corr_sparsity = correlative_sparsity(pop, relaxation_spec, solver_config.cs_algo)

    # Compute partial objectives for each clique
    cliques_objective = map(c -> project_to_clique(pop.objective, c), corr_sparsity.cliques)

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    if !isnothing(solver_config.moment_basis)
        total_basis = _state_total_basis(pop, corr_sparsity, cliques_term_sparsities)
        _validate_state_relaxation_support(pop, total_basis; source="The supplied `moment_basis`")
    end

    return SparsityResult(corr_sparsity, initial_activated_supps, cliques_term_sparsities)
end

"""
    cs_nctssos(pop::PolyOpt{P}, solver_config::SolverConfig; dualize::Bool=true) where {P}

Solve a polynomial optimization problem using the CS-NCTSSOS method with correlative sparsity and term sparsity exploitation.

# Arguments
- `pop::PolyOpt{P}`: The polynomial optimization problem to solve
- `solver_config::SolverConfig`: Configuration containing optimizer, moment order, and sparsity algorithms

# Keyword Arguments
- `dualize::Bool=true`: Whether to dualize the moment relaxation to a sum-of-squares problem
- `sos_hermitian_representation::Symbol=:real_lift`: Hermitian SOS dual cone representation for ordinary complex moment problems; use `:native` to keep native Hermitian cones
- `direct_linear::Bool=false`: Explicitly request directly emitted `MomentLinearData`; supported ordinary no-symmetry polynomial problems (`MonoidAlgebra`, `PauliAlgebra`, and `PBWAlgebra`) select this path automatically, and Pauli-chain translation fast paths still use it when requested
- Pauli translation fast-path keywords: `momenta`, `real_moment_matrix`, `phase_atol`, `contiguous_rdm_k`, `contiguous_rdm_decomposition`, `contiguous_rdm_support`, `u1_symmetry`, `su2_symmetry`, `base_su2_extend_rdm`, `su2_moment_quotient`, `su2_moment_quotient_atol`, `su2_moment_quotient_condition_limit`, `qmbcertify_base_construct`, `qmbcertify_base_extra`, `qmbcertify_base_three_type`, `axis_rotation_equalities`, `axis_rotation_quotient`, `singlet_channel_equalities`, `singlet_channel_atol`, `linear_state_opt_width`, `linear_state_opt_mode`, and `psd_state_opt_width`

# Returns
- `PolyOptResult`: Result containing the objective value, correlative sparsity structure, and term sparsity information

# Description
This function solves a polynomial optimization problem by:
1. Computing correlative sparsity to decompose the problem into smaller cliques
2. Computing term sparsity for each clique to further reduce problem size
3. Formulating and solving either the moment relaxation or its SOS dual
4. Returning the optimal objective value and sparsity information

The moment order is automatically determined from the polynomial degrees if not specified in `solver_config`.
"""
function cs_nctssos(
    pop::OP,
    solver_config::SolverConfig;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
    sos_hermitian_representation::Symbol=:real_lift,
    direct_linear::Bool=false,
    momenta=nothing,
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition=nothing,
    contiguous_rdm_support=nothing,
    u1_symmetry::Bool=false,
    su2_symmetry::Bool=false,
    base_su2_extend_rdm::Bool=false,
    su2_moment_quotient::Bool=false,
    su2_moment_quotient_atol::Real=1e-11,
    su2_moment_quotient_condition_limit::Real=1e10,
    qmbcertify_base_construct::Bool=false,
    qmbcertify_base_extra=nothing,
    qmbcertify_base_three_type::Tuple{<:Integer,<:Integer}=(1, 1),
    axis_rotation_equalities::Bool=false,
    axis_rotation_quotient::Bool=false,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
    linear_state_opt_width=nothing,
    linear_state_opt_mode=nothing,
    psd_state_opt_width=nothing,
) where {A<:AlgebraType, P, OP<:OptimizationProblem{A,P}}
    pauli_fast_path_options = _pauli_translation_fast_path_options(
        ;
        direct_linear,
        momenta,
        real_moment_matrix,
        phase_atol,
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        contiguous_rdm_support,
        u1_symmetry,
        su2_symmetry,
        base_su2_extend_rdm,
        su2_moment_quotient,
        su2_moment_quotient_atol,
        su2_moment_quotient_condition_limit,
        qmbcertify_base_construct,
        qmbcertify_base_extra,
        qmbcertify_base_three_type,
        axis_rotation_equalities,
        axis_rotation_quotient,
        singlet_channel_equalities,
        singlet_channel_atol,
        linear_state_opt_width,
        linear_state_opt_mode,
        psd_state_opt_width,
    )
    fast_path_result = _maybe_pauli_translation_fast_path(
        pop,
        solver_config;
        dualize,
        formulation,
        representation,
        orphan_policy,
        sos_hermitian_representation,
        pauli_fast_path_options...,
    )
    isnothing(fast_path_result) || return fast_path_result
    _check_no_pauli_translation_fast_path_options(pauli_fast_path_options; context=(
        "`cs_nctssos` did not match an ordinary Pauli-chain `SolverConfig` supported by the translation backend"
    ))
    trivial_finite_symmetry =
        !isnothing(solver_config.symmetry) &&
        _is_trivial_finite_symmetry(solver_config.symmetry)
    use_direct_linear_no_symmetry = _use_no_symmetry_direct_linear(
        pop,
        solver_config;
        direct_linear,
        trivial_finite_symmetry,
    )

    sparsity = compute_sparsity(pop, solver_config)
    if !trivial_finite_symmetry
        _check_symmetry_mvp_support(pop, solver_config, sparsity)
    end

    if isnothing(solver_config.symmetry) || trivial_finite_symmetry
        moment_problem = if use_direct_linear_no_symmetry
            _moment_relax_linear(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        else
            moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        end
        result = solve_sdp(
            moment_problem,
            solver_config.optimizer;
            dualize=dualize,
            formulation=formulation,
            representation=representation,
            orphan_policy=orphan_policy,
            sos_hermitian_representation=sos_hermitian_representation,
        )
        n_unique_elements = use_direct_linear_no_symmetry ?
            _moment_matrix_element_count(
                A,
                _moment_matrix_basis(sparsity.cliques_term_sparsities),
            ) :
            result.n_unique_elements
        return PolyOptResult(result.objective, sparsity, result.model, n_unique_elements)
    end

    moment_problem, symmetry_report = moment_relax_symmetric(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
        solver_config.symmetry,
    )
    result = solve_sdp(
        moment_problem,
        solver_config.optimizer;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
        sos_hermitian_representation=sos_hermitian_representation,
    )
    return PolyOptResult(
        result.objective,
        sparsity,
        result.model,
        result.n_unique_elements;
        symmetry=symmetry_report,
    )
end

function _higher_step_graph_support(
    activated_supp::Vector{M}, poly::P, basis::Vector{M}
) where {A<:AlgebraType, T<:Integer, C<:Number, P<:Polynomial{A,T,C}, M<:NormalMonomial{A,T}}
    graph = get_term_sparsity_graph(last.(poly.terms), activated_supp, basis)
    return term_sparsity_graph_supp(graph, basis, poly)
end

function _higher_step_graph_support(
    activated_supp::Vector{M}, poly::P, basis::Vector{M}
) where {ST<:StateType, A<:AlgebraType, T<:Integer, C<:Number, P<:NCStatePolynomial{C,ST,A,T}, M<:NCStateWord{ST,A,T}}
    graph = get_term_sparsity_graph(monomials(poly), activated_supp, basis)
    return term_sparsity_graph_supp(graph, basis, poly)
end

"""
    cs_nctssos_higher(pop::PolyOpt{T}, prev_res::PolyOptResult, solver_config::SolverConfig; dualize::Bool=true) where {T}

Solve a polynomial optimization problem using another raw term-sparsity iteration based on a previous result.

# Arguments
- `pop::PolyOpt{T}`: The polynomial optimization problem to solve
- `prev_res::PolyOptResult`: Previous optimization result containing sparsity information to build upon
- `solver_config::SolverConfig`: Configuration containing optimizer and sparsity algorithms

# Keyword Arguments
- `dualize::Bool=true`: Whether to dualize the moment relaxation to a sum-of-squares problem
- `sos_hermitian_representation::Symbol=:real_lift`: Hermitian SOS dual cone representation for ordinary complex moment problems
- `direct_linear::Bool=false`: Explicitly request directly emitted `MomentLinearData` for supported ordinary no-symmetry polynomial problems; these problems also select this path automatically

# Returns
- `PolyOptResult`: Result containing the objective value, correlative sparsity structure, and updated term sparsity information

# Description
For supported ordinary no-symmetry polynomial problems, the updated relaxation
is emitted as `MomentLinearData` directly, matching the default `cs_nctssos`
route.

This function performs another term-sparsity iteration of the CS-NCTSSOS method by:
1. Reusing the correlative sparsity structure and relaxation order from the previous result
2. Recomputing the next activated supports from the previous iteration's raw term-sparsity graphs
3. Formulating and solving either the moment relaxation or its SOS dual with the updated sparsity
4. Returning the optimal objective value and updated sparsity information

Importantly, this does **not** increase the relaxation order or switch to a denser basis. It only updates the raw term-sparsity support at the same order. If that raw support is already at a fixed point, `cs_nctssos_higher` may return exactly the same relaxation and objective value as the previous solve. To obtain a tighter bound in that situation, increase `order` or use a less aggressive sparsity configuration.
"""
function cs_nctssos_higher(
    pop::OP,
    prev_res::PolyOptResult,
    solver_config::SolverConfig;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
    sos_hermitian_representation::Symbol=:real_lift,
    direct_linear::Bool=false,
) where {A<:AlgebraType, P, OP<:OptimizationProblem{A,P}}
    isnothing(solver_config.moment_basis) ||
        throw(ArgumentError("`cs_nctssos_higher` reuses the basis from `prev_res`; do not pass `moment_basis`."))
    trivial_finite_symmetry =
        !isnothing(solver_config.symmetry) &&
        _is_trivial_finite_symmetry(solver_config.symmetry)
    (isnothing(solver_config.symmetry) || trivial_finite_symmetry) ||
        throw(ArgumentError("`cs_nctssos_higher` does not yet accept nontrivial `symmetry` in `solver_config`."))
    isnothing(prev_res.symmetry) ||
        throw(ArgumentError(
            "`cs_nctssos_higher` cannot continue from a symmetry-reduced `prev_res`; rerun from `cs_nctssos` without symmetry."
        ))
    use_direct_linear_no_symmetry = _use_no_symmetry_direct_linear(
        pop,
        solver_config;
        direct_linear,
        trivial_finite_symmetry,
    )

    prev_sparsity = prev_res.sparsity
    prev_corr_sparsity = prev_sparsity.corr_sparsity

    # Recompute the previous raw graph supports before clique fill.
    # Using the stored post-fill support here can widen the next iteration too aggressively.
    initial_activated_supps_nm = map(enumerate(prev_sparsity.initial_activated_supps)) do (clique_idx, activated_supp)
        cons_idx = prev_corr_sparsity.clq_cons[clique_idx]
        clique_cons = prev_corr_sparsity.cons[cons_idx]
        clique_mom_basis = prev_corr_sparsity.clq_mom_mtx_bases[clique_idx]
        clique_localizing_bases = prev_corr_sparsity.clq_localizing_mtx_bases[clique_idx]

        supports = Vector{eltype(activated_supp)}[
            _higher_step_graph_support(activated_supp, one(eltype(clique_cons)), clique_mom_basis),
        ]
        for (poly, basis) in zip(clique_cons, clique_localizing_bases)
            push!(supports, _higher_step_graph_support(activated_supp, poly, basis))
        end
        reduce(sorted_union, supports)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps_nm, prev_corr_sparsity.clq_cons, prev_corr_sparsity.clq_mom_mtx_bases, prev_corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, prev_corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    # Create new sparsity result with updated term sparsities
    sparsity = SparsityResult(prev_corr_sparsity, initial_activated_supps_nm, cliques_term_sparsities)

    moment_problem = if use_direct_linear_no_symmetry
        _moment_relax_linear(pop, prev_corr_sparsity, cliques_term_sparsities)
    else
        moment_relax(pop, prev_corr_sparsity, cliques_term_sparsities)
    end

    result = solve_sdp(
        moment_problem,
        solver_config.optimizer;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
        sos_hermitian_representation=sos_hermitian_representation,
    )
    n_unique_elements = use_direct_linear_no_symmetry ?
        _moment_matrix_element_count(
            A,
            _moment_matrix_basis(sparsity.cliques_term_sparsities),
        ) :
        result.n_unique_elements
    return PolyOptResult(result.objective, sparsity, result.model, n_unique_elements)
end


# =============================================================================
# State Polynomial Optimization Solver
# =============================================================================

"""
    cs_nctssos(pop::PolyOpt{A,T,P}, solver_config::SolverConfig; dualize::Bool=true)

Solve a state polynomial optimization problem using the CS-NCTSSOS method.

# Arguments
- `pop::PolyOpt{A,T,P}`: The polynomial optimization problem with NCStatePolynomial objective
- `solver_config::SolverConfig`: Configuration containing optimizer, moment order, and sparsity algorithms

# Keyword Arguments
- `dualize::Bool=true`: Whether to dualize the moment relaxation to a sum-of-squares problem

# Returns
- `PolyOptResult`: Result containing the objective value and sparsity information

# Description
This function solves a state polynomial optimization problem by:
1. Computing correlative sparsity to decompose the problem into smaller cliques
2. Computing term sparsity for each clique to further reduce problem size
3. Formulating and solving the SOS dual of the moment relaxation
4. Returning the optimal objective value and sparsity information

# Note on order selection
The relaxation order should be at least `ceil(max_objective_degree / 2)` to properly
capture the objective. For state polynomials with degree-2 terms like ⟨x₁y₁⟩,
use `order >= 1`. If `order=0` is specified, it will be automatically computed
from the maximum polynomial degree.
"""
function cs_nctssos(
    pop::PolyOpt{A,T,P},
    solver_config::SolverConfig;
    dualize::Bool=true,
    formulation::Symbol=:moment_variables,
    representation::Symbol=:real,
    orphan_policy::Symbol=:error,
    sos_hermitian_representation::Symbol=:real_lift,
    direct_linear::Bool=false,
    momenta=nothing,
    real_moment_matrix::Bool=true,
    phase_atol::Real=1e-12,
    contiguous_rdm_k=nothing,
    contiguous_rdm_decomposition=nothing,
    contiguous_rdm_support=nothing,
    u1_symmetry::Bool=false,
    su2_symmetry::Bool=false,
    base_su2_extend_rdm::Bool=false,
    su2_moment_quotient::Bool=false,
    su2_moment_quotient_atol::Real=1e-11,
    su2_moment_quotient_condition_limit::Real=1e10,
    qmbcertify_base_construct::Bool=false,
    qmbcertify_base_extra=nothing,
    qmbcertify_base_three_type::Tuple{<:Integer,<:Integer}=(1, 1),
    axis_rotation_equalities::Bool=false,
    axis_rotation_quotient::Bool=false,
    singlet_channel_equalities::Bool=false,
    singlet_channel_atol::Real=1e-12,
    linear_state_opt_width=nothing,
    linear_state_opt_mode=nothing,
    psd_state_opt_width=nothing,
) where {A<:AlgebraType,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    pauli_fast_path_options = _pauli_translation_fast_path_options(
        ;
        direct_linear,
        momenta,
        real_moment_matrix,
        phase_atol,
        contiguous_rdm_k,
        contiguous_rdm_decomposition,
        contiguous_rdm_support,
        u1_symmetry,
        su2_symmetry,
        base_su2_extend_rdm,
        su2_moment_quotient,
        su2_moment_quotient_atol,
        su2_moment_quotient_condition_limit,
        qmbcertify_base_construct,
        qmbcertify_base_extra,
        qmbcertify_base_three_type,
        axis_rotation_equalities,
        axis_rotation_quotient,
        singlet_channel_equalities,
        singlet_channel_atol,
        linear_state_opt_width,
        linear_state_opt_mode,
        psd_state_opt_width,
    )
    if direct_linear
        _check_no_symmetry_direct_linear_support(pop)
    end
    _check_no_pauli_translation_fast_path_options(pauli_fast_path_options; context=(
        "state/trace polynomial optimization cannot use them"
    ))
    isnothing(solver_config.symmetry) || throw(ArgumentError(
        "Symmetry reduction MVP does not yet support state/trace polynomial optimization."
    ))
    sparsity = compute_sparsity(pop, solver_config)
    moment_problem = moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    result = solve_sdp(
        moment_problem,
        solver_config.optimizer;
        dualize=dualize,
        formulation=formulation,
        representation=representation,
        orphan_policy=orphan_policy,
        sos_hermitian_representation=sos_hermitian_representation,
    )
    return PolyOptResult(result.objective, sparsity, result.model, result.n_unique_elements)
end
