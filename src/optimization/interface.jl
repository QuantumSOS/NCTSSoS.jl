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
- `moment_matrix_sizes::Vector{Vector{Int}}`: Per-clique vector of term sparsity block sizes for the moment matrix
- `n_unique_moment_matrix_elements::Int`: Number of unique moment variables in all moment matrices (after canonicalization)
"""
struct PolyOptResult{T, A<:AlgebraType, TI<:Integer, P, M, ST}
    objective::T
    sparsity::SparsityResult{A,TI,P,M,ST}
    model::GenericModel{T}
    moment_matrix_sizes::Vector{Vector{Int}}
    n_unique_moment_matrix_elements::Int
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
    n_unique_elements::Int
) where {T, A<:AlgebraType, TI<:Integer, P, M, ST}
    moment_matrix_sizes = _compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)

    return PolyOptResult{T,A,TI,P,M,ST}(
        objective,
        sparsity,
        model,
        moment_matrix_sizes,
        n_unique_elements
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
        println(io, "   Moment Matrix Block Sizes: ", result.moment_matrix_sizes[i])
        println(io, "   Moment Matrix:")
        println(io, sparsities[1])
        println(io, "   Localizing Matrix:")
        for sparsity in sparsities[2:end]
            show(io, sparsity)
        end
    end
    println(io, "Unique Moment Matrix Elements: ", result.n_unique_moment_matrix_elements)
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
    all_polys = [pop.objective; pop.eq_constraints; pop.ineq_constraints]
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

    error("Solver failed: termination=$status, primal=$primal, dual=$dual")
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
    solve_sdp(moment_problem, optimizer; dualize::Bool=true)

Solve the SDP relaxation, either via SOS dualization or directly as moment problem.

Returns a named tuple `(objective, model, n_unique_elements, status)`.

Throws an error if the solver fails (infeasible, unbounded, numerical error).
"""
function solve_sdp(moment_problem, optimizer; dualize::Bool=true)
    if dualize
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, optimizer)
        optimize!(sos_problem.model)
        status = _check_solver_status(sos_problem.model)
        return (
            objective = objective_value(sos_problem.model),
            model = sos_problem.model,
            n_unique_elements = moment_problem.n_unique_moment_matrix_elements,
            status = status
        )
    else
        result = solve_moment_problem(moment_problem, optimizer)
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
    SolverConfig(; optimizer, order, cs_algo=NoElimination(), ts_algo=NoElimination())

Configuration for solving polynomial optimization problems.

# Keyword Arguments
- `optimizer` (required): The optimizer to use for solving the SDP problem (e.g. Clarabel.Optimizer)
- `order::Int`: The order of the moment relaxation (default: 0)
- `cs_algo::EliminationAlgorithm`: Algorithm for correlative sparsity exploitation (default: NoElimination())
- `ts_algo::EliminationAlgorithm`: Algorithm for term sparsity exploitation (default: NoElimination())


# Examples
```jldoctest; setup=:(using NCTSSoS, Clarabel)
julia> solver_config = SolverConfig(optimizer=Clarabel.Optimizer, order=2) # default elimination algorithms
SolverConfig(Clarabel.MOIwrapper.Optimizer, 2, NoElimination(), NoElimination())
```

"""
@kwdef struct SolverConfig
    optimizer
    order::Int = 0
    cs_algo::EliminationAlgorithm = NoElimination()
    ts_algo::EliminationAlgorithm = NoElimination()
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
    order = compute_relaxation_order(pop, solver_config.order)
    corr_sparsity = correlative_sparsity(pop, order, solver_config.cs_algo)

    # Compute partial objectives for each clique
    cliques_objective = map(c -> project_to_clique(pop.objective, c), corr_sparsity.cliques)

    initial_activated_supps_nm = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps_nm, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    return SparsityResult(corr_sparsity, initial_activated_supps_nm, cliques_term_sparsities)
end

"""
    compute_sparsity(pop::PolyOpt, solver_config::SolverConfig) -> SparsityResult

Compute correlative and term sparsity for a state polynomial optimization problem.
"""
function compute_sparsity(pop::PolyOpt{A,T,P}, solver_config::SolverConfig) where {A<:AlgebraType,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    order = compute_relaxation_order(pop, solver_config.order)
    corr_sparsity = correlative_sparsity(pop, order, solver_config.cs_algo)

    # Compute partial objectives for each clique
    cliques_objective = map(c -> project_to_clique(pop.objective, c), corr_sparsity.cliques)

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
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
function cs_nctssos(pop::OP, solver_config::SolverConfig; dualize::Bool=true) where {A<:AlgebraType, P, OP<:OptimizationProblem{A,P}}
    sparsity = compute_sparsity(pop, solver_config)
    moment_problem = moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    result = solve_sdp(moment_problem, solver_config.optimizer; dualize)
    return PolyOptResult(result.objective, sparsity, result.model, result.n_unique_elements)
end

"""
    cs_nctssos_higher(pop::PolyOpt{T}, prev_res::PolyOptResult, solver_config::SolverConfig; dualize::Bool=true) where {T}

Solve a polynomial optimization problem using higher-order term sparsity based on a previous result.

# Arguments
- `pop::PolyOpt{T}`: The polynomial optimization problem to solve
- `prev_res::PolyOptResult`: Previous optimization result containing sparsity information to build upon
- `solver_config::SolverConfig`: Configuration containing optimizer and sparsity algorithms

# Keyword Arguments
- `dualize::Bool=true`: Whether to dualize the moment relaxation to a sum-of-squares problem

# Returns
- `PolyOptResult`: Result containing the objective value, correlative sparsity structure, and updated term sparsity information

# Description
This function performs a higher-order iteration of the CS-NCTSSOS method by:
1. Using the correlative sparsity structure from the previous result
2. Computing new term sparsity based on the union of previously activated supports
3. Formulating and solving either the moment relaxation or its SOS dual with the refined sparsity
4. Returning the optimal objective value and updated sparsity information

This is typically used when the previous relaxation was not tight enough and a higher-order relaxation is needed.
"""
function cs_nctssos_higher(pop::OP, prev_res::PolyOptResult, solver_config::SolverConfig; dualize::Bool=true) where {A<:AlgebraType, P, OP<:OptimizationProblem{A,P}}
    prev_sparsity = prev_res.sparsity

    # Compute new initial activated supports from union of previous term sparsity graph supports
    initial_activated_supps_nm = [
        reduce(sorted_union, [poly_term_sparsity.term_sparse_graph_supp for poly_term_sparsity in term_sparsities_vec])
        for term_sparsities_vec in prev_sparsity.cliques_term_sparsities
    ]

    prev_corr_sparsity = prev_sparsity.corr_sparsity

    cliques_term_sparsities = map(zip(initial_activated_supps_nm, prev_corr_sparsity.clq_cons, prev_corr_sparsity.clq_mom_mtx_bases, prev_corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, prev_corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    # Create new sparsity result with updated term sparsities
    sparsity = SparsityResult(prev_corr_sparsity, initial_activated_supps_nm, cliques_term_sparsities)

    moment_problem = moment_relax(pop, prev_corr_sparsity, cliques_term_sparsities)

    result = solve_sdp(moment_problem, solver_config.optimizer; dualize)
    return PolyOptResult(result.objective, sparsity, result.model, result.n_unique_elements)
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
function cs_nctssos(pop::PolyOpt{A,T,P}, solver_config::SolverConfig; dualize::Bool=true) where {A<:AlgebraType,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    sparsity = compute_sparsity(pop, solver_config)
    moment_problem = moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    result = solve_sdp(moment_problem, solver_config.optimizer; dualize)
    return PolyOptResult(result.objective, sparsity, result.model, result.n_unique_elements)
end
