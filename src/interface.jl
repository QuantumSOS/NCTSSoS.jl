"""
    PolyOptResult{T, A<:AlgebraType, TI<:Integer, P<:Polynomial{A,TI}, M<:Monomial{A,TI}}

Result of a polynomial optimization problem solution.

# Type Parameters
- `T`: Coefficient type for the objective value (supports BigFloat for high precision)
- `A`: Algebra type
- `TI`: Index type
- `P`: Polynomial type
- `M`: Monomial type

# Fields
- `objective::T`: Optimal objective value
- `corr_sparsity::CorrelativeSparsity{A,TI,P,M}`: Correlative sparsity structure
- `cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}`: Term sparsity for each clique
- `model::GenericModel{T}`: JuMP model used for solving
"""
struct PolyOptResult{T, A<:AlgebraType, TI<:Integer, P<:Polynomial{A,TI}, M<:Monomial{A,TI}}
    objective::T
    corr_sparsity::CorrelativeSparsity{A,TI,P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
end

function Base.show(io::IO, result::PolyOptResult)
    println(io, "Objective: ", result.objective)
    show(io, result.corr_sparsity)
    println(io, "Term Sparsity:")
    for (i, sparsities) in enumerate(result.cliques_term_sparsities)
        println(io, "Clique $i:")
        println(io, "   Moment Matrix:")
        println(io, sparsities[1])
        println(io, "   Localizing Matrix:")
        for sparsity in sparsities[2:end]
            show(io, sparsity)
        end
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

    order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order

    corr_sparsity = correlative_sparsity(pop, order, solver_config.cs_algo)

    # Compute partial objectives for each clique
    # cliques are now Vector{T} (indices), use variable_indices() for comparison
    cliques_objective = map(corr_sparsity.cliques) do clique_indices
        clique_set = Set(clique_indices)
        reduce(+, [
            issubset(variable_indices(mono), clique_set) ? coef * mono : zero(coef) * one(mono)
            for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))
        ])
    end

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    # Create unified symbolic moment problem (handles both real and complex via algebra traits)
    moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

    if dualize
        # Dualize to SOS problem and solve
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities, sos_problem.model)
    else
        # Solve moment problem directly
        result = solve_moment_problem(moment_problem, solver_config.optimizer)
        return PolyOptResult(result.objective, corr_sparsity, cliques_term_sparsities, result.model)
    end
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
    initial_activated_supps = [sorted_union([poly_term_sparsity.term_sparse_graph_supp for poly_term_sparsity in term_sparsities_vec]...)
                               for term_sparsities_vec in prev_res.cliques_term_sparsities]

    prev_corr_sparsity = prev_res.corr_sparsity

    cliques_term_sparsities = map(zip(initial_activated_supps, prev_corr_sparsity.clq_cons, prev_corr_sparsity.clq_mom_mtx_bases, prev_corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, prev_corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end


    moment_problem = moment_relax(pop, prev_res.corr_sparsity, cliques_term_sparsities)

    if dualize
        # Dualize to SOS problem and solve
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), prev_res.corr_sparsity, cliques_term_sparsities, sos_problem.model)
    else
        # Solve moment problem directly
        result = solve_moment_problem(moment_problem, solver_config.optimizer)
        return PolyOptResult(result.objective, prev_res.corr_sparsity, cliques_term_sparsities, result.model)
    end
end

# =============================================================================
# State Polynomial Optimization Result
# =============================================================================

"""
    StatePolyOptResult{T, A<:AlgebraType, ST<:StateType, TI<:Integer, P<:NCStatePolynomial, M<:NCStateWord}

Result of a state polynomial optimization problem solution.

# Type Parameters
- `T`: Coefficient type for the objective value
- `A`: Algebra type
- `ST`: State type
- `TI`: Index type
- `P`: NCStatePolynomial type
- `M`: NCStateWord type

# Fields
- `objective::T`: Optimal objective value
- `corr_sparsity::StateCorrelativeSparsity{A,ST,TI,P,M}`: State correlative sparsity structure
- `cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}`: Term sparsity for each clique
- `model::GenericModel{T}`: JuMP model used for solving
"""
struct StatePolyOptResult{T, A<:AlgebraType, ST<:StateType, TI<:Integer, P<:NCStatePolynomial, M<:NCStateWord}
    objective::T
    corr_sparsity::StateCorrelativeSparsity{A,ST,TI,P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
end

function Base.show(io::IO, result::StatePolyOptResult)
    println(io, "State Optimization Result")
    println(io, "Objective: ", result.objective)
    show(io, result.corr_sparsity)
end

# =============================================================================
# State Polynomial Optimization Solver
# =============================================================================

"""
    cs_nctssos(pop::StatePolyOpt{A,ST,P}, solver_config::SolverConfig; dualize::Bool=true)

Solve a state polynomial optimization problem using the CS-NCTSSOS method.

# Arguments
- `pop::StatePolyOpt{A,ST,P}`: The state polynomial optimization problem to solve
- `solver_config::SolverConfig`: Configuration containing optimizer, moment order, and sparsity algorithms

# Keyword Arguments
- `dualize::Bool=true`: Whether to dualize the moment relaxation to a sum-of-squares problem

# Returns
- `StatePolyOptResult`: Result containing the objective value and sparsity information

# Description
This function solves a state polynomial optimization problem by:
1. Computing correlative sparsity to decompose the problem into smaller cliques
2. Computing term sparsity for each clique to further reduce problem size
3. Formulating and solving the SOS dual of the moment relaxation
4. Returning the optimal objective value and sparsity information
"""
function cs_nctssos(pop::StatePolyOpt{A,ST,P}, solver_config::SolverConfig; dualize::Bool=true) where {A<:AlgebraType, ST<:StateType, C<:Number, T<:Integer, P<:NCStatePolynomial{C,ST,A,T}}

    order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order

    corr_sparsity = correlative_sparsity(pop, order, solver_config.cs_algo)

    # Compute partial objectives for each clique
    # For NCStatePolynomial, use variable_indices on NCStateWord
    cliques_objective = map(corr_sparsity.cliques) do clique_indices
        clique_set = Set(clique_indices)
        result_coeffs = C[]
        result_ncsws = NCStateWord{ST,A,T}[]
        for (coef, ncsw) in zip(coefficients(pop.objective), monomials(pop.objective))
            if issubset(variable_indices(ncsw), clique_set)
                push!(result_coeffs, coef)
                push!(result_ncsws, ncsw)
            end
        end
        isempty(result_coeffs) ? zero(P) : NCStatePolynomial(result_coeffs, result_ncsws)
    end

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    # Create symbolic state moment problem
    moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

    M = NCStateWord{ST,A,T}

    if dualize
        # Dualize to SOS problem and solve
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return StatePolyOptResult{C,A,ST,T,P,M}(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities, sos_problem.model)
    else
        # Solve moment problem directly
        result = solve_moment_problem(moment_problem, solver_config.optimizer)
        return StatePolyOptResult{C,A,ST,T,P,M}(result.objective, corr_sparsity, cliques_term_sparsities, result.model)
    end
end
