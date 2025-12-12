struct PolyOptResult{T,P,M}
    objective::T # support for high precision solution
    corr_sparsity::CorrelativeSparsity{P,M}
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

    cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
    end

    # Dispatch between real and complex moment relaxation based on algebra traits
    moment_problem = if _is_complex_problem(A)
        !dualize && error("Solving Moment Problem for complex algebra is not supported (dualize=false)")
        complex_moment_relax(pop, corr_sparsity, cliques_term_sparsities)
    else
        moment_relax(pop, corr_sparsity, cliques_term_sparsities)
    end
    problem_to_solve = !dualize ? moment_problem : sos_dualize(moment_problem)

    set_optimizer(problem_to_solve.model, solver_config.optimizer)
    optimize!(problem_to_solve.model)
    return PolyOptResult(objective_value(problem_to_solve.model), corr_sparsity, cliques_term_sparsities, problem_to_solve.model)
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

    problem_to_solve = !dualize ? moment_problem : sos_dualize(moment_problem)

    set_optimizer(problem_to_solve.model, solver_config.optimizer)
    optimize!(problem_to_solve.model)
    return PolyOptResult(objective_value(problem_to_solve.model), prev_res.corr_sparsity, cliques_term_sparsities, problem_to_solve.model)
end
