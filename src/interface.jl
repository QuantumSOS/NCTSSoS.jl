struct PolyOptResult{RT,T,P,M}
    objective::RT # support for high precision solution
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{RT}
    moment_matrices::Vector{Vector{Matrix{T}}}
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
function cs_nctssos(pop::OP, solver_config::SolverConfig; dualize::Bool=true) where {P,OP<:OptimizationProblem{P}}
    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)
    order = iszero(solver_config.order) ? maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order

    corr_sparsity = correlative_sparsity(pop, order, solver_config.cs_algo)

    cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
    end

    cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
    end

    moment_problem = moment_relax(pop, corr_sparsity, cliques_term_sparsities)

    (pop isa ComplexPolyOpt{P} && !dualize) && error("Solving Moment Problem for Complex Poly Opt is not supported")

    if dualize
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities, sos_problem.model, get_moment_matrices(sos_problem))
    else
        set_optimizer(moment_problem.model, solver_config.optimizer)
        optimize!(moment_problem.model)
        return PolyOptResult(objective_value(moment_problem.model), corr_sparsity, cliques_term_sparsities, moment_problem.model, get_moment_matrices(moment_problem.model, cliques_term_sparsities))
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
function cs_nctssos_higher(pop::OP, prev_res::PolyOptResult, solver_config::SolverConfig; dualize::Bool=true) where {T,OP<:OptimizationProblem{T}}
    sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)

    initial_activated_supps = [sorted_union([poly_term_sparsity.term_sparse_graph_supp for poly_term_sparsity in term_sparsities]...)
                               for term_sparsities in prev_res.cliques_term_sparsities]

    prev_corr_sparsity = prev_res.corr_sparsity

    cliques_term_sparsities = map(zip(initial_activated_supps, prev_corr_sparsity.clq_cons, prev_corr_sparsity.clq_mom_mtx_bases, prev_corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        term_sparsities(init_act_supp, prev_corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
    end


    moment_problem = moment_relax(pop, prev_res.corr_sparsity, cliques_term_sparsities)

    if dualize
        sos_problem = sos_dualize(moment_problem)
        set_optimizer(sos_problem.model, solver_config.optimizer)
        optimize!(sos_problem.model)
        return PolyOptResult(objective_value(sos_problem.model), prev_res.corr_sparsity, cliques_term_sparsities, sos_problem.model, get_moment_matrices(sos_problem))
    else
        set_optimizer(moment_problem.model, solver_config.optimizer)
        optimize!(moment_problem.model)
        return PolyOptResult(objective_value(moment_problem.model), prev_res.corr_sparsity, cliques_term_sparsities, moment_problem.model, get_moment_matrices(moment_problem.model, cliques_term_sparsities))
    end
end

"""
    get_moment_matrices(result::PolyOptResult)

Extract moment matrices from a solved polynomial optimization problem.

For primal formulation, retrieves the moment matrices directly from the named
PSD constraints in the JuMP model. Each moment matrix is stored with the name
`mom_mtx_clique_i_block_j` where i is the clique index and j is the block index.

# Arguments
- `result::PolyOptResult`: The solved optimization result

# Returns
- `Vector{Vector{Matrix{T}}}`: Hierarchical structure of moment matrices:
  - Outer vector: One entry per clique
  - Inner vector: One entry per block within the clique
  - Matrix: The moment matrix values

# Example
```julia
result = cs_nctssos(pop, solver_config, dualize=false)
moments = get_moment_matrices(result)
M11 = moments[1][1]  # Moment matrix for clique 1, block 1
```
"""
function get_moment_matrices(model::GenericModel{T}, cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}) where {T,M}
    num_cliques = length(cliques_term_sparsities)

    moment_mats = Vector{Vector{Matrix{T}}}(undef, num_cliques)

    for clq_idx in 1:num_cliques
        # Get number of blocks in first term sparsity (moment matrix blocks)
        num_blocks = length(cliques_term_sparsities[clq_idx][1].block_bases)

        clique_mats = Vector{Matrix{T}}(undef, num_blocks)

        for blk_idx in 1:num_blocks
            constraint_name = Symbol("mom_mtx_clique_$(clq_idx)_block_$(blk_idx)")

            clique_mats[blk_idx] = value(model[constraint_name])
        end

        moment_mats[clq_idx] = clique_mats
    end

    return moment_mats
end
