using NCTSSoS
using MosekTools
using JuMP

# Setup the CHSH problem
@ncpolyvar x[1:2]
@ncpolyvar y[1:2]
f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
pop = polyopt(f; comm_gps=[x, y], is_unipotent=true)

# Configure solver
solver_config = SolverConfig(optimizer=Mosek.Optimizer, order=4)

# Build the problem
sa = NCTSSoS.SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)
actual_order = 4

corr_sparsity = NCTSSoS.correlative_sparsity(pop, actual_order, solver_config.cs_algo)

cliques_objective = [reduce(+, [issubset(sort!(NCTSSoS.variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(NCTSSoS.coefficients(pop.objective), NCTSSoS.monomials(pop.objective))]) for clique in corr_sparsity.cliques]

initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
    NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
end

cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
    NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
end

moment_problem = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)
model = moment_problem.model

println("\n=== Model Statistics ===")
println("Number of variables: ", num_variables(model))
println("\n=== Constraint Types ===")
for (F, S) in list_of_constraint_types(model)
    num = num_constraints(model, F, S)
    println("$F in $S: $num constraints")
    
    # Try to extract PSD matrix dimensions
    if S <: Union{MOI.PositiveSemidefiniteConeTriangle, MOI.PositiveSemidefiniteConeSquare}
        println("  This is a PSD constraint type!")
        constraints = all_constraints(model, F, S)
        for (i, con) in enumerate(constraints)
            set = MOI.get(model, MOI.ConstraintSet(), con)
            println("    Constraint $i: set = $set, type = $(typeof(set))")
            if hasfield(typeof(set), :side_dimension)
                println("      side_dimension = $(set.side_dimension)")
            elseif hasfield(typeof(set), :dimension)
                println("      dimension = $(set.dimension)")
            end
        end
    end
end

println("\n=== Total Constraints ===")
println("Total (excluding variable bounds): ", num_constraints(model; count_variable_in_set_constraints=false))
