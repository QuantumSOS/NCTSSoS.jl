"""
Benchmark for CHSH Inequality using NCTSSoS.jl
Collects: Matrix size, Variables, Constraints, Set-up time
"""

using NCTSSoS
using NCTSSoS: polyopt, PolyOptResult, SolverConfig
using MosekTools
using JuMP
using Printf

function get_problem_stats(model::JuMP.GenericModel)
    """Extract problem statistics from the JuMP model
    
    Returns:
    - num_vars: Number of scalar decision variables (excluding the fixed variable for 1)
    - num_cons: Number of constraints (excluding variable bounds)
    - max_matrix_size: Size of the largest PSD matrix block
    - total_matrix_size: Sum of all PSD matrix block sizes
    - matrix_blocks: Number of PSD matrix blocks
    """
    
    # Get number of variables (excluding the fixed variable for 1)
    # The first variable in moment relaxation is always fixed to 1
    num_vars = num_variables(model) - 1
    
    # Get number of constraints (excluding variable bounds)
    num_cons = num_constraints(model; count_variable_in_set_constraints=false)
    
    # Get matrix sizes from PSD constraints
    matrix_sizes = Int[]
    
    # Iterate through all constraint types
    for (F, S) in list_of_constraint_types(model)
        # Check for PSD constraints (both Triangle and Square representations)
        if S <: Union{MOI.PositiveSemidefiniteConeTriangle, MOI.PositiveSemidefiniteConeSquare}
            constraints = all_constraints(model, F, S)
            for con in constraints
                # Get the dimension from the constraint set
                set = constraint_object(con).set
                if hasfield(typeof(set), :side_dimension)
                    push!(matrix_sizes, set.side_dimension)
                elseif hasfield(typeof(set), :dimension)
                    # For square form, dimension is n^2, so we need sqrt
                    dim = isqrt(set.dimension)
                    push!(matrix_sizes, dim)
                end
            end
        end
    end
    
    total_matrix_size = sum(matrix_sizes)
    max_matrix_size = isempty(matrix_sizes) ? 0 : maximum(matrix_sizes)
    
    return (
        num_vars = num_vars,
        num_cons = num_cons,
        max_matrix_size = max_matrix_size,
        total_matrix_size = total_matrix_size,
        matrix_blocks = length(matrix_sizes),
        all_matrix_sizes = matrix_sizes  # Keep all sizes for detailed analysis
    )
end

function benchmark_chsh(order::Int)
    """Benchmark CHSH inequality for a given relaxation order"""
    
    println("\n" * "="^70)
    println("Benchmarking CHSH Inequality - Relaxation Order: $order")
    println("="^70)
    
    # Setup the CHSH problem
    @ncpolyvar x[1:2]
    @ncpolyvar y[1:2]
    
    # CHSH objective: x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    
    # Create polynomial optimization problem
    # Variables are unipotent (squared to 1)
    # x and y are in different commuting groups
    pop = polyopt(f; comm_gps=[x, y], is_unipotent=true)
    
    # Configure solver
    solver_config = SolverConfig(
        optimizer=Mosek.Optimizer, 
        order=order
    )
    
    # Measure ONLY setup time (not solving time)
    result = nothing
    setup_time = @elapsed begin
        # This is the setup phase - building the problem
        sa = NCTSSoS.SimplifyAlgorithm(comm_gps=pop.comm_gps, is_projective=pop.is_projective, is_unipotent=pop.is_unipotent)
        actual_order = iszero(solver_config.order) ? maximum([ceil(Int, NCTSSoS.maxdegree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order
        
        corr_sparsity = NCTSSoS.correlative_sparsity(pop, actual_order, solver_config.cs_algo)
        
        cliques_objective = [reduce(+, [issubset(sort!(NCTSSoS.variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(NCTSSoS.coefficients(pop.objective), NCTSSoS.monomials(pop.objective))]) for clique in corr_sparsity.cliques]
        
        initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
            NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base, sa)
        end
        
        cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
            NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo, sa)
        end
        
        moment_problem = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)
        
        set_optimizer(moment_problem.model, solver_config.optimizer)
        
        # Store the model for statistics extraction
        result_model = moment_problem.model
    end
    
    # Extract problem statistics (WITHOUT solving)
    stats = get_problem_stats(result_model)
    
    # Print results - focus on setup time
    println("\nResults:")
    println("  Relaxation Order:         $order")
    println("  Matrix Size (largest):    $(stats.max_matrix_size)")
    println("  Total Matrix Size:        $(stats.total_matrix_size)")
    println("  Number of PSD Blocks:     $(stats.matrix_blocks)")
    println("  Scalar Variables:         $(stats.num_vars)")
    println("  Constraints:              $(stats.num_cons)")
    println("  Set-up Time:              $(@sprintf("%.4f", setup_time)) seconds")
    
    return (
        order = order,
        matrix_size = stats.max_matrix_size,
        total_matrix_size = stats.total_matrix_size,
        num_blocks = stats.matrix_blocks,
        num_vars = stats.num_vars,
        num_cons = stats.num_cons,
        setup_time = setup_time
    )
end

function print_summary_table(results)
    """Print a formatted summary table of all results"""
    
    println("\n" * "="^80)
    println("CHSH Inequality Benchmark Summary")
    println("="^80)
    println()
    
    # Header
    println(@sprintf("%-5s | %-12s | %-10s | %-10s | %-12s", 
        "Order", "Matrix Size", "Variables", "Constraints", "Set-up Time"))
    println("-"^80)
    
    # Data rows
    for r in results
        println(@sprintf("%-5d | %-12d | %-10d | %-10d | %-12.4f", 
            r.order, r.matrix_size, r.num_vars, r.num_cons, r.setup_time))
    end
    
    println("="^80)
    println()
end

# Main execution
function main()
    println("CHSH Inequality Benchmark using NCTSSoS.jl")
    println("Testing relaxation orders: 4, 5, 6")
    
    orders = [4, 5]
    results = []
    
    for order in orders
        try
            result = benchmark_chsh(order)
            push!(results, result)
        catch e
            println("Error benchmarking order $order:")
            println(e)
            println(stacktrace(catch_backtrace()))
        end
    end
    
    # Print summary table
    if !isempty(results)
        print_summary_table(results)
    end
    
    return results
end

main()
