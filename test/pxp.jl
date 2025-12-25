using NCTSSoS, NCTSSoS.FastPolynomials
using MosekTools
using JuMP

# TODO: do 3x3 case with various operators' expectation value

function cs_nctssos_with_blockade(pop::OP, solver_config::SolverConfig, blockade_constraints::Vector{Polynomial{T}}, eigen_state_constraints::Vector{Polynomial{T}}; dualize::Bool=true) where {T,P<:Polynomial{T},OP<:NCTSSoS.OptimizationProblem{P}}

   order = iszero(solver_config.order) ? maximum([ceil(Int, degree(poly) / 2) for poly in [pop.objective; pop.eq_constraints; pop.ineq_constraints]]) : solver_config.order

   corr_sparsity = NCTSSoS.correlative_sparsity(pop, order, solver_config.cs_algo)

   cliques_objective = [reduce(+, [issubset(sort!(variables(mono)), clique) ? coef * mono : zero(coef) * one(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

   initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
        NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
   end

   cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
        NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, solver_config.ts_algo)
   end

   moment_problem = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)

   for c in blockade_constraints
     push!(moment_problem.constraints,(:Zero, [c;;]))
   end

   for c in eigen_state_constraints
     push!(moment_problem.constraints,(:HPSD, [c;;]))
   end

   # recover blockade constraints
    for (type, cons) in moment_problem.constraints
        type == :HPSD && continue
        if cons[1, 1] in blockade_constraints
            cons[2:end, 2:end] .*= zero(T)
            # why can we ignore those zeroing?
            # cons[2:end, 1] .*= zero(T)
            # cons[1, 2:end] .*= zero(T)
        end
        # @show cons
    end

    for (type, cons) in moment_problem.constraints
        type == :Zero && continue
        if cons[1, 1] in eigen_state_constraints
            cons[2:end, 2:end] .*= zero(T)
            cons[2:end, 1] .*= zero(T)
            cons[1, 2:end] .*= zero(T)
        end
        # @show cons
    end

   # Note: All algebra types now support moment relaxation
   problem_to_solve = !dualize ? moment_problem : NCTSSoS.sos_dualize(moment_problem)

   set_optimizer(problem_to_solve.model, solver_config.optimizer)
   optimize!(problem_to_solve.model)
   return NCTSSoS.PolyOptResult(objective_value(problem_to_solve.model), corr_sparsity, cliques_term_sparsities, problem_to_solve.model)
end

function solve_n1n2_bounds(target_energy, cur_spreading)
    # hard-code problem into the function
    data = Dict("Detuning" => [12.786857270237986, 15.486419419279784, 13.363757638514809, 26.836874327070852], "Exact GSE" => -39.2040626133127, "Vanderwaals" => [[1, 4, 81.59613025943814], [2, 3, 81.59613025943814]], "Rabi" => [25.132741228718345, 25.132741228718345, 25.132741228718345, 25.132741228718345], "PXP" => [[1, 2], [1, 3], [2, 4], [3, 4]])

    T = ComplexF64
    Lx = 2
    Ly = 2
    N = Lx * Ly
    reg, (σx, σy, σz) = create_pauli_variables(1:N)

    H = -sum(data["Detuning"] ./ 2 .* (ones(T, N) .- σz)) + sum(data["Rabi"] ./ 2 .* σx) + sum(map(V -> V[3] / 4 * (1 - σz[Int(V[1])]) * (1 - σz[Int(V[2])]), data["Vanderwaals"]); init=zero(T))

    Pauli_algebra =
        reduce(vcat, [[σx[i] * σy[i] - im * σz[i], σy[i] * σx[i] + im * σz[i], σy[i] * σz[i] - im * σx[i], σz[i] * σy[i] + im * σx[i], σz[i] * σx[i] - im * σy[i], σx[i] * σz[i] + im * σy[i]] for i in 1:N])

    blockade_constraints = [one(T) - σz[e[1]] - σz[e[2]] + σz[e[1]] * σz[e[2]] for e in data["PXP"]]

    n1n2 = ((one(T) - σz[1]) / 2) * ((one(T) - σz[2]) / 2)


    pop_lower = polyopt(n1n2, reg)

    SOLVER = optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8, "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8, "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)))

    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    energy_cons = [-(H - target_energy) * (H - target_energy) + cur_spreading]

    res_lower = cs_nctssos_with_blockade(pop_lower, solver_config, blockade_constraints, energy_cons)

    pop_upper = polyopt(-n1n2, reg)

    res_upper = cs_nctssos_with_blockade(pop_upper, solver_config, blockade_constraints, energy_cons)

    return res_lower, res_upper
end

# get energy spectrum
energy_spectrum = [-39.20406261331271,
    -24.472055568089665,
    -15.90113602877199,
    -14.262691405978654,
    15.245634541416734,
    47.298087863292544,
    57.54066642011262]

res_lower, res_upper = solve_n1n2_bounds(energy_spectrum[1], 0.1)

function does_solve(target_energy, spreading_upper)
    spreading_lower = 0.0
    cur_spreading = (spreading_upper + spreading_lower) / 2
    is_solved = false
    while !is_solved
        res_lower, res_upper = solve_n1n2_bounds(target_energy, cur_spreading)
        if is_solved_and_feasible(res_upper.model) && is_solved_and_feasible(res_lower.model)
            if (spreading_upper - spreading_lower) < 1e-3
                @show "Both models are solved and feasible"
                return cur_spreading
            else
                @info "Lowering current bound $(cur_spreading)"
                spreading_upper, cur_spreading = cur_spreading, (cur_spreading + spreading_lower) / 2
            end
        else
            @info "Raising current bound $(cur_spreading)"
            spreading_lower, cur_spreading = cur_spreading, (spreading_upper + spreading_lower) / 2
        end
    end
end

tight_spreading = [does_solve(energy_spectrum[i], 0.001) for i in 1:length(energy_spectrum)]

tight_spreadings = [277.8101921081543, 47.6008415222168, 66.45956039428711, 82.20663070678711, 156.52642250061035, 309.25512313842773, 312.9763603210449]
