# Compare graph-selected joint bases with arXiv:2604.01555v2, Table 3.
# Usage: julia --project=<licensed Mosek environment> perf/pauli_physics_basis_gap.jl N RANGE [ORDER] [moment|sos]
# This is a bound-quality diagnostic, not a repeated performance benchmark.
using NCTSSoS, JuMP, MosekTools, LinearAlgebra

function main()
    n = parse(Int,ARGS[1])
    distance = parse(Int,ARGS[2])
    order = length(ARGS) >= 3 ? parse(Int,ARGS[3]) : 4
    form = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :moment
    form in (:moment,:sos) || error("Form must be moment or sos.")
    paper = Dict(10=>-0.4515446,20=>-0.4452196,100=>-0.4432378)
    reg, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops),reg)
    started = time()
    basis = physics_moment_basis(pop,order;pair_distance=distance)
    println("INPUT n=$n order=$order pair_distance=$distance rows=$(length(basis)) degree_counts=$([count(m->degree(m)==d,basis) for d in 0:order])")
    println("SYMMETRIES translation=true sign=true reflection=true conjugation=true axis_permutations=true continuous=false")
    println("STRENGTHENING rdm=[] state_optimality=none form=$form")
    flush(stdout)
    mp, report = pauli_translation_invariant_moment_relaxation(pop,ops,order;
        basis,axis_permutation_symmetry=true)
    println("BUILT seconds=$(time()-started) blocks=$(length(report.psd_block_sizes)) max_block=$(maximum(report.psd_block_sizes)) moments=$(length(mp.linear.moments))")
    flush(stdout)
    model = form == :moment ? first(build_jump_model(mp)) : NCTSSoS.sos_dualize(mp).model
    set_optimizer(model,optimizer_with_attributes(MosekTools.Optimizer,
        "MSK_IPAR_NUM_THREADS"=>1,"MSK_IPAR_LOG"=>0,
        "MSK_DPAR_OPTIMIZER_MAX_TIME"=>600.0,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-9,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-9,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-9))
    set_silent(model)
    println("MODELED seconds=$(time()-started)")
    flush(stdout)
    optimize!(model)
    primal, dual = objective_value(model),dual_objective_value(model)
    violation = maximum(values(primal_feasibility_report(model;atol=0.0));init=0.0)
    gap = abs(primal-dual)/max(1.0,abs(primal),abs(dual))
    accepted = termination_status(model)==JuMP.MOI.OPTIMAL &&
        primal_status(model)==JuMP.MOI.FEASIBLE_POINT && dual_status(model)==JuMP.MOI.FEASIBLE_POINT &&
        violation<=1e-7 && gap<=1e-7
    target = get(paper,n,NaN)
    lower_estimate = form == :moment ? dual : primal
    println("RESULT n=$n order=$order range=$distance primal_per_spin=$(primal/n) dual_per_spin=$(dual/n) lower_estimate_per_spin=$(lower_estimate/n) paper=$target shortfall=$(target-lower_estimate/n) reduced_violation=$violation relative_gap=$gap status=$(termination_status(model)) primal_status=$(primal_status(model)) dual_status=$(dual_status(model)) accepted=$accepted elapsed=$(time()-started)")
    accepted || error("Numerical quality gates failed; do not count this as an accepted bound comparison.")
end

main()
