# Bound-quality comparison, not a repeated timing benchmark.
# Use a local environment with NCTSSoS, Graphs, JuMP and licensed MosekTools.
using NCTSSoS, Graphs, JuMP, MosekTools, LinearAlgebra, SHA
import NCTSSoS: degree
const N = NCTSSoS
const MOI = JuMP.MOI
const PERTURBATION_ORDER = parse(Int,get(ENV,"PERTURBATION_ORDER","2"))
const SOLVER_TOLERANCE = parse(Float64,get(ENV,"SOLVER_TOLERANCE","1e-9"))
@assert PERTURBATION_ORDER in (1,2)

function screen(mp, solver; formulation=Symbol(get(ENV,"FORMULATION","moment_variables")))
    started = time()
    model = formulation == :sos ? N.sos_dualize(mp).model : first(build_jump_model(mp;formulation))
    set_optimizer(model,solver)
    set_silent(model)
    built = time()
    optimize!(model)
    solved = time()
    if !has_values(model) || !has_duals(model)
        return (;accepted=false,status=termination_status(model),lower=NaN,primal=NaN,
            gap=Inf,violation=Inf,build_seconds=built-started,solve_seconds=solved-built)
    end
    primal,dual = objective_value(model),dual_objective_value(model)
    lower = objective_sense(model) == MOI.MIN_SENSE ? dual : primal
    gap = abs(primal-dual)/max(1.0,abs(primal),abs(dual))
    violation = maximum(values(primal_feasibility_report(model;atol=0.0));init=0.0)
    accepted = termination_status(model) == MOI.OPTIMAL &&
        primal_status(model) == dual_status(model) == MOI.FEASIBLE_POINT &&
        gap <= 1e-7 && violation <= 1e-7
    return (;accepted,status=termination_status(model),lower,primal,gap,violation,
        build_seconds=built-started,solve_seconds=solved-built)
end

function chain_comparison(n, solver; control=false)
    started = time()
    reg,ops = create_pauli_variables(1:n)
    h = heisenberg_chain_hamiltonian(ops;coupling=control ? -1/4 : 1/4)
    pop = polyopt(h,reg)
    radius = control ? 0 : 1
    basis = physics_moment_basis(reg,cycle_graph(n);radius,max_rows=2000)
    config = SolverConfig(optimizer=solver,moment_basis=basis)
    # Exact original matrices check moment coverage BEFORE symmetry projection.
    original = N._continuous_exact_problem(pop,compute_sparsity(pop,config))
    k = control ? 2 : 4
    patches = [sort(mod1.(i:i+k-1,n)) for i in 1:n]
    qs = reduce(vcat,ops)
    if PERTURBATION_ORDER == 2
        append!(qs,[only(monomials(a[i]*b[mod1(i+1,n)])) for i in 1:n for a in ops for b in ops])
    end
    extra = N._pauli_strengthening_constraints(pop,original;subsystems=patches,
        perturbations=qs,stationarity=true,max_block_size=256)
    N._check_translation_basis_closure(extra.qs,n)
    N._check_reflection_basis_closure(extra.qs,n)
    N._check_pauli_axis_permutation_basis_closure(extra.qs)
    main,report = pauli_translation_invariant_moment_relaxation(pop,ops,maximum(degree,basis);
        basis,axis_permutation_symmetry=true)
    # The chain interface validates the Hamiltonian and basis under each action.
    # Reuse those SAME moment identifications, including sign averaging, for the
    # exact generic strengthening matrices. Do not use the chain backend's
    # different built-in optimality basis or introduce any extra moments.
    reducer = N._PauliMomentReducer{eltype(basis)}(n,true,true,true)
    function project(p)
        even = Polynomial([(c,m) for (c,m) in terms(p) if iszero(N._pauli_sign_signature(m))])
        return N._reduce_moment_polynomial(even,reducer)
    end
    projected = [(label.kind,map(project,mat)) for (label,(_,mat)) in zip(extra.labels,extra.constraints)]
    densities = [mat for (kind,mat) in projected if kind == :subsystem]
    # Sorted site labels reorder tensor factors on windows crossing the ring
    # boundary. Verify the exact basis permutation before keeping one window.
    for (label,mat) in zip(extra.labels,densities)
        start = findfirst(==(label.sites),patches)
        cyclic = mod1.(start:start+k-1,n)
        positions = indexin(cyclic,label.sites)
        permutation = [1+sum(((b >> (k-j)) & 1) << (k-positions[j]) for j in 1:k)
            for b in 0:(2^k-1)]
        @assert mat[permutation,permutation] == first(densities)
    end
    perturbative = only([mat for (kind,mat) in projected if kind == :perturbative])
    equalities = filter(!iszero,project.(extra.equalities))
    domain = Set(main.total_basis)
    @assert all(m in domain for (_,mat) in projected for p in mat for (_,m) in terms(p))
    @assert all(m in domain for p in equalities for (_,m) in terms(p))
    FP = typeof(main.objective)
    density_real,_ = N._continuous_real_cone(first(densities))
    # Split only exact structural zeros. This preserves the full perturbative
    # matrix condition and does not discard cross-component nonzero entries.
    graph = SimpleGraph(size(perturbative,1))
    for j in axes(perturbative,2), i in 1:j-1
        (!iszero(perturbative[i,j]) || !iszero(perturbative[j,i])) && add_edge!(graph,i,j)
    end
    perturbative_blocks = [first(N._continuous_real_cone(perturbative[c,c])) for c in connected_components(graph)]
    println("SETUP n=$n radius=$radius rows=$(length(basis)) original_moments=$(length(original.words)) reduced_moments=$(length(main.total_basis)) main_blocks=$(length(main.constraints)) main_max=$(maximum(report.psd_block_sizes)) subsystem_sites=$k subsystem_block=$(size(density_real,1)) perturbation_rows=$(length(extra.qs)) perturbative_blocks=$([size(b,1) for b in perturbative_blocks]) nonzero_stationarity=$(length(equalities)) seconds=$(time()-started)")
    flush(stdout)
    results = Dict{Symbol,Any}()
    for mode in (:none,:subsystem,:perturbative,:both)
        constraints = copy(main.constraints)
        if mode in (:subsystem,:both)
            push!(constraints,(:PSD,map(p->convert(FP,p),density_real)))
        end
        if mode in (:perturbative,:both)
            append!(constraints,[(:PSD,map(p->convert(FP,p),b)) for b in perturbative_blocks])
            append!(constraints,[(:Zero,reshape([convert(FP,p)],1,1)) for p in equalities])
        end
        mp = N.MomentProblem(main.objective,constraints,copy(main.total_basis),length(main.total_basis);real_moments=true)
        result = screen(mp,solver)
        results[mode] = result
        println("RESULT n=$n mode=$mode lower=$(result.lower) lower_per_spin=$(result.lower/n) primal=$(result.primal) gap=$(result.gap) violation=$(result.violation) status=$(result.status) accepted=$(result.accepted) build_seconds=$(result.build_seconds) solve_seconds=$(result.solve_seconds)")
        flush(stdout)
        if control
            generic = pauli_strengthened_relaxation(pop,config;
                subsystems=mode in (:subsystem,:both) ? patches : Vector{Int}[],
                perturbations=mode in (:perturbative,:both) ? qs : eltype(basis)[])
            reference = screen(generic.moment_problem,solver)
            @assert result.accepted && reference.accepted
            @assert abs(result.lower-reference.lower) <= 1e-7
            println("CONTROL mode=$mode generic=$(reference.lower) fast=$(result.lower) accepted=true")
        end
        GC.gc()
    end
    monotone = results[:subsystem].lower >= results[:none].lower-1e-7 &&
        results[:perturbative].lower >= results[:none].lower-1e-7 &&
        results[:both].lower >= max(results[:subsystem].lower,results[:perturbative].lower)-1e-7
    accepted = all(r.accepted for r in values(results)) && monotone
    println("COMPLETE n=$n accepted=$accepted monotone=$monotone elapsed=$(time()-started)")
    flush(stdout)
    return accepted
end

function main()
    println("ENV julia=$VERSION cpu=$(first(Sys.cpu_info()).model) julia_threads=$(Threads.nthreads()) blas_threads=$(BLAS.get_num_threads()) JuMP=$(pkgversion(JuMP)) MosekTools=$(pkgversion(MosekTools))")
    println("FORMULATION $(get(ENV,"FORMULATION","moment_variables"))")
    paths = ["src/optimization/"*file*".jl" for file in ("pauli_strengthening","physics_sparsity","pauli_chains")]
    println("SOURCE sha256=$(bytes2hex(sha256(join(read.(paths,String)))))")
    println("CONTRACT geometry_radius=1 perturbation_order=$PERTURBATION_ORDER continuous=false translation=true reflection=true conjugation=true sign=true axis_permutations=true acceptance_gap=1e-7 acceptance_violation=1e-7 solver_tolerance=$SOLVER_TOLERANCE solver_time_limit=120 single_run_diagnostic=true")
    solver = optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_NUM_THREADS"=>1,
        "MSK_IPAR_LOG"=>0,"MSK_DPAR_OPTIMIZER_MAX_TIME"=>120.0,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>SOLVER_TOLERANCE,"MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>SOLVER_TOLERANCE,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>SOLVER_TOLERANCE)
    license = Model(solver)
    set_silent(license)
    @variable(license,v>=1)
    @objective(license,Min,v)
    optimize!(license)
    @assert termination_status(license) == MOI.OPTIMAL
    @assert chain_comparison(2,solver;control=true)
    accepted = true
    for n in (isempty(ARGS) ? [10,20] : parse.(Int,ARGS))
        accepted = chain_comparison(n,solver) && accepted
    end
    accepted || error("At least one solve failed the frozen acceptance gates; inspect retained results.")
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
