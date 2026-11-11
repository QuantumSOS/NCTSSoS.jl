# Same Hamiltonians, geometry, patches and perturbations as the discrete-only
# comparison. This adds full continuous spin invariance and matrix reduction.
include("pauli_strengthening_comparison.jl")

function complete_comparison()
    formulation = Symbol(get(ENV,"FORMULATION","sos"))
    solver = optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_NUM_THREADS"=>1,
        "MSK_IPAR_LOG"=>0,"MSK_DPAR_OPTIMIZER_MAX_TIME"=>120.0,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>SOLVER_TOLERANCE,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>SOLVER_TOLERANCE,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>SOLVER_TOLERANCE)
    license = Model(solver)
    set_silent(license)
    @variable(license,v>=1)
    @objective(license,Min,v)
    optimize!(license)
    @assert termination_status(license) == MOI.OPTIMAL
    println("ENV julia=$VERSION cpu=$(first(Sys.cpu_info()).model) julia_threads=$(Threads.nthreads()) blas_threads=$(BLAS.get_num_threads()) JuMP=$(pkgversion(JuMP)) MosekTools=$(pkgversion(MosekTools))")
    println("CONTRACT continuous_spin=true translation=true reflection=true conjugation=true global_spin_cliffords=contained_in_continuous_group exhaustive_clifford_discovery=false radius=1 subsystem_sites=4 perturbations=single_and_adjacent_pair unchanged_acceptance=1e-7 solver_tolerance=$SOLVER_TOLERANCE")
    println("SOURCE sha256=$(bytes2hex(sha256(read("src/optimization/pauli_strengthening.jl"))))")
    println("FORMULATION $formulation")
    baseline = Dict{Tuple{Int,String},Float64}()
    if haskey(ENV,"BASELINE_CSV")
        for row in Iterators.drop(eachline(ENV["BASELINE_CSV"]),1)
            fields = split(row,',')
            baseline[parse(Int,fields[1]),fields[2]] = parse(Float64,fields[3])
        end
    end
    accepted_all = true
    for n in (isempty(ARGS) ? [10,20] : parse.(Int,ARGS))
        started = time()
        reg,ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops),reg)
        basis = physics_moment_basis(reg,cycle_graph(n);radius=1,max_rows=2000)
        patches = [sort(mod1.(i:i+3,n)) for i in 1:n]
        qs = [reduce(vcat,ops); [only(monomials(a[i]*b[mod1(i+1,n)])) for i in 1:n for a in ops for b in ops]]
        full = pauli_chain_strengthened_relaxation(pop,ops;basis,subsystems=patches,perturbations=qs)
        sizes(kind) = [size(mat,1) for ((_,mat),label) in zip(full.moment_problem.constraints,full.provenance) if label==kind]
        println("SETUP n=$n rows=$(length(basis)) discrete_moments=$(length(full.spin.index)) spin_invariant_moments=$(length(full.spin.words)) coefficient_residual=$(full.coefficient_error) main_blocks=$(sizes(:moment)) subsystem_blocks=$(sizes(:subsystem)) perturbative_blocks=$(sizes(:perturbative)) seconds=$(time()-started)")
        flush(stdout)
        results = Dict{Symbol,Any}()
        for mode in (:none,:subsystem,:perturbative,:both)
            kinds = [:moment]
            mode in (:subsystem,:both) && push!(kinds,:subsystem)
            mode in (:perturbative,:both) && append!(kinds,[:perturbative,:stationarity])
            constraints = [c for (c,kind) in zip(full.moment_problem.constraints,full.provenance) if kind in kinds]
            mp = N.MomentProblem(full.moment_problem.objective,constraints,
                copy(full.moment_problem.total_basis),length(full.moment_problem.total_basis);real_moments=true)
            result = screen(mp,solver;formulation)
            results[mode] = result
            old = get(baseline,(n,String(mode)),NaN)
            difference = result.lower/n-old
            equivalent = isnan(old) || abs(difference*n) <= 1e-7
            println("RESULT n=$n mode=$mode lower=$(result.lower) lower_per_spin=$(result.lower/n) gap=$(result.gap) violation=$(result.violation) status=$(result.status) accepted=$(result.accepted) previous_per_spin=$old difference_per_spin=$difference equivalent=$equivalent build_seconds=$(result.build_seconds) solve_seconds=$(result.solve_seconds)")
            flush(stdout)
            accepted_all &= result.accepted && equivalent
        end
        monotone = results[:subsystem].lower >= results[:none].lower-1e-7 &&
            results[:perturbative].lower >= results[:none].lower-1e-7 &&
            results[:both].lower >= max(results[:subsystem].lower,results[:perturbative].lower)-1e-7
        accepted_all &= monotone
        println("COMPLETE n=$n monotone=$monotone elapsed=$(time()-started)")
        flush(stdout)
    end
    accepted_all || error("A solve or same-relaxation equivalence gate failed; do not accept the campaign.")
end
complete_comparison()
