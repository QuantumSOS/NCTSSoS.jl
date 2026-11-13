# Heisenberg chain construction and optional solve, up to the paper's N=100.
# This is a configurable reproduction attempt, NOT a reproduced paper recipe.
# julia --project=<NCTSSoS/JuMP environment> perf/wang_chain_reproduction.jl 20
# Add --solve in an environment with licensed MosekTools. No solve by default.
# FORMULATION=sos (historical default) or moment_variables; compare on the SAME input.
using NCTSSoS, JuMP, LinearAlgebra, TOML, SHA
const N = NCTSSoS
const MOI = JuMP.MOI

# Historical b18830a density blocks, recovered and entry-validated in the frozen
# campaign. Keep this recipe distinct from the public magnetization blocks.
function historical_density8(::Type{P},reduce_word) where {T,P<:Polynomial{PauliAlgebra,T}}
    M = NormalMonomial{PauliAlgebra,T}
    C = N._coefficient_type(P)
    options = [Tuple{Int,C,M}[] for _ in 1:256]
    for initial in 0:4^8-1
        code=initial; word=T[]; counts=zeros(Int,3); flip=0; yz=0
        for site in 1:8
            axis=code&3; code>>=2
            axis==0 && continue
            counts[axis]+=1; push!(word,T(3*(site-1)+axis))
            bit=1<<(8-site)
            axis in (1,2) && (flip|=bit)
            axis in (2,3) && (yz|=bit)
        end
        all(iseven,counts) || continue
        rep=reduce_word(M(word)); isnothing(rep) && continue
        push!(options[flip+1],(yz,(-1)^(counts[2]÷2),rep))
    end
    matrices=Matrix{P}[]
    for residue in (0,1,2)
        states=[s for s in 0:255 if mod(count_ones(s),4)==residue]
        mat=Matrix{P}(undef,length(states),length(states))
        for (j,b) in enumerate(states), (i,a) in enumerate(states)
            mat[i,j]=P([(isodd(count_ones(b&yz)) ? -phase : phase,m)
                for (yz,phase,m) in options[xor(a,b)+1]])
        end
        push!(matrices,mat)
    end
    return matrices
end

function wang_chain_problem(n; pair_distance=min(10,n÷2), density_sites=min(8,n),
    five_site_window=n÷2, seven_site_window=min(7,n), perturbative=false,
    density_kind=density_sites == 8 ? :historical : :magnetization,
    normalize_density=false)
    n >= 4 || throw(ArgumentError("This runner requires at least four spins."))
    0 <= density_sites <= n || throw(ArgumentError("density_sites must lie in 0:nsites."))
    density_kind in (:historical,:magnetization) || throw(ArgumentError("Unknown density_kind"))
    density_kind == :historical && density_sites != 8 && throw(ArgumentError("Historical density blocks require eight sites."))
    reg,ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops;coupling=1//4),reg)
    basis = physics_moment_basis(pop,4;pair_distance)
    patches = density_sites == 0 ? Vector{Int}[] : [collect(1:density_sites)]
    t = time()
    selected = pauli_chain_stationarity_operators(pop,ops;basis,subsystems=patches,
        five_site_window,seven_site_window)
    selection_seconds = time()-t
    t = time()
    main,report = pauli_translation_invariant_moment_relaxation(pop,ops,4;
        basis,axis_permutation_symmetry=true,
        rdm_levels=density_sites == 0 || density_kind == :historical ? Int[] : [density_sites],
        state_optimality=perturbative ? :linear_psd : :none)
    FP = typeof(main.objective)
    constraints = copy(main.constraints)
    if density_kind == :historical
        reducer = N._PauliMomentReducer{eltype(basis)}(n,true,true,true)
        reduce_word(m) = iszero(N._pauli_sign_signature(m)) ? N._reduce_moment(reducer,m) : nothing
        densities = historical_density8(FP,reduce_word)
        append!(constraints,[(:PSD,normalize_density ? mat/256 : mat) for mat in densities])
    elseif normalize_density
        for (i,label) in enumerate(report.block_labels)
            haskey(label,:rdm) || continue
            cone,mat = constraints[i]
            constraints[i] = (cone,mat/2^label.rdm)
        end
    end
    matrix_seconds = time()-t
    # Reuse the exact pre-cancellation domain and selected equations. Matrix
    # support is audited by SET equality, never taken as the original domain.
    domain = Set(selected.moments)
    matrix_words = Set(m for (_,mat) in constraints for p in mat for m in monomials(p))
    perturbative || matrix_words == domain || error(
        "Pre-cancellation domain differs from matrix support; inspect exact sets before comparing formulations.")
    words = sort!(collect(union(domain,matrix_words)))
    append!(constraints,[(:Zero,reshape([convert(FP,p)],1,1)) for p in selected.equations])
    mp = N.MomentProblem(main.objective/n,constraints,words,length(words);real_moments=true)
    return (;moment_problem=mp,report,selection=selected,selection_seconds,matrix_seconds)
end

function wang_chain(n; solve=false, output="", pair_distance=min(10,n÷2),
    density_sites=min(8,n), five_site_window=n÷2, seven_site_window=min(7,n),
    perturbative=false, density_kind=density_sites == 8 ? :historical : :magnetization,
    normalize_density=false, formulation=:sos, solver_threads=4,
    solver_seconds=900.0, solver_tolerance=1e-8)
    formulation in (:moment_variables,:sos) || throw(ArgumentError("Unknown formulation $formulation"))
    BLAS.set_num_threads(1)
    record = Dict{String,Any}("n_sites"=>n,"degree"=>4,"pair_distance"=>pair_distance,
        "density_sites"=>density_sites,"five_site_window"=>five_site_window,
        "seven_site_window"=>seven_site_window,"perturbative"=>perturbative,
        "formulation"=>string(formulation),"density_kind"=>string(density_kind),
        "normalized_density"=>normalize_density,"normalization"=>"energy per spin",
        "recipe"=>"configurable attempt, not established paper table recipe",
        "acceptance_contract"=>"numerical screen only; not an exact certificate or paper-quality match",
        "peak_scope"=>"entire Julia process, including previous runs in the same process",
        "julia_version"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(),"phase"=>"construction")
    save() = isempty(output) ? nothing : open(io->TOML.print(io,record;sorted=true),output,"w")
    save()
    started = time()
    timed = @timed wang_chain_problem(n;pair_distance,density_sites,five_site_window,
        seven_site_window,perturbative,density_kind,normalize_density)
    result = timed.value
    record["cumulative_allocations_bytes"] = timed.bytes
    record["process_peak_rss_bytes"] = Sys.maxrss()
    mp = result.moment_problem
    record["selection_seconds"] = result.selection_seconds
    record["matrix_seconds"] = result.matrix_seconds
    record["construction_seconds"] = time()-started
    record["rows"] = result.report.basis_size
    record["moments"] = length(mp.total_basis)
    record["candidates"] = result.selection.candidate_count
    record["stationarity_equations"] = length(result.selection.equations)
    record["block_sides"] = [size(mat,1) for (cone,mat) in mp.constraints if cone == :PSD]
    record["matrix_nonzero_terms"] = sum(length(terms(p)) for (_,mat) in mp.constraints for p in mat)
    record["moments_sha256"] = bytes2hex(sha256(join((join(m.word,",") for m in result.selection.moments),"\n")))
    record["equations_sha256"] = bytes2hex(sha256(join((string(terms(p)) for p in result.selection.equations),"\n")))
    record["phase"] = "constructed"
    save(); println("CONSTRUCTION ",record); flush(stdout)
    solve || return record

    # Verify the license with a trivial solve before constructing the large model.
    isdefined(@__MODULE__,:MosekTools) || throw(ArgumentError("Load MosekTools before calling wang_chain with solve=true."))
    optimizer = optimizer_with_attributes(MosekTools.Optimizer,
        "MSK_IPAR_NUM_THREADS"=>solver_threads,
        "MSK_DPAR_OPTIMIZER_MAX_TIME"=>solver_seconds,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>solver_tolerance,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>solver_tolerance,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>solver_tolerance)
    license = Model(optimizer)
    set_silent(license)
    @variable(license,v>=1)
    @objective(license,Min,v)
    optimize!(license)
    termination_status(license) == MOI.OPTIMAL || error("MOSEK license check failed")
    record["solver_threads"] = solver_threads
    record["solver_tolerance"] = solver_tolerance
    record["solver_time_limit"] = solver_seconds
    record["phase"] = "model_construction"; save()
    GC.gc()
    t = time()
    model = formulation == :sos ? N.sos_dualize(mp).model : first(build_jump_model(mp))
    record["model_seconds"] = time()-t
    set_optimizer(model,optimizer)
    record["phase"] = "solve"; save()
    t = time(); optimize!(model)
    record["solve_call_seconds"] = time()-t
    record["solver_seconds"] = solve_time(model)
    record["termination"] = string(termination_status(model))
    record["raw_status"] = raw_status(model)
    record["primal_status"] = string(primal_status(model))
    record["dual_status"] = string(dual_status(model))
    record["accepted"] = false
    if has_values(model) && has_duals(model)
        primal,dual = objective_value(model),dual_objective_value(model)
        record["primal_per_spin"],record["dual_per_spin"] = primal,dual
        record["objective_sense"] = string(objective_sense(model))
        record["lower_objective_per_spin"] = objective_sense(model) == MOI.MIN_SENSE ? dual : primal
        record["relative_gap"] = abs(primal-dual)/max(1.0,abs(primal),abs(dual))
        violations = primal_feasibility_report(model;atol=0.0)
        eqmax,psdmax = 0.0,0.0
        for (F,S) in list_of_constraint_types(model), c in all_constraints(model,F,S)
            object = constraint_object(c)
            if object.set isa MOI.EqualTo
                eqmax = max(eqmax,abs(value(object.func)-object.set.value))
            elseif object.set isa MOI.Zeros
                eqmax = max(eqmax,maximum(abs,value.(object.func);init=0.0))
            elseif object.set isa Union{MOI.PositiveSemidefiniteConeTriangle,MOI.PositiveSemidefiniteConeSquare}
                psdmax = max(psdmax,get(violations,c,0.0))
            else
                error("Unscreened constraint type $S")
            end
        end
        record["equality_coordinate_max"] = eqmax
        record["psd_distance_max"] = psdmax
        record["accepted"] = termination_status(model) == MOI.OPTIMAL &&
            primal_status(model) == dual_status(model) == MOI.FEASIBLE_POINT &&
            record["relative_gap"] <= 1e-7 && eqmax <= 1e-7 && psdmax <= 1e-7
    end
    record["phase"] = "complete"
    record["total_seconds"] = time()-started
    record["process_peak_rss_bytes"] = Sys.maxrss()
    save(); println("RESULT ",record)
    return record
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    "--solve" in ARGS && (@eval using MosekTools)
    n = parse(Int,ARGS[1])
    k = parse(Int,get(ENV,"DENSITY_SITES",string(min(8,n))))
    wang_chain(n;solve="--solve" in ARGS,
        output=get(ENV,"OUTPUT",""),
        pair_distance=parse(Int,get(ENV,"PAIR_DISTANCE",string(min(10,n÷2)))),
        density_sites=k,
        density_kind=Symbol(get(ENV,"DENSITY_KIND",k == 8 ? "historical" : "magnetization")),
        normalize_density=get(ENV,"NORMALIZE_DENSITY","0") == "1",
        five_site_window=parse(Int,get(ENV,"FIVE_WINDOW",string(parse(Int,ARGS[1])÷2))),
        perturbative=get(ENV,"PERTURBATIVE","0") == "1",
        formulation=Symbol(get(ENV,"FORMULATION","sos")),
        solver_threads=parse(Int,get(ENV,"SOLVER_THREADS","4")),
        solver_seconds=parse(Float64,get(ENV,"SOLVER_SECONDS","900")))
end
