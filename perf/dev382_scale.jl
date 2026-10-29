# Run in the external benchmark environment; see dev382_scale_plan.md.
using NCTSSoS, JuMP, MosekTools, LinearAlgebra, TOML

const MODE, FAMILY, OUTPUT = ARGS[1], ARGS[2], ARGS[5]
const N, REPEATS = parse(Int, ARGS[3]), parse(Int, ARGS[4])
MODE in ("complex", "antiunitary", "complex-ts", "antiunitary-ts", "specialized") || error("Unknown mode")
FAMILY in ("heisenberg", "inhomogeneous") || error("Unknown family")
BLAS.set_num_threads(1)
const OPTIMIZER = optimizer_with_attributes(MosekTools.Optimizer,
    "MSK_IPAR_NUM_THREADS" => 4, "MSK_DPAR_OPTIMIZER_MAX_TIME" => 120.0,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8)

function construct()
    registry, ops = create_pauli_variables(1:N)
    x, y, z = ops
    h = if FAMILY == "heisenberg"
        heisenberg_chain_hamiltonian(ops; coupling=1/4, periodic=true)
    else
        # Real, non-translation-invariant nearest-neighbor Hamiltonian.
        sum((1 + (i % 5)/10) * (x[i]*x[mod1(i+1,N)] +
            0.7y[i]*y[mod1(i+1,N)] + 1.2z[i]*z[mod1(i+1,N)])/4 +
            0.13x[i] + (0.07 + (i % 3)/100)*z[i] for i in 1:N)
    end
    pop = polyopt(h, registry)
    basis = pauli_contiguous_chain_basis(ops, 2; periodic=true)
    reduction = nothing
    mp = if MODE == "specialized"
        FAMILY == "heisenberg" || error("Specialized comparator requires translation invariance")
        first(pauli_translation_invariant_moment_relaxation(pop, ops, 2; basis))
    else
        config = SolverConfig(optimizer=nothing, moment_basis=basis,
            ts_algo=endswith(MODE,"-ts") ? MF() : NoElimination())
        sparsity = compute_sparsity(pop, config)
        if startswith(MODE,"antiunitary")
            reduction = antiunitary_moment_relaxation(pop, sparsity)
            reduction.moment_problem
        else
            NCTSSoS.moment_relax(pop,sparsity.corr_sparsity,sparsity.cliques_term_sparsities)
        end
    end
    model = NCTSSoS.sos_dualize(mp).model
    set_optimizer(model, OPTIMIZER)
    set_silent(model)
    return (;mp, model, reduction, basis_words=length(basis))
end

function validate(problem)
    (;mp,model,reduction) = problem
    cones = [c for (F,S) in list_of_constraint_types(model) if S <: MOI.PositiveSemidefiniteConeTriangle
        for c in all_constraints(model,F,S)]
    # Accumulate coefficients once: repeated polynomial addition is quadratic.
    residual = Dict(m => ComplexF64(c) for (c,m) in terms(mp.objective-objective_value(model)))
    mineig, psd_penalty = Inf, 0.0
    sides = Int[]
    for (b,cone) in enumerate(cones)
        object = constraint_object(cone)
        Z = Matrix(JuMP.reshape_vector(value.(object.func), JuMP.shape(object)))
        push!(sides,size(Z,1))
        λ = eigmin(Symmetric(Z))
        mineig = min(mineig,λ)
        mat = mp.constraints[mp.linear.psd_block_constraint_idx[b]][2]
        m = size(mat,1)
        Q = if !isnothing(reduction)
            mat = reduction.original.constraints[reduction.original.linear.psd_block_constraint_idx[b]][2]
            transpose(recover_antiunitary_gram(reduction,b,Z))
        elseif size(Z,1) == 2m
            Z[1:m,1:m]+Z[m+1:end,m+1:end]+im*(Z[m+1:end,1:m]-Z[1:m,m+1:end])
        else
            Z
        end
        # A diagonal PSD repair costs the diagonal polynomial's coefficient l1
        # norm; the doubled embedding contributes twice that diagonal.
        diagonal_l1 = sum(sum(abs,coefficients(mat[i,i]);init=0.0) for i in 1:m)
        psd_penalty += max(0,-λ)*diagonal_l1*(size(Z,1)==2m ? 2 : 1)
        for i in 1:m, j in 1:m, (c,mono) in terms(mat[i,j])
            residual[mono] = get(residual,mono,0.0im) - Q[i,j]*c
        end
    end
    l1 = sum(abs,values(residual);init=0.0)
    maxres = maximum(abs,values(residual);init=0.0)
    gap = abs(objective_value(model)-dual_objective_value(model))
    tolerance = 1e-5*N
    valid = maxres <= 1e-6 && l1+psd_penalty <= tolerance && gap <= tolerance && mineig >= -1e-6
    return Dict("valid"=>valid,"residual_max"=>maxres,"residual_l1"=>l1,
        "gram_min_eigenvalue"=>mineig,"psd_repair_penalty"=>psd_penalty,
        "gap"=>gap,"bound"=>objective_value(model),
        "numerically_repaired_lower"=>objective_value(model)-l1-psd_penalty,
        "cone_sides"=>sides)
end

data = Dict{String,Any}("mode"=>MODE,"family"=>FAMILY,"n"=>N,"julia"=>string(VERSION),
    "cpu"=>Sys.cpu_info()[1].model,"julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
    "solver_threads"=>4,"solver_limit_s"=>120,"tolerance"=>1e-8,
    "git_commit"=>readchomp(`git rev-parse HEAD`),"samples"=>Any[])
function save()
    data["peak_rss_bytes"] = Sys.maxrss()
    open(OUTPUT,"w") do io
        TOML.print(io,data;sorted=true)
    end
end
save()
for sample in 0:REPEATS
    GC.gc(true)
    row = Dict{String,Any}("sample"=>sample,"cold"=>sample==0,"stage"=>"construction")
    push!(data["samples"],row); save()
    build = @timed construct()
    problem = build.value
    model = problem.model
    merge!(row,Dict("construction_s"=>build.time,"construction_bytes"=>build.bytes,
        "basis_words"=>problem.basis_words,"moments"=>length(problem.mp.linear.moments),
        "variables"=>num_variables(model),"stage"=>"solve")); save()
    row["solve_wall_s"] = @elapsed optimize!(model)
    row["solver_s"] = solve_time(model)
    row["status"] = string(termination_status(model))
    row["stage"] = "validation"; save()
    if termination_status(model) != MOI.OPTIMAL
        row["valid"] = false
        row["stage"] = "solver_failed"; save()
        println(row); break
    end
    row["validation_s"] = @elapsed merge!(row,validate(problem))
    row["stage"] = "complete"
    data["solver_version"] = MOI.get(backend(model),MOI.SolverVersion())
    save()
    println(row); flush(stdout)
    row["valid"] || error("Certificate validation failed")
    # Do not keep the previous solver/model alive in the next sample.
    problem = nothing; model = nothing; build = nothing
end
