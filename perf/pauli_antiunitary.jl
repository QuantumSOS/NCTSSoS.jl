# Run each mode in a fresh process. Requires MosekTools in the active environment.
# julia --project=<benchmark-env> perf/pauli_antiunitary.jl MODE N REPEATS
# MODE: complex, antiunitary, complex-ts, antiunitary-ts, specialized
using NCTSSoS, JuMP, MosekTools, LinearAlgebra, Statistics, Printf

mode = get(ARGS, 1, "antiunitary")
n = parse(Int, get(ARGS, 2, "4"))
repeats = parse(Int, get(ARGS, 3, "3"))
mode in ("complex", "antiunitary", "complex-ts", "antiunitary-ts", "specialized") || error("Unknown mode")
BLAS.set_num_threads(1)
registry, ops = create_pauli_variables(1:n)
pop = polyopt(heisenberg_chain_hamiltonian(ops; coupling=1/4, periodic=true), registry)
basis = pauli_contiguous_chain_basis(ops, 2; periodic=true)
optimizer = optimizer_with_attributes(Mosek.Optimizer, "MSK_IPAR_NUM_THREADS"=>1,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-9, "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-9,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-9)

# Independent dense spin-matrix oracle, not another relaxation.
paulis = ([0 1; 1 0], [0 -im; im 0], [1 0; 0 -1])
identity2 = Matrix{ComplexF64}(I,2,2)
hamiltonian = zeros(ComplexF64,2^n,2^n)
for i in 1:n, p in paulis
    j = mod1(i+1,n)
    hamiltonian .+= reduce(kron,[k == i || k == j ? p : identity2 for k in 1:n]) / 4
end
exact = eigmin(Hermitian(hamiltonian))

function construct()
    reduction = nothing
    mp = if mode == "specialized"
        first(pauli_translation_invariant_moment_relaxation(pop,ops,2; basis=basis))
    else
        config = SolverConfig(optimizer=nothing,moment_basis=basis,
            ts_algo=endswith(mode,"-ts") ? MF() : NoElimination())
        sparsity = compute_sparsity(pop,config)
        if startswith(mode,"antiunitary")
            reduction = antiunitary_moment_relaxation(pop,sparsity)
            reduction.moment_problem
        else
            NCTSSoS.moment_relax(pop,sparsity.corr_sparsity,sparsity.cliques_term_sparsities)
        end
    end
    model = NCTSSoS.sos_dualize(mp).model
    set_optimizer(model,optimizer)
    set_silent(model)
    return mp,model,reduction
end

function validate(mp,model,reduction)
    termination_status(model) == MOI.OPTIMAL || error("Nonoptimal solve: $(termination_status(model))")
    cones = [c for (F,S) in list_of_constraint_types(model) if S <: MOI.PositiveSemidefiniteConeTriangle
        for c in all_constraints(model,F,S)]
    certificate = zero(mp.objective)
    mineig = Inf
    sides = Int[]
    for (b,cone) in enumerate(cones)
        Z = Matrix(JuMP.reshape_vector(value.(constraint_object(cone).func),JuMP.shape(constraint_object(cone))))
        push!(sides,size(Z,1))
        mineig = min(mineig,eigmin(Symmetric(Z)))
        idx = mp.linear.psd_block_constraint_idx[b]
        mat = mp.constraints[idx][2]
        m = size(mat,1)
        if !isnothing(reduction)
            Q = recover_antiunitary_gram(reduction,b,Z)
            mat = reduction.original.constraints[reduction.original.linear.psd_block_constraint_idx[b]][2]
            certificate += sum(Q[j,i]*mat[i,j] for i in 1:m,j in 1:m)
        elseif size(Z,1) == 2m
            # Same lifted-dual convention as sos_dualize's coefficient equations.
            Q = Z[1:m,1:m]+Z[m+1:end,m+1:end]+im*(Z[m+1:end,1:m]-Z[1:m,m+1:end])
            certificate += sum(Q[i,j]*mat[i,j] for i in 1:m,j in 1:m)
        else
            certificate += sum(Z[i,j]*mat[i,j] for i in 1:m,j in 1:m)
        end
    end
    residual = maximum(abs,coefficients(mp.objective-objective_value(model)-certificate);init=0.0)
    gap = abs(objective_value(model)-dual_objective_value(model))
    residual < 2e-7 || error("Certificate residual $residual")
    mineig >= -2e-7 || error("Non-PSD certificate $mineig")
    objective_value(model) <= exact+2e-7 || error("Invalid lower bound")
    gap < 2e-7 || error("Duality gap $gap")
    return (;residual,mineig,gap,sides,bound=objective_value(model))
end

function sample()
    GC.gc(true)
    build = @timed construct()
    mp,model,reduction = build.value
    elapsed = @elapsed optimize!(model)
    proof = validate(mp,model,reduction)
    return (;build=build.time,allocated=build.bytes,solve_wall=elapsed,
        solver=solve_time(model),proof)
end

warm = sample() # compile and validate the exact workload before measuring
samples = [sample() for _ in 1:repeats]
println("Julia=$VERSION CPU=$(Sys.CPU_NAME) threads=$(Threads.nthreads()) BLAS_threads=$(BLAS.get_num_threads())")
println("mode=$mode n=$n order=2 basis=$(length(basis)) repeats=$repeats exact=$exact")
println("cold_build=$(warm.build) cold_solve_wall=$(warm.solve_wall)")
for (i,s) in enumerate(samples)
    println("sample=$i build=$(s.build) allocated=$(s.allocated) solver=$(s.solver) solve_wall=$(s.solve_wall) bound=$(s.proof.bound) residual=$(s.proof.residual) mineig=$(s.proof.mineig) gap=$(s.proof.gap)")
end
for field in (:build,:allocated,:solver,:solve_wall)
    values = [getproperty(s,field) for s in samples]
    println("$field median=$(median(values)) min=$(minimum(values)) max=$(maximum(values))")
end
println("cone_sides=$(samples[end].proof.sides) peak_rss_bytes=$(Sys.maxrss())")
