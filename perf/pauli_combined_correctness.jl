# Small same-relaxation checks, not a timing benchmark or a paper-scale claim.
# Run in an environment containing NCTSSoS, JuMP and a licensed MosekTools.
using NCTSSoS, JuMP, MosekTools, LinearAlgebra, SHA, Test

function rotations(reg, ops)
    vectors = [[axis[i] for axis in ops] for i in eachindex(first(ops))]
    zero_poly = zero(Polynomial(first(first(ops))))
    return [SuppliedDerivation(reg, Dict(v[i].word[1] =>
        (i == a ? Polynomial(v[b]) : i == b ? -Polynomial(v[a]) : zero_poly)
        for v in vectors for i in 1:3)) for a in 1:3 for b in a+1:3]
end

function exact_energy(n, couplings)
    identity = Matrix{ComplexF64}(I, 2, 2)
    paulis = ([0 1; 1 0], [0 -im; im 0], [1 0; 0 -1])
    h = zeros(ComplexF64, 2^n, 2^n)
    for (i,c) in enumerate(couplings), p in paulis
        j = mod1(i+1,n)
        h += c * reduce(kron, [k == i || k == j ? p : identity for k in 1:n])
    end
    return eigmin(Hermitian(h))
end

solver = optimizer_with_attributes(MosekTools.Optimizer,
    "MSK_IPAR_NUM_THREADS"=>1, "MSK_IPAR_LOG"=>0,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-10,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-10,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-10)

files = ["src/optimization/"*name*".jl" for name in
    ("continuous_symmetry", "continuous_reduction", "physics_sparsity", "interface", "moment_linear", "sos")]
println("# Julia=", VERSION, " CPU=", first(Sys.cpu_info()).model,
    " Julia_threads=", Threads.nthreads(), " BLAS_threads=", BLAS.get_num_threads(),
    " MosekTools=", pkgversion(MosekTools), " JuMP=", pkgversion(JuMP))
println("# source_sha256=", bytes2hex(sha256(join(read.(files, String)))))
println("case,n,region_sites,path,objective,dual_objective,exact_energy,original_psd_violation,relative_gap,blocks,max_block,moments,elapsed_s,termination,primal_status,dual_status")

@testset "Small correctness: $label" for (label,n,couplings) in [("ring4",4,fill(1/4,4)), ("ring6",6,fill(1/4,6)),
    ("ring8",8,fill(1/4,8)), ("unequal4",4,[1/4,1/2,3/4,1/2])]
    reg, ops = create_pauli_variables(1:n)
    h = sum(couplings[i]*sum(axis[i]*axis[mod1(i+1,n)] for axis in ops) for i in 1:n)
    pop = polyopt(h,reg)
    actions = rotations(reg,ops)
    zero_actions = [SuppliedDerivation(reg,Dict(i=>zero(h) for i in indices(reg)))]
    finite = label == "unequal4" ? [first(pauli_global_axis_rotation_generators(ops))] :
        [pauli_site_permutation([2:n;1]), pauli_site_permutation(reverse(1:n))]
    metric = Matrix{Int}(I,3n,3n)
    config = SolverConfig(optimizer=solver,order=2)
    exact = exact_energy(n,couplings)
    previous = -Inf
    for sites in (2,3)
        physics = PhysicsSparsity(strategy=:connected,max_region_sites=sites,max_block_size=37)
        # Every mode uses these identical selected blocks. Selection is outside
        # the times below, so they must not be presented as end-to-end benchmarks.
        sparsity = compute_sparsity(pop,config,physics;actions,clifford_generators=finite)
        original = NCTSSoS.moment_relax(pop,sparsity.corr_sparsity,sparsity.cliques_term_sparsities)
        reference = NaN
        for mode in (:plain,:finite,:continuous,:combined)
            started = time_ns()
            if mode == :plain
                # Use the existing matrix-variable lowering uniformly. Scalar
                # moment lowering stalled on the unequal-coupling control.
                model, extract = build_jump_model(original;
                    formulation=:psd_blocks,representation=:complex)
                set_optimizer(model,solver)
                set_silent(model)
                optimize!(model)
                moment_values = extract()
                key_type = typeof(original.linear.identity)
                matrices = [map(p -> sum(c*moment_values[NCTSSoS._moment_key(key_type,m)]
                    for (c,m) in terms(p); init=0.0), mat) for (cone,mat) in original.constraints if cone in (:PSD,:HPSD)]
                block_sizes = [b.size for b in original.linear.psd_blocks_lin]
                nmoments = length(original.linear.moments)
            else
                result = continuous_nctssos(pop,config,mode == :finite ? zero_actions : actions;
                    compact_metric=metric,physics,
                    clifford_generators=mode == :continuous ? CliffordSymmetry[] : finite)
                @assert result.reduction.sparsity.corr_sparsity.clq_mom_mtx_bases == sparsity.corr_sparsity.clq_mom_mtx_bases
                model = result.model
                matrices = collect(values(result.matrices))
                block_sizes = [b.size for b in result.reduction.moment_problem.linear.psd_blocks_lin]
                nmoments = length(result.reduction.moment_problem.linear.moments)
            end
            elapsed = (time_ns()-started)/1e9
            objective, bound = objective_value(model), dual_objective_value(model)
            violation = max(0.0,-minimum(eigmin(Hermitian(ComplexF64.(a))) for a in matrices))
            gap = abs(objective-bound)/max(1.0,abs(objective),abs(bound))
            println(join((label,n,sites,mode,objective,bound,exact,violation,gap,
                length(block_sizes),maximum(block_sizes),nmoments,elapsed,
                termination_status(model),primal_status(model),dual_status(model)),','))
            flush(stdout)
            @test termination_status(model) == JuMP.MOI.OPTIMAL
            @test primal_status(model) == JuMP.MOI.FEASIBLE_POINT && dual_status(model) == JuMP.MOI.FEASIBLE_POINT
            @test violation <= 1e-7 && gap <= 1e-7
            @test objective <= exact+1e-7
            if mode == :plain
                reference = objective
                @test objective >= previous-1e-7
                previous = objective
            else
                @test abs(objective-reference) <= 1e-7
            end
        end
    end
end
