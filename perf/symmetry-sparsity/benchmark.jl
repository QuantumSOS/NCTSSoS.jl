using NCTSSoS, LinearAlgebra, JuMP, MosekTools, TOML
include("experiment.jl")
const E=SymmetrySparsityExperiment
const N=NCTSSoS
BLAS.set_num_threads(4)
const route=first(ARGS)
@assert route in ("A","B","C","D","E","F")
const n=length(ARGS)>1 ? parse(Int,ARGS[2]) : 20
const repetitions=length(ARGS)>2 ? parse(Int,ARGS[3]) : 1
const window_width=parse(Int,get(ENV,"WINDOW_WIDTH","4"))
const moment_limit=parse(Int,get(ENV,"MOMENT_LIMIT","100000"))
const numeric_projection=get(ENV,"NUMERIC_PROJECTION","false")=="true"
const report_suffix=get(ENV,"REPORT_SUFFIX","")
const OPT=optimizer_with_attributes(Mosek.Optimizer,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-8,"MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-8,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-8,"MSK_IPAR_NUM_THREADS"=>4,
    "MSK_DPAR_OPTIMIZER_MAX_TIME"=>300.0)
const report=Dict{String,Any}("route"=>route,"n"=>n,"julia"=>string(VERSION),
    "cpu"=>Sys.cpu_info()[1].model,"julia_threads"=>Threads.nthreads(),
    "blas_threads"=>BLAS.get_num_threads(),"form"=>"explicit dual, realified Hermitian PSD",
    "group_enumeration"=>false,"moment_limit"=>moment_limit,
    "numeric_projection"=>numeric_projection,"window_width"=>window_width,"samples"=>Any[])
function checkpoint()
    open(joinpath(@__DIR__,"benchmark-$route-$n$report_suffix.toml"),"w") do io
        TOML.print(io,report)
    end
end
function stage(f,label,sample)
    sample["active_stage"]=label; checkpoint()
    println("BEGIN ",label);flush(stdout)
    result=@timed f()
    push!(sample["stages"],Dict("name"=>label,"seconds"=>result.time,
        "allocated_bytes"=>result.bytes,"gc_seconds"=>result.gctime))
    sample["active_stage"]="between stages";checkpoint()
    println("END ",label," seconds=",result.time," allocated_bytes=",result.bytes);flush(stdout)
    result.value
end
function run(sample)
    reg,ops,H,seed=stage("input",sample) do
        reg,ops=create_pauli_variables(1:n)
        H=heisenberg_chain_hamiltonian(ops;coupling=1/4,periodic=true)
        (reg,ops,H,pauli_contiguous_chain_basis(reg,4))
    end
    gs=stage("symmetry",sample) do
        if route=="A"
            CliffordSymmetry[]
        elseif route=="B"
            x,y,z=ops
            h=CliffordSymmetry(Dict(w=>image for i in 1:n for (w,image) in
                (x[i]=>(1,z[i]),y[i]=>(-1,y[i]),z[i]=>(1,x[i])));nqubits=n)
            s=CliffordSymmetry(Dict(w=>image for i in 1:n for (w,image) in
                (x[i]=>(1,y[i]),y[i]=>(-1,x[i]),z[i]=>(1,z[i])));nqubits=n)
            T=eltype(first(seed).word)
            [pauli_site_permutation([2:n;1];integer_type=T),
             pauli_site_permutation(reverse(1:n);integer_type=T),h,s]
        else
            sympleq_symmetry_spec(H).clifford_generators
        end
    end
    @assert all(g->N._act_polynomial(g,H)==H,gs)
    rows=stage("row_closure",sample) do
        E.close_rows(seed,gs)
    end
    actions=stage("signed_row_actions",sample) do
        E.actions_on(rows,gs)
    end
    sample["rows"]=length(rows);sample["generators"]=length(gs)
    groups=stage("sparsity",sample) do
        if route in ("A","B","C")
            [collect(eachindex(rows))]
        else
            @assert 1<=window_width<=n
            seeds=E.window_groups(seed,rows,n,window_width)
            sample["seed_block_sizes"]=sort!(unique(length.(seeds)))
            parents=route=="D" ? [collect(eachindex(rows))] :
                E.close_groups(seeds,actions)
            if route=="E" && window_width>=4
                # Closure preserves inclusion: check the seed inclusion, then
                # block_orbits checks invariance of the complete collection.
                smaller=E.window_groups(seed,rows,n,max(4,window_width-1))
                @assert all(C->any(D->issubset(C,D),seeds),smaller)
                @assert Set(vcat(parents...))==Set(eachindex(rows))
            end
            sample["parent_blocks"]=length(parents)
            if route=="E"
                parents
            else
                result=E.term_blocks(parents,rows,actions,Set([one(first(rows));monomials(H)]))
                sample["completion_rounds"]=result.rounds
                sample["parent_edge_count"]=result.edge_count
                result.blocks
            end
        end
    end
    sample["blocks_before_reduction"]=length(groups)
    sample["largest_before_reduction"]=maximum(length,groups)
    # Preflight moment storage before paying for numerical decomposition.
    coordinates=stage("moment_coordinates",sample) do
        E.moment_coordinates(H,rows,groups,gs;budget=E.Budget(max_items=moment_limit))
    end
    sample["moment_words"]=length(coordinates.coordinates)
    sample["moment_variables"]=coordinates.nvars
    sample["zero_moment_orbits"]=coordinates.zero_orbits
    orbits=stage("block_orbits_stabilizers",sample) do
        E.block_orbits(groups,actions,length(rows))
    end
    sample["block_orbits"]=length(orbits)
    sample["stabilizer_generators"]=[length(o.stabilizers) for o in orbits]
    plans=stage("decomposition",sample) do
        [(o.representative,E.decompose_signed(o.stabilizers,length(o.representative))) for o in orbits]
    end
    built=stage("symbolic_assembly",sample) do
        E.assemble(H,rows,groups,gs,OPT;reduced=route!="A",dual_form=true,
            coordinates,block_plans=plans,numeric_projection)
    end
    sample["hermitian_block_sizes"]=sort(built.sizes)
    stage("solve",sample) do
        optimize!(built.model)
    end
    sample["termination_status"]=string(termination_status(built.model))
    sample["primal_status"]=string(primal_status(built.model))
    sample["dual_status"]=string(dual_status(built.model));checkpoint()
    optimal=termination_status(built.model)==MOI.OPTIMAL
    allow_feasible=get(ENV,"ALLOW_FEASIBLE","false")=="true"
    (optimal || (allow_feasible && primal_status(built.model)==MOI.FEASIBLE_POINT &&
        dual_status(built.model)==MOI.FEASIBLE_POINT)) || error("solver did not report OPTIMAL")
    stage("verification",sample) do
        mv=built.moment_values()
        valueword(w)=begin
            k,s=coordinates.coordinates[w]; k==0 ? 0.0 : s*mv[k]
        end
        energy=sum(real(c)*valueword(w) for (c,w) in terms(H))
        bound=objective_value(built.model)
        gap=abs(energy-bound)
        sample["energy_from_moments"]=energy;sample["energy_lower_bound"]=bound
        sample["primal_dual_gap"]=gap;checkpoint()
        @assert gap<=2e-6*max(1,abs(energy))
        @assert abs(valueword(one(first(rows)))-1)<1e-7
        residual=0.0
        for (w,_) in coordinates.coordinates, g in gs
            s,v=N._act_monomial(g,w)
            residual=max(residual,abs(valueword(w)-s*valueword(v)))
        end
        @assert residual<1e-8
        # Validate the full original representative matrices, not only their
        # reduced blocks. The other blocks follow by exact signed transport.
        mineig=Inf
        for o in orbits
            C=o.representative
            M=[begin
                phase,w=E.product(rows[i],rows[j]);phase*valueword(w)
            end for i in C,j in C]
            @assert norm(M-M')<1e-8
            mineig=min(mineig,eigmin(Hermitian(M)))
        end
        sample["min_original_block_eigenvalue"]=mineig
        sample["symmetry_residual"]=residual
        violations=primal_feasibility_report(built.model)
        sample["max_conic_dual_violation"]=maximum(values(violations);init=0.0)
        checkpoint()
        @assert mineig>=-2e-6
        @assert sample["max_conic_dual_violation"]<=2e-6
    end
    sample["status"]=optimal ? "verified" : "residual_checks_passed_not_optimal"
    sample["sum_stage_seconds"]=sum(s["seconds"] for s in sample["stages"])
    checkpoint()
end

# The license is checked before relying on the solver, outside route timings.
license=Model(OPT);set_silent(license)
@variable(license,t>=1);@objective(license,Min,t);optimize!(license)
@assert termination_status(license)==MOI.OPTIMAL && abs(value(t)-1)<1e-8
for repetition in 0:repetitions-1
    GC.gc(true)
    sample=Dict{String,Any}("sample"=>repetition,"status"=>"running","stages"=>Any[])
    push!(report["samples"],sample);checkpoint()
    try
        run(sample)
        println("RESULT status=",sample["status"]," sample=",repetition," route=",route," seconds=",sample["sum_stage_seconds"],
            " bound=",sample["energy_lower_bound"]);flush(stdout)
    catch err
        sample["status"]=err isa E.LimitExceeded ? "capped" : "error"
        sample["error"]=sprint(showerror,err);checkpoint();rethrow()
    end
end
