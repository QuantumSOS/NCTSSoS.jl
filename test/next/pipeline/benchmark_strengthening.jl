# Run each kind in a fresh process: julia --project benchmark_strengthening.jl density|soc
# Requires MosekTools and a working license. Construction and package loading
# are outside the timer; the first build includes JIT. No solves are timed.
using NCTSSoS, JuMP, MosekTools, LinearAlgebra, Statistics
const Nx=NCTSSoS.Next

function strengthening_measurement(name)
    name in ("density","soc") || error("Choose density or soc.")
    registry,(x,y,z)=Nx.pauli_variables(3)
    H=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
    problem=Nx.Problem(H)
    W=Nx.word_type(registry,4)
    symmetry=Nx.FamilySymmetry([Nx.SitePermutation([1,3,2])];moment_only=(Nx.AntiunitaryPauli(Nx.pauli_letter(W,1,1)),))
    kind=name=="density" ? Nx.DensityMatrix(1:3) : Nx.InvariantCone(1:3,"YII")
    description=Nx.BlockRelaxation(order=1,strengthen=(kind,),symmetry=symmetry)
    plain=Nx.BlockRelaxation(order=1,symmetry=symmetry)
    println("ENV Julia=",VERSION," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads()," CPU=",Sys.cpu_info()[1].model)
    cold=@timed Base.invokelatest(Nx.build,problem,description)
    println("COLD ",name," seconds=",cold.time," bytes=",cold.bytes)
    report=cold.value.sizes.closure.strengthening
    println("ADDED ",name," blocks=",report.added_blocks," equalities=",report.added_equalities,
        " sides=",[(b.kind,b.side) for b in cold.value.program.blocks[report.first_block:end]],
        " reconstruction_error=",report.reconstruction_error)
    Nx.build(problem,plain)
    samples=[Tuple{Float64,Int}[] for _ in 1:2]
    for pair in 1:7,i in (isodd(pair) ? (1,2) : (2,1))
        GC.gc(); m=@timed Nx.build(problem,i==1 ? plain : description)
        push!(samples[i],(m.time,m.bytes))
    end
    for (label,values) in zip(("plain",name),samples)
        println("SAMPLES ",label," ",values)
        println("WARM ",label," seconds=",median(first.(values))," IQR=",quantile(first.(values),[0.25,0.75]),
            " bytes=",median(last.(values)))
    end
    optimizer=optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_LOG"=>0)
    check=Model(optimizer)
    @variable(check,v>=1); @objective(check,Min,v); optimize!(check)
    @assert termination_status(check)==JuMP.MOI.OPTIMAL && abs(objective_value(check)-1)<1e-8
    for (label,r) in (("plain",plain),(name,description))
        result=Nx.solve(problem,r;solver=Nx.JuMPBackend(optimizer))
        println("BOUND ",label," value=",result.bound," status=",result.status)
    end
    paulis=Dict('I'=>ComplexF64[1 0;0 1],'X'=>ComplexF64[0 1;1 0],
        'Y'=>ComplexF64[0 -im;im 0],'Z'=>ComplexF64[1 0;0 -1])
    exact=eigmin(Hermitian(-sum(foldl(kron,(paulis[c] for c in label)) for label in ("XYI","XIY","IXX","IYY","IZZ"))))
    println("EXACT value=",exact)

    registry,(x,y,z)=create_pauli_variables(1:3)
    h=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
    pop=polyopt(h,registry)
    config=SolverConfig(order=1,optimizer=optimizer)
    if name=="density"
        old=pauli_strengthened_relaxation(pop,config;subsystems=[[1,2,3]],stationarity=false,moment_domain=:extend)
        model,_=build_jump_model(old.moment_problem)
        set_optimizer(model,optimizer)
    else
        old=pauli_invariant_soc_nctssos(pop,config;max_rounds=0)
        model=old.model
        add_pauli_invariant_soc!(model,old.moment_problem,[1,2,3],"YII")
    end
    optimize!(model)
    println("LEGACY ",name," value=",objective_value(model)," status=",termination_status(model))
end

if abspath(PROGRAM_FILE)==@__FILE__
    strengthening_measurement(only(ARGS))
end
