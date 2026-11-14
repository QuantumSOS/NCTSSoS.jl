# Fresh-process first build after package precompilation; no solver or loading
# inside the timer. Optional first argument saves/checks exact program data.
using NCTSSoS, LinearAlgebra, Serialization, Statistics
const Nx=NCTSSoS.Next

function block_key_measurement(path=nothing)
    n=12
    _,(x,y,z)=Nx.pauli_variables(n)
    H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.5z[i]*z[mod1(i+1,n)] for i in 1:n)
    problem=Nx.Problem(H)
    translation=Nx.SitePermutation([2:n;1])
    reflection=Nx.SitePermutation([mod1(2-i,n) for i in 1:n])
    description=Nx.BlockRelaxation(order=2,grouping=Nx.Correlative(),split=Nx.TermSplit(),
        symmetry=Nx.FamilySymmetry([translation,reflection]))
    println("ENV Julia=",VERSION," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads()," CPU=",Sys.cpu_info()[1].model)
    cold=@timed Base.invokelatest(Nx.build,problem,description)
    println("COLD seconds=",cold.time," compile=",cold.compile_time," recompile=",cold.recompile_time," bytes=",cold.bytes)
    built=cold.value; p=built.program; closure=built.sizes.closure
    data=(nvars=p.nvars,blocks=p.blocks,objective=(p.objective_ids,p.objective_coefs),
        forms=(p.forms.offsets,p.forms.ids,p.forms.coefs),labels=built.table.labels,
        family=[(b.source,b.rows) for b in closure.family],orbits=closure.block_orbits)
    if !isnothing(path)
        if isfile(path)
            @assert data==deserialize(path) "Program, labels or closure changed"
            println("ORACLE exact match")
        else
            serialize(path,data)
            println("ORACLE saved")
        end
    end
    println("SHAPE ",built.sizes," row_counts=",sort!(unique(length(b.rows) for b in closure.family)),
        " key_types=",unique(typeof(Nx._block_key(b)) for b in closure.family))
    warm=map(1:7) do _
        GC.gc()
        m=@timed Base.invokelatest(Nx.build,problem,description)
        (m.time,m.bytes)
    end
    println("WARM median=",median(first.(warm))," IQR=",quantile(first.(warm),[0.25,0.75])," bytes=",median(last.(warm))," samples=",warm)
end

if abspath(PROGRAM_FILE)==@__FILE__
    block_key_measurement(isempty(ARGS) ? nothing : only(ARGS))
end
