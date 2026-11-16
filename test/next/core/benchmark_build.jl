# Standalone packed build measurements. No solver is invoked.
using NCTSSoS, Statistics
const Nx = NCTSSoS

function measure_build(name,f::F,args::Vararg{Any,N}) where {F,N}
    f(args...)
    times = Float64[]
    for _ in 1:7
        GC.gc()
        push!(times,@elapsed f(args...))
    end
    GC.gc()
    result = @timed f(args...)
    println(name," median_ms=",round(1000median(times);digits=3),
        " q25_ms=",round(1000quantile(times,.25);digits=3)," q75_ms=",round(1000quantile(times,.75);digits=3),
        " bytes=",result.bytes," allocations=",Base.gc_alloc_count(result.gcstats))
    flush(stdout)
end
function new_build(algebra,W,n,order,objective)
    rows = Nx.basis_rows(W,algebra,1:n,order)
    return Nx.plain_program(rows,objective)
end
function new_standard(algebra,W,n,order,objective)
    p,t = new_build(algebra,W,n,order,objective)
    return Nx.moment_standard_form(p),t
end
function benchmark_build(name,H,order)
    algebra=H.registry.algebra
    W=eltype(H.data.words)
    n=H.registry.nsites
    objective=collect(Nx.terms(H))
    p,t = new_build(algebra,W,n,order,objective)
    println(name," rows=",p.blocks[1].side," keys=",length(t.labels),
        " real_variables=",p.nvars," entries=",length(p.forms)," terms=",length(p.forms.ids))
    measure_build("$name packed rows+program",new_build,algebra,W,n,order,objective)
    measure_build("$name packed rows+program+standard",new_standard,algebra,W,n,order,objective)
end
function main()
    println("Julia=",VERSION," CPU=",Sys.CPU_NAME," threads=",Threads.nthreads())
    _,(x,y,z) = Nx.pauli_variables(10)
    ham = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:9)
    benchmark_build("Heisenberg10/d2",ham,2)
    _,(a,ad) = Nx.fermionic_variables(4)
    ham = sum(-(ad[i]*a[i+1]+ad[i+1]*a[i]) for i in 1:3)+1.3ad[1]*a[1]*ad[2]*a[2]
    benchmark_build("Fermion4/d2",ham,2)
    _,(x,y) = Nx.site_variables(Nx.Unipotent(),2,2)
    ham = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    benchmark_build("CHSH/d2",ham,2)
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
