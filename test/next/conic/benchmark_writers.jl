# Run one solver/route/model per process. /usr/bin/time -v records whole-process
# peak RSS (including load/JIT); @timed records warmed build-plus-solve only.
using NCTSSoS, JuMP, SCS, Mosek, MosekTools, BPSDP, LinearAlgebra, Statistics
const N = NCTSSoS.Next
BLAS.set_num_threads(1)

function inputs(name)
    if name=="heisenberg"
        r,(x,y,z) = create_pauli_variables(1:4)
        h = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
        nr,(x,y,z) = N.pauli_variables(4)
        nh = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
        return h,r,nh,1,-9.
    elseif name=="hubbard"
        r,(a,ad) = create_fermionic_variables(1:4)
        h = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
        nr,(a,ad) = N.fermionic_variables(4)
        nh = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
        return h,r,nh,2,1-sqrt(5)
    elseif name=="chsh"
        r,(x,y) = create_unipotent_variables([("x",1:2),("y",1:2)])
        h = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
        nr,(x,y) = N.site_variables(N.Unipotent(),2,2)
        nh = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
        return h,r,nh,1,-2sqrt(2)
    end
    error("Unknown model $name")
end

function benchmark(solver,route,name)
    old,r,h,order,expected = inputs(name)
    if solver=="scs"
        options = (;verbose=0,eps_abs=1e-7,eps_rel=1e-7)
        direct = N.SCSBackend(;options...)
        factory = optimizer_with_attributes(SCS.Optimizer,(String(k)=>v for (k,v) in pairs(options))...)
    elseif solver=="mosek"
        direct = N.MosekBackend(MSK_IPAR_NUM_THREADS=1)
        factory = optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_NUM_THREADS"=>1)
    elseif solver=="bpsdp"
        options = (;print_level=0,sdp_objective_convergence_abs=1e-7,sdp_objective_convergence_rel=1e-7,
            sdp_error_convergence_abs=1e-7,sdp_error_convergence_rel=1e-7)
        direct = N.BPSDPBackend(;options...)
        factory = optimizer_with_attributes(BPSDP.Optimizer,(String(k)=>v for (k,v) in pairs(options))...)
    else
        error("Unknown solver $solver")
    end
    call = route=="legacy" ? ()->cs_nctssos(polyopt(old,r),SolverConfig(optimizer=factory,order=order);dualize=solver=="bpsdp").objective :
        ()->N.solve(N.Problem(h),N.Relaxation(order=order);solver=direct).bound
    println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
    first = @timed call()
    @assert first.value!==nothing && isapprox(first.value,expected;atol=2e-5) (first.value,expected)
    samples = NamedTuple[]
    for _ in 1:5
        GC.gc()
        t = @timed call()
        @assert t.value!==nothing && isapprox(t.value,expected;atol=2e-5) (t.value,expected)
        push!(samples,(seconds=t.time,bytes=t.bytes,allocations=Base.gc_alloc_count(t.gcstats),bound=t.value))
    end
    seconds = [s.seconds for s in samples]
    println("RESULT solver=$solver route=$route model=$name cold_s=",first.time,
        " median_s=",median(seconds)," q25_s=",quantile(seconds,.25)," q75_s=",quantile(seconds,.75),
        " median_bytes=",median([s.bytes for s in samples])," median_allocations=",median([s.allocations for s in samples]),
        " bound=",samples[end].bound," process_peak_rss=",Sys.maxrss())
    println("SAMPLES ",samples)
end
benchmark(ARGS...)
