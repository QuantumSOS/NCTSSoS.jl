# Paired warmed build-plus-solve measurements, with GC and model input creation
# outside the timer. Run with no arguments for the three small models.
using NCTSSoS, JuMP, SCS, LinearAlgebra, Statistics
const N = NCTSSoS.Next
BLAS.set_num_threads(1)

function inputs(name)
    if name in ("heisenberg","heisenberg10")
        n = name=="heisenberg" ? 4 : 10
        r,(x,y,z) = create_pauli_variables(1:n)
        h = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:n-1)
        nr,(x,y,z) = N.pauli_variables(n)
        nh = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:n-1)
        return h,r,nh,n==4 ? 1 : 2,n==4 ? -9. : nothing
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

function benchmark(mode,name;samples=5,tol=1e-7)
    old,r,h,order,expected = inputs(name)
    mode in ("direct","indirect") || error("Unknown linear solver $mode")
    linear_solver = mode=="direct" ? SCS.DirectSolver : SCS.IndirectSolver
    options = (;verbose=0,eps_abs=tol,eps_rel=tol)
    backend = N.SCSBackend(;linear_solver,options...)
    factory = optimizer_with_attributes(SCS.Optimizer,"linear_solver"=>linear_solver,(String(k)=>v for (k,v) in pairs(options))...)
    calls = (legacy=()->cs_nctssos(polyopt(old,r),SolverConfig(optimizer=factory,order=order);dualize=false).objective,
        native=()->N.solve(N.Problem(h),N.Relaxation(order=order);solver=backend))
    println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
    reference = calls.legacy()
    isnothing(expected) || @assert isapprox(reference,expected;atol=2e-5)
    warm = calls.native()
    @assert warm.bound!==nothing && isapprox(warm.bound,reference;atol=5tol,rtol=5tol)
    println("ORACLE mode=$mode model=$name tol=$tol legacy=$reference native=",warm.bound," diagnostics=",warm.solution.diagnostics)
    measurements = (legacy=NamedTuple[],native=NamedTuple[])
    for i in 1:samples
        for route in (isodd(i) ? (:legacy,:native) : (:native,:legacy))
            GC.gc()
            t = @timed getproperty(calls,route)()
            bound = route==:legacy ? t.value : t.value.bound
            @assert bound!==nothing && isapprox(bound,reference;atol=5tol,rtol=5tol)
            push!(getproperty(measurements,route),(seconds=t.time,bytes=t.bytes,allocations=Base.gc_alloc_count(t.gcstats),bound=bound))
        end
    end
    for route in (:legacy,:native)
        data = getproperty(measurements,route)
        isempty(data) && continue
        seconds = [s.seconds for s in data]
        println("RESULT mode=$mode route=$route model=$name median_s=",median(seconds)," q25_s=",quantile(seconds,.25)," q75_s=",quantile(seconds,.75),
            " median_bytes=",median([s.bytes for s in data])," median_allocations=",median([s.allocations for s in data])," bound=",data[end].bound)
        println("SAMPLES ",data)
    end
    flush(stdout)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    for mode in ("direct","indirect"), name in ("heisenberg","hubbard","chsh")
        benchmark(mode,name)
    end
end
