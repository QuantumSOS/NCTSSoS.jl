# Warmed native build-plus-solve measurements, with GC and model input creation
# outside the timer. Run with no arguments for the three small models.
using NCTSSoS, JuMP, SCS, LinearAlgebra, Statistics
const N = NCTSSoS.Next
BLAS.set_num_threads(1)

function inputs(name)
    if name in ("heisenberg","heisenberg10")
        n = name=="heisenberg" ? 4 : 10
        nr,(x,y,z) = N.pauli_variables(n)
        nh = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:n-1)
        return nh,n==4 ? 1 : 2,n==4 ? -9. : nothing
    elseif name=="hubbard"
        nr,(a,ad) = N.fermionic_variables(4)
        nh = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
        return nh,2,1-sqrt(5)
    elseif name=="chsh"
        nr,(x,y) = N.site_variables(N.Unipotent(),2,2)
        nh = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
        return nh,1,-2sqrt(2)
    end
    error("Unknown model $name")
end

function benchmark(mode,name;samples=5,tol=1e-7)
    h,order,expected = inputs(name)
    mode in ("direct","indirect") || error("Unknown linear solver $mode")
    linear_solver = mode=="direct" ? SCS.DirectSolver : SCS.IndirectSolver
    options = (;verbose=0,eps_abs=tol,eps_rel=tol)
    backend = N.SCSBackend(;linear_solver,options...)
    call = ()->N.solve(N.Problem(h),N.Relaxation(order=order);solver=backend)
    println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
    warm = call()
    @assert warm.bound!==nothing
    isnothing(expected) || @assert isapprox(warm.bound,expected;atol=5tol,rtol=5tol)
    println("CHECK mode=$mode model=$name tol=$tol expected=$expected native=",warm.bound," diagnostics=",warm.solution.diagnostics)
    measurements = NamedTuple[]
    for i in 1:samples
        GC.gc()
        t = @timed call()
        bound = t.value.bound
        @assert bound!==nothing && isapprox(bound,warm.bound;atol=5tol,rtol=5tol)
        push!(measurements,(seconds=t.time,bytes=t.bytes,allocations=Base.gc_alloc_count(t.gcstats),bound=bound))
    end
    if !isempty(measurements)
        seconds = [s.seconds for s in measurements]
        println("RESULT mode=$mode route=native model=$name median_s=",median(seconds)," q25_s=",quantile(seconds,.25)," q75_s=",quantile(seconds,.75),
            " median_bytes=",median([s.bytes for s in measurements])," median_allocations=",median([s.allocations for s in measurements])," bound=",measurements[end].bound)
        println("SAMPLES ",measurements)
    end
    flush(stdout)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    for mode in ("direct","indirect"), name in ("heisenberg","hubbard","chsh")
        benchmark(mode,name)
    end
end
