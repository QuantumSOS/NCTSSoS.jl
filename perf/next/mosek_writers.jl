# Warm build + native-write comparisons; inputs and optimization are not timed.
# Run with the same environment as test/next/conic/mosek.jl (plus SCS for the
# shared writer-benchmark inputs). Julia allocations exclude Mosek's C heap.
include(joinpath(@__DIR__,"../../test/next/conic/benchmark_writers.jl"))
import Mosek, MosekTools

function write_case(p,r,route)
    f=N.moment_standard_form(N.build(p,r).program)
    options=(MSK_IPAR_NUM_THREADS=1,MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-10,
        MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-10,MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-10)
    if route==:native
        return N.prepare_standard(f,N.MosekBackend(;options...))
    end
    model,_=N.to_jump(f)
    JuMP.set_optimizer(model,JuMP.optimizer_with_attributes(MosekTools.Optimizer,
        (String(k)=>v for (k,v) in pairs(options))...))
    JuMP.MOI.Utilities.attach_optimizer(JuMP.backend(model))
    return model
end

function measure_writers(name;samples=5)
    _,_,h,order,expected=inputs(name)
    p,r=N.Problem(h),N.Relaxation(;order)
    for route in (:native,:jump)
        handle=write_case(p,r,route)
        if expected!==nothing
            bound=if route==:native
                result=N.solve_standard(handle)
                @assert result.status==JuMP.MOI.OPTIMAL
                result.bound
            else
                JuMP.optimize!(handle)
                @assert JuMP.termination_status(handle)==JuMP.MOI.OPTIMAL
                JuMP.objective_value(handle)
            end
            @assert isapprox(bound,expected;atol=1e-7)
        end
        route==:native && Mosek.deletetask(handle.task)
    end
    data=Dict(route=>NamedTuple[] for route in (:native,:jump))
    for i in 1:samples,route in (isodd(i) ? (:native,:jump) : (:jump,:native))
        GC.gc()
        t=@timed write_case(p,r,route)
        push!(data[route],(seconds=t.time,bytes=t.bytes,allocations=Base.gc_alloc_count(t.gcstats)))
        route==:native && Mosek.deletetask(t.value.task)
    end
    for route in (:native,:jump)
        s=getproperty.(data[route],:seconds)
        println("RESULT model=$name route=$route median_s=",median(s)," q25_s=",quantile(s,.25),
            " q75_s=",quantile(s,.75)," median_bytes=",median(getproperty.(data[route],:bytes)),
            " median_allocations=",median(getproperty.(data[route],:allocations)))
        println("SAMPLES ",data[route])
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads())
    for name in (isempty(ARGS) ? ["heisenberg","hubbard","chsh"] : ARGS)
        measure_writers(name)
    end
end
