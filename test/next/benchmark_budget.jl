# Standalone, build-only budget report. Uses the writer-benchmark inputs.
# Input construction and GC are outside the timer; no solver is called.
include("conic/benchmark_writers.jl")

native_budget_build(h,order)=N.build(N.Problem(h),N.Relaxation(order=order))

function budget_report()
    println("ENV Julia=",VERSION," CPU=",Sys.cpu_info()[1].model," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads())
    println("Seven warm samples; GC outside timer; no solves. Real scalar slots include normalization.")
    for name in ("heisenberg","hubbard","chsh","heisenberg10")
        h,order,_=inputs(name)
        built=native_budget_build(h,order)
        sizes=(largest=built.sizes.largest_block,moments=built.sizes.moments,
            real_slots=built.program.nvars,cones=unique(b.kind for b in built.program.blocks))
        samples=NamedTuple[]
        for _ in 1:7
            GC.gc()
            t=@timed native_budget_build(h,order)
            push!(samples,(seconds=t.time,bytes=t.bytes,allocations=Base.gc_alloc_count(t.gcstats)))
        end
        println("BUDGET model=",name," order=",order," route=native",
            " allocations=",median(getproperty.(samples,:allocations))," bytes=",median(getproperty.(samples,:bytes)),
            " seconds=",median(getproperty.(samples,:seconds))," IQR=",quantile(getproperty.(samples,:seconds),[.25,.75])," sizes=",sizes)
        println("SAMPLES ",samples)
        flush(stdout)
    end
end
abspath(PROGRAM_FILE)==(@__FILE__) && budget_report()
