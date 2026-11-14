# Standalone, build-only budget report. Uses the exact writer-benchmark inputs.
# JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=<SCS environment> test/next/benchmark_budget.jl
# Input polynomials are outside the timer; checked problem creation, row generation
# and moment-program assembly are inside. No solver is attached or called.
include("conic/benchmark_writers.jl")

function legacy_budget_build(h,registry,order)
    problem=polyopt(h,registry)
    sparsity=compute_sparsity(problem,SolverConfig(optimizer=SCS.Optimizer,order=order))
    return NCTSSoS.moment_relax(problem,sparsity.corr_sparsity,sparsity.cliques_term_sparsities)
end
native_budget_build(h,order)=N.build(N.Problem(h),N.Relaxation(order=order))

function budget_report()
    println("ENV Julia=",VERSION," CPU=",Sys.cpu_info()[1].model," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads())
    println("Seven alternating warm pairs; GC outside timer; no solves. Largest block is its native side, before any real embedding.")
    println("Real scalar slots include normalization; legacy count is from untimed JuMP lowering, native count from the moment program.")
    for name in ("heisenberg","hubbard","chsh","heisenberg10")
        old,registry,h,order,_=inputs(name)
        calls=(legacy=()->legacy_budget_build(old,registry,order),native=()->native_budget_build(h,order))
        reference=calls.legacy()
        built=calls.native()
        W=eltype(built.table.keys)
        function packed(m)
            letters=name=="chsh" ? [(mod1(NCTSSoS.decode_operator_id(l),2)<<8)|NCTSSoS.decode_site(l) for l in m.word] : m.word
            return N.word_from_legacy(W,letters)
        end
        @assert Set(N.moment_key(packed(m))[1] for m in values(reference.linear.key_to_monomial))==Set(label.word for label in built.table.labels)
        # Inspect the existing legacy lowering rather than estimating its real
        # variable count from complex moment keys. This is outside all timings.
        model,_=build_jump_model(reference)
        oldsize=(largest=maximum(b.size for b in reference.linear.psd_blocks_lin),
            moments=length(reference.linear.moments),real_slots=JuMP.num_variables(model),
            cones=unique(b.meta.cone for b in reference.linear.psd_blocks_lin))
        newsize=(largest=built.sizes.largest_block,moments=built.sizes.moments,
            real_slots=built.program.nvars,cones=unique(b.kind for b in built.program.blocks))
        samples=(legacy=NamedTuple[],native=NamedTuple[])
        for i in 1:7,route in (isodd(i) ? (:legacy,:native) : (:native,:legacy))
            GC.gc()
            t=@timed getproperty(calls,route)()
            push!(getproperty(samples,route),(seconds=t.time,bytes=t.bytes,allocations=Base.gc_alloc_count(t.gcstats)))
        end
        for route in (:legacy,:native)
            data=getproperty(samples,route)
            size=route==:legacy ? oldsize : newsize
            println("BUDGET model=",name," order=",order," route=",route,
                " allocations=",median(getproperty.(data,:allocations))," bytes=",median(getproperty.(data,:bytes)),
                " seconds=",median(getproperty.(data,:seconds))," IQR=",quantile(getproperty.(data,:seconds),[.25,.75])," sizes=",size)
            println("SAMPLES ",data)
        end
        println("RATIO model=",name," allocation_reduction=",median(getproperty.(samples.legacy,:allocations))/median(getproperty.(samples.native,:allocations)),
            " time_reduction=",median(getproperty.(samples.legacy,:seconds))/median(getproperty.(samples.native,:seconds)))
        flush(stdout)
    end
end
abspath(PROGRAM_FILE)==(@__FILE__) && budget_report()
