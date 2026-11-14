# Full order-two Heisenberg rows, translation and reflection. No solver in the
# Next timed region. Three fresh warm builds; each phase owns its input and
# exact elimination may extend its table. First-call compilation is separate.
# OPENBLAS_NUM_THREADS=1 julia --project test/next/pipeline/benchmark_continuous.jl 8 16
using NCTSSoS, LinearAlgebra, Statistics, Profile
const Nx=NCTSSoS.Next

function continuous_case(n)
    registry,(x,y,z)=Nx.pauli_variables(n)
    H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    generators=(Nx.SitePermutation([2:n;1]),Nx.SitePermutation([1;n:-1:2]))
    ds=Nx.spin_rotations(registry)
    plain=Nx.Problem(H;symmetry=Nx.Symmetry(generators...))
    invariant=Nx.Problem(H;symmetry=Nx.Symmetry(generators...;continuous=ds))
    return plain,invariant,Nx.Relaxation(order=2),ds
end

function measure_continuous(n)
    plain,invariant,description,ds=continuous_case(n)
    cold=@timed Nx.build(invariant,description)
    report=cold.value.sizes.closure.continuous.blocks
    println("FIRST n=$n seconds=$(cold.time) bytes=$(cold.bytes) variables=$(cold.value.program.nvars-1) largest_before=$(report.largest_before) largest_after=$(report.largest_after) reconstruction=$(report.reconstruction_error) orthogonality=$(report.orthogonality_error) coefficient_bytes=$(report.max_coefficient_bytes)")
    flush(stdout)
    # The empty-continuous declaration has a different concrete type. Warm
    # that exact family call too, rather than timing its first specialization.
    warm=Nx.build(plain,description)
    p,_,_=Nx._continuous_program(warm.program,warm.table,ds)
    Nx._continuous_blocks(p)
    samples=NamedTuple[]
    for sample in 1:3
        GC.gc()
        family=@timed Nx.build(plain,description)
        input=family.value
        elimination=@timed Nx._continuous_program(input.program,input.table,ds)
        split=@timed Nx._continuous_blocks(first(elimination.value))
        @assert last(split.value).largest_after==report.largest_after
        @assert first(split.value).nvars==cold.value.program.nvars
        push!(samples,(family_seconds=family.time,elimination_seconds=elimination.time,split_seconds=split.time,
            family_bytes=family.bytes,elimination_bytes=elimination.bytes,split_bytes=split.bytes))
        println("SAMPLE n=$n i=$sample ",last(samples)); flush(stdout)
    end
    for key in keys(first(samples))
        values=getproperty.(samples,key)
        println("SUMMARY n=$n $key median=$(median(values)) range=$(extrema(values))")
    end
    if n==16
        input=Nx.build(plain,description)
        Profile.clear()
        Profile.@profile Nx._continuous_program(input.program,input.table,ds)
        println("ELIMINATION PROFILE n=$n")
        Profile.print(;format=:flat,sortedby=:count,mincount=100,C=true)
    end
    flush(stdout)
end

if abspath(PROGRAM_FILE)==@__FILE__
    println("ENV Julia=",VERSION," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads()," CPU=",Sys.cpu_info()[1].model)
    for n in (isempty(ARGS) ? (8,16) : parse.(Int,ARGS))
        measure_continuous(n)
    end
end
