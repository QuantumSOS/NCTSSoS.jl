# julia --project perf/next/functional.jl cost
# julia --project=<solver environment> perf/next/functional.jl first-scs
# Run first-scs in three fresh processes after precompilation.
using NCTSSoS, LinearAlgebra, Test
const N=NCTSSoS.Next
const SIGNS=[1 1 -1;1 1 1;-1 1 0]
cov(A,B,s)=-sum(SIGNS[i,j]*(s(A[i]*B[j])-s(A[i])*s(B[j])) for i in 1:3,j in 1:3 if SIGNS[i,j]!=0)
function model(trace)
    if trace
        _,(x,)=N.site_variables(N.Unipotent(),1,6;functional=N.RealTrace())
        return N.Problem(cov(x[1:3],x[4:6],N.trace_symbol))
    end
    _,(a,b)=N.site_variables(N.Unipotent(),2,3;functional=N.RealState())
    N.Problem(cov(a,b,N.state_symbol))
end
function samples(label,f)
    value=f()
    for _ in 1:3
        GC.gc()
        t=@timed f()
        println("SAMPLE ",label," seconds=",t.time," bytes=",t.bytes," allocations=",Base.gc_alloc_count(t.gcstats))
    end
    value
end
function product_kernel!(sink,words,repeats)
    checksum=UInt(0)
    for _ in 1:repeats,i in eachindex(words)
        empty!(sink); N.mul!(sink,words[i],words[mod1(i+1,length(words))],1.0)
        checksum⊻=hash(sink.words[1])
    end
    checksum
end
function closure_kernel!(out,words,repeats)
    for _ in 1:repeats,i in eachindex(words)
        out[i]=first(N.moment_key(words[i]))
    end
    out
end
function kernel_words(F,limbs)
    algebra=N.FunctionalAlgebra(F(),N.SiteAlgebra(N.Free(),3))
    W=N.word_type(algebra,1,limbs==1 ? 1 : N._functional_capacity(3,limbs-1)+1)
    words=W[]
    for offset in 1:8
        w=one(W)
        for j in 1:N._fcap(W)÷3
            w=N.mul(w,N.site_letter(W,1,mod1(j+offset,3)))[2]
        end
        scalar=first(N.moment_key(N.site_letter(W,1,mod1(offset,3))))
        push!(words,N.mul(w,scalar)[2])
    end
    words
end
function cost()
    for trace in (false,true)
        label=trace ? "trace" : "state"
        p=samples(label*"_new_model",()->model(trace))
        built=samples(label*"_new_assembly",()->N.build(p,N.Relaxation(order=2)))
        form=samples(label*"_new_standard",()->N.moment_standard_form(built.program))
        samples(label*"_new_sos_standard",()->N.sos_standard_form(built.program))
        println("DIMENSIONS ",label," ",built.sizes," standard=",size(form.A))
    end
    for F in (N.RealState,N.RealTrace),k in (1,2,4,8)
        words=kernel_words(F,k); W=eltype(words)
        sink=N.PolynomialBuffer{W,Float64}(1); out=similar(words)
        product_kernel!(sink,words,1); closure_kernel!(out,words,1)
        @test (@allocated product_kernel!(sink,words,100))==0
        @test (@allocated closure_kernel!(out,words,100))==0
        label=string(F,"_limbs",k)
        println("KERNEL ",label," word_bytes=",sizeof(W)," input_degree=",N.degree(first(words))," calls=8000")
        samples(label*"_product",()->product_kernel!(sink,words,1000))
        samples(label*"_closure",()->closure_kernel!(out,words,1000))
    end
end
function first_calls(backend)
    for trace in (false,true)
        # Keep inference/code generation inside each timed call, rather than
        # allowing compilation of this driver to compile its callees first.
        input=@timed Base.invokelatest(model,trace)
        p=input.value
        println("INPUT trace=",trace," seconds=",input.time)
        for i in 1:2
            t=@timed Base.invokelatest(N.solve,p,N.Relaxation(order=2);solver=backend)
            bound=t.value.bound
            @test bound ≈ -5 atol=1e-3
            println("SOLVE trace=",trace," call=",i," seconds=",t.time," bytes=",t.bytes," bound=",bound)
        end
    end
end

BLAS.set_num_threads(1)
println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
mode=isempty(ARGS) ? "cost" : only(ARGS)
if mode=="cost"
    cost()
elseif mode=="first-scs"
    using SCS
    first_calls(N.SCSBackend(verbose=0,eps_abs=1e-6,eps_rel=1e-6))
elseif mode=="first-mosek"
    using Mosek
    first_calls(N.MosekBackend())
else
    error("Use cost, first-scs or first-mosek")
end
