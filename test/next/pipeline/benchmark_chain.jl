# Standalone build benchmark; run the Next correctness suite first.
# OPENBLAS_NUM_THREADS=1 julia --project test/next/pipeline/benchmark_chain.jl
# Contiguous degree-4 rows, periodic Pauli Heisenberg objective. Inputs and
# solves are outside timing. Each sample builds a fresh program; seven warm
# pairs alternate route order to limit drift. No cross-build cache is used.
using NCTSSoS, LinearAlgebra, Statistics
const Nx=NCTSSoS.Next
function chain_case(n,d=4)
    W=Nx.word_type(Nx.Pauli(),n,2d)
    words=Set([one(W)])
    for width in 1:d,start in 1:n,code in 0:3^width-1
        w,v=one(W),code
        for offset in 0:width-1
            w=Nx.mul(w,Nx.pauli_letter(W,mod1(start+offset,n),v%3))[2]; v÷=3
        end
        push!(words,w)
    end
    rows=sort!(collect(words))
    objective=[(1.,Nx.mul(Nx.pauli_letter(W,i,a),Nx.pauli_letter(W,mod1(i+1,n),a))[2]) for i in 1:n for a in 0:2]
    translation=Nx.SitePermutation([2:n;1])
    point=Nx.AbstractAction[Nx.SitePermutation(collect(n:-1:1)),
        Nx.CliffordMap([Nx.pauli_letter(W,i,a) for a in (1,0) for i in 1:n]),
        Nx.CliffordMap([Nx.pauli_letter(W,i,a) for a in (1,2) for i in 1:n])]
    anti=Nx.AntiunitaryPauli(one(W))
    rows,objective,translation,point,anti
end
function comparison(n,builds)
    for (name,f) in builds
        cold=@timed f()
        println("COLD ",name," n=",n," seconds=",cold.time," bytes=",cold.bytes)
        if cold.value isa Tuple && first(cold.value) isa Nx.MomentProgram
            p=first(cold.value)
            println("SIZE ",name," n=",n," vars=",p.nvars," blocks=",length(p.blocks)," maxside=",maximum(b.side for b in p.blocks))
        end
        flush(stdout)
    end
    samples=[NamedTuple{(:seconds,:bytes,:allocations),Tuple{Float64,Int,Int}}[] for _ in builds]
    for i in 1:7,j in (isodd(i) ? (1,2) : (2,1))
        GC.gc(); m=@timed builds[j][2]()
        push!(samples[j],(;seconds=m.time,bytes=m.bytes,allocations=Base.gc_alloc_count(m.gcstats)))
    end
    for ((name,_),s) in zip(builds,samples)
        println("SAMPLES ",name," n=",n," ",s)
        println("MEDIAN ",name," n=",n," seconds=",median(getproperty.(s,:seconds)),
            " IQR=",quantile(getproperty.(s,:seconds),[.25,.75])," bytes=",median(getproperty.(s,:bytes)))
    end
    flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
    println("ENV Julia=",VERSION," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads()," CPU=",Sys.cpu_info()[1].model)
    for n in (16,32)
        rows,h,t,point,anti=chain_case(n)
        println("ROWS n=",n," count=",length(rows)); flush(stdout)
        reg,ops=create_pauli_variables(1:n)
        pop=polyopt(heisenberg_chain_hamiltonian(ops;periodic=true,coupling=1),reg)
        comparison(n,(("lattice",()->Nx.lattice_program(rows,h,t,n,point;moment_generators=(anti,))),
            ("legacy",()->pauli_translation_invariant_moment_relaxation(pop,ops,4;axis_permutation_symmetry=true))))
    end
end
