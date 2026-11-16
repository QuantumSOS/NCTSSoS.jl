# Standalone, solver-free kernel measurements. Run the NCTSSoS correctness suite first.
# Input construction, compilation, and scratch allocation are outside measured regions.
using NCTSSoS, Random, Statistics, Test
const Nx = NCTSSoS

mutable struct DigestSink
    checksum::UInt
    count::Int
end
function Nx.push_term!(s::DigestSink,c::Number,w::Nx.AbstractWord)
    s.checksum ⊻= hash(w,hash(c)); s.count += 1
    return s
end
function packed_products!(sink,words)
    sink.checksum = 0; sink.count = 0
    for a in words, b in words
        Nx.mul!(sink,a,b,1.0+0im)
    end
    return sink.checksum
end
function packed_inspection(words)
    digest = UInt(0)
    for a in words, b in words
        digest ⊻= hash(Nx.star(a)[2],hash(a)) ⊻ UInt(isless(a,b))
        for site in Nx.sites(a); digest ⊻= UInt(site); end
    end
    return digest
end
function measure(label,f::F,args::Vararg{Any,N}) where {F,N}
    for _ in 1:3; f(args...); end
    times = Float64[]
    checksum = UInt(0)
    for _ in 1:31
        GC.gc(false)
        t = @elapsed for _ in 1:5
            checksum ⊻= hash(f(args...))
        end
        push!(times,t/5)
    end
    stats = @timed f(args...)
    println(label," median_us=",round(median(times)*1e6;digits=3),
        " q25_us=",round(quantile(times,.25)*1e6;digits=3),
        " q75_us=",round(quantile(times,.75)*1e6;digits=3),
        " bytes=",stats.bytes," allocations=",Base.gc_alloc_count(stats.gcstats)," checksum=",checksum)
end
function benchmark_words(label,words)
    println(label," words=",length(words)," pairs=",length(words)^2," sizeof_word=",sizeof(eltype(words)))
    measure("$label packed product",packed_products!,DigestSink(0,0),words)
    measure("$label packed inspection",packed_inspection,words)
end
function packed_polynomial!(out,a,b,scratch)
    Nx.mul!(out,a,b,scratch)
    return hash(out.n,hash(out.coefficients[1]))
end
function main()
    println("Julia=",VERSION," threads=",Threads.nthreads()," CPU=",Sys.CPU_NAME)
    rng = MersenneTwister(20261114)
    W = Nx.PauliWord{1}
    pauli = [Nx.PauliWord((rand(rng,UInt64)&0xfff,),(rand(rng,UInt64)&0xfff,)) for _ in 1:64]
    benchmark_words("Pauli",pauli)
    p,q,out,scratch = (Nx.PolynomialBuffer{W,ComplexF64}(4096) for _ in 1:4)
    for i in 1:32
        Nx.push_term!(p,ComplexF64(i%7-3),pauli[i])
        Nx.push_term!(q,ComplexF64(i%5-2),pauli[i+32])
    end
    Nx.canonicalize!(p); Nx.canonicalize!(q)
    measure("Pauli packed polynomial",packed_polynomial!,out,p,q,scratch)
    fermi = [Nx.word_from_legacy(Nx.FermiWord{1},sort!(unique(rand(rng,[-12:-1;1:12],3)))) for _ in 1:64]
    benchmark_words("Fermionic",fermi)
    for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent())
        W = Nx.word_type(Nx.SiteAlgebra(rule,4),4,8)
        words = [foldl((a,b)->Nx.mul(a,b)[2],(Nx.site_letter(W,rand(rng,1:4),rand(rng,1:4)) for _ in 1:4);init=one(W)) for _ in 1:64]
        benchmark_words(string(typeof(rule)),words)
    end
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
