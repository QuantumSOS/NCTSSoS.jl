# Standalone, solver-free M1 comparison. Run with julia --project test/next/core/benchmark.jl.
# Input construction, compilation, and scratch allocation are outside measured regions.
using NCTSSoS, Random, Statistics, Test
const Nx = NCTSSoS.Next

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
function old_products!(sink,words)
    sink.checksum = 0; sink.count = 0
    for a in words, b in words
        p = a*b
        for (c,w) in terms(p)
            sink.checksum ⊻= hash(w,hash(c)); sink.count += 1
        end
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
old_star(a) = adjoint(a)
old_star(a::NormalMonomial{FermionicAlgebra}) = adjoint(1.0 * a)
function old_inspection(words)
    digest = UInt(0)
    for a in words, b in words
        digest ⊻= hash(old_star(a),hash(a)) ⊻ UInt(isless(a,b))
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
function compare(label,A,T,words)
    old = [NormalMonomial{A,T}(T.(Nx.legacy_letters(w))) for w in words]
    sink = Nx.PolynomialBuffer{eltype(words),ComplexF64}(256)
    for (a,b) in zip(words,reverse(words))
        empty!(sink); Nx.mul!(sink,a,b,1.0+0im); Nx.canonicalize!(sink)
        got = Dict(Nx.legacy_letters(sink.words[i])=>sink.coefficients[i] for i in 1:length(sink))
        p = NormalMonomial{A,T}(T.(Nx.legacy_letters(a))) * NormalMonomial{A,T}(T.(Nx.legacy_letters(b)))
        @test got == Dict(Int.(m.word)=>ComplexF64(c) for (c,m) in terms(p) if !iszero(c))
    end
    println(label," words=",length(words)," pairs=",length(words)^2," sizeof_word=",sizeof(eltype(words)))
    measure("$label packed product",packed_products!,DigestSink(0,0),words)
    measure("$label legacy product",old_products!,DigestSink(0,0),old)
    # Packed inspection also visits sites; legacy inspection does not (conservative comparison).
    measure("$label packed inspection",packed_inspection,words)
    measure("$label legacy inspection",old_inspection,old)
    return old
end
function packed_polynomial!(out,a,b,scratch)
    Nx.mul!(out,a,b,scratch)
    return hash(out.n,hash(out.coefficients[1]))
end
old_polynomial(a,b) = hash(a*b)

function main()
    println("Julia=",VERSION," threads=",Threads.nthreads()," CPU=",Sys.CPU_NAME)
    rng = MersenneTwister(20261114)
    W = Nx.PauliWord{1}
    pauli = [Nx.PauliWord((rand(rng,UInt64)&0xfff,),(rand(rng,UInt64)&0xfff,)) for _ in 1:64]
    old = compare("Pauli",PauliAlgebra,UInt16,pauli)
    p,q,out,scratch = (Nx.PolynomialBuffer{W,ComplexF64}(4096) for _ in 1:4)
    for i in 1:32
        Nx.push_term!(p,ComplexF64(i%7-3),pauli[i])
        Nx.push_term!(q,ComplexF64(i%5-2),pauli[i+32])
    end
    Nx.canonicalize!(p); Nx.canonicalize!(q)
    oldp = sum((i%7-3)*old[i] for i in 1:32)
    oldq = sum((i%5-2)*old[i+32] for i in 1:32)
    Nx.mul!(out,p,q,scratch)
    @test Dict(Nx.legacy_letters(out.words[i])=>out.coefficients[i] for i in 1:length(out)) == Dict(Int.(m.word)=>ComplexF64(c) for (c,m) in terms(oldp*oldq) if !iszero(c))
    measure("Pauli packed polynomial",packed_polynomial!,out,p,q,scratch)
    measure("Pauli legacy polynomial",old_polynomial,oldp,oldq)
    fermi = [Nx.word_from_legacy(Nx.FermiWord{1},sort!(unique(rand(rng,[-12:-1;1:12],3)))) for _ in 1:64]
    compare("Fermionic",FermionicAlgebra,Int,fermi)
    for (rule,A) in ((Nx.Free(),NonCommutativeAlgebra),(Nx.Projector(),ProjectorAlgebra),(Nx.Unipotent(),UnipotentAlgebra))
        W = Nx.word_type(Nx.SiteAlgebra(rule,4),4,8)
        words = [foldl((a,b)->Nx.mul(a,b)[2],(Nx.site_letter(W,rand(rng,1:4),rand(rng,1:4)) for _ in 1:4);init=one(W)) for _ in 1:64]
        compare(string(typeof(rule)),A,UInt32,words)
    end
end
main()
