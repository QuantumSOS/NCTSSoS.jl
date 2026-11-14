# Tests of the core layer of NCTSSoS.Next. Include the test files of this layer here.
using Test, Random, NCTSSoS
isdefined(@__MODULE__, :check_word_contract) || include("../contract_checks.jl")
isdefined(@__MODULE__, :legacy_product_terms) || include("../reference.jl")

function bridge_terms(s::TermDict)
    out = Dict{Vector{Int},ComplexF64}()
    for (w, c) in s.terms, (d, letters) in Nx.legacy_terms(w)
        out[letters] = get(out, letters, 0.0 + 0im) + c*d
    end
    return filter(kv -> !iszero(kv.second), out)
end

function inspect_into!(out, words)
    for (i,w) in enumerate(words)
        phase,s = Nx.star(w)
        key,p,conjugated = Nx.moment_key(w)
        h = hash(s,hash(w,hash(key))) ⊻ UInt(Nx.degree(w)) ⊻ UInt(phase.k) ⊻ UInt(p.k)
        h ⊻= UInt(isless(w,words[mod1(i+1,length(words))])) ⊻ UInt(conjugated)
        for site in Nx.sites(w)
            h = hash(site,h)
        end
        out[i] = h
    end
    return out
end
function inspection_allocations(words)
    out = zeros(UInt,length(words))
    inspect_into!(out,words)
    return @allocated inspect_into!(out,words)
end

function compare_words(A, T, words)
    # Write all results to caller storage so the compiler cannot discard inspection.
    @test inspection_allocations(words) == 0
    for a in words, b in words
        expected = Dict{Vector{Int},ComplexF64}()
        for (ca, la) in Nx.legacy_terms(a), (cb, lb) in Nx.legacy_terms(b)
            for (l,c) in legacy_product_terms(A,T,la,lb)
                expected[l] = get(expected,l,0.0+0im) + ca*cb*c
            end
        end
        @test bridge_terms(product_terms(a,b,ComplexF64)) == filter(kv->!iszero(kv.second),expected)
    end
    for a in words
        p,b = Nx.star(a)
        expected = Dict{Vector{Int},ComplexF64}()
        for (c,l) in Nx.legacy_terms(a)
            old = adjoint(1.0 * NormalMonomial{A,T}(T.(l)))
            for (d,m) in terms(old)
                key = Int.(m.word)
                expected[key] = get(expected,key,0.0+0im) + conj(c)*d
            end
        end
        got = Nx.mul!(TermDict{typeof(a),ComplexF64}(),one(typeof(a)),b,Nx.scale(p,1.0+0im))
        @test bridge_terms(got) == expected
        sink = TermDict{typeof(a),ComplexF64}()
        for (c,l) in Nx.legacy_terms(a)
            Nx.terms_from_legacy!(sink,typeof(a),l,c)
        end
        @test clean(sink) == Dict(a=>1.0+0im)
    end
end

@testset "Packed Pauli" begin
    rng = MersenneTwister(25)
    for K in (1,2,3)
        W = Nx.PauliWord{K}
        words = [one(W); [Nx.PauliWord(ntuple(_->rand(rng,UInt64),K),ntuple(_->rand(rng,UInt64),K)) for _ in 1:12]]
        check_word_contract(W,words)
        compare_words(PauliAlgebra,UInt16,words)
    end
    @test Nx.word_type(Nx.Pauli(),64,1000) === Nx.PauliWord{1}
    @test Nx.word_type(Nx.Pauli(),65,0) === Nx.PauliWord{2}
    @test_throws ArgumentError Nx.word_from_legacy(Nx.PauliWord{1},[1,2])
    @test_throws ArgumentError Nx.pauli_letter(Nx.PauliWord{1},65,0)
end

@testset "Packed site words" begin
    rng = MersenneTwister(71)
    for (rule,A) in ((Nx.Free(),NonCommutativeAlgebra),(Nx.Projector(),ProjectorAlgebra),(Nx.Unipotent(),UnipotentAlgebra))
        W = Nx.word_type(Nx.SiteAlgebra(rule,3),3,12)
        words = [one(W)]
        for _ in 1:20
            raw = [Nx.site_letter(W,rand(rng,1:3),rand(rng,1:3)) for _ in 1:rand(rng,1:4)]
            word = foldl((a,b)->Nx.mul(a,b)[2], raw; init=one(W))
            oldraw = UInt32[only(Nx.legacy_letters(w)) for w in raw]
            @test Nx.legacy_letters(word) == Int.(NCTSSoS.simplify(A,oldraw))
            push!(words,word)
        end
        check_word_contract(W,words; C=Float64)
        compare_words(A,UInt32,words)
        # Adjacent cancellation can expose another matching pair.
        a,b = Nx.site_letter(W,1,1),Nx.site_letter(W,1,2)
        ab,ba = Nx.mul(a,b)[2],Nx.mul(b,a)[2]
        rule isa Nx.Unipotent && @test isone(Nx.mul(ab,ba)[2])
        @test_throws ArgumentError Nx.word_from_legacy(W,[258,257])
        narrow = Nx.word_type(Nx.SiteAlgebra(rule,3),3,2)
        x,y = Nx.site_letter(narrow,1,1),Nx.site_letter(narrow,1,2)
        xy = Nx.mul(x,y)[2]
        @test_throws ArgumentError Nx.mul(xy,xy)
        wide = Nx.word_type(Nx.SiteAlgebra(rule,3),3,32)
        raw = [Nx.site_letter(wide,1,isodd(i) ? 1 : 2) for i in 1:19]
        w = foldl((a,b)->Nx.mul(a,b)[2],raw;init=one(wide))
        @test sizeof(wide) == 16
        @test Nx.degree(w) == 19
        @test Nx.star(Nx.star(w)[2])[2] == w
        @test Nx.word_from_legacy(wide,Nx.legacy_letters(w)) == w
        @test Nx.mul(w,one(wide))[2] == w
    end
    @test_throws ArgumentError Nx.word_type(Nx.SiteAlgebra(Nx.Free(),3),3,33)
    @test_throws ArgumentError Nx.SiteWord{Nx.Free,UInt8,3,1}(0xc0)
end

function buffer_dict(p)
    Dict(p.words[i] => p.coefficients[i] for i in 1:length(p))
end
function polynomial_allocations(out,a,b,scratch)
    Nx.add!(out,a,b)
    Nx.mul!(out,a,b,scratch)
    Nx.star!(out,a)
    return @allocated begin
        Nx.add!(out,a,b)
        Nx.mul!(out,a,b,scratch)
        Nx.star!(out,a)
    end
end
@testset "Caller-owned polynomial arithmetic" begin
    W = Nx.PauliWord{1}
    for C in (ComplexF64,Complex{Rational{Int}})
        a,b,out,scratch = (Nx.PolynomialBuffer{W,C}(32) for _ in 1:4)
        x,y,z = (Nx.pauli_letter(W,1,k) for k in 0:2)
        Nx.push_term!(a,C(2),x); Nx.push_term!(a,C(3),y); Nx.canonicalize!(a)
        Nx.push_term!(b,C(5),y); Nx.push_term!(b,C(-7),z); Nx.canonicalize!(b)
        Nx.add!(out,a,b,2,-1)
        @test buffer_dict(out) == Dict(x=>C(4),y=>C(1),z=>C(7))
        Nx.mul!(out,a,b,scratch)
        @test buffer_dict(out) == Dict(one(W)=>C(15),x=>C(-21im),y=>C(14im),z=>C(10im))
        Nx.star!(out,out)
        @test buffer_dict(out) == Dict(one(W)=>C(15),x=>C(21im),y=>C(-14im),z=>C(-10im))
        @test polynomial_allocations(out,a,b,scratch) == 0
        Nx.mul!(a,a,b,scratch)
        @test buffer_dict(a) == Dict(one(W)=>C(15),x=>C(-21im),y=>C(14im),z=>C(10im))
        @test_throws ArgumentError Nx.add!(a,a,b)
        @test_throws ArgumentError Nx.mul!(out,a,b,a)
        empty!(out)
        for (c,w) in ((2,z),(3,x),(-2,z),(-3,x),(5,y))
            Nx.push_term!(out,C(c),w)
        end
        Nx.canonicalize!(out)
        @test buffer_dict(out) == Dict(y=>C(5))
    end
end

@testset "Packed fermions" begin
    rng = MersenneTwister(123)
    for K in (1,2,3)
        W = Nx.FermiWord{K}
        modes = K == 1 ? [1,2,3,4] : [1,64,65,64K]
        words = [one(W)]
        for _ in 1:20
            creators = [-s for s in modes if rand(rng,Bool)]
            annihilators = [s for s in modes if rand(rng,Bool)]
            push!(words,Nx.word_from_legacy(W,sort!([creators;annihilators])))
        end
        check_word_contract(W,words; C=Float64)
        compare_words(FermionicAlgebra,Int,words)
        a,ad = Nx.fermi_letter(W,1),Nx.fermi_letter(W,-1)
        @test isempty(clean(product_terms(a,a,Float64)))
        @test clean(product_terms(a,ad,Float64)) == Dict(one(W)=>1.0,Nx.FermiWord(ad.c,a.a)=>-1.0)
        @test_throws ArgumentError Nx.word_from_legacy(W,[-1,-1])
        @test_throws ArgumentError Nx.word_from_legacy(W,[1,-1])
        @test_throws ArgumentError Nx.fermi_letter(W,0)
        # Raw products include contractions and nilpotency, not just input rejection.
        for _ in 1:100
            raw = rand(rng,[-modes; modes],rand(rng,0:7))
            new = TermDict{W,ComplexF64}(); Nx.push_term!(new,1.0+0im,one(W))
            old = 1.0 * one(NormalMonomial{FermionicAlgebra,Int})
            for l in raw
                new = times(new,Nx.fermi_letter(W,l))
                old *= NormalMonomial{FermionicAlgebra,Int}([l])
            end
            @test bridge_terms(new) == Dict(Int.(m.word)=>ComplexF64(c) for (c,m) in terms(old) if !iszero(c))
        end
        p,q,out,scratch = (Nx.PolynomialBuffer{W,Float64}(128) for _ in 1:4)
        Nx.push_term!(p,2.0,a); Nx.push_term!(p,3.0,ad); Nx.canonicalize!(p)
        Nx.push_term!(q,5.0,a); Nx.push_term!(q,7.0,ad); Nx.canonicalize!(q)
        Nx.mul!(out,p,q,scratch)
        @test buffer_dict(out) == Dict(one(W)=>14.0,Nx.FermiWord(ad.c,a.a)=>1.0)
        @test polynomial_allocations(out,p,q,scratch) == 0
    end
end

@testset "Registry guardrails" begin
    r = Nx.Registry(Nx.Pauli(),65,[:left_x=>(1,0),:right_y=>(65,1)])
    @test length(r) == 2
    @test r[:right_y] == (65,1)
    @test Nx.word_type(r,4) === Nx.PauliWord{2}
    @test Nx.word(r,:right_y,Nx.word_type(r,4)) == Nx.pauli_letter(Nx.PauliWord{2},65,1)
    @test_throws ArgumentError r[:absent]
    @test_throws ArgumentError Nx.Registry(Nx.Pauli(),1,[:x=>(1,0),:x=>(1,1)])
    @test_throws ArgumentError Nx.Registry(Nx.Pauli(),1,[:x=>(1,0),:y=>(1,0)])
    @test_throws ArgumentError Nx.Registry(Nx.Pauli(),1,[:x=>(2,0)])
    @test_throws ArgumentError Nx.Registry(Nx.Pauli(),1,[:x=>(1,3)])
    r = Nx.Registry(Nx.Fermionic(),2,[:create=>(2,-1),:destroy=>(2,1)])
    @test Nx.word(r,:create,Nx.word_type(r,6)) == Nx.fermi_letter(Nx.FermiWord{1},-2)
    r = Nx.Registry(Nx.SiteAlgebra(Nx.Projector(),4),3,[:p=>(3,4)])
    W = Nx.word_type(r,6)
    @test Nx.word(r,:p,W) == Nx.site_letter(W,3,4)
end

include("actions.jl")
include("moments.jl")
