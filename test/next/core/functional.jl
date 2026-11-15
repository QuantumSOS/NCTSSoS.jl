module FunctionalCoreChecks
using Test, Random, LinearAlgebra, NCTSSoS
import NCTSSoS: Next as N
include("../contract_checks.jl")

legacy_rule(::Type{N.Free})=NonCommutativeAlgebra
legacy_rule(::Type{N.Projector})=ProjectorAlgebra
legacy_rule(::Type{N.Unipotent})=UnipotentAlgebra
function legacy(w::N.FunctionalWord{F,A,B,S}) where {F,A,B,S}
    LA=legacy_rule(A); ST=F===N.RealState ? Arbitrary : MaxEntangled
    letters(i,j)=UInt32[(UInt32(N._ftoken(w,k)>>S)<<8)|UInt32(N._fsite(N._ftoken(w,k),S)) for k in i:j]
    tail=NormalMonomial{LA,UInt32}(letters(1,N._ftail(w)))
    syms=NCTSSoS.StateSymbol{ST,LA,UInt32}[]
    i=N._ftail(w)+1
    while i<=N.degree(w)
        j=N._fstop(w,i)
        push!(syms,NCTSSoS.StateSymbol{ST}(NormalMonomial{LA,UInt32}(letters(i,j))))
        i=j+1
    end
    NCTSSoS.NCStateWord(StateWord{ST,LA,UInt32}(syms),tail)
end

function sample(rng,::Type{W},nsites) where {W}
    w=one(W)
    for _ in 1:rand(rng,0:3)
        v=one(W)
        for _ in 1:rand(rng,1:3)
            v=N.mul(v,N.site_letter(W,rand(rng,1:nsites),rand(rng,1:3)))[2]
        end
        rand(rng,Bool) && (v=first(N.moment_key(v)))
        w=N.mul(w,v)[2]
    end
    w
end
function exercise!(out,a,b)
    for i in eachindex(out)
        p,w=N.mul(a,b)
        w=first(N.moment_key(w))
        _,w=N.star(w)
        h=hash(w)
        for l in N.support_letters(w); h=hash(l,h); end
        for s in N.sites(w); h=hash(s,h); end
        out[i]=h
    end
end
function allocations(a,b)
    out=zeros(UInt,4)
    exercise!(out,a,b)
    @allocated exercise!(out,a,b)
end

# Independent orbit enumeration: rotate each site substring independently,
# and reverse the entire word. No production canonicalizer is used.
function trace_orbit(parts)
    candidates=Tuple[]
    for blocks in (parts,reverse.(parts))
        rotations=[[circshift(p,k) for k in 0:length(p)-1] for p in blocks]
        for choice in Iterators.product(rotations...)
            push!(candidates,Tuple(vcat(choice...)))
        end
    end
    minimum(candidates)
end
function pureword(::Type{W},parts) where {W}
    w=one(W)
    for (site,labels) in enumerate(parts),label in labels
        w=N.mul(w,N.site_letter(W,site,label))[2]
    end
    w
end

@testset "Packed real functional words" begin
    rng=MersenneTwister(8263)
    for F in (N.RealState,N.RealTrace), A in (N.Free,N.Projector,N.Unipotent)
        algebra=N.FunctionalAlgebra(F(),N.SiteAlgebra(A(),3))
        # Trace legacy comparison has one group: its multi-group oracle is wrong.
        nsites=F===N.RealTrace ? 1 : 2
        W=N.word_type(algebra,nsites,30)
        samples=[sample(rng,W,nsites) for _ in 1:50]
        check_word_contract(W,samples[1:8])
        for w in samples
            @test W(w.bits)==w
            @test isbitstype(W)
            @test legacy(N.star(w)[2])==adjoint(legacy(w))
            @test legacy(first(N.moment_key(w))).sw==NCTSSoS.expval(legacy(w))
            support=collect(N.support_letters(w))
            @test length(support)==N.degree(w)
            @test all(l->N.degree(l)==1,support)
            @test sort!(unique!(reduce(vcat,collect.(N.sites.(support));init=Int[])))==collect(N.sites(w))
            if N._ftail(w)==N.degree(w)
                @test collect(N.word_letters(w))==support
            else
                @test_throws r"support_letters" N.word_letters(w)
            end
            @test N.moment_key(first(N.moment_key(w)))==(first(N.moment_key(w)),one(N.Phase),false)
        end
        for _ in 1:150
            a,b,c=rand(rng,samples,3)
            @test legacy(N.mul(a,b)[2])==legacy(a)*legacy(b)
            @test N.mul(N.mul(a,b)[2],c)[2]==N.mul(a,N.mul(b,c)[2])[2]
            @test N.star(N.mul(a,b)[2])[2]==N.mul(N.star(b)[2],N.star(a)[2])[2]
        end
        @test allocations(samples[3],samples[5])==0
    end

    @testset "Independent two-group trace oracle" begin
        algebra=N.FunctionalAlgebra(N.RealTrace(),N.SiteAlgebra(N.Free(),3))
        W=N.word_type(algebra,2,20)
        matrices=[let a=randn(rng,3,3); a+a' end for _ in 1:6]
        id=Matrix{Float64}(I,3,3)
        ops=[[kron(matrices[i],id) for i in 1:3],[kron(id,matrices[i]) for i in 4:6]]
        for _ in 1:100
            parts=[rand(rng,1:3,rand(rng,1:4)) for _ in 1:2]
            key=first(N.moment_key(pureword(W,parts)))
            canonical=[[Int(N._ftoken(key,i)>>2) for i in 1:N.degree(key) if N._fsite(N._ftoken(key,i),2)==s] for s in 1:2]
            # Oracle tokens include the site so ordering agrees without losing boundaries.
            tagged=[[(label<<2)|s for label in parts[s]] for s in 1:2]
            @test Tuple(Int(N._ftoken(key,i)) for i in 1:N.degree(key))==trace_orbit(tagged)
            evaluate(ps)=prod(prod(ops[s][l] for l in ps[s]) for s in 1:2)
            @test LinearAlgebra.tr(evaluate(parts))/9 ≈ LinearAlgebra.tr(evaluate(canonical))/9 atol=1e-10
            @test W(key.bits)==key
        end
        # Evaluate mixed words independently as scalar factors times matrices.
        # The vector state is not tracial; both domains share the same operators.
        ψ=randn(rng,9); ψ/=norm(ψ)
        matrix_segment(w,indices)=foldl((m,i)->m*ops[N._fsite(N._ftoken(w,i),2)][Int(N._ftoken(w,i)>>2)],
            indices;init=Matrix{Float64}(I,9,9))
        function matrix_word(w,phi)
            value=matrix_segment(w,1:N._ftail(w)); i=N._ftail(w)+1
            while i<=N.degree(w)
                j=N._fstop(w,i)
                value*=phi(matrix_segment(w,i:j)); i=j+1
            end
            value
        end
        for F in (N.RealState,N.RealTrace)
            W=N.word_type(N.FunctionalAlgebra(F(),N.SiteAlgebra(N.Free(),3)),2,30)
            phi=F===N.RealState ? m->dot(ψ,m*ψ) : m->LinearAlgebra.tr(m)/9
            for _ in 1:50
                a,b=sample(rng,W,2),sample(rng,W,2)
                @test matrix_word(N.mul(a,b)[2],phi) ≈ matrix_word(a,phi)*matrix_word(b,phi)
                @test matrix_word(N.star(a)[2],phi) ≈ matrix_word(a,phi)'
                @test matrix_word(first(N.moment_key(a)),phi) ≈ phi(matrix_word(a,phi))*Matrix{Float64}(I,9,9)
            end
        end
    end

    @testset "Widths and capacity" begin
        for F in (N.RealState,N.RealTrace), limbs in 1:8
            algebra=N.FunctionalAlgebra(F(),N.SiteAlgebra(N.Free(),3))
            required=limbs==1 ? 1 : N._functional_capacity(3,limbs-1)+1
            W=N.word_type(algebra,1,required)
            @test sizeof(W)==8limbs
            a=pureword(W,[[mod1(i,3) for i in 1:required÷2]])
            b=first(N.moment_key(pureword(W,[[3,1,2]])))
            @test W(N.mul(a,b)[2].bits)==N.mul(a,b)[2]
            @test allocations(a,b)==0
            maxword=pureword(W,[fill(1,N._fcap(W))])
            @test N.degree(maxword)==N._fcap(W)
            @test_throws r"requires degree.*token.*capacity" N.mul(maxword,N.site_letter(W,1,1))
            @test_throws r"requires degree.*token.*capacity" N.word_type(algebra,1,N._functional_capacity(3,8)+1)
        end
        algebra=N.FunctionalAlgebra(N.RealState(),N.SiteAlgebra(N.Unipotent(),3))
        W=N.word_type(algebra,1,1)
        a=pureword(W,[[mod1(i,3) for i in 1:N._fcap(W)]])
        @test isone(N.mul(a,N.star(a)[2])[2])
        @test_throws ArgumentError W((typemax(UInt64),))
    end

    @testset "Caller-owned polynomial storage" begin
        algebra=N.FunctionalAlgebra(N.RealState(),N.SiteAlgebra(N.Free(),3))
        W=N.word_type(algebra,2,8)
        a=N.site_letter(W,1,1); b=N.site_letter(W,2,2)
        sa=first(N.moment_key(a)); sb=first(N.moment_key(b))
        p,q,out,scratch=(N.PolynomialBuffer{W,Float64}(4) for _ in 1:4)
        N.push_term!(p,2.0,a); N.push_term!(p,-3.0,sb); N.canonicalize!(p)
        N.push_term!(q,5.0,sa); N.push_term!(q,7.0,b); N.canonicalize!(q)
        N.mul!(out,p,q,scratch)
        expected=Dict(legacy(a)*legacy(sa)=>10.0,legacy(a)*legacy(b)=>14.0,
                      legacy(sb)*legacy(sa)=>-15.0,legacy(sb)*legacy(b)=>-21.0)
        @test Dict(legacy(out.words[i])=>out.coefficients[i] for i in 1:length(out))==expected
        @test (@allocated N.mul!(out,p,q,scratch))==0
    end
end
end
