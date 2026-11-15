using Test, Random, NCTSSoS
isdefined(@__MODULE__,:check_word_contract) || include("../contract_checks.jl")
isdefined(@__MODULE__,:legacy_product_terms) || include("../reference.jl")

function boson_product_bytes(s,a,b,c)
    empty!(s); Nx.mul!(s,a,b,c)
    empty!(s)
    return @allocated Nx.mul!(s,a,b,c)
end

@testset "Packed bosonic words" begin
    rng=MersenneTwister(843)
    for n in (2,23)
        W=Nx.word_type(Nx.Bosonic(),n,18)
        words=[Nx.word_from_legacy(W,sort!(rand(rng,[-n,-1,1,n],rand(rng,0:3)))) for _ in 1:12]
        check_word_contract(W,unique([one(W);words]);triples=50)
        registry,_=create_bosonic_variables(1:n)
        T=eltype(indices(registry))
        for _ in 1:100
            a=Nx.word_from_legacy(W,sort!(rand(rng,[-n,-1,1,n],rand(rng,0:6))))
            b=Nx.word_from_legacy(W,sort!(rand(rng,[-n,-1,1,n],rand(rng,0:6))))
            new=Dict(Nx.legacy_letters(w)=>c for (w,c) in clean(product_terms(a,b,ComplexF64)))
            @test new==legacy_product_terms(BosonicAlgebra,T,Nx.legacy_letters(a),Nx.legacy_letters(b))
            letters=collect(Nx.word_letters(a))
            @test length(letters)==Nx.degree(a)
            @test all(l->Nx.degree(l)==1,letters)
            @test reduce(vcat,Nx.legacy_letters.(letters);init=Int[])==Nx.legacy_letters(a)
            _,adj=Nx.star(a)
            oldadj=NCTSSoS._moment_linear_adjoint_monomial(NormalMonomial{BosonicAlgebra,T}(T.(Nx.legacy_letters(a))))
            @test Nx.legacy_letters(adj)==Int.(oldadj.word)
            folded=TermDict{W,ComplexF64}(); Nx.push_term!(folded,1.0+0im,one(W))
            for l in letters; folded=times(folded,l); end
            @test clean(folded)==Dict(a=>1.0+0im)
            raw=shuffle(rng,Nx.legacy_letters(a))
            ordered=TermDict{W,ComplexF64}(); Nx.push_term!(ordered,1.0+0im,one(W))
            for i in raw; ordered=times(ordered,Nx.boson_letter(W,i)); end
            old=foldl(*,(NormalMonomial{BosonicAlgebra,T}(T[i]) for i in raw);
                init=1.0*NormalMonomial{BosonicAlgebra,T}(T[]))
            @test Dict(Nx.legacy_letters(w)=>c for (w,c) in clean(ordered))==
                Dict(Int.(m.word)=>c for (c,m) in terms(old))
        end
        a=Nx.word_from_legacy(W,[1,1,n])
        b=Nx.word_from_legacy(W,[-n,-1,-1])
        for C in (Int,Int128,Complex{Int},Float64,ComplexF64)
            sink=Nx.PolynomialBuffer{W,C}(6)
            @test boson_product_bytes(sink,a,b,one(C))==0
            @test length(sink)==6
        end
        permutation=Nx.SitePermutation(reverse(1:n))
        _,image=Nx.act(permutation,a)
        @test Nx.legacy_letters(image)==[1,n,n]
        @test Nx.act(permutation,image)==(one(Nx.Phase),a)
    end
    W=Nx.word_type(Nx.Bosonic(),1,255)
    high=Nx.word_from_legacy(W,fill(1,255))
    @test Nx.degree(high)==255
    @test clean(product_terms(high,one(W),Int8))==Dict(high=>Int8(1))
    @test_throws r"exponent 256.*capacity" Nx.mul!(CountSink(0),high,Nx.boson_letter(W,1),1.0)
    @test_throws r"degree at most 255" Nx.word_type(Nx.Bosonic(),1,256)
    @test_throws ArgumentError Nx.word_from_legacy(W,[1,-1])
    @test_throws ArgumentError Nx.boson_letter(W,0)
    @test_throws ArgumentError Nx.BosonWord{3,1}((UInt64(1)<<63,),(UInt64(0),))
    creators=Nx.word_from_legacy(W,fill(-1,21))
    annihilators=Nx.word_from_legacy(W,fill(1,21))
    @test_throws r"coefficient exceeds" Nx.mul!(CountSink(0),annihilators,creators,1)
    @test_throws r"coefficient exceeds" Nx.mul!(CountSink(0),annihilators,creators,floatmax(Float64))
    for C in (Int128,BigInt,Complex{Int128})
        terms=clean(product_terms(annihilators,creators,C))
        @test terms[one(W)]==factorial(big(21))
        @test terms[Nx.word_from_legacy(W,[-1,1])]==21factorial(big(21))
    end
    registry,(b,bd)=Nx.bosonic_variables(2)
    @test bd==adjoint.(b)
    @test b[1]*bd[1]-bd[1]*b[1]==one(b[1])
    @test b[1]*bd[2]==bd[2]*b[1]
    @test b[1]^2*bd[1]^2==bd[1]^2*b[1]^2+4bd[1]*b[1]+2
    @test sprint(show,bd[2]^2*b[1])=="bd2 bd2 b1"
    W=Nx.word_type(registry,6)
    w=Nx.word_from_legacy(W,[-2,-2,1])
    sink=TermDict{W,ComplexF64}()
    Nx.derive!(sink,Nx.NumberRotation(registry),w)
    @test clean(sink)==Dict(w=>1im)
    empty!(sink.terms)
    Nx.derive!(sink,Nx.NumberRotation(registry;modes=[2]),w)
    @test clean(sink)==Dict(w=>2im)
end
