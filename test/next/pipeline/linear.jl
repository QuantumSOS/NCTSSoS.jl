struct TestSpinDerivation <: Nx.AbstractDerivation
    n::Int
    axis::Int
end
function Nx.derive!(sink,d::TestSpinDerivation,w::ToyPauli,c::Number=1)
    for site in 1:d.n
        p,u = Nx.mul(toy_letter(site,d.axis),w)
        q,v = Nx.mul(w,toy_letter(site,d.axis))
        a = Nx.scale(p,complex(c)*im/2)
        b = -Nx.scale(q,complex(c)*im/2)
        # Each Pauli string is Hermitian; the commutator's imaginary
        # coefficient vanishes termwise after combining equal products.
        u == v || error("Pauli commutator words must coincide")
        iszero(a+b) || Nx.push_term!(sink,real(a+b),u)
    end
    sink
end
@testset "Linear generator equations" begin
    words=toy_words(2,2)
    index=Dict(w=>i for (i,w) in enumerate(words))
    A=Nx.derivation_constraints(words,[TestSpinDerivation(2,a) for a in 0:2]; C=Rational{Int})
    Z=zeros(Rational{Int},length(words),2)
    Z[index[one(ToyPauli)],1]=1
    for a in 0:2
        Z[index[Nx.mul(toy_letter(1,a),toy_letter(2,a))[2]],2]=1
    end
    @test iszero(A*Z)
    @test rank(Matrix{Float64}(A)) == 14
    @test_throws ArgumentError Nx.derivation_constraints([toy_letter(1,0)],[TestSpinDerivation(2,2)])
    W=Nx.word_type(Nx.Pauli(),2,2)
    packed=[Nx.word_from_legacy(W,Nx.legacy_letters(w)) for w in words]
    core=Nx.derivation_constraints(packed,[Nx.SpinRotation(a,(UInt64(3),)) for a in 0:2];C=Rational{Int})
    @test core == -A # the fixture uses +i/2, core uses the specified -i/2
end
