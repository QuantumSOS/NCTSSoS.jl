# Compare tensor matrices by their Frobenius norm without materializing a large
# Kronecker product: <tensor(A),tensor(B)> = product(<A_i,B_i>).
mixed_raw(w)=map(Nx.legacy_letters,w.components)
function mixed_equal(oracles,left,right)
    terms=[left;[(-c,w) for (c,w) in right]]
    values=[map(operator_matrix,oracles,raw) for (_,raw) in terms]
    inner(a,b)=prod(dot(x,y) for (x,y) in zip(a,b))
    error=sum(conj(terms[i][1])*terms[j][1]*inner(values[i],values[j])
        for i in eachindex(terms),j in eachindex(terms))
    scale=sum(abs(c)*sqrt(real(inner(v,v))) for ((c,_),v) in zip(terms,values))
    return abs(error)<=1e-11*max(1,scale^2)
end

function mixed_oracles(algebras,rows)
    map(eachindex(algebras)) do i
        OperatorOracle(algebras[i],[Nx.legacy_letters(w.components[i]) for w in rows])
    end
end

function mixed_inspect!(out,words)
    for (i,w) in enumerate(words)
        h=zero(UInt)
        for s in Nx.sites(w); h=hash(s,h); end
        for l in Nx.word_letters(w); h=hash(l,h); end
        out[i]=h
    end
    out
end
function mixed_inspection_bytes(words)
    out=zeros(UInt,length(words))
    mixed_inspect!(out,words)
    @allocated mixed_inspect!(out,words)
end
function mixed_buffer_products!(out,words)
    for a in words,b in words
        empty!(out)
        Nx.mul!(out,a,b,1.0+1im)
    end
    out
end
function mixed_buffer_bytes(words::Vector{W}) where {W}
    out=Nx.PolynomialBuffer{W,ComplexF64}(128)
    mixed_buffer_products!(out,words)
    @allocated mixed_buffer_products!(out,words)
end

@testset "Packed tensor words" begin
    rng=MersenneTwister(460)
    cases=((Nx.Pauli(),Nx.Fermionic()),(Nx.Fermionic(),Nx.Bosonic()),
        (Nx.Pauli(),Nx.Pauli()),(Nx.Bosonic(),Nx.Bosonic()),
        (Nx.SiteAlgebra(Nx.Free(),2),Nx.SiteAlgebra(Nx.Unipotent(),2)),
        (Nx.SiteAlgebra(Nx.Projector(),2),Nx.Pauli()),
        (Nx.Pauli(),Nx.Fermionic(),Nx.Bosonic()))
    for algebras in cases
        names=ntuple(i->Symbol(:part,i),length(algebras))
        specs=NamedTuple{names}(map(a->(a,2),algebras))
        registry,_=Nx.mixed_variables(;specs...)
        W=Nx.word_type(registry,6)
        rows=Nx.basis_rows(W,registry.algebra,1:registry.nsites,2)
        oracles=mixed_oracles(algebras,rows)
        samples=unique!([one(W);rand(rng,rows,18)])
        check_word_contract(W,samples;triples=100)
        @test mixed_inspection_bytes(samples)==0
        @test mixed_buffer_bytes(samples)==0
        for w in samples
            folded=TermDict{W,ComplexF64}()
            Nx.push_term!(folded,1.0+0im,one(W))
            for letter in Nx.word_letters(w); folded=times(folded,letter); end
            @test clean(folded)==Dict(w=>1.0+0im)
            expected=map(raw_adjoint,algebras,mixed_raw(w))
            phase,s=Nx.star(w)
            @test mixed_equal(oracles,[(Nx.scale(phase,1.0+0im),mixed_raw(s))],[(1,expected)])
        end
        for _ in 1:60
            a,b=rand(rng,rows,2)
            got=[(c,mixed_raw(w)) for (w,c) in clean(product_terms(a,b,ComplexF64))]
            expected=map(vcat,mixed_raw(a),mixed_raw(b))
            @test mixed_equal(oracles,got,[(1,expected)])
        end
        # Enumerating each component then filtering is an independent basis oracle.
        parts=map((a,n,t)->Nx.basis_rows(t,a,1:n,2),algebras,ntuple(_->2,length(algebras)),fieldtypes(W.parameters[1]))
        @test Set(rows)==Set(W(ps) for ps in Iterators.product(parts...) if sum(Nx.degree,ps)<=2)
    end
    r,p=Nx.mixed_variables(spin=(Nx.Pauli(),65),matter=(Nx.Fermionic(),66),osc=(Nx.Bosonic(),17))
    x,_,_=p.spin; a,ad=p.matter; b,bd=p.osc
    samples=[last(only(Nx.terms(poly))) for poly in
        (one(x[1]),x[65],ad[66],b[17]^2,x[65]*a[66]*bd[17]^2)]
    check_word_contract(eltype(samples),samples;triples=40)
    @test mixed_inspection_bytes(samples)==0
    @test mixed_buffer_bytes(samples)==0
    @test collect(Nx.sites(last(samples)))==[65,131,148]
    W=Nx.word_type(Nx.mixed_variables(spin=(Nx.Pauli(),1),matter=(Nx.Fermionic(),1))[1],2)
    @test_throws r"outside" W((Nx.pauli_letter(Nx.PauliWord{1},2,0),one(Nx.FermiWord{1})))

    # Direct SingleTerm multiplication must retain the Pauli phase while the
    # unipotent factor squares to identity. Independent Kronecker matrices.
    _,p=Nx.mixed_variables(spin=(Nx.Pauli(),1),sign=(Nx.SiteAlgebra(Nx.Unipotent(),1),1))
    x,y,z=p.spin; u=p.sign[1][1]
    a,b=last(only(Nx.terms(x[1]*u))),last(only(Nx.terms(y[1]*u)))
    sx=[0 1;1 0]; sy=[0 -im;im 0]; sz=[1 0;0 -1]; id=[1 0;0 1]
    for (left,right,expected) in ((a,b,kron(sx,sz)*kron(sy,sz)),(b,a,kron(sy,sz)*kron(sx,sz)))
        phase,word=Nx.mul(left,right)
        @test word==last(only(Nx.terms(z[1])))
        @test Nx.scale(phase,1.0+0im)*kron(sz,id)==expected
    end
end
