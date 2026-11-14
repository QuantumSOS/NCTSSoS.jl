using LinearAlgebra

function permutation_matrix(images)
    n = length(images)
    U = zeros(ComplexF64,2^n,2^n)
    for col in 0:2^n-1
        row = sum(((col>>(i-1))&1)<<(images[i]-1) for i in 1:n)
        U[row+1,col+1] = 1
    end
    return U
end

function pauli_matrix(w,n)
    X = ComplexF64[0 1;1 0]; Y = ComplexF64[0 -im;im 0]; Z = ComplexF64[1 0;0 -1]
    letters = Dict((l-1)÷3+1 => (l-1)%3 for l in Nx.legacy_letters(w))
    out = ones(ComplexF64,1,1)
    for site in n:-1:1
        local_matrix = !haskey(letters,site) ? Matrix{ComplexF64}(I,2,2) : (X,Y,Z)[letters[site]+1]
        out = kron(out,local_matrix)
    end
    return out
end
function gate_matrix(gate,n,a,b=0)
    U = zeros(ComplexF64,2^n,2^n)
    for col in 0:2^n-1
        if gate==:H
            bit = (col>>(a-1))&1
            U[col+1,col+1] = (bit==0 ? 1 : -1)/sqrt(2)
            U[(col ⊻ (1<<(a-1)))+1,col+1] = 1/sqrt(2)
        elseif gate==:S
            U[col+1,col+1] = ((col>>(a-1))&1)==0 ? 1 : im
        elseif gate==:CNOT
            row = ((col>>(a-1))&1)==0 ? col : col ⊻ (1<<(b-1))
            U[row+1,col+1] = 1
        else
            different = ((col>>(a-1)) ⊻ (col>>(b-1)))&1
            row = different==0 ? col : col ⊻ (1<<(a-1)) ⊻ (1<<(b-1))
            U[row+1,col+1] = 1
        end
    end
    return U
end
function action_allocations(g,w)
    sink = CountSink(0)
    Nx.act!(sink,g,w,1.0+0im)
    return @allocated Nx.act!(sink,g,w,1.0+0im)
end
function derivation_allocations(d,w)
    sink = CountSink(0)
    Nx.derive!(sink,d,w,1.0)
    return @allocated Nx.derive!(sink,d,w,1.0)
end

@testset "Canonical mixed Pauli actions" begin
    W = Nx.PauliWord{1}
    for n in 1:4
        rows = Nx.basis_rows(W,Nx.Pauli(),1:n,n)
        images = [mod1(i+1,n) for i in 1:n]
        permutation = Nx.SitePermutation(images)
        anti = Nx.AntiunitaryPauli(Nx.pauli_letter(W,n,1))
        canonical_anti = Nx.CliffordMap(anti)
        cases = [(Nx.clifford_h(W,n,1),gate_matrix(:H,n,1),false),
                 (Nx.clifford_s(W,n,n),gate_matrix(:S,n,n),false),
                 (permutation,permutation_matrix(images),false),
                 (anti,pauli_matrix(anti.word,n),true),
                 (canonical_anti,pauli_matrix(anti.word,n),true)]
        n>1 && push!(cases,(Nx.clifford_cnot(W,n,1,n),gate_matrix(:CNOT,n,1,n),false))
        for (g,U,ag) in cases,(h,V,ah) in cases
            combined = Nx.compose(g,h)
            @test Nx.is_antiunitary(combined)==xor(ag,ah)
            C = combined isa Nx.SitePermutation ? Nx.CliffordMap(W,n,combined) : combined
            @test C isa Nx.CliffordMap
            @test isone(Nx.compose(C,inv(C))) && isone(Nx.compose(inv(C),C))
            @test Nx.compose(C,one(C))==C
            @test hash(Nx.compose(C,one(C)))==hash(C)
            @test action_allocations(C,rows[end])==0
            for w in rows
                # Complex scalar detects a lost conjugation flag, even on real Pauli matrices.
                input = (0.7+1.3im)*pauli_matrix(w,n)
                inner = V*(ah ? conj(input) : input)*V'
                expected = U*(ag ? conj(inner) : inner)*U'
                coefficient = Nx.is_antiunitary(C) ? 0.7-1.3im : 0.7+1.3im
                sink = Nx.act!(TermDict{W,ComplexF64}(),C,w,coefficient)
                actual = sum(c*pauli_matrix(v,n) for (v,c) in sink.terms)
                @test actual ≈ expected atol=2e-13
            end
        end
        point = Nx.clifford_h(W,n,1)
        canonical_permutation = Nx.CliffordMap(W,n,permutation)
        @test Nx.compose(inv(point),permutation,point)==Nx.compose(inv(point),canonical_permutation,point)
        @test hash(Nx.compose(inv(point),permutation,point))==hash(Nx.compose(inv(point),canonical_permutation,point))
        @test Nx.compose(anti,anti)==one(canonical_anti)
        @test !Nx.is_antiunitary(Nx.compose(anti,canonical_anti))
        @test Nx.CliffordMap(W,n,permutation)==Nx.CliffordMap(permutation,W,n)
    end
    # Ambient identity padding cannot change equality or hashing.
    h1,h4 = Nx.clifford_h(W,1,1),Nx.clifford_h(W,4,1)
    @test h1==h4 && hash(h1)==hash(h4)
    K = Nx.CliffordMap(Nx.AntiunitaryPauli(one(W)))
    @test !isone(K) && Nx.is_antiunitary(K)
    @test Nx.act(K,Nx.pauli_letter(W,64,1))==(Nx.Phase(2),Nx.pauli_letter(W,64,1))
    @test isone(Nx.compose(K,K))
    @test_throws ArgumentError Nx.CliffordMap(W,1,Nx.SitePermutation([2,1]))
end

@testset "Clifford matrix and legacy oracles" begin
    W = Nx.PauliWord{1}
    for n in 1:4
        rows = Nx.basis_rows(W,Nx.Pauli(),1:n,n)
        gates = [(:H,Nx.clifford_h(W,n,1),1,0),(:S,Nx.clifford_s(W,n,n),n,0)]
        n>1 && append!(gates,[(:CNOT,Nx.clifford_cnot(W,n,n,1),n,1),(:SWAP,Nx.clifford_swap(W,n,1,n),1,n)])
        for (name,g,a,b) in gates
            U = gate_matrix(name,n,a,b)
            oldg = b==0 ? CliffordSymmetry(name,a;nqubits=n,integer_type=UInt16) : CliffordSymmetry(name,a,b;nqubits=n,integer_type=UInt16)
            @test isone(Nx.compose(g,inv(g)))
            @test isone(Nx.compose(inv(g),g))
            @test Nx.compose(g,one(g))==g
            @test hash(Nx.compose(g,one(g)))==hash(g)
            for w in rows
                phase,image = Nx.act(g,w)
                @test Nx.scale(phase,1.0+0im)*pauli_matrix(image,n) ≈ U*pauli_matrix(w,n)*U' atol=1e-13
                sign,oldimage = NCTSSoS._act_monomial(oldg,NormalMonomial{PauliAlgebra,UInt16}(UInt16.(Nx.legacy_letters(w))))
                @test Nx.scale(phase,1.0)==sign && Nx.legacy_letters(image)==Int.(oldimage.word)
            end
            @test action_allocations(g,rows[end])==0
        end
        h,s = Nx.clifford_h(W,n,1),Nx.clifford_s(W,n,1)
        for w in rows
            p,u = Nx.act(h,w); q,v = Nx.act(s,u)
            @test Nx.act(Nx.compose(s,h),w)==(p*q,v)
        end
    end
    @test Nx.degree(Nx.act(Nx.clifford_cnot(W,2,1,2),Nx.pauli_letter(W,1,0))[2])==2
    @test_throws ArgumentError Nx.CliffordMap([Nx.pauli_letter(W,1,0),Nx.pauli_letter(W,1,0)])
    wide = Nx.PauliWord{2}
    g = Nx.clifford_cnot(wide,65,64,65)
    @test Nx.act(g,Nx.pauli_letter(wide,64,0))[2] == Nx.mul(Nx.pauli_letter(wide,64,0),Nx.pauli_letter(wide,65,0))[2]
    @test action_allocations(g,Nx.pauli_letter(wide,64,1))==0
end

@testset "Site and fermionic mode permutations" begin
    g = Nx.SitePermutation([3,1,4,2])
    @test isone(Nx.compose(g,inv(g)))
    @test Nx.compose(g,one(g))==g
    @test hash(Nx.compose(g,one(g)))==hash(g)
    @test_throws ArgumentError Nx.SitePermutation([2,2])
    for (algebra,W) in ((Nx.Pauli(),Nx.PauliWord{1}),(Nx.Fermionic(),Nx.FermiWord{1}),
            ((Nx.SiteAlgebra(rule,2),Nx.word_type(Nx.SiteAlgebra(rule,2),4,6)) for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent()))...)
        rows = Nx.basis_rows(W,algebra,1:4,2)
        for w in rows
            p,v = Nx.act(g,w); q,back = Nx.act(inv(g),v)
            @test p*q==one(Nx.Phase) && back==w
            if W <: Nx.FermiWord
                oldg = FermionicModePermutation(Dict(i=>g.images[i] for i in eachindex(g.images)))
                s,m = NCTSSoS._act_monomial(oldg,NormalMonomial{FermionicAlgebra,Int}(Nx.legacy_letters(w)))
                @test Nx.scale(p,1.0)==s && Nx.legacy_letters(v)==m.word
            end
        end
        @test action_allocations(g,rows[end])==0
        for a in rows[1:8],b in rows[1:8]
            lhs = TermDict{W,ComplexF64}()
            for (w,c) in product_terms(a,b,ComplexF64).terms; Nx.act!(lhs,g,w,c); end
            pa,aa = Nx.act(g,a); pb,bb = Nx.act(g,b)
            rhs = Nx.mul!(TermDict{W,ComplexF64}(),aa,bb,Nx.scale(pa*pb,1.0+0im))
            @test clean(lhs)==clean(rhs)
        end
    end
end

@testset "Antiunitary Pauli action and spin derivations" begin
    W = Nx.PauliWord{1}
    rows = Nx.basis_rows(W,Nx.Pauli(),1:3,3)
    g = Nx.AntiunitaryPauli(Nx.mul(Nx.pauli_letter(W,1,1),Nx.pauli_letter(W,3,0))[2])
    U = pauli_matrix(g.word,3)
    @test Nx.is_antiunitary(g)
    @test isone(Nx.compose(g,inv(g)))
    for w in rows
        p,v = Nx.act(g,w)
        @test Nx.scale(p,1.0+0im)*pauli_matrix(v,3) ≈ U*conj(pauli_matrix(w,3))*U'
    end
    @test action_allocations(g,rows[end])==0
    reg,(x,y,z) = create_pauli_variables(1:3)
    letters = vcat(x,y,z)
    images = (vcat(zero.(Polynomial.(x)),Polynomial.(z),-Polynomial.(y)),
        vcat(-Polynomial.(z),zero.(Polynomial.(y)),Polynomial.(x)),
        vcat(Polynomial.(y),-Polynomial.(x),zero.(Polynomial.(z))))
    for axis in 0:2
        d = Nx.SpinRotation(axis,(UInt64(7),))
        oldd = SuppliedDerivation(reg,Dict(v.word[1]=>image for (v,image) in zip(letters,images[axis+1])))
        @test d==Nx.SpinRotation(axis,(UInt64(7),))
        @test hash(d)==hash(Nx.SpinRotation(axis,(UInt64(7),)))
        for w in rows
            sink = Nx.derive!(TermDict{W,ComplexF64}(),d,w,1.0+0im)
            old = differentiate(oldd,NormalMonomial{PauliAlgebra,eltype(indices(reg))}(eltype(indices(reg)).(Nx.legacy_letters(w))))
            @test bridge_terms(sink)==Dict(Int.(m.word)=>ComplexF64(c) for (c,m) in terms(old) if !iszero(c))
        end
        @test derivation_allocations(d,rows[end])==0
    end
end

function derived_terms(d,p::TermDict{W,C}) where {W,C}
    out=TermDict{W,C}()
    for (w,c) in p.terms
        Nx.derive!(out,d,w,c)
    end
    return out
end
function derived_polynomial(d,p::Nx.Polynomial{W}) where {W}
    out=TermDict{W,ComplexF64}()
    for (c,w) in Nx.terms(p)
        Nx.derive!(out,d,w,c)
    end
    return clean(out)
end

@testset "User derivations and contract 1" begin
    rng=MersenneTwister(423)
    r,(x,y,z)=Nx.pauli_variables(3)
    rotations=Nx.spin_rotations(r)
    @test rotations==Tuple(Nx.SpinRotation(r,axis) for axis in (:x,:y,:z))
    h=sum(x[i]*x[mod1(i+1,3)]+y[i]*y[mod1(i+1,3)]+z[i]*z[mod1(i+1,3)] for i in 1:3)
    xxz=h+sum(z[i]*z[mod1(i+1,3)] for i in 1:3)
    @test all(isempty(derived_polynomial(d,h)) for d in rotations)
    @test map(d->isempty(derived_polynomial(d,xxz)),rotations)==(false,false,true)
    @test !isempty(derived_polynomial(Nx.SpinRotation(r,:z;sites=[1]),h))
    @test isempty(derived_polynomial(Nx.SpinRotation(r,:x;sites=Int[]),h))
    @test_throws r"named spin axis" Nx.SpinRotation(r,:bad)
    @test_throws r"distinct" Nx.SpinRotation(r,:x;sites=[1,1])
    @test_throws r"integer sites" Nx.SpinRotation(r,:x;sites=[0])
    @test_throws r"integer sites" Nx.SpinRotation(r,:x;sites=[1.5])
    @test_throws r"fermionic registry" Nx.NumberRotation(r)
    f,(a,ad)=Nx.fermionic_variables(3)
    @test_throws r"Pauli registry" Nx.SpinRotation(f,:z)
    number=Nx.NumberRotation(f)
    hopping=ad[1]*a[2]+ad[2]*a[1]+ad[2]*a[3]+ad[3]*a[2]
    pair=ad[1]*ad[2]+a[2]*a[1]
    @test isempty(derived_polynomial(number,hopping))
    @test derived_polynomial(number,pair)==Dict(w=>c for (c,w) in Nx.terms(2im*(ad[1]*ad[2]-a[2]*a[1])))
    @test !isempty(derived_polynomial(Nx.NumberRotation(f;modes=[1]),hopping))
    @test isempty(derived_polynomial(Nx.NumberRotation(f;modes=Int[]),pair))

    pauli_rows=Nx.basis_rows(Nx.PauliWord{1},Nx.Pauli(),1:3,3)
    cases=[(d,pauli_rows) for d in Nx.spin_rotations(r;sites=[1,3])]
    for (n,active,selected) in ((4,[1,2,4],[1,4]),(100,[1,64,65,100],[1,65,100]))
        registry,_=Nx.fermionic_variables(n)
        d=Nx.NumberRotation(registry;modes=selected)
        @test d==Nx.NumberRotation(registry;modes=reverse(selected))
        @test hash(d)==hash(Nx.NumberRotation(registry;modes=reverse(selected)))
        @test d!=Nx.NumberRotation(registry;modes=Int[])
        W=Nx.word_type(registry,1)
        K=cld(n,64)
        masks=[ntuple(k->sum((UInt64(1)<<((site-1)%64) for site in active
            if (site-1)÷64+1==k && rand(rng,Bool));init=UInt64(0)),K) for _ in 1:80]
        words=[Nx.FermiWord(rand(rng,masks),rand(rng,masks)) for _ in 1:100]
        push!(words,Nx.fermi_letter(W,-last(selected)))
        for w in words
            charge=-sum((sign(l) for l in Nx.legacy_letters(w) if abs(l) in selected);init=0)
            got=Nx.derive!(TermDict{W,ComplexF64}(),d,w,2-3im)
            @test clean(got)==(charge==0 ? Dict{W,ComplexF64}() : Dict(w=>im*charge*(2-3im)))
            @test derivation_allocations(d,w)==0
        end
        # Keep the cases heterogeneous: the tests must exercise both mask widths.
        for _ in 1:100
            u,v=rand(rng,words),rand(rng,words)
            du=Nx.derive!(TermDict{W,ComplexF64}(),d,u,1)
            dv=Nx.derive!(TermDict{W,ComplexF64}(),d,v,1)
            rhs=times(du,v)
            for (w,c) in times(u,dv).terms; Nx.push_term!(rhs,c,w); end
            @test clean(derived_terms(d,product_terms(u,v,ComplexF64)))==clean(rhs)
            single=TermDict{W,ComplexF64}(Dict(u=>2-3im))
            @test clean(derived_terms(d,star_terms(single)))==clean(star_terms(derived_terms(d,single)))
        end
    end
    for (d,words) in cases, _ in 1:100
        u,v=rand(rng,words),rand(rng,words)
        W=typeof(u)
        du=Nx.derive!(TermDict{W,ComplexF64}(),d,u,1)
        dv=Nx.derive!(TermDict{W,ComplexF64}(),d,v,1)
        rhs=times(du,v)
        for (w,c) in times(u,dv).terms; Nx.push_term!(rhs,c,w); end
        @test clean(derived_terms(d,product_terms(u,v,ComplexF64)))==clean(rhs)
        single=TermDict{W,ComplexF64}(Dict(u=>2-3im))
        @test clean(derived_terms(d,star_terms(single)))==clean(star_terms(derived_terms(d,single)))
        @test derivation_allocations(d,u)==0
    end
end
