function torus_actions(nx,ny)
    site(x,y)=mod(x,nx)+nx*mod(y,ny)+1
    perm(f)=[site(f(x,y)...) for y in 0:(ny-1) for x in 0:(nx-1)]
    tx=test_map(nx*ny;sites=perm((x,y)->(x+1,y)))
    ty=test_map(nx*ny;sites=perm((x,y)->(x,y+1)))
    reflect=test_map(nx*ny;sites=perm((x,y)->(x,-y)))
    point=[reflect,test_map(nx*ny;axes=(1,2,0)),test_map(nx*ny;axes=(1,0,2),signs=(1,-1,1))]
    nx==ny && push!(point,test_map(nx*ny;sites=perm((x,y)->(-y,x))))
    tx,ty,point
end
function torus_objective(nx,ny)
    site(x,y)=mod(x,nx)+nx*mod(y,ny)+1
    [(1.,Nx.mul(toy_letter(site(x,y),a),toy_letter(site(x+dx,y+dy),a))[2])
        for y in 0:(ny-1) for x in 0:(nx-1) for (dx,dy) in ((1,0),(0,1)) for a in 0:2]
end

@testset "Ambiguous antiunitary fixed basis preserves the Hermitian sector" begin
    # Ordinary K on the real X,Z operator span. This change of basis is
    # unitary, but its first real/imaginary candidates are nearly parallel.
    θ=3e-9
    Q=cis(0.7)*[cos(θ) im*sin(θ); im*sin(θ) cos(θ)]
    @test Q'*Q ≈ LinearAlgebra.I atol=1e-15
    @test isnothing(Nx._fixed_real_columns(Q,conj(Q)))

    W=Nx.PauliWord{1}
    rows=sort([Nx.pauli_letter(W,1,0),Nx.pauli_letter(W,1,2)])
    K=Nx.CliffordMap(Nx.AntiunitaryPauli(one(W)))
    g=one(K)
    translations=Nx.TranslationGroup([g],[1])
    orbits=Nx.AbelianRows(rows,translations)
    index=Dict(w=>(a,1) for (a,w) in enumerate(rows))
    sectors=[[[(Q[i,j],i) for i in 1:2] for j in 1:2]]
    retained,real_blocks=Nx._real_character_sectors(sectors,orbits,index,1,K,[g])
    @test retained==sectors
    @test real_blocks==[false]

    # A well-conditioned basis still succeeds and has the same projector.
    U=ComplexF64[0 im; 1 0; 0 0]
    F=hcat(Nx._fixed_real_columns(U,conj(U))...)
    @test F==[0 1; 1 0; 0 0]
    @test F*F' ≈ U*U' atol=1e-15
end

@testset "Fourier assembly retries a nonreal fixed basis" begin
    W=Nx.PauliWord{1}
    identity,x=one(W),Nx.pauli_letter(W,1,0)
    K=Nx.CliffordMap(Nx.AntiunitaryPauli(identity))
    translations=Nx.TranslationGroup([one(K)],[1])
    orbits=Nx.AbelianRows([identity,x],translations)
    original=[[[ (1.0+0im,1)],[(0.0+1im,2)]], [[(1.0+0im,1)]]]
    candidate=deepcopy(original)
    # These approximate K-fixed rows have overlap -2e-10im on identity:
    # below the fixed-basis norm budget, above the affine coefficient budget.
    candidate[1][2]=[(1.0+0im,2),(0.0+2e-10im,1)]
    table=Nx.OrbitMoments(W,(K,))
    p=Nx.MomentProgram{ComplexF64}()
    labels=Tuple{Int,Int}[]
    Nx._assemble_character!(p,Nx.OrbitFormSink(table,Dict{Int32,ComplexF64}()),
        orbits,ComplexF64[1],1,candidate,labels,Dict();real_blocks=[true,true],original_sectors=original)
    p.nvars=length(table.labels)
    @test Nx.validate(p)===p
    @test labels==[(1,1),(1,2)]
    @test candidate==original
    @test [(b.kind,b.side) for b in p.blocks]==[(Nx.HPSD,2),(Nx.PSD,1)]
    # Direct Pauli products: the original lower entry is (iX)'I = -iX.
    @test evaluated_blocks(p,[1.0,0.375])==[ComplexF64[1 0.375im; -0.375im 1],ones(ComplexF64,1,1)]
end

@testset "Product translations and point group" begin
    for (nx,ny,d) in ((2,2,1),(2,3,1),(2,2,2))
        tx,ty,point=torus_actions(nx,ny)
        translations=Nx.TranslationGroup([tx,ty],[nx,ny])
        rows=toy_words(nx*ny,d); objective=torus_objective(nx,ny)
        p,table,report=Nx.abelian_lattice_program(rows,objective,translations,point)
        dense,dt=dense_program(rows,objective)
        @test report.translation_order==nx*ny
        @test solve_reference(p)[1]≈solve_reference(dense)[1] atol=1e-6
        real_program,_,_=Nx.abelian_lattice_program(rows,objective,translations,point;moment_generators=(test_conjugate(nx*ny),))
        @test [b.side for b in real_program.blocks]==[b.side for b in p.blocks]
        @test all(b.kind==Nx.PSD for b in real_program.blocks)
        @test solve_reference(real_program)[1]≈solve_reference(dense)[1] atol=1e-6
        relations=[Nx.orbit_coordinates!(table,l.word) for l in dt.labels]
        y=randn(MersenneTwister(911),length(table.labels)); y[1]=1
        yd=[real(coordinate_value(m,y)) for m in relations]
        @test eigmin(Hermitian(only(evaluated_blocks(dense,yd))))≈minimum(eigmin(Hermitian(A)) for A in evaluated_blocks(p,y)) atol=1e-9
    end
    tx,ty,point=torus_actions(3,3)
    @test_throws ArgumentError Nx.TranslationGroup([tx,tx],[3,3])
    @test_throws ArgumentError Nx.TranslationGroup([ty,point[1]],[3,2])
    nx,ny=2,3; n=nx*ny; W=Nx.word_type(Nx.Pauli(),n,2)
    tx,ty,point=torus_actions(nx,ny)
    permutation(g)=Nx.SitePermutation([only(Nx.sites(Nx.act(g,toy_letter(i,0))[2])) for i in 1:n])
    translations=Nx.TranslationGroup([permutation(tx),permutation(ty)],[nx,ny])
    points=Nx.AbstractAction[permutation(first(point));[core_action(g,W) for g in point[2:end]]]
    rows=Nx.basis_rows(W,Nx.Pauli(),1:n,1)
    objective=[(c,Nx.word_from_legacy(W,Nx.legacy_letters(w))) for (c,w) in torus_objective(nx,ny)]
    p,_,_=Nx.abelian_lattice_program(rows,objective,translations,points;
        moment_generators=(Nx.AntiunitaryPauli(one(W)),))
    plain,_=Nx.plain_program(rows,objective)
    @test all(b.kind==Nx.PSD for b in p.blocks)
    @test solve_reference(p)[1]≈solve_reference(plain)[1] atol=1e-6
end
