function sector_matrix(p,n)
    sum(c*physical_pauli(w,1:n) for (c,w) in Nx.terms(p))
end

@testset "Sector states, not operator identities" begin
    _,(x,y,z)=Nx.pauli_variables(2)
    N=(2-z[1]-z[2])/2
    S2=((x[1]+x[2])^2+(y[1]+y[2])^2+(z[1]+z[2])^2)/4
    H=x[1]*x[2]+y[1]*y[2]+0.5*z[1]*z[2]
    for (Q,q) in ((N,1),(S2,0))
        built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(Q,q),)))
        start=built.sizes.closure.strengthening.first_block
        M=sector_matrix(Q,2)
        decomposition=eigen(Hermitian(M))
        V=decomposition.vectors[:,findall(v->abs(v-q)<1e-10,decomposition.values)]
        @test size(V,2)==(Q===N ? 2 : 1)
        for seed in 1:4
            A=randn(MersenneTwister(seed),ComplexF64,size(V,2),size(V,2))
            rho=V*A*A'*V'; rho/=LinearAlgebra.tr(rho)
            coords=physical_coordinates(built.table,rho,2)
            @test maximum(abs,equality_values(built.program,coords,start))<1e-12
        end
        outside=ComplexF64[1,0,0,0]
        @test maximum(abs,equality_values(built.program,physical_coordinates(built.table,outside*outside',2),start))>0.1
        @test strengthening_bound(built.program)≈minimum(eigvals(Hermitian(V'*sector_matrix(H,2)*V))) atol=1e-6
        @test occursin("eigenspace",sprint(show,built.valid_for))
    end
    # At order zero the mean number alone accepts the vacuum/full mixture.
    b=Nx.build(Nx.Problem(zero(H)),Nx.BlockRelaxation(order=0,strengthen=(Nx.Sector(N,1),)))
    check_used_coordinates(b.program)
    rho=diagm(ComplexF64[0.5,0,0,0.5])
    @test maximum(abs,equality_values(b.program,physical_coordinates(b.table,rho,2),b.sizes.closure.strengthening.first_block))>0.1
    @test_throws r"not Hermitian" Nx.Sector(im*x[1],0)
    @test_throws r"finite" Nx.Sector(N,Inf)
    _,(_,_,other)=Nx.pauli_variables(2)
    @test_throws r"registry" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(other[1],1),)))
    @test_throws r"exactly zero commutator" Nx.build(Nx.Problem(H+1e-14*x[1]),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(N,1),Nx.Stationarity())))
    combined=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(N,1),Nx.Stationarity())))
    @test strengthening_bound(combined.program)≈-2.5 atol=1e-6
    @test_throws r"StateOptimality with a sector" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(N,1),Nx.StateOptimality())))
end

@testset "Four-site sector symmetry guard" begin
    _,(x,y,z)=Nx.pauli_variables(4)
    N=(4-sum(z))/2
    H=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+0.5*z[i]*z[mod1(i+1,4)] for i in 1:4)
    W=typeof(first(x[1].data.words))
    flip=Nx.PauliConjugation(Nx.word_from_legacy(W,[1,4,7,10]),false)
    for symmetry in (Nx.FamilySymmetry([flip]),Nx.FamilySymmetry(;moment_only=(flip,)))
        @test_throws r"remove the generator or choose an invariant sector" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,symmetry=symmetry,strengthen=(Nx.Sector(N,1),)))
        built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,symmetry=symmetry,strengthen=(Nx.Sector(N,2),)))
        check_used_coordinates(built.program)
        @test built.sizes.closure.strengthening.added_equalities>0
    end
end

@testset "Fermionic particle-number sector" begin
    _,(a,ad)=Nx.fermionic_variables(2)
    N=sum(ad[i]*a[i] for i in 1:2)
    H=ad[1]*a[2]+ad[2]*a[1]+0.3*N
    registry,(old_a,old_ad)=create_fermionic_variables(1:2)
    old_H=old_ad[1]*old_a[2]+old_ad[2]*old_a[1]+0.3sum(old_ad[i]*old_a[i] for i in 1:2)
    for number in (0,1,2)
        built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(N,number),)))
        # Exact Fock-space energies: vacuum, one-particle hopping, full.
        bound=strengthening_bound(built.program)
        @test bound≈(0.0,-0.7,0.6)[number+1] atol=1e-6
        old=cs_nctssos(polyopt(old_H,registry;moment_eq_constraints=particle_number_constraint(registry,number)),
            SolverConfig(optimizer=SOLVER,order=1);dualize=false)
        @test bound≈old.objective atol=1e-6
    end
end

# Independent matrix products, expanded in the complete orthogonal Pauli
# basis. Both real and imaginary equations are needed for a one-sided product.
function independent_sector!(reference,Q,q,n,order)
    W=eltype(reference.table.representatives.keys)
    words=Nx.basis_rows(W,Nx.Pauli(),1:n,n)
    matrices=[physical_pauli(w,1:n) for w in words]
    shifted=sector_matrix(Q,n)-q*I
    products=[physical_pauli(a,1:n)'*sum(c*physical_pauli(w,1:n) for (c,w) in block.polynomial)*physical_pauli(b,1:n)
        for block in reference.sizes.closure.family for a in block.rows for b in block.rows]
    domain=[P for (w,P) in zip(words,matrices) if Nx.degree(w)<=2order &&
        any(abs(LinearAlgebra.tr(P*M))>1e-10 for M in products)]
    equations=Matrix{ComplexF64}[shifted^2]
    append!(equations,[shifted*P for P in domain])
    append!(equations,[P*shifted for P in domain])
    for M in equations
        physical_form!(reference.program,reference.table,(M+M')/2,words,matrices)
        physical_form!(reference.program,reference.table,(M-M')/(2im),words,matrices)
    end
    Nx.push_block!(reference.program,Nx.ZERO,2length(equations))
    reference.program.nvars=length(reference.table.labels)
end

@testset "Sector reduced and independent equality spaces" begin
    _,(x,y,z)=Nx.pauli_variables(4)
    N=(4-sum(z))/2
    H=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+0.5*z[i]*z[mod1(i+1,4)] for i in 1:4)
    problem=Nx.Problem(H)
    generators=[Nx.SitePermutation([2,3,4,1]),Nx.SitePermutation([1,4,3,2])]
    for symmetry in (Nx.FamilySymmetry(),Nx.FamilySymmetry(generators)),
        selection in ((;),(;grouping=Nx.Correlative(),split=Nx.TermSplit()))
        built=Nx.build(problem,Nx.BlockRelaxation(;order=1,symmetry,selection...,strengthen=(Nx.Sector(N,1),)))
        reference=Nx.build(problem,Nx.BlockRelaxation(;order=1,symmetry,selection...))
        start=length(reference.program.blocks)+1
        independent_sector!(reference,N,1,4,1)
        coordinates=[Nx.orbit_coordinates!(built.table,label.word) for label in reference.table.labels]
        built.program.nvars=length(built.table.labels)
        columns=[Float64[i==j for i in 1:length(built.table.labels)] for j in 1:length(built.table.labels)]
        A=hcat([equality_values(built.program,v,built.sizes.closure.strengthening.first_block) for v in columns]...)
        B=hcat([equality_values(reference.program,real.([coordinate_value(c,v) for c in coordinates]),start) for v in columns]...)
        @test rank(A;atol=1e-9)==rank(B;atol=1e-9)==rank(vcat(A,B);atol=1e-9)
        @test strengthening_bound(built.program)≈strengthening_bound(reference.program) atol=1e-6
    end
end

@testset "Sector includes moments supplied by localizers" begin
    _,(x,y,z)=Nx.pauli_variables(4)
    N=(4-sum(z))/2
    parity=prod(z)
    problem=Nx.Problem(sum(z);constraints=(Nx.Positive(1+parity),))
    # These row products do not contain parity. Only its localizer does.
    family=[[one(z[1]);z;x[1]*x[2]]]
    built=Nx.build(problem,Nx.BlockRelaxation(family=family,strengthen=(Nx.Sector(N,1),)))
    check_used_coordinates(built.program)
    reference=Nx.build(problem,Nx.BlockRelaxation(family=family))
    start=length(reference.program.blocks)+1
    independent_sector!(reference,N,1,4,2)
    coordinates=[Nx.orbit_coordinates!(built.table,label.word) for label in reference.table.labels]
    built.program.nvars=length(built.table.labels)
    columns=[Float64[i==j for i in 1:length(built.table.labels)] for j in 1:length(built.table.labels)]
    A=hcat([equality_values(built.program,v,built.sizes.closure.strengthening.first_block) for v in columns]...)
    B=hcat([equality_values(reference.program,real.([coordinate_value(c,v) for c in coordinates]),start) for v in columns]...)
    @test rank(A;atol=1e-9)==rank(B;atol=1e-9)==rank(vcat(A,B);atol=1e-9)
end

@testset "Fermionic singlet states from Jordan-Wigner matrices" begin
    n=4
    _,(a,ad)=Nx.fermionic_variables(n)
    N=sum(ad[i]*a[i] for i in 1:n)
    Sz=(ad[1]*a[1]-ad[2]*a[2]+ad[3]*a[3]-ad[4]*a[4])/2
    Sp=ad[1]*a[2]+ad[3]*a[4]
    S2=Sz^2+(Sp*Sp'+Sp'*Sp)/2
    H=-(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])
    ladder=[foldl(kron,(j<i ? ComplexF64[1 0;0 -1] : j==i ? ComplexF64[0 1;0 0] : Matrix{ComplexF64}(I,2,2) for j in 1:n)) for i in 1:n]
    wordmatrix(w)=foldl(*,(l<0 ? ladder[-l]' : ladder[l] for l in Nx.legacy_letters(w));init=Matrix{ComplexF64}(I,2^n,2^n))
    matrix(p)=sum(c*wordmatrix(w) for (c,w) in Nx.terms(p))
    built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(N,2),Nx.Sector(S2,0))))
    number=matrix(N); spin=matrix(S2)
    # The sum of squares selects their joint kernel independently of the
    # normal-order product and the sector implementation.
    e=eigen(Hermitian((number-2I)^2+spin^2))
    V=e.vectors[:,findall(v->abs(v)<1e-10,e.values)]
    @test size(V,2)==3
    for seed in 1:3
        A=randn(MersenneTwister(seed),ComplexF64,3,3)
        rho=V*A*A'*V'; rho/=LinearAlgebra.tr(rho)
        coords=[label.part==2 ? imag(LinearAlgebra.tr(rho*wordmatrix(label.word))) : real(LinearAlgebra.tr(rho*wordmatrix(label.word))) for label in built.table.labels]
        @test maximum(abs,equality_values(built.program,coords,built.sizes.closure.strengthening.first_block))<1e-12
    end
    @test strengthening_bound(built.program)≈-2 atol=1e-6
end
