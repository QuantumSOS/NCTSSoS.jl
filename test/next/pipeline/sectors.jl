function sector_matrix(p,n)
    sum(c*physical_pauli(w,1:n) for (c,w) in Nx.terms(p))
end

@testset "Continuous sector keeps the base equations" begin
    registry,(a,ad)=Nx.fermionic_variables(1)
    number=ad[1]*a[1]
    built=Nx.build(Nx.Problem(number;symmetry=(Nx.NumberRotation(registry),)),
        Nx.Relaxation(order=1,strengthen=Nx.Sector(number,1)))
    form=Nx.moment_standard_form(built.program)
    # Frozen from redesign 5a794cb; exact support splitting only permutes
    # the occupation, vacancy and normalization cones. The sector stays explicit.
    @test sort(collect(zip(vec(Matrix(form.A))[1:3],form.b[1:3])))==
        sort([(-1.0,0.0),(1.0,1.0),(0.0,1.0)])
    @test (form.A[4,1],form.b[4])==(1.0,1.0)
    @test form.c==[1.0]
    @test form.offset==0.0
    @test [(c.kind,c.side,c.first,c.length,c.label) for c in form.cones]==
        [(:hpsd,1,1,1,0),(:hpsd,1,2,1,0),(:hpsd,1,3,1,0),(:zero,3,4,1,1)]
    report=built.sizes.closure.sector_elimination
    @test !report.accepted && report.reason==:continuous_symmetry
    @test report.rank===nothing && report.coefficients_before===nothing
    @test report.variables==report.initial_variables==built.sizes.real_variables
    eliminated=Nx.build(Nx.Problem(number),Nx.Relaxation(order=1,strengthen=Nx.Sector(number,1)))
    accepted=eliminated.sizes.closure.sector_elimination
    @test accepted.accepted && accepted.reason==:accepted
    @test keys(report)==keys(accepted)==(:accepted,:reason,:rank,:coefficients_before,
        :coefficients_after,:arithmetic,:initial_variables,:variables)
end

# Check the same physical equations whether they remain explicit or have
# become the affine moment map. Undetermined moments are not queried.
function sector_state_residual(built,physical)
    part(label,z)=label.part==2 ? imag(z) : real(z)
    y=[part(label,physical(label.word)) for label in built.table.labels]
    residual=maximum(abs,equality_values(built.program,y,1);init=0.0)
    if built.table isa Nx.SectorMoments
        for label in built.table.base.labels
            coefficients=Dict{Int,ComplexF64}()
            for (id,c) in Nx.moment_terms(built.table,label.word)
                coefficients[id]=get(coefficients,id,0im)+c
            end
            all(iszero(c) || Nx.moment_covered(built.table,id) for (id,c) in coefficients) || continue
            value=sum((c*y[id] for (id,c) in coefficients if !iszero(c));init=0im)
            residual=max(residual,abs(part(label,value)-part(label,physical(label.word))))
        end
    end
    residual
end

function sector_reference_map(table,labels)
    pairs=[collect(Nx.moment_terms(table,label.word)) for label in labels]
    width=maximum((id for row in pairs for (id,c) in row);init=length(table.labels))
    M=zeros(Float64,length(labels),width)
    for (i,(label,row)) in enumerate(zip(labels,pairs)),(id,c) in row
        M[i,id]+=label.part==2 ? imag(c) : real(c)
    end
    M
end

@testset "Sector overflow preserves the exact program" begin
    _,(b,bd)=Nx.bosonic_variables(1)
    number=bd[1]*b[1]
    Q=bd[1]^12*b[1]^12
    # The last contraction in Q*(bd^60*b^60) exceeds Int64. Scaling Q
    # delays overflow until normalization; a still smaller scale cannot fit
    # the variance's denominator. Every retry must preserve the equations,
    # including the already completed number-sector block.
    @test prod(BigInt.(49:60))>typemax(Int64)
    family=[[one(b[1]),b[1]^60]]
    forms=map((1,2.0^-16,2.0^-40)) do factor
        built=Nx.build(Nx.Problem(zero(number)),Nx.BlockRelaxation(family=family,
            strengthen=(Nx.Sector(number,0),Nx.Sector(factor*Q,0))))
        @test built.sizes.closure.sector_elimination.arithmetic==Rational{BigInt}
        Nx.moment_standard_form(built.program)
    end
    for form in forms[2:end]
        @test isequal(forms[1].A,form.A)
        @test isequal(forms[1].b,form.b)
        @test isequal(forms[1].c,form.c)
        @test isequal(forms[1].offset,form.offset)
        @test isequal(forms[1].cones,form.cones)
    end
end

@testset "Sector elimination does not increase cone fill" begin
    _,(x,y,z)=Nx.pauli_variables(4)
    H=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+z[i]*z[mod1(i+1,4)] for i in 1:4)
    plain=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,split=Nx.TermSplit()))
    built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,split=Nx.TermSplit(),
        strengthen=(Nx.Sector(sum(z),0),)))
    report=built.sizes.closure.sector_elimination
    @test report.coefficients_after>report.coefficients_before
    @test !report.accepted && report.reason==:coefficient_growth
    @test keys(report)==(:accepted,:reason,:rank,:coefficients_before,
        :coefficients_after,:arithmetic,:initial_variables,:variables)
    # Keep the original cone entries and explicit sector equations, rather
    # than accepting a denser substitution with the same bound.
    count=length(plain.program.forms.ids)
    @test built.program.forms.ids[1:count]==plain.program.forms.ids
    @test built.program.forms.coefs[1:count]==plain.program.forms.coefs
    @test any(b.kind==Nx.ZERO for b in built.program.blocks)
end

@testset "Eliminated sector moments retain coverage" begin
    _,(b,bd)=Nx.bosonic_variables(2)
    n1,n2=bd.*b
    number=n1+n2
    problem=Nx.Problem(n1)
    relaxation=Nx.Relaxation(order=1,strengthen=Nx.Sector(number,1))
    built=Nx.build(problem,relaxation)
    result=Nx.solve(built;solver=Nx.JuMPBackend(SOLVER))
    @test result.bound≈0 atol=1e-6
    @test Nx.expectation(result,n2)≈1 atol=1e-6
    @test Nx.expectation(result,number)≈1 atol=1e-6
    # These degree-four moments have no PSD constraint at order one. The
    # sector relation is readable only after their free parts cancel.
    @test Nx.expectation(result,n1*(number-1))≈0 atol=1e-6
    @test_throws r"not among the moments of the selected blocks" Nx.expectation(result,n1*n2)
    pairs=collect(Nx.moment_terms(built.table,argmax(Nx.degree,(w for (c,w) in Nx.terms(n1*n2)))))
    @test any(!Nx.moment_covered(built.table,id) for (id,c) in pairs if !iszero(c))
    @test built.sizes.real_variables==built.program.nvars-1
    @test built.sizes.blocks==length(built.program.blocks)
    @test_throws r"Sector equations contradict the normalized identity" Nx.build(
        Nx.Problem(n1),Nx.Relaxation(order=1,strengthen=Nx.Sector(one(n1),2)))
end

@testset "DensityMatrix uses the contained sector face" begin
    _,(x,y,z)=Nx.pauli_variables(2)
    H=z[2]+0.2y[1]*y[2]
    # A non-diagonal sector whose exact kernel basis need not be orthonormal.
    Q=3x[1]+4z[1]
    kinds=(Nx.DensityMatrix(1:2),Nx.Sector(Q,5))
    built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=kinds))
    report=built.sizes.closure.strengthening
    @test only(report.density_faces).before==4
    @test only(report.density_faces).after==2
    @test occursin("density sector face 4 → 2",sprint(show,built.sizes))
    @test strengthening_bound(built.program)≈-1 atol=1e-6
    e=eigen(Hermitian(sector_matrix(Q,2)))
    V=e.vectors[:,findall(v->abs(v-5)<1e-10,e.values)]
    outside=e.vectors[:,findfirst(v->abs(v+5)<1e-10,e.values)]
    for seed in 1:3
        A=randn(MersenneTwister(seed),ComplexF64,2,2)
        rho=V*A*A'*V'; rho/=LinearAlgebra.tr(rho)
        coords=physical_coordinates(built.table,rho,2)
        @test sector_state_residual(built,w->LinearAlgebra.tr(rho*physical_pauli(w,1:2)))<1e-12
        @test minimum(extra_values(built.program,coords,report.first_block))>=-1e-10
    end
    # Removing off-face entries without equations would accept this matrix:
    # it has the same compressed density as the valid mixed-sector state.
    rho=V*V'/2+0.1(outside*V[:,1]'+V[:,1]*outside')
    @test sector_state_residual(built,w->LinearAlgebra.tr(rho*physical_pauli(w,1:2)))>0.01
    local_density=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.DensityMatrix([1]),Nx.Sector(z[1]+z[2],0))))
    @test isempty(local_density.sizes.closure.strengthening.density_faces)
    trivial=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.DensityMatrix(1:2),Nx.Sector(one(H),1))))
    @test only(trivial.sizes.closure.strengthening.density_faces).after==4
    joint=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(kinds...,Nx.Sector(z[2],1))))
    @test only(joint.sizes.closure.strengthening.density_faces).after==1
    @test_throws r"Sector 2.*not an eigenvalue.*choose" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.DensityMatrix(1:2),Nx.Sector(Q,4))))
    @test_throws r"empty common kernel" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(kinds...,Nx.Sector(Q,-5))))
    _,(b,bd)=Nx.bosonic_variables(1)
    @test_throws r"bosonic variables use Sector" Nx.build(Nx.Problem(bd[1]*b[1]),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.DensityMatrix([1]),)))
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
            @test sector_state_residual(built,w->LinearAlgebra.tr(rho*physical_pauli(w,1:2)))<1e-12
        end
        outside=ComplexF64[1,0,0,0]
        @test sector_state_residual(built,w->LinearAlgebra.tr(outside*outside'*physical_pauli(w,1:2)))>0.1
        @test strengthening_bound(built.program)≈minimum(eigvals(Hermitian(V'*sector_matrix(H,2)*V))) atol=1e-6
        @test occursin("eigenspace",sprint(show,built.valid_for))
    end
    # At order zero the mean number alone accepts the vacuum/full mixture.
    b=Nx.build(Nx.Problem(zero(H)),Nx.BlockRelaxation(order=0,strengthen=(Nx.Sector(N,1),)))
    check_used_coordinates(b.program)
    rho=diagm(ComplexF64[0.5,0,0,0.5])
    @test sector_state_residual(b,w->LinearAlgebra.tr(rho*physical_pauli(w,1:2)))>0.1
    @test_throws r"not Hermitian" Nx.Sector(im*x[1],0)
    @test_throws r"finite" Nx.Sector(N,Inf)
    _,(_,_,other)=Nx.pauli_variables(2)
    @test_throws r"registry" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(other[1],1),)))
    @test_throws r"exactly zero commutator.*remove Stationarity or choose a commuting sector operator" Nx.build(Nx.Problem(H+1e-14*x[1]),Nx.BlockRelaxation(order=1,strengthen=(Nx.Sector(N,1),Nx.Stationarity())))
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
        @test built.sizes.closure.sector_elimination.rank>0
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
        coordinates=sector_reference_map(built.table,reference.table.labels)
        columns=[Float64[i==j for i in 1:size(coordinates,2)] for j in 1:size(coordinates,2)]
        A=hcat([equality_values(built.program,v,built.sizes.closure.strengthening.first_block) for v in columns]...)
        B=hcat([equality_values(reference.program,coordinates*v,start) for v in columns]...)
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
    coordinates=sector_reference_map(built.table,reference.table.labels)
    columns=[Float64[i==j for i in 1:size(coordinates,2)] for j in 1:size(coordinates,2)]
    A=hcat([equality_values(built.program,v,built.sizes.closure.strengthening.first_block) for v in columns]...)
    B=hcat([equality_values(reference.program,coordinates*v,start) for v in columns]...)
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
        @test sector_state_residual(built,w->LinearAlgebra.tr(rho*wordmatrix(w)))<1e-12
    end
    @test strengthening_bound(built.program)≈-2 atol=1e-6
end
