# Explicit Jordan-Wigner matrices, independent of normal-order multiplication.
function fermi_matrices(n)
    identity=Matrix{ComplexF64}(I,2^n,2^n)
    ladders=[foldl(kron,(j<i ? ComplexF64[1 0;0 -1] : j==i ? ComplexF64[0 1;0 0] : Matrix{ComplexF64}(I,2,2)
        for j in 1:n)) for i in 1:n]
    word(w)=foldl(*,(l<0 ? ladders[-l]' : ladders[l] for l in Nx.legacy_letters(w));init=identity)
    polynomial(p)=sum((c*word(w) for (c,w) in Nx.terms(p));init=zeros(ComplexF64,size(identity)))
    (;identity,ladders,word,polynomial)
end
function fermi_coordinates(table,rho,matrices)
    [label.part==2 ? imag(LinearAlgebra.tr(rho*matrices.word(label.word))) :
        real(LinearAlgebra.tr(rho*matrices.word(label.word))) for label in table.labels]
end

function fermi_condition_matrix(condition,modes,rho,matrices)
    a=matrices.ladders
    pairs=[a[p]*a[q] for p in modes for q in modes if p<q]
    cubic=condition in (:T1,:T2,:T2prime)
    rows=condition==:P ? pairs : condition==:Q ? adjoint.(pairs) :
        condition==:G ? [[matrices.identity];[a[p]'*a[q] for p in modes for q in modes]] :
        condition==:T1 ? [a[p]*a[q]*a[r] for p in modes for q in modes for r in modes if p<q<r] :
        [a[r]'*a[p]*a[q] for p in modes for q in modes for r in modes if p<q]
    ntriples=length(rows)
    condition==:T2prime && append!(rows,a[modes])
    [LinearAlgebra.tr(rho*(rows[i]'*rows[j]+(cubic && i<=ntriples && j<=ntriples ? rows[j]*rows[i]' : zero(matrices.identity))))
        for i in eachindex(rows),j in eachindex(rows)]
end

@testset "Fermionic positivity: independent pure and mixed states" begin
    n=4
    _,(a,ad)=Nx.fermionic_variables(n)
    matrices=fermi_matrices(n)
    for condition in (:P,:Q,:G,:T1,:T2,:T2prime)
        built=Nx.build(Nx.Problem(zero(a[1])),Nx.BlockRelaxation(order=1,
            strengthen=(Nx.FermionicConditions(1:n;conditions=(condition,)),)))
        report=built.sizes.closure.strengthening
        @test report.added_blocks==1
        @test occursin("valid for each state",sprint(show,built.valid_for))
        @test maximum(Nx.degree(label.word) for label in built.table.labels)<=4
        check_used_coordinates(built.program)
        for rank in (1,5),seed in 1:3
            A=randn(MersenneTwister(seed),ComplexF64,2^n,rank)
            rho=A*A'; rho/=LinearAlgebra.tr(rho)
            values=fermi_coordinates(built.table,rho,matrices)
            spectrum=extra_values(built.program,values,report.first_block)
            reference=fermi_condition_matrix(condition,1:n,rho,matrices)
            @test spectrum≈eigvals(Hermitian(reference)) atol=2e-12
            @test minimum(spectrum)>=-2e-12
        end
    end
end

@testset "Each fermionic condition strengthens a finite lower bound" begin
    _,(a,ad)=Nx.fermionic_variables(3)
    numbers=[ad[i]*a[i] for i in 1:3]; N=sum(numbers)
    n1,n2,n3=numbers
    t2=n1*n2+n3-n1*n3-n2*n3
    cases=((:P,n1*n2,(:Q,:G)),(:Q,(1-n1)*(1-n2),(:P,:G)),
        (:G,(N-1.5)^2,(:P,:Q)),(:T1,1-N+n1*n2+n1*n3+n2*n3,(:P,:Q,:G)),
        (:T2,t2,(:P,:Q,:G)),(:T2prime,t2,(:P,:Q,:G)))
    matrices=fermi_matrices(3)
    # Number-one subspace, obtained independently from its explicit matrix.
    sector=findall(v->abs(v-1)<1e-12,diag(matrices.polynomial(N)))
    for (condition,H,others) in cases
        problem=Nx.Problem(H)
        before=Nx.build(problem,Nx.BlockRelaxation(order=1,
            strengthen=(Nx.FermionicConditions(1:3;conditions=others),)))
        after=Nx.build(problem,Nx.BlockRelaxation(order=1,
            strengthen=(Nx.FermionicConditions(1:3;conditions=(others...,condition)),)))
        plain,bound=strengthening_bound(before.program),strengthening_bound(after.program)
        matrix=matrices.polynomial(H)
        @test bound>plain+0.05
        @test bound<=eigmin(Hermitian(matrix))+1e-6
        restricted=Nx.build(problem,Nx.BlockRelaxation(order=1,
            strengthen=(Nx.FermionicConditions(1:3;conditions=(others...,condition)),Nx.Sector(N,1))))
        @test strengthening_bound(restricted.program)<=eigmin(Hermitian(matrix[sector,sector]))+1e-6
        # The same blocks and the sector equations accept pure and mixed states.
        for rank in (1,3)
            A=randn(MersenneTwister(rank),ComplexF64,length(sector),rank)
            rho=zeros(ComplexF64,8,8); rho[sector,sector]=A*A'; rho/=LinearAlgebra.tr(rho)
            values=fermi_coordinates(restricted.table,rho,matrices)
            @test maximum(abs,equality_values(restricted.program,values,1))<1e-12
            for block in restricted.program.blocks
                block.kind==Nx.ZERO && continue
                forms=[sum((restricted.program.forms.coefs[t]*values[restricted.program.forms.ids[t]]
                    for t in Nx.term_range(restricted.program.forms,f));init=0im) for f in Nx.form_range(block)]
                M=zeros(ComplexF64,block.side,block.side)
                for j in 1:block.side,i in j:block.side
                    M[i,j]=forms[Nx.tri_offset(block.side,i,j)+1]; M[j,i]=conj(M[i,j])
                end
                @test eigmin(Hermitian(M))>=-1e-12
            end
        end
    end
end

@testset "Fermionic condition guards and empty row sets" begin
    @test_throws r"distinct conditions" Nx.FermionicConditions(1:3;conditions=(:D,))
    @test_throws r"distinct conditions" Nx.FermionicConditions(1:3;conditions=(:P,:P))
    @test_throws r"distinct conditions" Nx.FermionicConditions(1:3;conditions=())
    @test_throws r"sorted and distinct" Nx.FermionicConditions([2,1])
    _,(x,_,_)=Nx.pauli_variables(2)
    @test_throws r"requires fermionic" Nx.build(Nx.Problem(x[1]),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.FermionicConditions(1:2),)))
    _,(a,ad)=Nx.fermionic_variables(3)
    @test_throws r"outside the registry" Nx.build(Nx.Problem(zero(a[1])),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.FermionicConditions(1:4),)))
    H=ad[1]*a[1]*ad[3]*a[3]
    @test_throws r"do not cover the objective" Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.FermionicConditions(1:2),)))
    for (modes,conditions) in ((1:1,(:P,:Q,:T1,:T2)),(1:2,(:T1,)))
        built=Nx.build(Nx.Problem(zero(a[1])),Nx.BlockRelaxation(order=0,
            strengthen=(Nx.FermionicConditions(modes;conditions),)))
        @test built.sizes.closure.strengthening.added_blocks==0
    end
end

function fermi_permutation(images)
    n=length(images)
    U=zeros(ComplexF64,2^n,2^n)
    for state in 0:2^n-1
        occupied=[i for i in 1:n if !iszero(state & (1<<(n-i)))]
        mapped=images[occupied]
        target=sum((1<<(n-i) for i in mapped);init=0)
        inversions=sum((mapped[i]>mapped[j] for i in eachindex(mapped) for j in i+1:length(mapped));init=0)
        U[target+1,state+1]=isodd(inversions) ? -1 : 1
    end
    U
end

@testset "Fermionic modes close under symmetry, not moment-only maps" begin
    _,(a,ad)=Nx.fermionic_variables(4)
    numbers=[ad[i]*a[i] for i in 1:4]; N=sum(numbers)
    H=sum(-(ad[i]*a[mod1(i+1,4)]+ad[mod1(i+1,4)]*a[i])+0.3numbers[i]*numbers[mod1(i+1,4)] for i in 1:4)
    problem=Nx.Problem(H)
    translation=TestFermiMap([2,3,4,1],0,false)
    reflection=TestFermiMap([1,4,3,2],0,false)
    anti=TestFermiMap([1,2,3,4],0,true)
    conditions=(:P,:Q,:G,:T1,:T2,:T2prime)
    matrices=fermi_matrices(4)
    U=fermi_permutation(translation.permutation); V=fermi_permutation(reflection.permutation)
    for i in 1:4
        @test U*matrices.ladders[i]*U'≈matrices.ladders[translation.permutation[i]] atol=1e-12
    end
    for modes in (1:3,1:4),real_form in (false,true)
        generators=length(modes)==3 ? [translation] : [translation,reflection]
        moment_only=Tuple([length(modes)==3 ? [reflection] : TestFermiMap[];real_form ? [anti] : TestFermiMap[]])
        symmetry=Nx.FamilySymmetry(generators;moment_only)
        extra=(Nx.FermionicConditions(modes;conditions),Nx.Sector(N,2))
        built=Nx.build(problem,Nx.BlockRelaxation(order=1,grouping=Nx.Correlative(),split=Nx.TermSplit(),symmetry=symmetry,strengthen=extra))
        report=built.sizes.closure.strengthening
        @test report.initial_cones==6
        @test report.closed_cones==(length(modes)==3 ? 24 : 6)
        check_used_coordinates(built.program)
        # An explicit closed family with all symmetries applied only to moments
        # is the unreduced reference; no symmetry-adapted row projections.
        mode_sets=length(modes)==3 ? [[i for i in 1:4 if i!=missing] for missing in 1:4] : [collect(1:4)]
        reference=Nx.build(problem,Nx.BlockRelaxation(order=1,grouping=Nx.Correlative(),split=Nx.TermSplit(),
            symmetry=Nx.FamilySymmetry(;moment_only=(generators...,moment_only...)),
            strengthen=(map(m->Nx.FermionicConditions(m;conditions),mode_sets)...,Nx.Sector(N,2))))
        value=strengthening_bound(built.program)
        @test value≈strengthening_bound(reference.program) atol=1e-6
        sector=findall(v->abs(v-2)<1e-12,diag(matrices.polynomial(N)))
        @test value<=eigmin(Hermitian(matrices.polynomial(H)[sector,sector]))+1e-6
        for rank in (1,3)
            A=randn(MersenneTwister(rank),ComplexF64,length(sector),rank)
            rho=zeros(ComplexF64,16,16); rho[sector,sector]=A*A'
            rho=sum(U^k*(rho+V*rho*V')*(U')^k for k in 0:3)
            real_form && (rho=real.(rho))
            rho/=LinearAlgebra.tr(rho)
            values=fermi_coordinates(built.table,rho,matrices)
            @test maximum(abs,equality_values(built.program,values,1))<1e-12
            # Extra PSD blocks precede the sector equality block.
            last_psd=findlast(b->b.kind!=Nx.ZERO,built.program.blocks)
            p=deepcopy(built.program); resize!(p.blocks,last_psd)
            @test minimum(extra_values(p,values,report.first_block))>=-1e-12
        end
    end
end

@testset "Fermionic T2prime agrees with the frozen structured route" begin
    old=LEGACY_EXPECTED["fermionic-t2prime"]
    _,(a,ad)=Nx.fermionic_variables(4)
    H=sum(complex(term["real"],term["imag"])*foldl(*,
        (l<0 ? ad[-l] : a[l] for l in term["word"]);init=one(a[1])) for term in old["terms"])
    N=sum(ad[i]*a[i] for i in 1:4)
    built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.FermionicConditions(1:4;conditions=(:P,:Q,:G,:T1,:T2prime)),Nx.Sector(N,2))))
    value=strengthening_bound(built.program)
    @test value≈old["objective"] atol=1e-6
    matrices=fermi_matrices(4)
    sector=findall(v->abs(v-2)<1e-12,diag(matrices.polynomial(N)))
    @test value≈eigmin(Hermitian(matrices.polynomial(H)[sector,sector])) atol=1e-6
end

@testset "T1 exposes unused imaginary coordinates" begin
    _,(a,ad)=Nx.fermionic_variables(4)
    H=-(ad[1]*a[2]+ad[2]*a[1]+ad[3]*a[4]+ad[4]*a[3])
    description=Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,1,4,3])]),
        strengthen=(Nx.FermionicConditions(1:4;conditions=(:T1,)),))
    built=Nx.build(Nx.Problem(H),description)
    missing=check_used_coordinates(built.program;expected=2)
    @test built.sizes.closure.unused_coordinates==2
    @test all(built.table.labels[i].part==2 for i in missing)
    solver=Nx.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "decompose"=>false,"eps_abs"=>1e-9,"eps_rel"=>1e-9))
    result=Nx.solve(Nx.Problem(H),description;solver)
    @test result.status==MOI.OPTIMAL
    @test result.bound≈-2 atol=1e-7
    w=ad[3]*ad[1]*a[1]*a[4]
    coordinates=Nx.moment_terms(result.table,last(only(Nx.terms(w))))
    @test length(coordinates)==2
    @test count(pair->Nx.moment_covered(result.table,first(pair)),coordinates)==1
    before=(length(result.table.labels),length(result.table.cache),copy(result.table.covered))
    @test_throws r"Moment .*selected blocks" Nx.expectation(result,w)
    @test_throws r"Moment .*selected blocks" Nx.expectation(result,im*(w-w'))
    @test isfinite(Nx.expectation(result,w+w'))
    # The relaxation need not determine the physical ground state's quartics.
    # Check reconstruction separately on independent, feasible true moments.
    matrices=fermi_matrices(4)
    ground=eigen(Hermitian(matrices.polynomial(H))).vectors[:,1]
    values=fermi_coordinates(result.table,ground*ground',matrices)[2:end]
    solution=Nx.ConicSolution(-2.0,MOI.OPTIMAL,"true-state fixture",values,nothing)
    physical=Nx.SolveResult(-2.0,MOI.OPTIMAL,result.sizes,result.valid_for,result.problem,
        result.relaxation,solution,result.table,result.solver)
    @test Nx.expectation(physical,w+w')≈real(dot(ground,matrices.polynomial(w+w')*ground)) atol=1e-12
    @test Nx.expectation(result,H)≈result.bound atol=1e-7
    @test before==(length(result.table.labels),length(result.table.cache),result.table.covered)
end
