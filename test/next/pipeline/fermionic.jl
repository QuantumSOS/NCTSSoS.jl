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

@testset "Fermionic T2prime agrees with the legacy structured route" begin
    h=ComplexF64[-0.4 -0.7+0.2im;-0.7-0.2im 0.2]
    eri=zeros(ComplexF64,2,2,2,2); eri[1,1,1,1]=4; eri[2,2,2,2]=3.5
    linear=build_pqg_moment_data(Dict(0=>h),Dict((0,0,0,0)=>eri);
        nk=1,norb=2,nelec_per_cell=2,blocking=:none,spin_resolved_trace=false,
        singlet_s2=false,include_one_d=true,three_index=:T1T2prime)
    M=typeof(first(values(linear.key_to_monomial))); T=eltype(linear.identity)
    P=Polynomial{FermionicAlgebra,T,ComplexF64}
    mp=NCTSSoS._moment_problem_with_linear(FermionicAlgebra,T,M,P,zero(P),Tuple{Symbol,Matrix{P}}[],M[],0,linear)
    old=NCTSSoS.sos_dualize(mp)
    set_optimizer(old.model,SOLVER); optimize!(old.model)
    @test termination_status(old.model)==MOI.OPTIMAL
    _,(a,ad)=Nx.fermionic_variables(4)
    H=sum(c*foldl(*,(l<0 ? ad[-l] : a[l] for l in word);init=one(a[1])) for (word,c) in linear.objective_lin)
    N=sum(ad[i]*a[i] for i in 1:4)
    built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        strengthen=(Nx.FermionicConditions(1:4;conditions=(:P,:Q,:G,:T1,:T2prime)),Nx.Sector(N,2))))
    value=strengthening_bound(built.program)
    @test value≈objective_value(old.model) atol=1e-6
    matrices=fermi_matrices(4)
    sector=findall(v->abs(v-2)<1e-12,diag(matrices.polynomial(N)))
    @test value≈eigmin(Hermitian(matrices.polynomial(H)[sector,sector])) atol=1e-6
end

@testset "T1 exposes unused imaginary coordinates" begin
    _,(a,ad)=Nx.fermionic_variables(4)
    H=-(ad[1]*a[2]+ad[2]*a[1]+ad[3]*a[4]+ad[4]*a[3])
    built=Nx.build(Nx.Problem(H),Nx.BlockRelaxation(order=1,
        symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,1,4,3])]),
        strengthen=(Nx.FermionicConditions(1:4;conditions=(:T1,)),)))
    missing=check_used_coordinates(built.program;expected=2)
    @test built.sizes.closure.unused_coordinates==2
    @test all(built.table.labels[i].part==2 for i in missing)
end
