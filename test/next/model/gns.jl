# The public reconstruction contract. Solver-specific accuracy sweeps live in
# perf/next/gns.jl; these cases also accept exact moments without a solver.
using Test, LinearAlgebra

function gns_cases()
    N=NCTSSoS
    _,(a,b)=N.site_variables(N.Unipotent(),2,2)
    chsh=N.Problem(-a[1]*b[1]-a[1]*b[2]-a[2]*b[1]+a[2]*b[2])
    _,(v,)=N.site_variables(N.Free(),1,2)
    x,y=v
    free=N.Problem(x*y*x;constraints=(N.Positive(1-x^2),N.Positive(1-y^2)))
    _,(x,y,z)=N.pauli_variables(4)
    h=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+.5z[i]*z[mod1(i+1,4)] for i in 1:4)
    return ((:chsh,chsh,4),(:free,free,3),(:xxz,N.Problem(h),2))
end

function check_gns_backend(backend)
    N=NCTSSoS
    @testset "GNS through $(typeof(backend))" begin
        for (name,p,d) in gns_cases()
            result=N.solve(p,N.Relaxation(order=d);solver=backend)
            @test result.bound!==nothing
            g=name===:chsh ? check_gns_chsh(result) : name===:free ? check_gns_free(result) : check_gns_xxz(result)
            @test g.verdict!==:withheld
        end
    end
end

function check_gns_coverage(backend)
    N=NCTSSoS
    @testset "GNS coverage and state restrictions" begin
        registry,(x,y,z)=N.pauli_variables(4)
        h=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+z[i]*z[mod1(i+1,4)] for i in 1:4)
        p=N.Problem(h;symmetry=(N.SitePermutation([2,3,4,1]),N.spin_rotations(registry)...))
        result=N.solve(p,N.Relaxation(order=2);solver=backend)
        g=N.gns(result)
        @test g.verdict===:moments_reproduced
        @test N.expectation(g,h) ≈ -8 atol=1e-6
        @test g.checks.program_moment_error<=1e-6

        _,(x,y,z)=N.pauli_variables(2)
        grouped=N.solve(N.Problem(-z[1]-z[2]),N.Relaxation(order=1,grouping=N.Correlative());solver=backend)
        @test_throws r"sparse clique.*grouping=nothing" N.gns(grouped)
        custom=N.solve(N.Problem(-z[1]),N.Relaxation(rows=[one(z[1]),z[1]]);solver=backend)
        @test_throws r"GNS needs moment.*custom rows.*complete rows" N.gns(custom)

        _,(a,ad)=N.fermionic_variables(2)
        number=ad[1]*a[1]+ad[2]*a[2]
        h=-(ad[1]*a[2]+ad[2]*a[1])
        sector=N.solve(N.Problem(h),N.Relaxation(order=2,strengthen=N.Sector(number,1));solver=backend)
        sg=N.gns(sector)
        @test sg.verdict!==:withheld
        @test N.expectation(sg,number) ≈ 1 atol=1e-6
        @test maximum(sg.checks.sector_errors)<=1e-6
    end
end

function check_gns_chsh(result)
    N=NCTSSoS
    g=N.gns(result)
    @test g.status===:reproduced
    @test g.rank==g.full_rank==4
    variables=N.variables(result.problem.objective.registry)
    A=[g[variables[Symbol("s1_",i)]] for i in 1:2]
    B=[g[variables[Symbol("s2_",i)]] for i in 1:2]
    for a in (A...,B...)
        @test a*a ≈ I atol=1e-6
        @test a ≈ a' atol=1e-6
    end
    for a in A,b in B
        @test a*b ≈ b*a atol=1e-6
    end
    @test N.expectation(g,result.problem.objective) ≈ -2sqrt(2) atol=1e-6
    @test g.checks.hankel_moment_error<=g.flatness.check_atol
    @test g.checks.program_moment_error<=g.flatness.check_atol
    @test occursin("moments reproduced",sprint(show,g))
    return g
end

function check_gns_xxz(result)
    N=NCTSSoS
    g=N.gns(result)
    @test g.status===:reproduced
    @test (g.full_degree,g.principal_degree,g.relaxation_order)==(3,2,2)
    @test g.rank==g.full_rank==16
    @test occursin("above relaxation order 2",sprint(show,g))
    v=N.variables(result.problem.objective.registry)
    for i in 1:4
        X,Y,Z=(g[v[Symbol(axis,i)]] for axis in (:X,:Y,:Z))
        @test X*Y ≈ im*Z atol=1e-6
        @test Y*Z ≈ im*X atol=1e-6
        @test Z*X ≈ im*Y atol=1e-6
        @test X*X ≈ I atol=1e-6
        @test Y*Y ≈ I atol=1e-6
        @test Z*Z ≈ I atol=1e-6
    end
    for i in 1:4,j in 1:i-1,a in (:X,:Y,:Z),b in (:X,:Y,:Z)
        A,B=g[v[Symbol(a,i)]],g[v[Symbol(b,j)]]
        @test A*B ≈ B*A atol=1e-6
    end
    @test N.expectation(g,result.problem.objective) ≈ result.bound atol=1e-6
    return g
end

function check_gns_free(result)
    N=NCTSSoS
    g=N.gns(result)
    @test g.verdict in (:moments_reproduced,:optimizer_certified)
    @test N.expectation(g,result.problem.objective) ≈ result.bound atol=1e-6
    for c in result.problem.constraints
        @test eigmin(Hermitian(N.evaluate(g,c.polynomial)))>=-1e-6
    end
    return g
end

# Populate a plain solved result from independently supplied matrix moments.
# The numerical algorithm sees the same public result that a solver returns.
function gns_fixture(problem,order,moment)
    N=NCTSSoS
    relaxation=N.Relaxation(order=order)
    built=N.build(problem,relaxation)
    values=[label.part==0x01 ? real(moment(label.word)) : imag(moment(label.word)) for label in built.table.labels[2:end]]
    bound=real(sum(c*moment(w) for (c,w) in N.terms(problem.objective);init=0.0im))
    solution=N.ConicSolution(bound,:OPTIMAL,"exact matrix fixture",values,nothing)
    return N.SolveResult(bound,:OPTIMAL,built.sizes,built.valid_for,problem,relaxation,solution,built.table,nothing)
end

function gns_matrix_word(w,matrices)
    dim=size(first(values(matrices)),1)
    return foldl((A,i)->A*matrices[i],NCTSSoS.legacy_letters(w);init=Matrix{ComplexF64}(I,dim,dim))
end

@testset "GNS from exact matrix moments" begin
    N=NCTSSoS
    id=Matrix{ComplexF64}(I,2,2)
    X=ComplexF64[0 1;1 0]; Y=ComplexF64[0 -im;im 0]; Z=ComplexF64[1 0;0 -1]
    _,(x,y,z)=N.pauli_variables(4)
    physical=Dict(3(i-1)+a=>reduce(kron,(j==i ? (X,Y,Z)[a] : id for j in 1:4)) for i in 1:4 for a in 1:3)
    h=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+.5z[i]*z[mod1(i+1,4)] for i in 1:4)
    H=sum(physical[3(i-1)+a]*physical[3(mod1(i+1,4)-1)+a]*(a==3 ? .5 : 1) for i in 1:4 for a in 1:3)
    psi=eigen(Hermitian(H)).vectors[:,1]
    result=gns_fixture(N.Problem(h),2,w->dot(psi,gns_matrix_word(w,physical)*psi))
    g=check_gns_xxz(result)
    @test g.flatness.minimum_eigenvalue>=-g.flatness.psd_tolerance
    @test !occursin("; reproduced;",sprint(show,g))
    @test occursin("max error $(round(g.checks.program_moment_error;sigdigits=3));",sprint(show,g))
    # Enumerate the tensor-product Pauli basis without a package word builder.
    words=[Int[3(i-1)+a for (i,a) in enumerate(axes) if a!=0]
        for axes in Iterators.product(ntuple(_->0:3,4)...)]
    represented(word)=foldl((A,i)->A*physical[i],word;init=Matrix{ComplexF64}(I,16,16))
    principal=hcat([represented(word)*psi for word in words if length(word)<=2]...)
    full=hcat([represented(word)*psi for word in words if length(word)<=3]...)
    @test (rank(principal;atol=1e-6),rank(full;atol=1e-6))==(g.rank,g.full_rank)
    @test svdvals(principal'*principal) ≈ g.singular_values atol=1e-6
    @test svdvals(full'*full) ≈ g.full_singular_values atol=1e-6
    variables=[v for i in 1:4 for v in (x[i],y[i],z[i])]
    for word in words
        p=foldl(*,(variables[i] for i in word);init=one(x[1]))
        @test N.expectation(g,p) ≈ dot(psi,represented(word)*psi) atol=1e-6
    end
    _,(other,_,_)=N.pauli_variables(4)
    @test_throws r"different registry" g[other[1]]
    @test_throws r"different registry" N.evaluate(g,other[1])
    @test_throws r"single variable" g[x[1]+x[2]]
    @test_throws r"nonnegative" N.gns(result;atol=-1)

    _,(a,b)=N.site_variables(N.Unipotent(),2,2)
    bell=N.Problem(-a[1]*b[1]-a[1]*b[2]-a[2]*b[1]+a[2]*b[2])
    matrices=Dict(257=>kron(Z,id),513=>kron(X,id),258=>kron(id,(Z+X)/sqrt(2)),514=>kron(id,(Z-X)/sqrt(2)))
    psi=ComplexF64[1,0,0,1]/sqrt(2)
    bell_result=gns_fixture(bell,4,w->dot(psi,gns_matrix_word(w,matrices)*psi))
    check_gns_chsh(bell_result)
    # A represented moment outside the selected degree-three Hankel still
    # matters. The same matrices remain an optimizer, not a reconstruction.
    changed=deepcopy(bell_result)
    slot=findfirst(l->N.degree(l.word)>6,changed.table.labels)
    @test !isnothing(slot)
    changed.solution.values[slot-1]+=0.125
    optimizer=N.gns(changed)
    @test optimizer.verdict===:optimizer_certified
    @test optimizer.checks.hankel_moment_error<1e-6
    @test optimizer.checks.program_moment_error>0.1
    @test occursin("moments NOT reproduced",sprint(show,optimizer))
    @test N.expectation(optimizer,changed.problem.objective) ≈ changed.bound atol=1e-6
    # A reserved coordinate is not a solved moment. The public coverage hook
    # must exclude its arbitrary value, while checking the covered component.
    coordinates=[[Int32(i)=>big(1)//big(1)] for i in eachindex(changed.table.labels)]
    covered=BitSet(eachindex(coordinates)); delete!(covered,slot)
    partial=N.ContinuousMoments(changed.table,changed.table.labels,coordinates,(),covered)
    partial_result=N.SolveResult(changed.bound,changed.status,changed.sizes,changed.valid_for,
        changed.problem,changed.relaxation,changed.solution,partial,changed.solver)
    partial_gns=N.gns(partial_result)
    @test partial_gns.verdict===:moments_reproduced
    @test partial_gns.checks.program_moment_error<1e-6
    @test N.expectation(partial_gns,changed.problem.objective) ≈ changed.bound atol=1e-6
    # A missing low-degree component cannot be filled from arbitrary solver
    # coordinates just because a flat extension would otherwise exist.
    low=findfirst(l->N.degree(l.word)==1,changed.table.labels)
    delete!(covered,low)
    @test_throws r"GNS needs moment.*complete rows" N.gns(partial_result)
    # Solver OPTIMAL alone is insufficient for the optimizer verdict.
    wrong_bound=N.SolveResult(changed.bound-0.132345,changed.status,changed.sizes,changed.valid_for,
        changed.problem,changed.relaxation,changed.solution,changed.table,changed.solver)
    withheld=N.gns(wrong_bound)
    @test withheld.verdict===:withheld
    @test occursin("check_atol",sprint(show,withheld))
    # Rounded advice must not recommend a tolerance below the measured error.
    recommended=parse(Float64,match(r"check_atol ≥ ([0-9.e+-]+)",withheld.advice)[1])
    @test recommended>=0.132345
    @test N.gns(bell_result;atol=1e6).verdict===:withheld
    # Even exact flat moments and an OPTIMAL status cannot certify matrices
    # that violate a problem constraint, including an operator equality.
    for constraint in (N.Positive(-one(a[1])),N.Zero(one(a[1])))
        problem=N.Problem(bell.objective;constraints=(constraint,))
        invalid=N.SolveResult(bell_result.bound,bell_result.status,bell_result.sizes,bell_result.valid_for,
            problem,bell_result.relaxation,bell_result.solution,bell_result.table,bell_result.solver)
        @test N.gns(invalid).verdict===:withheld
    end

    # A trace-one functional can be positive on the order-one rows but fail
    # positivity on the higher Hankel matrix built from those same moments.
    _,(x,y,z)=N.pauli_variables(2)
    physical=Dict(3(i-1)+a=>kron(i==1 ? (X,Y,Z)[a] : id,i==2 ? (X,Y,Z)[a] : id) for i in 1:2 for a in 1:3)
    rho=Diagonal([-0.01,1.01/3,1.01/3,1.01/3])
    bad=gns_fixture(N.Problem(zero(x[1])),1,w->LinearAlgebra.tr(rho*gns_matrix_word(w,physical)))
    matrices=[Matrix{ComplexF64}(I,4,4),[physical[i] for i in 1:6]...]
    gram=[LinearAlgebra.tr(rho*A'*B) for A in matrices,B in matrices]
    @test eigmin(Hermitian(gram))>0
    rejected=N.gns(bad)
    @test rejected.full_degree>rejected.relaxation_order
    @test rejected.status===:not_positive
    @test rejected.flatness.minimum_eigenvalue < -rejected.flatness.psd_tolerance
    @test rejected.operators===nothing && rejected.state===nothing
    @test occursin("higher relaxation order",sprint(show,rejected))
    @test !occursin("psd_rtol",sprint(show,rejected))
    @test_throws r"withheld" rejected[x[1]]
    @test_throws r"withheld" N.expectation(rejected,one(x[1]))

    # Two-mode CAR, including complex moments and creation adjoints.
    _,(a,ad)=N.fermionic_variables(2)
    annih=ComplexF64[0 1;0 0]
    f=Dict(1=>kron(annih,id),2=>kron(Z,annih),-1=>kron(annih',id),-2=>kron(Z,annih'))
    psi=ComplexF64[1,2im,3,4im]; psi/=norm(psi)
    fermion=gns_fixture(N.Problem(ad[1]*a[1]),2,w->dot(psi,gns_matrix_word(w,f)*psi))
    fg=N.gns(fermion)
    @test fg.status===:reproduced
    for i in 1:2,j in 1:2
        @test fg[a[i]]*fg[a[j]]+fg[a[j]]*fg[a[i]] ≈ zeros(fg.rank,fg.rank) atol=1e-6
        @test fg[a[i]]*fg[ad[j]]+fg[ad[j]]*fg[a[i]] ≈ (i==j ? Matrix{ComplexF64}(I,fg.rank,fg.rank) : zeros(fg.rank,fg.rank)) atol=1e-6
        @test fg[ad[i]] ≈ fg[a[i]]' atol=1e-6
    end
end
