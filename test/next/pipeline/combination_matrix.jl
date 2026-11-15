module CombinationMatrix

import NCTSSoS
using LinearAlgebra, Random, Test
const Nx=NCTSSoS.Next

const PROBLEMS=(:heisenberg,:xxz,:ising,:fermions,:unipotent,:positive)
const FINITE=(:none,:translation,:dihedral,:moment_only,:antiunitary,:clifford,:detected,:translation_antiunitary)
const CONTINUOUS=(:none,:u1,:su2,:number)
const FAMILIES=[(order=o,rows=r,group=g,split=s) for o in 1:2
    for r in (:complete,:contiguous) for g in (:none,:correlative,:physics)
    for s in (:none,:term,:physics)]
append!(FAMILIES,[(order=o,rows=:custom,group=:none,split=:none) for o in 1:2])

function strengthening_choices(problem)
    kinds=problem==:unipotent ? Symbol[] : problem==:fermions ? [:sector,:fermionic] :
        [:density,:invariant,:stationarity,:optimality,:sector]
    [();[(k,) for k in kinds];[(kinds[i],kinds[j]) for i in eachindex(kinds) for j in i+1:length(kinds)]]
end

six_site_controls()=[(;n=6,problem,family=first(FAMILIES),finite,continuous,strengthening=())
    for (problem,finite,continuous) in ((:heisenberg,:dihedral,:su2),(:fermions,:translation,:number))]

detected_continuous_control()=(;n=5,problem=:xxz,family=first(FAMILIES),finite=:detected,
    continuous=:u1,strengthening=())

bosonic_cases()=[(;n=2,problem=:bosons,family,finite,continuous,strengthening)
    for family in FAMILIES for finite in (:none,:translation,:dihedral,:moment_only)
    for continuous in (:none,:number) for strengthening in ((),(:sector,))]

combinations()=[[(;n=4,problem,family,finite,continuous,strengthening) for problem in PROBLEMS
    for family in FAMILIES for finite in FINITE for continuous in CONTINUOUS
    for strengthening in strengthening_choices(problem)];six_site_controls();detected_continuous_control();bosonic_cases()]

# Total-particle cutoff, not a per-mode truncation. N+6 covers every word
# in these order-two families, including the sector's extra degree-two factor.
function boson_matrices(n,cutoff)
    occupations=filter(v->sum(v)<=cutoff,vec(collect(Iterators.product(ntuple(_->0:cutoff,n)...))))
    index=Dict(v=>i for (i,v) in enumerate(occupations))
    [begin
        a=zeros(ComplexF64,length(occupations),length(occupations))
        for (j,v) in enumerate(occupations)
            v[mode]==0 && continue
            u=ntuple(i->v[i]-(i==mode),n)
            a[index[u],j]=sqrt(v[mode])
        end
        a
    end for mode in 1:n]
end

function fixture(name,n=4;anisotropy=0.4)
    I2=Matrix{ComplexF64}(I,2,2)
    X=ComplexF64[0 1;1 0]; Y=ComplexF64[0 -im;im 0]; Z=ComplexF64[1 0;0 -1]
    embed(A,i)=foldl(kron,(j==i ? A : I2 for j in 1:n))
    matrices=Dict{Int,Matrix{ComplexF64}}()
    if name==:bosons
        registry,(a,ad)=Nx.bosonic_variables(n)
        edges=n==2 ? [(1,2)] : [(i,mod1(i+1,n)) for i in 1:n]
        H=-sum(ad[i]*a[j]+ad[j]*a[i] for (i,j) in edges)+
            sum(ad[i]^2*a[i]^2-0.5ad[i]*a[i] for i in 1:n)
        Q=sum(ad.*a); q=n; operators=[a;ad]
        for (i,A) in enumerate(boson_matrices(n,n+6))
            matrices[i]=A; matrices[-i]=A'
        end
    elseif name==:fermions
        registry,(a,ad)=Nx.fermionic_variables(n)
        H=sum((1+0.3im)*ad[i]*a[mod1(i+1,n)]+(1-0.3im)*ad[mod1(i+1,n)]*a[i] for i in 1:n)
        Q=sum(ad[i]*a[i] for i in 1:n); q=n÷2
        operators=[a;ad]
        for i in 1:n
            A=foldl(kron,(j<i ? Z : j==i ? ComplexF64[0 1;0 0] : I2 for j in 1:n))
            matrices[i]=A; matrices[-i]=A'
        end
    elseif name==:unipotent
        registry,ops=Nx.site_variables(Nx.Unipotent(),n,2)
        a=[ops[i][1] for i in 1:n]; b=[ops[i][2] for i in 1:n]
        H=-sum(a[i]*a[mod1(i+1,n)]+0.7b[i] for i in 1:n)
        Q=zero(H); q=0; operators=[a;b]
        for (ps,Ms) in ((a,X),(b,Z)),i in 1:n
            matrices[only(Nx.legacy_letters(first(ps[i].data.words)))]=embed(Ms,i)
        end
    else
        registry,(x,y,z)=Nx.pauli_variables(n)
        H=name==:ising ? -sum(z[i]*z[mod1(i+1,n)]+0.7x[i] for i in 1:n) :
            sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+(name==:xxz ? anisotropy : 1)*z[i]*z[mod1(i+1,n)] for i in 1:n)
        Q=sum(z); q=0; operators=[x;y;z]
        for (a,M) in enumerate((X,Y,Z)),i in 1:n
            matrices[3(i-1)+a]=embed(M,i)
        end
    end
    dimension=size(first(values(matrices)),1)
    identity=Matrix{ComplexF64}(I,dimension,dimension)
    cache=Dict{Any,Matrix{ComplexF64}}()
    word(w)=get!(cache,w) do
        @assert name!=:bosons || Nx.degree(w)<=6 "Increase the Fock cutoff for this word degree"
        foldl(*,(matrices[Int(l)] for l in Nx.legacy_letters(w));init=identity)
    end
    polynomial(p)=sum((c*word(w) for (c,w) in Nx.terms(p));init=zeros(ComplexF64,size(identity)))
    letters=Dict(only(Nx.legacy_letters(first(p.data.words)))=>p for p in operators)
    row(w)=foldl(*,(letters[l] for l in Nx.legacy_letters(w));init=one(H))
    constraints=name==:positive ? (Nx.Positive(3n+H),) : ()
    # Full ground-space projector, not one arbitrarily chosen eigenvector.
    states=Dict{Bool,Matrix{ComplexF64}}()
    function state(sector)
        get!(states,sector) do
            U=if sector || name==:bosons
                E=eigen(Hermitian(polynomial(Q)))
                E.vectors[:,findall(v->abs(v-q)<1e-9,E.values)]
            else
                identity
            end
            E=eigen(Hermitian(U'*polynomial(H)*U))
            V=U*E.vectors[:,findall(v->abs(v-first(E.values))<1e-9,E.values)]
            V*V'/size(V,2)
        end
    end
    (;name,n,registry,H,Q,q,operators,row,word,polynomial,constraints,state,
        declarations=Dict(),moments=Dict{Tuple{Bool,Any},ComplexF64}())
end

function declaration(f,finite,continuous)
    get!(f.declarations,(finite,continuous)) do
        make_declaration(f,finite,continuous)
    end
end
function make_declaration(f,finite,continuous)
    ds=continuous==:none ? () : continuous==:u1 ? (Nx.SpinRotation(f.registry,:z),) :
        continuous==:su2 ? Nx.spin_rotations(f.registry) : (Nx.NumberRotation(f.registry),)
    t=Nx.SitePermutation([2:f.n;1]); r=Nx.SitePermutation([1;f.n:-1:2])
    if finite==:detected
        detected=Nx.detect_symmetry(Nx.Problem(f.H;constraints=f.constraints);max_automorphisms=16)
        return Nx.Symmetry(detected.generators...;moment_only=detected.moment_only,continuous=ds)
    end
    # Wrong-algebra declarations pass through the public domain guard too.
    _,(x,_,_)=Nx.pauli_variables(f.n)
    anti=Nx.AntiunitaryPauli(first(one(x[1]).data.words))
    flip=Nx.CliffordMap(Nx.PauliConjugation(first(prod(x).data.words),false))
    gs=finite in (:translation,:translation_antiunitary) ? (t,) : finite==:dihedral ? (t,r) : finite==:clifford ? (flip,) : ()
    ms=finite==:moment_only ? (t,) : finite in (:antiunitary,:translation_antiunitary) ? (anti,) : ()
    Nx.Symmetry(gs...;moment_only=ms,continuous=ds)
end

function strengthening(f,names)
    Tuple(k==:density ? Nx.DensityMatrix(1:f.n) : k==:invariant ? Nx.InvariantCone(1:f.n,"Y"*"I"^(f.n-1)) :
        k==:stationarity ? Nx.Stationarity() : k==:optimality ? Nx.StateOptimality() :
        k==:sector ? Nx.Sector(f.Q,f.q) : Nx.FermionicConditions(1:f.n;conditions=(:P,:Q,:G,:T1,:T2,:T2prime)) for k in names)
end

function relaxation(f,choice,kinds)
    if choice.rows==:custom
        W=Nx.word_type(f.registry,2choice.order)
        rows=f.row.(Nx.basis_rows(W,f.registry.algebra,1:f.n,choice.order))
        return Nx.Relaxation(family=[rows];strengthen=kinds)
    end
    grouping=choice.group==:none ? nothing : choice.group==:correlative ? Nx.Correlative() : Nx.PhysicsInformed(level=0)
    split=choice.split==:none ? nothing : choice.split==:term ? Nx.TermSplit() : Nx.PhysicsInformed(level=1)
    row_rule=choice.rows==:contiguous ? Nx.Contiguous(periodic=true) : nothing
    Nx.Relaxation(order=choice.order;grouping,split,row_rule,strengthen=kinds)
end

function reference(f,problem,description,built)
    if !hasproperty(built.sizes,:closure)
        # Dense no-symmetry front-door route is already the reference.
        return built
    end
    closed=built.sizes.closure
    W=eltype(first(closed.family).rows)
    objective=[(ComplexF64(c),Nx.word_from_legacy(W,Nx.legacy_letters(w))) for (c,w) in Nx.terms(f.H)]
    p,t,c=Nx.family_program(closed.family,objective;real_form=Nx._real_form(problem),check_coverage=false)
    kinds=description.strengthen
    initial=[Nx._pauli_extra(k,W,length(problem.constraints)+i) for (i,k) in enumerate(kinds)
        if k isa Union{Nx.DensityMatrix,Nx.InvariantCone}]
    if !isempty(initial)
        extras,_=Nx._extra_orbits(initial,closed.block_generators)
        for extra in extras; Nx._append_extra!(p,t,extra); end
    end
    Nx._hamiltonian_strengthen!(p,t,problem,kinds,c)
    Nx._fermionic_strengthen!(p,t,problem,kinds,c)
    sectors=Nx._checked_sectors(problem,kinds,Nx.Symmetry(),W)
    Nx._sector_equalities!(p,t,c,sectors)
    p.nvars=length(t.labels)
    Nx.validate(p)
    (;program=p,table=t)
end

function coordinate_map(target,source)
    # A reference stationarity equation can contain a whole orbit that cancels
    # after symmetry, so the reduced program need not register that orbit.
    # Keep it free in the comparison space, rather than assigning zero or
    # dropping reference equations. A lost constraint will still be detected.
    if target isa Nx.OrbitMoments && any(l->isnothing(Nx.moment_terms(target,l.word)),source.labels)
        target=deepcopy(target)
        for label in source.labels
            Nx.orbit_coordinates!(target,label.word)
        end
    end
    M=zeros(Float64,length(source.labels),length(target.labels))
    for (i,label) in enumerate(source.labels)
        pairs=Nx.moment_terms(target,label.word)
        isnothing(pairs) && error("Reference word is absent from reduced moment domain: $(label.word)")
        for (j,c) in pairs
            M[i,j]+=label.part==2 ? imag(c) : real(c)
        end
    end
    M,target
end

function cone_values(program,y)
    spectra=Float64[]; soc=Float64[]; equalities=ComplexF64[]
    value(f)=sum((program.forms.coefs[t]*y[program.forms.ids[t]] for t in Nx.term_range(program.forms,f));init=0im)
    for b in program.blocks
        if b.kind==Nx.ZERO
            append!(equalities,value.(Nx.form_range(b)))
        elseif b.kind==Nx.SOC
            v=value.(Nx.form_range(b))
            append!(soc,[real(v[1])-norm(v[2:end]),real(v[1])+norm(v[2:end])])
        else
            A=zeros(ComplexF64,b.side,b.side)
            for j in 1:b.side,i in j:b.side
                A[i,j]=value(b.first_form+Nx.tri_offset(b.side,i,j)); A[j,i]=conj(A[i,j])
            end
            append!(spectra,eigvals(Hermitian(A)))
        end
    end
    (;spectra=sort!(spectra),soc=sort!(soc),equalities)
end

function distinct(values;tol=1e-7)
    out=Float64[]
    for v in values
        (isempty(out) || abs(v-last(out))>tol*max(1,abs(v))) && push!(out,v)
    end
    out
end

function equality_basis(program,M)
    rows=Vector{Float64}[]
    for b in program.blocks
        b.kind==Nx.ZERO || continue
        for f in Nx.form_range(b)
            row=zeros(ComplexF64,size(M,2))
            for t in Nx.term_range(program.forms,f)
                row .+= program.forms.coefs[t].*view(M,program.forms.ids[t],:)
            end
            for part in (real,imag)
                r=part.(row); norm(r)>1e-10 && push!(rows,r/norm(r))
            end
        end
    end
    isempty(rows) && return zeros(0,size(M,2))
    S=svd(reduce(hcat,rows))
    S.U[:,findall(v->v>1e-8*first(S.S),S.S)]'
end

table_size(t)=(length(t.labels),hasproperty(t,:cache) ? length(t.cache) : 0,
    hasproperty(t,:base) ? table_size(t.base) : nothing)

function check_program(f,p,r,built;trials=2,solver=nothing)
    ref=reference(f,p,r,built)
    M,sampling=coordinate_map(built.table,ref.table)
    a=equality_basis(ref.program,M)
    b=equality_basis(built.program,Matrix{Float64}(I,built.program.nvars,length(sampling.labels)))
    equality_error=maximum(abs,a'*a-b'*b;init=0.0)
    @assert isapprox(a'*a,b'*b;atol=1e-7,rtol=1e-7) "Equality subspaces differ"
    rng=MersenneTwister(731)
    worst=0.0
    for _ in 1:trials
        y=randn(rng,length(sampling.labels)); y[1]=1
        va=cone_values(ref.program,M*y); vb=cone_values(built.program,y)
        for key in (:spectra,:soc)
            x=distinct(getproperty(va,key)); z=distinct(getproperty(vb,key))
            @assert length(x)==length(z) "Distinct $key counts differ: $(length(x)) / $(length(z))"
            err=maximum(abs.(x-z);init=0.0); worst=max(worst,err)
            @assert isapprox(x,z;atol=1e-7,rtol=1e-7) "Cone spectra differ: $err"
        end
        equality_error=max(equality_error,abs(norm(a*y)-norm(b*y)))
        @assert isapprox(norm(a*y),norm(b*y);atol=1e-7,rtol=1e-7) "Equality residuals differ"
    end
    sector=any(k->k isa Nx.Sector,r.strengthen)
    rho=f.state(sector)
    value(w)=get!(f.moments,(sector,w)) do
        tr(rho*f.word(w))
    end
    y=[label.part==2 ? imag(value(label.word)) : real(value(label.word)) for label in sampling.labels]
    @assert isapprox(M*y,[label.part==2 ? imag(value(label.word)) : real(value(label.word)) for label in ref.table.labels];atol=1e-8,rtol=1e-8) "True state is not invariant under coordinate reduction"
    for g in (p.symmetry.generators...,p.symmetry.moment_only...),label in ref.table.labels
        phase,w=Nx.act(g,label.word)
        z=Nx.scale(phase,value(w))
        @assert isapprox(value(label.word),Nx.is_antiunitary(g) ? conj(z) : z;atol=1e-8,rtol=1e-8) "True state violates a declared finite symmetry"
    end
    for d in p.symmetry.continuous,label in ref.table.labels
        derivative=Nx._derivative(d,label.word)
        @assert abs(sum((c*value(w) for (w,c) in derivative);init=0im))<1e-8 "True state violates a declared derivation"
    end
    v=cone_values(built.program,y)
    @assert minimum([v.spectra;v.soc];init=0.0)>=-1e-8 "True state violates a cone"
    @assert maximum(abs,v.equalities;init=0.0)<1e-8 "True state violates an equality"
    queried=blocked=0
    bound_error=NaN
    solve_status="not_selected"
    if !isnothing(solver)
        result=Nx.solve(p,r;solver)
        baseline=Nx.solve_standard(Nx.moment_standard_form(ref.program),solver)
        solve_status="$(result.status)/$(baseline.status)"
        if isnothing(result.bound) || isnothing(baseline.bound)
            return (;spectral_error=worst,equality_error,state=:pass,bound_error,
                solve_status,queries="not_checked")
        end
        bound_error=abs(result.bound-baseline.bound)
        @assert bound_error<2e-6 "Reduced and unreduced bounds differ: $bound_error"
        @assert isapprox(Nx.expectation(result,f.H),result.bound;atol=1e-7,rtol=1e-7)
        before=table_size(result.table)
        queries=[f.row(w) for w in unique(label.word for label in ref.table.labels)]
        # A degree-three query is outside many order-one families; a symmetry
        # may nevertheless force it to zero. Neither answer may add coordinates.
        push!(queries,prod(f.operators[1:3]))
        W=typeof(first(ref.table.labels).word)
        for query in queries
            w=Nx._pack_word(W,first(query.data.words))
            pairs=Nx.moment_terms(result.table,w)
            coefficients=Dict{Int,ComplexF64}()
            for (id,c) in something(pairs,()); coefficients[id]=get(coefficients,id,0im)+c; end
            covered=!isnothing(pairs) && all(iszero(c) || Nx.moment_covered(result.table,id) for (id,c) in coefficients)
            if covered
                Nx.expectation(result,query); queried+=1
            else
                failure=try
                    Nx.expectation(result,query); nothing
                catch e
                    e
                end
                @assert failure isa ArgumentError "An uncovered query did not raise ArgumentError"
                blocked+=1
            end
        end
        @assert before==table_size(result.table) "An expectation query changed the moment domain"
    end
    (;spectral_error=worst,equality_error,state=:pass,bound_error,solve_status,
        queries=isnothing(solver) ? "not_selected" : "$queried readable; $blocked blocked")
end

function run(c;f=fixture(c.problem,c.n),solver=nothing)
    built=nothing
    try
        symmetry=declaration(f,c.finite,c.continuous)
        p=Nx.Problem(f.H;constraints=f.constraints,symmetry)
        r=relaxation(f,c.family,strengthening(f,c.strengthening))
        built=Nx.build(p,r)
        checks=check_program(f,p,r,built;solver)
        return (;c...,outcome="program",variables=built.program.nvars-1,
            blocks=length(built.program.blocks),largest=maximum(b.side for b in built.program.blocks),checks...,message="")
    catch e
        named=isnothing(built) && e isa ArgumentError &&
            occursin(r"(?i)remove|use |supply|select|increase|add |construct|choose|put |requires? .*variables|supports? .*only",sprint(showerror,e))
        return (;c...,outcome=named ? "ArgumentError" : "DEFECT",variables=isnothing(built) ? 0 : built.program.nvars-1,
            blocks=isnothing(built) ? 0 : length(built.program.blocks),largest=0,spectral_error=NaN,state=:not_checked,
            equality_error=NaN,bound_error=NaN,solve_status="not_checked",queries="not_checked",message=sprint(showerror,e))
    end
end

function solved_cases()
    selected=NamedTuple[(;n=4,problem,family=first(FAMILIES),finite=:translation,continuous=:none,
        strengthening=()) for problem in PROBLEMS]
    for (problem,order,rows,group,split,finite,continuous,strengthening) in (
        (:heisenberg,1,:complete,:none,:none,:dihedral,:su2,()),
        (:fermions,1,:complete,:none,:none,:translation,:number,()),
        (:heisenberg,1,:complete,:correlative,:term,:moment_only,:none,()),
        (:xxz,1,:complete,:physics,:term,:clifford,:u1,(:density,:stationarity)),
        (:ising,1,:complete,:correlative,:term,:antiunitary,:none,(:density,:sector)),
        (:heisenberg,1,:complete,:physics,:physics,:detected,:none,(:density,:optimality)),
        (:fermions,1,:complete,:correlative,:term,:translation,:number,(:sector,:fermionic)),
        (:positive,2,:custom,:none,:none,:translation,:su2,()))
        push!(selected,(;n=4,problem,family=(;order,rows,group,split),finite,continuous,strengthening))
    end
    append!(selected,six_site_controls())
    append!(selected,[(;n=4,problem=:xxz,family=(order=o,rows=:complete,group=:none,split=:none),
        finite=:translation_antiunitary,continuous=:none,strengthening=()) for o in 1:2])
    push!(selected,detected_continuous_control())
    selected
end

function smoke_cases()
    all=filter(c->c.problem!=:bosons && c.finite!=:translation_antiunitary && c!=detected_continuous_control(),combinations())
    # Fixed, stratified indices across the complete product, plus each named
    # axis in a valid dense control. No RNG chooses the acceptance subset.
    selected=all[unique(round.(Int,range(1,length(all);length=80)))]
    for problem in PROBLEMS,finite in filter(!=(:translation_antiunitary),FINITE)
        push!(selected,(;n=4,problem,family=first(FAMILIES),finite,continuous=:none,strengthening=()))
    end
    for problem in (:heisenberg,:fermions),continuous in CONTINUOUS,order in 1:2
        push!(selected,(;n=4,problem,family=(order=order,rows=:complete,group=:none,split=:none),
            finite=:translation,continuous,strengthening=()))
    end
    append!(selected,six_site_controls())
    push!(selected,(;n=4,problem=:heisenberg,family=first(FAMILIES),finite=:dihedral,continuous=:su2,strengthening=()))
    push!(selected,(;n=4,problem=:xxz,family=(order=2,rows=:contiguous,group=:none,split=:physics),
        finite=:detected,continuous=:none,strengthening=(:stationarity,)))
    append!(selected,[(;n=4,problem,family=first(FAMILIES),finite=:translation_antiunitary,
        continuous=:none,strengthening=()) for problem in PROBLEMS])
    push!(selected,detected_continuous_control())
    append!(selected,filter(c->c.family==(order=2,rows=:complete,group=:none,split=:none),bosonic_cases()))
    unique(selected)
end

function write_table(io,results;header=true)
    header && println(io,"problem\tn\torder\trows\tgroup\tsplit\tfinite\tcontinuous\tstrengthening\toutcome\tvariables\tblocks\tlargest\tspectral_error\tequality_error\tstate\tbound_error\tsolve_status\tqueries\tmessage")
    for r in results
        cells=(r.problem,r.n,r.family.order,r.family.rows,r.family.group,r.family.split,r.finite,r.continuous,
            join(r.strengthening,"+"),r.outcome,r.variables,r.blocks,r.largest,r.spectral_error,r.equality_error,r.state,r.bound_error,r.solve_status,r.queries,r.message)
        println(io,join((replace(string(c),'\t'=>' ','\n'=>' ') for c in cells),'\t'))
    end
end

function smoke_test()
    # Recorded build/error outcomes, not merely "no uncaught exception": a new
    # rejection of a previously supported combination must fail this test.
    builds=Set([1,2,3,4,5,6,7,9,11,15,17,19,21,23,26,28,30,32,34,49,51,53,55,59,60,61,
        67,69,71,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,
        104,108,109,110,111,115,116,117,118,119,120,122,123,124,125,126,129,134,135,136,137,138,
        139,140,141,144,145])
    cases=smoke_cases()
    fixtures=Dict((c.problem,c.n)=>fixture(c.problem,c.n) for c in unique(c->(c.problem,c.n),cases))
    @testset "Combination product: fixed acceptance slice" begin
        @test length(combinations())==84_515
        @test length(cases)==161
        for (i,c) in enumerate(cases)
            @testset "$i: $c" begin
                r=run(c;f=fixtures[(c.problem,c.n)])
                expected=i in builds || c.problem==:bosons ? "program" : "ArgumentError"
                r.outcome==expected || @info "Combination mismatch" case=c expected result=r
                @test r.outcome==expected
            end
        end
    end
end

end # module
