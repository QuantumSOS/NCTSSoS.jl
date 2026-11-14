# Independent Hilbert-space expansion. This reference never calls the
# strengthening products, row projections, or commutator implementation.
function physical_form!(program,table,A,words,matrices)
    for (w,P) in zip(words,matrices)
        c=LinearAlgebra.tr(P*A)/size(A,1)
        abs(c)<1e-12 && continue
        m=Nx.orbit_coordinates!(table,w)
        m.re!=0 && Nx.push_term!(program.forms,c*m.real_coefficient,m.re)
        m.im!=0 && Nx.push_term!(program.forms,c*m.imag_coefficient,m.im)
    end
    Nx.close_form!(program.forms)
end

function independent_hamiltonian!(reference,problem,kind;family=reference.sizes.closure.family,operators=nothing)
    n=problem.objective.registry.nsites
    W=eltype(reference.table.representatives.keys)
    words=Nx.basis_rows(W,Nx.Pauli(),1:n,n)
    matrices=[physical_pauli(w,1:n) for w in words]
    H=sum(c*physical_pauli(w,1:n) for (c,w) in Nx.terms(problem.objective))
    p,t=reference.program,reference.table
    if kind isa Nx.StateOptimality
        for block in family
            rows=[physical_pauli(w,1:n) for w in block.rows]
            for j in eachindex(rows),i in j:length(rows)
                A,B=rows[i],rows[j]
                physical_form!(p,t,A'*H*B-(H*A'*B+A'*B*H)/2,words,matrices)
            end
            Nx.push_block!(p,Nx.HPSD,length(rows))
        end
    else
        if isnothing(operators)
            # Derive the operator domain from physical row products, rather
            # than the pipeline's packed multiplication or orbit selection.
            operators=filter(words) do w
                P=physical_pauli(w,1:n)
                any(family) do block
                    rows=[physical_pauli(v,1:n) for v in block.rows]
                    any(abs(LinearAlgebra.tr(P*A'*B))>1e-10 for A in rows,B in rows)
                end
            end
        end
        for w in operators
            U=physical_pauli(w,1:n)
            physical_form!(p,t,-im*(H*U-U*H),words,matrices)
        end
        isempty(operators) || Nx.push_block!(p,Nx.ZERO,length(operators))
    end
    p.nvars=length(t.labels)
    Nx.validate(p)
end

function equality_values(p,y,start)
    Float64[real(sum((p.forms.coefs[t]*y[p.forms.ids[t]] for t in Nx.term_range(p.forms,f));init=0im))
        for b in p.blocks[start:end] if b.kind==Nx.ZERO for f in Nx.form_range(b)]
end
function hamiltonian_slacks(p,y,start)
    psd=deepcopy(p)
    empty!(psd.blocks)
    append!(psd.blocks,filter(b->b.kind!=Nx.ZERO,p.blocks[start:end]))
    Float64[extra_values(psd,y,1)...]
end
function physical_coordinates(table,rho,n)
    [label.part==2 ? imag(LinearAlgebra.tr(rho*physical_pauli(label.word,1:n))) :
        real(LinearAlgebra.tr(rho*physical_pauli(label.word,1:n))) for label in table.labels]
end

@testset "Hamiltonian strengthening: independent reduced data" begin
    problem,registry=strengthening_fixture(3)
    W=Nx.word_type(registry,5)
    translation=Nx.SitePermutation([2,3,1])
    for kind in (Nx.Stationarity(),Nx.StateOptimality()),mode in (:plain,:translation,:reflection,:anti,:split,:dihedral)
        symmetry=mode==:plain ? Nx.FamilySymmetry() : mode==:reflection ?
            Nx.FamilySymmetry([Nx.SitePermutation([1,3,2])]) : mode==:anti ?
            Nx.FamilySymmetry([Nx.SitePermutation([1,3,2])];moment_only=(Nx.AntiunitaryPauli(one(W)),)) :
            mode==:dihedral ? Nx.FamilySymmetry([translation,Nx.SitePermutation([1,3,2])]) : Nx.FamilySymmetry([translation])
        built=Nx.build(problem,Nx.BlockRelaxation(;order=1,symmetry,
            split=mode==:split ? Nx.TermSplit() : nothing,strengthen=(kind,)))
        report=built.sizes.closure.strengthening
        reference=Nx.build(problem,Nx.BlockRelaxation(order=1))
        start=length(reference.program.blocks)+1
        family=built.sizes.closure.family[built.sizes.closure.representatives]
        independent_hamiltonian!(reference,problem,kind;family=kind isa Nx.Stationarity ? built.sizes.closure.family : family)
        coordinates=[Nx.orbit_coordinates!(built.table,label.word) for label in reference.table.labels]
        if kind isa Nx.Stationarity
            # Compare the complete equality spaces, including operators
            # omitted by orbit reduction. No implementation-selected domain.
            columns=[Float64[i==j for i in 1:length(built.table.labels)] for j in 1:length(built.table.labels)]
            A=hcat([equality_values(built.program,v,report.first_block) for v in columns]...)
            B=hcat([equality_values(reference.program,real.([coordinate_value(c,v) for c in coordinates]),start) for v in columns]...)
            @test rank(A;atol=1e-9)==rank(B;atol=1e-9)==rank(vcat(A,B);atol=1e-9)
        end
        for seed in 1:3
            values=randn(MersenneTwister(seed+19),length(built.table.labels)); values[1]=1
            other=[real(coordinate_value(c,values)) for c in coordinates]
            if kind isa Nx.StateOptimality
                reduced=hamiltonian_slacks(built.program,values,report.first_block)
                full=hamiltonian_slacks(reference.program,other,start)
                # Nonabelian minimal projections omit repeated eigenvalues.
                # Check both inclusions of the spectrum, not just its minimum.
                @test all(v->minimum(abs.(full.-v))<1e-9,reduced)
                @test all(v->minimum(abs.(reduced.-v))<1e-9,full)
            end
        end
        @test occursin(kind isa Nx.Stationarity ? "eigenstates of the objective" : "ground space of the objective",sprint(show,built.valid_for))
        @test occursin("equalities",sprint(show,Nx.inspect(problem,Nx.BlockRelaxation(order=1,strengthen=(kind,)))))
    end
end

@testset "Hamiltonian strengthening: physical validity and guards" begin
    registry,(x,y,z)=Nx.pauli_variables(2)
    H=z[1]*z[2] # A two-dimensional ground space, not a unique vector.
    problem=Nx.Problem(H)
    E=eigen(Hermitian(sum(c*physical_pauli(w,1:2) for (c,w) in Nx.terms(H))))
    @test count(v->v≈minimum(E.values),E.values)==2
    swap=ComplexF64[1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1]
    for kind in (Nx.Stationarity(),Nx.StateOptimality()),symmetric in (false,true)
        symmetry=Nx.FamilySymmetry(symmetric ? [Nx.SitePermutation([2,1])] : Nx.AbstractAction[])
        built=Nx.build(problem,Nx.BlockRelaxation(;order=1,strengthen=(kind,),symmetry))
        start=built.sizes.closure.strengthening.first_block
        # A nontrivial mixture and a coherent superposition in the ground
        # space. Both must work, as must their symmetry averages.
        v=(E.vectors[:,1]+im*E.vectors[:,2])/sqrt(2)
        for rho in (v*v',0.3E.vectors[:,1]*E.vectors[:,1]'+0.7E.vectors[:,2]*E.vectors[:,2]')
            symmetric && (rho=(rho+swap*rho*swap')/2)
            values=physical_coordinates(built.table,rho,2)
            @test maximum(abs,equality_values(built.program,values,start);init=0.)<1e-10
            @test minimum(hamiltonian_slacks(built.program,values,start);init=0.)>=-1e-10
        end
        # An excited eigenstate violates optimality but obeys stationarity.
        excited=E.vectors[:,end]; values=physical_coordinates(built.table,excited*excited',2)
        if kind isa Nx.StateOptimality
            @test minimum(hamiltonian_slacks(built.program,values,start)) < -1
        elseif !symmetric
            @test maximum(abs,equality_values(built.program,values,start);init=0.)<1e-10
            v=(E.vectors[:,1]+im*E.vectors[:,end])/sqrt(2)
            values=physical_coordinates(built.table,v*v',2)
            @test maximum(abs,equality_values(built.program,values,start))>1
        end
        @test_throws ArgumentError Nx.build(Nx.Problem(H;constraints=(Nx.Positive(z[1]),)),Nx.BlockRelaxation(order=1,strengthen=(kind,)))
        _,(a,ad)=Nx.fermionic_variables(2)
        @test_throws ArgumentError Nx.build(Nx.Problem(ad[1]*a[1]),Nx.BlockRelaxation(order=1,strengthen=(kind,)))
    end
    # Antiunitary odd u can have zero moment while -i[H,u] is nonzero.
    registry,(x,y,z)=Nx.pauli_variables(1)
    W=Nx.word_type(registry,3)
    built=Nx.build(Nx.Problem(z[1]),Nx.BlockRelaxation(order=1,strengthen=(Nx.Stationarity(),),
        symmetry=Nx.FamilySymmetry(;moment_only=(Nx.AntiunitaryPauli(one(W)),))))
    values=physical_coordinates(built.table,ComplexF64[1 1;1 1]/2,1)
    @test maximum(abs,equality_values(built.program,values,built.sizes.closure.strengthening.first_block))≈2
end

@testset "Hamiltonian strengthening: independent bounds and combinations" begin
    registry,(x,y,z)=Nx.pauli_variables(3)
    H=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
    problem=Nx.Problem(H)
    exact=eigmin(Hermitian(sum(c*physical_pauli(w,1:3) for (c,w) in Nx.terms(H))))
    for kind in (Nx.Stationarity(),Nx.StateOptimality())
        built=Nx.build(problem,Nx.BlockRelaxation(order=1,strengthen=(kind,),symmetry=Nx.FamilySymmetry([Nx.SitePermutation([1,3,2])])))
        reference=Nx.build(problem,Nx.BlockRelaxation(order=1))
        independent_hamiltonian!(reference,problem,kind)
        bound=strengthening_bound(built.program)
        @test bound≈strengthening_bound(reference.program) atol=1e-6
        @test bound<=exact+1e-6
    end
    combined=Nx.build(problem,Nx.BlockRelaxation(order=1,split=Nx.TermSplit(),
        strengthen=(Nx.DensityMatrix(1:3),Nx.InvariantCone(1:3,"YII"),Nx.Stationarity(),Nx.StateOptimality())))
    @test strengthening_bound(combined.program)≈exact atol=1e-6

    # A degree-changing moment-only action must not add transformed row
    # blocks. Reference: plain family PLUS independently built moment
    # equalities, never the plain unreduced bound without those equalities.
    registry,(x,y,z)=Nx.pauli_variables(2)
    problem=Nx.Problem(z[2]+z[1]*z[2])
    W=Nx.word_type(registry,4)
    cnot=Nx.clifford_cnot(W,2,1,2)
    U=ComplexF64[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
    words=Nx.basis_rows(W,Nx.Pauli(),1:2,2)
    matrices=[physical_pauli(w,1:2) for w in words]
    for kind in (Nx.Stationarity(),Nx.StateOptimality())
        built=Nx.build(problem,Nx.BlockRelaxation(order=1,strengthen=(kind,),
            symmetry=Nx.FamilySymmetry(;moment_only=(cnot,))))
        reference=Nx.build(problem,Nx.BlockRelaxation(order=1))
        independent_hamiltonian!(reference,problem,kind)
        for P in matrices
            physical_form!(reference.program,reference.table,P-U*P*U',words,matrices)
        end
        Nx.push_block!(reference.program,Nx.ZERO,length(words))
        reference.program.nvars=length(reference.table.labels)
        @test built.sizes.closure.closed_blocks==1
        @test built.sizes.closure.strengthening.closed_cones==(kind isa Nx.StateOptimality ? 1 : 0)
        @test strengthening_bound(built.program)≈strengthening_bound(reference.program) atol=1e-6
    end
end
