# Independent Hilbert-space oracle: tensor the physical Pauli matrices, not
# the packed product/action implementation or the density builder.
function physical_pauli(w,sites)
    letters=Nx.legacy_letters(w)
    axes=Dict((Int(l)-1)÷3+1 => (Int(l)-1)%3+1 for l in letters)
    matrices=([1 0;0 1],[0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    foldl(kron,(ComplexF64.(matrices[get(axes,s,0)+1]) for s in sites))
end

function independent_extra!(program,table,W,sites;r=nothing)
    words=Nx.basis_rows(W,Nx.Pauli(),sites,length(sites))
    matrices=[physical_pauli(w,sites) for w in words]
    if isnothing(r)
        d=2^length(sites)
        for j in 1:d,i in j:d
            for (w,P) in zip(words,matrices)
                c=P[i,j]/d
                iszero(c) && continue
                id=Nx.orbit_coordinates!(table,w).re
                Nx.push_term!(program.forms,c,id)
            end
            Nx.close_form!(program.forms)
        end
        Nx.push_block!(program,Nx.HPSD,d)
    else
        paulis=Dict('I'=>[1 0;0 1],'X'=>[0 1;1 0],'Y'=>[0 -im;im 0],'Z'=>[1 0;0 -1])
        R=foldl(kron,(ComplexF64.(paulis[c]) for c in r))
        Nx.push_term!(program.forms,sqrt(2.0^(length(sites)-1)-1),1)
        Nx.close_form!(program.forms)
        count=1
        for (w,P) in zip(words,matrices)
            (isone(w) || R*conj(P)*R!=P) && continue
            Nx.push_term!(program.forms,1.,Nx.orbit_coordinates!(table,w).re)
            Nx.close_form!(program.forms); count+=1
        end
        Nx.push_block!(program,Nx.SOC,count)
    end
    program.nvars=length(table.labels)
    program
end

function extra_values(p,y,start)
    result=[]
    for b in p.blocks[start:end]
        values=[sum((p.forms.coefs[t]*y[p.forms.ids[t]] for t in Nx.term_range(p.forms,f));init=0im)
            for f in Nx.form_range(b)]
        if b.kind==Nx.SOC
            push!(result,real(values[1])-norm(values[2:end]))
        else
            A=zeros(ComplexF64,b.side,b.side)
            for j in 1:b.side,i in j:b.side
                A[i,j]=values[Nx.tri_offset(b.side,i,j)+1]; A[j,i]=conj(A[i,j])
            end
            append!(result,eigvals(Hermitian(A)))
        end
    end
    sort!(Float64.(result))
end

function strengthening_fixture(n=3)
    registry,(x,y,z)=Nx.pauli_variables(n)
    H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    Nx.Problem(H),registry
end

function strengthening_bound(program)
    # COSMO's complex chordal decomposition crashes on a reduced density
    # cone. Keep the cone unchanged and solve without that preprocessing.
    model,_=Nx.to_jump(Nx.moment_standard_form(program))
    set_optimizer(model,SOLVER)
    set_optimizer_attribute(model,"decompose",false)
    optimize!(model)
    @test termination_status(model)==JuMP.MOI.OPTIMAL
    objective_value(model)
end

@testset "Pauli strengthening: independent data and physical states" begin
    cases=[(n,kind) for n in (3,4) for kind in (Nx.DensityMatrix(1:3),Nx.InvariantCone(1:3,"YII"))]
    append!(cases,[(4,Nx.DensityMatrix(1:4)),(4,Nx.InvariantCone(1:4,"YIII"))])
    for (n,kind) in cases,symmetric in (false,true)
        problem,registry=strengthening_fixture(n)
        W=Nx.word_type(registry,4)
        g=Nx.SitePermutation([2:n;1])
        symmetry=Nx.FamilySymmetry(symmetric ? [g] : Nx.AbstractAction[])
        description=Nx.BlockRelaxation(order=1,split=Nx.TermSplit(),symmetry=symmetry,strengthen=(kind,))
        built=Nx.build(problem,description)
        report=built.sizes.closure.strengthening
        @test length(report.representatives)==1
        expected=symmetric ? (kind isa Nx.InvariantCone || length(kind.sites)<n ? n : 1) : 1
        @test report.closed_cones==expected
        @test report.added_equalities==0
        @test occursin("valid for each state",sprint(show,built.valid_for))
        @test occursin(string(nameof(typeof(kind))),sprint(show,Nx.inspect(problem,description)))

        # Reference uses an independent full subsystem matrix, with no
        # coefficient compression. At invariant moments all orbit images
        # are equivalent, so compare one representative's complete spectrum.
        reference=Nx.build(problem,Nx.BlockRelaxation(order=1,split=Nx.TermSplit()))
        start=length(reference.program.blocks)+1
        independent_extra!(reference.program,reference.table,W,kind.sites;r=kind isa Nx.InvariantCone ? kind.r : nothing)
        for seed in 1:3
            coordinates=[Nx.orbit_coordinates!(built.table,label.word) for label in reference.table.labels]
            values=randn(MersenneTwister(78+seed),length(built.table.labels)); values[1]=1
            previous=[label.part==2 ? imag(coordinate_value(c,values)) : real(coordinate_value(c,values))
                for (label,c) in zip(reference.table.labels,coordinates)]
            @test extra_values(built.program,values,report.first_block)≈extra_values(reference.program,previous,start) atol=1e-9
        end

        H=sum(c*physical_pauli(w,1:n) for (c,w) in Nx.terms(problem.objective))
        states=eigen(Hermitian(H))
        for rho in (states.vectors[:,1]*states.vectors[:,1]',states.vectors[:,end]*states.vectors[:,end]',
            0.3states.vectors[:,1]*states.vectors[:,1]'+0.7states.vectors[:,2]*states.vectors[:,2]')
            moments=map(built.table.labels) do label
                w=label.word
                value=0im
                for shift in 0:(symmetric ? n-1 : 0)
                    value+=LinearAlgebra.tr(rho*physical_pauli(w,1:n)); _,w=Nx.act(g,w)
                end
                value/=symmetric ? n : 1
                label.part==2 ? imag(value) : real(value)
            end
            @test minimum(extra_values(built.program,moments,report.first_block))>=-1e-9
        end
    end
end

@testset "Strengthening bounds, closure and guards" begin
    problem,registry=strengthening_fixture()
    W=Nx.word_type(registry,4)
    exact=eigmin(Hermitian(sum(c*physical_pauli(w,1:3) for (c,w) in Nx.terms(problem.objective))))
    for kind in (Nx.DensityMatrix(1:3),Nx.InvariantCone(1:3,"YII"))
        # Three explicitly named cones in the independent reference versus
        # one generator-closed representative. No invariant-moment shortcut
        # is applied to the reference solve.
        built=Nx.build(problem,Nx.BlockRelaxation(order=1,strengthen=(kind,),
            symmetry=Nx.FamilySymmetry([Nx.SitePermutation([2,3,1])])) )
        reference=Nx.build(problem,Nx.BlockRelaxation(order=1))
        for r in (kind isa Nx.InvariantCone ? ("YII","IYI","IIY") : (nothing,))
            independent_extra!(reference.program,reference.table,W,1:3;r)
        end
        bound=strengthening_bound(built.program)
        @test bound≈strengthening_bound(reference.program) atol=1e-6
        @test bound<=exact+1e-6
        kind isa Nx.DensityMatrix && @test bound≈exact atol=1e-6
    end
    # CNOT turns a one-site algebra into a non-site-local algebra. Closure
    # must retain its positivity, not replace it by a larger site density.
    registry,(x,y,z)=Nx.pauli_variables(2)
    W=Nx.word_type(registry,2)
    cnot=Nx.clifford_cnot(W,2,1,2)
    problem=Nx.Problem(z[1]+x[2])
    built=Nx.build(problem,Nx.BlockRelaxation(order=2,strengthen=(Nx.DensityMatrix([1]),),
        symmetry=Nx.FamilySymmetry([cnot])))
    @test built.sizes.closure.strengthening.closed_cones==2
    @test sum(b.side for b in built.program.blocks[built.sizes.closure.strengthening.first_block:end])==2
    @test strengthening_bound(built.program)≈-2 atol=1e-6
    @test_throws ArgumentError Nx.DensityMatrix([2,1])
    @test_throws ArgumentError Nx.InvariantCone(1:3,"III")
    @test_throws ArgumentError Nx.build(problem,Nx.BlockRelaxation(order=1,strengthen=(Nx.DensityMatrix([3]),)))
    r,(a,ad)=Nx.fermionic_variables(2)
    @test_throws r"require Pauli" Nx.build(Nx.Problem(ad[1]*a[1]),Nx.BlockRelaxation(order=1,strengthen=(Nx.DensityMatrix([1]),)))
end

@testset "Nonredundant named SOC and mixed antiunitary symmetry" begin
    registry,(x,y,z)=Nx.pauli_variables(3)
    H=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
    problem=Nx.Problem(H)
    W=Nx.word_type(registry,4)
    symmetry=Nx.FamilySymmetry([Nx.SitePermutation([1,3,2])];moment_only=(Nx.AntiunitaryPauli(Nx.pauli_letter(W,1,1)),))
    @test strengthening_bound(Nx.build(problem,Nx.BlockRelaxation(order=1)).program)≈-5 atol=1e-6
    states=eigen(Hermitian(sum(c*physical_pauli(w,1:3) for (c,w) in Nx.terms(H))))
    swap=zeros(ComplexF64,8,8)
    for i in 0:7
        j=(i&4)|((i&1)<<1)|((i&2)>>1)
        swap[j+1,i+1]=1
    end
    R=kron(ComplexF64[0 -im;im 0],Matrix{ComplexF64}(I,4,4))
    for (kind,expected) in ((Nx.InvariantCone(1:3,"YII"),-sqrt(15)),(Nx.DensityMatrix(1:3),-3.0))
        plain=Nx.build(problem,Nx.BlockRelaxation(order=1,strengthen=(kind,)))
        reduced=Nx.build(problem,Nx.BlockRelaxation(order=1,strengthen=(kind,),symmetry=symmetry))
        @test strengthening_bound(plain.program)≈expected atol=1e-6
        @test strengthening_bound(reduced.program)≈expected atol=1e-6
        coordinates=[Nx.orbit_coordinates!(reduced.table,l.word) for l in plain.table.labels]
        values=randn(MersenneTwister(91),length(reduced.table.labels)); values[1]=1
        previous=[l.part==2 ? imag(coordinate_value(c,values)) : real(coordinate_value(c,values)) for (l,c) in zip(plain.table.labels,coordinates)]
        @test extra_values(reduced.program,values,reduced.sizes.closure.strengthening.first_block)≈
            extra_values(plain.program,previous,plain.sizes.closure.strengthening.first_block) atol=1e-9
        for v in (states.vectors[:,1],states.vectors[:,end])
            rho=v*v'
            rho=(rho+swap*rho*swap'+R*conj(rho)*R'+swap*R*conj(rho)*R'*swap')/4
            moments=[real(LinearAlgebra.tr(rho*physical_pauli(l.word,1:3))) for l in reduced.table.labels]
            @test minimum(extra_values(reduced.program,moments,reduced.sizes.closure.strengthening.first_block))>=-1e-9
        end
    end
    optimizer=optimizer_with_attributes(COSMO.Optimizer,"decompose"=>false,"verbose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8)
    result=Nx.solve(problem,Nx.BlockRelaxation(order=1,symmetry=symmetry,
        strengthen=(Nx.DensityMatrix(1:3),Nx.InvariantCone(1:3,"YII")));solver=Nx.JuMPBackend(optimizer))
    @test result.status==JuMP.MOI.OPTIMAL
    @test result.bound≈-3 atol=1e-6
    @test length(result.valid_for.kinds)==2
end
