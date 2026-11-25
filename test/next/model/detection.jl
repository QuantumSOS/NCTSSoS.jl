function detection_ring(n,kind)
    _,(x,y,z)=Nx.pauli_variables(n)
    h=kind==:heisenberg ? sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n) :
        sum(-z[i]*z[mod1(i+1,n)]-0.7*x[i] for i in 1:n)
    return Nx.Problem(h)
end

@testset "Detected symmetry without nontrivial generators" begin
    @testset "Manual provenance" begin
        manual=Nx.FamilySymmetry([Nx.SitePermutation([2,1])])
        @test sprint(show,manual)=="Symmetry: 1 block generators, 0 moment-only generators, 0 continuous generators"
        @test sprint(show,MIME"text/plain"(),Nx.FamilySymmetry())==
            "Symmetry: 0 block generators, 0 moment-only generators, 0 continuous generators"
        registry,_=Nx.pauli_variables(4)
        continuous=Nx.Symmetry(Nx.SitePermutation([2,3,4,1]),Nx.SitePermutation([4,3,2,1]);
            continuous=Nx.spin_rotations(registry))
        for render in (repr,s->repr(MIME"text/plain"(),s))
            @test render(continuous)=="Symmetry: 2 block generators, 0 moment-only generators, 3 continuous generators"
        end
    end
    @testset "Asymmetric Hamiltonian" begin
        rng=MersenneTwister(414)
        _,ops=Nx.pauli_variables(4)
        h=sum(randn(rng)*ops[a][i]*ops[b][j] for i in 1:4 for j in i+1:4 for a in 1:3 for b in 1:3)+
            sum(randn(rng)*ops[a][i] for i in 1:4 for a in 1:3)
        problem=Nx.Problem(h)
        symmetry=Nx.detect_symmetry(problem)
        @test isempty(symmetry.generators) && isempty(symmetry.moment_only)
        @test occursin("group order not reported by the detector",sprint(show,symmetry))
        detected=Nx.solve(problem,Nx.BlockRelaxation(;order=1,symmetry);solver=Nx.JuMPBackend(SOLVER))
        plain=Nx.solve(problem,Nx.BlockRelaxation(order=1);solver=Nx.JuMPBackend(SOLVER))
        @test detected.bound ≈ plain.bound atol=1e-6
        @test detected.sizes.program==plain.sizes.program
    end
end

@testset "Detected symmetry in user problems" begin
    @testset "Independent local signs exceed the finite group limit" begin
        _,(_,_,z)=Nx.pauli_variables(16)
        H=sum(i*z[i] for i in eachindex(z))
        symmetry=Nx.detect_symmetry(H)
        built=Nx.build(Nx.Problem(H;symmetry),Nx.Relaxation(order=1))
        # Identity and all Z rows share a character; X_i and Y_i share another.
        @test sort([b.side for b in built.program.blocks])==[fill(2,16);17]
    end
    @testset "$kind ($n sites)" for (n,kind) in ((8,:heisenberg),(16,:heisenberg),(8,:ising))
        problem=detection_ring(n,kind)
        symmetry=Nx.detect_symmetry(problem)
        W=eltype(problem.objective.data.words)
        translation=Nx.SitePermutation([mod1(i+1,n) for i in 1:n])
        @test any(all(Nx.act(g,Nx.pauli_letter(W,i,a))==Nx.act(t,Nx.pauli_letter(W,i,a))
            for i in 1:n for a in 0:2) for g in symmetry.generators for t in (translation,inv(translation)))
        @test all(Nx.degree(Nx.act(g,Nx.pauli_letter(W,i,a))[2])==1
            for g in symmetry.generators for i in 1:n for a in 0:2)
        @test all(any(Nx.degree(Nx.act(g,Nx.pauli_letter(W,i,a))[2])>1
            for i in 1:n for a in 0:2) for g in symmetry.moment_only)
        description=Nx.BlockRelaxation(;order=1,symmetry)
        reduced=Nx.solve(problem,description;solver=Nx.JuMPBackend(SOLVER))
        unreduced=Nx.solve(problem,Nx.BlockRelaxation(order=1);solver=Nx.JuMPBackend(SOLVER))
        @test reduced.bound ≈ unreduced.bound atol=1e-6
        @test reduced.sizes.program.real_variables < unreduced.sizes.program.real_variables
        @test occursin("$(length(symmetry.generators)) block generators",sprint(show,symmetry))
        @test occursin("$(length(symmetry.moment_only)) moment-only generators",sprint(show,symmetry))
        @test occursin("group order not reported by the detector; the result can be a subgroup",sprint(show,MIME"text/plain"(),symmetry))
    end

    # Independent reference: unreduced PSD matrix plus explicit symmetry
    # equations on all 4-site Pauli words, including words outside its rows.
    problem=detection_ring(4,:heisenberg)
    symmetry=Nx.detect_symmetry(problem)
    @test !isempty(symmetry.moment_only)
    plain=Nx.build(problem,Nx.Relaxation(order=1))
    model,y=Nx.to_jump(plain.program)
    W=eltype(problem.objective.data.words)
    moments=Dict{W,AffExpr}(one(W)=>AffExpr(1.0))
    for (id,label) in enumerate(plain.table.labels)
        id==1 && continue
        moments[label.word]=1.0*y[id]
    end
    domain=Nx.basis_rows(W,Nx.Pauli(),1:4,4)
    for w in domain
        haskey(moments,w) || (moments[w]=1.0*@variable(model))
    end
    for g in (symmetry.generators...,symmetry.moment_only...),w in domain
        phase,image=Nx.act(g,w)
        @constraint(model,moments[w]==real(Nx.scale(phase,1.0+0im))*moments[image])
    end
    set_optimizer(model,SOLVER)
    optimize!(model)
    reduced=Nx.solve(problem,Nx.BlockRelaxation(;order=1,symmetry);solver=Nx.JuMPBackend(SOLVER))
    @test reduced.bound ≈ objective_value(model) atol=1e-6

    # A cap reports its warning and still supplies a verified subgroup.
    partial=@test_logs (:warn,r"hit max_automorphisms=1") match_mode=:any Nx.detect_symmetry(problem;max_automorphisms=1)
    @test partial isa Nx.FamilySymmetry
    @test_throws ArgumentError Nx.detect_symmetry(problem;ga_backend=:unsupported)
    _,ops=Nx.site_variables(Nx.Free(),1,2)
    @test_throws r"Pauli problems only" Nx.detect_symmetry(Nx.Problem(ops[1][1]))
    @test_throws r"with constraints" Nx.detect_symmetry(Nx.Problem(problem.objective;constraints=(Nx.Positive(one(problem.objective)),)))
    for value in (0,1)
        constant=Nx.Problem(value*one(problem.objective))
        symmetry=Nx.detect_symmetry(constant)
        @test isempty(symmetry.generators) && isempty(symmetry.moment_only)
        result=Nx.solve(constant,Nx.BlockRelaxation(;order=0,symmetry);solver=Nx.JuMPBackend(SOLVER))
        @test result.bound==value
    end
    @test_throws r"cannot yet return expectations" Nx.moment_value(nothing,one(W),Float64[])

    # Unused registry sites stay fixed; unequal coefficients and local fields
    # prevent an axis/site permutation from passing by accidental symmetry.
    _,(x,y,z)=Nx.pauli_variables(3)
    p=2*x[1]*x[2]+3*y[1]*y[2]+0.7*z[1]*z[2]+0.31*z[1]+0.17*z[2]
    sym=Nx.detect_symmetry(Nx.Problem(p))
    expected=Dict(w=>c for (c,w) in Nx.terms(p))
    for g in (sym.generators...,sym.moment_only...)
        image=Dict(Nx.act(g,w)[2]=>Nx.scale(Nx.act(g,w)[1],c) for (c,w) in Nx.terms(p))
        @test image==expected
        @test all(Nx.act(g,op[3].data.words[1])==(one(Nx.Phase),op[3].data.words[1]) for op in (x,y,z))
    end

    # Detection must not require a registry name for every axis. These X/Z
    # products contain Y on two sites, neither of which has a registered name.
    partial=Nx.Registry(Nx.Pauli(),3,[:left=>(1,0),:right=>(1,2),:other=>(2,0),:axis=>(2,2)])
    v=Nx.variables(partial)
    h=v[:left]*v[:other]+v[:right]*v[:axis]-v[:left]*v[:right]*v[:other]*v[:axis]
    full=x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    sparse_sym,full_sym=Nx.detect_symmetry(h),Nx.detect_symmetry(full)
    @test sparse_sym.generators==full_sym.generators
    @test sparse_sym.moment_only==full_sym.moment_only
    forms=map(((h,sparse_sym),(full,full_sym))) do (objective,symmetry)
        built=Nx.build(Nx.Problem(objective;symmetry),Nx.Relaxation(order=1))
        Nx.moment_standard_form(built.program)
    end
    f,g=forms
    @test (f.A,f.b,f.c,f.offset,f.cones)==(g.A,g.b,g.c,g.offset,g.cones)
end
