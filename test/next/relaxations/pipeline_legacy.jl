module PipelineLegacyChecks
using Test, LinearAlgebra, TOML
using NCTSSoS.Next
import JuMP, COSMO

word(p)=last(only(Next.terms(p)))
backends(solver,sos_solver)=isnothing(sos_solver) ? (solver,) : (solver,sos_solver)

function sizes(id,built)
    data=TOML.parsefile(joinpath(@__DIR__,"../../data/expectations/next_pipeline_relaxations.toml"))
    expected=only(filter(c->c["id"]==id,data["cases"]))["expected"]
    @test sum(length,Next.block_rows(built))==expected["rows"]
    @test built.program.nvars==expected["real_variables"]
    @test sort(Next.block_sides(built))==expected["sides"]
end

function antiunitary(solver,sos_solver)
    @testset "Explicit antiunitary legacy controls" begin
        _,(x,y,z)=pauli_variables(2)
        h=x[1]+2z[1]+y[1]*y[2]
        rows=[one(h),x[1],y[1],z[1],y[2]]
        conjugation=AntiunitaryPauli(word(one(h)))
        built=Next.build(Problem(h;symmetry=(conjugation,)),Relaxation(rows=rows))
        sizes("antiunitary_custom",built)
        for backend in backends(solver,sos_solver)
            plain=solve(Problem(h),Relaxation(rows=rows);solver=backend)
            reduced=solve(Problem(h;symmetry=(conjugation,)),Relaxation(rows=rows);solver=backend)
            @test reduced.bound ≈ -sqrt(6) atol=2e-6
            @test reduced.bound ≈ plain.bound atol=2e-6
        end
        # Odd-site time reversal has theta^2=-I, but mixed invariant states exist.
        _,(x,y,z)=pauli_variables(3)
        h=sum(a[i]*a[i+1] for a in (x,y,z) for i in 1:2)
        reversal=AntiunitaryPauli(word(prod(y)))
        plain=solve(Problem(h),Relaxation(order=1);solver)
        reduced=solve(Problem(h;symmetry=(reversal,)),Relaxation(order=1);solver)
        @test reduced.bound ≈ plain.bound atol=3e-6

        _,(x,y,z)=pauli_variables(2)
        h=x[1]+0.3z[2]
        constraints=(Positive(1+z[1]),Positive(2+y[1]*y[2]))
        conjugation=AntiunitaryPauli(word(one(h)))
        for (grouping,split) in ((nothing,nothing),(Correlative(algorithm=MF()),TermSplit(algorithm=MF())))
            relaxation=Relaxation(order=2;grouping,split,strengthen=Sector(z[2],1))
            plain=solve(Problem(h;constraints),relaxation;solver)
            reduced=solve(Problem(h;constraints,symmetry=(conjugation,)),relaxation;solver)
            @test reduced.bound ≈ plain.bound atol=3e-6
        end
    end
end

function continuous(solver,sos_solver)
    @testset "Built-in continuous legacy controls" begin
        r,(x,y,z)=pauli_variables(2)
        h=x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
        rotations=spin_rotations(r)
        for (objective,constraints,expected) in ((h,(),-3),(-h,(),-3),(-h,(Zero(h-1),),-1))
            for backend in backends(solver,sos_solver)
                plain=solve(Problem(objective;constraints),Relaxation(order=1);solver=backend)
                reduced=solve(Problem(objective;constraints,symmetry=rotations),Relaxation(order=1);solver=backend)
                @test plain.bound ≈ expected atol=2e-5
                @test reduced.bound ≈ expected atol=2e-5
            end
        end
        plain=solve(Problem(h),Relaxation(order=2);solver)
        reduced=solve(Problem(h;symmetry=rotations),Relaxation(order=2);solver)
        @test plain.bound ≈ -3 atol=2e-5
        @test reduced.bound ≈ plain.bound atol=2e-5

        r,(a,ad)=fermionic_variables(2)
        h=ad[1]*a[2]+ad[2]*a[1]
        for backend in backends(solver,sos_solver)
            result=solve(Problem(h;symmetry=(NumberRotation(r),)),Relaxation(order=1);solver=backend)
            @test result.bound ≈ -1 atol=2e-5
        end
    end
end

function interface(solver,sos_solver)
    @testset "Complex equalities and affine squares" begin
        _,(x,y,z)=pauli_variables(1)
        @test_throws ArgumentError Problem(x[1];constraints=(Positive(x[1]+im*y[1]),))
        @test_throws ArgumentError Problem(x[1];constraints=(Zero(x[1]+im*y[1]),))
        for backend in backends(solver,sos_solver)
            result=solve(Problem((x[1]+y[1]+z[1])/2),Relaxation(order=1);solver=backend)
            @test result.bound ≈ -sqrt(3)/2 atol=1e-6
            # Next requires the two Hermitian parts explicitly.
            result=solve(Problem(x[1];constraints=(Zero(x[1]),Zero(y[1]))),Relaxation(order=1);solver=backend)
            @test result.bound ≈ 0 atol=1e-6
            _,(a,)=site_variables(Free(),1,2)
            result=solve(Problem((a[1]+a[2])^2+3),Relaxation(order=2);solver=backend)
            @test result.bound ≈ 3 atol=1e-6
            result=solve(Problem(a[1]^2+a[2]^2+1),Relaxation(order=1);solver=backend)
            @test result.bound ≈ 1 atol=1e-4
        end
    end
    @testset "COSMO early stop gives no bound" begin
        _,(xs,)=site_variables(Free(),1,1); x=only(xs)
        p=Problem(10x^2+2x;constraints=(Positive(1-x^2),))
        optimizer=JuMP.optimizer_with_attributes(COSMO.Optimizer,
            "verbose"=>false,"max_iter"=>30,"check_termination"=>1)
        stopped=solve(p,Relaxation(order=1);solver=JuMPBackend(optimizer))
        @test stopped.bound === nothing
        @test stopped.status==JuMP.MOI.ITERATION_LIMIT
        optimizer=JuMP.optimizer_with_attributes(COSMO.Optimizer,
            "verbose"=>false,"max_iter"=>10000,"check_termination"=>1)
        control=solve(p,Relaxation(order=1);solver=JuMPBackend(optimizer))
        @test control.status==JuMP.MOI.OPTIMAL
        @test control.bound ≈ -0.1 atol=1e-4
    end
end

function finite(solver,sos_solver)
    @testset "Finite and continuous non-tight controls" begin
        r,(x,y,z)=pauli_variables(3)
        h=sum(a[i]*a[j] for a in (x,y,z) for i in 1:3 for j in i+1:3)
        plain=solve(Problem(h),Relaxation(order=1);solver)
        reduced=solve(Problem(h;symmetry=(SitePermutation([2,3,1]),spin_rotations(r)...)),Relaxation(order=1);solver)
        @test plain.bound ≈ -4.5 atol=2e-5
        @test reduced.bound ≈ plain.bound atol=2e-5
        _,(x,y,z)=pauli_variables(2)
        h=(x[1]*x[2]+y[1]*y[2]+z[1]*z[2])/4
        rows=[one(h);x;y;z]
        result=solve(Problem(h;symmetry=(SitePermutation([2,1]),)),Relaxation(rows=rows);solver)
        @test result.bound ≈ -0.75 atol=1e-4
        _,(x,)=site_variables(Free(),1,2)
        h=x[1]*x[2]+x[2]*x[1]
        constraints=(Positive(1-x[1]^2-x[2]^2),)
        swap=LetterPermutation(x[1]=>x[2],x[2]=>x[1])
        plain=solve(Problem(h;constraints),Relaxation(order=1);solver)
        reduced=solve(Problem(h;constraints,symmetry=(swap,)),Relaxation(order=1);solver)
        @test plain.bound ≈ -1 atol=1e-5
        @test reduced.bound ≈ plain.bound atol=1e-5
        result=solve(Problem((x[1]*x[2]-x[2]*x[1])^2;
            constraints=(Zero(x[1]^2+x[2]^2-1),)),Relaxation(order=2);solver)
        @test result.bound ≈ -1 atol=2e-5
    end
    @testset "Signed localizers and the site-algebra sector limit" begin
        _,(x,)=site_variables(Unipotent(),1,3)
        @test_throws ArgumentError Next.build(Problem(x[1]),Relaxation(order=1,
            grouping=Correlative(),strengthen=Sector(x[2],-1)))
        _,(x,)=site_variables(Free(),1,2)
        h=x[1]^4+x[2]^4-x[1]^2-2x[2]^2
        constraints=(Positive(3-x[1]^2-x[2]^2),)
        flip=LetterPermutation(x[1]=>-x[1])
        relaxation=Relaxation(order=2,split=TermSplit())
        plain=solve(Problem(h;constraints),relaxation;solver)
        reduced=solve(Problem(h;constraints,symmetry=(flip,)),relaxation;solver)
        @test reduced.bound ≈ -1.25 atol=1e-6
        @test reduced.bound ≈ plain.bound atol=1e-6
    end
end

function custom_rows(solver)
    @testset "Complete custom rows preserve constrained and functional programs" begin
        r,(x,)=site_variables(Free(),1,2)
        h=x[1]^2+x[2]^2+1
        problems=Any[(r,Problem(h;constraints=(Positive(1-x[1]*x[2]-x[2]*x[1]),)),2)]
        r,(x,)=site_variables(Unipotent(),1,3;functional=Next.RealState())
        push!(problems,(r,Problem(state_symbol(x[1])*state_symbol(x[3])+state_symbol(x[2])),1))
        r,(x,)=site_variables(Free(),1,2;functional=Next.RealTrace())
        push!(problems,(r,Problem(trace_symbol(x[1]^2+x[2]*x[1]^2*x[2]+1)),2))
        for (registry,problem,order) in problems
            regular=Next.build(problem,Relaxation(;order))
            custom=Next.build(problem,Relaxation(rows=collect(Next.basis(registry,order))))
            @test Set.(Next.block_rows(custom))==Set.(Next.block_rows(regular))
            @test solve(custom;solver).bound ≈ solve(regular;solver).bound atol=1e-6
        end
    end
end

function bosonic_sector(solver)
    @testset "Bosonic fixed-number contact model without signed reflection" begin
        _,(a,ad)=bosonic_variables(2)
        s=inv(sqrt(2pi))
        h=0.5ad[1]*a[1]+1.5ad[2]*a[2]
        for i in 1:2,j in 1:2,k in 1:2,l in 1:2
            odd=count(==(2),(i,j,k,l))
            isodd(odd) && continue
            h+=0.5s*(odd==0 ? 1.0 : odd==2 ? 0.5 : 0.75)*ad[i]*ad[j]*a[l]*a[k]
        end
        rows=[one(h);a;[a[i]*a[j] for i in 1:2 for j in i:2];[ad[i]*a[j] for i in 1:2 for j in 1:2]]
        number=sum(ad[i]*a[i] for i in 1:2)
        # Independent occupation matrix on |20>, |11>, |02>.
        expected=eigmin(Symmetric([1+s 0 s/2;0 2+s 0;s/2 0 3+3s/4]))
        for (name,split) in (("bosonic_sector",nothing),("bosonic_sector_split",TermSplit()))
            built=Next.build(Problem(h),Relaxation(;rows,split,strengthen=Sector(number,2)))
            sizes(name,built)
            result=solve(built;solver)
            @test result.bound ≈ expected atol=1e-6
        end
    end
end

function selected_rows(solver)
    @testset "Contiguous row selection and density controls" begin
        _,(x,y,z)=pauli_variables(4)
        for (degree,count) in ((0,1),(1,13),(2,49),(3,157),(4,238))
            built=Next.build(Problem(zero(x[1])),Relaxation(order=degree,row_rule=Contiguous(periodic=true)))
            rows=only(collect(Next.block_rows(built)))
            @test length(rows)==count
            @test any(isone,rows)
            @test all(w->Next.degree(w)<=degree,rows)
        end
        built=Next.build(Problem(zero(x[1])),Relaxation(order=2,row_rule=Contiguous()))
        @test length(only(collect(Next.block_rows(built))))==40
        translation=SitePermutation([2,3,4,1])
        for h in (sum(x),sum(y),sum(x[i]*z[mod1(i+1,4)] for i in 1:4))
            plain=solve(Problem(h),Relaxation(order=1);solver)
            reduced=solve(Problem(h;symmetry=(translation,)),Relaxation(order=1);solver)
            @test reduced.bound ≈ plain.bound atol=1e-6
        end
        r,(x,y,z)=pauli_variables(2)
        h=-(x[1]*x[2]+y[1]*y[2]+z[1]*z[2])/4
        plain=solve(Problem(h),Relaxation(order=1);solver)
        density=solve(Problem(h),Relaxation(order=1,strengthen=DensityMatrix(1:2));solver)
        reduced=solve(Problem(h;symmetry=(SitePermutation([2,1]),spin_rotations(r)...)),
            Relaxation(order=1,strengthen=DensityMatrix(1:2));solver)
        @test plain.bound ≈ -0.75 atol=2e-5
        @test density.bound ≈ -0.25 atol=2e-5
        @test reduced.bound ≈ -0.25 atol=2e-5
        r,ops=pauli_variables(3)
        h=sum(a[i]*a[mod1(i+1,3)] for a in ops for i in 1:3)/4
        for symmetry in ((),spin_rotations(r))
            result=solve(Problem(h;symmetry),Relaxation(order=1,strengthen=DensityMatrix(1:3));solver)
            @test result.bound ≈ -0.75 atol=2e-5
        end
    end
    @testset "Detected small signed and phase symmetries" begin
        _,(x,y,z)=pauli_variables(2)
        for h in (x[1]+z[1],x[1],x[1]-z[1],x[1]*x[2]+y[1]*y[2])
            plain=solve(Problem(h),Relaxation(order=1);solver)
            reduced=solve(Problem(h;symmetry=detect_symmetry(h)),Relaxation(order=1);solver)
            @test reduced.bound ≈ plain.bound atol=1e-6
        end
    end
    @testset "Sparse involution swap" begin
        _,(x,)=site_variables(Unipotent(),1,2)
        h=-sum(x)
        symmetry=(LetterPermutation(x[1]=>x[2],x[2]=>x[1]),)
        for relaxation in (Relaxation(order=1,grouping=Correlative()),
                Relaxation(rows=[one(h);x],split=TermSplit(algorithm=MMD())))
            @test solve(Problem(h;symmetry),relaxation;solver).bound ≈ -2 atol=1e-6
        end
    end
end

function fermionic_controls(solver,sos_solver)
    @testset "Number rotations with sparse complex hopping and custom rows" begin
        r,(a,ad)=fermionic_variables(4)
        h=-(ad[1]*a[2]+ad[2]*a[1])-2im*ad[3]*a[4]+2im*ad[4]*a[3]+0.3ad[3]*a[3]
        relaxation=Relaxation(order=1,grouping=Correlative(algorithm=MMD()),split=TermSplit(algorithm=MMD()))
        for backend in backends(solver,sos_solver)
            result=solve(Problem(h;symmetry=(NumberRotation(r),)),relaxation;solver=backend)
            @test result.bound ≈ -1+eigmin(Hermitian([0.3 -2im;2im 0])) atol=2e-6
        end
        number=ad[1]*a[1]+ad[2]*a[2]
        result=solve(Problem(number;symmetry=(NumberRotation(r),)),
            Relaxation(rows=[one(h);a[1:2]],strengthen=Sector(number,1));solver)
        @test result.bound ≈ 1 atol=2e-6
        r,(a,ad)=fermionic_variables(2)
        numbers=[ad[i]*a[i] for i in 1:2]
        h=-sum(numbers)
        plain=solve(Problem(h),Relaxation(rows=[one(h);numbers]);solver)
        reduced=solve(Problem(h;symmetry=(NumberRotation(r),)),Relaxation(rows=[one(h);numbers]);solver)
        @test reduced.bound ≈ plain.bound atol=1e-6
        @test reduced.bound ≈ -2 atol=1e-6
        _,(_,_,z)=pauli_variables(2)
        for h in (sum(z),h)
            relaxation=Relaxation(order=1,grouping=Correlative(),split=TermSplit())
            plain=solve(Problem(h),relaxation;solver)
            reduced=solve(Problem(h;symmetry=(SitePermutation([2,1]),)),relaxation;solver)
            @test reduced.bound ≈ -2 atol=1e-6
            @test reduced.bound ≈ plain.bound atol=1e-6
        end
    end
end

function short_orbits()
    @testset "Short translation orbits retain the thermal spectrum" begin
        r,ops=pauli_variables(4)
        paulis=([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
        identity=Matrix{ComplexF64}(I,16,16)
        operators=Dict(word(axis[i])=>reduce(kron,[j==i ? p : Matrix{ComplexF64}(I,2,2) for j in 1:4])
            for (axis,p) in zip(ops,paulis) for i in 1:4)
        matrix(w)=foldl(*,(operators[l] for l in Next.word_letters(w));init=identity)
        thermal=sum(weight*operators[word(axis[i])]*operators[word(axis[mod1(i+1,4)])]
            for (axis,weight) in zip(ops,(0.7,1.3,0.2)) for i in 1:4)
        rho=exp(-thermal); rho/=tr(rho)
        h=sum(axis[i]*axis[mod1(i+1,4)] for axis in ops for i in 1:4)/4
        rows=collect(Next.basis(r,2)) # includes antipodal pairs, with orbit length two
        built=Next.build(Problem(h;symmetry=(SitePermutation([2,3,4,1]),)),Relaxation(;rows))
        values=[real(tr(rho*matrix(label.word))) for label in built.table.labels]
        spectrum=Float64[]
        # Read contract-2 affine blocks, not an implementation's row transform.
        for block in built.program.blocks
            block.kind in (Next.PSD,Next.HPSD) || continue
            m=zeros(ComplexF64,block.side,block.side)
            for j in 1:block.side,i in j:block.side
                f=block.first_form+Next.tri_offset(block.side,i,j)
                m[i,j]=sum((built.program.forms.coefs[t]*values[built.program.forms.ids[t]]
                    for t in Next.term_range(built.program.forms,f));init=0im)
                m[j,i]=conj(m[i,j])
            end
            append!(spectrum,eigvals(Hermitian(m)))
        end
        matrices=matrix.(word.(rows))
        dense=[tr(rho*(a'*b)) for a in matrices,b in matrices]
        @test sort(spectrum) ≈ eigvals(Hermitian(dense)) atol=1e-10
    end
end

function guards()
    @testset "Public input guards retained from internal ledger" begin
        r,(x,y,z)=pauli_variables(2)
        @test_throws ArgumentError CliffordMap(word.([one(x[1]),x[2],z[1],z[2]]))
        @test_throws ArgumentError CliffordMap(word.([x[1],z[1]]),[false])
        @test_throws ArgumentError Next.build(Problem(z[1];symmetry=(
            CliffordMap(word.([z[1],x[2],x[1],z[2]])),SpinRotation(r,:z))),Relaxation(order=1))
        @test_throws ArgumentError DensityMatrix([1,1])
        @test_throws ArgumentError DensityMatrix(Int[])
        @test_throws ArgumentError Next.build(Problem(z[1]),Relaxation(order=1,strengthen=DensityMatrix(1:3)))
        @test_throws ArgumentError Problem(z[1];symmetry=(:reflection,))
        @test_throws ArgumentError Relaxation(rows=[one(x[1]),(1,2)])
        free,_=site_variables(Free(),1,2)
        @test_throws ArgumentError SpinRotation(free,:z)
        for functional in (Next.RealState(),Next.RealTrace())
            _,(x,)=site_variables(Unipotent(),1,1;functional)
            symbol=functional isa Next.RealState ? state_symbol : trace_symbol
            @test_throws ArgumentError Relaxation(rows=[one(x[1]),2symbol(x[1])])
            @test_throws ArgumentError Relaxation(rows=[one(x[1]),symbol(x[1])+one(x[1])])
            @test_throws ArgumentError Next.build(Problem(symbol(x[1])),Relaxation(rows=[one(x[1])]))
        end
    end
end

function chain_references(solver)
    @testset "Different programs: nested chain density strengthening" begin
        n=6
        r,(x,y,z)=pauli_variables(n)
        h=sum(a[i]*a[mod1(i+1,n)] for a in (x,y,z) for i in 1:n)/4
        translation=SitePermutation([mod1(i+1,n) for i in 1:n])
        reflection=SitePermutation([mod1(2-i,n) for i in 1:n])
        flipx=CliffordMap(word.([x;z]),[falses(n);trues(n)])
        flipz=CliffordMap(word.([x;z]),[trues(n);falses(n)])
        basic=(translation,flipx,flipz)
        mirrored=(basic...,reflection,AntiunitaryPauli(word(one(h))))
        relaxation=Relaxation(order=2,row_rule=Contiguous(periodic=true))
        base=solve(Problem(h;symmetry=basic),relaxation;solver).bound
        mirror=solve(Problem(h;symmetry=mirrored),relaxation;solver).bound
        pair=solve(Problem(h;symmetry=mirrored),Relaxation(order=2,row_rule=Contiguous(periodic=true),
            strengthen=DensityMatrix(1:2));solver).bound
        four=solve(Problem(h;symmetry=mirrored),Relaxation(order=2,row_rule=Contiguous(periodic=true),
            strengthen=(DensityMatrix(1:2),DensityMatrix(1:4)));solver).bound
        spin=solve(Problem(h;symmetry=(mirrored...,spin_rotations(r)...)),
            Relaxation(order=2,row_rule=Contiguous(periodic=true),
                strengthen=(DensityMatrix(1:2),DensityMatrix(1:4)));solver).bound
        @test mirror>=base-1e-6
        @test mirror<=-2.802775637+1e-4
        @test pair>=mirror-1e-6
        @test four>=pair-1e-6
        @test four>mirror+1e-3
        @test four<=-2.802775637+1e-5
        @test spin>=four-1e-6
        @test spin<=-2.802775637+1e-5
    end
    @testset "Different programs: odd-ring density positivity retains all spins" begin
        n=5
        r,ops=pauli_variables(n)
        h=sum(a[i]*a[mod1(i+1,n)] for a in ops for i in 1:n)/4
        patches=[sort(mod1.(i:i+3,n)) for i in 1:n]
        triples=[Set(mod1.(i:i+2,n)) for i in 1:n]
        rows=filter(collect(Next.basis(r,3))) do p
            support=Set(Next.sites(word(p)))
            any(t->issubset(support,t),triples)
        end
        paulis=([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
        matrix=sum(reduce(kron,[j in (i,mod1(i+1,n)) ? p : Matrix{ComplexF64}(I,2,2)
            for j in 1:n]) for i in 1:n for p in paulis)/4
        symmetry=(SitePermutation([mod1(i+1,n) for i in 1:n]),
            SitePermutation([mod1(2-i,n) for i in 1:n]),spin_rotations(r)...)
        for sign in (1,-1)
            result=solve(Problem(sign*h;symmetry),Relaxation(;rows,
                strengthen=Tuple(DensityMatrix(p) for p in patches));solver)
            @test result.bound ≈ eigmin(Hermitian(sign*matrix)) atol=2e-5
        end
    end
end

function fixed_point_reference(solver,sos_solver)
    isnothing(sos_solver) && return
    @testset "Fixed-point inequalities on both solver forms" begin
        _,(x,)=site_variables(Free(),1,2)
        p=Problem((x[1]+x[2])^2+3;constraints=(Positive(-10-x[1]),
            Positive(-10-x[2]),Positive(x[1]+10),Positive(x[2]+10)))
        primal=solve(p,Relaxation(order=2);solver)
        dual=solve(p,Relaxation(order=2);solver=sos_solver)
        @test primal.bound ≈ dual.bound atol=1e-3
    end
end

function check(;solver,sos_solver=nothing,references=false)
    guards()
    antiunitary(solver,sos_solver)
    continuous(solver,sos_solver)
    interface(solver,sos_solver)
    finite(solver,sos_solver)
    custom_rows(solver)
    bosonic_sector(solver)
    selected_rows(solver)
    fermionic_controls(solver,sos_solver)
    short_orbits()
    references && chain_references(solver)
    references && fixed_point_reference(solver,sos_solver)
end
end
