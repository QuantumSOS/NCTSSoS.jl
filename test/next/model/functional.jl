module FunctionalModelChecks
using Test, LinearAlgebra, NCTSSoS
import NCTSSoS: Next as N

const DOC_SIGNS=[1 1 1; 1 1 -1; 1 -1 0]
const BENCH_SIGNS=[1 1 -1; 1 1 1; -1 1 0]
covariance(A,B,symbol,signs)=-sum(signs[i,j]*(symbol(A[i]*B[j])-symbol(A[i])*symbol(B[j])) for i in 1:3,j in 1:3 if signs[i,j]!=0)
function covariance_model(trace=false;signs=DOC_SIGNS)
    if trace
        r,(x,)=N.site_variables(N.Unipotent(),1,6;functional=N.RealTrace())
        return r,covariance(x[1:3],x[4:6],N.trace_symbol,signs)
    end
    r,(A,B)=N.site_variables(N.Unipotent(),2,3;functional=N.RealState())
    r,covariance(A,B,N.state_symbol,signs)
end
function old_covariance(trace=false;signs=DOC_SIGNS)
    if trace
        r,(x,)=create_unipotent_variables([("x",1:6)])
        A,B=x[1:3],x[4:6]; symbol=NCTSSoS.tr
    else
        r,(A,B)=create_unipotent_variables([("A",1:3),("B",1:3)])
        symbol=ς
    end
    polyopt(covariance(A,B,symbol,signs)*one(A[1]),r)
end

function check(backend;factory=nothing,references=false)
    @testset "Functional modeling and outer moments" begin
        r,(x,)=N.site_variables(N.Unipotent(),1,3;functional=N.RealState())
        a,b,c=x; symbol=N.state_symbol
        @test symbol(a*b*c)==symbol(c*b*a)
        @test symbol(a*b*c)!=symbol(b*c*a)
        @test symbol(a*b)!=symbol(a)*symbol(b)
        @test symbol(symbol(a)*b)==symbol(a)*symbol(b)
        @test symbol(one(a))==one(a)
        @test_throws r"matching functional" N.trace_symbol(a)
        rt,(y,)=N.site_variables(N.Unipotent(),1,3;functional=N.RealTrace())
        @test_throws r"different registries" symbol(a)+N.trace_symbol(y[1])
        for constructor in (N.pauli_variables,N.fermionic_variables,N.bosonic_variables)
            _,ops=constructor(1)
            @test_throws r"unsupported" symbol(first(ops[1]))
        end
        @test occursin("⟨",sprint(show,symbol(a)))
        @test occursin("tr(",sprint(show,N.trace_symbol(y[1])))
        mixed_text=sprint(show,symbol(a)*b)
        @test occursin("⟨s1_1⟩",mixed_text) && occursin("s1_2",mixed_text)
        p=N.Problem(symbol(a)); relaxation=N.Relaxation(order=2)
        @test occursin("normalized real states",sprint(show,p))
        @test occursin("normalized real traces",sprint(show,N.Problem(N.trace_symbol(y[1]))))
        @test_throws r"declared symmetry.*remove symmetry" N.Problem(symbol(a);symmetry=(N.SitePermutation([1]),))
        built=N.build(p,relaxation)
        # An ensemble of the scalar evaluations a=b=c=+1 and a=b=c=-1.
        # Its outer moment of <a>^2 is one while the first moment is zero.
        vals=[isodd(N.degree(label.word)) ? 0.0 : 1.0 for label in built.table.labels[2:end]]
        result=N.SolveResult(0.0,:test,built.sizes,built.valid_for,p,relaxation,(values=vals,),built.table,nothing)
        @test N.expectation(result,symbol(a))==0
        @test N.expectation(result,symbol(a)^2)==1
        @test N.expectation(result,symbol(a)*symbol(b))==1
        @test_throws r"scalar-factor decomposition.*expectation" N.gns(result)
        for options in ((strengthen=(N.DensityMatrix(1:1),),),(strengthen=(N.InvariantCone(1:3,"YII"),),),
                        (strengthen=(N.Stationarity(),),),(strengthen=(N.StateOptimality(),),),
                        (strengthen=(N.Sector(a,1),),),(strengthen=(N.FermionicConditions(1:1),),))
            @test_throws r"strengthening.*remove strengthen.*Zero/Positive" N.build(p,N.Relaxation(;order=2,options...))
        end
        for constructor in (N.ExpectationZero,N.ExpectationPositive)
            @test_throws r"outer averaged functional conditions are not supported" constructor(symbol(a))
        end
        @test_throws r"continuous symmetry.*remove continuous generators" N.Problem(symbol(a);symmetry=(N.SpinRotation(0,(UInt64(1),)),))
        # Closing before row multiplication would turn this unit entry into <a>^2.
        W=N.word_type(r,4); wa=N.word(r,:s1_1,W)
        program,table=N.plain_program([one(W),wa],[(1.0,one(W))];C=Float64)
        @test program.forms.ids[last(N.term_range(program.forms,3))]==1
        @test length(table)==2

        # A central equality is localized, not only one equation on its mean.
        _,(q,)=N.site_variables(N.Free(),1,1;functional=N.RealState())
        s=N.state_symbol(q[1])
        constrained=N.Problem(-s^2;constraints=(N.Zero(s),N.Positive(1-q[1]^2)))
        constrained_build=N.build(constrained,N.Relaxation(order=2))
        @test constrained_build.sizes.localizing_rows==[3,3]
        solved=N.solve(constrained,N.Relaxation(order=2);solver=backend)
        @test solved.bound ≈ 0 atol=1e-6
        @test N.expectation(solved,s^2) ≈ 0 atol=1e-6
        @test_throws r"scalar-factor decomposition" N.gns(solved)
        positive=N.Problem(s;constraints=(N.Positive(s-0.25),N.Positive(1-q[1]^2)))
        @test N.solve(positive,N.Relaxation(order=2);solver=backend).bound ≈ 0.25 atol=1e-6

        _,(z,)=N.site_variables(N.Unipotent(),1,4;functional=N.RealTrace())
        chsh=-N.trace_symbol(z[1]*z[3]+z[1]*z[4]+z[2]*z[3]-z[2]*z[4])
        solved=N.solve(N.Problem(chsh),N.Relaxation(order=1);solver=backend)
        @test solved.bound ≈ -2sqrt(2) atol=1e-5
        @test N.expectation(solved,chsh) ≈ solved.bound atol=1e-6
    end
    if references
        @testset "Functional acceptance references" begin
            for trace in (false,true),signs in (DOC_SIGNS,BENCH_SIGNS)
                _,h=covariance_model(trace;signs)
                result=N.solve(N.Problem(h),N.Relaxation(order=2);solver=backend)
                @test result.bound ≈ -5 atol=1e-3
                old=cs_nctssos(old_covariance(trace;signs),SolverConfig(optimizer=factory,order=2,
                    cs_algo=NoElimination(),ts_algo=NoElimination());dualize=false)
                @test result.bound ≈ old.objective atol=1e-6
                @info "Functional covariance" trace signs result.bound legacy=old.objective difference=result.bound-old.objective
                if trace
                    @test result.sizes.rows==115
                    @test result.sizes.moments==1010
                end
            end
            _,(x,)=N.site_variables(N.Projector(),1,3;functional=N.RealTrace())
            h=N.trace_symbol(x[1]*x[2]*x[3])+N.trace_symbol(x[1]*x[2])*N.trace_symbol(x[3])
            for (order,expected) in ((2,-0.046717378455438933),(3,-0.03124998978001017))
                result=N.solve(N.Problem(h),N.Relaxation(;order);solver=backend)
                @test result.bound ≈ expected atol=1e-6
                @info "Projector trace" order result.bound
            end
        end
    end
end

function check_blocks(backend)
    @testset "Functional sparse block route" begin
        for F in (N.RealState(),N.RealTrace())
            symbol=F isa N.RealState ? N.state_symbol : N.trace_symbol
            r,ops=N.site_variables(N.Unipotent(),F isa N.RealState ? 2 : 1,F isa N.RealState ? 2 : 4;functional=F)
            a,b=F isa N.RealState ? ops : (ops[1][1:2],ops[1][3:4])
            h=-symbol(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])
            p=N.Problem(h)
            plain=N.build(p,N.Relaxation(order=1))
            block=N.build(p,N.BlockRelaxation(order=1))
            @test N.block_sides(plain)==N.block_sides(block)==[9]
            @test plain.sizes.real_variables==block.sizes.real_variables==20
            @test all(c->c.kind==N.PSD,block.program.blocks)
            @test N.solve(block;solver=backend).bound ≈ -2sqrt(2) atol=1e-6
            rows=collect(N.basis(r,1))
            split=N.TermSplit(algorithm=MMD())
            automatic=N.build(p,N.Relaxation(order=1,split=split))
            starting=N.build(p,N.Relaxation(rows=rows,split=split))
            f,g=N.moment_standard_form(automatic.program),N.moment_standard_form(starting.program)
            @test (f.A,f.b,f.c,f.offset,f.cones)==(g.A,g.b,g.c,g.offset,g.cones)
            custom=N.build(p,N.Relaxation(family=[rows]))
            @test N.solve(custom;solver=backend).bound ≈ -2sqrt(2) atol=1e-6

            for rule in (N.Free(),N.Projector(),N.Unipotent())
                registry,(v,)=N.site_variables(rule,1,3;functional=F)
                a,b,c=v
                objective=-symbol(a*b)^2+0.25symbol(a)^2-0.125symbol(c)^2
                constraints=rule isa N.Free ? Tuple(N.Positive(1-q^2) for q in v) : ()
                problem=N.Problem(objective;constraints)
                dense=N.solve(problem,N.Relaxation(order=2);solver=backend).bound
                for grouping in (nothing,N.Correlative())
                    previous=-Inf
                    for level in 1:3
                        built=N.build(problem,N.Relaxation(;order=2,grouping,split=N.TermSplit(;level)))
                        bound=N.solve(built;solver=backend).bound
                        @test bound<=dense+1e-6
                        @test bound>=previous-1e-6
                        previous=bound
                    end
                    @test previous ≈ dense atol=1e-6
                end
                grouped=N.build(problem,N.Relaxation(order=2,grouping=N.Correlative()))
                dependencies=Set(Set(N.operator_name(registry,l) for w in rows for l in N.support_letters(w))
                    for rows in Iterators.take(N.block_rows(grouped),2))
                @test dependencies==Set((Set((:s1_1,:s1_2)),Set((:s1_3,))))
            end

            registry,(v,)=N.site_variables(N.Free(),1,3;functional=F)
            a,b,c=v
            constrained=N.Problem(symbol(a)^2+symbol(c)^2;constraints=(N.Zero(symbol(a)-symbol(c)),))
            grouped=N.build(constrained,N.Relaxation(order=1,grouping=N.Correlative()))
            dependencies=Set(Set(N.operator_name(registry,l) for w in rows for l in N.support_letters(w))
                for rows in Iterators.take(N.block_rows(grouped),2))
            @test dependencies==Set((Set((:s1_1,:s1_3)),Set((:s1_2,))))

            _,(q,)=N.site_variables(N.Free(),1,1;functional=F)
            s=symbol(q[1])
            for (problem,expected) in ((N.Problem(-s^2;constraints=(N.Zero(s),N.Positive(1-q[1]^2))),0.0),
                (N.Problem(s;constraints=(N.Positive(s-0.25),N.Positive(1-q[1]^2))),0.25))
                for options in ((grouping=N.Correlative(),),(split=N.TermSplit(),),
                    (grouping=N.Correlative(),split=N.TermSplit()))
                    @test N.solve(problem,N.Relaxation(;order=2,options...);solver=backend).bound ≈ expected atol=1e-6
                end
            end

            _,(u,)=N.site_variables(N.Unipotent(),1,1;functional=F)
            a=only(u); s=symbol(a); rows=[one(a),a,s]
            problem=N.Problem(-s^2)
            # Normalized traces are also states. The balanced 2x2 involution
            # has inner mean zero. The ensemble of scalar +1/-1 has outer
            # mean zero but outer second moment one. Both have <a^2>=1.
            for matrices in ([Float64[1 0;0 -1]], [fill(1.0,1,1),fill(-1.0,1,1)])
                means=[tr(m)/size(m,1) for m in matrices]
                first_moment=sum(means)/length(means)
                second_moment=sum(abs2,means)/length(means)
                @test first_moment==0
                for relaxation in (N.Relaxation(family=[rows]),N.Relaxation(rows=rows,split=N.TermSplit()))
                    built=N.build(problem,relaxation)
                    values=[N.degree(label.word)==2 ? second_moment : first_moment for label in built.table.labels[2:end]]
                    result=N.SolveResult(-second_moment,:test,built.sizes,built.valid_for,problem,relaxation,(;values),built.table,nothing)
                    @test N.expectation(result,s)==first_moment
                    @test N.expectation(result,s^2)==second_moment
                    @test N.expectation(result,a^2)==1
                    y=[1.0;values]; forms=built.program.forms
                    for block in built.program.blocks
                        matrix=zeros(block.side,block.side)
                        for j in 1:block.side,i in j:block.side
                            value=sum(real(forms.coefs[t])*y[forms.ids[t]] for t in
                                N.term_range(forms,block.first_form+N.tri_offset(block.side,i,j));init=0.0)
                            matrix[i,j]=matrix[j,i]=value
                        end
                        @test eigmin(Symmetric(matrix))>=-1e-12
                    end
                end
            end
        end
    end
end
end
