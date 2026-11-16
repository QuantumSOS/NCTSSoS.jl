module LocalReferenceChecks
using Test, LinearAlgebra, TOML
using NCTSSoS
import JuMP

# Local-only legacy cases. Bell equations/builds are routine; every solve
# requires references=true and the explicit native Mosek gate below.
# The old real moment functional identified f and f'; retain that objective
# explicitly when the old expression is not Hermitian as an operator.
hermitian_part(f)=(f+f')/2

function timed_solve(id,problem,relaxation,solver)
    println("START_LOCAL_REFERENCE ",id); flush(stdout)
    measured=@timed solve(problem,relaxation;solver)
    result=measured.value
    @info "Local reference" id seconds=measured.time bytes=measured.bytes bound=result.bound status=result.status
    println("LOCAL_DIAGNOSTICS ",id," ",result.solution.diagnostics," raw_status=",result.solution.raw_status)
    println("LOCAL_SIZES ",id," ",result.sizes)
    flush(stdout); flush(stderr)
    return result
end

function motzkin(sos_solver)
    @testset "T4 operator-only trace SOHS failure" begin
        registry,(v,)=site_variables(Free(),1,2;functional=NCTSSoS.RealTrace())
        x,y=v
        f=x*y^4*x+y*x^4*y-3x*y^2*x+1
        problem=Problem(trace_symbol(f))
        # Operator rows, not products of scalar trace symbols.
        rows=[one(x)]
        layer=[one(x)]
        for _ in 1:3
            layer=[w*v for w in layer for v in (x,y)]
            append!(rows,layer)
        end
        result=timed_solve("T4",problem,Relaxation(rows=rows),sos_solver)
        @test isnothing(result.bound)
        @test result.bound_kind===:unavailable
        @test result.status==JuMP.MOI.OTHER_ERROR
    end
end

function heisenberg(solver)
    @testset "Local Heisenberg rings and torus" begin
        _,ops=pauli_variables(6)
        h=sum(op[i]*op[mod1(i+1,6)] for op in ops for i in 1:6)/4
        result=timed_solve("XXX6",Problem(h),Relaxation(order=2),solver)
        @test result.bound/6 ≈ -0.467129 atol=1e-6
        for (n,couplings) in ((4,0.1:0.1:2.1),(6,(0.2,)))
            _,ops=pauli_variables(n)
            for j2 in couplings
                h=sum((op[i]*op[mod1(i+1,n)]+j2*op[i]*op[mod1(i+2,n)])/4 for op in ops for i in 1:n)
                problem=Problem(h)
                lower=timed_solve("J1J2-$n-$j2-level1",problem,Relaxation(order=2,split=TermSplit(algorithm=MMD(),level=1)),solver)
                higher=timed_solve("J1J2-$n-$j2-level2",problem,Relaxation(order=2,split=TermSplit(algorithm=MMD(),level=2)),solver)
                @test lower.bound!==nothing && higher.bound!==nothing
                @test higher.bound>=lower.bound-1e-6
                if n==4
                    @test higher.bound ≈ min(-2+j2,-3j2) atol=1e-6
                else
                    @test higher.bound/n ≈ -0.4270083225302217 atol=1e-6
                end
            end
        end
        _,ops=pauli_variables(9)
        index=LinearIndices((1:3,1:3))
        h=sum((op[index[i,j]]*op[index[i,mod1(j+1,3)]]+
               op[index[i,j]]*op[index[mod1(i+1,3),j]])/4 for op in ops for i in 1:3 for j in 1:3)
        previous=-Inf
        for level in 1:5
            result=timed_solve("torus3x3-level$level",Problem(h),
                Relaxation(order=3,grouping=Correlative(algorithm=MF()),split=TermSplit(algorithm=MMD(),level=level)),solver)
            @test result.bound!==nothing
            @test -4.5-1e-6<=result.bound<=9*(-0.44100019443650207)+1e-6
            @test result.bound>=previous-1e-6
            previous=result.bound
        end
    end
end

function benchmark_problem(kind)
    n=kind==:banded ? 4 : kind in (:singular,:wood) ? 8 : 6
    _,(x,)=site_variables(Free(),1,n)
    f=zero(x[1])
    if kind==:rosenbrock
        f+=n
        for i in 2:n
            f+=100x[i-1]^4-200x[i-1]^2*x[i]-2x[i]+101x[i]^2
        end
    elseif kind==:banded
        for i in 1:n
            neighbors=setdiff(max(1,i-5):min(n,i+1),i)
            residual=2x[i]+5x[i]^3+1-sum(x[j]+x[j]^2 for j in neighbors)
            f+=residual'*residual
        end
    elseif kind==:tridiagonal
        for i in 1:n
            residual=3x[i]-2x[i]^2+1
            i>1 && (residual-=x[i-1])
            i<n && (residual-=2x[i+1])
            f+=residual'*residual
        end
    elseif kind==:singular
        # This file uses NC fourth powers, unlike the E2 chip benchmark.
        for i in 1:2:n-3
            f+=(x[i]+10x[i+1])^2+5(x[i+2]-x[i+3])^2+
                (x[i+1]-2x[i+2])^4+10(x[i]-10x[i+3])^4
        end
    else
        f+=21n-41
        for i in 1:2:n-3
            f+=-2x[i]+x[i]^2+100x[i]^4-200x[i]^2*x[i+1]-40x[i+1]+
                110.1x[i+1]^2+19.8x[i+1]*x[i+3]-2x[i+2]+x[i+2]^2+
                90x[i+2]^4-180x[i+2]^2*x[i+3]-40x[i+3]+100.1x[i+3]^2
        end
    end
    return Problem(hermitian_part(f))
end

function benchmarks(solver,sos_solver,chained_solver)
    # Exact old Mosek references; keep the old 1e-6 tolerance.
    expected=(rosenbrock=(0.999999995930163,0.999999973478842,0.9999997821660428),
        banded=(-2.0628417282628906e-8,-2.0628417282628906e-8,-1.7480268537791176e-9),
        tridiagonal=(6.580588236100922e-11,-7.683336739295744e-9,-1.3851129121677345e-7),
        singular=(-6.428097254420157e-8,-1.2251940085171628e-7,-2.4177846987055753e-7),
        wood=(0.9999998013387236,0.9999998300594172,0.9999998815901131))
    # (side => multiplicity, moments). Dense/CS records match the old file.
    # TS records use NCTSSoS's row order for chordal fill and remove duplicate
    # blocks across groups; the old moment counts remain unchanged.
    shapes=(rosenbrock=(([43=>1],820),([7=>5],90),([1=>5,2=>14,3=>5],33)),
        banded=(([85=>1],2815),([85=>1],2815),([1=>37,2=>1,3=>17,4=>17,5=>1,7=>1,9=>2],134)),
        tridiagonal=(([43=>1],820),([13=>4],226),([2=>8,3=>10,4=>8,5=>4],62)),
        singular=(([73=>1],2413),([13=>6],328),([1=>11,2=>17,3=>6,4=>6],83)),
        wood=(([73=>1],2413),([7=>7],124),([1=>4,2=>18,3=>7],46)))
    @testset "Local NC polynomial benchmarks" begin
        for kind in keys(expected), route in 1:3
            problem=benchmark_problem(kind)
            relaxation=Relaxation(order=kind==:banded ? 3 : 2,
                grouping=route==1 ? nothing : Correlative(algorithm=MF()),
                split=route==3 ? TermSplit(algorithm=MMD(),level=1) : nothing)
            built=NCTSSoS.build(problem,relaxation)
            sides,moments=shapes[kind][route]
            @test sort(NCTSSoS.block_sides(built))==reduce(vcat,[fill(side,count) for (side,count) in sides])
            @test built.sizes.moments==moments
            # Old tests use SOS except Chained Singular CS+TS, which uses
            # moments and one-decade looser interior-point tolerances.
            backend=kind==:singular ? chained_solver(route==3 ? :moment : :sos) : sos_solver
            result=timed_solve("$kind-route$route",problem,relaxation,backend)
            if kind==:singular && route==3
                # Same original 1e-7 Mosek settings: old -9.9092722e-7,
                # native NCTSSoS -4.2915344e-6 (difference < 1e-5).
                # Exact minimum is zero. Keep the validity-side tolerance.
                @test -1e-5<=result.bound<=1e-6
            else
                @test result.bound ≈ expected[kind][route] atol=1e-6
            end
            if kind==:singular && route==3
                stopped=timed_solve("singular-route3-SOS",problem,relaxation,sos_solver)
                if stopped.status in (JuMP.MOI.SLOW_PROGRESS,JuMP.MOI.ITERATION_LIMIT)
                    @test isnothing(stopped.bound)
                else
                    @test stopped.bound!==nothing && stopped.bound<=0
                end
            end
        end
    end
end

function large_scale(solver,sos_solver)
    @testset "Local constrained Broyden n=10" begin
        n=10
        _,(x,)=site_variables(Free(),1,n)
        f=zero(x[1])
        for i in 1:n
            neighbors=setdiff(max(1,i-5):min(n,i+1),i)
            residual=2x[i]+5x[i]^3+1-sum(x[j]+x[j]^2 for j in neighbors)
            f+=residual'*residual
        end
        problem=Problem(f;constraints=Positive.(vcat([1-v^2 for v in x],[v-1/3 for v in x])))
        relaxation=Relaxation(order=3,grouping=Correlative(algorithm=MF()),split=TermSplit(algorithm=MMD(),level=1))
        moment=timed_solve("large10-moment",problem,relaxation,solver)
        @test moment.bound ≈ 3.011288353315061 atol=1e-4
        # Old 1544 counts moment blocks only. NCTSSoS counts the full program,
        # including localizers, with its row-order-dependent chordal fill.
        @test moment.sizes.blocks==2268
        @test moment.sizes.moments==5748
        sos=timed_solve("large10-SOS",problem,relaxation,sos_solver)
        @test sos.bound ≈ 3.011288353315061 atol=1e-6
    end
end

function bell_build(case)
    m,n=case["m"],case["n"]
    registry,(x,y)=site_variables(Projector(),2,max(m,n))
    equation=case["equation"]
    startswith(equation,'-') || (equation="+"*equation)
    matches=collect(eachmatch(r"([+-])(\d+)\*([XY])\[(\d+)\](?:\*([XY])\[(\d+)\])?",equation))
    @test join(match.match for match in matches)==equation
    variable(axis,index)=(axis=="X" ? x : y)[parse(Int,index)]
    h=zero(x[1])
    for match in matches
        sign,coefficient,axis,index,axis2,index2=match.captures
        term=variable(axis,index)
        isnothing(axis2) || (term*=variable(axis2,index2))
        h+=(sign=="-" ? -1 : 1)*parse(Int,coefficient)*term
    end
    # The site constructor has uniform letter counts. Remove unused labels
    # from rows to preserve the original (m,n) scenario.
    excluded=Set([Symbol(:s1_,i) for i in m+1:max(m,n)])
    union!(excluded,[Symbol(:s2_,i) for i in n+1:max(m,n)])
    rows=[p for p in NCTSSoS.basis(registry,case["order"]) if
        all(l->NCTSSoS.operator_name(registry,l) ∉ excluded,
            NCTSSoS.support_letters(last(only(NCTSSoS.terms(p)))))]
    NCTSSoS.build(Problem(h),Relaxation(rows=rows))
end

function bell(solver=nothing;references=false)
    data=TOML.parsefile(joinpath(@__DIR__,"bell_inequalities/local_bell_table.toml"))
    @testset "Local Bell equations and questioned references" begin
        for case in data["cases"]
            id=case["id"]
            @testset "$id" begin
                built=bell_build(case)
                @test built isa NCTSSoS.BuiltRelaxation
                # The old file's skip_instances is empty, despite its A17
                # investigation comment. Keep every equation and order.
                selected=case["order"]<=2 || max(case["m"],case["n"])<=4 || id in ("A9","A32")
                if !references || isnothing(solver) || !selected
                    @test false skip=true
                    continue
                end
                measured=@timed solve(built;solver)
                result=measured.value
                @info "Bell reference" id seconds=measured.time bound=result.bound reference=case["reference"] diagnostics=result.solution.diagnostics
                flush(stderr)
                @test result.bound!==nothing && isfinite(result.bound)
                if id=="A17"
                    # Old expected -0.3754473; NCTSSoS -0.37544733205129416.
                    # The old comment's -1.115 result is stale.
                    @test result.bound ≈ case["reference"] atol=1e-6
                end
                if !isnothing(result.bound) && !isapprox(result.bound,case["reference"];atol=1e-6)
                    # These questioned values are provenance, not assertions.
                    @info "Questioned Bell reference" id bound=result.bound reference=case["reference"] difference=result.bound-case["reference"]
                end
            end
        end
    end
end

function check(;references=false,mosek=false)
    licensed=references && mosek
    bell(licensed ? MosekBackend(MSK_IPAR_NUM_THREADS=8) : nothing;references=licensed)
    @testset "Mosek local-only removal references" begin
        if !licensed
            @test false skip=true
        else
            moment=MosekBackend(form=:moment,MSK_IPAR_NUM_THREADS=1)
            sos=MosekBackend(form=:sos,MSK_IPAR_NUM_THREADS=1)
            chained(form)=MosekBackend(;form,MSK_IPAR_NUM_THREADS=1,
                MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-7,MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-7,
                MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-7,MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-7)
            motzkin(sos)
            benchmarks(moment,sos,chained)
            heisenberg(sos)
            large_scale(moment,sos)
        end
    end
end
end
