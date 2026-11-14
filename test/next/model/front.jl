# Import only the declarative front in a fresh namespace, as a user does.
module FrontDoorChecks
using Test
using NCTSSoS.Next
import NCTSSoS: Next as N

function same_program(old_problem,old_relaxation,problem,relaxation)
    old = old_relaxation isa N.BlockRelaxation ?
        invoke(N.build,Tuple{N.Problem,N.BlockRelaxation},old_problem,old_relaxation) :
        N.build(old_problem,old_relaxation)
    new = N.build(problem,relaxation)
    a,b=N.moment_standard_form(old.program),N.moment_standard_form(new.program)
    for field in (:A,:b,:c,:offset,:cones)
        @test getproperty(a,field)==getproperty(b,field)
    end
    return new
end

function check(backend)
    @testset "One declarative front" begin
        registry,(x,y,z)=pauli_variables(4)
        h=sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+0.7z[i]*z[mod1(i+1,4)] for i in 1:4)
        symmetry=Symmetry(SitePermutation([2,3,4,1]))
        p=Problem(h)
        symmetric=Problem(h;symmetry)
        same_program(p,Relaxation(order=1),Problem(h;symmetry=Symmetry()),Relaxation(order=1))
        for options in ((;order=1),(;order=2,grouping=Correlative(),split=TermSplit()),
            (;order=1,strengthen=(DensityMatrix(1:3),Stationarity())),
            (;order=2,row_rule=Contiguous(periodic=true)))
            same_program(p,N.BlockRelaxation(;options...,symmetry),symmetric,Relaxation(;options...))
        end
        rows=[one(x[1]);x;y;z]
        same_program(p,Relaxation(rows=rows),Problem(h;symmetry=Symmetry()),Relaxation(rows=rows))
        same_program(p,N.BlockRelaxation(family=[rows],symmetry=symmetry),symmetric,Relaxation(rows=rows))
        same_program(p,N.BlockRelaxation(family=[rows],symmetry=symmetry),symmetric,Relaxation(family=[rows]))
        _,(a,b)=site_variables(Unipotent(),2,2)
        bell=Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2]))
        real_form=same_program(bell,Relaxation(order=1),Problem(bell.objective;symmetry=Symmetry()),Relaxation(order=1))
        @test eltype(real_form.program.forms.coefs)===Float64
        _,(c,cd)=fermionic_variables(2)
        fermion=Problem(-(cd[1]*c[2]+cd[2]*c[1]))
        sector=Sector(cd[1]*c[1]+cd[2]*c[2],1)
        same_program(fermion,N.BlockRelaxation(order=1,strengthen=(sector,)),fermion,
            Relaxation(order=1,strengthen=sector))

        # Empty declarations are neutral in either position, independent of provenance.
        for empty in (Symmetry(),N.FamilySymmetry(;detected=true))
            same_program(p,N.BlockRelaxation(order=1,symmetry=symmetry),
                Problem(h;symmetry=empty),N.BlockRelaxation(order=1,symmetry=symmetry))
            same_program(p,N.BlockRelaxation(order=1,symmetry=symmetry),
                symmetric,N.BlockRelaxation(order=1,symmetry=empty))
            same_program(p,Relaxation(order=1),Problem(h;symmetry=empty),Relaxation(order=1))
        end
        @test_throws r"both.*keep it on Problem" N.build(symmetric,N.BlockRelaxation(order=1,symmetry=symmetry))
        @test_throws r"Declare symmetry on Problem.*remove" Relaxation(order=2,symmetry=symmetry)
        @test_throws r"rows.*grouping.*remove" Relaxation(rows=rows,grouping=Correlative())
        @test_throws r"exactly one.*rows.*family" Relaxation(rows=rows,family=[rows])
        foreign=N.PauliWord{2}
        wrong=CliffordMap([N.pauli_letter(foreign,1,2),N.pauli_letter(foreign,1,0)])
        @test_throws r"generator 1.*word type.*construct" Problem(h;symmetry=Symmetry(wrong))
        @test_throws r"generator 1.*algebra.*construct" Problem(fermion.objective;symmetry=Symmetry(wrong))
        @test_throws r"generator 1.*outside sites.*rebuild" Problem(h;symmetry=Symmetry(SitePermutation([2,3,4,5,1])))
        constrained=Problem(h;constraints=(Positive(z[1]),),symmetry)
        @test_throws r"generator 1.*constraint 1.*Z1.*remove" N.build(constrained,Relaxation(order=1))
        @test_throws r"generator 1.*constraint 1.*remove" N.build(constrained,N.BlockRelaxation(order=1))
        @test_throws r"generator 1.*objective.*remove" N.build(Problem(z[1];symmetry),Relaxation(order=1))
        @test N.add_constraint(symmetric,Positive(sum(z))).symmetry === symmetry
        @test detect_symmetry(zero(h)).detected
        @test isempty(detect_symmetry(zero(h)).generators)
        same_program(p,N.BlockRelaxation(order=1,strengthen=(DensityMatrix(1:3),)),p,
            Relaxation(order=1,strengthen=DensityMatrix(1:3)))
        report=inspect(symmetric,Relaxation(order=2,grouping=Correlative(),split=TermSplit()))
        shown=sprint(show,report)
        @test all(occursin(text,shown) for text in ("block generators","order=2","grouping=","split=","Sizes:"))

        # Exercise the public spelling and the single result type with a real solve.
        n=6
        _,(x,y,z)=pauli_variables(n)
        h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.7z[i]*z[mod1(i+1,n)] for i in 1:n)
        symmetry=Symmetry(SitePermutation([2:n;1]),SitePermutation([mod1(2-i,n) for i in 1:n]))
        options=(;order=2,grouping=Correlative(),split=TermSplit(),strengthen=(DensityMatrix(1:3),Stationarity()))
        old=solve(Problem(h),N.BlockRelaxation(;options...,symmetry);solver=backend)
        result=solve(Problem(h;symmetry),Relaxation(;options...);solver=backend)
        @test result isa N.SolveResult
        @test result.bound ≈ old.bound atol=1e-6
        @test expectation(result,h) ≈ result.bound atol=1e-6
        @test isfinite(expectation(result,x[1]*x[2]))
        @test occursin("strengthen=",sprint(show,result))
        @info "Front-door XXZ6" backend=typeof(backend) old=old.bound new=result.bound

        sparse=solve(Problem(h;symmetry),Relaxation(order=2,grouping=Correlative(),split=TermSplit());solver=backend)
        @test_throws r"Moment Z2.*selected blocks.*no split.*site group" expectation(sparse,z[2])
        _,(_,_,z)=pauli_variables(2)
        p=Problem(sum(z);symmetry=Symmetry(SitePermutation([2,1])))
        initial=solve(p,Relaxation(order=1);solver=backend)
        appended=solve(initial,Positive(4+sum(z)))
        @test appended.problem.symmetry === p.symmetry
        @test appended.bound ≈ initial.bound atol=1e-6
        @test_throws r"generator 1.*constraint 1.*remove" solve(initial,Positive(z[1]))
    end
end
end
