module GNSLegacyChecks
using Test, LinearAlgebra
using NCTSSoS.Next
import JuMP, COSMO

function check(solver)
    @testset "GNS physical legacy predictions" begin
        registry,(v,)=site_variables(Free(),1,2)
        x,y=v
        cases=(
            (:ball,2-x^2+x*y^2*x-y^2,(Positive(1-x^2-y^2),),3,1.0,false),
            (:polydisc,x*y*x,(Positive(1-x^2),Positive(1-y^2)),3,-1.0,false),
            (:mixed,2-x^2+x*y^2*x-y^2,(Positive(4-x^2-y^2),Zero(x*y+y*x-2)),3,-1.0,true),
            (:pironio,x*y+y*x,(Positive(-y^2+y+0.5),Zero(x^2-x)),2,-0.75,true),
            (:motzkin,x*y^4*x+y*x^4*y-3x*y^2*x+1,(Positive(1-x^4-y^4),),3,-0.160813,false))
        for (name,f,constraints,order,bound,flat) in cases
            @testset "$name" begin
                result=solve(Problem(f;constraints),Relaxation(;order);solver)
                @test result.bound ≈ bound atol=(name==:motzkin ? 1e-3 : 1e-5)
                g=gns(result)
                @test g isa Next.GNSResult
                @test g.flatness.is_flat==flat
                @test g.verdict==(flat ? :moments_reproduced : :optimizer_certified)
                @test norm(g.state) ≈ 1 atol=1e-8
                @test norm(g[x]-g[x]')<=1e-8
                @test norm(g[y]-g[y]')<=1e-8
                # Old verify_gns loops over full_basis: degree <= order,
                # not all the moments of degree <= 2order in the program.
                for row in Next.basis(registry,order)
                    @test expectation(g,row) ≈ expectation(result,row) atol=5e-5
                end
                @test expectation(g,f) ≈ result.bound atol=1e-5
                for constraint in constraints
                    matrix=evaluate(g,constraint.polynomial)
                    if constraint isa Zero
                        @test norm(matrix)<=1e-5
                    else
                        @test eigmin(Hermitian(matrix))>=-(name==:ball ? 1e-6 : 5e-5)
                    end
                end
                if name in (:ball,:polydisc)
                    @test eigmin(Hermitian(evaluate(g,f))) ≈ bound atol=1e-5
                end
            end
        end
        # The old unconstrained polynomial is unbounded, not a fake optimizer.
        result=solve(Problem(2-x^2+x*y^2*x-y^2),Relaxation(order=2);solver)
        @test isnothing(result.bound)
        @test result.status in (JuMP.MOI.DUAL_INFEASIBLE,JuMP.MOI.INFEASIBLE,
            JuMP.MOI.INFEASIBLE_OR_UNBOUNDED,JuMP.MOI.SLOW_PROGRESS,JuMP.MOI.ITERATION_LIMIT)
        @test_throws ArgumentError gns(result)

        _,(z,)=site_variables(Free(),1,1)
        result=solve(Problem(z[1]^2+1),Relaxation(order=1);solver)
        @test expectation(result,one(z[1])) ≈ 1 atol=1e-8
    end
end

function noise()
    @testset "Default COSMO GNS noise contract" begin
        _,(v,)=site_variables(Free(),1,2)
        x,y=v
        problem=Problem(x*y*x;constraints=(Positive(1-x^2),Positive(1-y^2)))
        solver=JuMPBackend(JuMP.optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false))
        result=solve(problem,Relaxation(order=3);solver)
        # Old raw-Hankel PSD allowance is 1e-4 times its norm; the old
        # state-normalization check also uses 1e-4 for this solver setting.
        g=gns(result;psd_rtol=1e-4,check_atol=1e-4)
        @test g.flatness.minimum_eigenvalue>=-g.flatness.psd_tolerance
        @test norm(g.state)^2 ≈ 1 atol=1e-4
    end
end

function vacuum()
    @testset "Exact vacuum creation adjoints" begin
        _,(a,ad)=fermionic_variables(1)
        problem=Problem(ad[1]*a[1])
        relaxation=Relaxation(order=2)
        built=Next.build(problem,relaxation)
        # Normal-ordered vacuum moments vanish except for the identity.
        solution=Next.ConicSolution(0.0,:OPTIMAL,"exact vacuum fixture",
            zeros(length(built.table.labels)-1),nothing)
        result=Next.SolveResult(0.0,:OPTIMAL,built.sizes,built.valid_for,
            problem,relaxation,solution,built.table,nothing)
        g=gns(result;atol=1e-10)
        A,D=g[a[1]],g[ad[1]]
        @test D ≈ A' atol=1e-12
        @test !isapprox(D,A;atol=1e-12)
        @test A^2 ≈ zeros(size(A)) atol=1e-12
        @test D^2 ≈ zeros(size(D)) atol=1e-12
    end
end
end
