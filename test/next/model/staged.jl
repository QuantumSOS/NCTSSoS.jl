module StagedSolveChecks
using Test, NCTSSoS
import NCTSSoS: Next as N

function preparation_bytes(built,solver)
    N.prepare(built;solver)
    @allocated N.prepare(built;solver)
end

function check(solver;reconstruct=true)
    @testset "Public staged solve" begin
        for name in (:build,:BuiltRelaxation,:prepare,:prepare_standard,:solve_standard)
            @test Base.ispublic(N,name) && !Base.isexported(N,name)
        end
        _,(_,_,Z)=N.pauli_variables(2)
        H=-Z[1]-Z[2]; relaxation=N.Relaxation(order=2)
        for symmetry in (nothing,(N.SitePermutation([2,1]),))
            problem=N.Problem(H;symmetry)
            built=N.build(problem,relaxation)
            @test built isa N.BuiltRelaxation
            @test occursin("Sizes: ",sprint(show,built))
            @test isnothing(built.sizes.closure)==isnothing(symmetry)
            @test_throws r"solver" N.prepare(built)
            prepared=N.prepare(built;solver)
            @test preparation_bytes(built,solver)==0
            result=N.solve(built;solver)
            direct=N.solve(problem,relaxation;solver)
            @test result.bound==direct.bound
            @test result.bound ≈ -2 atol=1e-6
            @test result.table===built.table && result.sizes===built.sizes
            @test result.relaxation===relaxation
            @test N.expectation(result,H) ≈ result.bound atol=1e-7
            again=N.solve(built;solver)
            @test again.bound ≈ result.bound atol=1e-6
            @test N.prepare(built;solver)===prepared
            if reconstruct
                g=N.gns(result;atol=1e-5,check_atol=1e-4)
                @test g.flatness.is_flat
                @test N.expectation(g,H) ≈ -2 atol=1e-4
            end
        end
        # An appended solve may consume its native handle, but must not change
        # the program subsequently solved through the original built value.
        _,(x,)=N.site_variables(N.Free(),1,1); a=only(x)
        built=N.build(N.Problem(a;constraints=(N.Positive(1-a^2),)),N.Relaxation(order=2))
        initial=N.solve(built;solver)
        changed=N.solve(initial,N.Positive(a-.25))
        @test changed.bound ≈ .25 atol=1e-5
        @test N.solve(built;solver).bound ≈ -1 atol=1e-5
    end
end
end
