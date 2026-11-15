module NewtonChipChecks
using Test, NCTSSoS
import NCTSSoS: Next as N
import JuMP

same_rows(actual,expected)=length(actual)==length(expected) && all(x->count(==(x),actual)==1,expected)

# Test the public helper's conservative behavior with unaccepted or invalid
# numerical certificates. Neither backend proves exclusion of any row.
struct UncertainBackend <: N.AbstractSolver
    claimed_optimal::Bool
end
function N.solve_standard(form::N.MomentStandardForm,backend::UncertainBackend)
    N.ConicSolution(backend.claimed_optimal ? -1.0 : nothing,
        backend.claimed_optimal ? JuMP.MOI.OPTIMAL : JuMP.MOI.ITERATION_LIMIT,
        "no separating certificate",zeros(length(form.c)),nothing)
end

function check(backend)
    @testset "Newton-chip custom rows" begin
        @test Base.ispublic(N,:newton_chip_basis)
        @test !Base.isexported(N,:newton_chip_basis)
        @test NCTSSoS.newton_chip_basis !== N.newton_chip_basis
        _,(x,)=N.site_variables(N.Free(),1,2)
        a,b=x
        p=N.Problem(1+a^2+2b*a^2*b)
        rows=N.newton_chip_basis(p,2)
        @test same_rows(N.newton_chip_basis(p),rows)
        @test same_rows(rows,[one(a),a,b,a*b])
        dense=N.solve(p,N.Relaxation(order=2);solver=backend)
        chip=N.solve(p,N.Relaxation(rows=rows);solver=backend)
        @test dense.bound ≈ 1 atol=1e-6
        @test chip.bound ≈ 1 atol=1e-6
        @test (dense.sizes.rows,chip.sizes.rows)==(7,4)
        @test N.newton_chip_basis(N.Problem(one(a)),0)==[one(a)]
        @test N.newton_chip_basis(N.Problem(zero(a)),0)==[one(a)]
        @test N.newton_chip_basis(N.Problem(one(a)))==[one(a)]
        @test N.newton_chip_basis(N.Problem(zero(a)))==[one(a)]
        @test_throws r"at least 2" N.newton_chip_basis(p,1)
        @test_throws r"nonnegative" N.newton_chip_basis(p,-1)
        @test_throws r"unconstrained" N.newton_chip_basis(N.Problem(a^2;constraints=(N.Positive(a),)),1)
        @test_throws r"symmetry" N.newton_chip_basis(N.Problem(a^2;symmetry=(N.SitePermutation([1]),)),1)
        # A non-square support term cannot supply an ordinary right chip.
        @test same_rows(N.newton_chip_basis(N.Problem(a^2+b^2+a*b+b*a+a),1),[one(a),a,b])

        _,(y,)=N.site_variables(N.Free(),1,2;functional=N.RealTrace())
        a,b=y; tr=N.trace_symbol
        p=N.Problem(tr(1+a^2+b*a^2*b))
        rows=N.newton_chip_basis(p,2;solver=backend)
        @test same_rows(N.newton_chip_basis(p;solver=backend),rows)
        @test same_rows(rows,[one(a),a,a*b,b*a])
        dense=N.solve(p,N.Relaxation(order=2);solver=backend)
        chip=N.solve(p,N.Relaxation(rows=rows);solver=backend)
        # Both nonconstant terms are traces of squares; a=0 attains one.
        @test dense.bound ≈ 1 atol=1e-6
        @test chip.bound ≈ 1 atol=1e-6
        @test (dense.sizes.rows,chip.sizes.rows)==(19,4)
        @test N.newton_chip_basis(N.Problem(tr(one(a))),0;solver=backend)==[one(a)]
        @test N.newton_chip_basis(N.Problem(tr(zero(a))),0;solver=backend)==[one(a)]
        @test_throws r"solver=" N.newton_chip_basis(p,2)
        for bad in (N.Problem(tr(a)^2),N.Problem(tr(a)*b))
            @test_throws r"scalar single-trace" N.newton_chip_basis(bad,1;solver=backend)
        end
        for claimed in (false,true)
            conservative=N.newton_chip_basis(p,2;solver=UncertainBackend(claimed))
            @test same_rows(conservative,[one(a),a,b,a*b,b*a])
            @test N.solve(p,N.Relaxation(rows=conservative);solver=backend).bound ≈ 1 atol=1e-6
        end
        # Boundary exponent (1,1) lies in half conv{(0,0),(4,0),(0,4)}.
        # A numerical LP result must not remove boundary rows.
        boundary=N.newton_chip_basis(N.Problem(tr(a^4+b^4)),2;solver=backend)
        @test same_rows(boundary,[one(a),a,b,a^2,a*b,b*a,b^2])
        # Only letters occurring in the objective enter the candidate set.
        @test same_rows(N.newton_chip_basis(N.Problem(tr(a^2)),2;solver=backend),[one(a),a])

        for functional in (nothing,N.RealTrace())
            _,(a,b)=N.site_variables(N.Free(),2,1;functional)
            h=a[1]^2+b[1]^2
            functional===nothing || (h=tr(h))
            @test_throws r"one noncommuting group" N.newton_chip_basis(N.Problem(h),1;solver=backend)
        end
        for (rule,functional) in ((N.Projector(),nothing),(N.Unipotent(),N.RealTrace()),(N.Free(),N.RealState()))
            _,(x,)=N.site_variables(rule,1,1;functional)
            @test_throws r"ordinary Free" N.newton_chip_basis(N.Problem(x[1]^2),1;solver=backend)
        end
    end
end
end
