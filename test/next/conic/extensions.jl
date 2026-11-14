# Standalone entry point: the SCS extension CI job installs its environment.
using Test, NCTSSoS, JuMP, LinearAlgebra, SparseArrays, SCS
const N = NCTSSoS.Next
const TOL = 1e-7

function extension_models()
    r,(x,y,z) = N.pauli_variables(4)
    h = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
    r,(a,ad) = N.fermionic_variables(4)
    hub = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    r,(x,y) = N.site_variables(N.Unipotent(),2,2)
    chsh = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    return (("Heisenberg4",N.Problem(h),N.Relaxation(order=1),-9.),
        ("Hubbard2",N.Problem(hub),N.Relaxation(order=2),1-sqrt(5)),
        ("CHSH",N.Problem(chsh),N.Relaxation(order=1),-2sqrt(2)))
end

@testset "SCS direct and indirect moment writers" begin
    @test Base.get_extension(NCTSSoS,:NCTSSoSSCSExt)!==nothing
    for linear_solver in (SCS.DirectSolver,SCS.IndirectSolver)
        backend = N.SCSBackend(;linear_solver,verbose=0,eps_abs=TOL,eps_rel=TOL)
        factory = optimizer_with_attributes(SCS.Optimizer,"linear_solver"=>linear_solver,"verbose"=>0,"eps_abs"=>TOL,"eps_rel"=>TOL)
        # Unequal real and imaginary objective coefficients expose HPSD sign/scaling errors.
        f = N.MomentStandardForm(sparse([2,3],[1,2],[-sqrt(2),-sqrt(2)],4,2),[1.,0,0,1],[1.,2],0.,[N.DataCone(:hpsd,2,1,4,0)])
        complex_result = N.solve_standard(f,backend)
        @test complex_result.bound ≈ -sqrt(5) atol=2e-6
        @test complex_result.values ≈ [-1,-2]/sqrt(5) atol=2e-6
        @test_throws ArgumentError N.prepare_standard(f,N.SCSBackend(linear_solver=:indirect))
        # Deliberately non-SCS cone order: PSD, nonnegative, then SOC.
        mixed = N.MomentStandardForm(sparse([1,4,5],[2,2,1],[-1.,-1,-1],7,2),
            [0.,sqrt(2),1,-2,0,3,4],[1.,1],0.,
            [N.DataCone(:psd,2,1,3,0),N.DataCone(:nonnegative,1,4,1,0),N.DataCone(:soc,3,5,3,0)])
        mixed_result = N.solve_standard(mixed,backend)
        @test mixed_result.bound ≈ 7 atol=2e-6
        @test mixed_result.values ≈ [5,2] atol=2e-6
        for (name,p,r,expected) in extension_models()
            ref = N.solve(p,r;solver=N.JuMPBackend(factory))
            sol = N.solve(p,r;solver=backend)
            @test sol.status==JuMP.MOI.OPTIMAL
            @test sol.bound ≈ expected atol=2e-5
            @test sol.bound ≈ ref.bound atol=2e-6
            f = N.moment_standard_form(N.build(p,r).program)
            @test dot(f.c,sol.solution.values)+f.offset ≈ sol.bound atol=2e-6
            @test !sol.solution.diagnostics.warm_started
            @test all(isfinite,(sol.solution.diagnostics.primal_residual,sol.solution.diagnostics.dual_residual,sol.solution.diagnostics.gap))
            @info "SCS writer" name linear_solver bound=sol.bound reference=ref.bound tolerance=TOL diagnostics=sol.solution.diagnostics
        end
        r,ops = N.site_variables(N.Free(),1,1)
        x = ops[1][1]
        sol = N.solve(N.Problem(x;constraints=(N.Positive(1-x^2),)),N.Relaxation(order=2);solver=backend)
        @test sol.bound ≈ -1 atol=2e-5
        limited = N.solve(sol,N.Positive(x-0.25))
        @test limited.bound ≈ 0.25 atol=2e-5
        @test limited.solution.handle===sol.solution.handle
        @test limited.solution.diagnostics.warm_started
        # ZERO coordinates must move ahead of HPSD coordinates in the SCS cone order.
        fixed = N.solve(limited,N.Zero(x-0.6))
        @test fixed.bound ≈ 0.6 atol=2e-5
        @test fixed.solution.diagnostics.warm_started
        again = N.solve_standard(fixed.solution.handle)
        @test again.bound ≈ fixed.bound atol=2e-6
        @test again.diagnostics.iterations==0
        failed = N.solve(fixed,N.Positive(x-0.8))
        @test failed.bound===nothing
        @test failed.status==JuMP.MOI.INFEASIBLE

        _,ops = N.site_variables(N.Free(),1,2)
        a,b = ops[1]
        ball = N.Problem(a*b+b*a;constraints=(N.Positive(1-a^2-b^2),))
        relaxation = N.Relaxation(order=2)
        for (name,cut) in ((:restricting,N.Positive(a+0.25)),(:redundant,N.Positive(2-a^2-b^2)))
            initial = N.solve(ball,relaxation;solver=backend)
            warm = N.solve(initial,cut)
            fresh = N.solve(N.add_constraint(ball,cut),relaxation;solver=backend)
            @test warm.status==fresh.status==JuMP.MOI.OPTIMAL
            @test warm.bound ≈ -1 atol=2e-6
            @test warm.bound ≈ fresh.bound atol=2e-6
            @test warm.sizes.blocks==initial.sizes.blocks+1
            @test warm.solution.diagnostics.warm_started
            # A feasible appended slack with zero dual preserves the old optimum.
            # A cut violated by the old moments need not reduce the iteration count.
            if name==:redundant
                @test warm.solution.diagnostics.iterations < fresh.solution.diagnostics.iterations
            end
            @info "Appended constraint: warm versus fresh" linear_solver name warm=warm.solution.diagnostics.iterations fresh=fresh.solution.diagnostics.iterations
        end
    end
    _,p,r,_ = last(extension_models())
    failed = N.solve(p,r;solver=N.SCSBackend(verbose=0,max_iters=1))
    @test failed.bound===nothing
    @test failed.bound_kind==:unavailable
end
