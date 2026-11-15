module BosonModelChecks
using Test, LinearAlgebra, NCTSSoS
import NCTSSoS: Next as N

bh(b,bd;mu=0.5)=-(bd[1]*b[2]+bd[2]*b[1])+sum(bd[i]^2*b[i]^2-mu*bd[i]*b[i] for i in 1:2)

function check(backend,factory)
    @testset "Bosonic solve and workflow guards" begin
        registry,(b,bd)=N.bosonic_variables(2)
        old_registry,(old_b,old_bd)=create_bosonic_variables(1:2)
        for order in (2,3)
            displaced=bd[1]*b[1]+0.7*(b[1]+bd[1])+bd[2]*b[2]
            r=N.solve(N.Problem(displaced),N.Relaxation(;order);solver=backend)
            @test r.bound ≈ -0.49 atol=1e-6
            @test N.expectation(r,displaced) ≈ r.bound atol=1e-6
            @test_throws r"No finite matrices.*expectation\(result, p\)" N.gns(r)
            coupled=bd[1]*b[1]+bd[2]*b[2]+0.3*(bd[1]*b[2]+bd[2]*b[1])
            r=N.solve(N.Problem(coupled),N.Relaxation(;order);solver=backend)
            @test r.bound ≈ 0 atol=1e-6
            h=bh(b,bd)
            r=N.solve(N.Problem(h),N.Relaxation(;order);solver=backend)
            old=cs_nctssos(polyopt(bh(old_b,old_bd),old_registry),SolverConfig(
                optimizer=factory,order=order,cs_algo=NoElimination(),ts_algo=NoElimination());dualize=false)
            @test r.bound ≈ old.objective atol=1e-6
            # Independent finite occupation truncation is a variational upper bound.
            a=diagm(1=>sqrt.(1.0:6)); id=Matrix{Float64}(I,7,7)
            matrices=[kron(a,id),kron(id,a)]
            exact=eigmin(Hermitian(bh(matrices,adjoint.(matrices))))
            @test r.bound<=exact+1e-6
            @info "Bosonic BH without sector" order new=r.bound legacy=old.objective truncated=exact
        end
        h=bh(b,bd)
        dense=N.solve(N.Problem(h),N.Relaxation(order=2);solver=backend)
        symmetric=N.solve(N.Problem(h;symmetry=(N.SitePermutation([2,1]),)),N.Relaxation(order=2);solver=backend)
        @test symmetric.bound ≈ dense.bound atol=1e-6
        for (options,cs,ts) in (((grouping=N.Correlative(),),MF(),NoElimination()),
                                ((split=N.TermSplit(),),NoElimination(),MaximalElimination()))
            current=N.solve(N.Problem(h),N.Relaxation(;order=2,options...);solver=backend)
            old=cs_nctssos(polyopt(bh(old_b,old_bd),old_registry),SolverConfig(
                optimizer=factory,order=2,cs_algo=cs,ts_algo=ts);dualize=false)
            @test current.bound ≈ old.objective atol=1e-6
        end
        for kind in (N.Stationarity(),N.StateOptimality(),N.DensityMatrix(1:2),N.FermionicConditions(1:2))
            @test_throws ArgumentError N.build(N.Problem(h),N.Relaxation(order=2,strengthen=kind))
        end
        h=bh(b,bd;mu=0)
        for order in (2,3)
            plain=N.build(N.Problem(h),N.Relaxation(;order))
            rotated=N.build(N.Problem(h;symmetry=(N.NumberRotation(registry),)),N.Relaxation(;order))
            @test rotated.sizes.real_variables < plain.sizes.real_variables
            relaxation=N.Relaxation(;order,strengthen=N.Sector(sum(bd.*b),2))
            sector=N.solve(N.Problem(h),relaxation;solver=backend)
            number=N.solve(N.Problem(h;symmetry=(N.NumberRotation(registry),)),relaxation;solver=backend)
            @test sector.bound ≈ 1-sqrt(5) atol=1e-6
            @test number.bound ≈ sector.bound atol=1e-6
            @test number.bound ≈ 1-sqrt(5) atol=1e-6
            # Continuous symmetry skips sector elimination; follow up with elimination before continuous reduction.
            @test_broken number.sizes.real_variables <= sector.sizes.real_variables
        end
    end
end
end
