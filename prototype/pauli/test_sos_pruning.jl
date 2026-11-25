using Test
import NCTSSoS as NC

if !isdefined(@__MODULE__, :PauliWorkflowRedesign)
    include("PauliWorkflowRedesign.jl")
end
import .PauliWorkflowRedesign as WF

@testset "Solver-free ambient SOS pruning" begin
    _, axes = NC.pauli_variables(2)
    x, y, z = map(axis -> WF._wf_exact.(axis), axes)
    A, B = x[1]*x[2], y[1]*x[2]
    a, b = 99//101, 20//101
    # A and B anticommute, and a²+b²=1. Left multiplication by Z₁
    # changes the factor phases but not f†f. The exact ground energy is -4.
    H = -4*(a*A + b*B)
    f = (1+im)*z[1]*(one(A)-a*A-b*B)
    residual = 2b^2*one(A)-4b*B
    cost = 2b^2 + 4b

    rejected = WF.prune_sos_factors(H,-4; squares=(f,),loss_budget=cost-1//10^12)
    @test (rejected.removed,rejected.rejected) == (0,3)
    @test only(rejected.squares) == f
    @test rejected.certificate.bound == rejected.baseline_bound == -4
    @test rejected.certificate.residual_l1 == rejected.bound_loss == 0

    pruned = WF.prune_sos_factors(H,-4; squares=(f,),loss_budget=cost)
    @test (pruned.removed,pruned.rejected) == (1,2)
    @test only(pruned.squares) == (1+im)*z[1]*(one(A)-a*A)
    @test pruned.certificate.residual == Tuple(NC.terms(residual))
    @test pruned.certificate.residual_l1 == pruned.bound_loss == cost
    @test pruned.certificate.bound == -4-cost
    replay = WF.certify_rational_sos(H,-4; squares=pruned.squares)
    @test (replay.bound,replay.residual) == (pruned.certificate.bound,pruned.certificate.residual)

    # The same trial cancels a coefficient in a nonzero starting residual.
    improved = WF.prune_sos_factors(-4a*A,-4; squares=(f,),loss_budget=0)
    @test improved.removed == 1 && improved.bound_loss == 0
    @test improved.baseline_bound == -4-4b
    @test improved.certificate.bound == -4-2b^2

    # Each complete-factor removal costs b², but the budget is shared.
    cumulative = WF.prune_sos_factors(2b^2*one(A),0;
        squares=(b*x[1],b*z[2]),loss_budget=b^2)
    @test (cumulative.removed,cumulative.rejected) == (1,1)
    @test length(cumulative.squares) == 1
    @test cumulative.certificate.bound == -b^2

    # On the conserved ∑Z=0 sector, H's exact energy is b². The multiplier
    # must still be present after removing the sole square monomial.
    sector_H = z[1]+z[2]+b^2
    sector = (; target_mz_sector=0, sector_multipliers=(one(A),))
    pruned = WF.prune_sos_factors(sector_H,0; squares=(b*x[1],),loss_budget=b^2,sector...)
    @test isempty(pruned.squares) && pruned.baseline_bound == 0
    replay = WF.certify_rational_sos(sector_H,0; squares=pruned.squares,sector...)
    @test replay.bound == pruned.certificate.bound == -b^2
    @test replay.residual == Tuple(NC.terms(b^2*one(A)))
    @test_throws ArgumentError WF.prune_sos_factors(H,-4; squares=(f,),loss_budget=-1)
end
