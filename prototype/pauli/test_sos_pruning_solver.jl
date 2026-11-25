# Run with Clarabel available, as for test_workflow.jl.
using Test
using Clarabel
import NCTSSoS as NC

if !isdefined(@__MODULE__, :PauliWorkflowRedesign)
    include("PauliWorkflowRedesign.jl")
end
import .PauliWorkflowRedesign as WF

@testset "Ambient SOS solver lift and exact replay" begin
    solver = NC.JuMP.optimizer_with_attributes(Clarabel.Optimizer,"verbose"=>false,
        "chordal_decomposition_enable"=>false)
    _, (x,y,z) = NC.pauli_variables(1)
    H = -x[1] + 0.5y[1] + 0.75z[1]
    no_symmetry = WF.PauliSymmetrySpec(time_reversal_gauge=false,anti_sign_realify=false,
        klein4_sign=false,xyz_permutation=false,rotated_lie_detection=false,
        cartan_u1_weyl=false,su2_mz0_tree=false,spatial_mode=:none)
    ambient = WF.PauliSparsitySpec(ts_algo=:None,rewrite_complete_local_rdm=false)
    budget = 1//100_000
    for form in (:dualized,:primalized)
        plan = WF.compile_relaxation(H; order=1,symmetry=no_symmetry,sparsity=ambient,moi_form=form)
        @test only(plan.program.blocks).side == 2length(only(plan.reports.original_rows))
        result = WF.solve_ground_state(plan,H; optimizer=solver,prune_budget=budget)
        @test result.bound ≈ -sqrt(29)/4 atol=1e-6
        @test !result.provenance.rigorous_lower_bound
        @test result.certificate_status === :verified_rational
        @test result.certified_bound == result.certificate.bound
        @test result.certified_bound < 0 && result.certified_bound^2 >= 29//16
        @test Float64(result.certified_bound) ≈ -sqrt(29)/4 atol=2e-5
        @test result.pruning.baseline_bound-result.certified_bound <= budget
        replay = WF.certify_rational_sos(H,result.certificate.level; squares=result.pruning.squares)
        @test replay.bound == result.certified_bound && replay.residual == result.certificate.residual
    end

    # The local-density substitution has no automatic adjoint in this hook.
    density = WF.PauliSparsitySpec(ts_algo=:None)
    reduced = WF.compile_relaxation(H; order=1,symmetry=no_symmetry,sparsity=density)
    @test !isempty(reduced.reports.support.density_rewrite.applied)
    result = WF.solve_ground_state(reduced,H; optimizer=solver,prune_budget=budget)
    @test result.bound ≈ -sqrt(29)/4 atol=1e-6
    @test result.certificate_status === :lift_unavailable
    @test isnothing(result.certified_bound) && isnothing(result.certificate) && isnothing(result.pruning)
end
