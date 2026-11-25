# Requires the square seed and symmetry compiler in addition to the verifier.
using Test
using Clarabel
import NCTSSoS as NC

if !isdefined(@__MODULE__, :PauliWorkflowRedesign)
    include("PauliWorkflowRedesign.jl")
end
import .PauliWorkflowRedesign as WF

@testset "Square lift mapping through realification and pruning" begin
    _, xyz = NC.pauli_variables(4)
    bonds = [(1,2),(1,3),(2,4),(3,4)]
    H = 1//8 + sum(a[i]*a[j] for a in xyz for (i,j) in bonds)
    seed = WF.build_square_xxx_seed(H,2;translation_quotient=true)
    before,table,report = WF.workflow_square_symmetry_program(H,seed.rows,2;
        strengthening=seed.strengthening, real_gauge=false,recover_gram_rows=true)
    real = WF._real_workflow_program(before)
    pruned, pruning = WF.workflow_term_sparsity(real; ts_algo=:Maximal,prune_constant_singletons=true)
    form = NC.moment_standard_form(pruned)
    lifts = WF.square_lifts(H,seed,before,real,pruned,pruning,report;form)
    @test !isempty(seed.stationarity_words)
    @test length(lifts) == length(form.cones)
    @test count(l -> l.kind === :rdm,lifts) == length(seed.plaquettes) == 1
    @test all(lifts[i].rows == [WF._wf_word(H,seed.stationarity_words[c.label])]
        for (i,c) in enumerate(form.cones) if c.kind === :zero)
    @test count(c -> c.kind === :zero,form.cones) < length(seed.stationarity_words)
    @test_throws ArgumentError WF.square_lifts(H,seed,before,real,pruned,
        merge(pruning,(;algorithm=:AggregateDualChordal)),report;form)

    model, _ = NC.to_jump(form)
    solver = NC.JuMP.optimizer_with_attributes(Clarabel.Optimizer,"verbose"=>false,
        "chordal_decomposition_enable"=>false)
    NC.JuMP.set_optimizer(model,solver)
    solve_seconds = @elapsed NC.JuMP.optimize!(model)
    dual = WF.pauli_dual_blocks(form,model)
    @test !isnothing(dual)
    finite_only = WF.verify_pauli_dual(H,form,table,dual,lifts;
        symmetry_generators=report.symmetry_generators)
    verification_seconds = @elapsed result = WF.verify_pauli_dual(H,form,table,dual,lifts;
        symmetry_generators=report.symmetry_generators,spin_average=true)
    estimate = NC.JuMP.objective_value(model)
    # The four bonds connect two pairs of spins, so H=4 S_(1,4)⋅S_(2,3)+1/8.
    @test estimate ≈ -63/8 atol=1e-6
    @test -63//8-1//100000 < result.bound <= -63//8
    @test finite_only.bound < result.bound
    @test result.max_lift_error <= 1e-12
    @info "Actual Fourier square SDP" estimate solve_seconds verification_seconds removed=pruning.pruned_singletons
    @info "Square certificate" bound=Float64(result.bound,RoundDown) raw_dual=Float64(result.raw_dual_objective) identity=Float64(result.residual_identity) nonidentity_l1=Float64(result.residual_nonidentity_l1) finite_only=Float64(finite_only.bound,RoundDown)
end

@testset "Removed singleton does not consume a cone" begin
    _, xyz = NC.pauli_variables(9)
    site(x,y) = mod(x,3)+3mod(y,3)+1
    H = sum(a[site(x,y)]*(a[site(x+1,y)]+a[site(x,y+1)]) for a in xyz for x in 0:2 for y in 0:2)
    rows = [one(H); reduce(vcat,xyz)]
    seed = (; stationarity_words=(),plaquettes=(),strengthening=(;equalities=(),inequalities=(),psd_blocks=()))
    before,table,report = WF.workflow_square_symmetry_program(H,rows,3;recover_gram_rows=true)
    real = WF._real_workflow_program(before)
    pruned, pruning = WF.workflow_term_sparsity(real; ts_algo=:Maximal,prune_constant_singletons=true)
    form = NC.moment_standard_form(pruned)
    @test !isempty(pruning.pruned_singletons)
    lifts = WF.square_lifts(H,seed,before,real,pruned,pruning,report;form)
    @test length(lifts) == length(form.cones)
    # Replay all surviving affine rows, even with zero multipliers.
    zero_dual = [zeros(c.side,c.side) for c in form.cones]
    result = WF.verify_pauli_dual(H,form,table,zero_dual,lifts;
        symmetry_generators=report.symmetry_generators)
    @test result.bound == -54
    @test result.max_lift_error <= 1e-12
end
