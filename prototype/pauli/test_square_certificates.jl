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

@testset "Folded replay equals the ambient 4×4 certificate" begin
    L = 4
    _, xyz = NC.pauli_variables(L^2)
    site(x,y) = mod(x,L) + L*mod(y,L) + 1
    H = (1//64)*sum(a[site(x,y)]*(a[site(x+1,y)]+a[site(x,y+1)])
        for a in xyz for x in 0:L-1 for y in 0:L-1)
    plan = WF.compile_minimal_square(H,L;recover_gram_rows=true,recover_folded_gram=true)
    form = plan.standard
    model,_ = NC.to_jump(form)
    optimizer = NC.JuMP.optimizer_with_attributes(Clarabel.Optimizer,"verbose"=>false,
        "chordal_decomposition_enable"=>false)
    NC.JuMP.set_optimizer(model,optimizer)
    NC.JuMP.optimize!(model)
    dual = WF.pauli_dual_blocks(form,model)
    folded = plan.lifts()
    grams = plan.symmetry.gram_rows()
    expanded = [l.kind === :fourier_gram ? (;kind=:gram,rows=grams[l.source],l.layout,l.indices) : l for l in folded]
    generators = plan.symmetry.symmetry_generators
    a = WF.verify_pauli_dual(H,form,plan.table,dual,expanded; symmetry_generators=generators,spin_average=true)
    b = WF.verify_pauli_dual(H,form,plan.table,dual,folded; symmetry_generators=generators,spin_average=true)
    @test a.residual_representation === :ambient_pauli && b.residual_representation === :square_orbit_mass
    @test a.bound == b.bound && a.residual_identity == b.residual_identity
    @test a.residual_nonidentity_l1 == b.residual_nonidentity_l1
    W = eltype(H.data.words)
    square = WF.SquareOrbitMoments(W,L,true)
    masses = Dict{W,WF._PauliCQ}()
    for (c,w) in a.residual
        NC._sector_add!(masses,WF._square_canonical(square,w),c)
    end
    @test masses == Dict(w=>c for (c,w) in b.residual)
end

@testset "Non-character Fourier samples and a nontrivial stabilizer" begin
    L = 4
    _, (x,_,z) = NC.pauli_variables(L^2)
    H = WF._wf_exact(one(x[1]))
    W = eltype(H.data.words)
    word(p) = only(NC.terms(p))[2]
    shifts = WF._workflow_translations(W,L^2,(L,L),:auto)
    tr = NC.TranslationGroup(shifts,[L,L])
    roots = word.([x[1],z[1]*z[3]])
    words = unique([last(NC.act(t,w)) for w in roots for t in tr.elements])
    orbits = NC.AbelianRows(words,tr)
    a = only(findall(ws -> length(ws)==L^2,orbits.words))
    b = only(findall(ws -> length(ws)==L^2÷2,orbits.words))
    # Deliberately not roots of unity. The samples are exactly constant on
    # integer phase classes, including the antipodal pair's stabilizer.
    chars = ComplexF64[iseven(t[1]) ? 1.0 : -1.0 + 0.017im for t in tr.coordinates]
    cols = [[(0.62+0.3im,a),(0.27-0.5im,b)],[(0.4-0.29im,a),(0.73+0.15im,b)]]
    fourier = (;orbits,characters=chars,momentum=(2,0),cols)
    lift = (;kind=:fourier_gram,fourier,layout=:realified)
    square = WF.SquareOrbitMoments(W,L,false)
    fold = WF._wf_square_fold_context(H,square,[lift],())
    rows = WF._square_fourier_lift(H.registry,orbits,chars,cols)
    full = WF._wf_dual_entries(H,(;kind=:gram,rows,layout=:realified),4)
    compact = WF._wf_dual_entries(H,lift,4,fold)
    for j in 1:4, i in j:4
        reference = Dict{W,WF._PauliCQ}()
        for (c,w) in NC.terms(full(i,j))
            NC._sector_add!(reference,WF._square_canonical(square,w),c)
        end
        @test reference == Dict(w=>c for (c,w) in NC.terms(compact(i,j)))
    end
    inconsistent = copy(chars)
    inconsistent[3] += 0.002im # Same integer class as the first sample.
    @test_throws ArgumentError WF._wf_folded_gram(H,merge(fourier,(;characters=inconsistent)),fold)
    # Samples pass the class test but this momentum is not constant under
    # translation by two, which fixes the antipodal pair.
    unsupported = ComplexF64[im^t[1] for t in tr.coordinates]
    @test_throws ArgumentError WF._wf_folded_gram(H,merge(fourier,(;momentum=(1,0),characters=unsupported)),fold)
end
