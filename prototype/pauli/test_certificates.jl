using Test
using Clarabel
using SparseArrays: sparse
import NCTSSoS as NC

if !isdefined(@__MODULE__, :PauliWorkflowRedesign)
    include("PauliWorkflowRedesign.jl")
end
import .PauliWorkflowRedesign as WF

@testset "Exact Pauli SDP dual replay" begin
    _, (_,_,z) = NC.pauli_variables(1)
    H = 1//4 + z[1]
    rows = [one(H), z[1]]
    program,table = NC.plain_program([only(NC.terms(p))[2] for p in rows],NC.terms(H); C=Float64)
    form = NC.moment_standard_form(program)
    lift = [(;kind=:gram, rows)]
    dual = [[1//2 1//2; 1//2 1//2]]
    result = WF.verify_pauli_dual(H,form,table,dual,lift)
    @test result.bound == result.level == -3//4
    @test iszero(result.residual_nonidentity_l1) && iszero(result.residual_identity)

    # Perturb the solver's identity coefficient, not the physical Gram.
    # The exact signed correction must work in both directions.
    for δ in (-1//2^38, 1//2^38)
        b = copy(form.b)
        b[1] += δ
        rounded = NC.MomentStandardForm(form.A,b,form.c,form.offset,form.cones)
        certificate = WF.verify_pauli_dual(H,rounded,table,dual,lift)
        @test certificate.level == -3//4-δ/2
        @test certificate.residual_identity == δ/2
        @test certificate.bound == -3//4
    end

    # An indefinite input with exact balance must still be weakened.
    indefinite = [[-1//8 1//2; 1//2 1//2]]
    repaired = WF.verify_pauli_dual(H,form,table,indefinite,lift)
    @test repaired.raw_dual_objective == -1//8 > -3//4
    @test repaired.bound <= -3//4 && only(repaired.shifts) > 0
    @test iszero(repaired.residual_nonidentity_l1)
    weak = WF.verify_pauli_dual(H,form,table,[[1//4 1//4; 1//4 1//4]],lift)
    @test weak.residual_nonidentity_l1 == 1//2 && weak.bound == -3//4
    @test_throws ArgumentError WF.verify_pauli_dual(H,form,table,dual,[(;kind=:kms_psd)])
    @test_throws ArgumentError WF.verify_pauli_dual(H,form,table,dual,[(;kind=:gram,rows=[one(H),-z[1]])])
end

@testset "Complex RDM compression and shift adjoint" begin
    _, (_,y,_) = NC.pauli_variables(1)
    H = 1//4 + y[1]
    rows = [one(H), y[1]]
    _,table = NC.plain_program([only(NC.terms(p))[2] for p in rows],NC.terms(H); C=Float64)
    # Contracting the local density matrix with (1,i) gives I+Y, not I-Y.
    form = NC.MomentStandardForm(sparse([1],[1],[-1.0]),[1.0],[1.0],0.25,
        [NC.DataCone(:psd,1,1,1,Int32(1))])
    lifts = [(;kind=:rdm,sites=[1],columns=reshape([1,im],2,1))]
    exact = WF.verify_pauli_dual(H,form,table,[fill(1,1,1)],lifts)
    @test exact.bound == -3//4 && isempty(exact.residual)
    repaired = WF.verify_pauli_dual(H,form,table,[fill(-1//2,1,1)],lifts)
    @test only(repaired.shifts) == 1//2
    @test repaired.level == 1//4 && repaired.residual_nonidentity_l1 == 1
    @test repaired.bound == -3//4
    W = eltype(H.data.words)
    K = NC.CliffordMap(NC.AntiunitaryPauli(one(W))) # K flips Y and is not a symmetry of H.
    @test_throws ArgumentError WF.verify_pauli_dual(H,form,table,[fill(1,1,1)],lifts; symmetry_generators=(K,))
    @test_throws ArgumentError WF.verify_pauli_dual(H,form,table,[fill(1,1,1)],lifts; spin_average=true)
end

@testset "Exact collective spin projection is not the Clifford fallback" begin
    _, (x,_,z) = NC.pauli_variables(4)
    H = one(z[1])
    rows = [one(H),prod(z)]
    word(p) = only(NC.terms(p))[2]
    program,table = NC.plain_program(word.(rows),NC.terms(H); C=Float64)
    form = NC.moment_standard_form(program)
    dual = [[1//2 1//2; 1//2 1//2]]
    lifts = [(;kind=:gram,rows)]
    exact = WF.verify_pauli_dual(H,form,table,dual,lifts;spin_average=true)
    finite = WF.verify_pauli_dual(H,form,table,dual,lifts;spin_average=true,spin_max_weight=3)
    @test exact.spin_projection === :su2_exact && finite.spin_projection === :clifford24
    # Sphere fourth moments: E z⁴=1/5 and E x²z²=1/15. For the proper
    # signed-axis orbit these are instead 1/3 and 0, respectively.
    coeffs(r) = Dict(w => c for (c,w) in r.residual)
    mixed = word(x[1]*x[2]*z[3]*z[4])
    @test coeffs(exact)[word(prod(z))] == -1//5
    @test coeffs(exact)[mixed] == -1//15
    @test coeffs(finite)[word(prod(z))] == -1//3
    @test !haskey(coeffs(finite),mixed)
end

@testset "Actual symmetry-reduced SDP with stationarity and RDM" begin
    _, (x,y,z) = NC.pauli_variables(3)
    dot(i,j) = x[i]*x[j]+y[i]*y[j]+z[i]*z[j]
    H = 1//8 + dot(1,2) + dot(2,3) + (1//2)*dot(1,3)
    # In total spin 3/2 the bond eigenvalues are +1. In spin 1/2,
    # sum of the three bonds is -3, so H = 1/8 - 3 - (1/2) σ₁⋅σ₃.
    ground = -27//8
    symmetry = WF.PauliSymmetrySpec(time_reversal_gauge=false,anti_sign_realify=false,
        rotated_lie_detection=false,cartan_u1_weyl=false,su2_mz0_tree=false,spatial_mode=:none)
    sparsity = WF.PauliSparsitySpec(ts_algo=:None,rewrite_complete_local_rdm=false)
    strengthening = WF.PauliStrengtheningSpec(commutator_linear_order=2,rdm_positivity_k=3)
    sites = [1,2,3]
    plan = WF.compile_relaxation(H; order=2,symmetry,sparsity,strengthening)
    @test plan.reports.symmetry.route === :pauli_grading
    @test isempty(plan.symmetry.continuous)
    generators = (plan.symmetry.generators...,plan.symmetry.moment_only...)
    @test !isempty(generators)
    original = only(plan.reports.original_rows)
    W = eltype(original)
    sectors = WF._workflow_grade_rows(original,plan.symmetry.moment_only,plan.symmetry.generators)
    grams = [[WF._wf_polynomial(H.registry,W,row) for row in sector] for sector in sectors]

    # The retained scalar cone label identifies the pre-quotient commutator;
    # moment_standard_form removes symmetry-forced zero rows.
    HQ = WF._wf_exact(H)
    tests = typeof(HQ)[]
    commutators = Set{eltype(tests)}()
    for w in NC.basis_rows(W,NC.Pauli(),sites,2)
        u = WF._wf_word(H,w)
        p = im*(HQ*u-u*HQ)
        (iszero(p) || p in commutators) && continue
        push!(tests,u)
        push!(commutators,p)
    end
    form = plan.standard
    g = 0
    lifts = map(form.cones) do c
        c.kind === :zero && return (;kind=:stationarity,rows=[tests[c.label]])
        if iszero(c.label)
            rows = grams[g+=1]
            return (;kind=:gram,rows,layout=c.side==length(rows) ? :real : :realified)
        end
        return (;kind=:rdm,sites,layout=c.side==2^length(sites) ? :real : :realified)
    end
    @test any(l -> l.kind === :stationarity,lifts) && any(l -> l.kind === :rdm,lifts)
    @test g == length(grams) >= 2
    @test any(l -> get(l,:layout,:none) === :realified,lifts)

    solver = NC.JuMP.optimizer_with_attributes(Clarabel.Optimizer,"verbose"=>false,
        "chordal_decomposition_enable"=>false)
    for mode in (:moment,:sos)
        model, variables = NC.to_jump(mode === :moment ? form : NC.sos_standard_form(form))
        NC.JuMP.set_optimizer(model,solver)
        elapsed = @elapsed NC.JuMP.optimize!(model)
        status = NC.JuMP.termination_status(model)
        @test status in (NC.JuMP.MOI.OPTIMAL,NC.JuMP.MOI.ALMOST_OPTIMAL)
        dual = WF.pauli_dual_blocks(form,model; sos_variables=mode === :sos ? variables : nothing)
        @test !isnothing(dual)
        certificate = WF.verify_pauli_dual(H,form,plan.table,dual,lifts; symmetry_generators=generators)
        estimate = NC.JuMP.objective_value(model)
        @test estimate ≈ Float64(ground) atol=1e-6
        @test certificate.status === :verified_rational
        @test ground-1//100000 < certificate.bound <= ground
        @test abs(estimate-Float64(certificate.raw_dual_objective)) <= 1e-6
        @test certificate.max_lift_error <= 1e-12
        @test certificate.raw_balance_inf <= 1e-6
        raw_replay = WF.verify_pauli_dual(H,form,plan.table,dual,lifts)
        @test raw_replay.bound < certificate.bound
        stationarity = findfirst(l -> l.kind === :stationarity,lifts)
        perturbed = deepcopy(dual)
        perturbed[stationarity][1] += 1//8
        replay = WF.verify_pauli_dual(H,form,plan.table,perturbed,lifts)
        u = only(lifts[stationarity].rows)
        p = im*(HQ*u-u*HQ)
        poly(r) = WF._wf_polynomial(H.registry,W,r.residual)
        δ = WF._PauliQ(perturbed[stationarity][1]) - WF._PauliQ(dual[stationarity][1])
        @test poly(raw_replay)-poly(replay) == δ*p
        @test replay.bound <= ground
        diagnostic = (; mode, status, elapsed, estimate, ground,
            raw_dual=Float64(certificate.raw_dual_objective), repaired_level=Float64(certificate.level),
            certified_bound=Float64(certificate.bound,RoundDown), raw_balance=Float64(certificate.raw_balance_inf),
            repaired_balance=Float64(certificate.repaired_balance_inf), max_shift=Float64(maximum(certificate.shifts)),
            identity=Float64(certificate.residual_identity), nonidentity_l1=Float64(certificate.residual_nonidentity_l1),
            unaveraged_bound=Float64(raw_replay.bound,RoundDown),
            stationarity_dual_max=maximum(maximum(abs,d) for (d,l) in zip(dual,lifts) if l.kind === :stationarity),
            rdm_dual_max=maximum(maximum(abs,d) for (d,l) in zip(dual,lifts) if l.kind === :rdm))
        @info "Actual Pauli dual replay (Pauli units)" diagnostic
        @info "Exact certified bound" mode bound=certificate.bound
    end
end
