# Run with Clarabel available in the Julia environment, e.g. the package test environment.
using Test
using Clarabel
import NCTSSoS as NC

if !isdefined(@__MODULE__, :PauliWorkflowRedesign)
    include("PauliWorkflowRedesign.jl")
end
using .PauliWorkflowRedesign

@testset "Pauli workflow integration" begin
    # Clarabel 0.11.1 can fail in its own PSD completion on this dependency set.
    solver = NC.JuMP.optimizer_with_attributes(Clarabel.Optimizer,"verbose"=>false,
        "chordal_decomposition_enable"=>false)
    _, (x,y,z) = NC.pauli_variables(1)
    no_symmetry = PauliSymmetrySpec(time_reversal_gauge=false,anti_sign_realify=false,
        klein4_sign=false,xyz_permutation=false,rotated_lie_detection=false,
        cartan_u1_weyl=false,su2_mz0_tree=false,spatial_mode=:none)
    ambient = PauliSparsitySpec(ts_algo=:None,rewrite_complete_local_rdm=false)

    @testset "Sweep guard and coefficient updates" begin
        H = 0.25 - x[1]
        for form in (:dualized,:primalized)
            plan = pauli_nctssos(H; order=1,preset=:ExactNPA,moi_form=form)
            @test plan.provenance.bound_status === :unavailable
            @test plan.provenance.exact_wrt_dense_npa
            first = solve_ground_state(plan,H; optimizer=solver)
            @test first.bound ≈ -0.75 atol=1e-6
            @test !first.solver_model_reused
            second = solve_ground_state(plan,0.5 - 2x[1]; optimizer=solver)
            @test second.bound ≈ -1.5 atol=1e-6
            @test second.plan === plan && second.reuse === :reused && second.solver_model_reused
            @test second.timings.rebuild_seconds == 0
            @test 0 <= second.timings.update_seconds <= second.timings.total_seconds
            @test 0 <= second.timings.optimizer_seconds <= second.timings.total_seconds
            @test second.provenance.bound_status === :numerical_estimate
            @test !second.provenance.rigorous_lower_bound
            maximum = solve_ground_state(plan,H; sense=:max,optimizer=solver)
            @test maximum.bound ≈ 1.25 atol=1e-6
            @test solve_ground_state(plan,H;optimizer=solver).bound ≈ -0.75 atol=1e-6
        end

        minus_x = compile_relaxation(-x[1];order=1,preset=:ExactNPA)
        changed = solve_ground_state(minus_x,-x[1]-z[1];optimizer=solver)
        @test changed.reuse === :support_changed && changed.plan !== minus_x
        @test !changed.solver_model_reused
        @test changed.bound ≈ -sqrt(2) atol=1e-6

        _, (x2,y2,z2) = NC.pauli_variables(2)
        H2 = x2[1]*x2[2]+y2[1]*y2[2]+z2[1]*z2[2]
        isotropic = compile_relaxation(H2;order=1,preset=:ExactNPA)
        broken = solve_ground_state(isotropic,H2+z2[1]*z2[2];optimizer=solver)
        @test broken.reuse === :symmetry_changed
        @test broken.bound ≈ -4 atol=1e-5
    end

    @testset "Hamiltonian-dependent cone refresh and sectors" begin
        plan = compile_relaxation(-x[1];order=1,symmetry=no_symmetry,sparsity=ambient,
            strengthening=PauliStrengtheningSpec(kms_psd_order=1))
        @test plan.h_dependent && :kms_psd in plan.provenance.active_strengthenings
        # A stale KMS cone for -X imposes ⟨X⟩≥0 and cannot attain the new minimum.
        changed = solve_ground_state(plan,x[1];optimizer=solver)
        @test changed.reuse === :hamiltonian_constraints_changed
        @test changed.bound ≈ -1 atol=1e-6
        @test_throws ArgumentError solve_ground_state(plan,-x[1];sense=:max,optimizer=solver)

        sector = solve_pauli_ground_state(-z[1];order=1,optimizer=solver,
            strengthening=PauliStrengtheningSpec(target_mz_sector=-1))
        @test sector.provenance.sector_targets == (:mz,)
        @test !sector.provenance.exact_wrt_dense_npa
        @test sector.bound ≈ 1 atol=1e-6

        # Ferromagnetic three-site triangle: H=9/2-2S², so the J=1/2
        # target has E=+3 although the unrestricted minimum is -3.
        _, (sx,sy,sz) = NC.pauli_variables(3)
        triangle = -sum(a[i]*a[j] for a in (sx,sy,sz) for i in 1:3 for j in i+1:3)
        spin = solve_pauli_ground_state(triangle;order=1,optimizer=solver,symmetry=no_symmetry,sparsity=ambient,
            strengthening=PauliStrengtheningSpec(target_spin_j=1//2))
        @test spin.provenance.sector_targets == (:spin_j,)
        @test spin.plan.reports.analysis.target.spin_j == 1//2
        @test spin.bound ≈ 3 atol=1e-6
    end

    @testset "Verified dual witness and numerical interval" begin
        H = -x[1] + 0.5y[1] + 0.75z[1]
        exact_energy = -sqrt(1+0.5^2+0.75^2)
        for form in (:dualized,:primalized)
            plan = compile_relaxation(H; order=1,symmetry=no_symmetry,sparsity=ambient,moi_form=form)
            result = solve_ground_state(plan,H;optimizer=solver,certify=true)
            @test result.certificate_status === :verified_rational
            @test result.certificate.status === :verified_rational
            @test Float64(result.certificate.bound) <= exact_energy + 1e-12
            @test Float64(result.certificate.bound) ≈ exact_energy atol=5e-6
            @test result.provenance.bound_status === :numerical_estimate
            @test result.certificate_provenance.rigorous_lower_bound
            x_id = only(NC.moment_terms(plan.table,only(NC.terms(x[1]))[2]))[1]
            @test result.moments[x_id-1] ≈ -1/exact_energy atol=1e-5
        end

        H = 0.25 - x[1]
        plan = compile_relaxation(H; order=1,preset=:ExactHierarchy)
        numerical = bound_observable_interval(plan,H,z[1],-0.75,-0.55;optimizer=solver)
        @test numerical.status === :numerical_estimate && numerical.window_status === :assumed_window
        @test numerical.lower ≈ -0.6 atol=1e-5
        @test numerical.upper ≈ 0.6 atol=1e-5
        # Two copies of the projector square sum to I-X exactly. The trial
        # state (1,2) gives ⟨X⟩=4/5, a nonzero-width verified window.
        projector = (one(x[1])-x[1])/2
        lower = certify_rational_sos(H,-3//4;squares=(projector,projector))
        window = certify_energy_window(H,lower,[1,2])
        @test window.lower == -3//4 && window.upper == -11//20
        bounded = bound_observable_interval(plan,H,z[1],window;optimizer=solver)
        @test bounded.status === :numerical_estimate && bounded.window_status === :verified_window
        @test bounded.lower ≈ -0.6 atol=1e-5
        @test bounded.upper ≈ 0.6 atol=1e-5
        window_plan = compile_relaxation(H;order=1,energy_window=window)
        @test :energy_window in window_plan.provenance.active_strengthenings
        @test solve_ground_state(window_plan,H;optimizer=solver).bound ≈ -0.75 atol=1e-6
        @test_throws ArgumentError bound_observable_interval(window_plan,H,z[1],-1.0,1.0;optimizer=solver)
        stationary = compile_relaxation(z[1];order=1,
            strengthening=PauliStrengtheningSpec(commutator_linear_order=1))
        @test_throws ArgumentError bound_observable_interval(stationary,z[1],x[1],-1.0,1.0;optimizer=solver)
        @test_throws ArgumentError solve_ground_state(stationary,z[1];sense=:max,optimizer=solver)
        free = bound_observable_interval(compile_relaxation(z[1];order=1),z[1],x[1],-1.0,1.0;optimizer=solver)
        @test free.lower ≈ -1 atol=1e-6
        @test free.upper ≈ 1 atol=1e-6
        no_lift = solve_ground_state(plan,H;optimizer=solver,certify=true)
        @test no_lift.certificate_status === :lift_unavailable && isnothing(no_lift.certificate)
        @test_throws ArgumentError bound_observable_interval(plan,-x[1],z[1],window;optimizer=solver)
    end

    @testset "Preset dispatch, bounded RDM selection, density rewrite" begin
        @test resolve_workflow_preset(:FastScalable) == resolve_workflow_preset(:GeneralGraphFast)
        @test_throws ArgumentError resolve_workflow_preset(:wrong)
        @test_throws ArgumentError compile_relaxation(-x[1];preset=:ExactNPA,sparsity=ambient)
        @test_throws ArgumentError compile_relaxation(-x[1];symmetry=PauliSymmetrySpec(drop_irreps=[0]))
        @test_throws ArgumentError compile_relaxation(-x[1];order=1,preset=:ScalableLattice)

        _, (rx,ry,rz) = NC.pauli_variables(4)
        ring = sum(rz[i]*rz[mod1(i+1,4)] for i in 1:4) + 0.5sum(rx)
        lattice = compile_relaxation(ring;order=1,preset=:ScalableLattice)
        @test lattice.reports.support.mode === :contiguous_lattice
        @test lattice.reports.symmetry.route === :lattice
        @test lattice.reports.analysis.rdm.selection === :local_regions
        @test lattice.reports.analysis.rdm.count == 1
        torus = rz[1]*rz[2]+rz[3]*rz[4]+2*(rz[1]*rz[3]+rz[2]*rz[4])+0.5sum(rx)
        named_torus = compile_relaxation(torus;order=1,preset=:ScalableLattice,lattice_shape=(2,2))
        @test named_torus.reports.symmetry.lattice_shape == (2,2)
        @test named_torus.reports.symmetry.reduction.translation_order == 4
        @test named_torus.reports.analysis.rdm.regions == 1
        generic = compile_relaxation(ring;order=1,preset=:GeneralGraphFast)
        @test generic.reports.support.mode === :correlative_cliques
        @test generic.reports.term_sparsity.algorithm === :MMD
        tight = compile_relaxation(-x[1];order=1,preset=:TightGroundState)
        @test :lightcone in tight.provenance.active_strengthenings
        @test tight.reports.analysis.rdm.count == 1

        rewritten = compile_relaxation(-x[1];order=1,symmetry=no_symmetry)
        audit = rewritten.reports.support.density_rewrite
        @test only(audit.applied).before == 4 && only(audit.applied).after == 2
        @test rewritten.provenance.exact_wrt_dense_npa
        @test solve_ground_state(rewritten,-x[1];optimizer=solver).bound ≈ -1 atol=1e-6
    end
end
