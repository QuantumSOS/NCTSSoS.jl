import NCTSSoS
realpath(pkgdir(NCTSSoS)) == realpath(joinpath(@__DIR__, "../../..")) ||
    error("the Julia environment does not load NCTSSoS from the reviewed checkout")
import PauliWorkflowRedesign
realpath(pathof(PauliWorkflowRedesign)) == realpath(joinpath(@__DIR__, "..", "PauliWorkflowRedesign.jl")) ||
    error("the Julia environment does not load the reviewed Pauli prototype")
using .PauliWorkflowRedesign: compile_relaxation, solve_ground_state,
    PauliSymmetrySpec, PauliSparsitySpec, PauliStrengtheningSpec

function prototype_run(result, path, lattice, L, pair_range, control, strengthened)
    n = lattice == "chain" ? L : L^2
    started = time_ns()
    _, (x, y, z) = NCTSSoS.pauli_variables(n)
    H = (1 / (4n)) * sum(x[i]*x[j] + y[i]*y[j] + z[i]*z[j] for (i, j) in bonds(lattice, L))
    full = control && !strengthened
    order = strengthened ? 3 : 2
    # Pin the equivalent moment form for the unstrengthened sparse campaign.
    # At chain-40, the auto-selected SOS form stalled at the same tolerances.
    moi_form = full || strengthened ? :auto : :dualized
    if full
        plan = compile_relaxation(H; order, preset=:ExactNPA, moi_form)
    else
        rdm_regions = strengthened ? unique([sort([mod1(i+j,n) for j in 0:3]) for i in 1:n]) : Vector{Int}[]
        strengthening_options = strengthened ? (;
            strengthening=PauliStrengtheningSpec(commutator_linear_order=3, rdm_positivity_k=4), rdm_regions) : (;)
        plan = compile_relaxation(H; order, moi_form,
            symmetry=PauliSymmetrySpec(spatial_mode=:circulant_dft),
            sparsity=PauliSparsitySpec(support_mode=:contiguous_lattice,
                two_body_range_r=pair_range, cs_algo=:NoElimination, ts_algo=:Maximal, ts_step=0),
            strengthening_options...)
    end
    result["build_seconds"] = (time_ns() - started) / 1e9
    if full
        plan.reports.support.mode === :full || error("ED control did not use full support")
    else
        plan.reports.support.mode === :contiguous_lattice && plan.reports.symmetry.route === :lattice &&
            :contiguous_lattice in plan.provenance.active_relaxations || error("lattice route not active")
    end
    result["prototype_settings"] = full ? "ExactNPA; moi_form=auto" :
        "moi_form=$moi_form; contiguous_lattice; range=$pair_range; circulant_dft; NoElimination; Maximal; ts_step=0"
    result["prototype_moi_form"] = string(plan.moi_form)
    result["prototype_order"] = order
    result["prototype_support_mode"] = string(plan.reports.support.mode)
    result["prototype_symmetry_route"] = string(plan.reports.symmetry.route)
    result["prototype_active_relaxations"] = string.(collect(plan.provenance.active_relaxations))
    result["prototype_active_strengthenings"] = string.(collect(plan.provenance.active_strengthenings))
    result["moment_variable_count"] = length(plan.standard.c)
    result["moment_conic_coordinate_count"] = length(plan.standard.b)
    result["moment_matrix_sides"] = [cone.side for cone in plan.standard.cones if cone.side > 0]
    if strengthened
        result["prototype_commutator_linear_order"] = 3
        result["prototype_rdm_positivity_k"] = 4
        result["prototype_requested_rdm_regions"] = rdm_regions
        result["prototype_rdm_region_count"] = plan.reports.analysis.rdm.regions
        result["rdm_psd_blocks"] = plan.reports.analysis.rdm.count
        result["stationarity_equalities"] = sum(c.length for c in plan.standard.cones
            if c.kind === :zero && c.label > 0; init=0)
    end
    save_result(path, result)
    if strengthened
        :stationarity in plan.provenance.active_strengthenings && result["stationarity_equalities"] > 0 &&
            result["rdm_psd_blocks"] > 0 || error("prototype strengthening missing")
    end

    # Capture the actual optimizer so status and gap are observed, not inferred
    # from the prototype returning a finite objective. No second solve is done.
    optimizer = Ref{Any}(nothing)
    factory = () -> (optimizer[] = MOI.instantiate(OPTIMIZER))
    result["solve_seconds"] = @elapsed r = solve_ground_state(plan, H; optimizer=factory, certify=false)
    model = optimizer[]
    model === nothing && error("prototype did not call the configured Mosek factory")
    result["prototype_timings_seconds"] = Dict(string(k) => v for (k,v) in pairs(r.timings))
    result["optimizer_reported_seconds"] = MOI.get(model, MOI.SolveTimeSec())
    result["prototype_reuse"] = string(r.reuse)
    result["prototype_provenance"] = repr(r.provenance)
    r.reuse === :reused && !r.solver_model_reused || error("expected original plan and a fresh solver model")
    accepted = record_status!(result, model)
    result["termination_status"] == string(r.status) || error("prototype/optimizer termination mismatch")
    return accepted ? r.bound : nothing
end
