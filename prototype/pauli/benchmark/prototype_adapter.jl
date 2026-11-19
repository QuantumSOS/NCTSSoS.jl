import NCTSSoS
realpath(pkgdir(NCTSSoS)) == realpath(joinpath(@__DIR__, "../../..")) ||
    error("the Julia environment does not load NCTSSoS from the reviewed checkout")
import PauliWorkflowRedesign
realpath(pathof(PauliWorkflowRedesign)) == realpath(joinpath(@__DIR__, "..", "PauliWorkflowRedesign.jl")) ||
    error("the Julia environment does not load the reviewed Pauli prototype")
using .PauliWorkflowRedesign: compile_relaxation, solve_ground_state,
    PauliSymmetrySpec, PauliSparsitySpec, PauliStrengtheningSpec

function prototype_run(result, path, lattice, L, strengthened)
    n = lattice == "chain" ? L : L^2
    started = time_ns()
    _, (x, y, z) = NCTSSoS.pauli_variables(n)
    H = (1 / (4n)) * sum(x[i]*x[j] + y[i]*y[j] + z[i]*z[j] for (i, j) in bonds(lattice, L))
    order = strengthened ? 3 : 2
    preset = strengthened ? :TightGroundState : :GeneralGraphFast
    plan = compile_relaxation(H; order, preset)
    result["build_seconds"] = (time_ns() - started) / 1e9
    result["prototype_settings"] = string(preset)
    result["prototype_moi_form"] = string(plan.moi_form)
    result["prototype_order"] = order
    result["prototype_support_mode"] = string(plan.reports.support.mode)
    result["prototype_symmetry_route"] = string(plan.reports.symmetry.route)
    result["prototype_active_relaxations"] = string.(collect(plan.provenance.active_relaxations))
    result["prototype_active_strengthenings"] = string.(collect(plan.provenance.active_strengthenings))
    result["moment_variable_count"] = length(plan.standard.c)
    result["moment_conic_coordinate_count"] = length(plan.standard.b)
    result["moment_matrix_sides"] = [cone.side for cone in plan.standard.cones if cone.side > 0]
    result["prototype_rdm_region_count"] = plan.reports.analysis.rdm.regions
    result["rdm_psd_blocks"] = plan.reports.analysis.rdm.count
    result["active_equality_rows"] = plan.reports.analysis.active_equality_rows
    save_result(path, result)
    if strengthened
        :stationarity in plan.provenance.active_strengthenings && result["active_equality_rows"] > 0 &&
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
