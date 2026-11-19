import NCTSSoS
realpath(pkgdir(NCTSSoS)) == realpath(joinpath(@__DIR__, "../../..")) ||
    error("the Julia environment does not load NCTSSoS from the reviewed checkout")
import PauliWorkflowRedesign
realpath(pathof(PauliWorkflowRedesign)) == realpath(joinpath(@__DIR__, "..", "PauliWorkflowRedesign.jl")) ||
    error("the Julia environment does not load the reviewed Pauli prototype")
using .PauliWorkflowRedesign: compile_relaxation, solve_ground_state,
    PauliSymmetrySpec, PauliSparsitySpec, PauliStrengtheningSpec

function prototype_inputs(lattice, L, extra, strengthened, mode)
    n = lattice == "chain" ? L : L^2
    _, (x, y, z) = NCTSSoS.pauli_variables(n)
    H = (1 / (4n)) * sum(x[i]*x[j] + y[i]*y[j] + z[i]*z[j] for (i, j) in bonds(lattice, L))
    order = 2
    radius = lattice == "chain" ? extra+1 : 4
    sparsity = PauliSparsitySpec(support_mode=mode, power_graph_radius=mode === :graph_distance ? radius : 2,
        ts_algo=:None,rewrite_complete_local_rdm=false)
    strengthening = strengthened ? PauliStrengtheningSpec(commutator_linear_order=lattice == "chain" ? 7 : 5,
        kms_psd_order=lattice == "chain" ? 2 : 0,rdm_positivity_k=8) : PauliStrengtheningSpec()
    return H,order,sparsity,strengthening
end

# Compare actual words before timing, expanding Wang's one triplet
# representative into its three axis-equivalent row families.
function prototype_audit(lattice,L,extra,strengthened,mode)
    NC = NCTSSoS
    wf = PauliWorkflowRedesign
    H,order,sparsity,_ = prototype_inputs(lattice,L,extra,strengthened,mode)
    W = eltype(H.data.words)
    n = H.registry.nsites
    family = wf.workflow_support(H,order;support_mode=mode,sparsity.power_graph_radius,
        sparsity.local_product_depth,sparsity.cs_algo).family
    rows = [W[last(only(NC.terms(p))) for p in block] for block in family]
    actual = Set(only(rows))
    function stock_word(w,shift)
        result = one(W)
        for a in w
            site,axis = divrem(Int(a)-1,3)
            if lattice == "chain"
                site += 1
            else
                x,y = QMBCertify.location(site+1)
                site = x+(y-1)*L
            end
            result = last(NC.mul(result,NC.pauli_letter(W,site,mod(axis+shift,3))))
        end
        return result
    end
    expected = Set([one(W)])
    union!(expected,(stock_word(w,0) for w in QMBCertify.get_basis(L,0,order;extra,lattice)))
    union!(expected,(stock_word(w,s) for w in QMBCertify.get_basis(L,1,order;extra,lattice) for s in 0:2))
    sym = wf.workflow_symmetry(H;family=rows,spin=false,uniform_parity_blocks=true)
    sym,layout = wf._wang_spatial_symmetry(H,nothing,rows,sym)
    letters = [[NC.pauli_letter(W,i,a) for a in 0:2] for i in 1:n]
    found = isempty(sym.generators) ? nothing : NC.discover_lattice(sym.generators,letters;required=false)
    orders = isnothing(found) ? Int[] : found.translations.orders
    point = isnothing(found) ? 0 : length(NC.OrbitGroup(found.point_generators))
    expected_orders,expected_point = lattice == "chain" ? ([L],48) : ([L,L],192)
    actions = (sym.generators...,sym.moment_only...)
    global_actions = [wf._workflow_axes(W,n);wf._wang_uniform_parities(H,[NC._sector_terms(H,W).terms]);
        NC.CliffordMap(NC.AntiunitaryPauli(one(W)))]
    retained = count(in(actions),global_actions)
    rows_match = actual == expected
    symmetry_match = !isnothing(layout) && sort(orders) == expected_orders && point == expected_point &&
        retained == length(global_actions)
    return Dict("comparison_lane" => rows_match && symmetry_match && mode != :local_products ? "matched_support" : "exploratory",
        "row_span_match"=>rows_match,"expected_rows"=>length(expected),"selected_rows"=>length(actual),
        "wang_only_rows"=>length(setdiff(expected,actual)),"candidate_only_rows"=>length(setdiff(actual,expected)),
        "retained_translation_orders"=>orders,"expected_translation_orders"=>expected_orders,
        "retained_point_order"=>point,"expected_point_order"=>expected_point,
        "retained_global_actions"=>retained,"expected_global_actions"=>length(global_actions),
        "continuous_spin_enabled"=>false,"verified_symmetry_match"=>symmetry_match,
        "selected_row_triangle_forms"=>length(actual)*(length(actual)+1)÷2)
end

function prototype_run(result, path, lattice, L, extra, strengthened, mode)
    started = time_ns()
    H,order,sparsity,strengthening = prototype_inputs(lattice,L,extra,strengthened,mode)
    plan = compile_relaxation(H; order, sparsity, strengthening, wang_matched=true, wang_lso_span=strengthened ? 4 : 0)
    result["build_seconds"] = (time_ns() - started) / 1e9
    result["prototype_settings"] = "wang_matched"
    result["graph_radius"] = sparsity.power_graph_radius
    result["local_product_depth"] = sparsity.local_product_depth
    result["prototype_moi_form"] = string(plan.moi_form)
    result["prototype_order"] = order
    result["prototype_support_mode"] = string(plan.reports.support.mode)
    result["prototype_symmetry_route"] = string(plan.reports.symmetry.route)
    reduction = get(plan.reports.symmetry,:reduction,(;))
    result["translation_orders"] = get(reduction,:translation_orders,Int[])
    result["point_order"] = get(reduction,:point_order,0)
    result["complex_cones_realified_by_doubling"] = plan.reports.symmetry.doubled
    result["prototype_row_counts"] = plan.reports.support.row_counts
    result["prototype_generation"] = repr(plan.reports.support.generation)
    result["prototype_active_relaxations"] = string.(collect(plan.provenance.active_relaxations))
    result["prototype_active_strengthenings"] = string.(collect(plan.provenance.active_strengthenings))
    result["moment_variable_count"] = length(plan.standard.c)
    result["constrained_moment_variable_count"] = length(plan.constraint_ids)-1
    result["affine_term_count"] = count(!iszero,plan.standard.A.nzval)
    result["moment_conic_coordinate_count"] = length(plan.standard.b)
    result["moment_matrix_sides"] = [cone.side for cone in plan.standard.cones if cone.side > 0]
    result["prototype_rdm_region_count"] = plan.reports.analysis.rdm.regions
    result["rdm_psd_blocks"] = plan.reports.analysis.rdm.count
    result["active_equality_rows"] = plan.reports.analysis.active_equality_rows
    wang = plan.reports.analysis.wang
    result["rdm_sites"] = wang.rdm.region
    result["rdm_psd_sides"] = wang.rdm.sides
    result["pso_input_rows"] = wang.pso.rows
    result["pso_psd_sides"] = wang.pso.sides
    if strengthened
        result["lso_pool"] = Dict(string(k)=>v for (k,v) in pairs(wang.lso))
    end
    save_result(path, result)
    if strengthened
        result["rdm_psd_blocks"] > 0 && (strengthening.kms_psd_order == 0 || !isempty(wang.pso.sides)) ||
            error("prototype PSD strengthening missing")
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
