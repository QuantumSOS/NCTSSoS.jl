# The prepared model is private to a plan and optimizer factory. Neither the
# cone restriction nor the symbolic coordinate map is rebuilt for an update.
function _prepare_numeric(plan, optimizer)
    moment = plan.standard
    if plan.moi_form === :dualized
        model, variables = NC.to_jump(moment)
        balance = nothing
    else
        sos = plan.sos
        model = JuMP.Model()
        JuMP.@variable(model, z[1:length(sos.c)])
        balance = JuMP.@constraint(model, sos.A*z .== sos.b)
        NC._jump_cones!(model, z, sos.cones)
        JuMP.@objective(model, Max, sum(sos.c[i]*z[i] for i in eachindex(z)) + sos.offset)
        variables = z
    end
    JuMP.set_optimizer(model, optimizer)
    JuMP.MOI.Utilities.attach_optimizer(model)
    return (; model, variables, balance, coefficients=copy(plan.objective_data))
end

function _update_numeric!(prepared, objective, sos)
    old = prepared.coefficients
    for id in union(keys(old), keys(objective))
        coefficient = get(objective, id, 0.0)
        get(old, id, 0.0) == coefficient && continue
        if id == 1
            # Set the constant on a copy; the solver sees the changed objective
            # even on backends that retain a copy of the previous expression.
            f = JuMP.objective_function(prepared.model, JuMP.AffExpr)
            f.constant = coefficient
            JuMP.set_objective_function(prepared.model, f)
        elseif sos
            JuMP.set_normalized_rhs(prepared.balance[id-1], -coefficient)
        else
            JuMP.set_objective_coefficient(prepared.model, prepared.variables[id-1], coefficient)
        end
    end
    empty!(old)
    merge!(old, objective)
    return nothing
end

function _solve_numeric(plan, data, sense, optimizer)
    sign = sense === :max ? -1.0 : 1.0
    sos = plan.moi_form === :primalized
    warm = haskey(plan.cache, optimizer)
    setup_seconds = @elapsed prepared = get!(plan.cache, optimizer) do
        _prepare_numeric(plan, optimizer)
    end
    update_seconds = @elapsed begin
        objective = sign == 1 ? data : Dict(id => -c for (id,c) in data)
        try
            _update_numeric!(prepared, objective, sos)
        catch
            # A failed incremental update could have modified some coefficients.
            delete!(plan.cache, optimizer)
            rethrow()
        end
    end
    optimizer_seconds = @elapsed JuMP.optimize!(prepared.model)
    model = prepared.model
    status = JuMP.termination_status(model)
    accepted = status == JuMP.MOI.OPTIMAL && JuMP.has_values(model)
    bound = accepted ? sign * Float64(JuMP.objective_value(model)) : nothing
    if !isnothing(bound) && !isfinite(bound)
        bound = nothing
    end
    moments = if isnothing(bound)
        Float64[]
    elseif sos
        # The SOS maximization is written A'z == -c; JuMP's equality duals
        # have the moment sign here (unlike raw native-solver shadow prices).
        JuMP.has_duals(model) ? JuMP.dual.(prepared.balance) : Float64[]
    else
        JuMP.value.(prepared.variables)
    end
    return (; bound, status, moments, raw_status=JuMP.raw_status(model),
        solver_model_reused=warm, setup_seconds, update_seconds, optimizer_seconds)
end

"""
    solve_ground_state(plan, H; sense=:min, optimizer, certify=..., prune_budget=nothing)

Verify a new Hamiltonian before reusing a compiled model. On a guard failure,
recompile with the same options. Timings separate guard, rebuild, model setup,
coefficient update and the optimizer call (which may include vendor presolve
or reload). Serialize calls on the same plan; its optimizer cache is mutable.

`bound` is the numerical solver objective. `certified_bound` is the exact
lower bound from `certificate`, or `nothing` if unavailable. A nonnegative
`prune_budget` requests certification and, for a supported lift, returns the
ambient factor pruning report in `pruning`. See `prune_sos_factors` for the
budget semantics. No certificate or pruning report is produced for a cone
whose adjoint to the original Pauli algebra is unavailable.
"""
function solve_ground_state(plan::CompiledSDPPlan, H;
    sense::Symbol=:min, optimizer=nothing, certify::Bool=plan.certify_default,
    prune_budget::Union{Nothing,Real}=nothing)
    optimizer = _optimizer_factory(optimizer)
    budget = isnothing(prune_budget) ? nothing : _wf_loss_budget(prune_budget)
    certify = certify || !isnothing(budget)
    sense in (:min, :max) || throw(ArgumentError("sense must be :min or :max."))
    sense === :max && plan.h_dependent && throw(ArgumentError("Maximization is not supported with Hamiltonian-dependent conditions."))
    sense === :max && certify && throw(ArgumentError("The lower-bound certificate hook requires sense=:min."))

    total_start = time_ns()
    guard_seconds = @elapsed begin
        _check_hamiltonian(H)
        H.registry === plan.hamiltonian.registry || throw(ArgumentError("Use the same Pauli registry as the plan."))
        reason, data = _guarded_objective(plan, H)
    end
    rebuild_seconds = 0.0
    if reason !== :reused
        rebuild_seconds = @elapsed plan = compile_relaxation(H; plan.options...)
        data = plan.objective_data
    end
    numeric = _solve_numeric(plan, data, sense, optimizer)
    floor = _pauli_floor(H)
    status = isnothing(numeric.bound) ? :unavailable : :numerical_estimate
    provenance = _bound_provenance(plan.provenance, status, floor)
    certificate_seconds = @elapsed begin
        certificate, certificate_status, pruning = certify ?
            _lift_ground_certificate(plan,H,numeric.bound,optimizer,budget) : (nothing,:not_requested,nothing)
    end
    certified_bound = isnothing(certificate) ? nothing : certificate.bound
    certificate_provenance = isnothing(certificate) ? nothing :
        BoundProvenance(:certified_lower_bound,true,provenance.exact_wrt_dense_npa,provenance.witnesses,
            provenance.active_relaxations,provenance.active_strengthenings,provenance.sector_targets,floor)
    timings = (; guard_seconds, rebuild_seconds, numeric.setup_seconds, numeric.update_seconds,
        numeric.optimizer_seconds, certificate_seconds, total_seconds=(time_ns() - total_start)/1e9)
    result = (; numeric.bound, numeric.status, numeric.raw_status, numeric.moments, provenance,
        plan, reuse=reason, numeric.solver_model_reused, sense, timings, certified_bound,
        certificate, certificate_status, certificate_provenance, pruning)
    return result
end

_bound_provenance(p, status, floor) = BoundProvenance(status, false, p.exact_wrt_dense_npa,
    p.witnesses, p.active_relaxations, p.active_strengthenings, p.sector_targets, floor)

function _optimizer_factory(optimizer)
    isnothing(optimizer) && throw(ArgumentError("Supply a JuMP optimizer factory or optimizer_with_attributes."))
    optimizer isa NC.JuMPBackend && return optimizer.optimizer
    optimizer isa NC.AbstractSolver && throw(ArgumentError("The prototype accepts JuMP optimizer factories, not native backends."))
    return optimizer
end

# Only lift an intact, unprojected moment cone. Other routes need their own
# adjoints (including class multiplicity and term-sparsity auxiliary terms).
function _lift_ground_certificate(plan,H,bound,optimizer,budget)
    isnothing(bound) && return nothing,:no_optimum,nothing
    length(plan.reports.original_rows)==1 && length(plan.program.blocks)==1 &&
        !NC._has_symmetry(plan.symmetry) && !plan.options.wang_matched &&
        plan.reports.term_sparsity.algorithm in (:None,:Maximal) || return nothing,:lift_unavailable,nothing
    rows = only(plan.reports.original_rows)
    n = length(rows)
    block = only(plan.program.blocks)
    block.label==0 && block.kind==NC.PSD && block.side in (n,2n) || return nothing,:lift_unavailable,nothing
    if plan.reports.term_sparsity.algorithm === :Maximal
        cliques = plan.reports.term_sparsity.cliques
        length(cliques)==1 && length(only(values(cliques)))==1 &&
            only(only(values(cliques))) == collect(1:block.side) || return nothing,:lift_unavailable,nothing
    end
    prepared = plan.cache[optimizer]
    model = prepared.model
    side = block.side
    if plan.moi_form === :dualized
        JuMP.has_duals(model) || return nothing,:no_dual_point,nothing
        constraints = JuMP.all_constraints(model; include_variable_in_set_constraints=false)
        length(constraints)==1 || return nothing,:lift_unavailable,nothing
        D = Matrix(JuMP.dual(only(constraints)))
    else
        # SOS variables are svec scaled, in lower-column order.
        z = JuMP.value.(prepared.variables)
        D = [z[1+NC.tri_offset(side,max(i,j),min(i,j))] / (i==j ? 1.0 : sqrt(2.0)) for i in 1:side, j in 1:side]
    end
    all(isfinite,D) || return nothing,:no_dual_point,nothing
    # The adjoint of [Re M -Im M; Im M Re M], conjugated to match
    # the ambient verifier's sum G[i,j]*P[i]'*P[j] convention.
    gram = side==n ? D :
        D[1:n,1:n] + D[n+1:end,n+1:end] + im*(D[1:n,n+1:end] - D[n+1:end,1:n])
    eig = eigen(Hermitian((gram+gram')/2))
    factors = sqrt.(max.(eig.values,0.0)) .* eig.vectors'
    squares = [_wf_polynomial(H.registry,eltype(rows),((factors[i,j],rows[j]) for j in 1:n))
        for i in eachindex(eig.values) if eig.values[i] > 0]
    try
        if isnothing(budget)
            return certify_rational_sos(H,bound; squares),:verified_rational,nothing
        end
        pruning = prune_sos_factors(H,bound; squares,loss_budget=budget)
        return pruning.certificate,:verified_rational,pruning
    catch error
        error isa ArgumentError || rethrow()
        return nothing,:lift_verification_failed,nothing
    end
end

"""
    bound_observable_interval(plan, H, O, energy_lb, energy_ub; optimizer)
    bound_observable_interval(plan, H, O, window::PauliEnergyWindow; optimizer)

Solve both ends of an observable relaxation over arbitrary states in the
given energy window and sector. A verified window is checked against H and
the target; raw endpoints are assumptions. The plan must not carry additional
ground-state or window restrictions. Solver endpoints remain numerical. Use
`certify_observable_interval` with ambient SOS witnesses to certify endpoints.
"""
function bound_observable_interval(plan::CompiledSDPPlan,H,O,window::PauliEnergyWindow; optimizer=nothing)
    _check_hamiltonian(H)
    lb,ub = _wf_window_bounds(window,_wf_exact(H),plan.reports.analysis.target)
    return _observable_interval(plan,H,O,Float64(lb,RoundDown),Float64(ub,RoundUp),optimizer,:verified_window)
end

function bound_observable_interval(plan::CompiledSDPPlan,H,O,lb::Real,ub::Real; optimizer=nothing)
    isfinite(lb) && isfinite(ub) && lb <= ub || throw(ArgumentError("Require finite ordered energy endpoints."))
    return _observable_interval(plan,H,O,Float64(lb,RoundDown),Float64(ub,RoundUp),optimizer,:assumed_window)
end

function _observable_interval(plan,H,O,lb,ub,optimizer,window_status)
    isfinite(lb) && isfinite(ub) && lb <= ub || throw(ArgumentError("Require finite ordered energy endpoints."))
    any(rule -> rule.kind in (:stationarity,:kms_psd,:energy_window,:centered_second_moment),
        plan.reports.analysis.rules) && throw(ArgumentError(
        "Observable intervals cover arbitrary states in the supplied window. Compile a plan without stationarity, KMS, a preexisting energy window, or variance cuts."))
    optimizer = _optimizer_factory(optimizer)
    _check_hamiltonian(H)
    _check_hamiltonian(O)
    H.registry === plan.hamiltonian.registry === O.registry || throw(ArgumentError("Use the plan's Pauli registry."))
    reason,energy = _guarded_objective(plan,H)
    # Constraining H leaves the eligible state class unchanged; changing O
    # can break averaging and also changes the Stage 4 support seed.
    if reason !== :reused || !(O == plan.options.observable)
        plan = compile_relaxation(H; merge(plan.options,(;observable=O))...)
        energy = plan.objective_data
    end
    data = _coordinate_data(O,plan.table)
    all(in(plan.constraint_ids),keys(data)) || throw(ArgumentError("The observable still needs unconstrained moments; increase order or change supports."))
    p = plan.standard
    entries = [(id-1,c) for (id,c) in energy if id != 1]
    offset = get(energy,1,0.0)
    extra = sparse(repeat([1,2];inner=length(entries)),repeat(first.(entries),2),
        vcat(-last.(entries),last.(entries)),2,length(p.c))
    cones = [p.cones; NC.DataCone(:nonnegative,2,length(p.b)+1,2,Int32(0))]
    standard = NC.MomentStandardForm([p.A;extra],[p.b;offset-lb;ub-offset],p.c,p.offset,cones)
    sos = plan.moi_form === :primalized ? NC.sos_standard_form(standard) : nothing
    numeric_plan = (; standard,sos,plan.moi_form,objective_data=plan.objective_data,cache=IdDict{Any,Any}())
    lower = _solve_numeric(numeric_plan,data,:min,optimizer)
    upper = _solve_numeric(numeric_plan,data,:max,optimizer)
    status = isnothing(lower.bound) || isnothing(upper.bound) ? :unavailable : :numerical_estimate
    return (; lower=lower.bound,upper=upper.bound,status,window_status,lower_status=lower.status,upper_status=upper.status)
end
