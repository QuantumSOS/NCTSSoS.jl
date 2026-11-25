using LinearAlgebra: Symmetric, dot, eigmin, norm

function square_psd_residual(v,cone)
    M = Matrix{Float64}(undef,cone.side,cone.side)
    k = cone.first
    for j in 1:cone.side, i in j:cone.side
        M[i,j] = M[j,i] = v[k] / (i==j ? 1.0 : sqrt(2.0))
        k += 1
    end
    max(0.0,-eigmin(Symmetric(M)))
end

function square_jump_points(form, model, variables)
    jump = NCTSSoS.JuMP
    x = jump.has_values(model) ? jump.value.(variables) : nothing
    z = nothing
    if jump.has_duals(model)
        # to_jump emits one constraint per PSD block and scalar equality, in
        # the same per-type order as the standard form. JuMP.dual restores
        # the symmetric shape; convert it back to Frobenius-isometric svec.
        psd = jump.all_constraints(model,Vector{jump.AffExpr},jump.MOI.PositiveSemidefiniteConeTriangle)
        equalities = jump.all_constraints(model,jump.AffExpr,jump.MOI.EqualTo{Float64})
        inequalities = jump.all_constraints(model,jump.AffExpr,jump.MOI.GreaterThan{Float64})
        z = zeros(length(form.b))
        p = e = n = 0
        for c in form.cones
            if c.kind == :zero
                for k in c.first:c.first+c.length-1
                    z[k] = jump.dual(equalities[e+=1])
                end
            elseif c.kind == :nonnegative
                for k in c.first:c.first+c.length-1
                    z[k] = jump.dual(inequalities[n+=1])
                end
            else
                D = jump.dual(psd[p+=1])
                k = c.first
                for j in 1:c.side, i in j:c.side
                    z[k] = D[i,j]*(i==j ? 1.0 : sqrt(2.0))
                    k += 1
                end
            end
        end
        @assert p == length(psd) && e == length(equalities) && n == length(inequalities)
    end
    return x,z
end

function square_original_diagnostics!(result, form, model, variables)
    return square_point_diagnostics!(result,form,square_jump_points(form,model,variables)...)
end

function square_point_diagnostics!(result, form, x, z)
    all(c -> c.kind in (:zero,:psd,:nonnegative),form.cones) || error("Unsupported square diagnostic cone")
    function publish(key,value)
        isfinite(value) && (result[key] = value)
    end
    if !isnothing(x) && all(isfinite,x)
        slack = form.b - form.A*x
        publish("original_primal_objective_per_site",form.offset + dot(form.c,x))
        publish("original_equality_inf",maximum((abs(slack[i]) for c in form.cones if c.kind==:zero
            for i in c.first:c.first+c.length-1);init=0.0))
        publish("original_primal_nonnegative_violation",maximum((-slack[i] for c in form.cones if c.kind==:nonnegative
            for i in c.first:c.first+c.length-1);init=0.0))
        all(isfinite,slack) && publish("original_primal_psd_violation",maximum(
            (square_psd_residual(slack,c) for c in form.cones if c.kind==:psd);init=0.0))
    end
    if !isnothing(z) && all(isfinite,z)
        publish("original_dual_objective_per_site",form.offset-dot(form.b,z))
        publish("original_balance_inf",norm(transpose(form.A)*z+form.c,Inf))
        publish("original_dual_nonnegative_violation",maximum((-z[i] for c in form.cones if c.kind==:nonnegative
            for i in c.first:c.first+c.length-1);init=0.0))
        publish("original_dual_psd_violation",maximum((square_psd_residual(z,c) for c in form.cones if c.kind==:psd);init=0.0))
    end
    if haskey(result,"original_primal_objective_per_site") && haskey(result,"original_dual_objective_per_site")
        p = result["original_primal_objective_per_site"]
        d = result["original_dual_objective_per_site"]
        publish("original_objective_relative_gap",abs(p-d)/max(1,abs(p),abs(d)))
    end
    return result
end

function square_moment_solve!(result, path, form)
    all(c -> c.kind in (:psd,:zero),form.cones) || error("Expected real square PSD and equality cones")
    result["prototype_numeric_form"] = "explicit_moment"
    save_result(path,result)
    result["solve_seconds"] = @elapsed begin
        result["jump_model_seconds"] = @elapsed model,variables = NCTSSoS.to_jump(form)
        optimizer = Ref{Any}(nothing)
        JuMP.set_optimizer(model,() -> (optimizer[] = MOI.instantiate(OPTIMIZER)))
        save_result(path,result)
        result["optimize_wall_seconds"] = @elapsed JuMP.optimize!(model)
        result["optimizer_reported_seconds"] = MOI.get(optimizer[],MOI.SolveTimeSec())
        accepted = record_status!(result,optimizer[])
        # Supplemental original-coordinate diagnostics, never substitutes
        # for the existing status/feasibility/relative-gap acceptance gate.
        result["original_diagnostics_seconds"] = @elapsed square_original_diagnostics!(result,form,model,variables)
    end
    return accepted ? JuMP.objective_value(model) : nothing
end
