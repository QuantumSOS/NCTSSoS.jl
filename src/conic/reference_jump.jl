# Reference lowering of a `MomentProgram` to JuMP, in moment form.
#
# This is the correctness reference for contract 2 and the JuMP export for
# solvers without a direct writer. It is not a fast path: it builds one JuMP
# expression for each form.

function _jump_form(model, y, forms::LinearForms{C}, k::Integer) where {C}
    expr = zero(JuMP.GenericAffExpr{C,JuMP.VariableRef})
    for j in term_range(forms, k)
        id = forms.ids[j]
        if id == 1
            JuMP.add_to_expression!(expr, forms.coefs[j])
        else
            JuMP.add_to_expression!(expr, forms.coefs[j], y[id])
        end
    end
    return expr
end

"""
    to_jump(program::MomentProgram) -> (model, y)

Build a JuMP model of `program` in moment form. `y[id]` is the JuMP variable of
moment `id` for `id >= 2`. The caller sets the optimizer and solves. After a
solve, `moment_values(program, y)` gives the value of each moment.
"""
function to_jump(p::MomentProgram{C}) where {C}
    validate(p)
    model = JuMP.Model()
    @variable(model, y[2:nvars(p)])
    f = p.forms
    for b in p.blocks
        r = form_range(b)
        if b.kind == ZERO
            for k in r
                e = _jump_form(model, y, f, k)
                @constraint(model, real(e) == 0)
                C <: Complex && @constraint(model, imag(e) == 0)
            end
        elseif b.kind == NONNEGATIVE
            for k in r
                @constraint(model, real(_jump_form(model, y, f, k)) >= 0)
            end
        elseif b.kind == SOC
            @constraint(model, [real(_jump_form(model, y, f, k)) for k in r] in JuMP.SecondOrderCone())
        else
            n = b.side
            E = typeof(_jump_form(model, y, f, first(r)))
            M = Matrix{E}(undef, n, n)
            for j in 1:n, i in j:n
                e = _jump_form(model, y, f, b.first_form + tri_offset(n, i, j))
                M[i, j] = e
                i == j || (M[j, i] = conj(e))
            end
            if b.kind == PSD
                @constraint(model, LinearAlgebra.Symmetric(real.(M)) in JuMP.PSDCone())
            else
                @constraint(model, LinearAlgebra.Hermitian(M, :L) in JuMP.HermitianPSDCone())
            end
        end
    end
    objective = zero(JuMP.AffExpr)
    for (id, c) in zip(p.objective_ids, p.objective_coefs)
        id == 1 ? JuMP.add_to_expression!(objective, c) : JuMP.add_to_expression!(objective, c, y[id])
    end
    @objective(model, Min, objective)
    return model, y
end

"""
    moment_values(program, y) -> Vector{Float64}

The solved value of each moment variable; entry 1 is the constant 1.
"""
function moment_values(p::MomentProgram, y)
    values = ones(Float64, nvars(p))
    for id in 2:nvars(p)
        values[id] = JuMP.value(y[id])
    end
    return values
end
