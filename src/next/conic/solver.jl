"Solver descriptions dispatch at the standard-form boundary; package extensions can add new subtypes."
abstract type AbstractSolver end
"Use any JuMP optimizer factory or optimizer_with_attributes value. No solver is loaded by this description."
struct JuMPBackend{O} <: AbstractSolver
    optimizer::O
end

mutable struct JuMPHandle{F}
    model::JuMP.Model
    variables::Vector{JuMP.VariableRef}
    form::F
end
"Numerical relaxation optimum and status, not a rigorous certificate. A failed solve has bound=nothing."
struct ConicSolution{H,S,D}
    bound::Union{Nothing,Float64}
    status::S
    raw_status::String
    values::Vector{Float64}
    handle::H
    diagnostics::D
end
ConicSolution(bound,status,raw_status,values,handle) = ConicSolution(bound,status,raw_status,values,handle,nothing)

"Solve standard data with a backend description, or re-solve an open backend handle."
function solve_standard end
function _optimize(handle::JuMPHandle)
    JuMP.optimize!(handle.model)
    status = JuMP.termination_status(handle.model)
    accepted = status == JuMP.MOI.OPTIMAL
    has_values = JuMP.has_values(handle.model)
    bound = accepted && has_values ? Float64(JuMP.objective_value(handle.model)) : nothing
    values = has_values ? JuMP.value.(handle.variables) : Float64[]
    return ConicSolution(bound,status,JuMP.raw_status(handle.model),values,handle)
end
function solve_standard(form::Union{MomentStandardForm,SOSStandardForm},solver::JuMPBackend)
    model,variables = to_jump(form)
    JuMP.set_optimizer(model,solver.optimizer)
    return _optimize(JuMPHandle(model,variables,form))
end
function solve_standard(form::MomentStandardForm,handle::JuMPHandle{MomentStandardForm})
    old = handle.form
    m,k = length(old.b),length(old.cones)
    unchanged = form.c==old.c && form.offset==old.offset && size(form.A,2)==size(old.A,2) &&
        length(form.b)>=m && length(form.cones)>=k && form.b[1:m]==old.b &&
        form.A[1:m,:]==old.A && form.cones[1:k]==old.cones
    unchanged || throw(ArgumentError("Re-solving an open handle only permits appended constraints with the same objective and variables; call solve_standard(form, backend) for a changed problem."))
    if length(form.b)>m
        coordinates = form.b[m+1:end]-form.A[m+1:end,:]*handle.variables
        cones = [DataCone(c.kind,c.side,c.first-m,c.length,c.label) for c in form.cones[k+1:end]]
        _jump_cones!(handle.model,coordinates,cones)
    end
    handle.form = form
    return _optimize(handle)
end
