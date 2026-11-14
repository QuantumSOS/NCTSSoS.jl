"""
Native SCS moment-form writer. Load SCS before solving. Set
`linear_solver=SCS.IndirectSolver` for the indirect linear system solver;
the default is SCS.DirectSolver. Other keywords are passed to SCS.
"""
struct SCSBackend{L,O<:NamedTuple} <: AbstractSolver
    linear_solver::L
    options::O
end
SCSBackend(;linear_solver=nothing,kwargs...) = SCSBackend(linear_solver,(;kwargs...))

"Construct native solver data without solving it."
function prepare_standard(form,backend::SCSBackend)
    throw(ArgumentError("No SCS writer for $(typeof(form)); load SCS with `using SCS` and pass a MomentStandardForm."))
end
solve_standard(form,backend::SCSBackend) = solve_standard(prepare_standard(form,backend))

function _check_append_only(old::MomentStandardForm,form::MomentStandardForm)
    m,k = length(old.b),length(old.cones)
    unchanged = form.c==old.c && form.offset==old.offset && size(form.A,2)==size(old.A,2) &&
        length(form.b)>=m && length(form.cones)>=k && form.b[1:m]==old.b &&
        form.A[1:m,:]==old.A && form.cones[1:k]==old.cones
    unchanged || throw(ArgumentError("An open solver handle only accepts appended constraints with the same objective and variables; start a fresh solve for other changes."))
    return nothing
end
