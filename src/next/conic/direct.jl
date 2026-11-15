public prepare_standard

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

"""
Native Mosek writer. Load `Mosek` before solving. `form=:auto` selects SOS
when there are no more moment variables than cone coordinates after Hermitian
realification, otherwise moment form. Override with `form=:moment` or `:sos`. Other keywords are
Mosek's `MSK_IPAR_`, `MSK_DPAR_`, or `MSK_SPAR_` parameters.
Use `MSK_DPAR_OPTIMIZER_MAX_TIME` for a time limit and
`MSK_IPAR_INTPNT_MAX_ITERATIONS` for an interior-point iteration limit.
Only Mosek's `OPTIMAL` solution status supplies a bound and moments; an
unaccepted solution at a limit reports `TIME_LIMIT` or `ITERATION_LIMIT`.
The raw termination code is retained in the solution diagnostics, including
`STALL` when Mosek nevertheless reports an optimal solution. Both forms
support `expectation`; SOS recovers moments from the equality duals.
Appended constraints rebuild the native task; no warm start is claimed.
"""
struct MosekBackend{O<:NamedTuple} <: AbstractSolver
    form::Symbol
    options::O
end
function MosekBackend(;form=:auto,kwargs...)
    form in (:auto,:moment,:sos) || throw(ArgumentError("Mosek form must be :auto, :moment, or :sos."))
    return MosekBackend(form,(;kwargs...))
end
function prepare_standard(form,backend::MosekBackend)
    throw(ArgumentError("Load Mosek with `using Mosek` and pass a MomentStandardForm to MosekBackend."))
end
solve_standard(form,backend::MosekBackend) = solve_standard(prepare_standard(form,backend))

function _check_append_only(old::MomentStandardForm,form::MomentStandardForm)
    m,k = length(old.b),length(old.cones)
    unchanged = form.c==old.c && form.offset==old.offset && size(form.A,2)==size(old.A,2) &&
        length(form.b)>=m && length(form.cones)>=k && form.b[1:m]==old.b &&
        form.A[1:m,:]==old.A && form.cones[1:k]==old.cones
    unchanged || throw(ArgumentError("An open solver handle only accepts appended constraints with the same objective and variables; start a fresh solve for other changes."))
    return nothing
end
