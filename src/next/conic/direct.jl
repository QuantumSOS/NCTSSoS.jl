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

`result.solution.diagnostics` includes `iterations`, `primal_residual` and
`dual_residual`: Mosek's INTPNT_ITER, INTPNT_PRIMAL_FEAS and INTPNT_DUAL_FEAS
information items. These describe the native formulation named by
`diagnostics.representation`; in SOS form primal/dual refer to the SOS task,
not the original moment problem. Mosek measures feasibility in its homogeneous
model, so these are not norms of the original problem's residuals. Native `-1`
values are retained when no interior-point measure is reported (for example,
after presolve); do not interpret them as small residuals. These are solver
reports, not certificates.

Keep `form=:auto` unless comparing formulations on a small problem; do not
force `:moment` on large PSD blocks. In the architect's XXZ4/order-2 review
(one Hermitian block of side 67, one thread), moment form hit a 20 s limit
or took 26 s, while SOS took 0.2–3 s.

An accurate bound does not imply accurate individual moments, even with a
unique ground state. On that dense XXZ ring, bound error was 7.5e-9 but the
largest checked moment error was 2.2e-5. Requested tolerances of 1e-10 gave
raw `STALL` and 3.8e-5 moment error; 1e-12 returned no bound. Native SCS had
2e-12 moment error in the earlier review. Validate observables separately;
tightening Mosek tolerances need not improve them. On the same case, default
JuMPBackend/Mosek bound error was 1.6e-6 versus about 8e-9 for native SOS,
so prefer native SOS as the reference there, not default JuMP/Mosek settings.
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

"""
    Next.task_sizes(prepared)

Read native Mosek task dimensions without solving, from the value returned by
Next.prepare or Next.prepare_standard with MosekBackend. Returns a named tuple:
`affine_rows` counts linear constraint rows plus affine conic expression rows;
`scalar_variables` excludes semidefinite bar-variable entries; `psd_sides` lists
real PSD matrix dimensions, first bar variables, then affine PSD cones.
These are loaded task sizes before presolve, including Hermitian realification,
not the word-family sizes or the pre-realification sides from Next.block_sides.
The returned vector is owned by the caller. Other backends are not supported.
"""
task_sizes(prepared) = throw(ArgumentError("task_sizes requires a native Mosek handle from Next.prepare or Next.prepare_standard with MosekBackend."))

function _check_append_only(old::MomentStandardForm,form::MomentStandardForm)
    m,k = length(old.b),length(old.cones)
    unchanged = form.c==old.c && form.offset==old.offset && size(form.A,2)==size(old.A,2) &&
        length(form.b)>=m && length(form.cones)>=k && form.b[1:m]==old.b &&
        form.A[1:m,:]==old.A && form.cones[1:k]==old.cones
    unchanged || throw(ArgumentError("An open solver handle only accepts appended constraints with the same objective and variables; start a fresh solve for other changes."))
    return nothing
end
