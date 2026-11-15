# [Solver preparation](@id next-solver-preparation)

The staged path separates formulation, backend preparation and optimization:

```julia
using NCTSSoS.Next, Mosek
import NCTSSoS.Next as Next

_, (_, _, z) = pauli_variables(1)
problem = Problem(z[1])
relaxation = Relaxation(order=1)
solver = MosekBackend(form=:auto)
built = Next.build(problem, relaxation)
prepared = Next.prepare(built; solver)
Next.task_sizes(prepared)
result = solve(built; solver)
@assert isapprox(result.bound, -1; atol=1e-6)
```

`Next.build` does not solve. `Next.prepare` prepares and caches the backend
task. `solve(built; solver)` uses that prepared task. See
[Reading a result](@ref next-result) for result fields.

```@docs
NCTSSoS.Next.build
NCTSSoS.Next.BuiltRelaxation
NCTSSoS.Next.prepare
```

## Backends and forms

The standard form is chosen by the backend option, not by a keyword on `solve`.
`MosekBackend(form=:auto)` selects a form automatically; `:moment` and `:sos`
override that choice. Both forms support `expectation`. `SCSBackend` is a
moment-form writer. `JuMPBackend` accepts an optimizer factory or an
`optimizer_with_attributes` value and loads no solver itself; it writes moment
form and has no `form` option.

```@docs
NCTSSoS.Next.MosekBackend
NCTSSoS.Next.SCSBackend
NCTSSoS.Next.JuMPBackend
```

## Native task sizes and diagnostics

`Next.task_sizes` reports the native Mosek task before presolve. It includes
Hermitian realification, unlike the polynomial-row and cone-side counts on
the built relaxation.

```@docs
NCTSSoS.Next.task_sizes
```

See [`solve`](@ref NCTSSoS.Next.solve) for the result fields, program sizes and backend diagnostics.
