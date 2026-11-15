# Loading Mosek invalidates compiled SCS paths

This is a diagnosis, not a performance fix. The effect also occurs on Pauli
problems, not only on functional words. No solver registration, dependency,
precompile workload or numerical tolerance was changed to hide it.

## Reproduction

Run `backend_first_call.jl` in a fresh Julia process for each domain and loading
mode, three times each, with no other Julia job running:

```sh
julia --startup-file=no --compiled-modules=existing --project=<solver-environment> perf/next/backend_first_call.jl state scs
julia --startup-file=no --compiled-modules=existing --project=<solver-environment> perf/next/backend_first_call.jl state both
julia --startup-file=no --compiled-modules=existing --project=<solver-environment> perf/next/backend_first_call.jl pauli scs
julia --startup-file=no --compiled-modules=existing --project=<solver-environment> perf/next/backend_first_call.jl pauli both
```

The timed region is `solve` through `Base.invokelatest`, after package loading
and model construction. It includes any compilation needed by that first
solve. The second call uses the same problem, relaxation and backend. The
script checks the covariance bound against −5 and the four-site periodic
Heisenberg bound against −8, within 1e-5. Both use order 2 and native SCS with
`eps_abs=eps_rel=1e-6`. Package-cache generation and import time are excluded.

Environment: Julia 1.12.7, SCS.jl 2.7.1, Mosek.jl 11.2.1; Linux x86-64 orb,
Intel Xeon at 2.60 GHz (Ice Lake), one Julia thread and one BLAS thread. The
host is virtualized, not an exclusive physical node. The source is based on
redesign 28fc773 with the Newton-chip helper added.

## Three fresh-process samples per case

Seconds below are medians, with minimum–maximum in parentheses. Compilation
time is measured by Julia's `@timed`; all second-call compilation times were
zero. These are diagnostic samples, not a claimed optimization speedup.

| Domain | Loaded | First solve s | First compilation s | Second solve s |
| --- | --- | ---: | ---: | ---: |
| State | SCS | 0.214 (0.213–0.243) | 0.00494 (0.00486–0.00561) | 0.238 (0.206–0.241) |
| State | SCS + Mosek | 0.342 (0.336–0.351) | 0.131 (0.131–0.145) | 0.240 (0.205–0.256) |
| Pauli | SCS | 0.394 (0.390–0.478) | 0.315 (0.309–0.378) | 0.0773 (0.0766–0.102) |
| Pauli | SCS + Mosek | 0.716 (0.635–0.735) | 0.625 (0.558–0.627) | 0.0938 (0.0762–0.102) |

Every solve returned −5.0000023487687795 (state) or −7.999999999927017
(Pauli), independent of the loading mode. The state first-call penalty is
smaller here than the architect's 0.8 s observation, but the compilation
difference is reproduced. The Pauli first call already needs compilation
with SCS alone; loading Mosek increases it. The second-call ranges overlap.
The full 12-process measurement took 75.24 s, within the 1–2 minute estimate.

## Invalidation evidence

With NCTSSoS and SCS already loaded, enable Julia's diagnostic invalidation
log around `using Mosek`:

```julia
events = ccall(:jl_debug_method_invalidation, Any, (Cint,), 1)
@eval using Mosek
ccall(:jl_debug_method_invalidation, Any, (Cint,), 0)
foreach(println, events)
```

The log includes the insertion of this method from Mosek's `msk_enums.jl:1505`:

```julia
abstract type MosekEnum <: Integer end
Base.convert(::Type{Int32}, c::T) where {T <: MosekEnum} = c.value
```

It invalidates the existing `convert(::Type{Int32}, ::Integer)` specialization
and its callers. The invalidation trees include standard-form lowering and
precompiled `Next.solve` methods for both `PauliWord` and `FunctionalWord`.
The log also records invalidations from enum iteration methods. This is
evidence of dependency-method invalidation, not merely inference from elapsed
times. The insertion occurs in Mosek.jl, not our native-writer registration.

The diagnostic uses a Julia internal API and is not a package dependency.
It was run separately from the measurements. No small, verified local fix was
found; do not infer that adding another warmup or precompile call removes the
underlying invalidation.

## Deferred functional-word kernel work

The architect recorded the degree-five, eight-byte functional-word product
cost as a later optimization item. The existing `functional.md` measurements
are 0.956 µs for state and 1.137 µs for trace, with zero allocated bytes.
Those are earlier measurements, not new numbers from this investigation.
No functional-word kernel optimization is included in the Newton-chip change.
