# Packed real state and trace polynomials

## Scope and reproduction

Core delivery (a): Free, Projector and Unipotent algebras, explicitly real
functionals, complete mixed rows, localizing constraints and outer-moment reads.
Newton-chip selection is delivery (b). Grouping, split, finite symmetry and
specialized strengthening remain named errors pending pipeline delivery (c).

Run `julia --project=<environment> perf/next/functional.jl cost`. No solver is
needed for the cost run. This measurement used Julia 1.12.7 on an icelake-server
Linux orb, one Julia thread and one BLAS thread. The old package is the reference
in the same checkout, based on redesign commit 79ac93d. There are three samples
after an untimed warmup, with GC before each sample. Times below are medians;
ranges are the minimum and maximum, not confidence intervals. Bytes and
allocation counts are identical across each set of three samples.

The build workload is the benchmark `covinn22`, n=3, order=2, with sign matrix
`[1 1 -1; 1 1 1; -1 1 0]`, minimizing the negative covariance expression. State
variables use two commuting groups of three generators; trace variables use
one group of six. Both use complete rows without symmetry or strengthening.

## Construction costs

| Domain / phase | New ms (range) | Old ms (range) | New bytes / allocations | Old bytes / allocations |
| --- | ---: | ---: | ---: | ---: |
| State model | 0.234 (0.200–0.284) | 0.145 (0.134–0.290) | 69,904 / 1,230 | 45,152 / 1,015 |
| State assembly | 7.815 (6.737–7.823) | 21.102 (18.885–22.429) | 268,256 / 124 | 42,566,728 / 1,100,258 |
| State SOS lowering | 0.439 (0.381–0.455) | 6.761 (5.982–6.766) | 819,832 / 101 | 8,236,144 / 128,972 |
| Trace model | 0.225 (0.224–0.267) | 0.161 (0.151–0.175) | 70,112 / 1,232 | 49,536 / 1,158 |
| Trace assembly | 8.963 (8.930–9.016) | 24.305 (24.240–27.167) | 289,504 / 125 | 55,531,744 / 1,470,408 |
| Trace SOS lowering | 0.504 (0.488–0.527) | 9.432 (8.578–10.817) | 1,320,664 / 109 | 13,995,496 / 254,284 |

Model construction alone regresses slightly; the larger assembly savings
outweigh it. These measurements do not include a numerical solve or peak RSS.
The new SOS lowering produces numeric sparse standard form; the old
`sos_dualize` produces the equivalent JuMP SOS model. The old state route has
no solver-free `build_jump_model` method for its primal form. Thus the lowering
comparison includes JuMP overhead and is not a comparison of identical output
containers. The new primal standard form separately takes 0.341 ms / 629,312
bytes (state), and 0.446 ms / 1,098,976 bytes (trace).

The old and new moment counts agree: state has 106 rows, 1,098 moments and
1,097 nonconstant real variables; trace has 115 rows, 1,010 moments and 1,009
nonconstant real variables. Both have one real PSD block. New primal matrix
sizes are 5,671 × 1,097 and 6,670 × 1,009, respectively.

## Width costs and the eight-limb limit

The payload contains operator tokens, scalar-factor endpoints and two lengths.
`word_type` selects the smallest limb count covering the required total degree;
the capacity for that limb count is derived from the layout. At four-bit tokens
the capacities of one, two, three and eight limbs are 11, 23, 36 and 99 letters.
An overflow is an `ArgumentError`; there is no vector fallback.

Each sample below is 8,000 calls over eight runtime inputs. A product writes
into a reserved polynomial buffer; closure writes into a preallocated word
vector. Every sample allocates **zero bytes and zero objects**. Inputs contain
one scalar factor and an operator tail, with three-bit tokens. Tail length grows
with capacity: this is a width-and-degree stress measurement, not an isolated
comparison of different widths on a fixed word. The allocation tests also
exercise every intermediate limb count.

| Limbs / bytes | Input degree | State product µs | State closure µs | Trace product µs | Trace closure µs |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 / 8 | 5 | 0.956 | 0.817 | 1.137 | 1.057 |
| 2 / 16 | 10 | 1.987 | 1.309 | 2.248 | 2.136 |
| 4 / 32 | 21 | 3.841 | 2.308 | 4.254 | 4.661 |
| 8 / 64 | 42 | 9.942 | 4.855 | 8.783 | 10.351 |

The eight-limb ceiling is the approved first implementation limit, not a
claim of an optimal cutoff. Trace closure includes cyclic normalization.
Times include buffer writes and loop overhead; product timing also consumes
the output hash to prevent dead-code elimination.

## Numerical acceptance

`test/next/model/functional.jl` exposes
`FunctionalModelChecks.check(backend; factory, references=true)` for the full
comparison. The native Mosek writer was checked against the old package using
MosekTools, after a trivial licensed LP succeeded. Default native settings and
the old solver factory's default tolerances were retained.

| Case | Native Mosek | Reference | Absolute difference |
| --- | ---: | ---: | ---: |
| Docs state covariance, order 2 | −5.000000000166828 | −4.999999121245736 | 8.79e-7 |
| Benchmark state covariance, n=3, order 2 | −5.000000000166828 | −4.999999116166029 | 8.84e-7 |
| Trace covariance, either sign matrix, order 2 | −4.9999999996301785 | −4.999999349100641 | 6.51e-7 |
| Projector trace, order 2 | −0.04671737627039185 | −0.046717378455438933 | 2.19e-9 |
| Projector trace, order 3 | −0.0312499934321917 | −0.03124998978001017 | 3.65e-9 |

The state comparison is checked at 1e-6 against the same old polynomial,
separately from the 1e-3 literature-value check. Projector references are the
stored relaxation values in `trace_poly_refs.json`, not asserted exact optima.
Trace CHSH is checked against −2√2 at 1e-5. Localizing scalar `Zero` and
`Positive` constraints are checked through solves; outer moments are also
checked against an explicit two-point ensemble whose first moment is zero
and whose second moment is one.

For two commuting trace groups, the old canonicalizer is not an oracle: it
can return different representatives for equivalent words and is not
idempotent. The new tests enumerate independent cyclic rotations of both
groups and a joint reversal, and check 100 random words against real matrix
traces. Another 100 random mixed-word pairs check products, adjoints and closure
against matrix products with a vector state and normalized trace. Single-group
products and adjoints and state products use term-by-term
legacy comparisons. All three base rules pass the word contract.

## First calls

The base precompile workload builds one state and one trace program. The SCS
extension solves both. The Mosek extension constructs both native tasks without
calling the optimizer; loading and precompilation do not require a license.
Fresh-process measurements use `perf/next/functional.jl first-scs`,
`first-mosek` and `first-old-mosek`. Each mode ran in three sequential fresh
processes after precompilation, without concurrent tests. `invokelatest` keeps
inference and code generation inside the timed calls rather than allowing
compilation of the benchmark driver to compile its callees first.

| Route / domain | Input construction s | First solve s (range) | Second solve s (range) | First solve bytes |
| --- | ---: | ---: | ---: | ---: |
| Native SCS / state | 0.203 | 0.267 (0.262–0.271) | 0.236 (0.226–0.251) | 1,322,880 |
| Native SCS / trace | 0.086 | 0.303 (0.295–0.374) | 0.452 (0.391–0.471) | 1,871,728 |
| Native Mosek / state | 0.173 | 0.297 (0.275–0.304) | 0.234 (0.222–0.269) | 6,503,808 |
| Native Mosek / trace | 0.095 | 0.276 (0.260–0.290) | 0.211 (0.189–0.255) | 3,289,904 |
| Legacy MosekTools / state | 6.001 | 18.887 (17.989–18.959) | 3.534 (3.263–3.746) | 2,121,801,880 |
| Legacy MosekTools / trace | 4.460 | 8.643 (8.299–8.768) | 4.807 (4.804–4.899) | 513,826,280 |

Values are medians. State runs before trace in each process: trace is its first
domain-specific call but can reuse common code compiled by state. The second
solve builds a fresh relaxation, rather than reusing a prepared solver handle.
The SCS trace timings are noisy; the slower second call is not a claimed
warmup benefit. Package startup, imports and driver compilation are outside
these solve timers. Full process wall time, including both domains and all
four solves, is 7.03 s (7.02–7.36) for native SCS, 6.89 s (6.73–7.05) for native
Mosek, and 51.78 s (51.45–51.97) for legacy MosekTools.

The solver environment uses SCS.jl 2.7.1, Mosek.jl 11.2.1 and MosekTools.jl
0.15.11. Mosek tolerance and thread settings are defaults on both routes.
Native SCS uses `eps_abs=eps_rel=1e-6`; its state and trace bounds are
−5.0000023488 and −4.9999993314. The native and legacy Mosek values agree with
the numerical acceptance table. These backend timings are not a solver
ranking: SCS and Mosek have different accuracy settings.

## Validation record

Commands run with `--startup-file=no` and an environment containing the local
package plus COSMO, SCS, Mosek, MosekTools, Aqua and ExplicitImports:

- `include("test/next/runtests.jl")`: all behavior testsets passed, followed by
  Aqua 11/11. The initial ExplicitImports run then rejected two non-public Base
  calls. Replaced them with the existing tuple update helper and iterator-length
  pattern, without changing the contract or adding an ignore entry.
- On that final implementation, `include("test/next/core/functional.jl")`:
  **9,782/9,782**. `FunctionalModelChecks.check` with native SCS and native
  Mosek: **33/33 each**. Aqua: **11/11**; ExplicitImports clean;
  `ExamplesStamp.verify_stamp()` passed. The full Next suite was not repeated
  after the two Base-call replacements; the changed kernels and both native
  model routes were rerun.
- `FunctionalModelChecks.check(MosekBackend(); factory,
  references=true)`: **14/14** reference checks, in addition to the model checks.
- `perf/next/functional.jl cost`: three samples per row above; all allocation
  assertions and legacy moment-count checks passed, 119.48 s.
- The nine final first-call processes took 197.28 s in total, within the
  3–5 minute estimate. The final targeted tests and quality run took 272.50 s,
  above its 2–4 minute estimate; the integrated run took 822.15 s, above its
  5–10 minute estimate. No claims are made about CI before it runs on the PR.
