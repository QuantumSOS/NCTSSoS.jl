# Direct Mosek writer

`using Mosek` activates the weak-dependency extension. Pass `solver=MosekBackend(; options...)` to the same `solve` and `expectation` front door as SCS. The writer constructs native Mosek tasks; it does not construct JuMP models or call MosekTools. MosekTools appears only in comparison tests and benchmarks. There is no Mosek precompile workload.

`form=:auto` chooses SOS when the number of moment variables is no larger than the number of cone coordinates **after Hermitian realification**; otherwise it chooses moment form. This size heuristic compares the SOS equality dimension with the moment cone dimension, not a promised solve-time optimum. Users can force `:sos` or `:moment`. SOS uses native PSD matrix variables and recovers moments as negative equality duals; moment form uses affine conic constraints and reads primal moment variables. Both use `[Re(M) -Im(M); Im(M) Re(M)]` for Hermitian PSD cones.

Only `MSK_SOL_STA_OPTIMAL` returns a bound and moments. Certificates map to infeasible/unbounded statuses with the SOS direction reversed. Unaccepted time/iteration limits return `TIME_LIMIT`/`ITERATION_LIMIT` and no bound. Set `MSK_DPAR_OPTIMIZER_MAX_TIME` or `MSK_IPAR_INTPNT_MAX_ITERATIONS`; raw termination remains in diagnostics even when an optimal solution has a `STALL` termination. `solve(result, constraint)` retains the model and backend settings but rebuilds the native task; it does not claim a warm start.

## Dense Bell I_nn22, n=5, order 3

Same relaxation on both routes: one PSD block of side 436, 27,555 nonconstant real moment variables, 95,265 standard-matrix nonzeros. The legacy reference uses its documented `sos_dualize` route, not a different hierarchy or a symmetry-reduced model. Objective is the negative Bell polynomial, as in the scale probe.

| Route | Status | Bound | Timed solve s | Whole process s | Peak RSS GB |
|---|---|---:|---:|---:|---:|
| Native `solve(Problem(...), Relaxation(order=3); solver=MosekBackend())`, auto/SOS | OPTIMAL | −0.3056390766662601 | 254.229 | 288.787 | 11.552 |
| Current package `sos_dualize`, MosekTools | OPTIMAL | −0.30563907748681035 | 248.970 | 290.020 | 11.683 |

The absolute bound difference is 8.21×10⁻¹⁰, below 10⁻⁶. Native recovered moments give objective −0.30563907252751077, within 4.14×10⁻⁹ of its bound. Native termination is `MSK_RES_OK`. These are single acceptance runs, not a speedup claim. Native timed solve includes the front-door rebuild; the legacy timer covers `JuMP.optimize!` after its symbolic build and write. Whole-process times include loading, compilation and diagnostic build/write calls.

Both runs used 16 Mosek threads, one Julia/BLAS thread, and primal/dual/gap tolerances 10⁻⁸, with a 1,800 s solver limit and a 20 GB RSS watchdog. Earlier native runs at 10⁻¹⁰ with one/eight threads hit a 600 s limit and returned no bound; they are not successful acceptance runs. The successful pair took 578.807 s (9.6 min), below the 20–60 minute estimate. Do not attribute the improvement to thread count alone: tolerance and available solve time changed too.

After rebasing on 0a48d3c, a second pair captured Mosek's interior-point feasibility measures and objectives explicitly. Both again returned `OPTIMAL`, and both bounds were bit-for-bit identical to the table above. The acceptance target remains an absolute bound difference below 10⁻⁶; it was not changed.

| Route | Solver primal feasibility | Solver dual feasibility | Absolute primal–dual objective gap |
|---|---:|---:|---:|
| Native SOS | 5.71268×10⁻⁸ | 9.23600×10⁻⁹ | 4.13875×10⁻⁹ |
| Legacy SOS | 1.78650×10⁻⁸ | 1.01615×10⁻⁸ | 4.14171×10⁻⁹ |

These are Mosek's reported measures for the solver's SOS representation, not independently normalized feasibility certificates; do not interpret the requested 10⁻⁸ tolerances as a guarantee that each absolute measure is at most 10⁻⁸. The second pair took 330.954 / 278.411 s whole process, 233.254 / 239.257 s timed solve, and 11.290 / 11.568 GB peak RSS (native / legacy). Native process time includes package/extension precompilation following the rebase. The initial paired timings above remain the comparison; neither pair establishes a speedup.

The previous moment-form `JuMPBackend` memory failure is investigated separately in [the lowering profile](mosek_lowering.md). No lowering fix was justified.

## The three tight-tolerance physics programs

All three front-door calls returned `OPTIMAL`, with one Mosek thread and requested primal/dual/gap tolerances of 10⁻¹⁰. These are Mosek-reported feasibility measures for the written SOS form, not independently normalized residual certificates. In particular, not every reported absolute measure is below 10⁻¹⁰.

| Program | Bound | Solver primal feasibility | Solver dual feasibility | Raw termination |
|---|---:|---:|---:|---|
| XXZ4, order 2, translation, DensityMatrix(1:4) | −6.744562646739805 | 7.63×10⁻¹² | 1.34×10⁻⁷ | OK |
| TFI6, order 2, site groups, term split, translation, StateOptimality | −6.842601652841712 | 5.38×10⁻¹¹ | 1.80×10⁻⁷ | OK |
| Heisenberg8, order 2, translation/reflection, no continuous symmetry | −14.60683242593556 | 3.38×10⁻¹⁰ | 3.90×10⁻⁹ | STALL |

Heisenberg's recovered moment objective is −14.606832394119962, differing by 3.18×10⁻⁸ from the bound. The direct writer therefore does **not** eliminate the stall at a requested 10⁻¹⁰ tolerance. It exposes both Mosek's `OPTIMAL` solution status and raw `STALL` termination. The three cases ran in 54.32 s (estimate roughly two minutes), including loading and compilation; their timed front-door calls after diagnostic preparation were 0.435, 0.165 and 1.251 s. The raw physics log is included with the measurements.

## Warm build plus write, not solve time

Five paired samples, alternating route order, after warmup; input construction and GC outside the timer. Both routes use one Mosek thread and identical 10⁻¹⁰ primal/dual/gap tolerances. Each small model was solved once per route to check its independently known bound before timing. Native uses automatic SOS form; JuMP writes moment form and attaches MosekTools. These are the backend choices under comparison, not identical representations. Julia allocations exclude native C memory.

| Model | Native median ms | JuMP median ms | Native bytes | JuMP bytes | Native allocations | JuMP allocations |
|---|---:|---:|---:|---:|---:|---:|
| Heisenberg4 | 0.596 | 6.207 | 154,464 | 914,776 | 311 | 17,806 |
| Hubbard2 (four modes) | 1.302 | 8.218 | 993,944 | 6,166,368 | 388 | 104,774 |
| CHSH | 0.491 | 4.613 | 15,248 | 172,752 | 218 | 5,418 |

[Raw samples, quartiles, solver output and Bell logs](mosek_results.json) preserve the measurements. The writer batch took 93.76 s, within the 1–3 minute estimate. A preceding trial stopped at its oracle check because the JuMP comparator used default solver settings; no failed-trial timing is included. The retained harness sets identical tolerances and thread counts on both routes.

## Validation

The combined Next, SCS writer, native Mosek, Aqua and ExplicitImports gate passed **39,211 assertions** in 705.07 s, within the 10–15 minute estimate. It used base d97ca7f before rebasing on the workload restoration. The standalone native writer suite contains 27 construction checks and 120 licensed checks, covering both forced forms and automatic selection, real/Hermitian PSD, SOC, linear and constant programs, equalities, appended constraints, infeasibility, time and iteration limits, and moment recovery. Small-model bounds match JuMPBackend with Mosek within 10⁻⁷. The standalone expanded suite and ExplicitImports took 109.29 s.

After rebase on 0a48d3c, the 147 native checks, Aqua 11/11 and ExplicitImports passed again in 81.45 s. With an empty `HOME` and no `MOSEKLM_LICENSE_FILE`, construction passed 27/27 and licensed solves printed the intended skip message. That process took 308.48 s because the isolated environment rebuilt dependency caches; the construction testset itself took 2.0 s. The post-rebase Bell pair plus these two checks took 999.30 s (16.7 min), exceeding the 10–15 minute estimate. The examples stamp and `git diff --check` also passed.

CI installs Mosek as a weak dependency, loads the extension and runs construction without a license. Licensed tests print an explicit skip message when neither the standard license file nor `MOSEKLM_LICENSE_FILE` is available; an invalid configured license is a failure, not a silent skip.

## Reproduction

Julia 1.12.7, Mosek 11.2.1, MosekTools 0.15.11, JuMP 1.32.1; Linux x86_64, Intel Xeon 2.60 GHz (`icelake-server`), shared orb, no pinning. The initial measurements use d97ca7f with this native-writer implementation; the explicitly labeled post-rebase pair uses 0a48d3c. Run numerical jobs serially, in a disposable environment developing this checkout and containing Mosek, MosekTools, SCS, JuMP, Aqua and ExplicitImports. Licensed solves verified the available license; construction checks also run without one.

```sh
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --project=PATH_TO_ENV --startup-file=no test/next/conic/mosek.jl
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  julia --project=PATH_TO_ENV --startup-file=no perf/next/mosek_writers.jl
MOSEK_THREADS=16 MOSEK_TOLERANCE=1e-8 MOSEK_MAX_TIME=1800 \
  OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=16 \
  julia --project=PATH_TO_ENV --startup-file=no perf/next/mosek_acceptance.jl bell auto
```

Repeat the last command with `bell legacy_sos` for the current-package reference. `xxz auto`, `tfi auto` and `heisenberg auto` select the three small acceptance programs; their default tolerance is 10⁻¹⁰. Use an external RSS watchdog for large runs: the solver time limit is not a memory limit. The measured watchdog sampled `/proc/PID/status` every 100 ms; peak RSS includes the whole process.
