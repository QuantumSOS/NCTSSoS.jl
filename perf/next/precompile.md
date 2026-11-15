# Precompile scope and first-call cost

**Selected by the owner on 2026-11-15: restore the original six workload blocks.** His decision was: "Go back to the 68 s workload; first-call time matters more than precompile". The package builds and lowers the original Pauli lattice/grouping/split, CHSH, Hubbard, free-ball, DensityMatrix and InvariantCone cases; the SCS extension solves the original three groups of cases. Symmetry declarations use tuples or qualified `Next.Symmetry`; the names remain public but unexported. The 45 s and 28 s workloads below remain measured alternatives, not the selected scope.

## Selected restoration and historical choices on the same five probes

Times are seconds, medians of three isolated processes per configuration. Precompile totals are the median of the per-repetition package-plus-extension sums; separately computed component medians need not sum to that total. Ranges are minimum–maximum. First-call times measure `solve` after input construction, not loading or input construction; those additional times and second calls are retained in [raw measurements](precompile.json).

| Workload | Package s | SCS extension s | Total s [min,max] | xxz6 | xxz12front | hubbard | ball | CHSH |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Old, all six blocks | 60.523 | 7.494 | 68.016 [65.490,74.946] | 0.954 | 0.022 | 0.103 | 0.002 | 0.001 |
| Restored, selected | 67.787 | 9.460 | 76.976 [75.334,78.745] | 0.972 | 0.028 | 0.089 | 0.002 | 0.002 |
| None | 11.235 | 4.764 | 15.999 [14.734,16.241] | 18.149 | 25.459 | 3.963 | 4.615 | 3.384 |
| Common-path alternative | 38.793 | 6.834 | 44.965 [44.125,48.209] | 0.910 | 6.319 | 1.682 | 2.290 | 0.001 |
| Smallest useful: plain Pauli and CHSH | 20.946 | 6.807 | 27.734 [26.312,29.904] | 19.327 | 26.343 | 2.043 | 2.168 | 0.001 |

The smallest useful point meets the former 30 s target, but leaves the two XXZ first calls at 19–26 s. It removes only the translated build/lowering and solve from the common-path alternative, retaining the same four-site plain Pauli problem and CHSH. The 45 s alternative keeps all five first calls below 8 s in these samples. The owner chose the original scope for its shorter first calls instead.

| First solve | Old s [min,max] | Common-path alternative s [min,max] |
|---|---:|---:|
| xxz6 | 0.954 [0.801,0.997] | 0.910 [0.874,1.173] |
| xxz12front | 0.022 [0.022,0.023] | 6.319 [6.058,7.027] |
| hubbard | 0.103 [0.092,0.106] | 1.682 [1.663,1.923] |
| ball | 0.002 [0.002,0.002] | 2.290 [2.272,2.492] |
| CHSH | 0.001 [0.001,0.001] | 0.001 [0.001,0.001] |

All recorded first solves report `OPTIMAL`; their bounds agree exactly across configurations for each input. First-call times are printed to 0.001 s resolution. The xxz6 timing difference is within the observed ranges and is not a claimed speedup.

The restored first-solve ranges are xxz6 0.870–1.049 s, xxz12front 0.026–0.042 s, Hubbard 0.086–0.125 s, ball 0.002–0.003 s and CHSH 0.001–0.002 s. The five bounds exactly match the original workload: respectively −10.841281370523458, −21.68256252593265, −1.2360679663714502, −0.9999999771195256 and −2.82842712459816. **The strict timing-range acceptance is not met:** xxz12front is outside its old range and CHSH has a higher median at the timer's 0.001 s resolution. The multi-second cold calls are removed, but 76.976 s is not the historical 68.016 s precompile total. The measurements do not isolate the cause of this increase across source revisions and shared-orb runs; no workload scope was removed to hide it.

## Four paths outside the restored workload

Measurement only: none of these paths was added to a workload. Each row is the median of three fresh processes. Input includes model-builder compilation; loading is excluded. Detection is timed separately and subtracted from the input time. The second solve is a warm-call reference, not the first-call result. All twelve first solves returned `OPTIMAL`, with identical bounds across repetitions.

| Front-door path | Input s | Detection s | First solve s | Second solve s | Bound |
|---|---:|---:|---:|---:|---:|
| Heisenberg8, order 2, translation/reflection/spin rotations | 0.844 | — | 18.392 | 5.7725 | −14.60683422165745 |
| Fermionic4 hopping ring, order 1, P/Q/G, Sector(N,2), NumberRotation | 0.969 | — | 17.442 | 0.0103 | −1.9999999739663321 |
| Heisenberg8, detect_symmetry, then order 1 | 4.632 | 8.959 | 1.250 | 0.0107 | −24.00000000077552 |
| XXZ6 (anisotropy 0.7), order 1, translation, DensityMatrix(1:3) | 0.172 | — | 1.947 | 0.0096 | −10.841281401965553 |

The exact inputs are in `first_call_extra.jl`. The original five-case `first_call.jl` is unchanged. Separate builders prevent adding these measurement branches from increasing compilation of the original input builder. Raw per-repetition input, detection, solve and loading times are retained in `precompile.json`.

## Isolated cost of the original blocks

Only the named block is retained in each row. Both modules still precompile. Subtract the no-workload row to obtain the additional compile cost and the first-call saving for a probe. These costs and savings **overlap and must not be added**: a block can compile methods that another block needs.

| Enabled old block | Package s | Extension s | Total s | xxz6 | xxz12front | hubbard | ball | CHSH |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P1: Pauli lattice/grouping/split and CHSH builds | 45.222 | 5.394 | 50.616 | 1.310 | 0.496 | 2.169 | 2.532 | 0.575 |
| P2: Hubbard and free-ball builds | 17.321 | 4.815 | 22.070 | 18.464 | 24.030 | 0.558 | 0.521 | 2.089 |
| P3: DensityMatrix and InvariantCone builds | 42.280 | 5.358 | 47.586 | 4.125 | 4.228 | 2.610 | 2.745 | 2.478 |
| E1: translated/split Pauli and CHSH solves | 11.968 | 35.807 | 48.069 | 0.854 | 0.106 | 2.173 | 2.131 | 0.001 |
| E2: Hubbard and free-ball solves | 12.094 | 14.216 | 26.245 | 20.037 | 25.935 | 0.108 | 0.002 | 1.838 |
| E3: DensityMatrix and InvariantCone solves | 11.630 | 37.197 | 48.550 | 4.167 | 3.730 | 2.134 | 2.355 | 1.883 |

For example, isolated P2 adds 6.1 s over the no-workload total, while saving 3.4 s on the Hubbard first solve and 4.1 s on the free ball. These paths already take less than 5 s without a workload. The retained translated family avoids the much longer uncovered xxz6 and xxz12front calls.

### Rejected trials are not monotonic cost evidence

Removing `Correlative()` from the proposed translation case gave totals of 52.051, 52.600 and 58.016 s. Package components were 43.728, 45.296 and 49.408 s; extension components were 8.324, 7.304 and 8.608 s. A two-site dense example took 52.957 s (45.738 package + 7.218 extension); putting that example only in the extension took 54.715 s (15.010 + 39.705). The two-site trials have only one sample each and are not used for the four-choice comparison.

**The timing increase is not explained by these measurements.** Removing grouping is not simply deleting executed calls: it changes the row family and the compiled pipeline paths. `pipeline/families.jl` has different paths for one block and an orbit of several blocks. This is a possible contributor, not a measured attribution. Shared-orb variation is also visible in the controls. No conclusion that less workload inherently costs more is drawn from these trials.

## Reproduction and environment

Julia 1.12.7, SCS 2.7.1, JuMP 1.32.1; full package versions are in `results.json`. Linux x86_64, Intel Xeon 2.60 GHz, shared orb, one Julia/BLAS/OpenMP thread, no CPU pinning. Numerical jobs ran serially. Measurements use the old source [a50c459](https://github.com/QuantumSOS/NCTSSoS.jl/commit/a50c45902aab81db2b5e5c9142161b4592f6e3b1) and the candidate on [7c079ec](https://github.com/QuantumSOS/NCTSSoS.jl/commit/7c079ecf3ecfd919d606b15ce278fbfc393bc57d), whose intervening changes only affect documentation. The PR then rebased on [3fb021a](https://github.com/QuantumSOS/NCTSSoS.jl/commit/3fb021a) (combination tests and one sector error text); the measured workloads did not change. The probe is the architect's five-case first-call test (with CHSH added), kept identical in every configuration.

`precompile_probe.py` changes workload blocks only in a disposable source copy under `.amp/in/orb-backup`. Its Julia environment must develop that copy and contain SCS. It retains dependency caches in a private depot, removes only the package and SCS-extension images before each repetition, then runs `Pkg.precompile(["NCTSSoS", "SCS"])`. Each probe runs in its own fresh Julia process. The script records source hashes, precompile component times, process times, first-call bounds/statuses and first/second-call times. The initial dependency compilation is setup, not part of the reported package/extension total.

The three restoration repetitions used base `aa06227` with the restored workloads. The PR then rebased onto `d97ca7f` (pipeline antiunitary fallback). These timings are not measurements of the post-rebase tree. On that tree, the SCS writer, Aqua, ExplicitImports and export-clash gate passed 154 checks in 196.94 s, within the 2–4 minute estimate. The isolated restoration batch took 587.110 s (9.8 minutes), within the 10–15 minute estimate to rounding. Its 27 first solves all returned `OPTIMAL`. An earlier trial with the extra branches in the original input builder was discarded because it changed the input-compilation measurement; those trial timings are not mixed into this table.

Example, after preparing a disposable source and its environment:

```sh
python3 perf/next/precompile_probe.py \
  --source=.amp/in/orb-backup/precompile/source \
  --environment=.amp/in/orb-backup/precompile/env \
  --output=.amp/in/orb-backup/precompile/before \
  --depot=.amp/in/orb-backup/precompile/depot \
  --variants none P1 P2 P3 E1 E2 E3 all --replicates 3
```

The restored source again has all six blocks; use `--variants all` to measure the selected workload. To reproduce historical alternatives, use their recorded source configuration. For the smallest point, remove only the translated build/lowering and solve from the common-path disposable copy, leaving the plain Pauli and CHSH calls intact. Use separate output directories for each source configuration: resumed measurements must match source hashes. The repository checkout is never rewritten by the harness. `--cases` selects fresh-process probes without adding any workload calls.

The restoration batch uses `--variants all --replicates 3 --cases xxz6 xxz12front hubbard ball CHSH heisenberg8continuous fermionic4number heisenberg8detect xxz6density`.

The old-block batch took about 37 minutes against a 30–60 minute estimate. The first proposed and dense-alternative batches each finished within their 3–6 minute estimates; the two single-sample two-site trials took 3.4 minutes against 2–4 minutes. The smallest useful point took 5.7 minutes against 5–8 minutes. All 175 recorded first solves reached `OPTIMAL` with identical bounds across configurations. Individual elapsed times are in `precompile.json`.
