# Precompile scope and first-call cost

**Selected by the owner on 2026-11-15: retain the restored workload and add continuous symmetry, fermionic conditions and detection.** His latest decision was: "I don't care about precompile"; first-call time decides. The original six blocks remain unchanged. The package also builds and lowers a four-site Heisenberg ring with translation/reflection/spin rotations, a four-mode hopping ring with P/Q/G, Sector(N,2) and NumberRotation, and an eight-site Heisenberg ring with detected symmetry. SCS solves these representatives. The 45 s and 28 s workloads below remain measured alternatives, not the selected scope.

The Mosek extension **only constructs tasks**, in both SOS and moment forms, for the same representatives. It calls `precompile` on the typed native solve entry point to compile solution extraction without executing it, then deletes the tasks. There is no optimizer call, license check or license-server dependency in the workload. The earlier solve-based variant was rejected for that dependency. These workloads do not change solver tolerances, accepted statuses or model semantics.

The sizes follow measurements, not a size-independent compilation assumption. Two fermionic modes left about 0.1 s of size-specialized `hcat` compilation on four modes. Four-site detection left tuple specializations such as `NTuple{24,Tuple{Int,NTuple{7,Int}}}` to compile on eight sites; a trace showed a 1.673 s first detection versus 0.068 s warm. The detector still uses the legacy bridge; this change precompiles the eight-site case without changing that algorithm.

## Expanded workload: measured first calls

Measurements use base `18a1476`. Times below are seconds, medians; brackets give minimum–maximum. For the selected construction-only Mosek workload, three isolated precompiles took **132.162 s [122.378,138.715]** in total; package-plus-SCS took **108.284 s [100.331,115.848]**. The package and SCS component medians are 98.199 and 10.085 s; the native extension takes about 23 s. Component medians need not sum to a total's median.

Each construction-only precompile was followed by three fresh-process native probes. That batch took **470.17 s** (146.40 s first trial, then 323.77 s for two more repetitions). The SCS rows come from the preceding three-precompile/33-process batch, which took **632.37 s**, within its 10–15 minute estimate. Its package/SCS source hashes are unchanged: SCS probes never load the Mosek extension. Numerical measurement jobs ran serially. The small late native/quality diagnostics are not timing samples.

| Backend / full-size probe | Input s | Detection s | First solve s [min,max] | Second solve s | Bound |
|---|---:|---:|---:|---:|---:|
| SCS / Heisenberg8 continuous, order 2 | 0.665 | — | 5.397 [4.184,6.179] | 5.1927 | −14.60683422165745 |
| SCS / fermionic4 P/Q/G + number sector | 0.789 | — | 0.042 [0.032,0.048] | 0.0112 | −1.9999999739663321 |
| SCS / Heisenberg8 detected, order 1 | 2.924 | 0.087 | 0.024 [0.022,0.025] | 0.0145 | −24.00000000077552 |
| Mosek / Heisenberg8 continuous, order 2 | 0.805 | — | 0.309 [0.233,0.363] | 0.2510 | −14.606832405653952 |
| Mosek / fermionic4 P/Q/G + number sector | 0.957 | — | 0.165 [0.158,0.187] | 0.0181 | −2.0000000053016898 |
| Mosek / Heisenberg8 detected, order 1 | 3.168 | 0.088 | 0.104 [0.101,0.121] | 0.0163 | −23.999999996467714 |

All 33 selected first solves returned `OPTIMAL`. Bounds were identical across repetitions within each backend. Every construction-only Mosek first solve is below 1 s, meeting the architect's threshold for removing licensed optimization from precompile. Mosek uses exactly `MosekBackend()`, including default tolerances and automatic solver thread selection; SCS uses `verbose=0, eps_abs=eps_rel=1e-7`. Julia, BLAS and OpenMP use one thread. **This is not a cross-solver speed comparison:** settings and achieved bounds differ. Input compilation and detection are excluded from solve times and reported separately; loading is in the raw data.

The continuous and detected SCS solves are within twice their warm calls in all three processes. **The fermionic case is not:** first/warm ratios are 3.75–4.32, despite falling from the historical 17.442 s first call to 0.042 s. A fresh-process `--trace-compile-timing` run identifies residual Base specializations: `copy(::Memory{UInt64})` (6.5 ms), reverse `issorted(::Vector{Int})` (5.4 ms), three-matrix complex `hcat` (4.4 ms), complex `iszero` (3.1 ms), and integer `reduce_first` (2.5 ms). The trace's solve was 0.032 s versus 0.0081 s warm. No pipeline algorithm was changed to remove these milliseconds. The Mosek fermionic and detected calls also retain a cold/warm gap; they are not claimed to meet the two-times target.

Separate traces of the solve-based Mosek trial found the same five Base specializations in the fermionic call. The detected call still compiled `Base.hashindex(::Tuple{Int},Int)`, `Iterators.product(::UnitRange{Int})` and `divrem(::Int,::Int)`. These traces identify remaining compilation in that trial, not the cause of every millisecond of cold/warm difference in the selected construction-only variant; native solver initialization is not isolated. They were diagnostic runs alongside correctness testing and are not included in the timing table. The architect accepted the 42 ms versus 11 ms SCS residual: the target was to remove seconds of compilation, not these milliseconds.

Mosek needed native compilation coverage, but not an optimizer run:

| Mosek workload | Continuous first s | Fermionic first s | Detected first s | Total precompile s |
|---|---:|---:|---:|---:|
| None (one control per case) | 7.053 | 5.764 | 5.025 | 109.450 |
| Real solves (three repetitions; rejected) | 0.252 | 0.073 | 0.036 | 118.224 |
| Both task forms + compile-only solve entry (three repetitions; selected) | 0.309 | 0.165 | 0.104 | 132.162 |

The no-native-workload control's warm solves were 0.1157, 0.0129 and 0.0123 s. It is a single control, not a timing distribution. All three variants have exactly the same bounds. The selected variant trades a small cold-call increase for precompile that cannot wait on a license server. It does not need exception handling or a non-waiting license parameter because it never calls the optimizer. The higher measured total is accepted; it is not evidence that construction is inherently slower than solving in precompile, given the shared orb and differing compiled paths.

The original five probes retain their bounds and show no new multi-second cold call:

| SCS probe | Input s | First solve s [min,max] | Second solve s |
|---|---:|---:|---:|
| xxz6 | 0.859 | 0.787 [0.743,0.851] | 0.6369 |
| xxz12front | 1.032 | 0.036 [0.030,0.038] | 0.0289 |
| hubbard | 0.819 | 0.114 [0.111,0.145] | 0.1087 |
| ball | 0.927 | 0.003 [0.002,0.003] | 0.0017 |
| CHSH | 0.896 | 0.001 [0.001,0.001] | 0.0007 |

Each original first/warm ratio is below 1.8. Their five bounds are bit-identical to the restoration values below. Absolute timings still vary in the shared orb: the 0.145 s Hubbard sample also took 0.1447 s warm, rather than showing extra first-call compilation. The four-site detection trial and two-mode fermionic trial remain in `precompile.json`, labeled as trials rather than mixed into this table.

Reproduce with the disposable source/environment setup below and `--variants all --replicates 3 --backend both --cases xxz6 xxz12front hubbard ball CHSH heisenberg8continuous fermionic4number heisenberg8detect`. The native-only comparison used `--backend mosek --cases heisenberg8continuous fermionic4number heisenberg8detect`. `first_call.jl` is unchanged. To inspect remaining compilation, run the extra probe in a new process with Julia's `--trace-compile=<file> --trace-compile-timing` flags; the optional second script argument `mosek` selects the default native backend.

Validation on the final import spelling: with an empty home and all `MOSEK*`/`MSK_*` environment variables removed, a freshly rebuilt Mosek extension precompiled and then loaded in a separate process. Both exited 0 with **0 stdout bytes and 0 stderr bytes**, in 29.45 s combined. The check runs `Pkg.precompile(["NCTSSoS","Mosek"]; strict=true, io=devnull)` so normal Pkg progress is silent, while capturing both output streams; it does not suppress workload output. Three final first-call smoke probes took 0.231, 0.177 and 0.122 s with identical bounds. The only change after the timed batch was replacing qualified private helper access with explicit imports of the same functions. The final native checks passed 147 assertions, Aqua passed 11, ExplicitImports was clean, and `julia docs/examples_stamp.jl` exited 0. Full Next and native SCS checks also passed; the COSMO suite and native/quality gate took 718.17 s. As agreed for the redesign, no docs page uses Next and only the stamp check is required, not example regeneration.

## Restoration and historical choices on the same five probes

Times are seconds, medians of three isolated processes per configuration. Precompile totals are the median of the per-repetition package-plus-extension sums; separately computed component medians need not sum to that total. Ranges are minimum–maximum. First-call times measure `solve` after input construction, not loading or input construction; those additional times and second calls are retained in [raw measurements](precompile.json).

| Workload | Package s | SCS extension s | Total s [min,max] | xxz6 | xxz12front | hubbard | ball | CHSH |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Old, all six blocks | 60.523 | 7.494 | 68.016 [65.490,74.946] | 0.954 | 0.022 | 0.103 | 0.002 | 0.001 |
| Restored | 67.787 | 9.460 | 76.976 [75.334,78.745] | 0.972 | 0.028 | 0.089 | 0.002 | 0.002 |
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

The architect accepted the restoration's timing deviations on 2026-11-15 after checking that the first and second calls were comparable and the bounds agreed. The intervening pipeline change and the 4–20 ms deviations were not grounds to remove workload paths. The later expansion below follows the owner's new decision, not a change to the old measurements.

## Baseline: four paths outside the restored workload

These historical measurements preceded the owner's decision to add the first three paths. Each row is the median of three fresh processes. Input includes model-builder compilation; loading is excluded. Detection is timed separately and subtracted from the input time. The second solve is a warm-call reference, not the first-call result. All twelve first solves returned `OPTIMAL`, with identical bounds across repetitions.

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

`precompile_probe.py` changes workload blocks only in a disposable source copy under `.amp/in/orb-backup`. Its Julia environment must develop that copy and contain SCS and Mosek. It retains dependency caches in a private depot, removes the package and both extension images before each repetition, then runs `Pkg.precompile(["NCTSSoS", "SCS", "Mosek"])`. Run `Pkg.resolve()` after refreshing a disposable source whose extension declarations changed. Each probe runs in its own fresh Julia process. The script records source hashes, precompile component times, process times, first-call bounds/statuses and first/second-call times. The initial dependency compilation is setup, not part of the reported package/extension total. Historical rows used the then-current package/SCS-only harness.

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

Use `--variants all` to measure the selected workload. To reproduce historical alternatives, use their recorded source configuration. For the smallest point, remove only the translated build/lowering and solve from the common-path disposable copy, leaving the plain Pauli and CHSH calls intact. Use separate output directories for each source configuration: resumed measurements must match source hashes. The repository checkout is never rewritten by the harness. `--cases` selects fresh-process probes without adding any workload calls. `--backend both` adds default-Mosek probes for the selected extra cases; it does not change the original five SCS probes.

The restoration batch uses `--variants all --replicates 3 --cases xxz6 xxz12front hubbard ball CHSH heisenberg8continuous fermionic4number heisenberg8detect xxz6density`.

The old-block batch took about 37 minutes against a 30–60 minute estimate. The first proposed and dense-alternative batches each finished within their 3–6 minute estimates; the two single-sample two-site trials took 3.4 minutes against 2–4 minutes. The smallest useful point took 5.7 minutes against 5–8 minutes. All 175 recorded first solves reached `OPTIMAL` with identical bounds across configurations. Individual elapsed times are in `precompile.json`.
