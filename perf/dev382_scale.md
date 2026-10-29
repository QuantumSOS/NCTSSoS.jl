# DEV-382 Phase 2: a repeated same-budget nonuniform-chain unlock

A relaxation is a tractable optimization whose optimum bounds the original
quantum problem. This campaign computes lower bounds from a fixed Pauli moment
relaxation; it does **not** solve the exact many-body ground state.

**Draft only: documentation regeneration and preview remain blocked.** Both
setup attempts stalled. Generated pages and the regeneration stamp are unchanged;
the stamp check fails. The full package suite passed, but release readiness is
not established. The preserved documentation evidence is linked below.

**The scale and unlock criteria are met against the named generic complex
baseline.** On a nonuniform 22-site chain, three fresh baseline processes hit
the precommitted 180-second external watchdog. All three Track A processes
finished with the same validated bound in 127.76–129.45 seconds, using
5.77–5.78 GiB rather than the baseline's observed 11.69–11.71 GiB before timeout.
The existing translation-invariant specialized-real method rejects this input.
This is not a claim that every possible existing algorithm fails.

Both pinned tracks are integrated. Track A supplies the measured direct-real
reduction. Track B's local second-order cones are present and tested, but this
campaign does not claim a scale gain from B or a composed A+B algorithm.

## The decisive confirmation used an independent watchdog

Machine: **a1.3xlarge orb**, 16 exposed logical CPUs (8 cores, 2 threads/core),
44 GB advertised memory, approximately 43 GiB exposed; Intel Xeon 2.60 GHz,
Linux 6.1.158+, Julia 1.12.7, Mosek 11.2.7, MosekTools 0.15.10. One job at a
time; Mosek 4 threads, Julia/BLAS 1 each. CPU affinity and clock frequency were
not fixed. No test, example, or other solver ran during performance measurements.
See [machine and merge provenance](results/dev382_scale/environment.txt) and the
captured [Project](results/dev382_scale/Project.toml) /
[Manifest](results/dev382_scale/Manifest.toml).

Budget: **180 seconds for the entire fresh process and 24 GiB resident memory**,
identical for both methods. Resident memory means memory held in physical RAM;
the process high-water mark includes compilation, construction and validation.
The watchdog samples every 0.1 seconds and kills the process group at a limit.
The internal solver limit was also 180 seconds, so the external limit expired
first. Observed timeout durations include process reaping, not extra solve time.
The [smoke test](results/dev382_scale/watchdog-smoke.txt) killed a five-second
sleep after 0.2015 seconds under a 0.2-second limit.

Execution order was complex-1, Track-A-1, complex-2, Track-A-2, complex-3,
Track-A-3. [Completion timestamps and frozen source hashes](results/dev382_scale/provenance.txt)
preserve that order. No confirmation protocol changed after any confirmation
outcome. The [protocol](dev382_scale_plan.md) explicitly separates exploratory
results, the one adaptive intermediate size, and held-back confirmation.

| Method, 22-site nonuniform chain | Cold construction s | Solver s | Whole-process s | Peak GiB | Outcome / bound |
|---|---:|---:|---:|---:|---|
| Complex, repeat 1 | 7.522 | interrupted | 180.121 | 11.691 | watchdog timeout; no validated bound |
| Track A, repeat 1 | 11.314 | 98.210 | 129.448 | 5.770 | OPTIMAL; −11.7409135349 |
| Complex, repeat 2 | 8.387 | interrupted | 180.094 | 11.709 | watchdog timeout; no validated bound |
| Track A, repeat 2 | 10.959 | 97.090 | 127.755 | 5.767 | OPTIMAL; −11.7409135349 |
| Complex, repeat 3 | 8.861 | interrupted | 180.173 | 11.706 | watchdog timeout; no validated bound |
| Track A, repeat 3 | 10.604 | 98.696 | 129.206 | 5.785 | OPTIMAL; −11.7409135349 |

The baseline timing is censored: it was interrupted before returning a solver
time or a bound. Do not treat 180 seconds as its completed solve time. Treatment
memory is about half the baseline's **observed pre-timeout** peak, not a measured
ratio of two completed solves. Cold construction is slower with Track A; its
overall benefit comes from the smaller solver representation.

## The input and certificate quality are frozen

The Hamiltonian is the deterministic periodic nearest-neighbor expression in
`perf/dev382_scale.jl:construct`: bond coefficient `(1 + (i % 5)/10)/4` multiplies
`XᵢXⱼ + 0.7YᵢYⱼ + 1.2ZᵢZⱼ`; fields are `0.13Xᵢ + (0.07 + (i % 3)/100)Zᵢ`,
with `j=mod1(i+1,N)`. No random-number generator or seed is used. The actual
Float64 coefficients, all **110 Pauli terms**, and all **265 ordered basis
words** are in [the input record](results/dev382_scale/input-inhomogeneous-22.toml).
The contiguous degree-two basis is identical between modes. The generic model
has 140,715 scalar solver variables; Track A has 35,245.

A positive-semidefinite (PSD) matrix has nonnegative eigenvalues. The audit
recovers PSD Gram matrices in the **original unprojected polynomial algebra**,
accumulates every coefficient of the objective-minus-certificate identity,
and checks the primal/dual objective gap. All three confirmation results have:

- bound: **−11.740913534928156**;
- largest coefficient residual: **7.17e-11**;
- sum of absolute coefficient residuals: **3.3789883e-8**;
- primal/dual gap: **2.5755222e-8**;
- smallest Gram eigenvalue: **+3.7925062e-12**;
- residual-adjusted numerical lower estimate: **−11.740913568718039**.

These satisfy the predeclared tolerances. This is **floating-point certificate
validation**, not rigorous interval or exact-rational certification. No dense
2^N-by-2^N diagonalization is hidden in the large-system timings. Validation is
included in the whole-process watchdog but timed separately from construction
and optimization in the raw records.

One-site translation changes the bond coefficients and fields, so the existing
specialized-real path does not apply. The [input audit](results/dev382_scale/input-audit.log)
actually called `pauli_translation_invariant_moment_relaxation`, with unrelated
sign/reflection options disabled, and received:

> ArgumentError: Translation-invariant Pauli relaxation requires a one-site translation-invariant objective.

The [audit helper](results/dev382_scale/capture_inputs.jl) reuses the measured
generator; it does not substitute another Hamiltonian.

## Speed controls and negative results remain visible

The [complete comparison table](dev382_scale_table.md) links **every one of the
24 process records**, its exact command/limits, per-sample measurements, bounds,
residuals and gaps. Warm timings below are medians of three independently rebuilt
samples after one exact-workload warmup; ranges are in the complete table.

| Matched workload | Complex build / solve s | Track A build / solve s | Bounds, complex / Track A |
|---|---:|---:|---|
| Nonuniform 12, dense | 0.07180 / 15.8343 | 0.05774 / 2.9775 | −6.3821651841 / −6.3821651714 |
| Heisenberg 12, dense | 0.07306 / 11.2651 | 0.04978 / 1.8047 | −5.4377852404 / −5.4377852136 |
| Nonuniform 22, term-sparse | 0.031 / 0.120 | 0.023 / 0.052 | −19.0675000141 / −19.0675000129 |

Nonuniform 12 gives a **5.32× solver-time reduction**, with matching bounds and
nonoverlapping timing ranges (15.793–16.076 versus 2.898–3.261 seconds). Heisenberg
12 gives 6.24× against the generic path. For the larger 22-site case, the decisive
evidence is the repeated whole-process deadline crossing and observed memory
reduction, not an invented completed-baseline speedup.

Important limitations and negative results:

1. **Uniform Heisenberg is not counted as an unlock.** At 24 sites, the existing
   specialized-real method takes 0.0209 solver seconds, versus Track A's 103.164,
   with matching bounds near −10.78547756. It wins decisively and should be used
   when its assumptions hold.
2. **Both dense methods failed on nonuniform 24** under the exploratory 120-second
   solver setting. Actual times were 175.51 and 126.27 seconds because Mosek
   checks that limit coarsely. Those failures remain; the confirmation watchdog
   was declared separately before collecting its data.
3. **Term sparsity changes bound quality here.** At nonuniform 22, the fast sparse
   bound is −19.0675, not the dense −11.7409 result. Sparse speed does not count
   as achieving the dense quality target. At uniform 24 its bound is −18, not
   −10.7855. No stronger sparse iteration or exhaustive alternative-method
   search is claimed.
4. No larger 48/64-site search was needed after locating this boundary. No SOC
   scale speedup or A+B composition is asserted. The earlier Track B report
   [retains its anticommutation and larger-basis counterexamples](pauli_invariant_soc.md).

## Reproduce without changing the recorded artifacts

Install the solver in an external environment, develop this checkout there,
and verify the Mosek license with a trivial solve. The saved environment pins
the measured dependency versions; its local developed-package path must point
to the reproduction checkout. Keep the repository on the load path for JuMP.

```sh
mkdir -p /tmp/dev382-repro
for repeat in 1 2 3; do
  for mode in complex antiunitary; do
    prefix="/tmp/dev382-repro/confirm-inhomogeneous-22-$mode-$repeat"
    DEV382_WALL_LIMIT_S=180 python3 perf/dev382_budget.py "$prefix" -- \
      env DEV382_SOLVER_LIMIT_S=180 \
      JULIA_LOAD_PATH="/tmp/dev382-env:$PWD:@stdlib" \
      OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
      julia --project=/tmp/dev382-env perf/dev382_scale.jl \
      "$mode" inhomogeneous 22 0 "$prefix.toml"
  done
done
```

For the original warm experiments, omit the two `DEV382_*_LIMIT_S` overrides
and use `REPEATS=3`; default limits remain 120 solver seconds / 900 process
seconds / 24 GiB. Modes are `complex`, `antiunitary`, their `-ts` variants, and
`specialized` (uniform Heisenberg only). Regenerate the table with
`python3 perf/dev382_scale_summary.py > perf/dev382_scale_table.md`.

## Delivery and validation record

Combined temporary base: `dev-382-scale-proof`, initially pushed at
[`7119fc1`](https://github.com/QuantumSOS/NCTSSoS.jl/commit/7119fc13a33c8c0289b04644d21b23890f19380a). Evidence branch:
`dev-382-scale-evidence`, targeting that temporary base, not main.

- Track B pinned head: [`4bf5ff2`](https://github.com/QuantumSOS/NCTSSoS.jl/commit/4bf5ff2d9808128c7a39d2f4b1309c6cc46307af)
  ([#380](https://github.com/QuantumSOS/NCTSSoS.jl/pull/380)), retained as a parent
  of merge [`6d9e752`](https://github.com/QuantumSOS/NCTSSoS.jl/commit/6d9e752b9be0e7f2e331ffd818caf87a754d8dd8).
- Track A pinned head: [`a2e1815`](https://github.com/QuantumSOS/NCTSSoS.jl/commit/a2e1815c2a59f8a9879bdcc00e3623b214da9d53)
  ([#381](https://github.com/QuantumSOS/NCTSSoS.jl/pull/381)), retained as a parent
  of merge [`1c5d326`](https://github.com/QuantumSOS/NCTSSoS.jl/commit/1c5d3264a2db942a61678bd8aaad4adca83f6f23).
- The only track-integration conflicts were adjacent source/test include lines;
  both includes were retained. The exact combined diff is in the provenance log.
- All measured intermediate commits remain reachable, including the amended
  exploratory protocol revision. No benchmark implementation or confirmation
  limit changed during confirmation runs. Main was never updated.
- Input audit passed in 28.90 seconds (estimate under one minute).
- [Record audit](results/dev382_scale/record-audit.txt): 24 processes, 59 validated
  solves, 7 retained failures; matched baseline/treatment bounds pass tolerance.
- `julia --project -e 'using Pkg; Pkg.test()'`: **7,637 passed, 2 marked broken**,
  19m02s including setup (estimate 15–25 minutes). Both tracks and doctests ran
  together; the pre-existing marked-broken tests were not changed.
  [Full test log](results/dev382_scale/full-tests.log).
- **Documentation regeneration/preview failed to complete.** First attempt:
  stopped after approximately 25m17s without new setup output; highest observed
  memory 25.3 GiB. [Command, timestamps, memory snapshots and created cache state](results/dev382_scale/docs-attempt1.txt);
  [first raw log](results/dev382_scale/examples-attempt1.log).
- Exactly one retry used the instantiated docs environment and an `8G` heap
  hint, with 600s total, 300s no-output and 12 GiB aggregate-memory limits.
  It stopped at the no-output limit after **301.84s**, peaking at **8.74 GiB**.
  [Retry command, timestamps and watchdog result](results/dev382_scale/docs-retry-budget.txt);
  [retry raw log](results/dev382_scale/examples-retry.log).
- Generated pages and stamp are unchanged. No Julia process remains.
  `make servedocs` was not reached. `julia docs/examples_stamp.jl` exits **1**
  with a stale source hash ([check output](results/dev382_scale/docs-stamp-check.txt),
  2.15s). No stamp was fabricated, source altered to bypass validation, or
  further retry attempted. The 25–40 minute docs estimate did not produce a
  successful regeneration/preview; this remains a delivery blocker.

Grant: `GRANT-dev-382-scale-proof-2026-10-29`, maximum $200. Accounting starts
2026-10-29 13:14:01 Asia/Shanghai, conservatively at $2.13 per running orb hour
plus actual agent/tool fees. At 15:05:12 Asia/Shanghai, elapsed time was 1.853h,
or **$3.95** standard-price equivalent. A two-hour envelope is **$4.26**, well
below the cap. Latest coordinating-thread telemetry reported agent/tool and
Web Extract fees of $0. Linked-subscription inference list-price estimates are
not Amp charges. No additional machine was rented.
