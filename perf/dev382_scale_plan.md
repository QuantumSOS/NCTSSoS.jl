# DEV-382 Phase 2 predeclared scale protocol

Recorded before the first scale solve on 2026-10-29. Both original PR heads
are merged on `dev-382-scale-proof`, never main. Grant:
`GRANT-dev-382-scale-proof-2026-10-29`.

## Question and acceptance

Does direct-real antiunitary lowering reduce solve time or memory at the same
Pauli relaxation, and allow a solve that the complex baseline cannot finish
within the same fixed budget? SOC means second-order cone; both tracks are
integrated, but this campaign does not assert that SOC and antiunitary APIs
compose. Track B's earlier negative comparisons remain in its original report.

Families: periodic Heisenberg chain and a fixed inhomogeneous real nearest-
neighbor chain (coefficients in the harness). Both use a contiguous degree-two
Pauli basis. The latter removes translation invariance, so the specialized
translation-invariant path is not applicable. No synthetic random search.
Initial sizes: 12, 24, 48; extend to 64 only if needed and economically feasible.
Stop increasing a family once both paths fail or a useful resource boundary
is established. Keep every attempted configuration, including failures.

Same machine: a1.3xlarge orb, 16 exposed CPUs, 44 GB advertised RAM, Linux x86-64.
One job at a time; Julia and BLAS one thread, Mosek four threads. No concurrent
tests during measurements. Tolerances 1e-8. All modes run in fresh processes.

**Fixed budget, before results:** 120 seconds of Mosek optimizer time per solve;
900 seconds total process time (startup, one cold run, three warmed samples,
validation); 24 GiB peak resident process memory. The external runner enforces
wall time and memory; Mosek enforces solve time. An unlock requires an OPTIMAL,
validated treatment and baseline timeout, memory limit, or solver failure with
these identical limits, not just a slower baseline. A repeat failure is recorded,
not silently dropped. Cold and warm times are separate. Ratios use three warmed
samples, median and range; a single success is only a feasibility observation.

Validation reconstructs the original polynomial from positive-semidefinite
(PSD: all eigenvalues nonnegative) Gram matrices, checks every coefficient,
the smallest eigenvalue, and the primal/dual objective gap. Each Pauli word
has operator norm one, so coefficient residual l1 gives a bound-error estimate.
Acceptance: maximum coefficient residual <=1e-6, residual l1 plus diagonal PSD
repair <=1e-5*N, eigenvalue >=-1e-6, gap <=1e-5*N. This is floating-point
certificate validation, **not** exact rational or interval certification.
Matching baseline/treatment bounds must agree within 2e-5*N. Sparse comparisons
are matched separately; a weaker sparse bound is not the dense result.

The specialized-real Heisenberg method is a required comparator, not omitted
when it wins. Dense baseline/antiunitary pairs are primary; term-sparse pairs
at the largest successful size test combination with existing methods.

## Reproduction

Create an external Julia environment with this checkout developed into it and
MosekTools installed; verify the Mosek license with a trivial solve first.

```
python3 perf/dev382_budget.py OUTPUT_PREFIX -- \
  env JULIA_LOAD_PATH='/tmp/dev382-env:/home/user/workspace/repo:@stdlib' \
  OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 \
  julia --project=/tmp/dev382-env perf/dev382_scale.jl \
  MODE FAMILY N 3 OUTPUT_PREFIX.toml
```

The runner saves command, wall time, peak resident memory, and termination reason.
The Julia record saves per-sample construction, solve, validation, certificate
quality, cone sizes, version and commit metadata. External process memory includes
compilation and validation; it is not isolated solver memory.

Budget accounting: $2.13 per running orb hour plus reported actual tool/agent
fees, upper limit $200. Initial time window 3–5 hours ($6.39–$10.65 orb equivalent).
Source Puck reported linked-subscription inference and no Amp agent charges at
13:19 Asia/Shanghai. Recheck before extended work; stop well before risking cap.

## Independent confirmation protocol (13:33 Asia/Shanghai)

The exploratory 24-site complex solve returned TIME_LIMIT after 170.47 seconds
despite the 120-second solver setting; the first treatment solve finished in
105.61 seconds. These are exploratory observations, not the confirmation data.
At Puck's request, freeze a second, externally enforced **180-second whole-
process** limit (120 seconds plus 60 seconds for startup/build/validation), with
the same 24 GiB memory limit, before collecting confirmation outcomes. Use one
cold solve per fresh process (`REPEATS=0`) and three independent repetitions per
mode, alternating complex then antiunitary. All processes keep the same 120-
second solver setting, threads, tolerances, and inputs. The Python watchdog
checks every 0.1 seconds and kills the isolated benchmark process group on limit;
reported elapsed includes process reaping. Set `DEV382_WALL_LIMIT_S=180` on the
runner command. Save under distinct `confirm-*` names; never replace exploratory
records. Require all three treatment runs to validate, and report every baseline
outcome. This confirms a resource boundary for the named baseline, not that all
existing methods fail: the specialized-real and sparse controls remain mandatory.
