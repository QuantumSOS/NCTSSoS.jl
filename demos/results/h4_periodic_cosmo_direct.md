# Probe 1 — bypass JuMP and build COSMO directly

**Date:** 2026-04-21  
**Script:** `demos/results/h4_periodic_cosmo_direct.jl`  
**Log:** `demos/results/h4_periodic_cosmo_direct.log`

## Goal

Remove JuMP and the MOI-to-COSMO copy path from the equation.

This probe builds the primal moment SDP directly as:

- `x = [y_re; y_im]` with `171,882` scalar decision variables,
- one normalization `ZeroSet` constraint,
- `313` PSD constraints encoded directly as `COSMO.PsdConeTriangle`,
- zero constraints encoded directly from the symbolic moment blocks,
- then times:
  1. symbolic-to-sparse direct data build,
  2. `COSMO.assemble!`,
  3. `COSMO.optimize!`.

If `assemble!` is slow, the pathology is data transfer. If `assemble!` is fast and `optimize!` still stalls before iteration 1, the pathology is inside COSMO's own setup/solve path.

## Observed timings

The invariant symbolic stages stayed in line with the prior diagnostic:

| Stage | Wall | Status |
| :--- | ---: | :--- |
| build polyopt | 46.31 s | ✅ |
| compute_sparsity (TS/MMD) | 135.12 s | ✅ |
| moment_relax (symbolic) | 7.51 s | ✅ |

Direct-COSMO-specific stages:

| Stage | Wall | Alloc | Status |
| :--- | ---: | ---: | :--- |
| direct data build (`MomentProblem` → sparse `A,q,constraints`) | 5.61 s | 22.94 GiB | ✅ |
| `COSMO.assemble!` | 7.89 s | 54.66 GiB | ✅ |
| `COSMO.optimize!` | **> 643 s** before tool timeout | — | ❌ no iteration output |

The `bash` wrapper had a 900 s wall limit. After subtracting the completed upstream stages, `COSMO.optimize!` had already spent roughly **643 s (~10.7 min)** with no banner, no timing table, and no iteration 1 output before the run was killed.

## Model footprint seen by direct COSMO

| Quantity | Value |
| :--- | ---: |
| direct variables | 171,882 |
| symbolic constraint blocks walked | 8,571 |
| COSMO constraint objects created before merge | 8,572 |
| final assembled cone count | 314 |
| assembled rows `m` | 1,997,801 |
| assembled columns `n` | 171,882 |
| assembled `A` nnz | 2,579,805 |
| assembled `b` nnz | 1 |

A notable detail: COSMO merged the many zero-set pieces down to a final cone count of `314`, i.e. one big zero cone plus the `313` PSD cones. That is sane. Nothing in the direct assembly path looked pathological.

## What this probe actually established

### 1. The JuMP/MOI transfer is **not** the main bottleneck

That hypothesis took a serious hit.

The whole direct sparse build plus `COSMO.assemble!` completed in about **13.5 seconds**. The earlier JuMP run spent **15–20 minutes** after entering `optimize!`. That giant wall is plainly **not** explained by "copying data into COSMO".

### 2. The stall survives after bypassing JuMP entirely

That matters more.

This run never touched JuMP's optimizer-attachment path, never relied on MOI bridges to construct COSMO constraints, and still stalled in `COSMO.optimize!` before iteration 1. So the core pathology sits **inside COSMO's own setup / solve path** on this assembled SDP.

The remaining suspects are now much narrower:

- COSMO's internal presolve / scaling,
- chordal-analysis-adjacent setup internal to `optimize!` even with `decompose=false`,
- initial KKT matrix formation / symbolic factorization / numeric factorization,
- or some other COSMO-internal pre-iteration step.

What it is **not** anymore:

- not `moment_relax`,
- not JuMP model construction,
- not MOI-to-COSMO affine-copy overhead as the dominant cost.

### 3. Probe 3 just got demoted

Typed `Matrix{AffExpr}` instead of `Matrix{Any}` might still shave some JuMP-side overhead. Fine. But after this direct-COSMO run, that is polishing the wrong piece first.

Even if typed expressions remove a minute or two from JuMP/MOI handling, they do not explain a direct-COSMO run that still hangs for 10+ minutes before iteration 1.

## Decision: is the next probe still warranted?

**Yes — Probe 2 is now more warranted than before.**

The right next question is no longer "is MOI transfer killing us?" The right next question is:

> Is this a COSMO-specific solver/setup pathology, or would another solver start immediately on the same JuMP primal model?

So the next probe should stay exactly as planned:

- **Probe 2:** keep the JuMP primal model, swap `COSMO.Optimizer` for `Clarabel.Optimizer`, cap at `max_iter = 10`, and see whether Clarabel reaches iteration 1 in reasonable time.

## Bottom line

**Bypassing JuMP did not fix the stall. `COSMO.assemble!` is fast (~8 s). The hang is inside `COSMO.optimize!` itself, before iteration 1.**

That is the kind of result you actually want from a diagnostic: one big hypothesis is dead now.
