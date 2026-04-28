# Probe 2 — keep the JuMP primal model, swap COSMO for Clarabel

**Date:** 2026-04-21  
**Script:** `demos/results/h4_periodic_clarabel_primal.jl`  
**Log:** `demos/results/h4_periodic_clarabel_primal.log`

## Goal

Hold the JuMP primal model fixed and change only the solver backend.

If Clarabel quickly reaches banner / presolve / iteration 1 while COSMO does not, then the pathology is likely COSMO-specific. If Clarabel also stalls before first solver output, then this problem structure is poisonous more generally, and JuMP-side typed-expression cleanup is unlikely to be the decisive fix.

## Observed timings

Invariant symbolic stages remained normal:

| Stage | Wall | Status |
| :--- | ---: | :--- |
| build polyopt | 49.39 s | ✅ |
| compute_sparsity (TS/MMD) | 139.13 s | ✅ |
| moment_relax (symbolic) | 7.60 s | ✅ |
| primal moment JuMP build | 12.64 s | ✅ |

The JuMP model size was the same as before:

- variables: `171,882`
- constraints: `16,831`

After that, the run entered:

```text
[stage] optimize! (Clarabel) : running SDP solver ...
---------------- begin Clarabel output ----------------
```

…and then produced **nothing else** before the 15-minute wall budget expired and the run was terminated.

## What this probe established

### 1. Clarabel did **not** provide an easy escape hatch

That was the whole point of this probe, and the answer is no.

Clarabel never printed its banner, settings, presolve summary, factorization summary, or iteration table. So this is **not** a story where "COSMO is uniquely broken but Clarabel sails through the same JuMP model." It doesn't.

### 2. The issue is not just "COSMO's MOI bridge"

That hypothesis was already badly wounded by Probe 1. This probe buried it.

- **Probe 1:** direct COSMO assembly was fast, but direct `optimize!` still stalled.
- **Probe 2:** JuMP + Clarabel also stalled before first visible solver progress.

So the bottleneck is not some cute little one-off in the COSMO bridge. The large H4 Nk=2 primal moment SDP is expensive in the exact part solvers care about most: the pre-iteration linear-algebra setup around these PSD blocks and equality structure.

### 3. Typed `Matrix{AffExpr}` is now a side quest, not the main road

Probe 3 was originally sensible because `Matrix{Any}` in the JuMP build smelled bad. It still smells bad. But bad smell is not the same as root cause.

With both:

- direct COSMO stalling after fast assembly, and
- JuMP + Clarabel stalling before iteration 1,

there is no serious evidence that replacing `Matrix{Any}` will turn this from "doesn't start" into "solves". At best it would shave some JuMP/MOI overhead from a workflow that still looks dominated by solver setup.

## Decision: is the next probe still warranted?

**Probe 3 is no longer strongly warranted.**

If the goal were upstream code hygiene, sure, typed expressions are still a good cleanup. If the goal is to get this H4 Nk=2 V2RDM problem actually moving, Probe 3 has been demoted.

The more useful remaining probe is now:

- **Probe 4 — shrink-and-extrapolate.**

That probe still matters because it answers a different question:

> Is the pipeline fundamentally fine on smaller fermionic SDPs, with this H4 case simply beyond the practical setup/factorization limit?

That is actionable. A typed-expression micro-optimization, right now, mostly is not.

## Bottom line

**Clarabel also failed to reach visible solver progress within the time budget.**

So the current evidence says this H4 Nk=2 order-2 primal moment SDP is bottlenecked by heavy solver setup / factorization cost for the assembled conic system itself, not by a COSMO-only transfer bug and not by JuMP model construction.
