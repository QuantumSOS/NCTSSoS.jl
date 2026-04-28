# Probe 3 — typed expression matrices

**Date:** 2026-04-21  
**Script:** `demos/results/h4_periodic_typed_exprs.jl`  
**Log:** `demos/results/h4_periodic_typed_exprs.log`

## Goal

Replace the local `Matrix{Any}` containers in the JuMP primal build with typed expression matrices, per the suspicion that boxing ~500k PSD scalar entries might be slowing optimizer attachment.

The intended experiment was:

- use the same primal moment + COSMO route,
- keep `max_iter = 10`,
- change only the matrix element container type,
- see whether iteration 1 becomes reachable.

## What actually happened

The run failed immediately during primal-moment model construction.

After the usual symbolic stages:

| Stage | Wall | Status |
| :--- | ---: | :--- |
| build polyopt | 46.86 s | ✅ |
| compute_sparsity (TS/MMD) | 149.46 s | ✅ |
| moment_relax (symbolic) | 9.69 s | ✅ |

…the typed build died on the first matrix fill:

```text
ERROR: LoadError: MethodError: Cannot `convert` an object of type QuadExpr to an object of type AffExpr
```

The failing assignment was into `Matrix{AffExpr}`.

## Follow-up inspection

A quick type inspection on the same code path showed that the current substitution helper is not returning affine expressions here:

- `typeof(re_expr) = QuadExpr`
- `typeof(im_expr) = QuadExpr`

So the naïve replacement

```julia
Matrix{JuMP.GenericAffExpr{Float64,JuMP.VariableRef}}
```

is simply **not valid** on this path as it stands.

That is a useful result, because it means the suggested cleanup was based on a false premise for this concrete model build.

## What this probe established

### 1. Probe 3, as originally phrased, is not runnable without deeper surgery

The local builder is not producing `AffExpr`s. It is producing `QuadExpr`s. So before benchmarking typed expression matrices, someone first has to answer the more basic question:

> Why is `_substitute_complex_poly` yielding quadratic JuMP expressions for what should be a linear moment substitution?

Until that is fixed, `Matrix{AffExpr}` is dead on arrival.

### 2. This further demotes the "typed expressions will save us" theory

Even if we wanted to chase this, it is now clearly a **separate upstream correctness / typing investigation**, not a quick performance tweak.

And after Probes 1 and 2, there is still no evidence that such a fix would materially change the real bottleneck for H4 Nk=2:

- direct COSMO assembly was already fast,
- direct COSMO still stalled in `optimize!`,
- Clarabel also failed to reach visible progress.

So this probe uncovered a real typing oddity, but not a likely escape route for the benchmark.

## Decision: is the next probe still warranted?

**Yes — Probe 4 still matters.**

Probe 4 answers a different, still-useful question: whether smaller in-repo fermionic periodic problems complete cleanly end-to-end under the same general primal-moment + COSMO route.

That helps distinguish:

- "the whole pipeline is fundamentally broken"

from

- "the pipeline is fine on smaller instances, but the H4 Nk=2 problem crosses a brutal setup/factorization threshold."

## Bottom line

**Typed `AffExpr` matrices are not a drop-in replacement here. The current substitution path returns `QuadExpr`, so the experiment failed before solver attachment.**

That is worth knowing, but it is not the main H4 bottleneck story.
