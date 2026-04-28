# H4 periodic Nk=2 V2RDM — COSMO investigation summary

**Date:** 2026-04-21  
**Worktree:** `/Users/exaclior/QuantumSOS/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`  
**Primary artifacts:**

- original diagnostic: `demos/results/h4_periodic_term_sparsity_cosmo_benchmark.md`
- probe 1: `demos/results/h4_periodic_cosmo_direct.{jl,md,log}`
- probe 2: `demos/results/h4_periodic_clarabel_primal.{jl,md,log}`
- probe 3: `demos/results/h4_periodic_typed_exprs.{jl,md,log}`
- probe 4: `demos/results/h4_periodic_smaller_system.{jl,md,log}`

## Executive verdict

**The current generic 32-mode spin-orbital order-2 primal-moment formulation is not a viable path to solving H4 periodic Nk=2 with COSMO or Clarabel.**

That is the answer.

Not because symbolic preprocessing is too slow. Not because JuMP model construction is too slow. Not because stdout was buffered. And not because COSMO's MOI bridge alone is secretly eating 20 minutes.

The problem is that, once assembled, this conic system lands squarely in the part of the solver stack that is brutally expensive: pre-iteration setup / factorization for a very large PSD-heavy SDP.

## What was observed, not guessed

### Baseline H4 diagnostic

From `h4_periodic_term_sparsity_cosmo_benchmark.md`:

| Stage | Wall | Status |
| :--- | ---: | :--- |
| build polyopt | ~44–47 s | ✅ |
| compute_sparsity (TS/MMD) | ~140 s | ✅ |
| moment_relax (symbolic) | ~7.5 s | ✅ |
| primal moment JuMP build | ~14 s | ✅ |
| `optimize!` with COSMO (`decompose=false`) | **>985 s** | ❌ no iteration output |
| `optimize!` with COSMO (`decompose=true`) | **>1302 s** | ❌ no iteration output |

Model size:

- TS PSD blocks: `313`
- largest block: `168 × 168`
- PSD upper-tri slots: `498,845`
- unique moments: `26,817`
- JuMP vars: `171,882`
- JuMP constraints: `16,831`

### Probe 1 — direct COSMO, no JuMP/MOI solver attachment

From `h4_periodic_cosmo_direct.md`:

- direct symbolic-to-sparse build: **5.61 s**
- `COSMO.assemble!`: **7.89 s**
- direct assembled rows/cols: `1,997,801 × 171,882`
- direct `A` nnz: `2,579,805`
- `COSMO.optimize!`: **stalled >643 s** before timeout, with no banner / no iter 1

This killed the most tempting bad theory:

> "Maybe the real bottleneck is just JuMP/MOI copying data into COSMO."

No. It isn't.

### Probe 2 — same JuMP primal model, Clarabel instead of COSMO

From `h4_periodic_clarabel_primal.md`:

- JuMP primal build: **12.64 s**
- Clarabel then produced **no banner, no presolve summary, no iteration output** before the wall budget expired.

So this is not a clean "COSMO is bad, Clarabel is fine" story either.

### Probe 3 — typed expression matrices

From `h4_periodic_typed_exprs.md`:

- the experiment did **not** reach solver attachment,
- the typed `Matrix{AffExpr}` build failed immediately,
- because the current substitution path returns `QuadExpr`, not `AffExpr`, on this route.

That is a real typing issue, but it is not the main H4 performance story.

### Probe 4 — smaller in-repo periodic fermionic stand-ins

From `h4_periodic_smaller_system.md`:

Because the repo only vendors the H4 Nk=2 asset, smaller H4/H2 chemistry assets were unavailable. So the probe used smaller periodic Hubbard stand-ins.

Those smaller cases **did** start iterating and solve end-to-end:

| Case | PSD slots | JuMP vars | COSMO setup | COSMO runtime | Status |
| :--- | ---: | ---: | ---: | ---: | :--- |
| periodic Hubbard N=2 | 83 | 50 | 0.38 s | 1.19 s | `OPTIMAL` |
| periodic Hubbard N=4 canonical | 569 | 1,982 | 2.84 s | 3.84 s | `OPTIMAL` |

This matters because it proves the pipeline is not universally broken. The H4 instance is the outlier because of size/structure, not because the whole primal path never works.

## What is ruled out

These hypotheses are now dead:

- **Not symbolic preprocessing.** That part is stable and finishes in a few minutes.
- **Not JuMP model construction.** That takes ~13–15 seconds on H4.
- **Not stdout buffering.** Flushed logs fixed that.
- **Not COSMO chordal decomposition.** `decompose=false` did not save the H4 run.
- **Not MOI→COSMO transfer as the dominant cost.** Direct `COSMO.assemble!` was fast.
- **Not a COSMO-only oddity.** Clarabel also failed to reach visible progress.

## What remains true

The only story left that fits all evidence is the boring one:

> the current H4 periodic Nk=2 formulation is too large / too ugly in solver setup terms for this solver stack.

More concretely, the bottleneck is in the solver-side pre-iteration regime:

- KKT assembly/factorization,
- cone/setup scaling,
- PSD-heavy linear algebra before iteration 1,
- or some combination of those.

For H4, the direct conic system is already around **2 million constraint rows** and **171k decision variables**. Compared with the smaller N=4 stand-in, the KKT dimension is larger by roughly **two orders of magnitude**. That is exactly where symbolic factorization stops being cute.

## Concrete recommendation for actually solving H4 Nk=2

### Recommendation 1 — stop polishing the current generic formulation

Do **not** spend more time on:

- `Matrix{Any}` micro-optimizations,
- JuMP bridge speculation,
- solver-attribute bikeshedding,
- or hoping another first-order solver will magically eat this raw formulation.

That road is not serious anymore.

### Recommendation 2 — reduce the formulation **before** solving

If you want this problem solved inside NCTSSoS, the next meaningful work is formulation reduction, not solver tweaking.

The high-leverage targets are:

1. **Momentum-aware / k-blocked formulation**  
   The current generic spin-orbital lift throws the whole Nk=2 chemistry problem into a broad fermionic moment basis and only later tries to recover structure via term sparsity. That is backwards for this case. The momentum conservation is physical structure; it should be baked into the basis from the start.

2. **Spin / particle-number sectoring structurally, not just by moment equalities**  
   Right now, particle number is enforced through moment-equality constraints after the basis is already enormous. Better is to build the relevant sector into the problem representation so the moment system never grows that large in the first place.

3. **Chemistry-specific 2-RDM formulation rather than generic NC-polynomial lift**  
   If the real goal is "get the H4 periodic V2RDM number," then a generic 32-mode NC moment hierarchy is the wrong hammer. Use a representation closer to the actual V2RDM object, with the known symmetries and selection rules encoded structurally.

Put bluntly: **earn your structure early**. Don't ask the solver to discover it after you've already blown the model up.

### Recommendation 3 — once the formulation is smaller, retest with a stronger solver path

After a real formulation reduction, then it is worth testing:

- direct COSMO again,
- Clarabel again,
- and especially Mosek if available.

But doing that **before** shrinking the model is just paying electricity to learn the same lesson again.

### Recommendation 4 — treat the `QuadExpr`/typed-matrix issue as a separate cleanup

Probe 3 found a real wart:

- the comment in `src/optimization/moment.jl` suggesting `Matrix{GenericAffExpr}` is misleading for this path,
- `_substitute_complex_poly` is yielding `QuadExpr` on the current generic-model route.

That deserves its own follow-up ticket because it is ugly and probably not what the code author intended.

But it is **not** the blocker for H4 Nk=2.

### Recommendation 5 — vendor a real smaller chemistry asset if future scaling studies matter

If you want a clean chemistry-specific shrink-and-extrapolate story later, vendor one of:

- H4 Nk=1,
- or H2 Nk=2,

into `test/data/assets/` and give it the same benchmark harness. Right now the repo simply does not contain that intermediate chemistry point.

## Should the dual route be revisited?

**Not as the next move on the current raw formulation.**

The probes do not support the fantasy that the same giant model, merely dualized, will suddenly become easy through COSMO/Clarabel. Maybe the dual behaves differently under a strong commercial solver or after symmetry reduction. Fine. But on today's evidence, revisiting the dual *before* shrinking the formulation is low-odds work.

## Recommended next implementation step

If this were my branch, I would do exactly this next:

1. open a focused issue/spec for a **k-blocked / symmetry-reduced H4 periodic formulation**,
2. keep the new probe scripts and reports as evidence,
3. stop spending time on generic primal-moment solver tuning for this case,
4. only come back to solver benchmarking after the reduced formulation exists.

If the immediate goal is simply to obtain the H4 Nk=2 V2RDM energy, the honest answer is even simpler:

> use a more specialized chemistry/V2RDM code path or the original research implementation, because the current generic NCTSSoS formulation is not the practical route to that number.

## Final bottom line

**H4 periodic Nk=2 is not failing because of preprocessing overhead or a dumb bridge bug. It is failing because the assembled conic problem is too large in exactly the solver-setup dimensions that matter.**

So the path forward is not another micro-optimization. It is a different formulation.
