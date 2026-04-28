# Probe 4 — shrink-and-extrapolate with smaller in-repo periodic fermionic systems

**Date:** 2026-04-21  
**Script:** `demos/results/h4_periodic_smaller_system.jl`  
**Log:** `demos/results/h4_periodic_smaller_system.log`

## Goal

Check whether the same general **fermionic primal-moment + COSMO** route completes end-to-end on smaller systems, and get a rough sense of how violently the setup cost grows with problem size.

## Asset reality check

Per `test/H4PeriodicAssets.jl`, the repo only vendors **one** H4 periodic asset: **Nk = 2**.

So the requested ideal proxies — H4 Nk=1 or H2 Nk=2 — are simply not present in this worktree. Rather than inventing fake chemistry data, this probe used the smallest honest **in-repo periodic fermionic stand-ins**:

1. periodic 2-site Hubbard ring, order 2,
2. periodic 4-site Hubbard ring at canonical half-filling, order 2.

These are not scientific substitutes for the H4 V2RDM instance. They are pipeline sanity checks.

## Important caveat about objective values

The logged objective mismatches versus the test oracles are **expected here**.

Why? Because the test oracles in `test/data/expectations/hubbard.toml` come from the package's standard `cs_nctssos(..., SolverConfig(order=2))` route, while this probe deliberately used the **same TS/MMD primal-moment construction style** as the H4 benchmark. That is a different relaxation, so the bound can be weaker. The question here is **completion and scaling**, not reproducing those dense-reference numbers.

## Observed results

### Case A — periodic Hubbard N=2 (grand-canonical)

| Quantity | Value |
| :--- | ---: |
| TS PSD blocks | 23 |
| largest block | 5 × 5 |
| PSD upper-tri slots | 83 |
| unique moments | 25 |
| JuMP vars / constraints | 50 / 25 |
| COSMO setup time | 0.38 s |
| COSMO runtime | 1.19 s |
| first iteration reached? | **Yes** |
| status | `OPTIMAL` |

### Case B — periodic Hubbard N=4 (canonical half-filling)

| Quantity | Value |
| :--- | ---: |
| TS PSD blocks | 87 |
| largest block | 9 × 9 |
| PSD upper-tri slots | 569 |
| unique moments | 239 |
| JuMP vars / constraints | 1,982 / 1,149 |
| COSMO setup time | 2.84 s |
| COSMO runtime | 3.84 s |
| first iteration reached? | **Yes** |
| status | `OPTIMAL` |

In both smaller cases, COSMO printed its banner immediately, reported setup time, reached iteration 1, and finished normally.

## Compare that with H4 periodic Nk=2

From the original H4 diagnostic plus Probe 1:

| Quantity | Smaller N=4 Hubbard | H4 periodic Nk=2 |
| :--- | ---: | ---: |
| PSD upper-tri slots | 569 | 498,845 |
| unique moments | 239 | 26,817 |
| JuMP vars | 1,982 | 171,882 |
| JuMP constraints | 1,149 | 16,831 |
| direct assembled rows `m` | 22,142 | 1,997,801 |
| direct `A` nnz | 77,657 | 2,579,805 |
| solver setup behavior | starts immediately | no banner / no iter 1 after >10 min |

That PSD-slot ratio alone is about **877×**. The direct H4 conic system also blows the KKT dimension up from roughly **27k** to well over **2 million** rows+cols. That is not a "small constant factor" problem. It is an entirely different regime.

## What this probe established

### 1. The pipeline is not fundamentally broken

This matters.

The same broad route — symbolic build, term sparsity, primal moment JuMP model, COSMO solve — **does work** on smaller periodic fermionic problems. So we are not looking at a universal bug where the primal path never reaches iteration 1.

### 2. H4 Nk=2 is crossing a size/setup cliff

Also matters.

The smaller stand-ins start iterating in seconds. H4 Nk=2 never even shows solver setup progress within the budget. That strongly supports the boring explanation — which is usually the right one:

> the H4 instance is simply too large and too nasty in solver setup / factorization terms for this current formulation.

Not mysterious. Just expensive in exactly the wrong place.

### 3. Probe 4 is enough; no further "maybe it just needs more patience" fantasy

The smaller cases prove the machinery can run. The H4 case proves this particular formulation is beyond the practical range of that machinery, at least with COSMO/Clarabel through the current path.

## Decision: are more probes still warranted?

**No more diagnostic probes are really needed to answer the original question.**

At this point the evidence is already overdetermined:

- symbolic stages are fine,
- JuMP build is fine,
- direct COSMO assembly is fine,
- direct COSMO still stalls,
- Clarabel also stalls,
- smaller periodic fermionic systems solve end-to-end.

That's enough to recommend a concrete path forward without pretending there is one more cute local tweak that will save H4 Nk=2.

## Bottom line

**Available smaller periodic fermionic systems solve cleanly and quickly; H4 periodic Nk=2 does not.**

So the present H4 formulation is not hitting a minor implementation wart. It is hitting a genuine solver-scale wall.
