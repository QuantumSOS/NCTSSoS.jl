# H4 periodic Nk=2 — PQG basis + COSMO benchmark

**Date:** 2026-04-22  
**Demo:** `demos/h4_periodic_pqg_cosmo_benchmark.jl`  
**Routes compared:**

- `--route=fat --direct=true` — single 2017×2017 PQG block
- `--route=paper --direct=true` — manual paper-style `(ΔN, K)` blocks `[513, 512, 256, 256, 240, 240]`

## What was changed relative to the earlier TS benchmark

This benchmark does four things the TS route did not:

1. replaces the dense order-2 basis with the explicit PQG / V2RDM basis,
2. bypasses generic TS fragmentation,
3. compares the **single-fat-block** and **manual paper-block** routes directly,
4. adds the paper-style scalar trace targets through `moment_eq_constraints`:
   - `N_up = 4`
   - `N_dn = 4`
   - `Tr(²D) = 28`
   - `Tr(²Q) = 276`
   - `Tr(²G) = 200`

Important correction: the CAR-derived PQG linear maps are **implicit** here because
`²D`, `²Q`, and `²G` are slices of one shared canonical moment map. They are not
modeled as separate free matrices, so they do not need to be re-added as explicit
linear constraints.

## Static block comparison

| Route | PSD blocks | largest block | Σ dim² | upper-tri slots |
| :--- | ---: | ---: | ---: | ---: |
| single PQG block | 1 | 2017 | 4,068,289 | 2,035,153 |
| paper-style `(ΔN, K)` blocks | 6 | 513 | 771,585 | 386,801 |

So the manual paper blocking cuts the PSD footprint by about:

- **3.93×** on max block size,
- **5.27×** on `Σ dim²`,
- **5.26×** on upper-tri PSD slots.

That is real structural improvement, not paperwork.

## Timed runs

### 1) Single-fat-block route

Command:

```bash
julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl --route=fat --direct=true
```

Observed stages before external timeout:

| Stage | Value |
| :--- | ---: |
| shared build | 208.63 s |
| PQG basis size | 2017 |
| route `moment_relax` | 28.19 s |
| route unique moments | 637,393 |
| route `total_basis` | 637,393 |
| direct variables | 1,274,786 |
| direct rows | 8,154,607 |
| direct A nnz | 8,401,311 |
| `COSMO.assemble!` | 40.79 s |
| assembled cone count | 2 |

Outcome:

- `moment_relax` completed.
- direct sparse conic assembly completed.
- `COSMO.optimize!` started.
- **No COSMO iteration banner / no iter 1 before timeout.**

### 2) Paper-style `(ΔN, K)` route

Command:

```bash
julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl --route=paper --direct=true
```

Observed stages before external timeout:

| Stage | Value |
| :--- | ---: |
| shared build | 207.98 s |
| PQG basis size | 2017 |
| route `moment_relax` | 6.03 s |
| route unique moments | 123,649 |
| route `total_basis` | 170,273 |
| direct variables | 340,546 |
| direct rows | 1,561,199 |
| direct A nnz | 2,061,919 |
| `COSMO.assemble!` | 8.05 s |
| assembled cone count | 7 |

Outcome:

- `moment_relax` completed.
- direct sparse conic assembly completed.
- `COSMO.optimize!` started.
- **No COSMO iteration banner / no iter 1 before timeout.**

## Direct comparison

| Quantity | fat block | paper blocks | reduction |
| :--- | ---: | ---: | ---: |
| `moment_relax` wall | 28.19 s | 6.03 s | 4.67× |
| direct variables | 1,274,786 | 340,546 | 3.74× |
| direct rows | 8,154,607 | 1,561,199 | 5.22× |
| direct A nnz | 8,401,311 | 2,061,919 | 4.07× |
| `COSMO.assemble!` | 40.79 s | 8.05 s | 5.07× |

So the paper-style block construction does exactly what it should do structurally.
It is much smaller than the fat-block route.

## Verdict

Two blunt conclusions:

1. **Manual paper-style blocking is the right route inside the current API.**
   The structural savings are large and immediate.
2. **COSMO still does not solve this H4 periodic PQG benchmark.**
   Even after the reduction, `COSMO.optimize!` never reached iteration 1 in the
   timed runs.

That means the new benchmark path is useful, reproducible, and better-shaped than
both the fat-block route and the earlier TS benchmark — but it is still not enough
for COSMO on this asset.

## Bottom line

This branch now has a concrete benchmark path for exactly the formulation we wanted:

- explicit PQG basis,
- explicit particle-number constraints,
- explicit `Tr(²D)`, `Tr(²Q)`, `Tr(²G)` targets,
- manual paper-style blocks,
- direct COSMO setup numbers.

What it does **not** have is a converged H4 periodic COSMO energy.

The good news is that we are no longer guessing where the pain lives.
The pain is squarely in solver-side pre-iteration work, even after the formulation
is cleaned up.
