# H4 periodic Nk=2 — spin-resolved native V2RDM Mosek probe

**Date:** 2026-04-22  
**Script:** `demos/h4_periodic_native_v2rdm_mosek_benchmark.jl`  
**Log:** `demos/results/h4_periodic_native_v2rdm_mosek_benchmark.log`

## Goal

Take the **spin-resolved native bridge** model for periodic H₄ `Nk = 2` and ask
one simple question:

> Does Mosek at least behave like a real conic solver on this model shape?

The comparison point is COSMO, which still stalls before even printing its
banner on both the K-only and spin-resolved bridge models.

## Model used

Selected refinement:

- `²D`, `²Q` blocked by `(K, 2S_z(pair))`
- `²G` blocked by `(Δk, Δ(2S_z))`

Block sizes:

- `²D`, `²Q`: `[56, 128, 56, 64, 128, 64]`
- `²G`: `[128, 256, 128, 128, 256, 128]`

Native / solver-facing counts:

- native Hermitian `²D` variables: `47,232`
- D-only moment-vector real variables: `94,464`
- Hermitian PQG blocks: `18`
- scalar equalities: `5`
- real-lift PSD rows: `584,160`
- total rows including equalities: `584,165`

Relative to the K-only native bridge, this spin refinement removes about:

- **`2.61×`** of the native free variables
- **`2.64×`** of the PSD rows

So if Mosek cannot even get traction here, the K-only model would be worse, not
better.

## Run configuration

Command used:

```bash
julia --project demos/h4_periodic_native_v2rdm_mosek_benchmark.jl \
    --refinement=spin_resolved --time-limit=30
```

Observed build stats from the log:

- outer JuMP build wall: `57.22 s`
- transient allocation during build: `212.06 GiB`
- internal native build total: `9.65 s`

Constraint summary from JuMP:

- variables: `47,232`
- constraints: `23`
  - `6` variable PSD blocks for `²D`
  - `12` affine PSD blocks for `²Q(D)` and `²G(D)`
  - `5` scalar equalities

## What Mosek printed

Unlike COSMO, Mosek immediately recognized the conic model and reported:

- objective sense: `minimize`
- problem type: `CONIC`
- constraints: `47,733`
- affine conic constraints: `12 (489,200 rows)`
- matrix variables: `6 (scalarized: 94,960)`

Then it proceeded to:

- `Optimizer started.`
- `Presolve started.`
- run the linear dependency checker,
- run the eliminator,
- and finish with:

```text
Presolve terminated. Time: 0.04
```

That is already materially better than COSMO, which never got as far as a
single solver-generated line.

## What Mosek did **not** print

Within the available outer wall budget, the run did **not** reach:

- a printed interior-point iteration line,
- a final termination summary,
- a solver-reported time-limit exit.

So the behavior is:

- **banner / model summary:** yes
- **presolve:** yes
- **iteration table:** no, not within the tested wall budget
- **clean timeout message:** no, not before the external timeout killed the run

The slightly annoying detail is that the internal `--time-limit=30` did not show
up as a clean solver stop before the external wall clock killed the process.
That strongly suggests the expensive work is still happening in the heavy
post-presolve setup before Mosek gets to a visible iteration checkpoint.

## Comparison against COSMO

| Solver | Banner | Presolve / setup progress | Iteration table | Outcome |
| :--- | :---: | :---: | :---: | :--- |
| COSMO on spin-resolved bridge | no | no | no | stalled before visible solver output |
| Mosek on spin-resolved bridge | yes | yes | no | parsed model and finished presolve before wall timeout |

That is not a win on solve time.
But it **is** a win on basic solver viability.

## Verdict

Bluntly:

1. **Mosek is the first solver here that actually engages with the model.**
   It prints the conic summary and gets through presolve.
2. **The spin-resolved bridge was worth doing.**
   Without that `2.6×` cone reduction, this probe would only be uglier.
3. **The benchmark is still hard.**
   Even Mosek did not reach a printed interior-point iteration under the tested
   wall budget.

So the ranking stays the same:

- COSMO: wrong tool
- Mosek: plausible tool, but still expensive
- the next improvement has to come from either more structure or more patience,
  not from pretending the solver problem is solved.
