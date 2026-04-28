# H4 periodic Nk=2 — D-only bridge / native V2RDM benchmark update

**Date:** 2026-04-22  
**COSMO script:** `demos/h4_periodic_native_v2rdm_jump_benchmark.jl`  
**Reviewed logs:**

- K-only native bridge: `demos/results/h4_periodic_native_v2rdm_jump_benchmark.log`
- spin-resolved native bridge: `demos/results/h4_periodic_native_v2rdm_jump_spin_resolved.log`

## What changed

The native bridge model now has two versions:

1. **K-only native bridge**
   - `²D`, `²Q` blocked only by pair momentum `K`
   - `²G` blocked only by momentum transfer `Δk`
2. **spin-resolved native bridge**
   - `²D`, `²Q` blocked by `(K, 2S_z(pair))`
   - `²G` blocked by `(Δk, Δ(2S_z))`

This is a **structural** refinement. We are not adding zero constraints after the
fact. We are shrinking the PSD blocks themselves.

## Block structure

### K-only bridge

- `²D`: `[240, 256]`
- `²Q`: `[240, 256]`
- `²G`: `[512, 512]`

Hermitian block layout:

- `[240, 256, 240, 256, 512, 512]`

### Spin-resolved bridge

Pair sectors for `²D` / `²Q`:

- `K = 0`: `2S_z = -2, 0, +2` → `[56, 128, 56]`
- `K = 1`: `2S_z = -2, 0, +2` → `[64, 128, 64]`

Particle-hole sectors for `²G`:

- `Δk = 0`: `Δ(2S_z) = -2, 0, +2` → `[128, 256, 128]`
- `Δk = 1`: `Δ(2S_z) = -2, 0, +2` → `[128, 256, 128]`

Hermitian block layout:

- `[56, 128, 56, 64, 128, 64, 56, 128, 56, 64, 128, 64, 128, 256, 128, 128, 256, 128]`

## Why this split is legitimate

For the pair spaces the split is straightforward: fixed `K` and fixed pair
`2S_z` sectors do not couple.

For `²G`, the right invariant is **spin transfer**, not naive spin label.
That is why we split into `Δ(2S_z) ∈ {-2, 0, +2}`.
We do **not** split the zero-transfer sector into separate `αα` and `ββ` blocks,
because they still couple through the `D` term in

```math
G_{pq,st} = δ_{qt} \, {}^{1}D_{ps} - {}^{2}D_{pt,sq}.
```

So the refinement is aggressive where it is mathematically justified, and no
further.

## Variable and row counts

| Formulation | Free-variable coordinates | Real vars | PSD rows | Equality rows | Total rows |
| :--- | ---: | ---: | ---: | ---: | ---: |
| current NCTSSoS paper-route | full PQG moment vector | `340,546` | `1,545,187` | `16,012` | `1,561,199` |
| K-only D-only bridge | complex `²D` entries only | `246,272` | — | `5` | — |
| K-only native Hermitian `²D` | Hermitian `²D` blocks | `123,136` | `1,543,136` | `5` | `1,543,141` |
| spin-resolved D-only bridge | complex `²D` entries only | `94,464` | — | `5` | — |
| spin-resolved native Hermitian `²D` | Hermitian `²D` blocks | `47,232` | `584,160` | `5` | `584,165` |

So the spin refinement buys, relative to the K-only native bridge:

- **`2.61×` fewer native variables** (`123,136 → 47,232`)
- **`2.64×` fewer PSD rows** (`1,543,136 → 584,160`)
- **`2.64×` fewer total rows** (`1,543,141 → 584,165`)

Relative to the current NCTSSoS paper-route, the spin-resolved native bridge buys:

- **`7.21×` fewer real variables** (`340,546 → 47,232`)
- **`2.67×` fewer total rows** (`1,561,199 → 584,165`)
- **`16,012 → 5` equality rows**

This is finally more than just deleting equality clutter.
The PSD burden drops for real.

## Actual JuMP model size for the spin-resolved bridge

Measured from the live JuMP build (`demos/results/h4_periodic_native_v2rdm_jump_spin_resolved.log`):

- **variables:** `47,232`
- **constraints:** `23`
  - `6` variable-in-cone Hermitian PSD blocks for `²D`
  - `12` affine Hermitian PSD constraints for `²Q(D)` and `²G(D)`
  - `5` scalar equality constraints

Constraint-type breakdown:

- `Vector{VariableRef} in HermitianPositiveSemidefiniteConeTriangle` → `6`
- `Vector{AffExpr} in HermitianPositiveSemidefiniteConeTriangle` → `12`
- `AffExpr in EqualTo{Float64}` → `5`

Internal model-assembly timings from that run:

- `¹D` contraction build: `0.087 s`
- `²Q` block build: `0.530 s`
- `²G` block build: `0.432 s`
- scalar constraints: `0.172 s`
- objective assembly: `5.004 s`
- internal native build total: `8.316 s`

The outer wall was still dominated by Julia compile / first-run overhead:

- **outer JuMP build wall:** `56.40 s`
- **transient allocation during build:** `210.37 GiB`

So the model construction is fine.
The solver handoff is still the sore spot.

## COSMO behavior

### K-only bridge

The earlier K-only run already showed the basic problem:

- no COSMO banner,
- no problem dimensions,
- no iteration table,
- no internal time-limit exit.

### Spin-resolved bridge

Command used:

```bash
julia --project demos/h4_periodic_native_v2rdm_jump_benchmark.jl \
    --refinement=spin_resolved --time-limit=30 --max-iter=10000
```

Observed outcome from `demos/results/h4_periodic_native_v2rdm_jump_spin_resolved.log`:

- the native model build finished,
- the script printed

```text
[stage] optimize! (COSMO) : running native JuMP model ...
---------------- begin COSMO output ----------------
```

- and then **nothing else** before the outer wall budget killed the process.

So even after reducing the solver-facing cone rows from `1.54e6` to `5.84e5`,
COSMO still does not reach visible setup progress.

## Verdict

Two things are now clear.

1. **Spin resolution is worth doing.**
   The bridge/native model gets materially smaller in the cones that matter.
2. **COSMO is still the wrong solver for this benchmark.**
   The cleaner model is better, but COSMO still stalls before visible progress.

That means the next serious probe is not “more optimism about COSMO”.
It is **Mosek on the spin-resolved bridge**.
