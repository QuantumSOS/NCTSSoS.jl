# H4 periodic (Nk=2) — order-2 V2RDM with term sparsity + COSMO — diagnostic

**Date:** 2026-04-21
**Goal:** Separate preprocessing wall time from actual ADMM wall time for the
primal-moment (`--dualize=false`) route, with `--max-iter=10`.

## Executive verdict

**The bottleneck is NOT the symbolic pipeline and NOT the JuMP model build.
It is entirely inside `optimize!(jump_model)` — specifically, in COSMO's
pre-iteration setup path (MOI → COSMO model transfer and/or initial KKT
symbolic factorization). COSMO never reached iteration 1 in either
`decompose=true` or `decompose=false` runs (wall times of 20+ and 17+ minutes
respectively, both aborted). With `--max-iter=10` we would have expected
10 iterations in seconds once the inner ADMM loop started.**

Primal-moment + COSMO, via JuMP, is **not viable as currently configured**
for this problem. The 3.2 min of symbolic work and 14 s of JuMP assembly
are fine — the hang is after that.

## Probe 1 — `decompose = true` (COSMO default)

Raw log: `demos/results/h4_periodic_term_sparsity_cosmo_benchmark_decomposeTRUE_STUCK.log`

| Stage | Wall (s) | Alloc | Status |
| :--- | ---: | ---: | :--- |
| build polyopt | 46.64 | 149.97 GiB | ✅ done |
| compute_sparsity (TS/MMD) | 139.71 | 11.16 GiB | ✅ done |
| moment_relax (symbolic) | 7.73 | 10.12 GiB | ✅ done |
| primal moment JuMP build | 13.67 | 106.42 GiB | ✅ done |
| optimize! (COSMO) | **>1302** | — | ⚠️ aborted before any output |

Symbolic pipeline total: **207.75 s ≈ 3.5 min**.
JuMP build: **13.67 s**.
COSMO `optimize!` reached: **nothing printed** — not even COSMO's banner.

## Probe 2 — `decompose = false`

Raw log: `demos/results/h4_periodic_term_sparsity_cosmo_benchmark_decomposeFALSE_STUCK.log`

| Stage | Wall (s) | Alloc | Status |
| :--- | ---: | ---: | :--- |
| build polyopt | 44.27 | 149.97 GiB | ✅ done |
| compute_sparsity (TS/MMD) | 142.22 | 11.16 GiB | ✅ done |
| moment_relax (symbolic) | 7.54 | 10.12 GiB | ✅ done |
| primal moment JuMP build | 14.55 | 106.42 GiB | ✅ done |
| optimize! (COSMO) | **>985** | — | ⚠️ aborted before any output |

Disabling COSMO's chordal decomposition (`decompose=false`) does **not** fix
the hang — it is not chordal analysis that is eating the time.

## Model size (invariant across probes)

| Quantity | Value |
| :--- | ---: |
| dense baseline PSD block | 2081 × 2081 |
| TS PSD blocks | 313 |
| TS largest block | 168 × 168 |
| sum(block_size²) | 983,593 |
| sum of upper-tri PSD scalar slots | 498,845 |
| unique moment variables (symbolic) | 26,817 |
| JuMP basis monomials (Hermitian embedding) | 85,941 |
| JuMP num_variables (= 2 × basis) | **171,882** |
| JuMP num_constraints (all types) | 16,831 |

Note: `basis length` (85,941) is larger than `unique_moment_variables`
(26,817) because the JuMP primal encodes a real/imag variable for every
`symmetric_canon(expval(m))` across the full `moment_problem.total_basis`,
while the unique-moment count is the distinct set of moment *values* that
appear inside the SDP constraint matrices. This is normal.

## Where is the time going?

After `@timed optimize!(jump_model)` starts, Julia spends 10+ minutes before
COSMO's verbose printer produces a single character. That rules out:

- ❌ Symbolic / NCTSSoS pipeline (already complete in 3.5 min, observed).
- ❌ JuMP model build (already complete in 14 s, observed).
- ❌ Per-iteration ADMM work (iteration 1 never started; max_iter=10 would
  have finished in seconds if we got there).
- ❌ COSMO chordal decomposition (`decompose=false` had same hang).

What remains, in decreasing likelihood:

1. **MOI → COSMO constraint copy.** JuMP's `optimize!` calls
   `MOI.Utilities.attach_optimizer`, which copies every affine expression in
   the 498,845 PSD-block scalar entries (plus the Zero-cone moment-equality
   entries) into COSMO's sparse constraint matrix. With a 2×2 Hermitian
   embedding, each entry is typically an affine combination of a handful of
   `y_re` / `y_im` variables, but the total nnz in `A` could be on the
   order of a few × 10⁶. Naive / bridge-heavy code paths in MOI can turn
   this into an O(nnz × something) process.

2. **KKT symbolic factorization.** QDLDL's initial symbolic analysis on the
   block-arrow KKT is usually fast, but for very wide PSD cones it can be
   slow. Worth checking after fixing (1).

3. **Constraint expansion in JuMP bridges.** The
   `@constraint(model, embedded in PSDCone())` lines feed dense `Matrix{Any}`
   into JuMP's PSD bridge. JuMP might be rebuilding typed expressions
   through `MOI.Utilities`. `Matrix{Any}` is already flagged in the source
   comment as suboptimal (see `_substitute_complex_poly` in
   `src/optimization/moment.jl`).

## What this run *did* establish, cleanly

- Flushed output works. Previous runs that produced 0-byte logs were a
  buffering artifact, not a silent crash. We now see every stage in real
  time; if the solver were iterating it would show up.
- Symbolic → JuMP build is fast and reproducible: **~3.5 min + 14 s = 3.7
  min total** before the hang.
- The primal-moment JuMP model is **171,882 scalar vars, 16,831
  constraints** — modestly sized for an SDP, but the PSD data it exposes
  (498k upper-tri scalars across 313 blocks, each with 2× Hermitian
  embedding) is where the transfer cost lives.
- The hang is deterministic and independent of COSMO's chordal-decomposition
  setting.

## Next experiments (ranked)

1. **Bypass JuMP and build COSMO's model directly.** Construct
   `COSMO.Constraint`s and `COSMO.Model` from the raw moment problem,
   `assemble!` it, and measure `assemble!` vs `optimize!` separately.
   COSMO's `assemble!` is the natural split: if `assemble!` is fast and
   `optimize!` still hangs, the pathology is in the solver; if `assemble!`
   hangs, the pathology is in the transfer.

2. **Try a different SDP solver with the same JuMP model.** Clarabel (in the
   docs env already) uses a very different data-ingestion path than COSMO.
   If Clarabel reaches iteration 1 quickly, the issue is specifically the
   COSMO/MOI bridge. Mosek (if a license is available) is the strongest
   sanity check.

3. **Typed expression matrices.** Replace the `Matrix{Any}` used for
   `mat_re`, `mat_im` with
   `Matrix{JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}` — both locally
   in the probe and upstream in `_solve_complex_moment_problem`. This
   removes one layer of boxing / `Any`-dispatch in the subsequent
   MOI-copy path.

4. **Shrink the problem as a proxy.** Run the same pipeline on a smaller
   problem (e.g. H₂ Nk=2 or H₄ Nk=1) that will certainly start iterating,
   and extrapolate wall time by problem size. This is cheap and separates
   "big but finite" from "hangs forever."

5. **Revisit the dual route with real awareness.** Originally rejected
   because the SOS dual ballooned to 2,006,058 JuMP vars / 172,195
   constraints — mostly lifted matrix multipliers, not moments. But if the
   pathology is in MOI bridging of PSD cones rather than raw problem size,
   the dual (which has fewer PSD cones but more scalar constraints) may
   still be more tractable through JuMP. Worth measuring once the MOI/COSMO
   split is instrumented.

## Artifacts

- `demos/h4_periodic_term_sparsity_cosmo_benchmark.jl` — updated with
  flushed logging, split JuMP build vs `optimize!` stages, and
  `--decompose` CLI flag.
- `demos/results/h4_periodic_term_sparsity_cosmo_benchmark.log` — latest
  run (decompose=false, aborted at 20 min).
- `demos/results/h4_periodic_term_sparsity_cosmo_benchmark_decomposeTRUE_STUCK.log`
  — saved from the earlier `decompose=true` probe.
- `demos/results/h4_periodic_term_sparsity_cosmo_benchmark_decomposeFALSE_STUCK.log`
  — saved from the `decompose=false` probe.
