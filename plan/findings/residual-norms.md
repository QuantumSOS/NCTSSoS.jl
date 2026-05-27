# Finding F1 — COSMO / SCS / BPSDP residuals are not on the same scale

**Status:** confirmed from source, 2026-05-14, while H₂/Nk=2 BPSDP "failed to converge".

**TL;DR.** The three iterative solvers in our sweep report primal/dual residuals using **different norms** and **different normalization conventions**. Treating their `eps_abs = 1e-5` tolerances as comparable is wrong: at the same numeric `1e-5` setting, **BPSDP is being asked for roughly √m × stricter precision than COSMO and SCS**, where `m ≈ 40 000` is the number of equality constraints in H₂/Nk=2.

The "BPSDP plateaus at 2.8 × 10⁻⁵" observation that drove our first retune attempt is largely an artifact of this mismatch.

---

## What each solver actually computes

### COSMO — `r_prim`, `r_dual`

Source: `COSMO/src/residuals.jl::calculate_residuals!`

```julia
# r_prim = A x + s - b   (then unscaled via Einv)
# r_dual = P x + q - A' μ   (then unscaled via Dinv, cinv)
return norm(r_prim, Inf), norm(r_dual, Inf)
```

So:

$$
r_{\text{prim}} = \|A x + s - b\|_{\infty}, \qquad
r_{\text{dual}} = \|P x + q - A^{\!\top} \mu\|_{\infty}
$$

both in **ℓ∞-norm**, in original problem units.

Termination is a mixed absolute + relative test (`residuals.jl::isapprox_feasible`):

$$
r_{\text{prim}} < \varepsilon_{\text{abs}} + \varepsilon_{\text{rel}} \cdot \max_{\text{prim}}
$$

with `eps_abs = eps_rel = 1e-5` in our config, and `max_norm_prim ≈ ‖b‖∞` for unit-magnitude data → effective threshold ≈ 2 × 10⁻⁵ in ∞-norm.

### SCS — `res_pri`, `res_dual`

Source: `scs/src/scs.c::compute_residuals` (and `glbopts.h`, `#define NORM SCS(norm_inf)`).

After dividing out the homogeneous-self-dual `tau`:

$$
\text{res\_pri} = \|A\hat x + \hat s - b\|_{\infty}, \qquad
\text{res\_dual} = \|P\hat x + A^{\!\top}\hat y + c\|_{\infty}
$$

In original problem units after un-normalization. **Same norm and same definition as COSMO**.

Termination test (`scs.c::has_converged`):

$$
\text{res\_pri} < \varepsilon_{\text{abs}} + \varepsilon_{\text{rel}} \cdot \max(\|A x\|_\infty, \|s\|_\infty, \|b\|_\infty)
$$

plus identical conditions on `res_dual` and `gap`. Default `eps_abs = eps_rel = 1e-4`; we set both to `1e-5`.

### BPSDP — `primal_error`, `dual_error`

Source: `BPSDP.jl/src/solver.jl::_update_metrics!`

```julia
problem.Au!(dual_work, state.x)
@. dual_work = dual_work - problem.b
state.primal_error = norm(dual_work)        # Julia LinearAlgebra.norm → 2-norm

problem.ATu!(primal_work, state.y)
@. primal_work = primal_work + state.z - problem.c
state.dual_error  = norm(primal_work)
```

`norm(v::AbstractVector)` in Julia defaults to **ℓ₂-norm**. So:

$$
\text{primal\_error} = \|A x - b\|_2, \qquad
\text{dual\_error} = \|A^{\!\top} y + z - c\|_2
$$

Termination test (`solver.jl::_has_converged`):

$$
\text{primal\_error} \le \texttt{sdp\_error\_convergence} + \texttt{sdp\_error\_convergence\_rel} \cdot (1 + \|b\|_2)
$$

Our config has `sdp_error_convergence_rel = 0.0` (default), so the threshold is **purely absolute** at `1e-5` — no data-magnitude scaling. Different from COSMO and SCS, which both use mixed abs+rel.

---

## Why the mismatch matters quantitatively

For any vector `v ∈ R^m`:

$$
\|v\|_{\infty} \;\le\; \|v\|_2 \;\le\; \sqrt{m}\, \|v\|_{\infty}.
$$

Near convergence, first-order methods (ADMM and boundary-point alike) tend to produce residuals that are roughly **uniformly distributed across constraints**. In that regime:

$$
\|v\|_2 \approx \sqrt{m} \cdot \|v\|_{\infty}.
$$

For H₂/Nk=2, the JuMP/MOI model has **40 128 complex equality constraints + 8 real equality constraints + 6 HPSD cones**. Real-valued residual dimension `m ≈ 40 136`, so `√m ≈ 200`.

### Implication 1 — same numeric tolerance asks for different precision

Setting `sdp_error_convergence = 1e-5` (BPSDP) and `eps_abs = 1e-5` (COSMO/SCS) does **not** ask for the same precision. The effective ∞-norm threshold each solver enforces is:

| Solver | Reported residual | Numeric tol | Effective ∞-norm threshold (uniform-residual estimate) |
|---|---|---|---|
| COSMO | ‖·‖∞ | 1e-5 | **~ 2 × 10⁻⁵** |
| SCS   | ‖·‖∞ | 1e-5 | **~ 2 × 10⁻⁵** |
| BPSDP | ‖·‖₂ | 1e-5 | **~ 1e-5 / √m ≈ 5 × 10⁻⁸** |

**At the "same" `1e-5`, BPSDP is being asked for ~400× tighter precision than COSMO/SCS.**

### Implication 2 — observed BPSDP "plateau" is actually high accuracy

Original H₂/Nk=2 BPSDP run (μ_freq=500, dynamic_cg=true, cg_tol=1e-10):

```
primal_error = 2.8 × 10⁻⁵   (2-norm)
dual_error   = 2.8 × 10⁻⁵   (2-norm)
```

Translating to ∞-norm via `‖v‖∞ ≈ ‖v‖₂ / √m`:

$$
\|A x - b\|_{\infty} \;\approx\; \frac{2.8 \times 10^{-5}}{\sqrt{40136}} \;\approx\; 1.4 \times 10^{-7}.
$$

**That is ~100× tighter than COSMO/SCS at "OPTIMAL" status.** BPSDP wasn't stuck — it had solved the problem to much higher precision than the other two, while being asked to keep going.

### Implication 3 — the retune story changes

Our first BPSDP retune attempt (BPSDP.jl reference config: μ_freq=25, dynamic_cg=false, cg_tol=1e-8) gave:

```
primal_error = 6.6 × 10⁻⁵   (2-norm)   →  ∞-norm est ≈ 3.3 × 10⁻⁷
```

Still tighter than COSMO/SCS, just less tight than the original. Both BPSDP configs are **already past the COSMO/SCS effective accuracy**; the difference between them is at a precision level COSMO/SCS never reach.

---

## How to put the three solvers on the same footing

Three sound options. Pick exactly one for the published comparison, and state the choice explicitly in the writeup.

### Option A — match effective ∞-norm precision per case

For each case, compute `m = num_equality_constraints` from `sdp_size.by_constraint_type` and set

$$
\texttt{sdp\_error\_convergence}_{\text{BPSDP}} \;=\; \sqrt{m}\,\varepsilon
$$

where `ε = 1e-5` is the target ∞-norm precision. For H₂/Nk=2, that is BPSDP tol ≈ 2 × 10⁻³.

Pros: rigorously apples-to-apples assuming uniform-residual heuristic.
Cons: depends on the heuristic; will look like a "weird" BPSDP tolerance to readers without context.

### Option B — use BPSDP's relative-residual path

Set `sdp_error_convergence_rel = 1.0` and `sdp_error_convergence = 0` (or a small abs floor). BPSDP threshold becomes ≈ `1 + ‖b‖₂`, which is data-dependent and somewhat parallels COSMO/SCS's `eps_rel × max_norm`.

Pros: uses BPSDP's own intended relative-error machinery.
Cons: still mixes 2-norm vs ∞-norm; "comparable" is qualitative.

### Option C — same numeric tolerance, report converted values

Keep `1e-5` for all three solvers. In the result.json and the published table, also report a converted column:

- For BPSDP: `primal_error_inf_est = primal_error / sqrt(m)`
- For COSMO/SCS: `primal_error_2_est = primal_error_inf × sqrt(m)`

Pros: minimal driver changes; no risk of mis-tuning a solver.
Cons: relies on readers trusting the conversion; doesn't change the actual stopping behavior, so BPSDP still runs to its (over-tight) tolerance and hits max_iter.

### Option D — record both norms by reading raw residual vectors

Patch the driver to extract `state.x`, `state.y`, `state.z` from BPSDP's MOI backend, compute `A x - b` and `A' y + z - c` directly, and store both 2-norm and ∞-norm. Mirror for COSMO/SCS.

Pros: fully rigorous; no heuristic.
Cons: ~30 min of driver work per solver; depends on stable internals.

---

## Decision (recorded after this finding lands)

> *(this section to be filled in after the user chooses; this is the placeholder for D8 — solver tolerance equalization)*

---

## Sources

- `COSMO/src/residuals.jl`, function `calculate_residuals!`
- `cvxgrp/scs` — `src/scs.c::compute_residuals`, `include/glbopts.h::NORM`, `docs/src/algorithm/index.rst` "Termination criteria"
- `BPSDP.jl/src/solver.jl` — `_update_metrics!` (lines ~298–308), `_has_converged` (lines ~117–127), `_error_convergence_tolerance` (lines ~104–107)
- Julia stdlib: `LinearAlgebra.norm(::AbstractVector)` defaults to 2-norm
- Empirical: `output/solver-evidence/h2/nk2/{cosmo,scs,bpsdp}/result.json` from the 2026-05-14 sweep runs (pre-decision, before tolerance equalization)
