# H₄ / Nk ∈ {2,…,8} sector pruning — exact post-processing of the JuMP/MOI/BPSDP pipeline

> **Status**: planning, not yet started. **Scope**: `demos/` and `probes/` only; no
> changes to `src/`, no changes to `BPSDP.jl`. **Claim**: removes the
> particle-number-non-conserving (and S_z, K, parity) orphan moment leakage from
> the lowering pipeline at zero relaxation. The post-pruned SDP and the
> as-built SDP have **identical optimal value** to machine precision.
>
> Companion to `PHASE2_PLAN.md` (H₂/Nk=2 calibration) and
> `MOMENT_PROBLEM_ENRICHMENT_PLAN.md` (Phase 3 lowering refactor). This plan
> fixes the *production* H₄/Nk≥2 pipeline today, without waiting for either.

---

## 0. The leak in one paragraph

The H₄ active-space Hamiltonian is U(1)_N × U(1)_{S_z} × ℤ_{Nk} symmetric. Its
ground state is a closed-shell singlet at fixed (N, S_z, K) =
(2·Nk, 0, 0). For any monomial m with ΔN(m) ≠ 0, ΔS_z(m) ≠ 0, or
ΔK(m) ≠ 0 (mod Nk), ⟨m⟩ = 0 is a **theorem**, not a relaxation.

NCTSSoS's `moment_relax` does not know this. It builds Schouten HPSD blocks
(²D, ²Q, ²G, ¹D) per momentum/spin sector — those are internally
ΔN-balanced and fine — but it also generates one-sided moment-equality
localizing rows ⟨b_i† · g⟩ = 0 with row basis
`moment_eq_row_bases = ⋃_b schouten_blocks[b]`. When b_i ∈ ²D (ΔN(b_i) =
−2), b_i† · g produces ΔN = +2 monomials. These canonical keys are
registered in `MomentLinearData.moments`, have **no HPSD pivot**, and fall
into `MomentLinearData.free_keys`. The current demo lowers them via
`orphan_policy = :aux_psd_free`, which packs them into 33×33 aux HPSD
blocks. The MOI bridge stack then explodes those into hundreds of thousands
of scalar 1×1 PSD cones plus tautological equality rows of the form
"linear combination of zeros = 0". BPSDP receives the bloat and stalls.

H₂/Nk=2 already exhibits this in miniature: 18,786 free moments → 406,706
scalar 1×1 PSD cones (numbers from `MOMENT_PROBLEM_ENRICHMENT_PLAN.md`).
H₄/Nk=2..8 scales worse.

**The fix here** is exact post-processing inside the demo: detect
symmetry-violating orphans, pin them to zero in JuMP, and delete the meq
constraints that become 0 = 0 after the pinning. No source change, no
relaxation.

---

## 1. Goal and non-goals

### Goal

For H₄ active-space PQG V2RDM at Nk ∈ {2, …, 8}, ship a `demos/`-only
pipeline that:

1. Builds the MomentProblem exactly as `demos/h4_periodic_moment_sos.jl`
   does today.
2. Lowers it to JuMP via the existing `build_jump_model(...; orphan_policy=:free_variables)`
   path (one knob change vs. the H₂/Nk=2 demo, no source mod).
3. **Post-prunes** the JuMP model by:
   - fixing every symmetry-violating orphan free variable to 0;
   - deleting every equality constraint that becomes tautologically
     `0 == 0` after the fixing.
4. Hands the pruned model to `BPSDP.Optimizer` on `HAI`.
5. Reports BPSDP-side metrics (variable count, equality rows, 1×1 cone
   count, A-matrix shape, nnz, wall-clock) before vs. after pruning, so
   the reduction is auditable.
6. Asserts `objective_value(pruned) ≈ objective_value(unpruned)` to ≤
   1e-9 absolute on the H₂/Nk=2 calibration case (Phase-2 instance, where
   unpruned still solves) before applying the same script to H₄/Nk≥2.

### Non-goals

- No modifications to `src/optimization/moment.jl`, `src/optimization/moment_linear.jl`,
  `src/optimization/lowering.jl`, `src/optimization/sos.jl`, or any other file
  under `src/`.
- No modifications to `BPSDP.jl`.
- No modifications to `demos/h4_periodic_moment_sos.jl`,
  `demos/h2_periodic_moment_sos.jl`, `demos/h2_periodic_nk2_moment_sos.jl`, or
  any other already-existing file under `demos/` (frozen list lives in
  `PHASE2_PLAN.md §12.1`).
- Not facial reduction. Facial reduction generalizes this; this plan
  exploits one specific, *known* face (the symmetry sector) by construction.
- Not symmetry block-diagonalization of the HPSD blocks. Those are already
  ΔN-homogeneous per Schouten block; the leak is in meq rows, not in
  HPSD-block off-diagonals.
- Not a substitute for the Phase 3 `MomentProblem` enrichment refactor. The
  refactor's `:psd_blocks` formulation kills MOI bridge clutter; this plan
  kills the *upstream* moment leak. Both are wins, on different layers.

---

## 2. What gets pruned, precisely

A canonical moment key `α` admits a **symmetry signature**

```
σ(α) = (ΔN(α), 2·ΔS_z(α), ΔK(α) mod Nk, parity(α))
```

computed from any representative monomial `m ∈ key_to_monomial[α]` via the
existing `basis_label(m; nk, norb)` helper in `demos/h4_periodic_moment_sos.jl`
(extended to also report fermion parity ΔN mod 2). The signature is
well-defined: σ does not depend on the choice of representative because
canonicalization within a fermionic CAR algebra preserves
(ΔN, ΔS_z, ΔK, parity).

A key `α` is **physically zero** iff σ(α) ≠ (0, 0, 0, even). The closed-shell
ground-state target sector is (0, 0, 0, even). For an excited-state
calculation in a different sector, replace the zero target with the
target sector signature; the pruning logic is identical.

The pruning operation is the pair (P1) + (P2):

**(P1) Variable fixing.** For every `α ∈ mp.linear.free_keys` with
σ(α) ≠ target sector, set the corresponding JuMP free variable(s) to 0:

```julia
JuMP.fix(orphan_re_var[α], 0; force = true)
JuMP.fix(orphan_im_var[α], 0; force = true)   # complex case only
```

**(P2) Tautological constraint deletion.** For every JuMP equality
constraint emitted by the lowering layer (zero constraints from
`mp.linear.zero_constraints`), check whether *every* variable referenced
by the constraint either (a) is fixed to 0 by (P1), or (b) appears with
zero coefficient in the constraint after substitution. If yes and the
constant term is zero, the constraint is `0 == 0`; `JuMP.delete(model, con)`.
If yes but the constant term is nonzero, that's a contradiction — the
symmetry classifier is wrong (or the user picked the wrong target sector).
Abort with a diagnostic; do not silently pass.

(P1) is exact because σ(α) ≠ target ⇒ ⟨α⟩ = 0 in the target sector.
(P2) is bookkeeping: a constraint over only zero-pinned variables
contributes nothing to the SDP.

> **Mathematical equivalence**: the unpruned SDP P_full and the pruned SDP
> P_pruned have the same feasible set restricted to the target sector and
> the same optimal value. P_full has additional feasible directions only in
> the orthogonal sectors, which the symmetry-violating moments parametrize;
> those directions cannot lower the objective because the objective lives
> in the (0,0,0,even) sector. P_full's nonzero solutions in those
> directions are equivalent to fictitious zero-energy degeneracies that
> P_pruned removes by construction.
>
> A two-page derivation belongs in the eventual probe artifact under
> `output/phase2/h4_nk_pruning/` (§9). It is not load-bearing for this
> plan's *implementation*; it is load-bearing for the *claim of no
> relaxation*. The smoke test in §6 is the empirical guard.

---

## 3. How orphans get from JuMP into a fixable handle

The current `_build_complex_psd_block_model` (with `orphan_policy=:free_variables`)
calls `_declare_free_orphan_moments!`, which creates two
`@variable(model, orphan_re[1:n_orphans])` /
`@variable(model, orphan_im[1:n_orphans])` arrays in **the same iteration order
as `L.free_keys`**. The returned dictionary `free_values[key] = orphan_re[idx] + im * orphan_im[idx]`
is *not* exposed by `build_jump_model`; only `model` and `extract_monomap`
escape.

So the demo cannot retrieve the orphan variables directly. Two robust
recovery paths, both source-free:

### 3.1 Path A (preferred): retrieve through `extract_monomap` indirection

After `build_jump_model`, call `extract_monomap()` to get
`Dict{key, JuMP.AffExpr or VariableRef}` (the resolver maps each
canonical key to its JuMP scalar). For a free key, the resolver returns the
free `VariableRef` (real algebra) or `aff = re + im*im_var` (complex).
Decompose `aff` via `JuMP.linear_terms(aff)` to recover both `re_var` and
`im_var`. This is public JuMP API.

```julia
model, extract_monomap = build_jump_model(mp;
    formulation     = :moment_variables,   # H₂/Nk=2 default; complex algebra path
    representation  = :real,               # current Hermitian-as-2n×2n lift
    orphan_policy   = :free_variables,
)

# extract_monomap is closed over (y_re, y_im) and the resolver
key_to_jump = extract_monomap()           # Dict{K, ComplexF64} *before* solve
# That returns *values*, not variables. Need to call before optimize!,
# but extract_monomap returns value(...) — useless pre-solve.
```

That fails because `extract_monomap` calls `JuMP.value(...)`. So:

### 3.2 Path B: re-derive orphan handles by symbolic substitution

The lowering layer's `AffineResolver` keys the free variable map by the
canonical key. We don't have that resolver, but we can rebuild the same
mapping deterministically from `mp.linear`:

```julia
# Pull cached linear data
L = mp.linear

# JuMP namespaces orphan_re and orphan_im as anonymous VariableRefs.
# Iterate model variables in declaration order; the orphan vars are the
# 2 * length(L.free_keys) trailing variables.
all_vars = JuMP.all_variables(model)
n_orph   = length(L.free_keys)

# In the complex path, free vars come AFTER y_re (length |moments|),
# AFTER y_im (length |moments|). So:
n_y      = length(L.moments)
@assert length(all_vars) == 2 * n_y + 2 * n_orph
orphan_re = all_vars[2*n_y + 1            : 2*n_y + n_orph]
orphan_im = all_vars[2*n_y + n_orph + 1   : 2*n_y + 2*n_orph]
free_handle = Dict(key => (orphan_re[i], orphan_im[i])
                   for (i, key) in enumerate(L.free_keys))
```

This is fragile if the lowering ever changes its declaration order. To
**lock the contract**, add a per-build assertion in the demo:

```julia
# spot-check three orphan keys: their JuMP variable's coefficient in the
# resolver should be 1 (real part) or i (imag part), nothing else.
for key in first(L.free_keys, min(3, length(L.free_keys)))
    handle = free_handle[key]
    # We can verify by comparing extract_monomap()'s post-solve output to
    # value(handle[1]) + im * value(handle[2]) on a trivial probe instance.
end
```

If the contract ever breaks, the assertion catches it before BPSDP runs.

### 3.3 Real algebra and `:psd_blocks` path

Real algebras: only `orphan_re`. Identical recipe with `n_orph` instead of
`2*n_orph`.

`:psd_blocks` path with `orphan_policy=:free_variables`: same recipe,
because `_build_complex_psd_block_model` also calls
`_declare_free_orphan_moments!`. The only path that differs is
`:aux_psd_free`, which constructs aux HPSD blocks instead of free
scalars; **avoid that path entirely** for this plan. (Use
`:free_variables` everywhere, including in the H₂/Nk=2 demo if we choose
to update it — but per §12.1 the H₂/Nk=2 demo is frozen, so this plan
only retargets the H₄ flow.)

---

## 4. Symmetry signature implementation

### 4.1 Per-monomial labels

`basis_label(m; nk, norb)` already returns `(ΔN, K, two_Sz)` for the H₄
demo's fermionic CAR algebra. Extend it (in the new H₄ pruning script,
not in the existing demo) with a fourth field `parity = mod(ΔN, 2)`. ΔN is
already computed; parity is `iseven(ΔN)`.

### 4.2 Per-key labels

For a canonical key `α ∈ mp.linear.moments`, pick the representative
monomial:

```julia
m = mp.linear.key_to_monomial[α]
σα = basis_label(m; nk, norb)
```

`symmetric_canon(expval(m))` is invariant under conjugation and basis
permutation but **does not change ΔN, ΔS_z, or ΔK**, so σ is well-defined
on canonical keys.

### 4.3 Target sector

Closed-shell ground state of H₄ at Nk: target = `(ΔN=0, K=0, two_Sz=0, parity=:even)`.
Expose this as a CLI flag:

```
--target-sector=N=0,K=0,Sz=0,P=even        (default)
--target-sector=N=2,K=0,Sz=0,P=even        (would target a charge-doublet sector)
```

Defining the sector at the script level — not inside NCTSSoS — keeps the
plan source-free.

### 4.4 Classifier output

`mp.linear.free_keys` is partitioned into:

- `expected_zero`: σ(α) ≠ target. Pin to 0.
- `physical_unexpected`: σ(α) == target but α has no HPSD pivot. **This
  shouldn't exist.** If `length(physical_unexpected) > 0`, that is the
  load-bearing diagnostic for an algebraic redundancy that NCTSSoS's
  current pivot-discovery missed (the Phase 2 "excess kernel" scaled up
  to H₄). Emit warnings; continue with `:free_variables` lowering for
  these specific keys (no pinning); record them in the run report.

This makes the script *also* a Phase-3-ready audit: an H₄/Nk=2 run with
zero `physical_unexpected` orphans is positive evidence that Phase 2's
bet on no excess kernel survives the production scale.

---

## 5. Execution shape

### 5.1 New file (the only one)

`demos/h4_periodic_moment_sos_sector_pruned.jl` — single new script.
Frozen-list discipline mirrors `PHASE2_PLAN.md §12.1`.

Top-level flow:

```
1. parse_options(argv)
   - all H₄ demo flags (--integrals, --blocking, --include-1d, …),
   - plus --target-sector=…,
   - plus the BPSDP knobs from PHASE2_PLAN.md §12.4,
   - plus --output-dir=output/h4_pruning/<tag>,
   - plus --pruning=on|off|both for A/B mode.

2. data = build_h4_pqg_moment_problem(options)
   - identical to existing helper in demos/h4_periodic_moment_sos.jl;
   - DO NOT clone the helper, `include(joinpath(@__DIR__, "h4_periodic_moment_sos.jl"))`
     and reuse `build_h4_pqg_moment_problem`. The frozen-file rule forbids
     edits, not reuse.

3. model, extract_monomap = build_jump_model(data.moment_problem;
       formulation    = :moment_variables,   # current H₄ default
       representation = :real,
       orphan_policy  = :free_variables,     # NOT :aux_psd_free
   )

4. orphan_handles = recover_orphan_handles(model, data.moment_problem)
   - Path B from §3.2 with the spot-check assertion.

5. classification = classify_orphans(data.moment_problem, options.target_sector;
                                     nk = data.nk, norb = data.norb)
   - Returns NamedTuple{(:expected_zero, :physical_unexpected), …}.

6. apply_sector_pruning!(model, orphan_handles, classification)
   - JuMP.fix(...; force=true) for every key in expected_zero.
   - JuMP.delete(model, con) for every all-zero zero-cone constraint.

7. solve with BPSDP (same recipe as PHASE2_PLAN §12.6):
   set_optimizer(model, BPSDP.Optimizer)
   set_attribute(model, "max_iter", options.bpsdp_max_iter)
   …
   optimize!(model)

8. write artifacts to options.output_dir:
   - report.json  (BPSDP termination, walltime, iterations, residuals)
   - cones.json   (cone counts before vs after pruning)
   - kernel.json  (free_keys split, expected_zero, physical_unexpected)
   - aastar.json  (optional: dense SVD of A * A^T for blocks at small Nk)

9. main always runs (no --no-solve). If --pruning=both, runs (3)→(8) twice
   with on/off; both must produce the same objective to ≤ 1e-9 abs, else
   abort with diagnostic.
```

### 5.2 Probe driver

`probes/h4_nk_pruning_diagnostic.jl` — consumes the per-Nk `output/h4_pruning/<tag>/`
directories produced by the demo and emits an aggregate
`output/h4_pruning/summary.md` with one row per (Nk, blocking, target_sector)
combo:

| Nk | blocking | n_moments | n_orphans | n_expected_zero | n_unexpected | bpsdp_iters_unpruned | walltime_unpruned | bpsdp_iters_pruned | walltime_pruned | obj_diff |
|----|----------|-----------|-----------|-----------------|--------------|----------------------|-------------------|--------------------|-----------------|----------|

Rows where `n_unexpected > 0` highlight in red.

### 5.3 Execution contract — same as Phase 2

Per `PHASE2_PLAN.md §12.6`: all BPSDP runs on `HAI` via `easy-ssh`; all
Python via `uv` (none needed for this plan). Pin `BPSDP.jl` to its
multithread branch checkout on `HAI`. Outputs land at
`/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark/output/h4_pruning/`,
mirrored to local via `easy-ssh sync` before write-up.

---

## 6. Acceptance gate

### 6.1 Equivalence proof, empirically

Run `--pruning=both` on H₂/Nk=2 (the smallest case where unpruned still
fits). Assert:

```
abs(obj_pruned - obj_unpruned) ≤ max(1e-9, 1e-9 * abs(obj_unpruned))
```

If this fails, the symmetry classifier is broken. Do not promote to H₄.

### 6.2 Cone reduction at H₂/Nk=2

After pruning vs. unpruned (with `:free_variables`, NOT `:aux_psd_free`):

```
n_variables_pruned  ≤  n_variables_unpruned - 2 * length(expected_zero)
n_eq_constraints_pruned  <  n_eq_constraints_unpruned
```

Concrete target — derived from the H₂/Nk=2 numbers cited in
`MOMENT_PROBLEM_ENRICHMENT_PLAN.md`:

- Unpruned, `:aux_psd_free`: 18,786 free moments → 406,706 scalar 1×1 cones.
- Unpruned, `:free_variables`: 18,786 free moments, no 1×1 cones, just
  18,786 free scalars.
- Pruned: number of free scalars equals `length(physical_unexpected)`. Aim
  for **≤ 100** at H₂/Nk=2. If physical_unexpected isn't empty, that's an
  excess-kernel signal feeding Phase 2.

### 6.3 BPSDP behavior at H₄/Nk=2

This is the production smoke test. Compare unpruned `:free_variables` vs
pruned at H₄/Nk=2:

- Both reach `OPTIMAL`, or
- Unpruned fails (timeout / OOM / `INFEASIBLE_OR_UNBOUNDED` from a stalling
  inner CG) and pruned succeeds.

The second case is the hypothesis. Either outcome is publishable: the
first is "leak doesn't matter at this scale, fix is moot at Nk=2"; the
second is "leak is what's killing BPSDP at Nk=2, this fix unblocks it."

### 6.4 Scaling probe

If 6.3 succeeds, repeat at Nk = 3, 4, …, 8. Tabulate. Expected pattern:
unpruned solve time grows superlinearly (orphan count × meq row count
both ∝ Nk^2 × |basis|^2), pruned solve time grows polynomially in problem
size. The crossover point is the empirical answer to "at what Nk does
this matter?"

---

## 7. Why this is the right cut

- **No source change.** Lives entirely in `demos/` + `probes/`. The
  frozen-file discipline of `PHASE2_PLAN.md §12.1` is preserved.
- **No relaxation.** Pinning σ ≠ target moments to zero is a theorem in a
  symmetry-pinned target sector; an a-priori-known face of the feasible
  cone, the kind facial reduction would discover by SVD. We just write
  it down explicitly because we know the symmetry.
- **Falsifiable.** Acceptance gate §6.1 catches symmetry classifier bugs.
  Acceptance gate §6.3 catches "the leak isn't actually the bottleneck"
  (we'd want to know that too).
- **Composes with Phase 2.** The H₂/Nk=2 calibration in `PHASE2_PLAN.md`
  measures `A * A^T` after lowering. Re-run Phase 2's §3 diagnostic on
  the pruned model and compare spectra: the pruned `A * A^T` should have
  a *much* smaller kernel and a cleaner spectral gap. That validation is
  free — Phase 2 already wrote `probes/h2_nk2_aastar_diagnostic.jl`.
- **Composes with Phase 3.** Once `MOMENT_PROBLEM_ENRICHMENT_PLAN.md`'s
  `:psd_blocks, representation = :complex` ships, this plan's pruning
  step layers on top with no rewrite — both layers compose by
  intersection of zero sets.
- **Reveals what facial reduction would buy.** `physical_unexpected`
  measures "leak that survives all known symmetries." Empty at H₄/Nk=2
  means facial reduction (Paper 5) is not on the critical path for
  production. Non-empty means it is, and gives a concrete catalog of
  algebraic identities to add.

---

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Lowering layer reorders free variables relative to `L.free_keys` | §3.2 spot-check assertion catches it before BPSDP runs. |
| Bridge translates `JuMP.fix(var, 0; force=true)` into a separate `EqualTo` constraint instead of removing the variable | Verify via `MOI.get(unsafe_backend(model), MOI.ListOfConstraintTypesPresent())` after fixing. If it does, follow with a manual MOI-level pass: detect the trivial constraints and `MOI.delete` them before BPSDP attaches. Acceptance gate §6.2 surfaces this. |
| Wrong target sector silently pinned (e.g., open-shell triplet target with closed-shell pinning) | (P2) catches it: a constraint over only zero-pinned variables with nonzero constant term is a contradiction. The script aborts with the constraint origin (which meq row, which sector). |
| `_close_adjoint_keys!` has already added adjoint partners — pinning α without pinning α* is inconsistent for complex algebras | σ is anti-equivariant under adjoint: σ(α*) = (−ΔN, …). Closing under adjoint means if α has σ ≠ target and target = (0,0,0,even), then α* also has σ ≠ target. So both get pinned. Verify by enumerating `mp.linear.adjoint_key`. |
| H₄/Nk=2 unpruned still solves under `:free_variables` and the leak is < 1% of variables | The plan is still correct (it's exact), but the production benefit becomes "scaling cushion" rather than "unblock." Document and move on. |
| Demo introspection breaks under future JuMP versions | Same fragility risk as `MOMENT_PROBLEM_ENRICHMENT_PLAN.md`'s smoke test. Pin Manifest. |
| The H₄ Hamiltonian has Re/Im noise that violates ΔN exactly by ε | Already handled in `build_h4_nk2_hamiltonian` via `0.5*(ham + adjoint(ham))`. ΔN-violating *coefficients* would survive this only at floating-point noise; their contribution to the SDP is at most 10⁻¹³ × |obj|. Drop them: any term with |coef| ≤ 10⁻¹² in the **canonicalized** Hamiltonian gets zeroed before `polyopt(...)`. (Optional — the demo already does this implicitly via numerical-zero registration.) |

---

## 9. Artifact layout

```
demos/
  h4_periodic_moment_sos_sector_pruned.jl   # ONE new file (§5.1)

probes/
  h4_nk_pruning_diagnostic.jl               # aggregator + summary writer (§5.2)

output/h4_pruning/
  summary.md                                # aggregated table across Nk

  h4_nk2/
    report.json                             # BPSDP termination, walltime, iters
    cones.json                              # cone counts: pruned vs unpruned
    kernel.json                             # free_keys split, expected_zero, unexpected
    aastar_pruned.json                      # optional: post-prune A*A^T spectrum
    aastar_unpruned.json                    # optional: comparison

  h4_nk3/...
  h4_nk4/...
  ...
  h4_nk8/...                                # may not finish unpruned; that's data
```

---

## 10. First-day list

1. Read `extract_monomap` closure construction in `_build_complex_psd_block_model`
   and `_build_complex_moment_variable_model` to confirm the orphan
   variable declaration order (§3.2). Locked? OK. Else, pick a different
   public hook and update §3.2.
2. Write `demos/h4_periodic_moment_sos_sector_pruned.jl` skeleton:
   parse_options, include the existing H₄ demo, build, `:free_variables`
   lowering, `recover_orphan_handles`, `classify_orphans`, no pruning yet,
   solve with BPSDP, dump `report.json` + `kernel.json`. Get an
   unpruned-but-`:free_variables` baseline at H₂/Nk=2 (since H₂/Nk=2 demo
   is frozen, run this script with H₂/Nk=2 integrals just to validate the
   baseline; H₂/Nk=2 production runs continue to use the existing demo).
3. Add `apply_sector_pruning!`. Run `--pruning=both` on H₂/Nk=2. Verify
   acceptance gate §6.1.
4. Run on H₄/Nk=2. Verify §6.3.
5. Stub `probes/h4_nk_pruning_diagnostic.jl`. Run on the H₄/Nk=2 outputs.
6. Scale up to Nk=3, 4, … as time permits. Each Nk adds one row to
   `summary.md`.

After day 1: we know whether the leak is the production bottleneck or a
red herring at this scale.

---

## 11. Bottom line

The leak is structural and known: the H₄ Hamiltonian has U(1)_N × U(1)_{S_z}
× ℤ_{Nk} × parity symmetry, and one-sided moment-equality localizing rows
generate canonical moments outside the target sector that lowering
treats as free orphans. The fix is to pin those moments to zero (theorem,
not relaxation) and delete the resulting trivial equalities — entirely in
the demo layer, post-build, with no source modification. The acceptance
gate is `obj_pruned == obj_unpruned` to 1e-9; everything else is
production engineering. If the gate passes and BPSDP starts converging
at H₄/Nk≥2 where it didn't, that's the unblock. If BPSDP still struggles,
the spectral gap of `A * A^T` on the pruned system tells us where the
*real* bottleneck lives, feeding directly into Phase 2's decision table
and Paper 5's facial-reduction targets.
