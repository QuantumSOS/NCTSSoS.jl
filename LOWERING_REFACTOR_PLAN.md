# Lowering refactor plan: `MomentProblem` → JuMP, two formulations

**Status**: proposal, not yet started.
**Companion doc**: `MOMENT_SOS_PIPELINE_ANALYSIS.md` (the why).
**This doc**: the how.

## Goal in one sentence

Factor the `MomentProblem → JuMP` lowering so that the same symbolic moment problem can be lowered into either the current free-moment-variable form (`:moment_variables`) or a PSD-block-variable form (`:psd_blocks`), with the second mode eliminating the bridge-induced 1×1 PSD-cone clutter that breaks BPSDP at H₂/Nk=2 and prevents H₄/Nk=2 from being solved at all.

## What this refactor is and isn't

| Is | Isn't |
|---|---|
| The infrastructure that lets BPSDP receive H₂/Nk=2 in a sensible shape. | A solution to H₄/Nk=2 by itself. |
| A clean cut between symbolic `MomentProblem` and JuMP construction. | A redesign of `MomentProblem` itself. `MomentProblem` stays as is. |
| A precondition for Phase 2 diagnostics on H₂/Nk=2. | The Phase 2 diagnostic. (See `PHASE2_PLAN.md`.) |
| A precondition for Phase 3 production runs at H₄/Nk=2. | Phase 3 itself (facial reduction, Mazziotti-CG, symmetry block-diagonalization). |
| Backward compatible — existing tests and COSMO path unchanged. | A new public API. Defaults preserve current behavior. |

## Scope

In scope:

- Extract current direct moment lowering into a builder function.
- Add `formulation = :moment_variables | :psd_blocks` and `representation = :real | :complex`.
- Implement `:psd_blocks, representation = :complex` for complex algebras (Pauli / Fermionic / Bosonic).
- Pivot discovery as a deterministic, library-internal procedure factored out of `demos/SDPALibsdpExport.jl`.
- An H₂/Nk=2 bridge-count smoke test on HAI as the acceptance gate.

Deferred to follow-ups:

- `:psd_blocks, representation = :real` (only needed if a backend can't take `HermitianPSDCone()` and we still want block-first).
- `SOSDualForm` symbolic dual + reusable SOS lowering. Tracked but not blocked on this PR.
- `ModelBuildResult` wrapper struct. Use a closure for extraction first; promote to a struct only when a regression test or external caller demands named-field access.

Out of scope (explicitly):

- Solving H₄/Nk=2.
- Symmetry-aware block-diagonal moment matrix construction.
- BPSDP tuning (warm start, ALM schedule, rank parameter).
- Facial reduction / Mazziotti-CG integration.

## Architectural cut: the Resolver pattern

Both formulations do the same operation on every polynomial they see:

> Walk the term list. For each `(c, m)`, replace `m` with some JuMP affine expression standing in for `⟨m⟩`. Sum with coefficients.

The only thing that differs between the two modes is **what stands in for `⟨m⟩`**. Call that function the resolver:

```julia
resolver(key::CanonicalKey) :: JuMP.AbstractJuMPScalar
```

Polynomial substitution is then generic:

```julia
substitute(poly, resolver) =
    sum(c_t * resolver(symmetric_canon(m_t)) for (c_t, m_t) in poly.terms)
```

The constraint-emission loop, objective handling, and identity normalization are shared. The two formulations differ in *one function plus block-variable declaration*.

### Two resolver shapes

**Affine resolver** (current behavior, `formulation = :moment_variables`):

```julia
# representation = :real (split re/im scalars)
resolver(key) = y_re[idx[key]] + im * y_im[idx[key]]
# representation = :complex (JuMP-native complex scalars)
resolver(key) = y[idx[key]]
```

**Pivot resolver** (new, `formulation = :psd_blocks`):

```julia
# representation = :complex (HermitianPSDCone block variables)
resolver(key) = X[pivot[key].block][pivot[key].i, pivot[key].j] / pivot[key].phase
# representation = :real (real-lift PSD blocks of size 2n)
# Same shape, with all four real-lift bindings emitted (see Trap 2 in companion doc).
```

## Minimum API surface

```julia
build_jump_model(mp::MomentProblem;
    formulation    = :moment_variables,        # | :psd_blocks
    representation = :real,                    # | :complex
    orphan_policy  = :error,                   # | :aux_psd_free
) :: Tuple{JuMP.GenericModel, Function}        # (model, extract_monomap)

solve_moment_problem(mp, optimizer;
    silent = true,
    formulation    = :moment_variables,
    representation = :real,
)
```

`extract_monomap` is a closure: call it after `optimize!(model)` to get `Dict{NormalMonomial, ComplexF64}`. No wrapper struct.

`representation` is irrelevant for real algebras; their resolver is always real-valued and `representation` defaults to `:real` silently.

## Open contracts (decide before coding)

Each contract has a **proposed default**. Override only with a stated reason.

### C1 — Pivot selection rule

**Default**: deterministic, total-order tie-break:

```
For each canonical key α appearing as a single-term entry phase * monomial,
prefer the pivot with the lexicographically smallest tuple
    (constraint_idx, i, j)
under the invariant that constraints, rows, and columns are iterated in
their stored order. First match wins; do not revisit.
```

This is reproducible across runs and across Julia versions because we never iterate `Dict` keys for pivot selection — we iterate `mp.constraints` (`Vector`) and matrix indices (`Int`).

### C2 — What counts as a pivot candidate

**Default**: `entry::Polynomial` qualifies as a pivot candidate iff after `simplify`:

- `length(entry.terms) == 1`,
- the single term is `(c, m)` with `c` exactly one of `±1, ±im` (no tolerance, no near-equality).

Anything else (multi-term, scaled, near-unity coefficient) does **not** qualify. This is conservative; we can loosen later if needed.

### C3 — Orphan definition

**Default**: an orphan is a canonical key that

- appears in `mp.objective` or in any `:Zero` matrix entry,
- but has no qualifying pivot candidate in any `:PSD` / `:HPSD` matrix.

Keys that appear *only* in PSD/HPSD entries but never as single-term pivot candidates are also orphans by this definition. The orphan_policy determines what happens:

- `:error` (default) — raise with the offending keys printed.
- `:aux_psd_free` — allocate a single auxiliary `1 × 1` Hermitian (or real) PSD block per orphan, mirroring `demos/SDPALibsdpExport.jl`. Logged at build time.

### C4 — Pivot-self skip rule

**Default**: when emitting bindings for block `b`, skip exactly the entry `(b, i*, j*)` whose pivot record is `(α, b, i*, j*, φ*)`. All other entries — including other entries that resolve to the same key α — emit a binding.

For native `HermitianPSDCone()` blocks, iterate the upper triangle only. The lower triangle is implied by the cone, not by binding equations.

### C5 — `:Zero` Hermiticity contract

**Default**: `MomentProblem.constraints` is a closed contract: every `:Zero` matrix is Hermitian by construction (this is what `_zero_constraint_components` in `moment_relax` already enforces). The builder asserts this on entry; non-Hermitian `:Zero` is a programmer error and raises.

If a future caller constructs `MomentProblem` by hand with non-Hermitian zero matrices, the assertion catches it before the lowering can quietly produce a wrong model.

### C6 — Identity normalization

**Default**: emit `resolver(key_of_identity) == 1` unconditionally as a regular constraint, *after* pivot discovery. Do not rely on the moment matrix `[1,1]` entry being identity (that's true today but it's algebra-dependent and not a load-bearing invariant). The constraint is one extra equality; cheap and explicit.

### C7 — Complex equality emission

**Default**: when emitting `X[i,j] == complex_affine_expr` against a `HermitianPSDCone` variable, let JuMP auto-split into real and imaginary equalities. Do not split manually. JuMP's split is correct and gives MOI a chance to combine equalities downstream.

If a debug regression later requires per-equation labels, switch to manual split locally.

## Implementation order

1. **Refactor without behavior change.**
   Extract the existing complex and real moment lowering into `build_jump_model(mp; formulation=:moment_variables, representation=:real)`. `solve_moment_problem` becomes a thin wrapper that calls the builder, attaches the optimizer, optimizes, and returns `(objective_value, extract(), status)`. Existing tests must pass with no edits.

2. **Resolver scaffolding.**
   Introduce internal `AffineResolver` and a generic `substitute` that uses it. Confirm the refactored `:moment_variables, :real` mode goes through the resolver. This is purely a code-shape change; outputs identical.

3. **Pivot discovery.**
   Factor pivot search out of `demos/SDPALibsdpExport.jl` into a library-internal `discover_pivots(mp) :: Dict{CanonicalKey, Pivot}` honoring contracts C1–C3. Unit-test against tiny hand-constructed `MomentProblem`s including the four phase cases (`+1, -1, +im, -im`) and an explicit orphan case.

4. **`:psd_blocks, representation = :complex`.**
   Wire `PivotResolver` into the builder. Implement block-variable declaration with `HermitianPSDCone()` for `:HPSD` and `PSDCone()` for `:PSD`. Honor C4–C7. Reuse `substitute` from step 2 unchanged.

5. **HAI smoke test.**
   Add the H₂/Nk=2 bridge-count regression test (predicate below). Run on HAI per server policy in `AGENTS.md`. This is the acceptance gate: the refactor is "done" when the smoke test passes.

6. **Optional follow-ups, separate PRs.**
   - `:psd_blocks, representation = :real` if a non-Hermitian-capable backend needs it.
   - `SOSDualForm` + `sos_dual_form` + matching builder.
   - Promote closure-based extraction to `ModelBuildResult` struct if regression tests demand named fields.

## Acceptance test (the smoke test predicate)

The refactor lands when this test passes on HAI:

```julia
# test/relaxations/h2_nk2_bridge_smoke.jl
using NCTSSoS, BPSDP, MathOptInterface
const MOI = MathOptInterface

@testset "H2/Nk=2 :psd_blocks delivers BPSDP-favorable cones" begin
    mp = build_h2_nk2_moment_problem()  # from demos/h2_periodic_nk2_moment_sos.jl
    model, _extract = build_jump_model(mp;
        formulation    = :psd_blocks,
        representation = :complex,
    )
    set_optimizer(model, BPSDP.Optimizer)
    MOI.Utilities.attach_optimizer(model)  # force bridges to materialize

    backend = MOI.get(unsafe_backend(model), MOI.ListOfConstraintTypesPresent())
    counts  = Dict(t => MOI.get(unsafe_backend(model), MOI.NumberOfConstraints{t...}())
                   for t in backend)

    n_hermitian = sum(v for ((F,S),v) in counts
                      if S <: MOI.HermitianPositiveSemidefiniteConeTriangle; init=0)
    n_psd_1x1   = sum(v for ((F,S),v) in counts
                      if S <: MOI.PositiveSemidefiniteConeTriangle &&
                         dimension_is_1x1(F, S); init=0)

    n_symbolic_blocks = count(c -> c[1] in (:PSD, :HPSD), mp.constraints)

    @test n_hermitian == n_symbolic_blocks_hermitian(mp)
    @test n_psd_1x1   <= length(orphans(mp))   # zero unless :aux_psd_free is in use
    @test n_psd_1x1   <  100                   # absolute ceiling vs current 406,706
end
```

Three assertions:

- **Cone count match.** Number of `HermitianPositiveSemidefiniteConeTriangle` instances equals the number of `:HPSD` matrices in `mp.constraints`. No silent bridging back to real-lift.
- **No 1×1 clutter.** Scalar 1×1 PSD cones at most equal the orphan count under `:error` policy this is zero. Hard ceiling of 100 catches any regression that reintroduces the pathology.
- **Recorded baseline.** The current direct path produces 406,706 scalar 1×1 cones. The new path must produce <100. Anything in between means a bridge is misbehaving and the diagnostic must run before declaring success.

If the third assertion fails, the bridge stack is the villain and the refactor is incomplete. If the first fails, `HermitianPSDCone()` is being silently realified by some MOI bridge and we need to investigate which one.

## What this enables downstream

| Phase | Without this refactor | With this refactor |
|---|---|---|
| Phase 2 diagnostic on H₂/Nk=2 | Diagnostic runs against the wrong shape of SDP — the 406,706-cone pathology dominates BPSDP behavior and obscures the real block structure. | Diagnostic sees the actual symbolic blocks. Per-block $\mathcal{A}\mathcal{A}^*$ analysis is meaningful. |
| Phase 3 production at H₄/Nk=2 | Cannot even build the JuMP model in a shape BPSDP can solve. | Build is clean. Remaining work is symmetry reduction, Phase-2-recommended solve strategy, BPSDP tuning. |
| COSMO compatibility | n/a | Unchanged — `:moment_variables, :real` is still the default. |

## Known traps (pointer)

The five implementation traps are documented in `MOMENT_SOS_PIPELINE_ANALYSIS.md` under "Trap 1–5":

1. `Symmetric(M) in PSDCone()` ≠ `M in PSDCone()`.
2. Real-lift PSD-first needs four bindings per entry.
3. Pivot phases include `±1, ±im`; phase must be tracked.
4. Native `HermitianPSDCone()` survival is bridge-dependent — verify, don't assume.
5. `symmetric_canon` does not collapse $w \leftrightarrow w^\dagger$ for complex algebras; do not introduce that collapse without a parallel rewrite of every consumer.

These should be linked from PR descriptions when the relevant code lands.

## Estimated change footprint

- **New code**: `src/optimization/lowering.jl` (resolvers, pivot discovery, `build_jump_model`). ~300–500 LOC.
- **Modified**: `src/optimization/moment.jl` (`solve_moment_problem` becomes a thin wrapper), `src/optimization/interface.jl` (kwargs plumbed through `solve_sdp` / `cs_nctssos`).
- **New tests**: `test/relaxations/lowering.jl` (unit tests for resolvers and pivot discovery), `test/relaxations/h2_nk2_bridge_smoke.jl` (HAI gate).
- **Untouched**: `MomentProblem` itself, sparsity layer, `polyopt`, every algebra-specific simplification.

## Bottom line

This refactor is the smallest cut that gets H₂/Nk=2 (and by extension H₄/Nk=2) into a JuMP shape BPSDP can consume without bridge-induced clutter. It is a precondition for Phase 2 diagnostics meaning anything and for Phase 3 production runs being possible. It does not solve H₄/Nk=2 by itself; symmetry reduction and solver tuning still come after.

Ship the resolver-based builder, the pivot discovery with strict orphan errors, the `:psd_blocks` mode for complex algebras, and the HAI smoke test. Defer the wrapper struct, the SOS dual form, and the real-lift PSD-block mode until they earn their keep.
