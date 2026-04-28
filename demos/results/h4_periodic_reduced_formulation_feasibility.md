# H4 periodic Nk=2 — reduced-formulation feasibility audit

**Date:** 2026-04-21  
**Purpose:** Continue the solver investigation by answering the only question that still matters:

> **Can the current NCTSSoS codebase support a materially smaller H4 periodic formulation without rewriting the whole package? If yes, where are the hooks, and what is the shortest honest path?**

This is a codebase audit, not another solver run.

## Short verdict

**Yes, the codebase already has most of the hooks needed for a reduced formulation.**

But there is one important catch:

**The real lever is not a naïve “pair basis” by itself. The real lever is symmetry-adapted basis construction and/or pre-TS sector partitioning.**

A dumb pair basis over the full 32-mode spin-orbital registry barely shrinks anything. If we want paper-like scale, we must encode momentum/sector structure *before* the solver sees the giant generic moment system.

## What already exists in the repo

### 1. Custom `moment_basis` is already first-class

The hook is real, not hypothetical.

Relevant code:

- `src/optimization/interface.jl`
  - `SolverConfig(; optimizer, order=0, moment_basis=nothing, ...)`
  - `_resolve_relaxation_spec(...)`
- `src/optimization/sparsity.jl`
  - `correlative_sparsity(pop, moment_basis::AbstractVector, elim_algo)`
  - `_normalize_moment_basis(...)`
  - `_filter_basis_to_clique(...)`

That means the package can already accept an explicit basis instead of the default full order-based basis.

This is the cleanest existing entry point for any H4-specific reduction.

### 2. `moment_eq_constraints` already work for ordinary polynomial problems

Relevant code:

- `src/optimization/problem.jl`
  - `polyopt(...; moment_eq_constraints=...)`
- `src/optimization/moment.jl`
  - `_add_moment_eq_constraints!(...)`

So particle-number fixing can stay exactly where it is for the H4 benchmark. No new sector-constraint API is needed just to keep `N_up = 4`, `N_dn = 4`.

### 3. `moment_relax` already consumes per-block bases independently

Relevant code:

- `src/optimization/moment.jl`
  - `moment_relax(...)`
  - inside it: `for ts_sub_basis in term_sparsity.block_bases`
  - `_build_constraint_matrix(poly, ts_sub_basis, cone)`

This matters because it means **symmetry partitioning can be inserted upstream of `moment_relax`** and the downstream symbolic assembly does not need architectural surgery. The symbolic builder already knows how to process one block basis at a time.

### 4. `dualize=false` is already a supported public solve path

Relevant code:

- `src/optimization/interface.jl`
  - `solve_sdp(moment_problem, optimizer; dualize::Bool=true)`
  - `cs_nctssos(pop, solver_config; dualize::Bool=true)`

So there is no missing primal-toggle problem. The benchmark already used the right public escape hatch.

## The useful correction: naïve pair basis is not enough

There is a tempting bad idea floating around:

> “Just use a pair basis and the H4 problem will collapse to paper scale.”

Not so fast.

For the current **32-mode** spin-orbital registry, the naïve full pair-style basis size is:

- identity: `1`
- annihilation pairs `a_p a_q`: `C(32,2) = 496`
- creation pairs `a†_p a†_q`: `496`
- number / particle-hole words `a†_p a_q`: `32 × 32 = 1024`

Total:

- **naïve pair basis = `2017` elements**

Compare with the current full dense order-2 basis:

- **full order-2 basis = `2081` elements**

That is a reduction of only **64 basis elements**. In other words: basically noise.

So a raw spin-orbital pair basis is **not** the paper's block-reduced V2RDM formulation. The real savings in the paper come from **symmetry-adapted pair sectors** and from handling `D/Q/G` objects in their own natural block structure, not from merely deleting single-operator rows.

## Why the codebase still looks promising

Because the architecture is already close to what we need.

### Existing capability A — explicit basis injection

This can support:

- momentum-labeled basis elements,
- sector-restricted bases,
- chemistry-specific reduced bases,
- hand-built benchmark bases for H4 periodic.

### Existing capability B — blockwise symbolic assembly

This means a sector partitioning step can happen before `term_sparsities(...)` or before `moment_relax(...)`, without rewriting the whole solver layer.

### Existing capability C — docs already describe the exact conceptual levers

There is already a relevant docs page:

- `docs/src/examples/literate/periodic_v2rdm_block_structure.jl`

That page is not fluff. It already lays out the core options:

- explicit `moment_basis`,
- symmetry-adapted partitioning before TS,
- primal moment instead of dual blow-up,
- long-term boundary-point solver story.

So the conceptual design work is largely done. What is missing is turning that roadmap into a small, concrete implementation sequence.

## What is *not* in place yet

### 1. No first-class symmetry partition hook before TS

Today `compute_sparsity(...)` ultimately does:

- build clique basis,
- call `term_sparsities(init_act_supp, cons, mom_mtx_bases, localizing_mtx_bases, ts_algo)`.

There is no public hook like:

- `sector_tag(monomial) -> quantum_numbers`
- `partition_basis_by = sector_tag`

before the TS graph is built.

That is the most useful missing extension point.

### 2. No out-of-the-box H4-specific reduced basis helper

There is no helper today that builds something like:

- momentum-sector pair operators,
- symmetry-adapted 2-RDM / 2-hole-RDM / particle-hole bases,
- or an H4 periodic benchmark-specific reduced moment basis.

So if we want a reduced H4 formulation, we must build that helper ourselves.

### 3. The typed-expression comment in `moment.jl` is misleading for this path

Probe 3 showed that the current substitution path yields `QuadExpr`, not `AffExpr`, on the H4 primal builder.

So any future cleanup around typed JuMP expression matrices needs to start by fixing or explaining that type behavior. That is separate from the H4 scaling problem, but it is real technical debt.

## The shortest honest path forward

### Step 1 — add a sector-partition hook before TS

**Best leverage, smallest meaningful change.**

Target area:

- `src/optimization/interface.jl`
- possibly a helper in `src/optimization/sparsity.jl`

Sketch:

- add an optional `basis_partition` or `sector_tag` field to `SolverConfig`,
- when present, partition each clique moment basis and localizing basis by that tag,
- run `term_sparsities(...)` inside each partition,
- concatenate the resulting block bases.

Why this is attractive:

- downstream moment assembly already supports per-block bases,
- no new solver backend work,
- no need to redesign `MomentProblem`.

### Step 2 — prototype an H4-specific momentum tagger in the benchmark layer first

Do **not** generalize everything blindly in one shot.

First prove the idea on the H4 benchmark side with a helper that tags monomials by whatever survives normal ordering and is physically meaningful for this asset, e.g. some combination of:

- `ΔN`,
- total crystal momentum `K`,
- maybe `S_z` if it survives cleanly through the operator encoding.

Start in demo/benchmark code, not package core. Earn the abstraction.

### Step 3 — add a structure-only benchmark before attempting a solve

Before asking a solver to do anything, measure whether the partition actually reduces:

- max block size,
- `Σ block_dim²`,
- PSD upper-tri slots,
- symbolic moment count,
- JuMP model size.

That benchmark is cheap and decisive.

### Step 4 — only then rerun primal COSMO / Clarabel

If the structural counts do not drop hard, stop. Don't waste more solver time.

If they *do* drop hard, rerun:

- primal `dualize=false` + COSMO,
- maybe Clarabel,
- maybe Mosek if available.

## Concrete implementation slices

Here is the order I would use.

### Slice 1 — structural API hook

Deliverable:

- optional basis-partition hook in `SolverConfig`
- tests on tiny fermionic toy cases showing partition purity and unchanged correctness on trivial cases

Likely files:

- `src/optimization/interface.jl`
- `src/optimization/sparsity.jl`
- `test/relaxations/interface.jl`
- likely a small new targeted test in `test/correlated_sparsity/` or `test/polynomials/`

### Slice 2 — H4 benchmark-side tagger + structure report

Deliverable:

- benchmark helper that tags H4 basis words by symmetry labels
- a new `demos/results/*.md` report comparing:
  - current TS blocks,
  - symmetry-partitioned TS blocks,
  - maybe raw sector sizes before TS

Likely files:

- `demos/` or `demos/results/` only
- maybe `test/H4PeriodicAssets.jl` for tiny shared helpers if they are genuinely reusable

### Slice 3 — rerun primal solve on reduced formulation

Deliverable:

- updated logs showing whether solver setup becomes practical
- honest stop/go decision

## Recommended acceptance criteria

A reduction attempt is worth keeping only if it materially moves at least one of these numbers on H4 Nk=2:

- max PSD block dimension,
- `Σ block_dim²`,
- PSD upper-tri scalar slots,
- direct assembled conic rows,
- or time to reach solver iteration 1.

If the change does not hit those, it is theater.

## Practical recommendation

**Do not chase more solver knobs on the current formulation. Implement the symmetry-partition hook first.**

That is the smallest change with a real chance of changing the problem class instead of just rearranging deck chairs.

If you want the blunt version:

- naive pair basis alone: not enough,
- typed expression cleanup: unrelated side quest,
- another solver probe on the same formulation: mostly wasted time,
- symmetry-aware reduction before TS: finally worth doing.

## Bottom line

The repo already contains the right architectural seams:

- explicit `moment_basis`,
- blockwise moment assembly,
- ordinary polynomial `moment_eq_constraints`,
- public primal solve path.

What it lacks is the one thing H4 actually needs:

- **a way to encode translation/symmetry sectors before the generic term-sparsity machinery shatters them.**

That is the next serious move.
