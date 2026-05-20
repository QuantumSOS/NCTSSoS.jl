# Symmetry-Adapted Basis Implementation Plan for `NCTSSoS.jl`

I read the Typst note and the repo pieces that matter.

## What’s already true in this repo

For the existing dense order-1 CHSH test, the package is already sitting exactly at the pre-symmetry point from the note:

- `Q = [1, x₁, x₂, y₁, y₂]`
- one dense `5×5` moment PSD block
- `mp.total_basis` has **13** words
- current solver uses **11** unique moment variables after `symmetric_canon`

That matches the note’s “big before symmetry” picture closely enough to be a solid implementation target.

Also: current **term sparsity is not symmetry reduction**.
The repo already gets CHSH to `[3,3,1]` blocks with TS, but that’s graph chopping, not representation theory. Different animal.

---

## My take

Don’t start by replacing the repo’s basis layer with linear-combination basis vectors.

That would be the fancy wrong move.

Why? Because the current pipeline is still quite hard-wired to:

- `NormalMonomial`
- `NCStateWord`
- clique-local ordinary word bases in `sparsity.jl`
- `MomentProblem{..., M<:NormalMonomial, ...}` in `moment.jl`

If you try to make the whole engine “symmetry-adapted basis vectors first”, you’ll touch half the package before you have one passing test.

**The minimal sane plan is:**
keep the existing word basis and moment construction,
then insert a **symmetry reduction pass that returns an ordinary reduced `MomentProblem`**.

That way:

- `compute_sparsity` mostly stays intact
- `solve_sdp`, `sos_dualize`, direct moment solving stay intact
- the new code lives in one place
- CHSH becomes a clean first acceptance test

---

# Recommended implementation plan

## Phase 1: MVP for dense monoid problems, with CHSH as the acceptance test

### Scope

Support only:

- `MonoidAlgebra` first, especially `UnipotentAlgebra`
- dense / no sparsity first:
  - `cs_algo = NoElimination()`
  - `ts_algo = NoElimination()`
- unitary actions only:
  - permutation or **signed permutation** actions on letters/words
- ordinary polynomial problems first
- no non-unitary correction matrices yet
- no state/trace path yet
- no PBW algebras yet

That is not cowardice. That is how you land the feature instead of writing a thesis in `git diff`.

---

## Phase 2: Add the symmetry API

### Files

- `src/optimization/interface.jl`
- new `src/optimization/symmetry.jl`
- `src/NCTSSoS.jl`

### Add to `SolverConfig`

Add something like:

```julia
symmetry::Union{Nothing,SymmetrySpec} = nothing
```

Where `SymmetrySpec` is explicit, not a bag of kwargs.

For MVP it can hold:

- `group`
- `action`
- maybe `check_invariance::Bool = true`

No need to design a cathedral.

---

## Phase 3: Implement word-level symmetry actions

### New module

`src/optimization/symmetry.jl`

### Core responsibilities

1. **Act on canonical words**
   - map encoded registry indices to new encoded indices
   - preserve order
   - re-canonicalize with existing `simplify!`

2. **Support signed permutations**
   For unipotent CHSH this is enough, and it avoids the affine/projector mess.

3. **Orbit normal form**
   Build a symmetry normal form for moment monomials:
   - quotient-reduce first
   - then orbit-min under the group
   - then still respect the repo’s existing moment convention (`symmetric_canon`)

4. **Check invariance**
   Fail fast unless:
   - objective is invariant
   - all active constraints are invariant
   - all moment equalities are invariant too, if present

5. **Check closure**
   For each block basis `Q`, verify `g(Qᵢ)` stays in `span(Q)`.

For the first slice, closure can be checked on the explicit CHSH `Q`.

---

## Phase 4: Use `SymbolicWedderburn`, don’t reinvent it

### New deps in `Project.toml`

You almost certainly want:

- `SymbolicWedderburn`
- `PermutationGroups`
- `AbstractPermutations`
- maybe `GroupsCore` as indirect plumbing

The note is right here: the representation theory backend already exists. Use it.

### Basis inputs

For CHSH order 1:

- `basis_half = Q = [1, x₁, x₂, y₁, y₂]`
- `basis_full =` canonical words appearing in `Q Q†`
  - here that’s the 13-word dense moment basis the repo already computes

Then:

- `WedderburnDecomposition(..., basis_full, basis_half)`
- use `invariants` to collapse moment variables
- use `Uπs` / block basis to block-diagonalize the PSD constraint

---

## Phase 5: Integrate symmetry **between sparsity and solving**

### Best hook

Not “replace basis generation everywhere”.

Instead:

- keep `compute_sparsity(...)`
- when `symmetry === nothing`, use current path
- when `symmetry !== nothing`, route to a new builder:

```julia
moment_problem, symmetry_report = moment_relax_symmetric(pop, sparsity, solver_config.symmetry)
```

This new path should:

1. build the ordinary symbolic constraint block(s)
2. apply symmetry reduction
3. emit a smaller, plain `MomentProblem`

Then `solve_sdp` stays untouched.

That’s the right seam.

### Why not reduce after `solve_sdp`?

Because then you’ve already paid the big SDP cost. Useless.

### Why not reduce before `moment_relax` with new basis types?

Because then you’re rewriting `sparsity.jl` and probably half the type signatures. Also useless for MVP.

---

## Phase 6: Result reporting

Current `PolyOptResult` only knows about sparsity block sizes, not symmetry blocks.

Add a small symmetry report, something like:

```julia
struct SymmetryReport
    group_order::Int
    invariant_moment_count::Int
    psd_block_sizes::Vector{Int}
    basis_half_size::Int
    basis_full_size::Int
end
```

Then extend `PolyOptResult` with:

```julia
symmetry::Union{Nothing,SymmetryReport}
```

Keep old fields for backward compatibility. Don’t overload “term sparsity block sizes” with symmetry meanings.

---

# CHSH acceptance test

## Where

Extend:

- `test/problems/bell_inequalities/chsh_simple.jl`

Don’t create a parallel near-copy unless you enjoy test sprawl.

## Formulation

Use the existing unipotent CHSH first, not the projector/affine version.

Reason:
- it is already in the repo
- it already gives the `5×5` dense block
- its symmetry can be expressed as **signed permutations**
- that avoids affine actions like `a_x ↦ 1 - a_x` in slice 1

## Make `Q` explicit

Use:

```julia
Q = [one(x[1]), x[1], x[2], y[1], y[2]]
```

via `moment_basis=Q` in the solver config.

That makes closure and debugging much cleaner.

## Candidate symmetry generators

For the unipotent CHSH operator
`x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂`,
a good signed-permutation generating set is:

- Alice setting swap + sign flip on `y₂`
- Bob setting swap + sign flip on `x₂`
- party swap `xᵢ ↔ yᵢ`

That should generate the standard size-16 CHSH symmetry on `span(Q)`.

## Assertions

The symmetry test should check at least:

- `result.objective ≈ -2sqrt(2)`
- `result.symmetry.group_order == 16`
- `result.symmetry.invariant_moment_count == 1`
- `result.symmetry.psd_block_sizes == [1, 1, 1]`

That’s the real success criterion. Anything weaker is just “we moved some matrices around”.

You can add a reviewed expectation entry to:
- `test/data/expectations/chsh_simple.toml`
with a new case id like `Symmetry_d1`.

---

# What to explicitly defer

Do **not** put these in the first implementation slice:

1. **Projector CHSH with affine actions**
   - useful follow-up
   - not MVP

2. **Non-unitary linear actions / correction matrices `Cᵢ`**
   - real feature
   - not needed for unipotent CHSH

3. **Composition with term sparsity**
   - current TS code assumes ordinary monomial bases
   - this needs its own refactor

4. **PBW / Pauli / fermionic symmetry**
   - separate complexity class

---

# Main correctness traps

1. **Using a group that does not preserve the objective/constraints**
   - reduction becomes wrong, quietly
   - must error, not warn

2. **Skipping closure checks**
   - then `ρ_g` is fake

3. **Confusing quotient canonicalization with orbit canonicalization**
   - they are different steps
   - both are needed

4. **Trying to expose symmetry-adapted basis vectors as first-class basis elements too early**
   - that’s the swamp

---

# Short PR sequence I’d use

## PR 1

- add deps
- add `SymmetrySpec`
- add word action + invariance + closure checks
- add low-level tests

## PR 2

- add `moment_relax_symmetric(...)`
- dense monoid/unipotent path only
- CHSH end-to-end symmetry test

## PR 3

- add docs snippet to `docs/src/examples/literate/bell.jl`
- update API docs for `SolverConfig`

## PR 4+

- per-clique correlative sparsity composition
- then term sparsity
- then projector/affine actions
- then non-unitary path

---

# Additional implementation notes preserved from the planning discussion

## Why this seam is the right one

The existing codebase already cleanly supports:

- symbolic moment relaxation construction
- canonical word rewriting
- user-supplied `moment_basis`
- plain dense and sparse SDP emission through `MomentProblem`

The missing layer is specifically the representation-theoretic reduction step.
That means the cleanest design is to add a symmetry-aware transformation layer that consumes the already-constructed symbolic data and emits a smaller ordinary symbolic problem.

In other words, symmetry should be a transformation on the relaxation, not a rewrite of the entire package’s basis abstraction.

## Why the repo is a good fit for this plan

The repo already provides the hard NC algebra plumbing:

- `simplify!` for quotient reduction
- `NormalMonomial` for canonical words
- `get_ncbasis` for explicit basis enumeration
- `moment_relax` for symbolic matrix assembly
- `SolverConfig(moment_basis=...)` for explicit `Q`

So the new feature should lean on those instead of replacing them.

## Why CHSH is the right first test

CHSH at order 1 is the exact sweet spot:

- small enough to inspect by hand
- already present in tests
- already produces the expected dense `5×5` block
- known symmetry group size: `16`
- known reduced answer: one invariant moment variable and three scalar positivity conditions
- known optimum: Tsirelson bound `2√2`

That makes it a real acceptance test instead of a vague “seems smaller” demo.

## What not to do first

Avoid these bad first moves:

- introducing a completely new basis element type for arbitrary linear combinations and threading it through `sparsity.jl`, `moment.jl`, `sos.jl`, and the test suite before proving the feature on CHSH
- coupling symmetry and term sparsity in the same first PR
- trying to support non-unitary actions before landing a correct unitary/signed-permutation path
- starting from projector CHSH if it forces affine rather than linear actions on the chosen basis

That’s how you get a giant diff and no trustworthy result.

## Minimal new data structures worth adding

The smallest useful additions are:

- `SymmetrySpec`
- `SymmetryReport`
- internal helpers for action/invariance/closure/orbit normal form
- a new `moment_relax_symmetric(...)` path

Anything beyond that should earn its keep.

## Expected first-slice behavior

For the initial implementation slice, symmetry support should probably be deliberately strict:

- error if the action is not supported
- error if closure fails
- error if invariance fails
- error if the reduction would require a not-yet-supported non-unitary path
- error if sparsity settings are active and the combination is not yet implemented

Strict is better than quietly wrong.

## Intended composition story, later

Once the dense path is correct, the intended long-term composition is:

1. build correlative sparsity cliques as today
2. for each clique-local moment block, apply symmetry reduction inside the clique
3. later, optionally combine that with term sparsity splitting inside each symmetry block

But that is a later step. First land the dense case and prove CHSH.

## Reporting expectations

When symmetry is enabled, callers should be able to inspect:

- the original half-basis size
- the original full moment-word basis size
- the number of invariant moments after orbit collapse
- the final reduced PSD block sizes
- the acting group order

That data belongs in a dedicated symmetry report, not jammed awkwardly into existing sparsity fields.

## Test philosophy for this feature

A real regression test for symmetry reduction should verify both:

- **correctness of the bound**
- **correctness of the reduced structure**

For CHSH, that means not merely checking `-2sqrt(2)`, but also checking the structural reduction:

- group order `16`
- invariant moment count `1`
- scalar PSD blocks `[1, 1, 1]`

Otherwise you can accidentally pass the bound while not actually doing symmetry reduction in the intended way.

## Practical first milestone

The first milestone should be modest and explicit:

> Solve dense order-1 unipotent CHSH with symmetry enabled, obtain the Tsirelson bound, and expose symmetry reduction metadata showing one invariant moment variable and three scalar PSD blocks.

That is enough to prove the design works.

## Practical second milestone

Only after the first milestone is stable:

- allow symmetry with clique-based correlative sparsity
- then revisit term sparsity composition
- then add projector/affine action support
- then tackle non-unitary linear actions and `Cᵢ`

That order matters.

## Final engineering stance

The simple sound solution is:

- keep the current monomial-based engine
- add a symmetry reduction layer at the relaxation seam
- make CHSH prove the design
- defer the clever parts until the dumb parts are airtight

That’s the plan worth implementing.