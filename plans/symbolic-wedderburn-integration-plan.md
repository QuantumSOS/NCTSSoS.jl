# SymbolicWedderburn Integration Plan for Symmetry Reduction

This plan replaces the current homemade representation-decomposition logic in `src/optimization/symmetry.jl` with `SymbolicWedderburn.jl`, while keeping the existing MVP symmetry API and the current architectural seam intact.

## Goal

Use `SymbolicWedderburn` for the part it is actually good at:

- symmetry-adapted block decomposition of the PSD constraint
- explicit change-of-basis matrices
- irrep / direct-summand extraction

Keep all noncommutative and quotient-sensitive logic local:

- word action on `NormalMonomial`
- quotient canonicalization with `simplify`
- strict invariance checks
- strict closure checks
- moment-variable orbit reduction into an ordinary reduced `MomentProblem`

In other words:

- keep the current public symmetry API
- keep the current "ordinary reduced `MomentProblem`" output contract
- swap only the internal representation-theory engine

## Scope

### In scope

- dense only
  - `cs_algo = NoElimination()`
  - `ts_algo = NoElimination()`
- one dense clique only
- ordinary polynomial optimization only
- monoid algebras only
- real signed-permutation actions only (`±1` coefficients)
- multiplicity-free symmetry reductions only
- CHSH order-1 dense unipotent symmetry as the acceptance target

### Out of scope

- multiple clique composition
- term sparsity composition
- projector affine actions
- non-unitary actions / correction matrices `Cᵢ`
- state / trace support
- PBW / Pauli / fermionic / bosonic symmetry support
- full replacement of monomial-based orbit reduction with invariant-vector coordinates

## Non-negotiable design choice

Do **not** rewrite the repo around symmetry-adapted linear-combination basis elements.

The existing package still wants to solve an ordinary `MomentProblem` whose symbolic entries are expressed in ordinary polynomial / monomial coordinates.

So the right split is:

- local code handles NC quotient semantics and reduced monomial variables
- `SymbolicWedderburn` handles block decomposition on the half-basis `Q`

That keeps the change small, reversible, and trustworthy.

## Public API: keep it unchanged

Keep these types and fields exactly as the symmetry-facing API:

- `SignedPermutation`
- `SymmetrySpec`
- `SymmetryReport`
- `SolverConfig.symmetry`
- `PolyOptResult.symmetry`

No user-facing redesign in this step.

## Dependency changes

### `Project.toml`

Add:

- `SymbolicWedderburn`

Only add extra group packages if needed by the actual integration path. Do not guess and spray dependencies around.

## Files to change

### 1. `src/optimization/symmetry.jl`

This file remains the symmetry integration layer.

#### Keep local

These functions / responsibilities stay in this file:

- finite group enumeration from `SymmetrySpec.generators`
- `SignedPermutation` validation
- signed action on registry indices
- signed action on `NormalMonomial`
- signed action on `Polynomial`
- invariance checks for:
  - objective
  - equality constraints
  - inequality constraints
  - moment equality constraints
- closure checks on basis blocks
- orbit reduction for full moment monomials
- emission of a reduced ordinary `MomentProblem`

#### Replace

Delete the current homemade representation-decomposition machinery and replace it with `SymbolicWedderburn` calls.

Candidate functions to remove or rewrite:

- `_IrrepSubspace`
- `_representation_matrices`
- `_commutant_basis`
- `_symmetric_commutant_basis`
- `_eigenvalue_groups`
- `_is_scalar_matrix`
- `_find_splitter`
- `_decompose_irreducible_subspaces`
- `_characters_match`
- `_symmetry_block_decomposition`

Those are exactly the "poor man’s Wedderburn" bits.

### 2. `src/optimization/interface.jl`

Keep the current hook point:

- when `solver_config.symmetry === nothing`: ordinary path
- when `solver_config.symmetry !== nothing`: call `moment_relax_symmetric(...)`

Keep the MVP support guards here:

- dense only
- one clique only
- no state/trace
- no term sparsity composition

### 3. `src/NCTSSoS.jl`

No architecture change.

Just ensure the dependency wiring and exports remain correct.

### 4. `test/problems/bell_inequalities/chsh_simple.jl`

Keep the current symmetry-enabled acceptance test.

It remains the main regression target.

### 5. `test/data/expectations/chsh_simple.toml`

Keep the `Symmetry_d1` fixture.

That is the acceptance oracle for this migration.

## Internal adapter design

Add a thin adapter from our word action to `SymbolicWedderburn`.

### New internal action type

Something like:

```julia
struct NCWordSignedPermutationAction{A,T} <: SymbolicWedderburn.BySignedPermutations
    algebra::Type{A}
end
```

Then define:

```julia
SymbolicWedderburn.action(act::NCWordSignedPermutationAction{A,T}, g, mono::NormalMonomial{A,T})
```

Behavior:

1. apply the signed letter map to `mono.word`
2. simplify using existing `simplify(A, word)`
3. return the image monomial together with the sign in the format expected by `SymbolicWedderburn`

This adapter must stay monomial-based. Do not generalize it to arbitrary linear combinations in this step.

## Representation-theory split

### Use `SymbolicWedderburn` for

- decomposition of the representation on the half-basis `Q`
- extraction of symmetry-adapted block bases
- block diagonalization of the PSD matrix

### Keep local for now

- orbit reduction of full moment monomials
- reduced symbolic objective / constraints expressed in ordinary monomial coordinates

This is the key engineering choice.

Do not try to make `total_basis` or the full reduced SDP variable layer depend directly on invariant vectors or abstract linear-combination basis elements yet.

That is a larger refactor and not needed for this slice.

## Integration inside `moment_relax_symmetric(...)`

The new flow should be:

1. **Enumerate the finite group** from `SymmetrySpec.generators`.
2. **Check strict invariance** of objective and constraints.
3. **Check strict closure** of each dense basis block under the action.
4. **Build the ordinary symbolic data** exactly as today.
5. **Construct a `SymbolicWedderburn` decomposition** on the half-basis `Q` for the PSD block.
6. **Use the decomposition to block-diagonalize the symbolic PSD matrix**.
7. **Apply the existing local orbit reducer** to the transformed polynomial entries.
8. **Verify off-block cancellation and scalar-block structure** under the MVP assumptions.
9. **Emit a smaller ordinary `MomentProblem`** consisting of scalar PSD blocks and reduced moment variables.
10. **Return a `SymmetryReport`** with the expected metadata.

## Recommended exact internal structure

### Step A — build the action adapter

Implement and test the adapter independently:

- action on `Q = [1, x₁, x₂, y₁, y₂]`
- closure holds
- images are correct signed monomials

### Step B — build Wedderburn decomposition on `Q`

Use either:

- `WedderburnDecomposition(...)`
- or the smaller `symmetry_adapted_basis(...)` path if that is sufficient for the half-basis-only decomposition

Choose the smallest API that gives us:

- direct summands
- change-of-basis matrices
- stable block extraction

Do not over-integrate package features we do not need.

### Step C — transform the symbolic PSD matrix

Given the ordinary symbolic moment matrix on `Q`, apply the package-provided basis change to obtain block-diagonal symbolic matrices.

The result should still be ordinary polynomial-valued matrices.

### Step D — reduce polynomial entries by moment orbits

Use the current local orbit reducer on the transformed symbolic entries.

This preserves the existing reduced-problem semantics:

- objective and constraints remain ordinary symbolic polynomials
- moment variables are still represented by reduced monomial orbit reps

### Step E — enforce MVP block shape

After reduction, require:

- off-block entries vanish exactly
- each surviving symmetry block is scalar
- multiplicity-free only

If not, error.

No silent fallback.

## Multiplicity policy

Even if `SymbolicWedderburn` can represent larger isotypic blocks, this migration should keep the current MVP strictness.

After decomposition, explicitly require:

- multiplicity `1` for every surviving direct summand relevant to the PSD block
- resulting reduced blocks match the scalar-block expectation of the current emission path

If repeated irreps occur, throw an error explaining that multiplicity-bearing reductions are not yet supported by the current reduced `MomentProblem` builder.

## CHSH acceptance target

Keep the exact first target:

- unipotent CHSH
- dense order-1 relaxation
- explicit

```julia
Q = [1, x₁, x₂, y₁, y₂]
```

- symmetry generators:
  - Alice setting swap + sign on `y₂`
  - Bob setting swap + sign on `x₂`
  - party swap

Expected assertions:

- `result.objective ≈ -2sqrt(2)`
- `result.symmetry.group_order == 16`
- `result.symmetry.invariant_moment_count == 1`
- `result.symmetry.psd_block_sizes == [1, 1, 1]`

The migration is not done until these still pass.

## Additional regression tests to add after the swap

### Low-level symmetry tests

Add narrow tests for:

1. **Invariance failure**
   - pass a symmetry that does not preserve the objective
   - expect a hard error

2. **Closure failure**
   - pass a basis `Q` not closed under the action
   - expect a hard error

3. **Multiplicity rejection**
   - if a cheap deterministic toy case exists, ensure repeated irreps currently error

These should be separate from the CHSH end-to-end test.

## Suggested implementation phases

### Phase 1 — dependency + adapter

- add `SymbolicWedderburn`
- add `NCWordSignedPermutationAction`
- add focused tests for action correctness on monomials / `Q`

### Phase 2 — swap decomposition internals

- remove homemade commutant / splitter logic
- replace with `SymbolicWedderburn` block decomposition
- keep local orbit reducer
- keep CHSH green

### Phase 3 — add strict regression tests

- invariance failure
- closure failure
- multiplicity rejection

### Phase 4 — optional docs follow-up

- add a small note in symmetry-facing docs once the new backend is stable

## What this migration intentionally does not do

This migration does **not** attempt to:

- support multiple cliques
- support term sparsity + symmetry together
- support projector affine symmetry actions
- support non-unitary `Cᵢ` correction paths
- support state / trace symmetry reduction
- replace the monomial-based reduced variable layer with invariant-vector coordinates

Those are future steps and should stay future steps.

## Definition of done

This migration is done when:

1. `SymbolicWedderburn` is the decomposition backend for dense MVP symmetry reduction.
2. The current public symmetry API is unchanged.
3. CHSH dense symmetry still passes with:
   - group order `16`
   - invariant moment count `1`
   - PSD block sizes `[1, 1, 1]`
   - objective `-2sqrt(2)`
4. Invariance and closure remain strict fail-fast checks.
5. Unsupported multiplicity / non-unitary / non-dense cases still error clearly instead of silently doing nonsense.

## Final engineering stance

The correct migration is not "replace the whole symmetry layer with SymbolicWedderburn".

The correct migration is:

- keep all NC quotient semantics local
- keep reduced ordinary `MomentProblem` output local
- use `SymbolicWedderburn` only for block decomposition and change of basis
- keep the slice dense, strict, and multiplicity-free

That is the smallest sound change.
