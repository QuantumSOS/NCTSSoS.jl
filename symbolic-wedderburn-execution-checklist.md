# SymbolicWedderburn Integration — Tight Execution Checklist

This is the implementation checklist version of `symbolic-wedderburn-integration-plan.md`.

The goal is simple:

- keep the current public symmetry API
- keep the current reduced ordinary `MomentProblem` output contract
- replace only the homemade representation-decomposition machinery with `SymbolicWedderburn`
- keep CHSH dense order-1 symmetry green the whole time

---

## Non-negotiable acceptance gate

Do not call this done unless the symmetry-enabled CHSH test still gives:

- `result.objective ≈ -2sqrt(2)`
- `result.symmetry.group_order == 16`
- `result.symmetry.invariant_moment_count == 1`
- `result.symmetry.psd_block_sizes == [1, 1, 1]`

And do not widen scope beyond:

- dense only
- one clique only
- ordinary polynomial problems only
- monoid algebras only
- real signed-permutation actions only
- multiplicity-free only

If anything tries to pull you beyond that, stop and error out.

---

## 0. Baseline capture

### Checklist

- [ ] Confirm the current repo still loads.
- [ ] Confirm the current CHSH symmetry regression passes before touching the backend.
- [ ] Record the current symmetry result values as the migration baseline.

### Commands

```bash
julia --project -e 'using NCTSSoS; println("ok")'
```

```bash
julia --project -e '
using NCTSSoS, Test
include("test/Expectations.jl")
using .TestExpectations: expectations_oracle
include("test/TestUtils.jl")
include("test/problems/bell_inequalities/chsh_simple.jl")
'
```

### Expected outcome

The `Dense + Symmetry` CHSH test passes with the known metadata.

---

## 1. Add the dependency, nothing else

### Files

- [ ] `Project.toml`
- [ ] `Manifest.toml` if the repo tracks it and the workflow expects it updated

### Checklist

- [ ] Add `SymbolicWedderburn` to `[deps]`.
- [ ] Only add extra group packages if `SymbolicWedderburn` actually requires them for the chosen integration path.
- [ ] Verify the package resolves and precompiles.

### Stop condition

If adding the dependency causes unrelated breakage, stop and isolate that first. Don’t pile changes on top of a broken environment.

---

## 2. Freeze the public API

### Files

- [ ] `src/optimization/interface.jl`
- [ ] `src/NCTSSoS.jl`

### Checklist

- [ ] Keep these user-facing names unchanged:
  - `SignedPermutation`
  - `SymmetrySpec`
  - `SymmetryReport`
  - `SolverConfig.symmetry`
  - `PolyOptResult.symmetry`
- [ ] Do **not** redesign the API while swapping the backend.
- [ ] Keep the current MVP support guards in place.

### Explicit non-goal

No new user-facing symmetry abstraction in this step.

---

## 3. Split the code into “keep local” vs “replace”

### File

- [ ] `src/optimization/symmetry.jl`

### Keep local

These stay owned by NCTSSoS:

- [ ] `SignedPermutation` validation and group enumeration
- [ ] action on registry indices
- [ ] action on `NormalMonomial`
- [ ] action on `Polynomial`
- [ ] strict invariance checks
- [ ] strict closure checks
- [ ] orbit reduction for moment monomials
- [ ] reduced ordinary `MomentProblem` emission

### Replace

These are the homemade representation-theory bits to remove or rewrite around `SymbolicWedderburn`:

- [ ] `_IrrepSubspace`
- [ ] `_representation_matrices`
- [ ] `_commutant_basis`
- [ ] `_symmetric_commutant_basis`
- [ ] `_eigenvalue_groups`
- [ ] `_is_scalar_matrix`
- [ ] `_find_splitter`
- [ ] `_decompose_irreducible_subspaces`
- [ ] `_characters_match`
- [ ] `_symmetry_block_decomposition`

### Rule

Do not delete the old machinery until the adapter path is in place and the file still loads.

---

## 4. Add the thinnest possible SymbolicWedderburn adapter

### File

- [ ] `src/optimization/symmetry.jl`

### Add

- [ ] `using SymbolicWedderburn`
- [ ] a thin internal adapter type, e.g.

```julia
struct NCWordSignedPermutationAction{A,T} <: SymbolicWedderburn.BySignedPermutations
    algebra::Type{A}
end
```

### Implement

- [ ] `SymbolicWedderburn.action(adapter, g, mono::NormalMonomial{A,T})`
- [ ] If needed, helper conversion from our enumerated group elements to the group element representation expected by the package.

### Adapter requirements

The adapter must:

- [ ] apply the signed letter map to a monomial word
- [ ] simplify using `simplify(A, word)`
- [ ] return the image plus sign in the package’s expected format
- [ ] remain monomial-based only

### Do not do

- [ ] Do not support arbitrary linear combinations here.
- [ ] Do not make the adapter understand full reduced moment-variable semantics.

---

## 5. Prove the adapter on CHSH before touching the solver path

### File

- [ ] temporary scratch code or a narrow internal test

### Checklist

On `Q = [1, x₁, x₂, y₁, y₂]`:

- [ ] verify every generator acts correctly
- [ ] verify the closure check still passes
- [ ] verify the generated group order is `16`
- [ ] verify the action is real signed-permutation / orthogonal as expected

### Required result

You should be able to say, with evidence, that the adapter represents the same CHSH symmetry we already use.

If not, stop. The rest is meaningless if the action is wrong.

---

## 6. Introduce a narrow Wedderburn wrapper helper

### File

- [ ] `src/optimization/symmetry.jl`

### Add helper(s)

Something like:

- [ ] `_sw_decompose_half_basis(...)`
- [ ] `_sw_block_change_of_basis(...)`
- [ ] `_sw_direct_summands(...)`

The exact names do not matter. The split does.

### Helper contract

Input:

- [ ] finite group
- [ ] signed-permutation adapter
- [ ] half-basis `Q`

Output:

- [ ] change-of-basis data for the PSD block
- [ ] direct-summand metadata
- [ ] enough information to transform the symbolic PSD matrix blockwise

### Rule

Use only the `SymbolicWedderburn` surface needed for the dense CHSH slice.
Do not integrate every API the package exposes just because it exists.

---

## 7. Enforce multiplicity-free immediately after decomposition

### File

- [ ] `src/optimization/symmetry.jl`

### Checklist

After obtaining the direct summands:

- [ ] inspect multiplicities
- [ ] inspect resulting block dimensions
- [ ] fail hard unless the decomposition is compatible with the current MVP emission path

### MVP policy

Require:

- [ ] multiplicity `1` for each relevant direct summand
- [ ] scalar reduced PSD blocks under the current CHSH-style path

### Error behavior

If repeated irreps show up, error with a message that says the current MVP only supports multiplicity-free symmetry reductions.

No silent fallback.

---

## 8. Swap `moment_relax_symmetric(...)` to use SymbolicWedderburn for block decomposition only

### File

- [ ] `src/optimization/symmetry.jl`

### Replace inside `moment_relax_symmetric(...)`

Old internal logic:

- homemade representation decomposition
- homemade basis-change extraction
- homemade scalar-block decomposition

New internal logic:

1. [ ] enumerate group
2. [ ] check invariance
3. [ ] check closure
4. [ ] build ordinary symbolic moment matrix on `Q`
5. [ ] call SymbolicWedderburn on `Q`
6. [ ] extract the symmetry-adapted block basis / change of basis
7. [ ] transform the symbolic PSD matrix
8. [ ] keep local orbit reduction on the transformed polynomial entries
9. [ ] verify off-block cancellation and scalar-block structure
10. [ ] emit reduced ordinary `MomentProblem`

### Important rule

Keep the full moment-variable reduction local.

Do **not** replace `total_basis` semantics with invariant vectors in this step.

That is a different refactor.

---

## 9. Keep local orbit reduction exactly where it buys us safety

### File

- [ ] `src/optimization/symmetry.jl`

### Keep alive

- [ ] `_build_orbit_reducer`
- [ ] `_orbit_reduce_polynomial`
- [ ] `_orbit_reduce_matrix`

### Rationale

This preserves:

- [ ] ordinary monomial-based reduced symbolic constraints
- [ ] ordinary reduced `MomentProblem`
- [ ] compatibility with the current solver path

### Rule

Use `SymbolicWedderburn` for block decomposition.
Use the local orbit reducer for reduced monomial variables.

That split is the whole point.

---

## 10. Delete the homemade decomposition code only after CHSH is green

### File

- [ ] `src/optimization/symmetry.jl`

### Checklist

- [ ] make the new path compile
- [ ] make the CHSH symmetry test pass
- [ ] then remove the dead homemade decomposition functions
- [ ] run formatting / sanity checks if the repo uses them

### Rule

No heroic “big bang” rewrite. Keep the file loadable at every intermediate step.

---

## 11. Re-run the CHSH acceptance test immediately

### File

- [ ] `test/problems/bell_inequalities/chsh_simple.jl`
- [ ] `test/data/expectations/chsh_simple.toml`

### Checklist

- [ ] run the CHSH file again
- [ ] confirm the symmetry case still passes
- [ ] confirm metadata is unchanged:
  - `group_order == 16`
  - `invariant_moment_count == 1`
  - `psd_block_sizes == [1, 1, 1]`
- [ ] confirm the objective is still `≈ -2sqrt(2)`

### Command

```bash
julia --project -e '
using NCTSSoS, Test
include("test/Expectations.jl")
using .TestExpectations: expectations_oracle
include("test/TestUtils.jl")
include("test/problems/bell_inequalities/chsh_simple.jl")
'
```

---

## 12. Add low-level regression tests after the backend swap

### Recommended new tests

Add focused tests for:

- [ ] invariance failure
- [ ] closure failure
- [ ] multiplicity rejection, if a cheap deterministic example exists

### Placement

Prefer a narrow symmetry-focused test file rather than bloating the CHSH problem test.

If no dedicated symmetry test file exists yet, create one under the relevant optimization / relaxation area.

### Rule

These tests should validate fail-fast behavior, not solver numerics.

---

## 13. Final sanity pass

### Checklist

- [ ] `using NCTSSoS` works cleanly
- [ ] CHSH symmetry regression passes
- [ ] no public API drift
- [ ] unsupported cases still error clearly
- [ ] no accidental broadening of supported scope

### Git diff smell test

The final diff should show:

- smaller / cleaner decomposition code in `src/optimization/symmetry.jl`
- no basis-system rewrite
- no solver rewrite
- no term-sparsity entanglement

If the diff starts looking like a new architecture, you took a wrong turn.

---

# Exact file touch list

## Must touch

- [ ] `Project.toml`
- [ ] `src/optimization/symmetry.jl`

## Probably touch

- [ ] `Manifest.toml`
- [ ] `test/problems/bell_inequalities/chsh_simple.jl`
- [ ] `test/data/expectations/chsh_simple.toml`

## Should not need major changes

- [ ] `src/optimization/interface.jl`
- [ ] `src/NCTSSoS.jl`

If those files need large edits, stop and check whether you are leaking the backend swap into the public API.

---

# Suggested commit order

## Commit 1

**chore(symmetry): add SymbolicWedderburn dependency**

- add dependency
- verify package loads

## Commit 2

**refactor(symmetry): add SymbolicWedderburn signed-word adapter**

- add adapter type
- add helper bridge code
- no solver-path change yet

## Commit 3

**refactor(symmetry): replace homemade block decomposition backend**

- wire `moment_relax_symmetric(...)` to `SymbolicWedderburn`
- keep orbit reducer local
- keep CHSH green

## Commit 4

**test(symmetry): add fail-fast regression tests**

- invariance failure
- closure failure
- multiplicity rejection if available

That order is sane. Anything more clever is just an invitation to waste time.

---

# Red flags — stop if you hit one

- [ ] You need to redesign `MomentProblem` to proceed.
- [ ] You need to thread linear-combination basis elements through `sparsity.jl`.
- [ ] You need to support non-unitary actions just to keep CHSH working.
- [ ] You need multiple-clique support for the first pass.
- [ ] You cannot explain whether a failure is in the NC action layer or the decomposition layer.

If any of these happen, the seam got broken.

---

# Definition of done

The migration is done when all of these are true:

- [ ] `SymbolicWedderburn` is the decomposition backend for dense MVP symmetry reduction
- [ ] the current public symmetry API is unchanged
- [ ] CHSH dense order-1 symmetry still passes
- [ ] metadata stays:
  - [ ] `group_order = 16`
  - [ ] `invariant_moment_count = 1`
  - [ ] `psd_block_sizes = [1, 1, 1]`
- [ ] strict fail-fast behavior remains for unsupported cases
- [ ] homemade decomposition code is gone
- [ ] the reduced problem is still emitted as an ordinary `MomentProblem`

That’s the target. Small, sharp, correct.
