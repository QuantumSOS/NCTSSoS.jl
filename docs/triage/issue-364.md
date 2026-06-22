# Triage: Issue #364 — Sparse chain basis + degree 3-4 charge/singlet + sign symmetry

**Issue**: https://github.com/QuantumSOS/NCTSSoS.jl/issues/364  
**Labels**: enhancement, SOTA-target  
**Author**: exAClior  
**Sub-issue of**: #356

---

## Issue Summary

Add four coordinated features that together enable the N=100 periodic Heisenberg XXX ground-state SDP at degree 4, matching block sizes from Wang et al. (arXiv:2604.01555):

1. **Sparse contiguous-support basis** — `pauli_contiguous_chain_basis(registry, degree; periodic=true)` producing O(N) elements instead of O(N²).
2. **Extend `PauliChargeSectorSpec` to degree 3-4** — lift the `max_degree ∈ 0:2` cap, with basis-driven charge-word generation to avoid combinatorial explosion.
3. **Sign symmetry support** — relax `_validate_pauli_spatial_group` to accept Clifford actions that act on charge words by scalar ±1 (not pure site permutations).
4. **Singlet constraints via variable elimination** — replace `:Zero` constraint-matrix emission with direct variable merging before SDP construction.

---

## User Impact / Severity

- **Blocking for SOTA target**: Without this, the Heisenberg XXX benchmark beyond N≈16 is infeasible (44,851 basis elements at N=100 full NPA).
- **Regression risk**: Low — new code paths for new features; existing `max_degree=2` flow is untouched when no custom basis is supplied.
- **API expansion**: Adds 2-3 new public functions and extends 2 existing specs.

---

## Evidence from Issue and Repository

### Current hard caps (from `src/optimization/symmetry.jl`)

| Location | Constraint |
|----------|-----------|
| Line ~1021 (SymmetrySpec constructor) | `pauli_charge.max_degree in 0:2` hard error |
| Line ~2776 (`_validate_pauli_charge_spec`) | Same `0:2` check |
| Line ~2826 (`_pauli_charge_words`) | `d in 0:2` explicit branch |
| Line ~2990 (`_validate_pauli_spatial_group`) | Rejects any CliffordSymmetry where `_pauli_spatial_permutation(g)` returns `nothing` |

### Basis-closure check (line ~1924)

```julia
function _check_basis_closure(label, basis, group)
    lookup = Set(basis)
    for g in group, mono in basis
        _, image = _act_monomial(g, mono)
        image in lookup || throw(...)
    end
end
```

This checks the monomial image is in the basis. For sign symmetry, `_act_monomial` returns `(sign, image)` where `image == mono` itself (σᵢˣ maps to itself with sign -1 under the Z-product conjugation). The monomial is still in the basis — the sign is carried by the orbit reducer's phase dict. **No change needed to `_check_basis_closure`.**

### Orbit reducer already handles phase (line ~1948)

The `_build_orbit_reducer` tracks `orbit_phase[image_sym]` and detects sign conflicts (`zero_orbit = true`). This means sign symmetry is already structurally supported by the reducer — it just needs the group to pass validation to reach this code.

### Charge-word generation is currently exhaustive (line ~2826)

`_pauli_charge_words(nqubits, max_degree)` generates ALL site combinations for degrees 0, 1, 2. At degree 4 with N=100, this would be C(100,4)·3⁴ ≈ 320M words. The issue proposes generating charge words **only for site-subsets present in the provided basis**, which is exactly 12,000 for contiguous d=4.

### `_pauli_charge_transform_groups` requires `length(charge_words) == length(basis)` (line ~3049)

This equality check assumes the charge basis and the monomial basis are the same size. For a sparse basis, this must be changed — only charge words corresponding to basis site-subsets should be generated, and the transform matrix becomes rectangular if the charge-word count differs from basis size.

Actually, re-reading: this is checking that the charge transformation is a **complete basis change**. If we only generate charge words for site-subsets in the basis, and the basis contains exactly those monomials, the sizes match again. The key insight is that for a contiguous-support basis, the set of relevant charge words is exactly the set of contiguous-support charge words — same count.

### Singlet constraints (line ~3165)

Currently emit `:Zero` constraints via `_append_unique_zero_constraint!`. The issue proposes replacing this with a variable-merging reducer (pre-solver elimination). This is an optimization, not a correctness fix.

---

## Likely Root Cause / Failure Mode

Not a bug — a feature gap. The existing engine is designed for complete bases through degree 2. The architecture (Wedderburn + orbit reducer + charge sectors) is extensible, but the validation guards reject the needed inputs.

---

## Recommended Best-Practice Resolution

Implement in **4 incremental sub-PRs** that can each be merged and tested independently:

### PR 1: Sparse contiguous-support basis constructor (~80 lines new code)

New file `src/optimization/pauli_chains.jl` with:
- `pauli_contiguous_chain_basis(registry, degree; periodic=true)` — returns `Vector{NormalMonomial{PauliAlgebra,T}}`
- Formula: `{1} ∪ { σᵢᵃ¹ ⋯ σᵢ₊ₗ₋₁ᵃˡ : 1 ≤ ℓ ≤ d, a ∈ {x,y,z}^ℓ, i ∈ 1:N }` with periodic wrapping

This is fully independent — it produces a basis vector compatible with the existing `moment_basis` keyword.

### PR 2: Sign symmetry support (~30 lines changed)

- New helper: `pauli_sign_symmetry(N)` — constructs a `CliffordSymmetry` for conjugation by ∏ᵢ σᵢᶻ
- Relax `_validate_pauli_spatial_group`: replace "must be pure site permutation" with "must map each charge word to a scalar (±1) multiple of another charge word with the same charge"
- Adjust `_pauli_charge_word_image` to return `(word, sign)` instead of just `word` when the action flips sign
- The `SymbolicWedderburn.action` for `NCPauliChargeSpatialAction` already returns `(word, sign)` — just need to compute the sign correctly for non-spatial Cliffords

### PR 3: Extend charge to degree 3-4 with basis-driven generation (~100 lines)

- Remove the `0:2` cap from `_validate_pauli_charge_spec` and `SymmetrySpec` constructor
- Add new method `_pauli_charge_words(basis, max_degree)` that:
  1. Extracts site-subsets from the provided basis monomials
  2. For each site-subset of size k ≤ max_degree, generates 3^k charge words
  3. Falls back to exhaustive generation for degree ≤ 2 when no basis is provided (backward compat)
- Thread the basis through `_pauli_charge_transform_groups` (it already receives it)

### PR 4: Singlet variable elimination (optional optimization, ~50 lines)

- New `_build_singlet_reducer(basis, nqubits)` that:
  - Applies 9 SU(2) rotations to each basis monomial
  - Groups equivalent monomials into orbits
  - Returns an `_OrbitReducer`-compatible map
- Integrate into `moment_relax_symmetric` as a second reduction pass
- Compose with the spatial orbit reducer

**Recommendation**: PRs 1-3 are required for the stated goal. PR 4 is a performance optimization that can follow.

---

## Step-by-Step Implementation Plan

### Step 1: `pauli_contiguous_chain_basis`

1. Create `src/optimization/pauli_chains.jl`
2. Implement the chain basis generator with periodic wrapping (`mod1`)
3. Add `include("optimization/pauli_chains.jl")` to `src/NCTSSoS.jl`
4. Export `pauli_contiguous_chain_basis`
5. Create `test/relaxations/pauli_chains.jl` with:
   - Size correctness for N=4,6,10 at various degrees
   - All elements are valid Pauli monomials with contiguous support
   - Identity is included

### Step 2: `pauli_sign_symmetry` + validation relaxation

1. Add `pauli_sign_symmetry(N; integer_type=Int)` to `src/optimization/pauli_chains.jl`
   - Constructs `CliffordSymmetry` where σᵢˣ → (-1, σᵢˣ), σᵢʸ → (-1, σᵢʸ), σᵢᶻ → (1, σᵢᶻ)
2. Replace `_validate_pauli_spatial_group` with `_validate_pauli_charge_compatible_group`:
   - For each generator `g`, check that for every charge word `w` in the active basis, the Clifford action produces `±w'` where `w'` has the same charge as `w`
3. Update `_pauli_charge_word_image(g, word)` to handle non-spatial Cliffords:
   - Factor through the monomial-level Clifford action
   - Decompose result back into charge words
   - Verify it's a scalar (±1) multiple of a single charge word
4. Export `pauli_sign_symmetry`
5. Add tests: basis closure under sign symmetry for contiguous chains

### Step 3: Charge degree 3-4 with basis-driven words

1. Remove `0:2` guards (2 places in validation, 1 in `_pauli_charge_words`)
2. Add `_pauli_charge_words_from_basis(basis, max_degree)`:
   - Extract site-subsets from each basis monomial
   - Generate charge words only for those subsets
   - Return sorted unique charge words
3. Modify `_pauli_charge_transform_groups`:
   - When `spec.max_degree > 2`: require `basis` argument, use `_pauli_charge_words_from_basis`
   - When `spec.max_degree ≤ 2`: use existing exhaustive path (backward compat)
4. Add tests: N=16 d=4 block sizes match 30

### Step 4 (optional): Singlet variable elimination

1. Add `_build_singlet_reducer(basis, nqubits, orbit_reducer)`:
   - For each basis monomial, apply SU(2) rotation group (S₃ × Z₂³ = 48 elements, or just the 9 generators)
   - If image is in basis, merge into same orbit
   - Compose with existing orbit reducer
2. Replace `_add_pauli_singlet_constraints!` call with pre-reducer pass
3. Test: same objective, fewer SDP variables

---

## Files Likely to Change

| File | Change |
|------|--------|
| `src/optimization/pauli_chains.jl` | **New** — chain basis constructor + sign symmetry helper |
| `src/NCTSSoS.jl` | `include` + 2 new exports |
| `src/optimization/symmetry.jl` | Lift degree cap, relax validation, basis-driven charge words, (optional) singlet reducer |
| `test/relaxations/pauli_chains.jl` | **New** — unit tests for chain basis |
| `test/problems/condensed_matter/heisenberg_symmetry.jl` | Add N=16/20/100 d=4 symbolic size assertions |

---

## Regression Tests and Validation Commands

### New tests to add

```
test/relaxations/pauli_chains.jl:
  - "chain basis size N=4 d=1..4"
  - "chain basis size N=10 d=2"
  - "chain basis closure under translation"
  - "chain basis closure under reflection"
  - "chain basis closure under sign symmetry"
  - "N=4 d=2 solve matches order-2 result" (COSMO)

test/problems/condensed_matter/heisenberg_symmetry.jl:
  - "N=16 d=4 max block size ≤ 30" (symbolic, no solver)
  - "N=20 d=4 max block size ≤ 30" (symbolic, no solver, if <60s)
  - "N=100 d=4 max block size ≤ 31" (symbolic, no solver, if <120s)
```

### Validation commands

```bash
# Unit tests only (fast, no solver)
julia --project -e 'using Pkg; Pkg.test(test_args=["relaxations/pauli_chains"])'

# Existing Heisenberg tests still pass
julia --project -e 'using Pkg; Pkg.test(test_args=["problems/condensed_matter"])'

# Full suite
make test
```

### Backward compatibility check

The existing 16-site order-2 test (`heisenberg_symmetry.jl` line 101) must continue passing unchanged — it exercises the full NPA basis with charge degree 2.

---

## Risks, Edge Cases, and Rollback Notes

### Risk 1: Periodic wrapping in charge-word generation

Site-subset `{99, 100, 1, 2}` for a 4-site window wrapping around a 100-site ring must be treated as a valid contiguous subset. The charge-word generator must not assume monotone site ordering equals spatial contiguity.

**Mitigation**: The chain basis generator produces monomials with site indices in `[i, i+1, ..., i+ℓ-1]` (mod N). The charge-word extractor should work from the monomial's site set, not from contiguity assumptions.

### Risk 2: `_pauli_charge_transform_groups` size equality check

Line ~3049: `length(charge_words) == length(basis)` asserts a square transform. This holds if and only if the basis-driven charge words and the basis have the same count. For a contiguous chain basis of degree d with N sites:
- Basis size: `1 + N·∑_{ℓ=1}^d 3^ℓ`
- Charge words generated: same formula (one charge word per basis monomial)

The identity element maps to the identity charge word. Each non-identity basis element maps to exactly one site-subset. Each site-subset generates 3^ℓ charge words = same as 3^ℓ basis monomials on that subset. **Sizes match.** But this must be verified by test.

### Risk 3: Sign symmetry composition with orbit reducer

The sign symmetry generator maps σᵢˣ → -σᵢˣ, σᵢʸ → -σᵢʸ. In a chain monomial with k occurrences of X or Y operators, the eigenvalue is (-1)^k. The orbit reducer will group monomials by their spatial orbit (translation/reflection) and accumulate the sign phase. Monomials with different X+Y counts are in different sign-symmetry eigenspaces and should NOT be merged.

**Mitigation**: The existing `_build_orbit_reducer` already handles this via `zero_orbit` detection when sign conflicts arise in the same orbit. If a monomial with even X+Y count and one with odd count are in the same spatial orbit (they aren't for contiguous chains), the orbit would be zeroed correctly.

### Risk 4: Singlet reducer + orbit reducer composition order

The singlet reducer (SU(2) axis permutations) and the spatial orbit reducer must compose. The simplest safe approach: apply both to find the canonical representative, choosing the lexicographically smallest. The issue proposes doing singlet reduction first as a "pre-reducer" — this avoids composition complexity.

### Rollback

Each sub-PR is independently revertable:
- PR 1 adds new code only (no existing code modified)
- PR 2 relaxes a validation check (revert = restore the `throw`)
- PR 3 removes a cap and adds a branch (revert = restore cap)
- PR 4 adds an optimization (revert = remove reducer call)

---

## Open Questions

### Q1: Should the chain basis support non-periodic (open) chains?

The issue specifies `periodic=true` as default with an option for open boundary conditions. For open chains, the boundary sites have fewer neighbors. The size formula changes slightly. **Recommendation**: Implement periodic first; open is a 5-line extension later.

### Q2: Does the `_pauli_charge_transform_groups` size-equality check need restructuring?

The current code assumes charge words and basis have the same count (square transform matrix). For a basis-driven approach, this is guaranteed if we generate exactly one charge word per basis monomial. But if the user provides an arbitrary sparse basis that doesn't correspond 1:1 to charge words, it breaks.

**Options**:
- (A) Require that custom bases for charge reduction are "charge-complete" (every monomial has a unique charge decomposition into the generated charge basis). Error otherwise. This is the safe conservative path.
- (B) Allow rectangular transforms (more charge words than basis elements or vice versa). This requires nontrivial changes to the Wedderburn machinery downstream.

**Recommendation**: Option A for this PR. The contiguous chain basis is charge-complete by construction.

### Q3: Performance budget for N=100 d=4 symbolic decomposition

The issue claims N=100 d=4 symbolic block decomposition should complete in <120s. The key bottleneck is the Wedderburn decomposition of 12,000 charge words under a group of order 400. SymbolicWedderburn's `symmetry_adapted_basis` scales with basis size × group order. 12,000 × 400 = 4.8M operations — should be feasible but needs profiling.

**Mitigation**: Add a timeout annotation to the test. If it exceeds 120s on CI, mark as `@test_skip` with a note.

### Q4: Should `PauliSingletConstraintSpec.max_degree` cap also be lifted?

The issue extends singlet to degree 4 via variable elimination (PR 4). The existing `_add_pauli_singlet_constraints!` only handles degree 1-2 cases. For degree 3-4, the singlet reducer approach is needed.

**Recommendation**: Lift the `PauliSingletConstraintSpec` cap to 4 only when PR 4 (singlet reducer) is implemented. Until then, keep cap at 2 for the constraint-emission path.
