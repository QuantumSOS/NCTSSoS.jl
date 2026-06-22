# Triage Explanation: Issue #364 — Sparse Chain Basis + Degree 3-4 Charge/Singlet + Sign Symmetry

## TL;DR

Issue #364 is a feature request (not a bug fix) that adds four coordinated capabilities to NCTSSoS.jl — a Julia package for computing quantum ground-state energy bounds using semidefinite programming (SDP) with symmetry reduction. The goal is to enable a benchmark computation (the N=100 periodic Heisenberg XXX model at degree 4) that is currently infeasible because the software's internal basis grows combinatorially (~44,851 elements at N=100 with the current "full" approach) instead of linearly (~12,000 elements with a sparse contiguous-support basis). The triage plan decomposes the work into four independent, incrementally mergeable pull requests that each unlock one capability, with PRs 1–3 required for the stated goal and PR 4 being an optional optimization.

---

## Decision Map

```
Issue #364: Enable N=100 Heisenberg XXX at degree 4
│
├─ Decision 1: Decompose into 4 sub-PRs (not a monolithic change)
│   ├─ Rationale: Each PR is independently testable, mergeable, and revertable
│   └─ Alternative rejected: Single large PR would be harder to review and riskier to revert
│
├─ Decision 2: PR 1 — New sparse basis constructor (no existing code modified)
│   ├─ Rationale: Pure addition; zero regression risk
│   └─ Creates `pauli_contiguous_chain_basis()` producing O(N) elements
│
├─ Decision 3: PR 2 — Relax validation to allow sign symmetry (~30 lines changed)
│   ├─ Rationale: The orbit reducer already handles sign phases internally;
│   │   the only blocker is a validation guard that rejects non-spatial Cliffords
│   └─ Key insight: `_check_basis_closure` needs NO change (sign carried in phase dict)
│
├─ Decision 4: PR 3 — Lift degree cap from 2 → 4, use basis-driven charge words (~100 lines)
│   ├─ Rationale: Exhaustive charge-word generation at degree 4 would produce ~320M words;
│   │   basis-driven generation produces ~12,000
│   └─ Backward compatible: degree ≤ 2 without custom basis still uses exhaustive path
│
└─ Decision 5: PR 4 — Singlet variable elimination (optional, ~50 lines)
    ├─ Rationale: Performance optimization only; not needed for correctness
    └─ Replaces constraint-matrix emission with pre-solver variable merging
```

---

## State Machine Trace

The **triage branch** (`triage/issue-364-feat-pauli-sparse-chain-basis-degree-3-4-charge`) is a Git branch created during the "triage" phase — the phase where the team investigates an issue, determines its scope, and writes a plan before any code changes are made. No code has been written yet; the branch contains only the triage document (`docs/triage/issue-364.md`).

The lifecycle for this issue follows:

1. **Issue opened** (#364, sub-issue of #356, labeled `enhancement` + `SOTA-target`)
2. **Triage** (current phase) — investigation complete, plan written, awaiting human review
3. **Resolve** (next phase) — code changes on a `resolve` branch, following the plan
4. **Validation** — tests pass, backward compatibility confirmed
5. **Merge** — PR(s) merged to main

The triage document is the deliverable of phase 2. It answers: *What needs to change, why, in what order, and what could go wrong?*

---

## Evidence

### Why the current code can't handle N=100 at degree 4

The triage document identifies four hard-coded barriers in `src/optimization/symmetry.jl`:

| Barrier | Location | What it blocks |
|---------|----------|----------------|
| `max_degree ∈ 0:2` validation | Lines ~1021, ~2776 | Any charge sector analysis above degree 2 |
| Exhaustive charge-word generation | Line ~2826 | At degree 4 with N=100 sites, this generates C(100,4)·3⁴ ≈ 320 million charge words — computationally infeasible |
| Spatial-only Clifford validation | Line ~2990 | Rejects sign symmetry (conjugation by ∏ᵢ σᵢᶻ), which is needed for full symmetry exploitation |
| No sparse basis support | No constructor exists | Users must use the full NPA basis, which has O(N⁴) elements at degree 4 |

### Why the architecture is ready for extension

The triage document provides concrete evidence that the *internal machinery* already supports the needed features — only the *validation guards* and *input constructors* are missing:

1. **Orbit reducer handles phases**: At line ~1948, `_build_orbit_reducer` tracks `orbit_phase[image_sym]` and detects sign conflicts via a `zero_orbit` flag. This means sign symmetry is structurally supported; it just can't reach this code because validation rejects sign-symmetry generators.

2. **Basis closure check is phase-aware**: At line ~1924, `_check_basis_closure` checks that the *monomial image* is in the basis set. For sign symmetry, the monomial maps to itself (with a sign), so the monomial is still in the basis. The sign is handled separately by the orbit reducer's phase dictionary.

3. **Charge-word/basis size equality holds for contiguous chains**: The document proves that for a contiguous chain basis of degree d with N sites, the number of basis monomials (`1 + N·∑_{ℓ=1}^d 3^ℓ`) exactly equals the number of basis-driven charge words. This means the existing square-transform assumption in `_pauli_charge_transform_groups` (line ~3049) is not violated.

### What "contiguous chain basis" means

In quantum spin systems, a **Pauli monomial** is a product of Pauli operators (σˣ, σʸ, σᶻ) acting on specific qubit sites. The "full NPA basis" at degree d includes *all* monomials on any subset of up to d sites out of N — yielding O(N^d) elements.

A **contiguous chain basis** restricts to monomials where the sites form a *contiguous block* on the lattice (with periodic wrapping for a ring geometry). For a 1D periodic chain, this yields only O(N·d) distinct site-subsets × 3^ℓ Pauli combinations per subset = O(N) total monomials for fixed d. At N=100, d=4: ~12,000 elements instead of ~44,851.

This is physically motivated: the Heisenberg Hamiltonian only has nearest-neighbor interactions, so a contiguous-support basis captures the relevant correlations without the combinatorial blowup.

---

## Implementation Plan

### PR 1: Sparse Contiguous-Support Basis Constructor (~80 lines, new code only)

**What**: Create a new file `src/optimization/pauli_chains.jl` containing `pauli_contiguous_chain_basis(registry, degree; periodic=true)`.

**Why this is safe**: This PR adds a new file and a new export. No existing code is modified. The function returns a `Vector{NormalMonomial}` that is compatible with the existing `moment_basis` keyword argument.

**Formula**: The basis contains the identity plus all monomials of the form σᵢᵃ¹ ⋯ σᵢ₊ₗ₋₁ᵃˡ where ℓ ranges from 1 to d, each aⱼ ∈ {x,y,z}, and i ranges from 1 to N with periodic wrapping via `mod1`.

**Files changed**:
- `src/optimization/pauli_chains.jl` (new)
- `src/NCTSSoS.jl` (add `include` + export)
- `test/relaxations/pauli_chains.jl` (new test file)

### PR 2: Sign Symmetry Support (~30 lines changed)

**What**: Add a `pauli_sign_symmetry(N)` constructor and relax the validation in `_validate_pauli_spatial_group` to accept Clifford actions that map charge words to ±1 multiples of other charge words (not only pure site permutations).

**Why this route**: The triage document identifies that the orbit reducer *already* handles phase accumulation and sign-conflict detection. The only change needed in the validation layer is to replace "must be a pure site permutation" with "must map each charge word to a scalar ±1 multiple of another charge word with the same charge." This is the minimal relaxation that enables sign symmetry without restructuring the Wedderburn decomposition pipeline.

**Key subtlety**: Under the sign symmetry generator (conjugation by ∏ᵢ σᵢᶻ), a Pauli X or Y operator picks up a factor of -1, while Z is unchanged. A monomial with k occurrences of X or Y gets eigenvalue (-1)^k. Monomials with different X+Y parities live in different eigenspaces and should not be merged — the orbit reducer handles this correctly via `zero_orbit` detection.

**Files changed**:
- `src/optimization/pauli_chains.jl` (add `pauli_sign_symmetry`)
- `src/optimization/symmetry.jl` (relax `_validate_pauli_spatial_group`, adjust `_pauli_charge_word_image`)

### PR 3: Extend Charge Sectors to Degree 3-4 (~100 lines)

**What**: Remove the `0:2` hard cap on `max_degree` and add a new method `_pauli_charge_words_from_basis(basis, max_degree)` that generates charge words only for site-subsets present in the provided basis.

**Why basis-driven generation**: At degree 4 with N=100, exhaustive generation produces C(100,4)·3⁴ ≈ 320M charge words. Basis-driven generation produces ~12,000 — the same count as the basis itself. This is possible because the contiguous chain basis has exactly one charge word per basis monomial (each monomial maps to a unique site-subset + Pauli-type combination).

**Backward compatibility**: When `max_degree ≤ 2` and no custom basis is provided, the existing exhaustive path is used unchanged. The degree-3/4 path is only activated when a basis is explicitly supplied.

**Files changed**:
- `src/optimization/symmetry.jl` (remove caps, add basis-driven generation, modify `_pauli_charge_transform_groups`)

### PR 4: Singlet Variable Elimination (Optional, ~50 lines)

**What**: Replace the current `:Zero` constraint-matrix emission for SU(2) singlet constraints with a pre-solver variable-merging pass. Instead of emitting constraints that the solver must enforce, equivalent variables are merged before SDP construction.

**Why optional**: This is a performance optimization. The existing singlet constraint mechanism is correct; it just adds constraints to the SDP instead of eliminating variables. For N=100, variable elimination reduces the SDP size, but the feature works without it.

**Files changed**:
- `src/optimization/symmetry.jl` (add `_build_singlet_reducer`, integrate as pre-reducer pass)

---

## Validation

### Regression Testing

The existing test at `test/problems/condensed_matter/heisenberg_symmetry.jl` (line ~101) exercises the full NPA basis with charge degree 2 for a 16-site system. This test must pass unchanged after all PRs. Since PR 1 adds new code only, PR 2 relaxes a validation check (existing inputs still pass the relaxed check), and PR 3 preserves the existing exhaustive path for degree ≤ 2, backward compatibility is expected.

### New Tests

| Test | What it validates |
|------|-------------------|
| Chain basis size for N=4,6,10 at degrees 1-4 | Correct element count formula: `1 + N·∑ 3^ℓ` |
| Chain basis closure under translation | Translating any element yields another element in the basis |
| Chain basis closure under reflection | Reflecting any element yields another element in the basis |
| Chain basis closure under sign symmetry | Sign symmetry maps each element to ±1 × (element in basis) |
| N=4 d=2 solve matches existing order-2 result | End-to-end correctness: sparse basis gives same energy bound as full basis |
| N=16 d=4 max block size ≤ 30 | Symbolic (no solver) check that Wedderburn block decomposition produces blocks matching Wang et al. |
| N=100 d=4 max block size ≤ 31 | Symbolic check with timeout annotation (skip if >120s on CI) |

### Validation Commands

```bash
# PR 1 tests only
julia --project -e 'using Pkg; Pkg.test(test_args=["relaxations/pauli_chains"])'

# Existing condensed matter tests (backward compat)
julia --project -e 'using Pkg; Pkg.test(test_args=["problems/condensed_matter"])'

# Full test suite
make test
```

---

## Caveats and Open Questions

### Preserved from Triage Document

**Q1: Should the chain basis support open (non-periodic) chains?**
The triage plan implements periodic chains first. Open boundary conditions change the size formula slightly (boundary sites have fewer neighbors) but are described as a ~5-line extension. No decision has been made on whether to include this in the initial PRs.

**Q2: What if a user provides an arbitrary sparse basis that isn't "charge-complete"?**
The current code assumes `length(charge_words) == length(basis)` (a square transform matrix). The triage plan chooses **Option A** (conservative): require that custom bases are "charge-complete" — i.e., every monomial maps to a unique charge word, and the counts match. Error otherwise. The contiguous chain basis satisfies this by construction. **Option B** (rectangular transforms) would require nontrivial changes to the Wedderburn machinery and is deferred.

**Q3: Will N=100 d=4 symbolic decomposition complete in <120s?**
The triage estimates 12,000 charge words × 400 group elements = 4.8M operations, which should be feasible but hasn't been profiled. A timeout annotation is recommended for the test. If it exceeds 120s on CI, it should be marked `@test_skip`.

**Q4: Should the `PauliSingletConstraintSpec.max_degree` cap also be lifted?**
The existing singlet constraint path (`_add_pauli_singlet_constraints!`) only handles degree 1-2. The triage recommends lifting this cap to 4 **only when PR 4 (singlet reducer) is implemented**, keeping the cap at 2 for the constraint-emission path until then.

### Risk: Periodic Wrapping in Charge-Word Generation

Site-subsets like `{99, 100, 1, 2}` on a 100-site ring must be treated as a valid contiguous-support subset. The charge-word extractor should use the monomial's actual site set (as produced by the chain basis generator), not infer contiguity from monotone ordering. This is called out as a specific risk requiring test coverage.

### Risk: Sign Symmetry Composition with Spatial Orbits

If two monomials with different X+Y parity happened to fall in the same spatial orbit (same site structure under translation/reflection), the orbit reducer would detect a sign conflict and zero the orbit. The triage document asserts this doesn't happen for contiguous chains but recommends defensive testing.

---

## Revision Targets

These are the specific decisions a human reviewer might reasonably want to annotate, challenge, or modify:

1. **4-PR decomposition vs. fewer PRs**: Is the overhead of 4 separate PRs justified, or should PRs 2 and 3 be combined since sign symmetry and degree extension are tightly coupled for the N=100 use case?

2. **Degree cap set to 4 vs. arbitrary**: The plan lifts the cap from 2 to 4. Should it be lifted to an arbitrary degree (with the basis-driven path handling any degree), or is 4 a deliberate ceiling?

3. **Option A (charge-complete requirement) vs. Option B (rectangular transforms)**: The conservative choice restricts custom bases to charge-complete ones. A reviewer may prefer the more general Option B if downstream Wedderburn changes are tractable.

4. **PR 4 marked optional**: The singlet variable elimination is deferred. If the N=100 SDP is too large without it, PR 4 may actually be required for the stated goal. Should it be re-prioritized?

5. **Open boundary condition deferral**: The `periodic=true` default is implemented first. If open chains are needed for other benchmarks in #356, this deferral may need revision.

6. **Timeout threshold (120s) for N=100 symbolic test**: This is an estimate. If CI hardware is significantly slower, this threshold may need adjustment. Should this test be tagged differently (e.g., `@test_broken` vs. `@test_skip`)?

7. **New file vs. extending existing file**: PR 1 creates a new file `src/optimization/pauli_chains.jl`. An alternative is adding the chain basis constructor directly to `symmetry.jl`. The new-file approach is cleaner for separation of concerns but adds a file to the include chain.

8. **Singlet reducer composition order**: PR 4 proposes applying singlet reduction *before* spatial orbit reduction ("pre-reducer"). The alternative — composing both into a single pass choosing the lexicographically smallest canonical representative — is more complex but potentially more optimal. The choice of composition order should be validated with profiling.