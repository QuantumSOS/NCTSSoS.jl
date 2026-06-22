# Triage: Issue #356 — Heisenberg 1D N=100 via Lattice Symmetry Exploitation

**Issue**: [#356](https://github.com/QuantumSOS/NCTSSoS.jl/issues/356)
**Filed**: 2025-06-15 by @exAClior
**Labels**: `enhancement`, `SOTA-target`
**Triage branch**: `triage/issue-356-sota-heisenberg-1d-reach-n-100-via-lattice-symme`

---

## Issue Summary

The goal is to match the scaling of Wang, Tavakoli et al. (arXiv:2604.01555, 2025), who solve the periodic 1D XXX Heisenberg chain at N=100, NPA degree d=4, with max SDP block size 31×31. Their bound is within <10⁻⁵ of DMRG.

NCTSSoS.jl currently works at N=8 dense (~50 min) and N=10 with term sparsity (~40 min). The issue reports that N=20, 40, 60, 100 symmetry runs all fail. The required reduction is from ~8 billion-dimensional moment matrices to 31×31 blocks — roughly 8 orders of magnitude.

## User Impact / Severity

- **Research impact**: High. This is a flagship benchmark for NC polynomial optimization applied to quantum many-body physics. The SOTA papers demonstrate that lattice symmetry exploitation makes industrial-scale spin chain problems tractable.
- **Package credibility**: This is the gap between "research prototype" and "competitive tool." The benchmark is well-defined with published reference values.
- **Scope**: Enhancement, not a bug. But the existing symmetry infrastructure was designed with this target in mind — the 16-site evidence test and Pauli charge/singlet pipeline are stepping stones toward this goal.

## Evidence from Issue and Repository

### What currently works

| N | Method | Basis size | Time | Status |
|---|--------|------------|------|--------|
| 2 | SWAP Clifford symmetry, order 2 | 7 | Fast | ✅ `test/problems/condensed_matter/heisenberg_symmetry.jl` |
| 4 | Charge + spatial + singlet, order 1 | Small | Fast | ✅ Same file |
| 8 | Dense, no symmetry | 3578 | ~50 min | ✅ Issue report |
| 10 | Term sparsity | 612 | ~40 min | ✅ Issue report |
| 16 | Size evidence only (no solve) | 1129 half-basis → 53 blocks, max 24 | N/A (symbolic only) | ✅ Same file; Mosek validation: 2147s total, 9.4s Mosek |
| 20+ | Symmetry adapter | — | — | ❌ "Adapter not working" |

### Existing symmetry infrastructure

The codebase has a mature symmetry reduction pipeline in `src/optimization/symmetry.jl` (3559 lines):

1. **`CliffordSymmetry`** — Represents finite Clifford conjugation actions on Pauli words. `pauli_site_permutation()` constructs spatial actions.
2. **`SymmetrySpec`** — Wraps generators + charge/singlet specs. Enforced: spatial Cliffords only when combined with charge/singlet.
3. **`PauliChargeSectorSpec`** — U(1) charge-sector block decomposition using `{Z, S⁺, S⁻}` basis. **Capped at `max_degree=2`.**
4. **`PauliSingletConstraintSpec`** — Order-2 SU(2) singlet moment equalities. **Capped at `max_degree=2`.**
5. **`CliffordSymmetryGroup`** — Full group enumeration from generators, with Wedderburn decomposition via SymbolicWedderburn.
6. **`_ConstraintMatrixEntryCache`** — Lazy entry computation to avoid materializing full dense matrices (used when `offblock_check=:off`).
7. **SympleQ pipeline** (`src/sympleq/`) — Automated Clifford symmetry detection from Hamiltonian structure.

### Documented bottleneck

From `docs/src/manual/pauli_charge_singlet_symmetry.md`:

> The remaining cost is not Mosek. In the N=16 run, Mosek took only about 9.4 seconds out of a 2146.5 second `cs_nctssos` call. The bottleneck is symbolic/JuMP construction: the current implementation still builds and transforms dense symbolic moment data before emitting the tiny reduced PSD blocks.

### Key code constraints

- `PauliChargeSectorSpec.max_degree` is validated to be in `0:2` (lines 847, 2777-2778 of `symmetry.jl`).
- `PauliSingletConstraintSpec.max_degree` is validated to be in `0:2` (lines 871, 2787-2788).
- `_validate_pauli_spatial_group()` (line 2990) rejects any CliffordSymmetry that is not a pure site permutation when combined with charge/singlet reductions. This blocks sign symmetry (σˣ→−σˣ, σʸ→−σʸ), which is axis-sign-changing but charge-compatible.
- The `moment_relax_symmetric` path at line 3199 builds the orbit reducer over `total_basis`, which for large N at degree 2+ grows as O(N²) in basis size and O(N⁴) in matrix entries.

### Research worktree

A research worktree exists at `NCTSSoS.jl-research-heisenberg-n100-tightening` (branched after the charge/singlet feature landed). This indicates active exploration.

## Likely Root Cause of Scaling Failure

The issue is not a single bug — it's a multi-layered capability gap:

### Layer 1: Dense intermediate representation (blocker for N>16 at order 2)

The `moment_relax_symmetric` function materializes dense polynomial constraint matrices of size `basis_dim × basis_dim` before extracting the tiny diagonal blocks. At N=20 order 2, the Pauli half-basis has ~6860 entries, yielding ~23.5M matrix entries. Each entry is a symbolic polynomial — memory and time explode. The `_ConstraintMatrixEntryCache` (offblock_check=:off mode) partially mitigates this by computing entries lazily, but the transformation loop still touches O(basis_dim² × group_order) entries.

### Layer 2: max_degree=2 cap (blocker for d=3 and d=4)

The charge/singlet infrastructure is hard-coded to degree ≤2. Wang et al. need d=3 and d=4 to get tight bounds. Extending `_pauli_charge_words` to degree 3+ requires generating 3-site and 4-site charge words, and the singlet constraints need higher-order analogues.

### Layer 3: Missing symmetries (gap vs. SOTA block sizes)

Wang et al. combine 7 symmetry reductions multiplicatively. NCTSSoS.jl currently implements 3 of them through the Pauli charge+spatial path:

| Symmetry | Wang et al. | NCTSSoS.jl | Gap |
|---|---|---|---|
| Translation invariance | DFT circulant decomposition | `pauli_site_permutation([2:N;1])` via Wedderburn | **Functional but not circulant-specialized** — generic Wedderburn over Z_N group is correct but may not exploit the diagonal structure as efficiently |
| Mirror (reflection) | Direct block decomposition | `pauli_site_permutation(reverse(1:N))` via Wedderburn | ✅ Works |
| U(1) magnetization | Sector decomposition | `PauliChargeSectorSpec` | ✅ Works at degree ≤2 |
| SU(2) singlet | Moment constraints | `PauliSingletConstraintSpec` | ✅ Works at degree ≤2 |
| Sign symmetry (σˣ→−σˣ, σʸ→−σʸ) | 4-fold factor | **Blocked** — `_validate_pauli_spatial_group` rejects non-site-permutation Cliffords combined with charge/singlet | ❌ Missing |
| xyz permutation (SU(2)) | 3-fold factor | Not implemented | ❌ Missing |
| Conjugate symmetry | Complex→real SDP, 2× reduction | Not implemented | ❌ Missing |
| ncKKT (state optimality) | Additional PSD constraints | Not implemented | ❌ Missing — from Araujo 2023, Fawzi 2024 |

### Layer 4: Group enumeration scaling

`_enumerate_symmetry_group` explicitly enumerates all group elements. For the dihedral group D_N (translation + reflection), the group order is 2N. At N=100, that's 200 elements, which is fine. But adding sign symmetry and xyz permutation would push to 200 × 4 × 3 = 2400 elements — still manageable for Wedderburn. The real cost is the basis transformation, not the group enumeration.

## Recommended Resolution: Phased Approach

This is too large for a single PR. The recommended path is a sequence of progressively harder, independently valuable milestones.

### Phase 0: Diagnose the N=20 failure (1-2 days)

**Goal**: Identify the exact failure mode when running the existing symmetry pipeline at N=20, order 2.

1. Run the perf benchmark script at N=20 with timing:
   ```bash
   NCTS_PERF_NS=20 julia --project perf/pauli_charge_singlet_prep.jl
   ```
2. If it OOMs, profile memory. If it hangs, profile time.
3. The likely bottleneck is either:
   - Group enumeration of D₂₀ (40 elements × 6860-dimensional domain)
   - `_build_orbit_reducer` over the full basis
   - `_build_constraint_matrix` materializing the dense matrix
   - Wedderburn decomposition inside `_pauli_charge_transform_groups`

**Deliverable**: Precise profiling data identifying which function consumes the most time/memory.

### Phase 1: Eliminate the dense intermediate (1-2 weeks)

**Goal**: Make N=20-40 work at order 2 with existing symmetries.

The `_ConstraintMatrixEntryCache` + `offblock_check=:off` path already avoids materializing the full dense matrix for the Pauli charge path (lines 3416-3421). The work needed:

1. **Make the lazy path the default for Pauli charge problems** — the offblock_check=:off path with `_diagonal_transformed_constraint_blocks` computes only the entries needed for the diagonal blocks. The orbit reducer should also be built lazily or in reduced coordinates.

2. **Build the orbit reducer in reduced coordinates** — `_build_orbit_reducer` currently iterates over `total_basis`, which includes unreduced monomials. For a translation-invariant system, the number of orbit representatives is O(basis_dim / group_order), dramatically smaller.

3. **Profile and optimize `_pauli_charge_transform_groups`** — this is where the charge-sector splitting meets the spatial Wedderburn decomposition. Ensure it doesn't materialize basis_dim² intermediate data.

**Files likely to change**:
- `src/optimization/symmetry.jl`: `_build_orbit_reducer`, `moment_relax_symmetric`, `_pauli_charge_transform_groups`, `_diagonal_transformed_constraint_blocks`

**Regression tests**:
- Extend the N=16 size evidence test to N=20, N=40 (symbolic only, no solve)
- Add wall-time regression guard: N=20 symbolic reduction should complete in <60s

### Phase 2: Extend charge/singlet to degree 3-4 (1-2 weeks)

**Goal**: Lift the `max_degree=2` cap.

1. **Extend `_pauli_charge_words`** to degree 3 and 4:
   - Degree 3: all ordered triples of sites with {Z, S⁺, S⁻} operators
   - Degree 4: all ordered 4-tuples
   - Count: degree-d charge words for N sites = ∑_{k=0}^{d} C(N,k)·3^k

2. **Extend `_add_pauli_singlet_constraints!`** to degree 3-4:
   - Higher-order SU(2) singlet equalities relate multi-point correlation functions
   - These are linear moment constraints, so the framework generalizes

3. **Update validation**: Change `0:2` range checks to `0:4` in `_validate_pauli_charge_spec` and `_validate_pauli_singlet_spec`.

**Files likely to change**:
- `src/optimization/symmetry.jl`: `_pauli_charge_words`, `_add_pauli_singlet_constraints!`, validation functions

**Regression tests**:
- Test charge basis construction at degree 3-4 for small N (N=4-6)
- Verify basis size matches theoretical count
- N=4, order 3-4 end-to-end solve with charge/singlet (small enough for COSMO)

### Phase 3: Sign symmetry (σˣ→−σˣ, σʸ→−σʸ) (3-5 days)

**Goal**: 4-fold additional reduction factor.

Sign symmetry is conjugation by ∏ᵢ σᵢᶻ. It maps:
- σˣᵢ → −σˣᵢ
- σʸᵢ → −σʸᵢ
- σᶻᵢ → σᶻᵢ

In the charge basis:
- Z → Z (charge 0, invariant)
- S⁺ → −S⁺ (charge +1, sign flip)
- S⁻ → −S⁻ (charge −1, sign flip)

This is **compatible** with the charge decomposition — it acts as multiplication by (−1)^(number of S⁺/S⁻ operators) within each sector. The current rejection by `_validate_pauli_spatial_group` is overly conservative.

1. **Relax `_validate_pauli_spatial_group`** to allow sign-compatible Clifford symmetries that act as scalar multiples within each charge sector. This requires a refined check: the group element must map each charge word to a scalar multiple of another charge word of the same charge.

2. **Add sign symmetry as a Clifford generator**: `CliffordSymmetry` with σˣᵢ→(−1,σˣᵢ) and σʸᵢ→(−1,σʸᵢ) for all i.

**Files likely to change**:
- `src/optimization/symmetry.jl`: `_validate_pauli_spatial_group`, `_pauli_charge_word_image`

**Regression tests**:
- N=4 Heisenberg with sign symmetry: verify additional block splitting
- Compare bounds with and without sign symmetry

### Phase 4: Conjugate symmetry (complex→real) (3-5 days)

**Goal**: 2× reduction by casting Hermitian PSD to real PSD.

The Heisenberg Hamiltonian has real coefficients in the charge basis. When the moment matrix is Hermitian with real structure, it can be decomposed into real and imaginary parts, halving the SDP size.

1. This may be implementable at the JuMP/solver level rather than the symmetry pipeline — detect when all PSD blocks are real-certifiable and emit real PSD constraints instead of Hermitian PSD.

**Files likely to change**:
- `src/optimization/lowering.jl` or `src/optimization/sos.jl`: real/Hermitian cone selection
- Possibly `src/optimization/symmetry.jl`: real-structure detection

### Phase 5: SU(2) xyz permutation (1-2 weeks)

**Goal**: 3-fold reduction from isotropic (XXX) structure.

For the XXX Heisenberg model, the Hamiltonian is invariant under global SU(2) rotations. The group of Pauli-axis permutations (x→y→z→x) is a subgroup of this. This is an axis-mixing Clifford that breaks the `{Z, S⁺, S⁻}` charge basis.

Two approaches:
1. **Non-Abelian SU(2) decomposition**: Replace the Abelian U(1) charge decomposition with a full SU(2) irrep decomposition. This is the mathematically clean approach but requires substantial new representation theory machinery.
2. **Post-charge xyz equivalence**: After charge+spatial reduction, identify blocks that are equivalent under the xyz permutation and drop redundant copies. This is simpler but less general.

**Files likely to change**: Major additions to `src/optimization/symmetry.jl` or a new `src/optimization/su2_symmetry.jl`.

### Phase 6: ncKKT (state optimality constraints) (2+ weeks)

**Goal**: Additional PSD constraints from state optimality conditions.

This implements the approach of Araujo (2023) and Fawzi (2024), adding first-order KKT conditions for the ground-state eigenvalue problem as additional PSD constraints. This is orthogonal to the symmetry pipeline — it adds constraints to the MomentProblem.

### Milestone targets

| Phase | N target | Order | Max block | Comparable to |
|-------|----------|-------|-----------|---------------|
| 0 | Diagnose N=20 failure | — | — | Understanding |
| 1 | N=20-40, order 2 | 2 | ~24 (N=16 reference) | Current capability, working |
| 2 | N=20, order 3-4 | 3-4 | TBD | Wang 2024 at small N |
| 1+2+3 | N=40-60, order 3-4 | 3-4 | TBD | Wang 2024 at medium N |
| 1+2+3+4+5 | N=100, order 4 | 4 | ~31 | Wang 2025 Table 1 |

## Step-by-Step Implementation Plan (Phase 0+1)

These are the concrete next steps that should be tackled first:

1. **Profile N=20 failure**: Run `perf/pauli_charge_singlet_prep.jl` at N=20. Capture the exact error or timeout.
2. **Identify the O(N⁴) bottleneck location**: Is it `_build_orbit_reducer`, `_build_constraint_matrix`, `_pauli_charge_transform_groups`, or the Wedderburn decomposition?
3. **Implement lazy orbit reduction**: Build the orbit representative map only for monomials that appear in the reduced blocks, not the full basis.
4. **Test `offblock_check=:off` + `_ConstraintMatrixEntryCache` at N=20**: Verify this path works end-to-end.
5. **Add N=20 symbolic size evidence test** to `test/problems/condensed_matter/heisenberg_symmetry.jl`.
6. **Add N=40 symbolic size evidence test** (same file).
7. **Benchmark end-to-end Mosek solve at N=20, order 2** with charge+spatial+singlet.

## Files Likely to Change

| File | Phase | Change |
|------|-------|--------|
| `src/optimization/symmetry.jl` | 1-5 | Core symmetry pipeline: orbit reducer, charge basis, validation |
| `perf/pauli_charge_singlet_prep.jl` | 0 | Extend benchmark to larger N |
| `test/problems/condensed_matter/heisenberg_symmetry.jl` | 1+ | Additional size evidence tests |
| `docs/src/manual/pauli_charge_singlet_symmetry.md` | 1+ | Update known limits and bottleneck section |
| `src/optimization/lowering.jl` | 4 | Real PSD cone detection |
| `src/optimization/moment.jl` | 1 | Lazy constraint matrix construction |

## Regression Tests and Validation Commands

```bash
# Current tests pass:
make test

# Specific symmetry tests:
julia --project -e 'using Pkg; Pkg.test()' -- --testset="heisenberg_symmetry"

# Perf benchmark (no solver):
NCTS_PERF_NS=8,12,16,20 julia --project perf/pauli_charge_singlet_prep.jl

# After Phase 1, add to CI:
# N=20 symbolic size evidence (no solver, fast)
# N=40 symbolic size evidence (no solver, <120s)
```

## Risks, Edge Cases, and Rollback Notes

- **Risk: Wedderburn numerical stability at large group order**. The dihedral group D₁₀₀ has order 200. The Wedderburn decomposition must be numerically stable at this scale. SymbolicWedderburn handles this generically, but the Pauli-specific charge sector path implements its own splitting.

- **Risk: Charge basis explosion at degree 4**. The number of degree-4 charge words is C(N,4)·3⁴ ≈ 3.2M for N=100. The charge transformation matrices will be large. The lazy evaluation path must extend to this regime.

- **Risk: Sign symmetry compatibility with charge basis**. The claim that sign symmetry is charge-compatible needs formal verification. Edge case: when the sign flip produces a phase factor in the charge-word expansion, the scalar multiplication within each sector must be verified.

- **Rollback**: Each phase is independently valuable and independently testable. Phase 0 is diagnostic only. Phases 1-2 are backward-compatible (existing tests continue to pass). Phases 3-5 add new SymmetrySpec options without changing default behavior.

## Open Questions

### Q1: What is the exact failure mode at N=20?

**What we know**: The issue says "adapter not working at these sizes." The N=16 evidence test runs successfully (symbolic only, 2147s with Mosek solve). The perf script (`perf/pauli_charge_singlet_prep.jl`) supports `NCTS_PERF_NS=8,12,16` by default.

**What is ambiguous**: Is the failure an OOM, a timeout, an error, or a numerical issue? The answer determines whether Phase 1 is a memory optimization, a time optimization, or a correctness fix.

**Options**:
- (a) Run the perf benchmark at N=20 and capture the failure → determines Phase 1 scope
- (b) Try `optimizer=nothing` at N=20 to see if the symbolic reduction completes → isolates solver vs. symbolic issue

### Q2: Does the research worktree contain relevant progress?

**What we know**: `NCTSSoS.jl-research-heisenberg-n100-tightening` exists, branched after the charge/singlet feature (#360). It has 10 commits visible.

**What is ambiguous**: Whether it contains partial implementations of sign symmetry, higher-degree charge, or the dense-intermediate fix.

**Action**: Check the research worktree for progress before starting Phase 1 implementation.

### Q3: Is the Wang et al. DFT-based translation invariance fundamentally different from Wedderburn on Z_N?

**What we know**: Wang et al. use DFT/circulant block diagonalization for translation invariance. NCTSSoS.jl uses generic Wedderburn decomposition of the cyclic group Z_N. Both produce the same block structure mathematically (irreducible representations of Z_N are one-dimensional, corresponding to the Fourier modes).

**What is ambiguous**: Whether the generic Wedderburn code path is asymptotically as efficient as a specialized DFT implementation, or whether it has overhead that matters at N=100.

**Recommendation**: The generic path should be fine since Z_N has only N irreps and the Wedderburn decomposition is exact for Abelian groups. Profile before specializing.

### Q4: Should higher-degree singlet constraints be moment equalities or full SU(2) Clebsch-Gordan decomposition?

**What we know**: The current order-2 singlet implementation adds linear moment equalities (〈σᵢᵃ〉=0, cross-component zeros, same-component equalities). Wang et al. use SU(2) in a more structural way.

**What is ambiguous**: Whether extending the linear-constraint approach to degree 3-4 is sufficient, or whether a representation-theoretic approach (Clebsch-Gordan coefficients for SU(2) tensor products) is needed to match the SOTA block sizes.

**Recommendation**: Start with the linear constraint extension (Phase 2) and benchmark the gap. If block sizes are still too large, implement the CG decomposition in Phase 5.
