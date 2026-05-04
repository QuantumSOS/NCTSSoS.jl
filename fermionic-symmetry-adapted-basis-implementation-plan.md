# Fermionic Symmetry-Adapted Basis Implementation Plan for `NCTSSoS.jl`

Reviewed against:

- note: `/Users/exaclior/MyBrain/YushengBrain/topics/nctssos-fermionic-symmetry/nctssos-fermionic-symmetry.typ`
- `src/optimization/symmetry.jl`
- `src/optimization/interface.jl`
- `src/optimization/moment.jl`
- `src/simplification/fermionic.jl`
- `src/types/monomial.jl`
- `test/relaxations/symmetry.jl`
- `test/problems/fermionic/fermionic.jl`
- `test/problems/condensed_matter/hubbard.jl`

---

## Implementation progress

- [x] Phase 0 — baseline capture and correctness blockers
- [x] Phase 1 — generic transformed-block symmetry pipeline
- [x] Phase 2 — fermionic low-level action support
- [x] Phase 3 — sector decomposition for parity / `U(1)` / `S_z`
- [x] Phase 4 — Abelian orbital symmetry composition
- [x] Phase 5 — SU(2) spin-adapted basis transform

### Low-cost benchmark ladder

- [x] Step 1 — two-orbital singlet/triplet structural SU(2) regression
- [x] Step 2 — H₂/STO-3G molecular-orbital regression
- [x] Step 3 — open two-site Hubbard dimer end-to-end benchmark

### Work log

- 2026-04-22: Started implementation on branch `feat/fermionic-symmetry-adapted-basis`.
- 2026-04-22: Loaded parent session excerpts and repo plan; beginning with symmetry-path correctness cleanup before widening support.
- 2026-04-22: Completed Phase 0 cleanup:
  - replaced the signed-permutation key with a collision-free representation;
  - added symmetry-support extraction that preserves signed fermionic operator indices;
  - routed symmetry-reduced fermionic relaxations through `_add_parity_constraints!`.
- 2026-04-22: Completed Phase 1 transformed-block refactor:
  - generalized block transforms to unitary/complex `U * M * U†`;
  - removed the scalar-only reduction path and now keep full reduced blocks;
  - extended `SymmetryReport` with block provenance/labels.
- 2026-04-22: Completed Phase 2 finite fermionic action support:
  - added `FermionicModePermutation` on physical modes;
  - generalized finite-group enumeration/Wedderburn plumbing beyond raw `SignedPermutation`;
  - added PBW-capable fermionic polynomial action support.
- 2026-04-22: Completed Phase 3 fermionic sector splitting:
  - added `FermionicModeLayout`, `FermionicSectorSpec`, and `FermionicSectorLabel`;
  - added exact sector-block splitting for parity / number / `S_z` with explicit zero constraints on forbidden cross-sector entries;
  - added Hubbard and spinless-fermion regressions comparing symmetric vs ordinary relaxations.
- 2026-04-22: Phase 4 / 5 remained intentionally loud-failing at the end of the first pass:
  - orbital-irrep composition (`split_irrep=true`) errored explicitly instead of pretending;
  - no automatic SU(2) CG / Be benchmark path was landed in that pass.
- 2026-04-23: Completed Phase 4 Abelian irrep composition:
  - added `AbelianIrrepTable` and active-metadata validation;
  - extended `FermionicModeLayout` with explicit irrep composition rules;
  - implemented exact irrep labels in `_fermionic_sector_label` and `SymmetryReport`.
- 2026-04-23: Completed Phase 5 SU(2) Casimir-based spin adaptation:
  - added `FermionicSpinAdaptationSpec` and deterministic `S^2` eigenspace handling;
  - composed spin transforms with the existing transformed-block pipeline;
  - converted non-scalar spin off-blocks into explicit zero constraints instead of pretending they were zero polynomials.
- 2026-04-23: Landed a CI-grade Be-derived benchmark slice:
  - reproduced the note’s Be 1s/2s worked-example `[2, 2]` spin blocks as a structural regression;
  - added a reviewed STO-6G 1s/2s one-body `N=2` numerical slice with expectation fixture and symmetric-vs-plain objective checks;
  - the full 4-electron / order-2 Be acceptance benchmark is still a larger follow-up, not something to fake.
- 2026-04-23: Expanded this plan with a fully specified low-cost benchmark ladder for:
  - H₂/STO-3G,
  - the two-orbital singlet/triplet SU(2) toy,
  - and the open two-site Hubbard dimer.
- 2026-04-22/23: Added regression coverage in:
  - `test/relaxations/symmetry.jl` for signed support extraction, collision-free keys, complex/unitary block transforms, fermionic mode actions, Abelian irrep labels, spin-Casimir transforms, and spin-adapted zero-constraint behavior;
  - `test/relaxations/interface.jl` for unsupported algebra/action combinations, including invalid spin-adaptation specs;
  - `test/problems/condensed_matter/hubbard.jl` for symmetry-adapted canonical two-site Hubbard;
  - `test/problems/fermionic/fermionic.jl` for parity + fermionic symmetry smoke coverage;
  - `test/problems/fermionic/fermionic_symmetry.jl` plus `test/data/expectations/fermionic_symmetry.toml` for end-to-end irrep sectors and the Be-derived symmetry benchmark slice.
- 2026-04-22: Validation run summary:
  - targeted symmetry/interface/CHSH/Hubbard subsets passed;
  - standalone `test/problems/fermionic/fermionic.jl` passed after fixing stale test assumptions;
  - full `make test` passed.
- 2026-04-23: Landed benchmark ladder step 1:
  - added the two-orbital, two-electron, `M_S = 0` direct SU(2) transform regression in `test/relaxations/symmetry.jl`;
  - verified the exact `[3, 1]` singlet/triplet split with `2S ∈ {0, 2}` using projector/span checks instead of brittle sign tests;
  - added reviewed structural fixture `two_orbital_pair_sector_spin_blocks`.
- 2026-04-23: Landed benchmark ladder step 2:
  - added H₂/STO-3G molecular-orbital metadata / irrep label regressions for `(σ_g, σ_u)` in `test/problems/fermionic/fermionic_symmetry.jl`;
  - added a solver-backed `N = 2` order-2 benchmark comparing exact diagonalization, ordinary dense relaxation, and the symmetry-adapted path;
  - committed reviewed fixtures `h2_sto3g_sector_labels` and `h2_sto3g_n2_spin_adapted`, explicitly documenting that the stored energy is electronic-only (no nuclear-repulsion constant).
- 2026-04-23: Landed benchmark ladder step 3:
  - added the open two-site Hubbard dimer benchmark in `test/problems/condensed_matter/hubbard.jl` without removing the older periodic-two-site coverage;
  - compared dense, sector-only, and sector + site-swap + spin-adapted order-2 paths against exact diagonalization and the analytic dimer energy;
  - committed reviewed fixtures `open_n2_u4_dimer_order2` and `open_n2_u4_dimer_order2_spin_symmetry`.
- 2026-04-23: Validation run summary for the benchmark ladder work:
  - `test/relaxations/symmetry.jl` passed with the new `[3, 1]` SU(2) structural regression;
  - `test/problems/fermionic/fermionic_symmetry.jl` passed with the new H₂/STO-3G checks plus the pre-existing Be slices;
  - `test/problems/condensed_matter/hubbard.jl` passed with the new open-dimer benchmark and the older periodic cases still intact;
  - full `make test` passed.

## Historical roadmap for the landed symmetry work

Slices A and B below are now landed, and the low-cost benchmark ladder added later
in this file is also landed. The only substantial acceptance target that still
remains intentionally open is the full Be minimal-basis / 4-electron / order-2
numerical benchmark from Slice C.

Anything that tries to jump straight to that final Be acceptance target without
first nailing the irrep and SU(2) plumbing is cargo-cult engineering.

### Slice A — Abelian orbital-irrep composition

#### Goal

Turn `split_irrep=true` from an explicit error into a real exact sector label.

#### What is missing today

`FermionicModeLayout` can store `irrep_of`, but the symmetry layer still has no
mathematics for how irreps compose across creation and annihilation operators.
Right now that field is just a bag of metadata with no semantics.

#### Required design decision

Introduce an explicit Abelian-irrep composition contract instead of hoping `Any`
will somehow become a group law by positive thinking.

Recommended addition:

```julia
struct AbelianIrrepTable{Γ}
    identity::Γ
    multiply::Function   # (Γ, Γ) -> Γ
    dual::Function       # Γ -> Γ
end
```

and extend `FermionicModeLayout` / `FermionicSectorSpec` to carry it.

Why this is the right seam:

- creation operators should contribute `Γ_p`;
- annihilation operators should contribute `dual(Γ_p)`;
- the label of a monomial is the ordered product in the Abelian irrep group;
- the neutral-selection rule is then exact and obvious.

For real Abelian point groups, `dual(Γ) == Γ` is often enough. But it must be a
first-class rule, not a hidden assumption.

#### Implementation steps

1. Add `AbelianIrrepTable` and validation helpers.
2. Tighten `FermionicModeLayout` / `FermionicSectorSpec` so `split_irrep=true`
   requires:
   - orbital ids for every active mode,
   - irreps for every active orbital,
   - an explicit irrep table.
3. Extend `_fermionic_sector_label` so the irrep part is computed exactly.
4. Remove the current loud-fail branch for `split_irrep=true`.
5. Keep the current sector-neutral constraint rule: for now, only neutral
   objective / localizing / moment-equality multipliers are accepted.
6. Extend `SymmetryReport.block_labels` to include the irrep label in testable form.

#### Tests

Add two cheap deterministic regressions:

1. **Pure label test**
   - tiny spin-orbital toy with orbitals tagged by Abelian irreps;
   - assert that monomials land in the expected irrep sectors.

2. **End-to-end symmetry test**
   - small fermionic Hamiltonian with a site-swap finite symmetry and orbital-irrep labels;
   - assert same optimum as ordinary path;
   - assert reduced blocks split by irrep as expected.

#### Files

- `src/optimization/symmetry.jl`
- probably new helper file once this grows: `src/optimization/fermionic_symmetry.jl`
- `test/relaxations/symmetry.jl`
- `test/problems/fermionic/fermionic.jl` or a new dedicated fermionic-symmetry file

#### Acceptance

- `split_irrep=true` no longer errors when metadata is complete;
- incomplete metadata still errors clearly;
- irrep labels show up in `SymmetryReport`;
- objective agrees with the unsymmetrized problem.

---

### Slice B — SU(2) spin-coupled basis transforms

#### Goal

Add a real spin-adapted basis transform on top of the existing `(parity, number,
S_z, irrep)` sector split.

#### Important choice

Do **not** start with Be-specific hard-coded CG tables.

The simplest sound first implementation for this repo is:

- keep the current fermionic sector split by `S_z`,
- build the **adjoint-action total-spin operators** on each fixed sector,
- diagonalize the spin Casimir inside that sector,
- use the resulting orthonormal eigenbasis as the transform `U`.

This fits the current architecture because the repo already knows how to:

- build symbolic polynomial matrices,
- keep complex/unitary transformed blocks,
- and report exact block provenance.

#### Why this is better than jumping straight to hand-coded CG tables

Because the half-basis here is a basis of **operators**, not an electron-configuration
basis in a quantum-chemistry package. The natural SU(2) action is the operator-space
adjoint action:

- `L_z(O) = [S_z, O]`
- `L_±(O) = [S_±, O]`
- `L² = L_z^2 + 1/2 (L_+ L_- + L_- L_+)`

That gives a generic, algebra-native path to spin adaptation.

#### Scope discipline

Do this in two substeps.

##### B1. `S^2`-resolved transforms inside fixed `S_z` sectors

This is the first honest milestone.

Implementation plan:

1. Add helpers to build total-spin operators from `FermionicModeLayout`:
   - `S_z`
   - `S_+`
   - `S_-`
2. For a fixed basis block and fixed `(parity, number, S_z, irrep)` label:
   - compute the action matrices of `L_z`, `L_+`, `L_-` on the basis;
   - build `L²`;
   - diagonalize `L²` numerically but deterministically;
   - partition eigenvectors by `S(S+1)` within tolerance.
3. Emit `BasisTransformBlock`s labelled by `(parity, number, 2S, 2M_S, irrep, multiplicity_tag)`.
4. Feed those transforms through the existing `U * M * U†` block path.

This lands useful SU(2)-aware structure without pretending we already have a
full reduced-tensor implementation.

##### B2. Multiplet-aware reduced basis across `M_S`

This is what you need for the nice compact Be-style block structure.

Implementation plan:

1. For each `S` multiplet, pick deterministic highest-weight vectors.
2. Generate the full multiplet by laddering with `L_-`.
3. Either:
   - keep the repeated `M_S` blocks but label them explicitly, or
   - compress to one representative block per `S` only when the symbolic entries
     are provably equivalent.

For the first useful Be target, it is acceptable to **keep multiplicity explicit**
and delay aggressive compression. What is not acceptable is silently identifying
blocks that are only approximately equal.

#### Determinism requirements

Eigen-decomposition is where things usually get stupid.

To keep the result testable:

- symmetrize action matrices before diagonalization;
- group eigenvalues with a fixed tolerance;
- within degenerate eigenspaces, choose a deterministic basis by:
  - sorting pivot support lexicographically,
  - QR/SVD with fixed column ordering,
  - and explicit sign/phase normalization.

If this sounds fussy, good. Numerical basis drift is exactly how you get tests
that pass on Tuesday and fail on Thursday.

#### Suggested API addition

Add a separate transform spec instead of cramming SU(2) into `FermionicSectorSpec`:

```julia
struct FermionicSpinAdaptationSpec
    mode_layout::FermionicModeLayout
    casimir_tol::Float64
    compress_multiplets::Bool
end
```

and let `SymmetrySpec` grow a field like `spin_adaptation::Union{Nothing,FermionicSpinAdaptationSpec}`.

That keeps:

- sector splitting,
- finite actions,
- and basis transforms

as three distinct ideas. Which they are.

#### Files

Strong recommendation: stop growing `src/optimization/symmetry.jl` into a kitchen sink.
Split out spin-specific code, for example:

- `src/optimization/fermionic_spin.jl`
- maybe also `src/optimization/fermionic_irreps.jl`

and include them from `src/NCTSSoS.jl`.

#### Tests

Add a dedicated spin-adaptation test file, e.g.

- `test/problems/fermionic/fermionic_symmetry.jl`

with low-cost milestones:

1. one-body spin doublet basis resolves into the expected `S = 1/2` sectors;
2. two-body singlet/triplet combinations emerge correctly in a tiny 2-orbital toy;
3. spin-adapted and ordinary relaxations give the same objective;
4. spin-adapted blocks are strictly smaller than the uncoupled `(N, S_z)` version.

#### Acceptance

- `S^2` labels are present in `SymmetryReport`;
- transformed blocks are deterministic across runs;
- ordinary vs spin-adapted objectives agree;
- no unsupported multiplicity case is silently compressed.

---

### Slice C — Be minimal-basis benchmark

#### Goal

Use the machinery from slices A and B to land the actual benchmark the note cares about.

#### Preconditions

Do **not** start this until both are true:

- Abelian irrep composition works,
- SU(2) transform blocks work in cheap toy regressions.

Otherwise the Be test is just a one-off stunt.

#### Benchmark plan

1. Encode the minimal spin-orbital layout explicitly:
   - spatial orbital ids,
   - spin labels,
   - any required orbital symmetry labels.
2. Build the benchmark Hamiltonian / moment constraints from reviewed literature data.
3. Add two acceptance levels:

##### C1. structural acceptance

- 1-RDM / order-1 block structure matches the expected singlet/triplet split;
- report labels show the expected `(N, 2S, 2M_S, Γ)` content.

##### C2. numerical acceptance

- order-2 spin-adapted relaxation matches the reviewed benchmark value within solver tolerance;
- unsymmetrized and adapted objectives agree;
- adapted blocks are strictly smaller.

#### Fixtures

Use reviewed expectation data, likely in:

- `test/data/expectations/fermionic_symmetry.toml`

Suggested case ids:

- `toy_two_orbital_singlet_triplet_order1`
- `be_minimal_order1_spin_blocks`
- `be_minimal_order2_spin_adapted`

#### Acceptance

The feature is only honestly “done” when Be passes without any hard-coded special-case
transform path.

---

## Detailed low-cost benchmark ladder (now landed) below the full Be target

These three examples are the right staircase below the full Be benchmark.
They are not interchangeable, and pretending otherwise is how you end up with
one giant “toy test” that proves nothing.

The intended order is:

1. **H₂/STO-3G** — smallest molecular literature example.
2. **Two-orbital singlet/triplet structural test** — smallest genuinely useful SU(2) test.
3. **Open two-site Hubbard dimer** — best small end-to-end benchmark for this codebase.

### General rules for all three examples

- Keep the output contract as an ordinary symmetry-reduced `MomentProblem`.
- Do **not** add a runtime PySCF/OpenFermion dependency just to run tests.
  - If molecular integrals are needed, generate/review them once offline,
    then commit the derived arrays / reviewed expectations.
- Be explicit about whether a reference value is:
  - **electronic energy only**, or
  - **electronic + nuclear/core constant**.
  Do not silently mix them.
- Reuse the existing `FermionicModeLayout` / `FermionicSectorSpec` /
  `FermionicSpinAdaptationSpec` seams.
- Keep finite discrete actions, Abelian irrep labels, and SU(2) basis transforms
  as separate layers even when a test uses more than one of them.
- Structural transform tests belong in `test/relaxations/symmetry.jl`.
- Solver-backed fermionic literature / benchmark tests belong in:
  - `test/problems/fermionic/fermionic_symmetry.jl` for H₂ and tiny abstract spin tests;
  - `test/problems/condensed_matter/hubbard.jl` for the Hubbard dimer.
- Reviewed numeric fixtures should live in:
  - `test/data/expectations/fermionic_symmetry.toml` for H₂ and the two-orbital SU(2) toy;
  - `test/data/expectations/hubbard.toml` for the Hubbard dimer.
- For each example, decide **CI / docs / both** up front instead of letting it drift:
  - H₂/STO-3G: **CI yes**, docs later if a public narrative is useful.
  - two-orbital singlet/triplet: **CI yes**, docs no; this is an internal structural test.
  - Hubbard dimer: **CI yes**, docs optional later if a symmetry tutorial page is added.

### Current status of this ladder

All three low-cost ladder examples below are now landed in CI as of 2026-04-23:

- H₂/STO-3G → landed in `test/problems/fermionic/fermionic_symmetry.jl` with
  fixtures `h2_sto3g_sector_labels` and `h2_sto3g_n2_spin_adapted`;
- two-orbital singlet/triplet → landed in `test/relaxations/symmetry.jl` with
  fixture `two_orbital_pair_sector_spin_blocks`;
- open two-site Hubbard dimer → landed in `test/problems/condensed_matter/hubbard.jl`
  with fixtures `open_n2_u4_dimer_order2` and `open_n2_u4_dimer_order2_spin_symmetry`.

What still remains open is the larger full Be minimal-basis / 4-electron / order-2
acceptance benchmark, not the ladder below.

### Example 1 — H₂/STO-3G (smallest molecular literature example; landed 2026-04-23)

#### Why this case exists

This is the smallest real molecular v2RDM / spin-adapted benchmark in the
literature. It is the right first “chemistry-shaped” test because it is tiny,
reviewable, and does not need Be-sized machinery.

#### What symmetry it should exercise here

Use it to test:

- fixed particle number `N = 2`;
- `S_z = 0`;
- `S^2 = 0` singlet adaptation;
- optional Abelian orbital irrep labels distinguishing `σ_g` and `σ_u`.

Do **not** force a finite permutation generator into this case unless there is
an exact orbital relabelling symmetry you actually use. H₂ is valuable even as a
sector-only / spin-adapted molecular example.

#### Literature anchor

Use DePrince’s H₂/STO-3G v2RDM tutorial as the primary literature source.
The implementation in this repo does not need to mimic the tutorial’s exact SDP
packing, but the integrals / reference energy must be traceable to a reviewed
source or reviewed offline derivation.

#### Concrete implementation choice

Implement H₂ in **molecular-orbital form**, not AO form and not via a runtime HF call.
That means:

1. store the reviewed one-electron integrals `h[p,q]` and two-electron integrals
   `(pq|rs)` for the chosen fixed geometry in a small helper or inline constant;
2. build the fermionic Hamiltonian directly as a polynomial over four spin orbitals;
3. use `moment_eq_constraints` to enforce `N = 2`;
4. compare against an exact diagonalization helper built in the test file.

That keeps the test deterministic and keeps the package out of the business of
being a quantum-chemistry front-end.

#### Required representation details

Use exactly:

- **2 spatial orbitals**: `σ_g`, `σ_u`;
- **4 spin orbitals**: `(σ_g↑, σ_g↓, σ_u↑, σ_u↓)`;
- orbital ids:
  - both `σ_g` spin orbitals → orbital id `1`;
  - both `σ_u` spin orbitals → orbital id `2`;
- spin labels:
  - up → `+1`, down → `-1` in `spin2_of`.

Optional irrep labels for `split_irrep=true`:

- `σ_g` → neutral/even irrep, e.g. `:Ag`;
- `σ_u` → odd irrep, e.g. `:B1u`.

Use an explicit two-element Abelian irrep table. Do **not** rely on string hacks.

#### Hamiltonian construction details

Build

```math
H = \sum_{pq,σ} h_{pq} a^†_{pσ} a_{qσ}
  + \frac12 \sum_{pqrs,στ} (pq|rs) a^†_{pσ} a^†_{qτ} a_{sτ} a_{rσ}
```

with the following implementation rules:

- skip coefficients below a fixed tiny threshold when constructing the polynomial;
- keep the coefficient arrays as reviewed constants in the test helper;
- if the reviewed source gives a **total** energy, either:
  - add the nuclear repulsion constant explicitly to the polynomial objective, or
  - keep the polynomial electronic-only and store the **electronic** reference value.
  Pick one and document it in the fixture notes.

#### Suggested basis / relaxation choices

Use two test layers.

##### H2-A. structural / label test

No solver required.

- Construct representative one-body monomials `a^†_{pσ} a_{qτ}`.
- Assert `_fermionic_sector_label` gives the expected
  `(parity, number, S_z, irrep)` labels.
- In particular, the diagonal densities `a^†_{gσ} a_{gσ}` and `a^†_{uσ} a_{uσ}`
  should be sector-neutral in irrep, while transition operators like
  `a^†_{gσ} a_{uσ}` should carry the non-neutral orbital irrep.

This is the smallest exact check that the Abelian irrep composition is not fake.

##### H2-B. solver-backed molecular test

Use the ordinary polynomial path with:

- `moment_basis = [1; all one-body number-conserving monomials]` as the first cheap pass;
- if the quartic Hamiltonian requires the next level for a stable lower bound,
  move to the smallest order / basis that matches the reviewed exact value reliably.

The symmetry spec should be:

```julia
SymmetrySpec(
    sector=FermionicSectorSpec(
        mode_layout=layout,
        split_parity=true,
        split_number=true,
        split_spin=true,
        split_irrep=true,
    ),
    spin_adaptation=FermionicSpinAdaptationSpec(mode_layout=layout),
)
```

#### Exact-reference strategy

Do **not** use literature energy text blindly.

Instead:

1. review the source integrals once;
2. build an exact finite-dimensional fermionic matrix in the test;
3. project to the `N = 2` sector;
4. take the minimum eigenvalue as the exact reference for the polynomial actually used.

That gives a deterministic oracle even if the literature reports the value in a
slightly different convention.

#### Expected test assertions

Minimum required assertions:

1. exact diagonalization matches the reviewed fixture;
2. symmetry-adapted objective matches the ordinary objective;
3. symmetry report includes `FermionicSpinBlockLabel`s with `2S = 0`;
4. if `split_irrep=true` is enabled, the report / direct-label checks show the
   expected `σ_g` vs `σ_u` irrep behavior;
5. adapted blocks are not larger than the unsymmetrized ones.

#### File / fixture plan

- add helper(s) in `test/problems/fermionic/fermionic_symmetry.jl`;
- add fixture ids to `test/data/expectations/fermionic_symmetry.toml`, for example:
  - `h2_sto3g_sector_labels`
  - `h2_sto3g_n2_spin_adapted`

If structural data needs richer fields than the shared oracle loader exposes,
parse the TOML directly in the test. Do not contort the loader.

---

### Example 2 — two-orbital singlet/triplet (smallest genuinely useful SU(2) test; landed 2026-04-23)

#### Why this case exists

This is the smallest case where `S^2` does real work instead of just producing
a cosmetic label. It is the first test that should fail if the Casimir-based
spin adaptation is wrong.

#### What this test should *not* be

Do **not** make this a “molecule” with borrowed integrals just to sound fancy.
It is a structural SU(2) test first.

#### Exact operator sector to use

Use a two-spatial-orbital, two-electron, `M_S = 0` pair-creation sector over
four spin orbitals `(1↑, 1↓, 2↑, 2↓)`.

Use the canonical monomials produced by the package for:

- `a^†_{1↑} a^†_{1↓}`
- `a^†_{2↑} a^†_{2↓}`
- `a^†_{1↑} a^†_{2↓}`
- `a^†_{2↑} a^†_{1↓}`

in whatever sign / order the canonical fermionic simplifier returns.

That basis has size `4` and should decompose under SU(2) as:

- one **triplet** `M_S = 0` channel of size `1`;
- three **singlet** channels of total size `3`.

So the transformed block sizes should be exactly:

- **`[3, 1]`**

with labels `2S = 0` and `2S = 2`.

That is the smallest nontrivial SU(2) block test worth having.

#### Why `[3, 1]` is the right answer

Because:

- double occupancies on a single orbital are forced singlets;
- the mixed-orbital `M_S = 0` subspace splits into one singlet and one triplet;
- total dimension is therefore `2 + 1 + 1 = 4`, grouped as singlet `3` plus triplet `1`.

#### Implementation details

This should be a **direct transform test**, not a solver-backed optimization.

Implementation steps:

1. build `FermionicModeLayout` with orbital ids `{1, 2}` and spin labels `±1`;
2. construct a `FermionicSectorSpec` with at least:
   - `split_parity=true`
   - `split_number=true`
   - `split_spin=true`
3. construct `FermionicSpinAdaptationSpec(mode_layout=layout)`;
4. call the internal transform helper on the four basis monomials;
5. assert the returned block labels and row-basis sizes.

#### Sign / phase handling rule

Do **not** hard-code analytic coefficient signs for the mixed singlet/triplet
vectors without respecting this package’s canonical fermionic ordering.

The test should instead assert one of the following:

- exact block sizes `[3, 1]` and exact `2S` labels; and/or
- projector / span equality for the singlet and triplet subspaces; and/or
- deterministic basis rows after the package’s own phase normalization rule.

What you should **not** do is write a brittle sign test that fails because the
canonical monomial order flips a global phase.

#### Optional second layer

If a tiny solver-backed follow-up is wanted, add a toy Hamiltonian whose
low-energy space is the singlet channel and whose first excited space is the triplet.
But that is optional. The structural `[3, 1]` transform test is the required part.

#### File / fixture plan

- add this primarily to `test/relaxations/symmetry.jl`;
- if structural fixture data is desired, use `test/data/expectations/fermionic_symmetry.toml` with ids like:
  - `two_orbital_pair_sector_spin_blocks`
  - `two_orbital_pair_sector_spin_labels`

---

### Example 3 — open two-site Hubbard dimer (best small benchmark for this codebase; landed 2026-04-23)

#### Why this case exists

This is the right small end-to-end physics benchmark for this repo because:

- it is naturally fermionic;
- it already fits the package’s existing Hubbard / lattice-test style;
- it has exact diagonalization and even closed-form energies;
- it exercises **finite discrete symmetry + particle number + `S_z` + `S^2`** in one honest example.

#### Important choice

Use the **open dimer**, not the current periodic two-site ring.

That means one bond `(1,2)`, not a periodic construction that effectively counts
that edge twice. If you accidentally reuse the periodic ring here, you are no
longer testing the textbook Hubbard dimer.

#### Hamiltonian

Use

```math
H = -t \sum_{σ \in \{↑,↓\}} (c^†_{1σ} c_{2σ} + c^†_{2σ} c_{1σ})
    + U (n_{1↑} n_{1↓} + n_{2↑} n_{2↓})
```

with the default benchmark choice:

- `t = 1.0`
- `U = 4.0`

#### Exact-reference value

At half-filling with `N = 2` and open boundary conditions, the singlet ground-state
energy is analytically:

```math
E_0 = \frac{U}{2} - \sqrt{\left(\frac{U}{2}\right)^2 + 4 t^2}
```

So for `t = 1`, `U = 4`:

- **`E0 = -0.8284271247461903`**

Use this closed form in the fixture notes and independently verify it by exact
matrix diagonalization in the test helper.

#### Symmetry content to exercise

This case should use all of:

- `N = 2` enforced by `moment_eq_constraints`;
- `S_z = 0` from sector splitting;
- `S^2` from `FermionicSpinAdaptationSpec`;
- site-swap `1 ↔ 2` as a finite `FermionicModePermutation` action.

This is precisely the kind of compositional test the codebase needs.

#### Required mode layout

For spin orbitals `(1↑, 2↑, 1↓, 2↓)` or the repo’s equivalent indexing:

- both spin orbitals on site 1 → orbital id `1`;
- both spin orbitals on site 2 → orbital id `2`;
- up → `+1`, down → `-1` in `spin2_of`.

Do **not** set `split_irrep=true` here unless you are deliberately adding an
additional Abelian orbital label on top. The finite site-swap action already
supplies the discrete symmetry worth testing.

#### Symmetry spec

Use:

```julia
SymmetrySpec(
    fermionic_generators=[site_swap],
    sector=FermionicSectorSpec(
        mode_layout=layout,
        split_parity=true,
        split_number=true,
        split_spin=true,
    ),
    spin_adaptation=FermionicSpinAdaptationSpec(mode_layout=layout),
)
```

where `site_swap` exchanges site 1 and site 2 separately for up and down modes.

#### Solver configuration

Use a solver-backed regression at the smallest honest order for a quartic Hamiltonian:

- default target: **order 2**;
- compare:
  1. ordinary dense path,
  2. sector-only path,
  3. sector + site-swap + spin-adapted path.

The test should not merely check “it runs”. It should compare the three paths.

#### Required assertions

Minimum required assertions:

1. exact diagonalization and closed-form analytic energy agree;
2. ordinary dense objective matches the exact ground-state energy within solver tolerance;
3. symmetry-adapted objective matches the ordinary objective;
4. symmetry-adapted maximum block size is strictly smaller than the unsymmetrized one;
5. symmetry report includes `FermionicSpinBlockLabel`s with at least the singlet `2S = 0` channel;
6. if the moment basis includes the triplet sector, the report should also expose `2S = 2` labels.

#### Where this belongs

This belongs in `test/problems/condensed_matter/hubbard.jl`, not in the general
fermionic-symmetry file.

Reason: it is a condensed-matter benchmark first and a symmetry test second.
The benchmark should live near the rest of the Hubbard coverage.

#### Fixture plan

Add reviewed fixture entries to `test/data/expectations/hubbard.toml`, for example:

- `open_n2_u4_dimer_order2`
- `open_n2_u4_dimer_order2_spin_symmetry`

If extra structural fields are needed beyond `objective`, `sides`, and `nuniq`,
store them in the TOML case and parse them directly in the test.

#### Relationship to existing Hubbard coverage

Do **not** replace the existing periodic-two-site ring test.
Keep both.

They protect different things:

- periodic two-site ring: existing package behavior / current benchmark lineage;
- open Hubbard dimer: smallest exact end-to-end symmetry benchmark.

---

## Completed implementation order for these three examples

The intended staircase was followed in practice:

1. **two-orbital singlet/triplet structural `[3, 1]` test**
   - fastest way to catch SU(2) mistakes,
   - no solver noise,
   - no external integrals.
2. **H₂/STO-3G**
   - first real molecular literature example,
   - proved the machinery survives contact with reviewed molecular coefficients.
3. **open Hubbard dimer**
   - best small end-to-end benchmark for this repo,
   - composes site-swap + sector splitting + spin adaptation cleanly.

That sequence gave a real benchmark ladder instead of one giant Be-shaped ball of mud.

---

## Current next implementation target

The low-cost ladder is done. The next honest acceptance target is now:

1. **full Be minimal-basis benchmark acceptance**
   - structural Be checks are already present;
   - the remaining work is the larger 4-electron / order-2 numerical acceptance case,
     with reviewed expectations and no hard-coded special-case transform path.

That is the remaining step worth caring about.

## Executive verdict

The note is directionally right, but its first implementation slice needs revision for this codebase.

`NCTSSoS.jl` already has a working symmetry MVP for dense **monoid** problems. The real job is **not** “add symmetry from scratch”; it is:

1. turn the current symmetry path into a **generic transformed-block pipeline**,
2. add **fermionic sector logic** for parity / particle number / `S_z`,
3. then add a **spin-coupled basis transform** for SU(2).

The clean implementation is **not**:

- jamming SU(2) into `SignedPermutation`,
- pretending a continuous symmetry is just a larger finite group,
- or hard-coding a Be-specific spin projector.

That would be a shortcut. Shortcuts are how you get a demo and then spend the next month undoing your own cleverness.

The correct north star is:

- keep the ordinary `MomentProblem` / solver pipeline,
- make symmetry a transformation on the symbolic relaxation,
- separate **sector splitting** from **basis transforms**,
- and only use `SymbolicWedderburn` where it actually fits: finite discrete actions on a basis that really is permuted linearly.

---

## What the note gets right

These parts of the note are solid and should stay as design anchors:

- the gap is real: `src/optimization/symmetry.jl` is monoid-only today;
- Jordan–Wigner is the wrong seam for spin adaptation;
- SU(2) and multiplicity are the real hard parts, not syntax;
- the Be minimal-basis benchmark is a good **final acceptance target**;
- v2RDM / QSpace / GUGA are the right literature lineages.

---

## What must be revised

### 1. “Lift the type gate” is necessary but nowhere near sufficient

The note’s step 1 is too optimistic for the current repo.

Today’s symmetry code is monoid-specific in several structural ways:

- `_act_monomial`, `_act_polynomial`, `_check_symmetry_invariance`, `_check_basis_closure`, `_build_orbit_reducer`, `_reduce_constraint_matrix_symmetric`, and `moment_relax_symmetric` are all restricted to `A<:MonoidAlgebra`.
- `_symmetry_domain` currently relies on `variable_indices`, which intentionally uses `abs(idx)` for PBW sparsity bookkeeping. That is correct for sparsity, but wrong for fermionic symmetry because it erases creation vs annihilation.
- `_signed_image_signature` collapses `(sign, dst)` to `sign * Int(dst)`, which is unsafe once signed fermionic operator indices matter.
- `moment_relax_symmetric` currently skips `_add_parity_constraints!`, so the fermionic parity rule would be silently lost even if the type gate were removed.

So no, “just widen `<:MonoidAlgebra` to `<:PBWAlgebra`” is not a plan. It is how you get incorrect results with a pleasant API.

### 2. `SignedPermutation` is not the long-term fermionic symmetry abstraction

It is still useful for:

- finite mode permutations,
- sign-carrying orbital relabellings,
- the current monoid MVP.

It is **not** enough for:

- `U(1)` particle number,
- `SU(2)` spin,
- complex/unitary basis changes,
- irreducible tensor operators,
- multiplicity-aware non-Abelian reductions.

So the repo should not contort fermionic symmetry into “more `SignedPermutation`”.

### 3. The note’s “step 4 needs a wider `MomentProblem`” is too pessimistic

`MomentProblem` already supports full matrix-valued PSD/HPSD constraints.

That means the first serious symmetry refactor should **remove the scalar-only assumption** in the symmetry path and allow reduced blocks of size `> 1`, instead of forcing everything through `1×1` scalar blocks.

That is a big deal:

- it makes spin-adapted fermionic blocks possible without rewriting the solver interface,
- it also improves the current monoid symmetry path,
- and it avoids a fake bottleneck that does not actually exist in the rest of the package.

### 4. SU(2) should enter through a basis transform, not through fake finite-group enumeration

For this codebase, the clean seam for spin adaptation is:

- construct a **spin-coupled half-basis**,
- transform the symbolic moment/localizing matrices by a unitary basis change,
- keep the resulting reduced block matrices as ordinary polynomial-valued PSD/HPSD constraints.

Do **not** start with a discrete-grid SU(2) projector wired into the current finite-group action API.

That approach is numerically awkward, conceptually wrong for the current architecture, and will be miserable to test deterministically.

---

## Current repository audit

### What is already in good shape

#### `src/optimization/moment.jl`

The ordinary moment builder is already fairly algebra-generic:

- `_build_constraint_matrix` already expands PBW products via `simplify(A, _neat_dot3(...))`;
- complex Hermitian problems are already supported;
- `MomentProblem` already stores full block matrices;
- `_add_parity_constraints!` already enforces fermionic parity superselection in the ordinary path.

This is good news. It means the symmetry feature should be implemented as a **reduction/transformation layer**, not by rewriting the basic moment machinery.

#### `test/problems/fermionic/fermionic.jl`

There is already a parity baseline and a small sanity net for fermionic moment relaxations.

That is useful and should be extended, not bypassed.

#### `test/problems/condensed_matter/hubbard.jl`

There are already realistic, number-conserving fermionic benchmarks with `moment_eq_constraints`.

That gives us a natural intermediate acceptance path before the Be benchmark.

### What is structurally wrong today for fermionic symmetry

#### `src/optimization/symmetry.jl`

This file currently mixes three different concerns that need to be teased apart:

1. **finite group enumeration and invariance checking**;
2. **basis decomposition / block transforms**;
3. **orbit reduction of symbolic polynomial entries**.

That was tolerable for the scalar monoid MVP. It is too cramped for fermionic symmetry.

Specific blockers:

- monoid-only dispatch throughout the symmetry path;
- scalar-only / multiplicity-free assumptions on reduced blocks;
- real-only transform helper `_transform_polynomial_block`;
- no conjugation-aware unitary transform path for complex fermionic blocks;
- no signed fermionic support extraction;
- no parity injection in the symmetric path.

---

## Design principles for the fermionic feature

### 1. Keep the ordinary symbolic output contract

The reduced problem should still be an ordinary symbolic `MomentProblem`.

That means:

- keep `solve_sdp` intact,
- keep `sos_dualize` intact,
- keep the core basis and polynomial types intact.

Symmetry should reduce the symbolic relaxation, not replace the package’s algebraic core.

### 2. Separate **sector splitting** from **basis transforms**

These are different things and should stay different:

- **sector splitting**: parity, particle number, `S_z`, Abelian irreps;
- **basis transforms**: Wedderburn blocks for finite groups, Clebsch–Gordan coupled bases for SU(2).

Mixing them into one fake action type is how designs rot.

### 3. Do not infer physics from variable names

Spin, orbital labels, and point-group irreps should come from explicit metadata, not from parsing `"c_up"` and hoping for the best.

Introduce explicit fermionic mode metadata.

### 4. Make the first fermionic slice useful **and** extensible

The staged path should be:

- parity / `U(1)` / `S_z` sectors first,
- Abelian orbital symmetry next,
- SU(2) coupled basis next,
- non-Abelian point groups / multiplicity compression after that.

That is not timidity. It is how you land the feature without painting yourself into a corner.

---

## Target architecture

The clean design is a three-layer symmetry system.

### Layer A: symmetry actions and metadata

Add explicit mode metadata for fermions.

Sketch:

```julia
struct FermionicModeLayout{T,Γ}
    orbital_of::Dict{T,Int}      # physical mode -> spatial orbital id
    spin2_of::Dict{T,Int}        # physical mode -> ±1 for α/β, or 0 if absent
    irrep_of::Dict{Int,Γ}        # spatial orbital id -> Abelian point-group irrep
end
```

This belongs in the symmetry spec, not the registry. The registry is for symbols; the symmetry layer is for physics.

For finite orbital relabellings, add a fermionic action that acts on **physical modes**, not raw signed operators:

```julia
struct FermionicModePermutation{T}
    mode_image::Dict{T,T}        # map on physical mode ids
end
```

It should induce the corresponding map on operator letters:

- `a_p      -> a_{g(p)}`
- `a†_p     -> a†_{g(p)}`

and then re-normal-order with `simplify(FermionicAlgebra, ...)`.

That keeps CAR semantics in one place.

### Layer B: sector labels

For the first serious fermionic slice, basis elements should carry additive quantum-number labels.

Sketch:

```julia
struct FermionicSectorLabel{Γ}
    parity::Bool      # odd/even operator parity
    number::Int       # creation count - annihilation count
    spin2::Int        # twice ΔS_z
    irrep::Γ          # Abelian orbital irrep product, optional
end
```

This is the right abstraction for:

- parity,
- `U(1)` particle number,
- `S_z`,
- Abelian point groups.

The selection rule for a matrix entry `⟨q_i† g q_j⟩` is then expressed in terms of sector compatibility, not ad-hoc casework.

### Layer C: transformed-block reduction

Generalize the current symmetry code to work with explicit basis transforms.

Sketch:

```julia
struct BasisTransformBlock{L}
    label::L
    U::Matrix{ComplexF64}        # rows = transformed basis vectors in original basis
end
```

Then a generic helper applies:

- `U_i * M * U_j†` for off-block checks,
- `U_i * M * U_i†` for diagonal blocks,
- and returns the resulting symbolic block matrices.

This becomes the common engine for:

- current `SymbolicWedderburn` discrete reductions,
- future fermionic sector blocks,
- future SU(2) spin-coupled basis transforms.

That is the right abstraction. One transform engine, multiple block providers.

---

## Public API direction

Do **not** break today’s CHSH API.

But do evolve `SymmetrySpec` into a composite spec with backward-compatible constructors.

### Backward-compatible direction

Keep:

```julia
SymmetrySpec(generators::SignedPermutation...)
```

as sugar for the current monoid/discrete path.

Extend toward something like:

```julia
SymmetrySpec(
    sectors = [...],
    finite_actions = [...],
    basis_transforms = [...],
    check_invariance = true,
)
```

No need to nail the exact field names today, but the direction matters: the spec must be able to describe more than finite signed permutations.

### Reporting

`SymmetryReport` should grow block metadata.

At minimum, each reduced block should know:

- size,
- label / sector,
- provenance (`:sector_split`, `:wedderburn`, `:su2_cg`, ...).

The current CHSH-oriented report is too small for fermionic work.

---

## Implementation roadmap

## Phase 0 — baseline capture and correctness blockers

### Goal

Lock down current behavior and fix the bugs that would make fermionic symmetry silently wrong.

### Changes

- add a symmetry-specific support extractor that preserves signed operator indices;
- replace `_signed_image_signature` / `_symmetry_key` with a collision-free representation;
- route fermionic symmetric reductions through `_add_parity_constraints!`;
- add a regression for the fact that the current symmetric path would otherwise miss parity zero constraints.

### Files

- `src/optimization/symmetry.jl`
- `src/optimization/moment.jl`
- `test/relaxations/symmetry.jl`
- `test/problems/fermionic/fermionic.jl`

### Acceptance

- all existing monoid symmetry tests still pass;
- a fermionic symmetry smoke test proves parity constraints are still present after the symmetric path runs;
- the new support extraction distinguishes `+i` from `-i`.

---

## Phase 1 — generic transformed-block symmetry pipeline

### Goal

Stop forcing symmetry reduction to produce only scalar blocks.

### Changes

Refactor the current symmetry path so that it can return **full reduced PSD/HPSD blocks**.

Concretely:

- generalize `_transform_polynomial_block` to complex/unitary transforms;
- use conjugation correctly (`U * M * U†`), not only real left/right multiplication;
- factor block transformation away from the `SymbolicWedderburn`-specific logic;
- replace `_reduce_to_scalar_blocks_sw` with a generic block-transform path;
- keep the current Wedderburn monoid path as one block-provider implementation.

### Files

- `src/optimization/symmetry.jl`
- `test/relaxations/symmetry.jl`
- `test/problems/bell_inequalities/chsh_simple.jl`

### Acceptance

- CHSH symmetry still returns the same optimum and report metadata;
- current scalar-block cases still work unchanged;
- a new test proves non-`1×1` transformed blocks can be carried through the symbolic pipeline.

This phase is non-negotiable. Fermionic spin-adapted blocks are impossible without it.

---

## Phase 2 — fermionic low-level action support

### Goal

Make the symmetry layer PBW-capable where the action is still monomial-preserving after normal ordering.

### Changes

- add a PBW-capable `_act_polynomial` path that can return a polynomial, not just a signed monomial;
- add fermionic action validation on physical modes;
- widen symmetry support from “all monoid algebras” to “supported algebra/action combinations”, explicitly including fermionic mode permutations;
- keep unsupported combinations failing loudly.

### Important constraint

This phase is still about **finite mode relabellings** and low-level correctness.

It is **not** SU(2) yet.

### Files

- `src/optimization/interface.jl`
- `src/optimization/symmetry.jl`
- `src/simplification/fermionic.jl` only if helper extraction is needed
- `test/relaxations/symmetry.jl`
- `test/problems/fermionic/fermionic.jl`

### Acceptance

- a tiny fermionic mode-swap symmetry test gives the same objective as the unsymmetrized path;
- the symmetric result keeps parity constraints;
- unsupported fermionic actions error cleanly.

Suggested first end-to-end problem: a 2-site number-conserving Hubbard slice or equally small even-parity fermionic Hamiltonian.

---

## Phase 3 — sector decomposition for parity, `U(1)`, and `S_z`

### Goal

Add the first genuinely useful fermionic symmetry reduction: basis blocks labelled by operator quantum numbers.

### Changes

- compute `FermionicSectorLabel` for each basis monomial;
- partition the half-basis and localizing bases by sector;
- use charge compatibility rules for `⟨q_i† g q_j⟩`;
- expose reduced block labels in the symmetry report;
- keep the reduction exact and symbolic.

### Why this phase matters

This is the right machinery for:

- parity adaptation,
- fixed particle number,
- fixed spin projection,
- later Abelian orbital symmetry.

And it does **not** require pretending these symmetries are finite signed-permutation groups.

### Files

- `src/optimization/symmetry.jl`
- `src/optimization/interface.jl`
- new tests in `test/relaxations/` and `test/problems/condensed_matter/hubbard.jl`
- expectation fixture, likely `test/data/expectations/fermionic_symmetry.toml`

### Acceptance

Use existing Hubbard coverage as the first serious regression target:

1. `test/problems/condensed_matter/hubbard.jl`
   - same objective as the ordinary path,
   - smaller moment/localizing block structure,
   - stable metadata for particle-number / `S_z` sectors.

2. canonical half-filled Hubbard case with `moment_eq_constraints`
   - sector labels consistent with the fixed-number problem,
   - no regression in solver behavior.

---

## Phase 4 — Abelian orbital symmetry composition

### Goal

Compose fermionic sector splitting with Abelian orbital irreps.

### Changes

- attach point-group irrep labels to spatial orbitals via `FermionicModeLayout`;
- extend sector labels with orbital irrep products;
- allow optional finite orbital permutations within each sector;
- use the existing discrete symmetry engine only where it actually applies.

### Why now

This gives the v2RDM-style “easy win” symmetries first:

- `U(1)` number,
- `S_z`,
- Abelian point groups.

That is the right bridge before SU(2).

### Acceptance

A small, deterministic toy problem with explicit orbital irrep labels should reduce into the expected Abelian blocks.

If a chemistry-style toy system is available cheaply, use it. If not, a hard-coded spin-orbital example is fine. What matters is the symmetry reduction, not fancy integrals.

---

## Phase 5 — SU(2) spin-adapted basis transform

### Goal

Land the actual spin-adapted basis feature for fermionic algebras.

### Design choice

Implement SU(2) through **coupled basis transforms** built from Clebsch–Gordan coefficients, not through fake finite-group enumeration.

### Changes

- add a generator for spin-coupled fermionic basis vectors using the explicit spin-orbital layout;
- start with 1-body operator sectors and then extend to the order-2 moment basis;
- transform symbolic blocks using the generic complex/unitary transform path from Phase 1;
- keep the resulting diagonal blocks as ordinary polynomial-valued HPSD constraints;
- record `(N, 2S, 2M_S, Γ)` or the appropriate reduced labels in the report.

### Important scope discipline

This phase should target **spin-adapted block structure**, not the full QSpace universe in one shot.

The point is to land a clean SU(2)-aware half-basis and reduced symbolic blocks.

### Deliverable benchmark

Use the note’s Be atom benchmark as the acceptance target, but do it in a way that fits this repo:

- hard-code the minimal spin-orbital metadata and published benchmark values;
- do **not** require a full molecular integral parser to claim success.

### Acceptance

Create a dedicated regression, likely in a new file such as:

- `test/problems/fermionic/fermionic_symmetry.jl`

with reviewed expectations in:

- `test/data/expectations/fermionic_symmetry.toml`

Minimum assertions:

1. the Be minimal-basis 1-RDM spin-adapted block structure matches the expected `[2, 2]` singlet/triplet split from the note;
2. the order-2 symmetry-adapted relaxation matches the published benchmark value within solver tolerance;
3. the unsymmetrized and symmetry-adapted formulations agree on the objective;
4. the symmetry-adapted formulation has strictly smaller PSD/HPSD blocks.

This is the first point where the feature is honestly “there”.

---

## Phase 6 — optional multiplicity compression and non-Abelian point groups

### Goal

Make the implementation elegant for the hard cases, not just correct.

### Changes

- compress repeated multiplicity blocks rather than merely keeping repeated transformed blocks;
- add cached CG / `6j` / `9j` data where needed;
- extend beyond Abelian point groups once the SU(2) machinery is settled.

### Status

This is explicitly **after** the first useful fermionic feature lands.

Do not hold the first spin-adapted Be benchmark hostage to every future non-Abelian ambition.

---

## File-by-file implementation map

### `src/optimization/symmetry.jl`

This remains the main seam.

It should be refactored to own:

- symmetry support extraction,
- invariance checks,
- sector partitioning,
- block transform application,
- optional finite-group/Wedderburn decomposition,
- reduced `MomentProblem` emission.

But the internal code should be separated conceptually into:

1. support / action helpers,
2. sector logic,
3. transform logic,
4. report construction.

If the file becomes unmanageable, split it into subfiles. Do **not** invent a mini framework just to feel clever.

### `src/optimization/interface.jl`

Replace today’s coarse “monoid only” guardrail with capability-based checks:

- monoid + finite signed permutations,
- fermionic + supported sector rules,
- fermionic + supported mode permutations,
- later fermionic + SU(2) coupled basis transforms.

Unsupported combinations must still error early.

### `src/optimization/moment.jl`

Keep the ordinary path intact.

Only shared hooks should be reused:

- polynomial matrix construction,
- parity constraints,
- cone handling,
- block constraint appending.

### `src/types/*`

Avoid deep type-hierarchy churn unless forced.

Do **not** rewrite `NormalMonomial` or `Polynomial` for phase 1–5. The whole point is that you should not need to.

---

## Test plan

### Low-level regression tests

Extend `test/relaxations/symmetry.jl` with:

- signed-support extraction for fermions;
- collision-free symmetry keying;
- complex/unitary block transform tests;
- PBW fermionic action tests;
- regression that the symmetric fermionic path still injects parity zero constraints.

### Existing end-to-end tests to strengthen

- `test/problems/fermionic/fermionic.jl`
- `test/problems/condensed_matter/hubbard.jl`

These should gain symmetry-enabled cases rather than spawning redundant near-copies.

### New end-to-end acceptance test

Add one dedicated fermionic symmetry test file for the Be benchmark and any tiny spin-adapted toy cases needed to get there.

### Expectation fixtures

Use a reviewed fixture file, likely:

- `test/data/expectations/fermionic_symmetry.toml`

Suggested case ids:

- `parity_sector_hubbard_n2_order2`
- `u1_sz_sector_hubbard_n4_order2`
- `be_minimal_1rdm_spin_blocks`
- `be_minimal_order2_spin_adapted`

---

## What not to do

Do **not** do any of these:

1. **Do not** widen the type gate and hope for the best.
2. **Do not** encode SU(2) as a `SignedPermutation` workaround.
3. **Do not** infer spin or orbital symmetry from symbol names.
4. **Do not** hard-code the Be transform without a general spin-orbital metadata layer.
5. **Do not** rebase the whole package around linear-combination basis elements before proving the feature.
6. **Do not** use Jordan–Wigner as the implementation seam for symmetry adaptation.
7. **Do not** keep the scalar-only symmetry assumption once full block matrices are already supported everywhere else.

---

## Recommended PR sequence

### PR 1 — symmetry correctness cleanup

- signed support extraction
- collision-free group keying
- parity constraints in symmetric path
- no public API expansion yet

### PR 2 — generic transformed-block support

- complex/unitary transform helper
- remove scalar-only restriction in the symmetry path
- keep CHSH green

### PR 3 — fermionic low-level action support

- PBW-capable symmetry action path
- fermionic mode-permutation support
- small fermionic symmetry smoke tests

### PR 4 — parity / `U(1)` / `S_z` sector blocks

- sector labels
- block partitioning
- Hubbard end-to-end regressions

### PR 5 — Abelian orbital symmetry

- orbital irrep metadata
- sector composition
- small deterministic orbital-symmetry regression

### PR 6 — SU(2) spin-adapted basis

- CG-based coupled basis transform
- Be benchmark acceptance test

### PR 7+ — multiplicity compression / non-Abelian point groups

- only after the first useful feature is done

---

## Final recommendation

Implement the fermionic feature as:

- **sector-aware symbolic reduction** for parity / `U(1)` / `S_z`, plus
- **basis-transform reduction** for SU(2) spin adaptation,

all while keeping `MomentProblem` as the output contract.

That is the graceful path.

Anything that tries to force the entire problem through today’s `SignedPermutation` monoid MVP is the wrong abstraction wearing a fake mustache.
