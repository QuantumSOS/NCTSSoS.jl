# Plan: Enhance Pauli Clifford Symmetry Documentation

**Branch:** `docs/enhance-pauli-clifford-symmetry`
**Date:** 2026-06-09

## Context

The current Pauli Clifford Symmetry docs consist of:
- **Example page** (`docs/src/examples/literate/pauli_clifford_symmetry.jl`): demonstrates manual `CliffordSymmetry(:SWAP,…)` and automatic `sympleq_symmetry_spec()` on a 2-site and 4-site Heisenberg model. Cute, but small — no timing, no performance comparison, no "feel the difference" moment.
- **Manual pages** (`manual/symmetry_adapted_basis.md`, `manual/extending_symmetry.md`): describe the SAB pipeline, what `SymbolicWedderburn` does, the MVP scope, and how fermionic extension would work. **Neither page explains *what algorithm* SympleQ uses to find Clifford symmetries.** A reader walks away knowing NCTSSoS has the feature but not *why* it works.

The implementation lives in `src/sympleq/` (7 files, ~1210 lines) and faithfully implements the find-side of Nation et al. (arXiv:2605.18966): symplectic tableau → coloured anticommutation graph → graph automorphisms → symplectic synthesis → phase verification → bridge to `SymmetrySpec`.

## Two deliverables

### Deliverable 1: New manual page — "Clifford Symmetry Detection"

**File:** `docs/src/manual/clifford_symmetry_detection.md`
**Location in nav:** Manual section, after "Extending Symmetry Support"

**Goal:** A reader who finishes this page understands *what algorithm* `sympleq_symmetry_spec` runs, *why* graph automorphisms correspond to Clifford symmetries, and where the theoretical limits are. No code required to follow it — that's the example page's job.

**Outline:**

1. **The problem statement** (~3 paragraphs)
   - Finding Hamiltonian symmetries by inspection works for textbook cases. For non-local Clifford symmetries on ≥20 qubits, it doesn't.
   - Prior methods only find *Pauli* symmetries (operators that commute term-by-term). Cliffords are strictly more general — they can permute and mix Pauli terms while preserving the Hamiltonian as a whole.
   - This page explains the SympleQ algorithm (Nation et al., 2026) that NCTSSoS uses to find Clifford symmetries automatically.

2. **Background: the symplectic representation** (~4 paragraphs + 1 equation block)
   - Each Pauli string on $n$ qubits ↔ a binary vector $\vec{p} \in \mathbb{Z}_2^{2n}$ plus a $\mathbb{Z}_4$ phase $\eta$.
   - A Clifford unitary ↔ a symplectic matrix $\underline{S} \in \text{Sp}(2n, \mathbb{F}_2)$ plus a phase vector $\vec{\phi}_S$.
   - Conjugation by a Clifford maps $\vec{p} \mapsto \vec{p}\,\underline{S}$, which preserves the symplectic product $\langle \vec{p}_i, \vec{p}_j \rangle = \vec{p}_i\,\Omega\,\vec{p}_j^T \mod 2$.
   - Reference to Aaronson-Gottesman (2004) and Hostens-Dehaene-De Moor (2005).

3. **The key insight: Clifford symmetries = graph automorphisms** (~3 paragraphs + 1 figure description)
   - Build a coloured graph from the Hamiltonian:
     - **Vertices:** one per Pauli term, coloured by coefficient.
     - **Edges:** connect anticommuting pairs ($\langle \vec{p}_i, \vec{p}_j \rangle = 1$).
     - **Auxiliary vertices:** one per GF(2)-linear dependency among the Pauli vectors (cycle vertices).
   - A colour-preserving graph automorphism $\Pi$ permutes Pauli terms while preserving (a) coefficients, (b) commutation relations, and (c) linear dependencies. These are exactly the three invariants of Clifford conjugation.
   - **Therefore:** every graph automorphism corresponds to a Clifford symmetry $\hat{S}$, and vice versa.

4. **The seven-step find algorithm** (numbered list matching `src/sympleq/` structure)
   - Step 1: Build the symplectic tableau from the Pauli polynomial. → `tableau.jl`
   - Step 2: Compute the symplectic product matrix ($O(M^2)$). → `tableau.jl`
   - Step 3: Colour vertices by coefficient, add anticommutation edges. → `graph.jl`
   - Step 4: Find GF(2) null space of the Pauli matrix; add auxiliary cycle vertices. → `cycles.jl`
   - Step 5: Run a graph-automorphism solver (nauty/bliss/igraph). → `automorphism.jl`
   - Step 6: Synthesize the symplectic matrix $\underline{S}$ from the permutation $\Pi$ using a basis of independent Pauli vectors. → `symplectic.jl`
   - Step 7: Solve for the phase vector $\vec{\phi}_S$ and verify it against the tableau. → `phase.jl`

5. **From SympleQ generators to the SAB pipeline** (~2 paragraphs)
   - The bridge (`bridge.jl`) converts each verified `SympleQGenerator` into a `CliffordSymmetry` and wraps the collection in a `SymmetrySpec`.
   - From there, the existing SAB machinery (group enumeration → `SymbolicWedderburn` decomposition → orbit reduction → smaller `MomentProblem`) takes over. See [Symmetry-Adapted Basis](@ref symmetry-adapted-basis).

6. **Scaling and limitations** (~3 paragraphs)
   - Time dominated by the symplectic product matrix: $O(M^2)$ where $M$ is the number of Pauli terms. In practice, systems up to $n = 1000$ qubits are tractable for the *find* step.
   - The *exploit* side (qubit-cost minimization, block decomposition into effective sub-Hamiltonians) from the SympleQ paper is **not implemented** — it requires algorithms from an unpublished companion paper. NCTSSoS takes the found Clifford generators and feeds them into the SAB pipeline instead.
   - Current scope: Hermitian Pauli Hamiltonians only (real ±1 phases on Clifford images). Generators with $\pm i$ phases are rejected with a clear error.

7. **References** — cite Nation et al. 2026, Aaronson-Gottesman 2004, Hostens et al. 2005, Rengaswamy et al. 2018, McKay-Piperno 2014 (nauty), Babai 2016.

8. **See also** links — to the SAB manual page, the extending-symmetry page, the Pauli Clifford example, and the API reference for `sympleq_symmetry_spec`, `CliffordSymmetry`, `SymmetrySpec`.

**Source material:**
- Implementation: `src/sympleq/*.jl` (7 files, the code comments are already step-labelled)
- Theory: `/Users/exaclior/MyBrain/YushengBrain/topics/sympleq-implementation-gaps/sympleq-implementation-gaps.typ` (the 11-step audit)
- Paper: `/Users/exaclior/MyBrain/YushengBrain/references/sources/2605.18966/PRL.tex` and `NOTES.md`
- Existing note: `/Users/exaclior/MyBrain/YushengBrain/topics/nctssos-symmetry-landscape/` (the three-community comparison)

---

### Deliverable 2: Enhanced example — larger system with timing

**File:** `docs/src/examples/literate/pauli_clifford_symmetry.jl` (modify in place)

**Goal:** Keep the existing 2-site and 4-site sections. Append a new section that:
1. **Builds a physically interesting Hamiltonian at $N = 8$ or $N = 10$ qubits** — a periodic Heisenberg chain or XXZ ladder. Large enough that the timing difference is obvious; small enough that COSMO/Mosek can handle the dense baseline within docs-regen time.
2. **Times the dense baseline vs. symmetry-reduced solve** using `@elapsed` (not `@benchmark` — we don't want BenchmarkTools as a docs dependency). Show wall-clock seconds for both.
3. **Prints a comparison table** summarizing: system size, PSD block sizes (dense vs. reduced), number of moment variables (dense vs. reduced), solve time (dense vs. reduced), group order, and the number of generators found.
4. **Makes the payoff visceral** — the reader should see something like "dense: 3.2s with a single 40×40 block; symmetric: 0.4s with 8 blocks of size ≤5".

**Sizing decision:**
- $N = 8$ periodic Heisenberg chain: 24 Pauli terms, order-1 basis has $3N + 1 = 25$ monomials. The dense moment matrix is $25 \times 25$. Symmetry group should be nontrivial (translation + reflection + spin-flip on a ring = dihedral × Z₂). This is the sweet spot — large enough to show clear speedup, small enough to solve in <30s for docs regen.
- If $N = 8$ is too fast to show timing differences, bump to $N = 10$.

**Structure of the new section:**

```
## Performance comparison: 8-site Heisenberg ring

### Build the Hamiltonian
### Dense baseline (timed)
### SympleQ detection + symmetry-reduced solve (timed)
### Comparison table
### Discussion: what makes the difference
```

**What NOT to change:**
- The existing 2-site manual SWAP and 4-site auto-detection sections stay intact. They serve a pedagogical purpose (introducing the API incrementally).
- The summary table at the bottom gets updated to include the new 8-site row.

---

## Nav / make.jl changes

1. Add the new manual page to `docs/make.jl`:
   ```julia
   "Manual" => Any[
       ...
       "Extending Symmetry Support" => "manual/extending_symmetry.md",
       "Clifford Symmetry Detection" => "manual/clifford_symmetry_detection.md",  # NEW
       "SDP Relaxation" => "manual/sdp_relaxation.md",
       ...
   ]
   ```

2. Add bib entries to `docs/src/refs.bib` for Nation et al. 2026, Aaronson-Gottesman 2004, Hostens et al. 2005, Rengaswamy et al. 2018, McKay-Piperno 2014 (nauty), Babai 2016 — if not already present.

---

## Verification steps

1. `make examples` — regenerate all Literate pages (the timing sections will execute).
2. `make servedocs` — preview the new manual page and the updated example.
3. Confirm the new manual page renders correctly: equations, step list, cross-links.
4. Confirm the timing comparison in the example shows a meaningful speedup.
5. `make test` — ensure no existing tests break.

---

## Ordering

1. Write `docs/src/manual/clifford_symmetry_detection.md` (Deliverable 1)
2. Update `docs/make.jl` nav + add bib entries
3. Update `docs/src/examples/literate/pauli_clifford_symmetry.jl` (Deliverable 2)
4. `make examples` + `make servedocs` + review
5. `make test` for safety
