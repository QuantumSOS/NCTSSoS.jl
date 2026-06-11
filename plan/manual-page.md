# Deliverable 1: New Manual Page — "Clifford Symmetry Detection"

**File:** `docs/src/manual/clifford_symmetry_detection.md`
**Location in nav:** Manual section, after "Extending Symmetry Support"

**Goal:** A reader who finishes this page understands *what algorithm* `sympleq_symmetry_spec` runs, *why* graph automorphisms correspond to Clifford symmetries, and where the theoretical limits are. No code required to follow it — that's the example page's job.

## Outline

### 1. The problem statement (~3 paragraphs)
- Finding Hamiltonian symmetries by inspection works for textbook cases. For non-local Clifford symmetries on ≥20 qubits, it doesn't.
- Prior methods only find *Pauli* symmetries (operators that commute term-by-term). Cliffords are strictly more general — they can permute and mix Pauli terms while preserving the Hamiltonian as a whole.
- This page explains the SympleQ algorithm (Nation et al., 2026) that NCTSSoS uses to find Clifford symmetries automatically.

### 2. Background: the symplectic representation (~4 paragraphs + 1 equation block)
- Each Pauli string on $n$ qubits ↔ a binary vector $\vec{p} \in \mathbb{Z}_2^{2n}$ plus a $\mathbb{Z}_4$ phase $\eta$.
- A Clifford unitary ↔ a symplectic matrix $\underline{S} \in \text{Sp}(2n, \mathbb{F}_2)$ plus a phase vector $\vec{\phi}_S$.
- Conjugation by a Clifford maps $\vec{p} \mapsto \vec{p}\,\underline{S}$, which preserves the symplectic product $\langle \vec{p}_i, \vec{p}_j \rangle = \vec{p}_i\,\Omega\,\vec{p}_j^T \mod 2$.
- Reference to Aaronson-Gottesman (2004) and Hostens-Dehaene-De Moor (2005).

### 3. The key insight: Clifford symmetries = graph automorphisms (~3 paragraphs + 1 figure description)
- Build a coloured graph from the Hamiltonian:
  - **Vertices:** one per Pauli term, coloured by coefficient.
  - **Edges:** connect anticommuting pairs ($\langle \vec{p}_i, \vec{p}_j \rangle = 1$).
  - **Auxiliary vertices:** one per GF(2)-linear dependency among the Pauli vectors (cycle vertices).
- A colour-preserving graph automorphism $\Pi$ permutes Pauli terms while preserving (a) coefficients, (b) commutation relations, and (c) linear dependencies. These are exactly the three invariants of Clifford conjugation.
- **Therefore:** every graph automorphism corresponds to a Clifford symmetry $\hat{S}$, and vice versa.

### 4. The seven-step find algorithm (numbered list matching `src/sympleq/` structure)
- Step 1: Build the symplectic tableau from the Pauli polynomial. → `tableau.jl`
- Step 2: Compute the symplectic product matrix ($O(M^2)$). → `tableau.jl`
- Step 3: Colour vertices by coefficient, add anticommutation edges. → `graph.jl`
- Step 4: Find GF(2) null space of the Pauli matrix; add auxiliary cycle vertices. → `cycles.jl`
- Step 5: Run a graph-automorphism solver (nauty/bliss/igraph). → `automorphism.jl`
- Step 6: Synthesize the symplectic matrix $\underline{S}$ from the permutation $\Pi$ using a basis of independent Pauli vectors. → `symplectic.jl`
- Step 7: Solve for the phase vector $\vec{\phi}_S$ and verify it against the tableau. → `phase.jl`

### 5. From SympleQ generators to the SAB pipeline (~2 paragraphs)
- The bridge (`bridge.jl`) converts each verified `SympleQGenerator` into a `CliffordSymmetry` and wraps the collection in a `SymmetrySpec`.
- From there, the existing SAB machinery (group enumeration → `SymbolicWedderburn` decomposition → orbit reduction → smaller `MomentProblem`) takes over. See [Symmetry-Adapted Basis](@ref symmetry-adapted-basis).

### 6. Scaling and limitations (~3 paragraphs)
- Time dominated by the symplectic product matrix: $O(M^2)$ where $M$ is the number of Pauli terms. In practice, systems up to $n = 1000$ qubits are tractable for the *find* step.
- The *exploit* side (qubit-cost minimization, block decomposition into effective sub-Hamiltonians) from the SympleQ paper is **not implemented** — it requires algorithms from an unpublished companion paper. NCTSSoS takes the found Clifford generators and feeds them into the SAB pipeline instead.
- Current scope: Hermitian Pauli Hamiltonians only (real ±1 phases on Clifford images). Generators with $\pm i$ phases are rejected with a clear error.

### 7. References
Cite Nation et al. 2026, Aaronson-Gottesman 2004, Hostens et al. 2005, Rengaswamy et al. 2018, McKay-Piperno 2014 (nauty), Babai 2016.

### 8. See also
Links to the SAB manual page, the extending-symmetry page, the Pauli Clifford example, and the API reference for `sympleq_symmetry_spec`, `CliffordSymmetry`, `SymmetrySpec`.
