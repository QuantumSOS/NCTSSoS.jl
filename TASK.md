# Current Task — SympleQ → `NCTSSoS.jl` symmetry bridge

## One-sentence goal

Add a *find*-side SympleQ engine to `NCTSSoS.jl` that takes a Pauli Hamiltonian,
discovers Clifford symmetries via graph automorphism on the symplectic tableau,
and feeds them into the existing `SymmetrySpec` / `SolverConfig.symmetry`
pipeline so the moment / SOS hierarchy block-diagonalizes automatically.

## Source of truth

The feature, its dependency graph, and its known holes are audited in:

- `~/MyBrain/YushengBrain/topics/sympleq-implementation-gaps/sympleq-implementation-gaps.typ`
  (rendered: `…/sympleq-implementation-gaps.pdf`)

Treat the audit's eleven-step decomposition and its **find/exploit split** as
the canonical task structure. This plan mirrors that numbering so every code
deliverable maps back to a step the audit either tagged ✓ (specified),
△ (fillable from cited refs), or ✗ (blocked on unposted companion paper).

## Verdict on what we can and cannot implement today

- *Find* half (Steps 1–7): ~60% directly from the PRL + supplemental + public
  refs. Two real gaps: cycle enumeration in Step 4 ("minimal" is ambiguous),
  and the Step 7 phase-vector recovery for the Clifford lift $\hat S$.
- *Exploit* half (Steps 8–11): blocked. Step 8 (qubit-cost minimisation
  $\hat B$) is described in one sentence and points at an unposted
  manuscript. Steps 9 and 11 inherit that block.

**Therefore: ship the find side now. The exploit side is _not_ on the critical
path for the bridge** — `NCTSSoS.jl` only needs a Clifford generator
$\hat S$ rewritten as a signed permutation on the Pauli word basis. That is
what the existing `SymmetrySpec(generators::Vector{SignedPermutation})` slot
already accepts; we just don't yet hook it up for the Pauli (twisted-group)
algebra and we don't yet have automatic discovery.

## What already exists in this repo (do **not** rewrite)

- `SymmetrySpec` — user-facing finite-action spec.
  - Built from `SignedPermutation` generators (and a separate fermionic path).
  - Exported from `src/NCTSSoS.jl`.
- `SolverConfig(; symmetry::Union{Nothing,SymmetrySpec} = nothing, …)` — the
  wiring into the moment / SOS pipeline.
- `NCWordSignedPermutationAction <: SymbolicWedderburn.BySignedPermutations`
  in `src/optimization/symmetry.jl` — the `SymbolicWedderburn` adapter that
  lifts a `SignedPermutation` on letters to a signed permutation on words
  for `MonoidAlgebra` problems. **This file is ~1.6k lines; treat it as
  load-bearing.**
- `NCFermionicModePermutationAction` — same adapter for fermionic mode
  permutations.

The SympleQ bridge does **not** touch this seam. It feeds it.

## Strict scope of this task

In scope:

- A new SympleQ subsystem under `src/sympleq/` (new directory) that operates
  on `TwistedGroupAlgebra` (Pauli) inputs.
- *Find* half only: Steps 1–7 from the audit.
- The Step 7 phase recovery delivered as **best-effort + clearly labelled**;
  when it fails we still return $\hat S$ as a symplectic generator and warn
  that the lift is unverified.
- Dense / no-sparsity Pauli problems first
  (`cs_algo = NoElimination()`, `ts_algo = NoElimination()`).
- One acceptance example: a small (≤ 6-qubit) Pauli Hamiltonian where
  SympleQ finds at least one non-trivial Clifford symmetry and the moment
  matrix visibly block-decomposes vs. the no-symmetry baseline.

Out of scope (deferred):

- Step 8 qubit-cost minimisation ($\hat B$). Wait for
  the companion paper / SympleQ code release. The audit's recommended action
  list ranks "clone SympleQ as the spec; email the authors" as the right
  unblocker — that work is **task-orthogonal** to the bridge.
- Steps 9–11 (tensor factorisation + per-block effective Hamiltonians).
- Term sparsity composition (`ts_algo`), multi-clique composition.
- Fermionic, bosonic, state-poly, trace-poly support inside SympleQ. The
  existing `NCFermionicModePermutationAction` covers fermionic symmetries
  the *user* knows; SympleQ-style automatic discovery for fermions is
  separate and not in this task.

## Phased implementation plan

Each phase = one PR, one new test file, one clear acceptance criterion.

### Phase 0 — Scaffolding

- Create `src/sympleq/` and stub `src/sympleq/SympleQ.jl` as a sub-include
  from `src/NCTSSoS.jl`.
- No new public exports yet.
- **Done when:** `using NCTSSoS` still works; `make test` is green.

### Phase 1 — Steps 1–3: tableau, symplectic product, anticommutation graph

Audit status: all ✓.

- `src/sympleq/tableau.jl`
  - `struct SymplecticTableau` holding $(\vec c, \underline p, \vec\eta)$
    per the PRL Eq. (1). $\vec p_i \in \mathbb F_2^{2n}$, $\eta_i \in \mathbb Z_4$.
  - Constructor `SymplecticTableau(H::Polynomial{<:TwistedGroupAlgebra,…})`.
- `src/sympleq/graph.jl`
  - Symplectic product matrix per PRL Eq. (3), $O(M^2)$ for $M$ Pauli terms.
  - Vertex colouring by coefficient $\vec c$; edges where
    $\langle \vec p_i, \vec p_j\rangle = 1 \pmod 2$.
- **Done when:** unit tests in `test/sympleq/tableau.jl` and
  `test/sympleq/graph.jl` round-trip a 2-qubit Pauli polynomial to the
  expected tableau and adjacency matrix.

### Phase 2 — Step 4: cycle enumeration (the first real gap)

Audit status: △/✗. "Minimal" is ambiguous; matroid circuits ≠ cycle basis.

- Pick **one** concrete strategy first and document the choice loudly:
  - default: a **cycle basis** of the GF(2) null space of $\underline p^\top$
    (size $M - \mathrm{rank}(\underline p)$, polynomial cost).
  - leave a `cycle_strategy::Symbol` knob (`:cycle_basis`, future
    `:matroid_circuits`) so we can swap when the SympleQ source clarifies.
- Add auxiliary green vertices, connect each to every $\vec p_i$ in its
  cycle.
- **Done when:** `test/sympleq/cycles.jl` documents the chosen
  convention with an explicit small example; tests **must** include a
  comment pointing at this gap in the audit.

### Phase 3 — Step 5: graph automorphism solver

Audit status: ✓.

- Wrap one of nauty / bliss / `Graphs.jl`'s automorphism backend.
  Recommendation: depend on
  [`Bliss.jl`](https://github.com/JuliaPlots/Graphs.jl) wrapper or
  call `nauty` via a thin shell-out if no clean Julia binding exists. Make
  this an extension (`ext/`) so non-SympleQ users do not pull the dep.
- Input: colored graph from Phase 1 + Phase 2. Output: a generating set of
  the colour-preserving automorphism group as `Vector{Permutation}` on the
  Pauli-vertex indices.
- **Done when:** on a hand-checked 4-Pauli toy graph the returned generators
  match the known stabiliser by direct enumeration in a unit test.

### Phase 4 — Step 6: symplectic $\underline S$ from permutation $\underline\Pi$

Audit status: △. Standard "synthesis from Paulis" problem.

- Implement
  $\underline S = \underline p_{\text{basis}}^{-1} (\underline\Pi\,\underline p)_{\text{basis}}$
  with the rank-deficient corner case patched per
  Rengaswamy–Calderbank–Newman 2018 §III (extend $\underline p$ with
  auxiliary independent symplectic vectors before inversion).
- **Done when:** synthesised $\underline S$ satisfies
  $\underline S^\top\,\underline\Omega\,\underline S = \underline\Omega$ on
  a fuzz battery of random tableaux, and
  $\underline\Pi\,\underline p = \underline p\,\underline S \pmod 2$ on the
  rows present in the basis.

### Phase 5 — Step 7: phase vector $\vec\phi_S$ (the second real gap)

Audit status: ✗. Companion-only.

- Implement Hostens–De Hænen–De Moor 2005 §3–4 phase update rules for
  Clifford generators (SUM, H, Phase, Pauli) operating on $\vec\eta$.
- Solve the resulting $\mathbb Z_4$ linear system over a Clifford
  decomposition of $\underline S$ to recover $\vec\phi_S$ such that
  $\underline\Pi\vec\eta \equiv \vec\eta + (\text{linear part of }\underline S\,\vec\eta) + \vec\phi_S \pmod 4$.
- If the solver fails for a given $\underline S$ (and it will, until we have
  the SympleQ source as ground truth), **return $\hat S$ as a
  symplectic-only generator with `phase_verified = false`** and surface a
  warning. Do not silently drop the symmetry.
- **Done when:** for the canonical SympleQ example in the PRL the
  recovered phase reproduces the published Clifford; for random tableaux
  we either succeed or fail loudly with `phase_verified = false`.

### Phase 6 — Bridge: Clifford $\hat S$ → `SignedPermutation` on Pauli words

The integration point. The existing `SymmetrySpec` accepts
`Vector{SignedPermutation}` over letter indices; a Clifford symmetry
$\hat S$ acts on the *Pauli letter basis* $\{I, X_i, Y_i, Z_i\}$ as
exactly a signed permutation (signs from the phase tracker). So:

- Add `NCPauliSignedPermutationAction <: SymbolicWedderburn.BySignedPermutations`
  in `src/optimization/symmetry.jl` (next to the existing
  `NCWordSignedPermutationAction`) that knows the `TwistedGroupAlgebra`
  word structure.
- Add a converter
  `clifford_to_signed_permutation(S::SymplecticMatrix, φ::PhaseVector,
   registry::VariableRegistry{<:TwistedGroupAlgebra}) ::
   SignedPermutation`.
- Add a one-shot helper:
  ```julia
  sympleq_symmetry_spec(H::Polynomial{<:TwistedGroupAlgebra,…};
                        cycle_strategy = :cycle_basis,
                        ga_backend     = :bliss) ::
      SymmetrySpec
  ```
  that runs Phases 1–5 and returns a ready-to-feed `SymmetrySpec`.
- **Done when:** for a 4-qubit Ising chain or transverse-field XY model,
  `solve_sdp(polyopt(H; symmetry = sympleq_symmetry_spec(H)), COSMO)` runs
  end-to-end and reports strictly smaller PSD blocks than the no-symmetry
  baseline at the same relaxation order. Reviewed fixture lives in
  `test/data/expectations/sympleq_bridge.json`.

### Phase 7 — Tests and a literate example

- `test/sympleq/` — unit tests for Phases 1–5, isolated from any solver.
- `test/relaxations/sympleq_bridge.jl` — COSMO-stable end-to-end test of
  Phase 6 against the reviewed fixture above.
- `docs/src/examples/literate/sympleq_symmetry.jl` — short narrative example
  ("from a Pauli Hamiltonian to a block-decomposed SOS"), routed via the
  Pauli / ground-state literate cluster per the repo's docs routing rules.
  Regenerate with `make examples`.

### Phase 8 (deferred) — Exploit side

Do not start until **either** the companion paper @NationCompanion2026
posts **or** the SympleQ source @SympleQcode goes public. When that happens,
the audit's Steps 8–11 become tractable; the work plugs into a per-block
relaxation loop on top of Phase 6 and does not invalidate any of the find
side. Track the unblock by re-reading the audit and updating this section.

## Recommended unblocking actions (not on the critical path)

From the audit, ranked by leverage and runnable in parallel with Phases 0–7:

1. Clone the SympleQ source @SympleQcode when it lands and read three
   modules: graph construction (Step 4 cycle convention), phase tracking
   (Step 7), Jordan-symplectic / qubit-cost (Step 8). Promote findings into
   this plan.
2. Email Charlie Nation `c.nation2@exeter.ac.uk` and Luca Dellantonio
   `l.dellantonio@exeter.ac.uk` for a draft of
   @NationSymmetryDiagonalization. They cite it four times as load-bearing.
3. Read Rengaswamy–Calderbank–Newman 2018 and Hostens et al. 2005 §3–4
   *before* Phases 4 and 5. The code then reads as confirmation rather
   than reverse engineering.

## Definition of done for the whole task

- Phases 0–7 merged with green CI on COSMO.
- One reviewed expectation fixture protects the bridge.
- One literate example shows the user-facing workflow.
- Audit gaps (Step 4 cycle convention, Step 7 phase recovery) are called
  out in code comments **and** in the literate example's narrative, so a
  future reader is not surprised by the limitations.
