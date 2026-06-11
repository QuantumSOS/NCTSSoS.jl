# Source Material for Clifford Symmetry Documentation

## Implementation code (in-repo)

The SympleQ find-side lives in `src/sympleq/` (7 files, ~1210 lines total). Code comments are already step-labelled to match the algorithm:

| File | SympleQ step | What it does |
|------|-------------|--------------|
| `tableau.jl` | Steps 1–2 | Pauli polynomial → binary symplectic tableau + product matrix |
| `graph.jl` | Step 3 | Coloured anticommutation graph construction |
| `cycles.jl` | Step 4 | GF(2) null-space cycle-basis augmentation |
| `automorphism.jl` | Step 5 | Colour-preserving graph automorphisms (igraph/nauty) |
| `symplectic.jl` | Step 6 | Synthesize symplectic $S$ from term permutation $\Pi$ |
| `phase.jl` | Step 7 | Best-effort $\mathbb{Z}_4$ phase verification |
| `bridge.jl` | — | Convert `SympleQGenerator` → `CliffordSymmetry` / `SymmetrySpec` |

The SAB pipeline: `src/optimization/symmetry.jl`, `src/optimization/interface.jl`.

## Theory notes (external — YushengBrain vault)

- **SympleQ algorithm audit (11 steps):** `/Users/exaclior/MyBrain/YushengBrain/topics/sympleq-implementation-gaps/sympleq-implementation-gaps.typ`
  - Score-card: find side 60% implementable from PRL alone; exploit side blocked by unpublished companion.
  - Gap analysis: Step 4 (cycle enumeration ambiguity), Step 7 (phase vector recovery), Step 8 (qubit-cost minimization — not implemented, not needed for NCTSSoS bridge).

- **Symmetry landscape comparison:** `/Users/exaclior/MyBrain/YushengBrain/topics/nctssos-symmetry-landscape/nctssos-symmetry-landscape.typ`
  - Three-community comparison: polyopt, tensor networks, NCTSSoS. Extracts three-axis pattern (group type × action target × SDP entry channel).

- **SAB tutorial:** `/Users/exaclior/MyBrain/YushengBrain/topics/symmetry-adapted-basis/symmetry-adapted-basis.typ`
  - Zero-assumption tutorial on SAB for NCPO. Includes CHSH worked example, implementation-notes appendix.

## Paper source

- **Nation et al. (arXiv:2605.18966):** `/Users/exaclior/MyBrain/YushengBrain/references/sources/2605.18966/`
  - `PRL.tex` — main manuscript
  - `NOTES.md` — reading notes with 78 citations grouped by role, recommended reading order
  - `cited.bib` — filtered bibliography (78 entries actually cited)
  - Key lines: 187 (motivation), 381–396 (graph construction + GA = Clifford symmetry theorem), 403 (symplectic synthesis), 564 (scaling benchmarks up to 1000 qubits)

## Bibliography entries needed in `docs/src/refs.bib`

Check whether these already exist; add missing ones:
- Nation et al. 2026 (arXiv:2605.18966) — the SympleQ PRL
- Aaronson & Gottesman 2004 — symplectic tableau
- Hostens, Dehaene & De Moor 2005 — stabilizer formalism + phase update rules
- Rengaswamy, Calderbank & Newman 2018 — symplectic synthesis from Paulis
- McKay & Piperno 2014 — nauty II (practical graph automorphism solver)
- Babai 2016 — quasi-polynomial graph isomorphism
