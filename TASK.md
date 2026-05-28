# Active Task: Clifford Symmetry for Pauli Algebra (Symmetry-Adapted Basis)

> Sibling worktree of `feat/fermionic-symmetry-adapted-basis`. Branched from it at commit `2fd6397` so the existing SAB scaffolding (`src/optimization/symmetry.jl`, `SymbolicWedderburn` integration, `fermionic_irreps.jl` patterns) is available as a template — **not** as a constraint.

## Goal
Extend the symmetry-adapted basis (SAB) machinery in NCTSSoS to support the **Clifford group action on the Pauli algebra** (`TwistedGroupAlgebra` / `PauliAlgebra`), so that Pauli-algebra moment/SOS relaxations can be block-diagonalized by Clifford symmetries supplied by the user (or by a small registered subgroup).

## Why this is non-trivial
The current SAB MVP in `src/optimization/symmetry.jl` is restricted to `SignedPermutation`s of registry indices. Clifford symmetries on Pauli words are strictly more general:

- A single-qubit Clifford permutes `{X, Y, Z}` up to a sign, e.g. the Hadamard exchanges X↔Z, the S gate sends X→Y, Y→−X.
- Multi-qubit Cliffords (CNOT, SWAP, etc.) entangle the action across sites — `IX → XX`, `IZ → ZZ`, etc. — so the action is **not** a site-local signed permutation of registry indices.
- The action on a Pauli word `P = c · σ_{i₁} ⊗ … ⊗ σ_{iₙ}` is: `U P U† = c · ϕ(U,P) · σ'_{j₁} ⊗ … ⊗ σ'_{jₙ}` with phase `ϕ ∈ {±1}` (Clifford preserves Hermitian Paulis up to sign).
- This is a **monomial action** on the Pauli basis (permutation + ±1 phase), which is exactly the regime `SymbolicWedderburn` handles via signed permutations on the **monomial basis** — but the underlying group elements are *Clifford operators*, not index permutations.

So the right abstraction is: a group homomorphism `ρ: G → SignedPerm(PauliWords)` induced by conjugation, where `G` is a finite subgroup of the n-qubit Clifford group (typically supplied by the user as a list of generators).

## Concrete deliverables
1. **Type:** `CliffordSymmetry` (or similar) representing a single n-qubit Clifford element acting on Pauli words by conjugation. Stored compactly — e.g. via its `(2n × 2n)` symplectic matrix over `𝔽₂` plus an `n`-vector of sign bits (the standard tableau representation). Action on a `NormalMonomial{PauliAlgebra}` returns `(sign, NormalMonomial)`.
2. **Group plumbing:** make a finite `CliffordSymmetryGroup` satisfy `GroupsCore.Group` so `SymbolicWedderburn` can iterate it. Closure check + a `closure(generators)` helper.
3. **SAB integration:** generalize the existing `symmetry.jl` reduction so it accepts a `CliffordSymmetryGroup` acting on a basis of `NormalMonomial{PauliAlgebra}` (currently it only consumes `SignedPermutation`s of registry indices). Reuse `SymbolicWedderburn.symmetry_adapted_basis` with the appropriate action.
4. **User-facing API:** ability to attach a `CliffordSymmetryGroup` to a Pauli `polyopt` problem the same way fermionic symmetries are attached.
5. **Tests:** at minimum, a 2-qubit and 3-qubit example where the SAB-reduced moment matrix matches the unreduced moment matrix's optimum to solver tolerance under COSMO. Pick at least one case where the block decomposition is non-trivial (multiple irreps).
6. **Docs:** one Literate example in `docs/src/examples/literate/` mirroring the style of `h2_fermionic_symmetry.jl` but for a small Pauli Hamiltonian with an obvious Clifford symmetry (e.g. 2-qubit Heisenberg with SWAP symmetry, or a TFIM segment with reflection × `Z_2` spin-flip).

## Out of scope (for this worktree)
- The full Clifford group as an automatic discovery problem. The user supplies generators.
- Non-Clifford unitary symmetries.
- Continuous symmetry groups.
- Re-touching fermionic SAB code beyond what is needed to factor shared utilities.

## Suggested order of work
1. Read `src/optimization/symmetry.jl` end-to-end. Understand the MVP path and the `SymbolicWedderburn` action contract it satisfies.
2. Read `src/simplification/pauli.jl` and the `TwistedGroupAlgebra` portions of `src/types/algebra.jl` to nail down how a Pauli word is stored and rewritten.
3. Prototype `CliffordSymmetry` and the conjugation action on a single Pauli word in isolation (unit tests against hand-computed examples first — Hadamard on X, S on Y, CNOT on `IX`, etc.).
4. Build the group closure + `GroupsCore.Group` interface.
5. Wire into `symmetry_adapted_basis` with a `SymbolicWedderburn.ByPermutations`-style action on the monomial basis.
6. Glue into `polyopt` for Pauli problems; add an end-to-end test.
7. Write the Literate doc; run `make examples` then `make servedocs`.

## Branch / worktree info
- Branch: `feat/clifford-pauli-symmetry-adapted-basis`
- Worktree path: `/Users/exaclior/QuantumSOS/NCTSSoS.jl-clifford-pauli-sab`
- Parent: `feat/fermionic-symmetry-adapted-basis` @ `2fd6397`

## House rules (from repo `AGENTS.md`)
- Surgical changes; add regression tests for bug fixes.
- Type-param order: algebra type `A` first.
- Tests stay deterministic and COSMO-stable.
- For solver-backed numeric expectations, prefer `test/data/expectations/*.json` fixtures over inline magic numbers.
- After any docs change: `make examples` then `make servedocs` before handoff.

## Open questions to resolve before serious coding
- **Action target:** should `CliffordSymmetry` act on the registry-level basis (Pauli letters per site) or on assembled `NormalMonomial`s? The fermionic path acts at the registry level; for general multi-qubit Cliffords we likely need monomial-level action. Confirm with a 2-qubit CNOT example.
- **Storage:** tableau (`F₂` symplectic + signs) vs. dict-of-images. Tableau is canonical and `O(n²)`; dict is easier to debug. Start with dict, switch to tableau if perf bites.
- **Group input format:** user-friendly generator list of named gates (`:H`, `:S`, `:CNOT`, …) vs. raw tableaus. Pick one for v1, document, and stop.
