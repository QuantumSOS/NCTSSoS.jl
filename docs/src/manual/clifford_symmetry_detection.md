# [Clifford Symmetry Detection](@id clifford-symmetry-detection)

Finding Hamiltonian symmetries by inspection works for textbook cases: the
two-site Heisenberg model is obviously invariant under qubit swap. For
Hamiltonians on 20 or more qubits, or Hamiltonians whose symmetries mix Pauli
types across sites, manual enumeration is tedious, error-prone, and
occasionally impossible without systematic search.

Prior automated methods typically find only **Pauli symmetries** — Pauli
operators that commute with every term individually. Clifford symmetries are
strictly more general: a Clifford unitary can permute and mix Pauli terms
while preserving the Hamiltonian *as a whole*, even when no single Pauli
operator commutes with every term. These are exactly the symmetries that
[`sympleq_symmetry_spec`](@ref) finds.

Detection is one half of a two-part story. **Finding** a symmetry group (this
page) and **exploiting** it (the
[Symmetry-Adapted Basis](@ref symmetry-adapted-basis) pipeline, which
block-diagonalizes the moment matrix and orbit-reduces moment variables) are
separate steps; each is useless without the other. The reason detection
earns its keep even on "obviously symmetric" models is that **the symmetry
visible by inspection is usually a small subgroup of what is actually
there** — and the SDP reduction scales with the group you give it, not the
group that exists. Concretely, on the 2-site Heisenberg model:

| symmetry source | group order | PSD blocks (from one ``7\times 7``) | nontrivial moments |
|:----------------|:------------|:------------------------------------|:-------------------|
| by inspection (qubit SWAP) | 2 | `[4, 3]` | 9 |
| SympleQ detection | 24 | `[1, 2]` | 1 |
| union of both | 48 | `[1, 1, 1]` | 1 |

The order-24 group SympleQ finds is the global octahedral rotation group of
the Pauli axes — the discrete Clifford remnant of the model's spin-rotation
invariance. It contains symmetries that permute the terms of ``H`` without
fixing any of them (e.g. the cyclic axis rotation ``x \to y \to z``), which
is precisely the kind of symmetry that humans rarely write down and that no
Pauli-symmetry search can represent. The runnable side-by-side
comparison is at [Pauli Symmetry Reduction](@ref pauli-clifford-symmetry).

This page explains the SympleQ algorithm [nation2026sympleq](@cite) that
NCTSSoS uses to detect Clifford symmetries automatically.

## [Background: the symplectic representation](@id sympleq-background)

Every Pauli string on ``n`` qubits can be encoded as a binary vector plus a
phase. Write the Pauli operator ``P = i^{\eta}\, X^{x_1} Z^{z_1} \otimes
\cdots \otimes X^{x_n} Z^{z_n}`` as a pair: a vector
``\vec{p} = (x_1, \ldots, x_n, z_1, \ldots, z_n) \in \mathbb{F}_2^{2n}``
and a phase ``\eta \in \mathbb{Z}_4``. For example, the Pauli ``Y`` on a
single qubit is ``iXZ``, so ``\vec{p} = (1, 1)`` and ``\eta = 1``.

A Clifford unitary ``U`` acts on Pauli strings by conjugation:
``U P U^\dagger = i^{\eta'} P'``. In the binary picture, this conjugation
becomes a matrix multiply:

```math
\vec{p} \;\mapsto\; \vec{p}\,\underline{S},
```

where ``\underline{S} \in \mathrm{Sp}(2n, \mathbb{F}_2)`` is a binary
symplectic matrix satisfying ``\underline{S}\,\Omega\,\underline{S}^T = \Omega``
and

```math
\Omega = \begin{pmatrix} 0 & I_n \\ I_n & 0 \end{pmatrix}
```

is the standard symplectic form. The symplecticity condition is exactly the
statement that conjugation by ``U`` preserves commutation relations: two Pauli
strings commute if and only if their symplectic product
``\langle \vec{p}_i, \vec{p}_j \rangle = \vec{p}_i\,\Omega\,\vec{p}_j^T
\bmod 2`` vanishes [aaronson2004improved](@cite) [hostens2005stabilizer](@cite).

The Clifford also carries a ``\mathbb{Z}_4`` phase vector
``\vec{\phi}_S`` that tracks how the signs of individual Pauli generators
transform. Finding ``\underline{S}`` alone tells you *which* Pauli string each
operator maps to; ``\vec{\phi}_S`` tells you *with what sign*.

## [The key insight: graph automorphisms are Clifford symmetries](@id sympleq-insight)

The SympleQ algorithm encodes the structure of a Pauli Hamiltonian
``H = \sum_{i=1}^{M} c_i P_i`` as a coloured graph:

- **Vertices:** one per Pauli term ``P_i``, coloured by ``|c_i|`` (the
  absolute value of the coefficient).
- **Edges:** connect every anticommuting pair
  (``\langle \vec{p}_i, \vec{p}_j \rangle = 1``).
- **Linear relations:** checked algebraically during the search, without
  auxiliary vertices. The graph alone is not a complete encoding.

The search accepts term permutations ``\Pi`` that preserve three invariants
simultaneously:

1. **Coefficients** (same colour → same ``|c_i|``).
2. **Commutation relations** (edges are preserved → symplectic products are
   preserved).
3. **Linear dependencies** (the algebraic checks preserve the kernel of the
   transposed Pauli matrix, including relations longer than three terms).

These are exactly the three invariants of Clifford conjugation. Therefore,
**a graph automorphism that also preserves all linear relations admits a
binary symplectic lift** [nation2026sympleq](@cite). The phase equations must
still admit a solution for that lift to preserve the signed Hamiltonian.
The graph construction translates
the algebraic problem of finding Clifford symmetries into the combinatorial
problem of finding graph automorphisms — a problem with mature, fast solvers
[mckay2014practical](@cite) [babai2016graph](@cite).

The correspondence is not one-to-one in the reverse direction: distinct
Clifford symmetries can induce the *same* permutation of the Hamiltonian's
terms (they differ by a Clifford that fixes every term individually). The
synthesis step below picks one symplectic representative per automorphism,
which has a practical consequence discussed in
[Scaling and limitations](@ref sympleq-limits).

## [The seven-step find algorithm](@id sympleq-algorithm)

The implementation in `src/sympleq/` follows these steps, each corresponding
to a source file:

1. **Build the symplectic tableau** from the Pauli polynomial. Each term
   becomes a binary row ``\vec{p}_i`` plus a ``\mathbb{Z}_4`` phase ``\eta_i``;
   polynomial coefficients are stored alongside.
   → [`tableau.jl`]

2. **Compute the symplectic product matrix** (``O(M^2)``). Entry ``(i,j)`` is
   ``\langle \vec{p}_i, \vec{p}_j \rangle \bmod 2``.
   → [`tableau.jl`]

3. **Colour vertices by coefficient; add anticommutation edges.** Two terms
   that anticommute get an edge. Vertex colours are ``|c_i|``, so only terms
   with equal coefficient magnitudes can be exchanged.
   → [`graph.jl`]

4. **Prepare exact relation checks, without relation vertices.**
   Gaussian elimination in search order finds a basis of the
   full relation space. Each equation is checked on the image rows as soon
   as its last term is assigned. Zero images of a basis imply zero images
   of its entire span; since the term map is a permutation, this preserves
   the full kernel. **Basis relations are equations, not graph vertices:**
   requiring a particular basis to permute into itself can discard valid
   translations. The chosen auxiliary-vertex set is empty, so it cannot
   discard any symmetry. The full-kernel checks enforce linear structure
   without enumerating all circuits or basis sums, or enlarging the graph
   with a redundant invariant family of short relations.
   The legacy keyword `cycle_strategy=:cycle_basis` selects this corrected
   implementation; `pauli_cycle_basis` itself remains an algebraic utility.
   → [`cycles.jl`]

5. **Run a graph-automorphism solver.** The current implementation uses a
   deterministic in-tree backtracking search (no external binary dependencies).
   The solver returns colour-preserving automorphisms satisfying the full
   relation checks, restricted to Pauli vertices. Enumeration stops with a
   warning at `max_automorphisms`; a truncated search is not complete.
   → [`automorphism.jl`]

6. **Synthesize the binary symplectic matrix** ``\underline{S}`` from each
   term permutation ``\Pi``. Pick a maximal ``\mathbb{F}_2``-independent subset
   of the Pauli rows as a "domain basis", read off the target rows under
   ``\Pi``, and extend to a full symplectic change-of-basis using a Witt
   extension procedure [rengaswamy2018synthesis](@cite). The result satisfies
   ``P \cdot \underline{S} = P[\Pi]`` over ``\mathbb{F}_2``, where ``P`` is
   the Pauli matrix and ``P[\Pi]`` is its row-permuted version.
   → [`symplectic.jl`]

7. **Solve for the phase vector** ``\vec{\phi}_S`` **and verify.** The binary
   ``\underline{S}`` determines *which* Pauli string each operator maps to;
   the phase vector determines *the sign*. The implementation sets up a
   ``\mathbb{F}_2`` linear system whose solution space gives all sign
   assignments compatible with ``H \mapsto H``, then checks each candidate
   against the full Hamiltonian.
   → [`phase.jl`]

Generators that pass all checks are collected as `SympleQGenerator` structs.
Candidates whose phase lift would require ``\pm i`` (incompatible with
Hermitian Pauli images) are rejected with an explicit error.

## [From SympleQ generators to the SAB pipeline](@id sympleq-bridge)

The bridge (`bridge.jl`) converts each verified `SympleQGenerator` into a
[`CliffordSymmetry`](@ref) — the action type understood by NCTSSoS's Pauli
symmetry path — by reading off how ``\underline{S}`` and ``\vec{\phi}_S``
transform the ``3n`` single-site Pauli generators (``X_i``, ``Y_i``, ``Z_i``
for each qubit ``i``). The ``Y`` image is derived from the ``X`` and ``Z``
images using ``Y = iXZ``, so only ``2n`` symplectic lookups are needed.

When the letters and their candidate images form a small closed set (at most
four times the number of letters), signed permutations on that set give a
faithful action. Permutation-group membership tests remove redundant
generators without enumerating the Clifford group. If the set is too large
or not closed, all candidates are retained. Reduction preserves the generated
group; it does not promise a minimum generating set.

The collection of `CliffordSymmetry` generators is wrapped in a
[`SymmetrySpec`](@ref) and handed to the existing SAB machinery:

```
SympleQ generators
  ↓  bridge.jl
CliffordSymmetry generators
  ↓  SymmetrySpec(...)
Group enumeration → SymbolicWedderburn decomposition → orbit reduction → smaller MomentProblem
```

From there, the pipeline described in [Symmetry-Adapted Basis](@ref sab-pipeline) takes over: group
enumeration, invariance checking, basis closure, symmetry-adapted decomposition
via [`SymbolicWedderburn`](https://github.com/kalmarek/SymbolicWedderburn.jl),
orbit reduction, and construction of a reduced `MomentProblem` with smaller PSD
blocks.

## [Scaling and limitations](@id sympleq-limits)

**Scaling.** With dense binary rows, computing symplectic products costs
``O(M^2 n)`` for ``M`` terms on ``n`` qubits. Relation preprocessing uses
Gaussian elimination and polynomial storage; no auxiliary vertices or
exponential set of circuits is constructed. Automorphism search and the number of term
permutations can still be exponential. Each accepted permutation is currently
synthesized before the final generator reduction, so compact output does not
remove this enumeration cost. The SympleQ paper reports systems up to
``n = 1000`` [nation2026sympleq](@cite); this is not a scaling guarantee for
this in-tree backend.

**Exploit side not implemented.** The SympleQ paper also describes an
"exploit" side — qubit-cost minimization and block decomposition into
effective sub-Hamiltonians — that requires algorithms from an unpublished
companion paper. NCTSSoS does not implement the exploit side. Instead, the
found Clifford generators are fed into the existing SAB pipeline, which
provides a different (but complementary) size reduction via moment-matrix
block diagonalization.

**Term-fixing Cliffords are not enumerated.** SympleQ sees a candidate
symmetry only through the permutation it induces on the Hamiltonian's terms.
For each such permutation, step 6 synthesizes a *single* symplectic
representative, and step 7 enumerates its compatible sign assignments. The
phase solve therefore does recover sign-flip refinements (Pauli symmetries
that negate operators while fixing every term), but **symplectically distinct
Cliffords inducing the same term permutation are missed**. The canonical
example is qubit SWAP on the 2-site Heisenberg model: it fixes all three
terms ``\sigma_1^a\sigma_2^a`` individually, so at the graph level it is
indistinguishable from the identity, and the detected group has order 24
rather than the full 48. The practical recipe is to **union the detected
generators with the structural generators you know**:

```julia
auto = sympleq_symmetry_spec(H)
spec = SymmetrySpec(auto.clifford_generators..., CliffordSymmetry(:SWAP, 1, 2))
```

Every manual generator is still verified for invariance before use, so a
wrong guess errors out instead of corrupting the relaxation. The
[Pauli Symmetry Reduction](@ref pauli-clifford-symmetry) example shows this
union recovering the order-48 group on the 2-site model. On Heisenberg rings,
the corrected relation checks also recover translation; the old order-16
detection result was an artifact of treating a relation basis as graph vertices.

**Detected maps need not preserve degree.** A Clifford can map a single-site
letter to a multi-site Pauli word while preserving the Hamiltonian. An
order-based moment basis need not be closed under it, and
`SolverConfig(symmetry=...)` then raises a closure error. Use a verified
basis-preserving subgroup or an explicitly closed basis. Detection does not
silently discard these valid symmetries to accommodate an SDP basis.

**Hermitian Pauli Hamiltonians only.** The current scope requires that every
Clifford image of a Pauli generator is a real ``(\pm 1)`` multiple of another
Pauli string. Generators whose phase lift would produce ``\pm i`` factors are
rejected with a clear error. This covers all physically relevant Hermitian
Hamiltonians but excludes certain non-Hermitian operators.

**Automorphism backend.** The current graph automorphism solver is a
deterministic in-tree backtracking search, adequate for moderate problem sizes.
A `backend=:bliss` keyword exists as a future extension point for the
[nauty/bliss](https://pallini.di.uniroma1.it/) family of dedicated graph
automorphism solvers, but currently falls back to the in-tree solver with a
warning.

## References

- [nation2026sympleq](@cite) — the SympleQ algorithm for Clifford symmetry
  detection in quantum many-body systems.
- [aaronson2004improved](@cite) — the Aaronson–Gottesman binary symplectic
  tableau for stabilizer simulation.
- [hostens2005stabilizer](@cite) — stabilizer formalism, phase update rules,
  and the symplectic representation for arbitrary dimensions.
- [rengaswamy2018synthesis](@cite) — symplectic synthesis of logical Clifford
  operators via Pauli isometries.
- [mckay2014practical](@cite) — nauty II, the practical graph automorphism
  solver underlying most implementations.
- [babai2016graph](@cite) — quasi-polynomial graph isomorphism (theoretical
  complexity bound).

## See also

- **Runnable examples:**
  [Pauli Symmetry Reduction](@ref pauli-clifford-symmetry) — manual Clifford
  gates, automatic SympleQ detection, the union of both, and a timing
  benchmark on a 10-site Heisenberg ring.
- **Architecture:**
  [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) — how NCTSSoS uses
  the detected symmetries to shrink the SDP.
- **API reference:**
  [`sympleq_symmetry_spec`](@ref), [`CliffordSymmetry`](@ref),
  [`SymmetrySpec`](@ref), [`SymmetryReport`](@ref).
