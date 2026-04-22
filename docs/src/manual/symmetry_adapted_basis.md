# [Symmetry-Adapted Basis](@id symmetry-adapted-basis)

Many noncommutative polynomial optimization problems carry a structural
symmetry: relabel the operators in some prescribed way, possibly with a sign
flip, and the objective and constraints come back unchanged. The CHSH Bell
inequality is the canonical example — it is invariant under swapping Alice's
two measurements (with a sign flip), under swapping Bob's two measurements
(with a sign flip), and under swapping the two parties.

When such a symmetry exists, the moment matrix that the relaxation builds is
not just positive semidefinite — it commutes with the corresponding group
action. By the standard representation-theoretic argument (the *symmetry-adapted
basis*), this large PSD constraint splits into several smaller, independent PSD
blocks, one per isotypic component. NCTSSoS exploits this to shrink the SDP
the user actually sends to the solver.

This page describes how that support is wired in **today**, what is and isn't
in scope, and how the pieces interact. It is intentionally aligned with the
implemented MVP, not a broader future redesign.

## [What "symmetry-adapted basis" means here](@id sab-meaning)

In the linear-algebra sense, a symmetry-adapted basis is a change of basis on
the half-basis ``B`` (the monomials indexing the rows and columns of the moment
matrix) that block-diagonalizes the action of the symmetry group ``G`` on
``B``:

```math
U^\top \, M_B \, U \;=\; \mathrm{blockdiag}\!\big(M^{(1)}, M^{(2)}, \ldots\big),
```

where each ``M^{(i)}`` corresponds to one isotypic component and the ``U`` is
built from ``G``-equivariant projectors. Because the moment matrix is
``G``-invariant, the off-diagonal coupling vanishes and the *single* large PSD
constraint ``M_B \succeq 0`` is equivalent to *several* smaller PSD constraints
``M^{(i)} \succeq 0``.

In NCTSSoS this idea is used **internally only**, to decide how to split one
PSD block into several smaller PSD blocks for the solver. Externally, the
relaxation continues to live in ordinary monomial coordinates: every PSD entry,
every constraint, every moment variable is still a monomial (or a polynomial in
monomials). The user does not have to think in terms of irreducible
representations to use the feature.

## [Key architectural decision: thin internal adapter](@id sab-decision)

There are two reasonable ways to bolt symmetry adaptation onto a moment-SOS
package:

1. **Rewrite-the-package approach.** Replace the canonical monomial basis with
   symmetry-adapted linear combinations everywhere — in the polynomial type,
   the simplification rules, the term-sparsity graph, the registries. Every
   moment variable becomes a coefficient on an irreducible component.

2. **Thin-adapter approach.** Keep monomials and polynomials exactly as they
   are. Build the moment problem normally. Use the symmetry-adapted basis only
   as a tool to decide *how to slice* one PSD constraint into several smaller
   ones, and emit those smaller blocks back in ordinary monomial coordinates
   (after orbit-reducing moment variables to one representative per orbit).

NCTSSoS deliberately takes approach **(2)**. The rest of the package — the
algebra hierarchy, the registries, `Polynomial`, `NormalMonomial`, the
correlative/term-sparsity machinery, the JuMP construction in
`MomentProblem`/`SOSProblem` — is not touched. The symmetry path lives in a
single file (`src/optimization/symmetry.jl`) and only takes over one job:
constructing a *reduced* `MomentProblem` from the same inputs the dense path
uses.

The output of the symmetry path is therefore an ordinary
[`MomentProblem`](@ref NCTSSoS.MomentProblem) just like the one produced by the
non-symmetric path — just with fewer moment variables and several small PSD
blocks instead of one big one.

## [Supported MVP scope](@id sab-scope)

The current MVP is intentionally narrow and is enforced by fail-fast checks in
[`SolverConfig`](@ref) handling. Symmetry reduction is supported only when
**all** of the following hold:

- **Dense** relaxation — `cs_algo = NoElimination()` and
  `ts_algo = NoElimination()`.
- **Single clique** — exactly one correlative-sparsity clique.
- **Ordinary polynomial problems** — `PolyOpt` over `Polynomial`. State and
  trace polynomial problems are not yet supported on the symmetry path.
- **Monoid algebras only** — `A <: MonoidAlgebra` (i.e.
  `NonCommutative`, `Projector`, `Unipotent`). PBW algebras (Bosonic,
  Fermionic) and `TwistedGroupAlgebra` (Pauli) are not yet supported.
- **Real signed-permutation actions only** — every group element acts on
  registry indices via a permutation possibly composed with a sign of ``\pm 1``
  per index. No general orthogonal action, no complex characters.
- **Multiplicity-free decompositions only** — every isotypic component returned
  by the internal decomposition must have multiplicity 1 and reduce to a scalar
  (1×1) block.

Anything outside this scope raises an `ArgumentError` instead of silently
constructing the wrong relaxation. This is the right default while the path
matures: the failure mode for a bug in symmetry reduction is a wrong but
plausible-looking objective, which is much worse than a clear error.

In particular, [`cs_nctssos_higher`](@ref) does not yet accept a
`SymmetrySpec`. The monomap-based GNS path is also not symmetry-aware today:
a symmetry-reduced `result.monomap` does not contain the dense unreduced
moments that `gns_reconstruct` needs, so it currently fails later with an
explicit missing-moments error rather than at the symmetry API boundary.

## [The pipeline](@id sab-pipeline)

When the user passes a [`SymmetrySpec`](@ref) via
`SolverConfig(; symmetry = ...)`, [`cs_nctssos`](@ref) routes the problem
through `moment_relax_symmetric` instead of the regular `moment_relax`. That
function performs the following steps in order:

1. **User supplies a `SymmetrySpec`.** A `SymmetrySpec` wraps a list of
   [`SignedPermutation`](@ref) generators acting on registry indices. Each
   generator records, per index ``i``, a target index ``j`` and a sign
   ``\sigma \in \{-1, +1\}`` such that the operator with index ``i`` is sent
   to ``\sigma \cdot x_j``. Indices not listed are fixed.

2. **The repo enumerates the signed-permutation group locally.** Starting from
   the supplied generators, NCTSSoS does its own breadth-first closure on
   signed permutations restricted to the *active variable domain* (the union
   of all registry indices that actually appear in the objective, the
   constraints, the moment-matrix basis, and any term-sparsity block bases).
   This produces the full finite group as a `Vector{SignedPermutation}`.

3. **Invariance of objective and constraints is checked.** Each enumerated
   group element is applied to the objective and to every equality, inequality,
   and moment-equality constraint. If any of these polynomials changes, the
   call errors out (this can be opted out of via
   `SymmetrySpec(...; check_invariance = false)`, but the default is to
   check). This is a hard prerequisite: if the *problem* is not invariant,
   block-diagonalizing its moment matrix is meaningless.

4. **Basis closure is checked.** The moment-matrix basis of the single clique,
   and the block bases of all term-sparsity blocks, must be mapped to
   themselves (up to sign) by every group element. If any element of the
   half-basis is sent outside the basis, the call errors out. This is what
   guarantees that the action descends to a well-defined linear action on the
   span of the half-basis.

5. **A local monomial action is built.** The function `_act_monomial` applies
   a `SignedPermutation` to a `NormalMonomial` letter-by-letter, multiplying
   the accumulated signs and re-running `simplify` for the relevant
   `MonoidAlgebra`. This is NCTSSoS's NC-aware action; nothing in
   `SymbolicWedderburn` knows about non-commutative words on its own.

6. **A thin adapter exposes that action to `SymbolicWedderburn`.**
   The `NCWordSignedPermutationAction <: SymbolicWedderburn.BySignedPermutations`
   struct, together with a tiny `GroupsCore`-compatible wrapper around the
   enumerated `SignedPermutation`s (`_SignedPermutationGroup` /
   `_SignedPermutationGroupElement`, with explicit multiplication, inverse,
   and order tables), is enough for `SymbolicWedderburn` to reason about the
   action without ever needing to look inside a `NormalMonomial`.

7. **`SymbolicWedderburn` computes the symmetry-adapted decomposition** of the
   half-basis via `SymbolicWedderburn.symmetry_adapted_basis(...)`. The output
   is a list of "blocks", one per isotypic component, each providing the
   change-of-basis rows used to project the moment matrix into that component.

8. **Unsupported decompositions are rejected.** Each returned block must (a)
   have multiplicity 1, (b) be simple, and (c) be 1×1. If any of those checks
   fails, the call errors out with a message indicating which assumption was
   violated. This is what enforces the "scalar / multiplicity-free only" part
   of the MVP scope.

9. **Moment variables are orbit-reduced and a smaller `MomentProblem` is
   emitted.** For each monomial appearing in the relaxation, NCTSSoS computes
   its orbit under the group, picks a canonical representative, and either
   maps the monomial to ``\pm`` that representative or — when the orbit has
   inconsistent signs and is therefore forced to zero by the symmetry — drops
   it entirely. The objective, every constraint matrix, and every constraint
   block is rewritten in terms of these representatives. The constraint
   matrices, after the change of basis induced by the
   `SymbolicWedderburn`-supplied projectors, are also checked to reduce to a
   single 1×1 scalar block per isotypic component (off-block entries must
   vanish numerically, within `_SYMMETRY_ATOL = 1e-10`); this catches
   inconsistencies between the declared symmetry and the actual algebraic
   structure.

The function returns both the reduced `MomentProblem` and a
[`SymmetryReport`](@ref) summarizing what happened: group order, number of
distinct invariant moment variables, list of PSD block sizes, and the original
half-basis size. On the final [`PolyOptResult`](@ref),
`moment_matrix_sizes` deliberately keeps its old meaning — the pre-symmetry
term-sparsity layout — while the actual reduced solver PSD block sizes live in
`result.symmetry.psd_block_sizes` and are printed separately.

## [What stays local vs. what `SymbolicWedderburn` does](@id sab-boundary)

Keeping the boundary explicit makes it easy to upgrade either side later.

**NCTSSoS-local responsibilities** (in `src/optimization/symmetry.jl`):

- The user-facing types `SignedPermutation`, `SymmetrySpec`, `SymmetryReport`.
- The active-domain computation `_symmetry_domain`.
- Group enumeration `_enumerate_symmetry_group` (BFS closure under composition,
  restricted to the active domain).
- The NC-word action `_act_monomial` / `_act_polynomial` — this is where
  algebra-specific simplification happens, via `simplify(A, ...)`.
- Invariance checks `_check_symmetry_invariance` and basis-closure checks
  `_check_basis_closure`.
- Orbit reduction `_build_orbit_reducer` and the rewriting helpers that map
  every polynomial in the relaxation to orbit representatives, with sign
  bookkeeping and zero-orbit handling.
- The construction of the final smaller `MomentProblem` and the populated
  `SymmetryReport`.

**`SymbolicWedderburn` responsibilities** (called via the thin adapter):

- The actual symmetry-adapted decomposition of the half-basis: computing
  characters, identifying isotypic components, building the change-of-basis
  rows.
- The `multiplicity`, `issimple`, and `degree` queries used to gate which
  decompositions the MVP accepts.

The adapter layer exposing the NC action to `SymbolicWedderburn` is
deliberately small: a `BySignedPermutations` action type, a `GroupsCore.Group`
wrapper (with multiplication, inverse, and element-order tables), and an
`action(::NCWordSignedPermutationAction, g, mono)` method that returns the
permuted monomial together with its accumulated sign. Everything else stays
in NCTSSoS.

## [Why the output is still an ordinary `MomentProblem`](@id sab-output)

After symmetry reduction, the JuMP/SDP layer downstream of
`moment_relax_symmetric` is unchanged: `solve_sdp` consumes a regular
`MomentProblem`, `sos_dualize` would consume a regular `SOSProblem`, the
moment variables are still keyed by `NormalMonomial`s. Only the *contents* of
the `MomentProblem` differ:

- The objective and every constraint are expressed using a small set of
  monomial **orbit representatives** rather than every monomial in the
  original basis.
- What used to be one large PSD block is now several small blocks (in the MVP
  scope, all 1×1).
- Off-block constraint entries that, after the change of basis and orbit
  reduction, must vanish are explicitly checked to be ``\le`` `_SYMMETRY_ATOL`
  in absolute value; if any leak through, the call errors out.

This is the central pragmatic gain of the architecture: **no downstream code
had to change**. `JuMP`, `Dualization`, `solve_sdp`, the result types, and the
GNS-side machinery all still speak ordinary monomial moments. The catch is
that monomap-based GNS reconstruction expects a *dense* unreduced moment table,
so a symmetry-reduced `result.monomap` is not directly consumable today. The
reduction is otherwise invisible past the boundary of `moment_relax_symmetric`.

## [Current limitations and fail-fast behavior](@id sab-limits)

!!! note "Where to go next"
    The companion page [Extending Symmetry Support](@ref extending-symmetry)
    audits each of these limitations and what it would take to lift them,
    using the fermionic case as the worked example.


The MVP errors out, rather than silently doing the wrong thing, on any of:

- non-monoid algebras (`PBWAlgebra`, `TwistedGroupAlgebra`);
- state or trace polynomial problems;
- any active correlative or term sparsity (`cs_algo`/`ts_algo` other than
  `NoElimination()`);
- more than one correlative-sparsity clique;
- term-sparsity block splitting inside a clique (more than one block basis per
  constraint);
- group elements that send a basis monomial outside the basis (basis closure
  fails);
- objective or constraint polynomials that are not invariant under the group
  (unless `check_invariance = false` is passed *and* the user has independently
  verified invariance);
- decompositions returned by `SymbolicWedderburn` whose blocks have
  multiplicity > 1, are non-simple, or have size > 1×1;
- constraint blocks whose off-isotypic coupling fails to vanish within
  `_SYMMETRY_ATOL = 1e-10`;
- use via [`cs_nctssos_higher`](@ref), or monomap-based GNS reconstruction
  from a symmetry-reduced solve.

These error messages name the violated assumption and the offending object
(generator, basis monomial, block index) so that misuse is debuggable rather
than mysterious.

## [Reference acceptance case: CHSH](@id sab-chsh)

The current acceptance case is the CHSH Bell inequality with the standard
half-basis ``Q = \{1, x_1, x_2, y_1, y_2\}``. The full setup lives in
`test/problems/bell_inequalities/chsh_simple.jl` and the corresponding
expectation fixture in `test/data/expectations/chsh_simple.toml`
(`Symmetry_d1`).

The symmetry group is generated by:

- an "Alice swap": ``x_1 \leftrightarrow x_2`` together with ``y_2 \to -y_2``;
- a "Bob swap": ``x_2 \to -x_2`` together with ``y_1 \leftrightarrow y_2``;
- a "party swap": ``x_i \leftrightrightarrows y_i`` for ``i = 1, 2``.

These three signed permutations generate a group of order 16. The symmetry
path produces the following acceptance values:

| quantity                  | value                |
|---------------------------|----------------------|
| objective                 | ``\approx -2\sqrt{2}``, observed ``-2.82842712455161`` |
| `group_order`             | 16                   |
| `invariant_moment_count`  | 1                    |
| `psd_block_sizes`         | `[1, 1, 1]`          |

In other words: one PSD block of side 5 (the dense moment matrix on ``Q``)
becomes three independent 1×1 PSD blocks, the eleven monomials in the dense
relaxation are collapsed to one nontrivial invariant moment plus the identity,
and the recovered Tsirelson bound matches ``-2\sqrt{2}`` to solver tolerance.
A representative call:

```julia
using NCTSSoS, COSMO

reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]
pop = polyopt(-f, reg)

alice_swap = SignedPermutation(
    x[1].word[1] => x[2].word[1],
    x[2].word[1] => x[1].word[1],
    y[2].word[1] => (-1, y[2].word[1]),
)
bob_swap = SignedPermutation(
    x[2].word[1] => (-1, x[2].word[1]),
    y[1].word[1] => y[2].word[1],
    y[2].word[1] => y[1].word[1],
)
party_swap = SignedPermutation(
    x[1].word[1] => y[1].word[1],
    x[2].word[1] => y[2].word[1],
    y[1].word[1] => x[1].word[1],
    y[2].word[1] => x[2].word[1],
)

config = SolverConfig(
    optimizer  = COSMO.Optimizer,
    moment_basis = [one(x[1]), x[1], x[2], y[1], y[2]],
    cs_algo    = NoElimination(),
    ts_algo    = NoElimination(),
    symmetry   = SymmetrySpec(alice_swap, bob_swap, party_swap),
)

result = cs_nctssos(pop, config)
result.objective              # ≈ -2sqrt(2)
result.symmetry.group_order   # 16
result.symmetry.psd_block_sizes  # [1, 1, 1]
```

Regression tests for the MVP failure modes (non-invariant objective,
basis-closure violation, multiplicity > 1) live in
`test/relaxations/symmetry.jl`.

## [Contrast with a full symmetry-adapted-coordinate redesign](@id sab-contrast)

It is worth saying explicitly what NCTSSoS is *not* doing, because the phrase
"symmetry-adapted basis" is sometimes used to describe a much more invasive
change.

A hypothetical full redesign would push symmetry adaptation all the way down
into the type system. Polynomials would be expressed natively as linear
combinations of *isotypic basis elements* rather than monomials; multiplication
would carry Clebsch–Gordan-like coefficients; simplification would be defined
on irreducible-representation labels; the `VariableRegistry` itself would be
indexed by ``(\text{irrep}, \text{multiplicity index}, \text{component})``;
and term-sparsity, correlative-sparsity, GNS, and reconstruction would all
have to be re-derived in those coordinates.

NCTSSoS deliberately does **none** of that. It keeps every existing data
structure and rewrite rule intact and uses the symmetry-adapted basis only as
an *internal* tool, late in the pipeline, to slice one PSD block into several.
The benefits of that conservative choice are:

- The non-symmetry path is exactly as before — no risk of regression in the
  much larger code surface that does not use symmetry.
- All downstream code (JuMP construction, dualization, GNS, reconstruction,
  result types) sees ordinary monomial moments and needs no changes.
- The symmetry feature is concentrated in one file and one boundary.

The cost is that the MVP is restricted to multiplicity-free, scalar-block
cases. Higher-dimensional irreducible blocks (multiplicity > 1, non-scalar
isotypic components) require a richer intermediate representation that the
package does not yet maintain. Lifting that restriction is a future extension
of this same boundary, not a rewrite of the rest of the package.

## See also

- The runnable walk-through:
  [CHSH with Symmetry Reduction](@ref chsh-symmetry).
- The contributor-facing roadmap for lifting MVP limitations (with the
  fermionic case as the worked example):
  [Extending Symmetry Support](@ref extending-symmetry).
- [`SignedPermutation`](@ref), [`SymmetrySpec`](@ref),
  [`SymmetryReport`](@ref), and the `symmetry` keyword of
  [`SolverConfig`](@ref).
- The reference acceptance case in
  `test/problems/bell_inequalities/chsh_simple.jl` (case `Symmetry_d1`).
- The MVP guard regression tests in `test/relaxations/symmetry.jl`.
- Related size-reduction strategies covered elsewhere in the manual:
  [Sparsities](@ref sparsities) and the
  [Moment Sum-of-Hermitian-Square Hierarchy](@ref moment-sohs-hierarchy).
