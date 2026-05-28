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
- **One symmetry family at a time**, chosen from:
  - `MonoidAlgebra` with [`SignedPermutation`](@ref) on registry indices;
  - `PauliAlgebra` with [`CliffordSymmetry`](@ref) on Pauli words;
  - `FermionicAlgebra` with [`FermionicModePermutation`](@ref),
    [`FermionicSectorSpec`](@ref), and optional
    [`FermionicSpinAdaptationSpec`](@ref).
- **Closed bases and invariant data** — every supplied generator must preserve
  the objective/constraints and map each moment/localizing half-basis back into
  itself.

Anything outside this scope raises an `ArgumentError` instead of silently
constructing the wrong relaxation. This is the right default while the path
matures: the failure mode for a bug in symmetry reduction is a wrong but
plausible-looking objective, which is much worse than a clear error.

In particular, [`cs_nctssos_higher`](@ref) does not yet accept a
`SymmetrySpec`. The GNS path is also not symmetry-aware today: it expects a
dense unreduced moment table, while a symmetry-reduced solve only builds the
reduced moment problem used by the SDP solver.

## [The pipeline](@id sab-pipeline)

When the user passes a [`SymmetrySpec`](@ref) via
`SolverConfig(; symmetry = ...)`, [`cs_nctssos`](@ref) routes the problem
through `moment_relax_symmetric` instead of the regular `moment_relax`. That
function performs the following steps in order:

1. **User supplies a `SymmetrySpec`.** A `SymmetrySpec` wraps one supported
   finite action family: [`SignedPermutation`](@ref) generators on monoid
   registry indices, [`CliffordSymmetry`](@ref) generators on Pauli words, or
   fermionic mode/sector/spin symmetry data. Indices or Pauli letters not
   listed by a finite generator are fixed.

2. **The repo enumerates the finite group locally.** Starting from the supplied
   generators, NCTSSoS does its own breadth-first closure on the relevant
   active domain. For signed permutations this is the active registry-index
   domain; for Clifford symmetry this is the action-closed Pauli-letter domain
   seeded by the problem support plus generator-touched sites; for fermionic
   mode permutations this is the active physical-mode domain.

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
   the selected finite action to a `NormalMonomial`, multiplying accumulated
   signs and re-running algebra-specific simplification where needed. The
   Clifford path composes signed Pauli-word images and tracks the Pauli phase;
   nothing in `SymbolicWedderburn` knows about these NC normal forms on its own.

6. **A thin adapter exposes that action to `SymbolicWedderburn`.**
   Small `SymbolicWedderburn.BySignedPermutations` adapters (`NCWordSignedPermutationAction`,
   `NCPauliCliffordAction`, and the fermionic mode adapter) plus
   `GroupsCore`-compatible finite-group wrappers are enough for
   `SymbolicWedderburn` to reason about the action without ever needing to look
   inside a `NormalMonomial`.

7. **`SymbolicWedderburn` computes the symmetry-adapted decomposition** of the
   half-basis via `SymbolicWedderburn.symmetry_adapted_basis(...)`. The output
   is a list of "blocks", one per isotypic component, each providing the
   change-of-basis rows used to project the moment matrix into that component.

8. **Transformed off-blocks are checked.** NCTSSoS keeps the diagonal blocks at
   the sizes returned by the decomposition, and checks that off-isotypic
   coupling vanishes after orbit reduction. If an off-block leaks through, the
   declared symmetry or the chosen basis is inconsistent and the call errors.

9. **Moment variables are orbit-reduced and a smaller `MomentProblem` is
   emitted.** For each monomial appearing in the relaxation, NCTSSoS computes
   its orbit under the group, picks a canonical representative, and either
   maps the monomial to ``\pm`` that representative or — when the orbit has
   inconsistent signs and is therefore forced to zero by the symmetry — drops
   it entirely. The objective, every constraint matrix, and every constraint
   block is rewritten in terms of these representatives. The constraint
   matrices, after the change of basis induced by the
   `SymbolicWedderburn`-supplied projectors, are also checked to have vanishing
   off-block entries within `_SYMMETRY_ATOL = 1e-10`; this catches
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

- The user-facing types `SignedPermutation`, `CliffordSymmetry`,
  `FermionicModePermutation`, `SymmetrySpec`, and `SymmetryReport`.
- The active-domain computation (`_symmetry_domain`, plus Pauli-letter and
  fermionic-mode variants).
- Group enumeration `_enumerate_symmetry_group` (BFS closure under composition,
  restricted to the relevant finite-action domain).
- The NC-word action `_act_monomial` / `_act_polynomial` — this is where
  algebra-specific simplification and Pauli phase bookkeeping happen.
- Invariance checks `_check_symmetry_invariance` and basis-closure checks
  `_check_basis_closure`.
- Orbit reduction `_build_orbit_reducer` and the rewriting helpers that map
  every polynomial in the relaxation to orbit representatives, with sign
  bookkeeping and zero-orbit handling.
- The construction of the final smaller `MomentProblem` and the populated
  `SymmetryReport`.

**`SymbolicWedderburn` responsibilities** (called via the thin adapter):

- The actual symmetry-adapted decomposition of the half-basis: computing
  characters, identifying isotypic components, and building the change-of-basis
  rows.

The adapter layer exposing the NC action to `SymbolicWedderburn` is
deliberately small: a `BySignedPermutations` action type, a `GroupsCore.Group`
wrapper (with multiplication, inverse, and element-order tables), and an
`action(..., g, mono)` method that returns the transformed monomial together
with its accumulated sign. Everything else stays in NCTSSoS.

## [Why the output is still an ordinary `MomentProblem`](@id sab-output)

After symmetry reduction, the JuMP/SDP layer downstream of
`moment_relax_symmetric` is unchanged: `solve_sdp` consumes a regular
`MomentProblem`, `sos_dualize` would consume a regular `SOSProblem`, the
moment variables are still keyed by `NormalMonomial`s. Only the *contents* of
the `MomentProblem` differ:

- The objective and every constraint are expressed using a small set of
  monomial **orbit representatives** rather than every monomial in the
  original basis.
- What used to be one large PSD block is now several smaller symmetry blocks.
- Off-block constraint entries that, after the change of basis and orbit
  reduction, must vanish are explicitly checked to be ``\le`` `_SYMMETRY_ATOL`
  in absolute value; if any leak through, the call errors out.

This is the central pragmatic gain of the architecture: **no downstream code
had to change**. `JuMP`, `Dualization`, `solve_sdp`, the result types, and most
solver-facing code still speak ordinary monomial moments. The catch is that
GNS reconstruction expects a *dense* unreduced moment table, so it is not
available from a symmetry-reduced solve today. The reduction is otherwise
invisible past the boundary of `moment_relax_symmetric`.

## [Current limitations and fail-fast behavior](@id sab-limits)

!!! note "Where to go next"
    The companion page [Extending Symmetry Support](@ref extending-symmetry)
    audits each of these limitations and what it would take to lift them,
    using the fermionic case as the worked example.


The MVP errors out, rather than silently doing the wrong thing, on any of:

- unsupported algebra/action combinations (for example, raw
  `SignedPermutation` on Pauli words instead of `CliffordSymmetry`);
- state or trace polynomial problems;
- any active correlative or term sparsity (`cs_algo`/`ts_algo` other than
  `NoElimination()`);
- more than one correlative-sparsity clique;
- term-sparsity block splitting inside a clique (more than one block basis per
  constraint);
- mixing finite action families in one `SymmetrySpec`;
- group elements that send a basis monomial outside the basis (basis closure
  fails);
- objective or constraint polynomials that are not invariant under the group
  (unless `check_invariance = false` is passed *and* the user has independently
  verified invariance);
- constraint blocks whose off-isotypic coupling fails to vanish within
  `_SYMMETRY_ATOL = 1e-10`;
- use via [`cs_nctssos_higher`](@ref), or GNS reconstruction from a
  symmetry-reduced solve.

These error messages name the violated assumption and the offending object
(generator, basis monomial, block index) so that misuse is debuggable rather
than mysterious.

## [Pauli Clifford symmetry](@id sab-pauli-clifford)

For Pauli problems, use [`CliffordSymmetry`](@ref) instead of raw
[`SignedPermutation`](@ref). Named constructors cover the common Clifford gate
actions on Pauli words: `:H`, `:S`, `:CNOT`, and `:SWAP`.

A two-site Heisenberg Hamiltonian is invariant under swapping the two qubits:

```julia
using NCTSSoS, COSMO

registry, (σx, σy, σz) = create_pauli_variables(1:2)
heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
pop = polyopt(heisenberg, registry)
basis = [one(σx[1]); σx; σy; σz]

spec = SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2))
config = SolverConfig(
    optimizer = COSMO.Optimizer,
    moment_basis = basis,
    cs_algo = NoElimination(),
    ts_algo = NoElimination(),
    symmetry = spec,
)

result = cs_nctssos(pop, config)
result.symmetry.group_order       # 2
result.symmetry.psd_block_sizes   # [4, 3]
```

Custom Clifford dictionaries are accepted too, but they are deliberately
validated as Pauli Clifford actions: sources must be single Pauli letters,
images must be non-identity Pauli words, the map must be injective up to sign,
and Pauli commutation plus local multiplication phases must be preserved. The
relaxation path supplies the sparse action domain automatically; direct
high-`nqubits` use of `CliffordSymmetryGroup` should pass `domain = ...` if it
only needs a local action.

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
basis-closure violation, malformed Clifford actions, and unsupported API
combinations) live in `test/relaxations/symmetry.jl`.

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
- Downstream solver code (JuMP construction, dualization, and result types)
  sees ordinary monomial moments and needs no changes.
- The symmetry feature is concentrated in one file and one boundary.

The cost is that the feature is still dense-only and ordinary-polynomial-only:
it does not yet compose with correlative/term sparsity, state/trace polynomial
optimization, higher-order continuation, or dense-moment GNS reconstruction.
Lifting those restrictions is a future extension of this same boundary, not a
rewrite of the rest of the package.

## See also

- The runnable walk-through:
  [CHSH with Symmetry Reduction](@ref chsh-symmetry).
- The contributor-facing roadmap for lifting MVP limitations (with the
  fermionic case as the worked example):
  [Extending Symmetry Support](@ref extending-symmetry).
- [`SignedPermutation`](@ref), [`CliffordSymmetry`](@ref),
  [`FermionicModePermutation`](@ref), [`SymmetrySpec`](@ref),
  [`SymmetryReport`](@ref), and the `symmetry` keyword of
  [`SolverConfig`](@ref).
- The reference acceptance case in
  `test/problems/bell_inequalities/chsh_simple.jl` (case `Symmetry_d1`).
- The MVP guard regression tests in `test/relaxations/symmetry.jl`, including
  Pauli Clifford validation and the 2-site Heisenberg SWAP regression.
- Related size-reduction strategies covered elsewhere in the manual:
  [Sparsities](@ref sparsities) and the
  [Moment Sum-of-Hermitian-Square Hierarchy](@ref moment-sohs-hierarchy).
