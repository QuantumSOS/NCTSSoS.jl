# [Symmetry-Adapted Basis](@id symmetry-adapted-basis)

Many noncommutative polynomial problems have a symmetry: relabel the operators,
possibly with signs, and the objective and constraints stay unchanged. CHSH is
one example. Swapping Alice's measurements with a sign flip on Bob's second
measurement preserves its objective.

| Input | Pipeline operation | Output |
|:------|:-------------------|:-------|
| `Problem(H; symmetry=...)` | Check invariance | Validated actions |
| `Relaxation(...)` | Select and close the row-block family | Block orbits |
| Closed family and actions | Project rows and identify invariant moments | Smaller conic program |
| Program and solver | Solve and map represented moments back | Bound and expectations |

This page covers **using** a known symmetry to reduce an SDP. The complementary
task, **finding** Clifford symmetries, is described in
[Clifford Symmetry Detection](@ref clifford-symmetry-detection).
See [CHSH with Symmetry Reduction](@ref chsh-symmetry) for a complete example.

## [What "symmetry-adapted basis" means here](@id sab-meaning)

Let ``B`` be the rows of a moment matrix and ``G`` a symmetry group that
preserves them. Averaging a feasible state over ``G`` preserves the objective
and feasibility of an invariant problem. Its moment matrix commutes with the
group action. A change of basis ``U`` separates the independent symmetry
components:

```math
U^\dagger M_B U = \operatorname{blockdiag}(M^{(1)}, M^{(2)}, \ldots).
```

Thus the invariant PSD constraint can be expressed with smaller PSD blocks.
For real actions, the adjoint above is a transpose. This is an exact symmetry
reduction of the invariant relaxation, not the removal of arbitrary matrix
entries. Adding a physical sector restriction is a different operation: it
changes which states are allowed.

## [Key architectural decision: internal projected rows](@id sab-decision)

The polynomial interface stays in the original operators. The pipeline uses
linear combinations of row words internally to split cones, and shares moment
coordinates across the resulting blocks. It does not replace the polynomial
algebra with irreducible-representation coordinates.

The same `Problem` and `Relaxation` declarations work for the direct conic
writers and the optional JuMP backend. The solver does not need to know the
group action. Use `expectation(result, observable)` to read an observable;
do not interpret a reduced coordinate as an original monomial by its index.

## [Supported scope](@id sab-scope)

- Use `LetterPermutation` for named letters in free, unipotent and projector
  site algebras. Negative targets are supported for free and unipotent letters;
  the negative of a projector is not a projector and is rejected.
- Use `SitePermutation` for site relabelling and `CliffordMap` for Pauli maps.
  `detect_symmetry` finds verified Clifford generators for an unconstrained
  Pauli objective. It can return a subgroup, not the full symmetry group.
- Declare actions on `Problem`. A tuple of actions is sufficient for ordinary
  block symmetry. `Next.Symmetry` also accepts `moment_only` and `continuous`
  generators. Antiunitary actions constrain moments; they are not unitary
  block generators.
- Symmetry works with selected block families, not only a single dense block.
  Automatically selected families are closed under block generators. A custom
  family must already be closed; the pipeline does not silently replace it.
- State and trace polynomials currently support the plain dense route only.
  Their functional block pipeline is not available yet.

## [The pipeline](@id sab-pipeline)

1. **Declare the invariant problem.** The pipeline checks the objective and
   each operator constraint under the supplied actions. A failed check names
   the polynomial and generator. There is no public flag to bypass it.
2. **Select row blocks.** An order, custom rows, or a custom family determines
   the starting blocks. Grouping and term splitting select smaller blocks when
   requested. They are distinct from symmetry reduction.
3. **Close the family.** Apply block generators until the family is closed.
   Blocks that map to one another form an orbit. One representative per orbit
   is enough when the shared moments obey the same symmetry.
4. **Reduce each representative.** Its stabilizer is the subgroup that maps
   that block to itself. The finite-group route builds shared group-algebra
   projectors with `SymbolicWedderburn` and applies them to row orbits. Lattice
   and continuous routes use their own structure where applicable.
5. **Identify moment coordinates.** Block and moment-only actions induce
   relations between moments, including sign-forced zeros. Projected entries
   remain linear forms in these shared coordinates.
6. **Apply requested strengthening and emit the conic program.** The result
   records the state domain for which the bound is valid. The conic writer
   receives the reduced blocks, objective and linear equations.

Finite block reduction is not limited to scalar or multiplicity-free blocks.
Continuous generators, such as `spin_rotations(registry)`, express spin
invariance without manually enumerating a large finite subgroup. They do not
by themselves impose a fixed spin sector or a singlet state.

## [What stays local vs. what `SymbolicWedderburn` does](@id sab-boundary)

**NCTSSoS owns** the operator algebra, the action on words, invariance checks,
block-family closure, block orbits and stabilizers, moment relations,
strengthening, and conic assembly. These operations live in `src/next/core/`,
`src/next/model/` and `src/next/pipeline/`.

**`SymbolicWedderburn` supplies** character and group-algebra projection tools
for the finite-group route. It does not implement noncommutative word
multiplication or decide which physical states a strengthening condition allows.
The local adapter applies its projectors to the algebra's row words.

## [Why the output keeps ordinary observables](@id sab-output)

The public staged path makes the boundary visible. Here `problem`,
`relaxation` and `solver` are declarations such as those in the complete
CHSH example below:

```julia
built = Next.build(problem, relaxation)
built.sizes
prepared = Next.prepare(built; solver)
result = solve(built; solver)
expectation(result, problem.objective)
```

`Next.build` does not solve. `Next.prepare` prepares the selected backend;
`solve(built; solver)` can also prepare it automatically. The result retains
the original problem and a map for represented moments. Not every observable
must be represented by a sparse relaxation. Reconstruction uses the checked
entry `gns(result)`; a symmetry-reduced solution alone does not certify flatness
or an optimizer.

`result.sizes` gives flat counts: `rows`, `moments`, `real_variables`, `blocks`
and `largest_block`. `moments` includes the identity; `real_variables` counts
free real coordinates. `sizes.closure` is `nothing` on the plain dense route
and reports family closure on the block route. To see individual cone sides,
use `Next.block_sides(built)`. It returns the final PSD/HPSD sides, with
Hermitian sides before realification. The generic closure report does not
provide the finite group order.

## [Current limitations and fail-fast behavior](@id sab-limits)

The pipeline rejects incompatible action/algebra pairs, non-invariant
objectives or operator constraints, and custom families that are not closed.
An action that maps rows outside a block must be handled by family closure or
used only for moment relations; it cannot be applied as a projection on an
unclosed row set. The generic finite-group route has a limit of 65,535 elements.
State and trace block routing remains unavailable.

The former check-bypass and off-block-check options are not part of this API.
Use the diagnostic and fix the declaration rather than disable its checks.

## [Reference acceptance case: CHSH](@id sab-chsh)

The rows are ``Q=\{1,A_1,A_2,B_1,B_2\}``. Alice swap, Bob swap and party swap
generate a group of order 16. Both the dense and reduced solves recover
``-2\sqrt{2}`` within `1e-6`. The reduced program has three scalar PSD blocks
and one free real coordinate. Its two moments include the identity.

```julia
using NCTSSoS.Next, COSMO
import NCTSSoS.Next as Next

reg, (A, B) = site_variables(Unipotent(), 2, 2)
S = A[1]*B[1] + A[1]*B[2] + A[2]*B[1] - A[2]*B[2]
alice = LetterPermutation(A[1]=>A[2], A[2]=>A[1], B[2]=>-B[2])
bob = LetterPermutation(A[2]=>-A[2], B[1]=>B[2], B[2]=>B[1])
party = LetterPermutation(A[1]=>B[1], A[2]=>B[2], B[1]=>A[1], B[2]=>A[2])
problem = Problem(-S; symmetry=(alice, bob, party))
relaxation = Relaxation(rows=[one(A[1]), A[1], A[2], B[1], B[2]])
solver = JuMPBackend(COSMO.Optimizer)
built = Next.build(problem, relaxation)
result = solve(built; solver)
@assert abs(result.bound + 2sqrt(2)) < 1e-6
@assert Next.block_sides(built) == [1, 1, 1]
@assert result.sizes.real_variables == 1
```

The [runnable CHSH page](@ref chsh-symmetry) uses the direct Mosek writer and
also checks the dense baseline and the expectation of the CHSH operator.

## [Contrast with symmetry-adapted polynomial coordinates](@id sab-contrast)

A different design could express every polynomial in irreducible-component
coordinates. Multiplication, algebraic reduction, sparsity and observable
queries would then all need that representation. Here only the internal row
projection uses symmetry coordinates. The user still writes the Hamiltonian
and observables in the original operators, and solver backends consume the
same conic contract whether or not symmetry is present.

## See also

- [CHSH with Symmetry Reduction](@ref chsh-symmetry).
- [Pauli Symmetry Reduction](@ref pauli-clifford-symmetry).
- [Clifford Symmetry Detection](@ref clifford-symmetry-detection).
- [Pauli Charge/Singlet Symmetry](@ref pauli-charge-singlet-symmetry).
- [Sparsities](@ref sparsities) and the
  [Moment Sum-of-Hermitian-Square Hierarchy](@ref moment-sohs-hierarchy).
