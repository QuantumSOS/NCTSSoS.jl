# [Clifford Symmetry Detection](@id clifford-symmetry-detection)

Given a nonconstant Pauli Hamiltonian, [`sympleq_symmetry_spec`](@ref) finds
generators of its Hamiltonian-preserving Clifford conjugation group. Global
unitary phases are irrelevant to conjugation. The site domain is `1:n`, where
`n` is the largest site index in the Hamiltonian.

Detection and moment-matrix reduction are separate operations. The default
detector searches the complete Clifford subgroup, subject to an explicit search
budget. The optional generator-only reduction backend avoids enumerating that
group, but requires a row basis preserved by its generators.

## [Background: the symplectic representation](@id sympleq-background)

A Pauli word has a binary label ``p=(x_1,\ldots,x_n,z_1,\ldots,z_n)``.
The labels of ``X``, ``Z``, and ``Y`` on one site are respectively ``(1,0)``,
``(0,1)``, and ``(1,1)``. Signs and the phase in ``Y=iXZ`` are tracked separately.
Clifford conjugation sends each Hermitian Pauli word to another such word with
a sign of ``+1`` or ``-1``. Its binary action is ``p\mapsto pS``, with

```math
S\Omega S^T=\Omega,\qquad
\Omega=\begin{pmatrix}0&I\\I&0\end{pmatrix}
```

over the field of two elements. This equation preserves commutation relations
[aaronson2004improved](@cite) [hostens2005stabilizer](@cite).

## [Graph candidates need algebraic checks](@id sympleq-insight)

The anticommutation graph has one vertex per Hamiltonian term, colored by
coefficient magnitude. Edges join anticommuting terms. Every Clifford symmetry
induces a color-preserving automorphism, but the converse needs two checks:

1. The term permutation must preserve **all binary linear relations**.
2. Its signs must satisfy the Hamiltonian's coefficient equations.

Checking only one chosen basis of relations is too restrictive: a translation
can preserve their span while changing the basis. Filtering the generators of
the unsigned graph group after the search is also insufficient: valid products
can be lost. The default search applies these checks inside the search for each
representative map.

## [Generator discovery](@id sympleq-algorithm)

The default `cycle_strategy=:relations` path uses a stabilizer chain:

1. Construct the tableau and anticommutation graph.
2. Choose independent term labels as base points.
3. For each base point, fix the earlier points and search for one valid map to
   each reachable image. Preserve linear relations during backtracking and
   solve the signed phase equations before accepting a map.
4. Extend each partial isometry to a symplectic map by Witt extension
   [rengaswamy2018synthesis](@cite).
5. Find the **term-fixing kernel**: maps invisible to term permutations. Fix the
   span of the Hamiltonian labels, extend it one independent vector at a time,
   and construct representative maps for the allowed images of each vector.
6. Add the Pauli sign characters that fix every Hamiltonian term.
7. Check every resulting generator directly against the original Hamiltonian.

The representatives at successive stabilizer levels generate the group; the
implementation does not form their Cartesian product. It can still visit
exponentially many candidates. Exceeding `max_search_nodes` raises an error,
not a partial result advertised as complete.

This kernel step recovers SWAP on the two-site Heisenberg Hamiltonian. Its full
Clifford conjugation stabilizer has order **192**, not 48. The order-48 group of
global axis rotations and site exchange is only a subgroup. Some additional
symmetries are entangling and do not preserve the seven degree-one Pauli rows.
Tests compare the detected subgroup with all 11,520 two-qubit Clifford actions
constructed independently by unitary conjugation.

## [Generator-only moment reduction](@id sympleq-bridge)

Here is a solver-free construction on two sites. The full Pauli row basis is
closed under every Clifford on those sites:

```julia
using NCTSSoS
reg, (x, y, z) = create_pauli_variables(1:2)
H = x[1]*x[2] + y[1]*y[2] + z[1]*z[2]
spec = sympleq_symmetry_spec(H; clifford_backend=:generators)
config = SolverConfig(
    optimizer=nothing,
    moment_basis=get_ncbasis(reg, 2),
    cs_algo=NoElimination(), ts_algo=NoElimination(), symmetry=spec,
)
pop = polyopt(H, reg)
sparsity = compute_sparsity(pop, config)
mp, report = NCTSSoS.moment_relax_symmetric(
    pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities, spec,
)
```

The backend performs these operations without a group list:

- Check generator invariance and row closure. Missing row images are errors;
  the row basis is never enlarged silently.
- Follow signed orbits of row pairs to obtain the invariant matrix space.
- Use two generic Hermitian samples to split irreducible copies and align
  equivalent copies, following the [RepLAB algorithm](https://arxiv.org/abs/1911.09154).
- Verify generator intertwining equations and compare the squared
  multiplicities with the exact dimension from the signed pair orbits.
- Build one multiplicity-size block per irreducible type. Repeated identical
  blocks are omitted.

These are floating-point checks, not rigorous interval certificates. Failure
raises an error; it does not silently drop blocks or switch to enumeration.
The complex decomposition handles real signed actions with complex or
quaternionic representation types without assuming all irreducibles are real.

`report.group_order` is `nothing`: the backend does not enumerate the group to
obtain a reporting statistic. Block provenance is `:generators`.

The existing `clifford_backend=:wedderburn` remains the default reduction
backend and a reference for small groups. It enumerates internally. The new
backend currently supports dense Pauli relaxations without charge/singlet
splitting. It accepts arbitrary supplied Clifford subgroups on closed rows.

## [Scaling and limitations](@id sympleq-limits)

- **Clifford symmetries only.** Continuous non-Clifford spin rotations are not
  detected by this method.
- **Row closure is separate from Hamiltonian invariance.** CNOT sends
  ``X_1`` to ``X_1X_2``. A degree-one row basis may therefore fail closure.
  Either choose a closed basis or explicitly supply a smaller subgroup.
- **Boundary conditions matter.** Cyclic translation preserves a periodic
  uniform ring, but generally not the corresponding open chain.
- **No universal speed guarantee.** Detection remains a bounded backtracking
  search. With ``r`` rows and ``k`` generators, signed pair traversal needs
  ``O(kr^2)`` work and ``O(r^2)`` storage; the dense eigendecomposition is
  ``O(r^3)``. Avoiding group enumeration does not remove row-dimension costs.
- **Legacy comparison.** `cycle_strategy=:cycle_basis` retains the old
  incomplete graph recognizer. Only that path uses `max_automorphisms` and the
  legacy `ga_backend=:bliss` fallback. Raw `sympleq_generators` also retains its
  legacy candidate-level contract; use `sympleq_symmetry_spec` for the complete
  checked path. The older [Pauli Symmetry Reduction](@ref pauli-clifford-symmetry)
  example selects the legacy strategy explicitly for its historical comparisons.

The SympleQ qubit-tapering algorithm is separate from this moment-matrix
reduction and is not implemented here [nation2026sympleq](@cite).

## References

- [nation2026sympleq](@cite): graph-based Clifford symmetry discovery.
- [aaronson2004improved](@cite), [hostens2005stabilizer](@cite): binary Clifford representations.
- [rengaswamy2018synthesis](@cite): extension of Pauli isometries.
- [RepLAB](https://arxiv.org/abs/1911.09154): numerical representation decomposition.

## See also

[Symmetry-Adapted Basis](@ref symmetry-adapted-basis),
[`CliffordSymmetry`](@ref), [`SymmetrySpec`](@ref), [`SymmetryReport`](@ref).
