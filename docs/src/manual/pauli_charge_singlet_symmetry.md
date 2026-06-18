# [Pauli Charge, Spatial, and Singlet Symmetry](@id pauli-charge-singlet-symmetry)

The first Pauli production symmetry path targets dense order-2 moment
relaxations for XXX-style spin Hamiltonians. It combines three reductions that
commute with each other when used carefully:

1. **U(1) charge sectors** in the local basis ``\{Z, S^+, S^-\}``, with
   charges ``q(Z)=0``, ``q(S^+)=+1``, and ``q(S^-)=-1``.
2. **Spatial site symmetries** represented as Pauli [`CliffordSymmetry`](@ref)
   actions that preserve the Pauli axis and only permute sites.
3. **Order-2 SU(2)-singlet moment equalities**, added as linear expectation
   constraints.

This is not a new solver pipeline. The feature is wired through
[`SymmetrySpec`](@ref) and the internal `moment_relax_symmetric` path just like
the older signed-permutation and Clifford reductions. The output is still an ordinary
[`MomentProblem`](@ref NCTSSoS.MomentProblem), with smaller PSD blocks and fewer
moment variables.

## [When to use it](@id pauli-charge-when)

Use this path for Pauli Hamiltonians whose objective and constraints are
charge-neutral and whose visible finite symmetries are lattice symmetries:
translations, reflections, rotations of sites, or any other site permutation
that preserves ``X``, ``Y``, and ``Z`` axes.

The intended first target is the periodic XXX Heisenberg chain

```math
H_N = \frac{1}{4}\sum_{i=1}^{N}\sum_{a\in\{x,y,z\}}
      \sigma_i^a\sigma_{i+1\; (\mathrm{mod}\; N)}^a.
```

Do **not** combine Pauli charge/singlet reduction with axis-mixing Cliffords
such as global Hadamards or phase rotations. Those are valid Clifford
symmetries for isotropic Hamiltonians, but they do not preserve the
``\{Z,S^+,S^-\}`` charge basis. NCTSSoS enforces this: when
`pauli_charge` or `pauli_singlet` is present, every supplied Clifford group
element must be a spatial site permutation.

## [API pattern](@id pauli-charge-api-pattern)

Use [`pauli_site_permutation`](@ref) for translations and reflections. It
constructs a [`CliffordSymmetry`](@ref) that maps ``\sigma_i^a`` to
``\sigma_{p(i)}^a`` for every Pauli axis ``a``.

```julia
using NCTSSoS, MosekTools

N = 16
registry, (σx, σy, σz) = create_pauli_variables(1:N)

H = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (σx, σy, σz) for i in 1:N
)

translation = pauli_site_permutation([2:N; 1])
reflection = pauli_site_permutation(reverse(1:N))

symmetry = SymmetrySpec(
    [translation, reflection];
    pauli_charge = PauliChargeSectorSpec(nqubits=N),
    pauli_singlet = PauliSingletConstraintSpec(nqubits=N),
)

config = SolverConfig(
    optimizer = Mosek.Optimizer,
    order = 2,
    cs_algo = NoElimination(),
    ts_algo = NoElimination(),
    symmetry = symmetry,
)

result = cs_nctssos(polyopt(H, registry), config)
```

The same `SymmetrySpec` can also be used without a solver to inspect the
symbolic reduction: set `optimizer = nothing` and call
`moment_relax_symmetric` through the normal `compute_sparsity` output, as the
regression tests do.

## [What each option does](@id pauli-charge-options)

[`PauliChargeSectorSpec`](@ref) changes the moment basis from Pauli words to
charge words through order 2. The PSD matrix is block-partitioned by total
charge. If spatial Clifford generators are also supplied, NCTSSoS runs the
finite Wedderburn split **inside each charge sector**.

[`PauliSingletConstraintSpec`](@ref) appends order-2 expectation equalities for
an SU(2)-singlet state:

- single-site expectations vanish: ``\langle \sigma_i^a\rangle = 0``;
- cross-component two-point correlators vanish:
  ``\langle \sigma_i^a\sigma_j^b\rangle = 0`` for ``a\ne b``;
- same-component two-point correlators are equal:
  ``\langle \sigma_i^x\sigma_j^x\rangle =
    \langle \sigma_i^y\sigma_j^y\rangle =
    \langle \sigma_i^z\sigma_j^z\rangle``.

These are linear moment constraints, not a full non-Abelian SU(2)
representation-theory implementation. That is the right level for the first
order-2 Pauli path: small, explicit, and easy to audit.

[`PauliChargeBlockLabel`](@ref) labels the produced PSD blocks. Its fields tell
you the charge sector, whether a further finite spatial block split was used,
the spatial group order, and the dimension of the unsplit charge sector.

## [Current scope and guardrails](@id pauli-charge-scope)

The implementation is intentionally narrow:

- only ordinary [`PolyOpt`](@ref) problems over `PauliAlgebra`;
- dense relaxations only: `cs_algo = NoElimination()` and
  `ts_algo = NoElimination()`;
- complete Pauli half-bases through degree 0, 1, or 2; `order = 2` is the
  production case;
- charge-neutral objectives and constraints when `PauliChargeSectorSpec` is
  enabled;
- optional finite symmetries must be spatial `CliffordSymmetry` actions that
  preserve Pauli axes; use [`pauli_site_permutation`](@ref);
- `PauliSingletConstraintSpec` currently emits only order-2 singlet moment
  equalities.

Violations raise `ArgumentError`. That is deliberate. A symmetry reduction bug
can give a plausible but wrong bound; failing early is cheaper.

## [N=16 XXX Heisenberg evidence](@id pauli-charge-n16-evidence)

The order-2 periodic XXX Heisenberg chain at ``N=16`` is the reference scale
case for this feature.

| quantity | dense order-2 relaxation | charge + spatial + singlet reduction |
|:--|--:|--:|
| half-basis dimension | 1129 | 1129 before reduction |
| scalar PSD variables | 637885 | 5108 |
| spatial group order | — | 32 |
| PSD blocks | one `1129×1129` block | 53 blocks |
| largest PSD block | 1129 | 24 |
| invariant moment count | — | 2669 |

The Mosek validation run gave:

| quantity | value |
|:--|--:|
| total `cs_nctssos` wall time | 2146.539 s |
| Mosek solve time | 9.449 s |
| relaxation objective | -7.153492909789651 |
| relaxation per site | -0.44709330686185317 |
| exact diagonalization energy | -7.142296360616788 |
| exact diagonalization per site | -0.44639352253854925 |
| relative lower-bound gap vs. `abs(ED)` | 0.157% |

The bound is below exact diagonalization, as a minimization lower bound should
be. The size reduction is the real win: the solver sees tiny blocks instead of
a dense million-entry PSD object.

## [Known bottleneck](@id pauli-charge-bottleneck)

The remaining cost is not Mosek. In the ``N=16`` run, Mosek took only about
9.4 seconds out of a 2146.5 second `cs_nctssos` call.

The bottleneck is symbolic/JuMP construction: the current implementation still
builds and transforms dense symbolic moment data before emitting the tiny
reduced PSD blocks. The next engineering step is to build directly in the
charge/spatial reduced basis so the package never materializes the dense
intermediate representation.

That is a construction optimization, not a change to the mathematical
relaxation.

## See also

- [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for the general
  reduction pipeline and report fields.
- [Pauli Symmetry Reduction](@ref pauli-clifford-symmetry) for manual and
  SympleQ-detected Clifford reductions without charge-sector splitting.
- [`SymmetrySpec`](@ref), [`PauliChargeSectorSpec`](@ref),
  [`PauliSingletConstraintSpec`](@ref), and [`pauli_site_permutation`](@ref).
