# 0.1.1 replaces the public API and is not source-compatible with 0.1.0

## Declarative API and packed-word implementation

Despite the patch version, this release breaks the 0.1.0 API. The redesign now lives in `NCTSSoS`, not
`NCTSSoS.Next`; use `using NCTSSoS`. There is no `Next` compatibility alias.
Public qualified names, such as `NCTSSoS.build` and `NCTSSoS.Symmetry`, remain
qualified; they are not added to the exports.

### Entry-point migration

| 0.1.0 entry point | 0.1.1 entry point |
|---|---|
| `polyopt` + `SolverConfig` + `cs_nctssos` | `Problem` + `Relaxation` + `solve`; choose `MosekBackend`, `SCSBackend` or `JuMPBackend` |
| `create_pauli_variables` | `pauli_variables` |
| `create_fermionic_variables` | `fermionic_variables` |
| `create_bosonic_variables` | `bosonic_variables` |
| Other `create_*_variables` for free, projector and unipotent algebras | `site_variables(Free(), ...)`, `site_variables(Projector(), ...)`, `site_variables(Unipotent(), ...)` |
| `moment_eq_constraints` for a fixed sector | `Sector(Q, q)` in the relaxation |
| `cs_nctssos_higher` | `TermSplit(level=...)`; set the algorithm explicitly to reproduce its choice |
| `gns_reconstruct`, `verify_gns`, `reconstruct`, `test_flatness`, `flat_extend`, `robustness_report` | `gns(result)` returns reconstruction and diagnostics; raw-Hankel, clique and Cholesky forms removed (D4) |
| `particle_number_constraint` | `Sector(N, q)` in `Relaxation(strengthen=(...,))` |
| `newton_chip_basis`; `get_ncbasis`, `get_state_basis` | `NCTSSoS.newton_chip_basis(problem)`; `NCTSSoS.basis(registry, degree)` |
| `ς` / `varsigma`, `tr` | `state_symbol`, `trace_symbol` on `site_variables(...; functional=NCTSSoS.RealState())` or `NCTSSoS.RealTrace()`; Pauli/fermionic state symbols removed (D12) |
| `sympleq_symmetry_spec`, `sympleq_clifford_symmetry` | `detect_symmetry`; declare with `Problem(...; symmetry=...)` |
| `SignedPermutation`, `FermionicModePermutation`, `CliffordSymmetry`, `CliffordSymmetryGroup`, `Symmetry`, `SymmetrySpec`, `check_symmetry`, `pauli_*` generator helpers | Declare `SitePermutation`, `LetterPermutation`, `CliffordMap`, or `spin_rotations` with `Problem(...; symmetry=...)`; combine with `NCTSSoS.Symmetry`; validation occurs during construction/build |
| `PauliChargeSectorSpec`, `PauliSingletConstraintSpec` | Fixed sectors: `Sector(Q, q)` (singlets: `Sector(S², 0)`); continuous spin invariance: `spin_rotations(registry)`. Spin invariance alone does not impose a singlet sector |
| `heisenberg_chain_hamiltonian`, `heisenberg_chain_symmetry_spec`, `pauli_translation_invariant_*` model helpers | No dedicated model helpers; construct the Hamiltonian with `pauli_variables`, then declare symmetry or use `detect_symmetry` |
| `compute_sparsity`, `moment_relax`, `sos_dualize`, `solve_sdp` | `NCTSSoS.build`, `inspect`, `NCTSSoS.prepare`, `solve`; choose the conic form through the backend |
| `monomials`, `coefficients`, `terms` | `NCTSSoS.terms(p)` iterates `(coefficient, word)` pairs |
| `ComposedMonomial`, old algebra/registry/polynomial types | `NCTSSoS.mixed_variables` for tensor products; use the variable constructors above rather than old concrete types (custom names/subregistries removed, D11) |
| Old graph and sparsity helpers | Not public; select `Correlative` / `TermSplit` in `Relaxation` and use `inspect` for supported diagnostics |

`Positive` and `Zero` are operator constraints. `ExpectationPositive` and
`ExpectationZero` constrain one scalar expectation. `Sector(Q, q)` imposes a
sector identity and is not a substitute for an arbitrary scalar constraint.
Term-sparsity block choices can differ from 0.1 because the word order changed.

`MosekBackend(form=:auto)` selects SOS when there are no more moment variables
than cone coordinates; use `form=:moment` or `form=:sos` to select either form.

### Removed without replacement in 0.1.1

- **D3:** user-supplied derivations, SO(5), and derivation search.
- **D4:** GNS on cliques, from a raw Hankel matrix, and the Cholesky route.
- **D5:** fermionic spin adaptation and irrep labels.
- **D6:** structured and paired V2RDM input.
- **D9:** antiunitary detection and Gram-certificate recovery.
- **D10:** rational polynomial coefficients and polynomial ordering.
- **D11:** user-chosen variable names and subregistries.
- **D12:** formal state symbols over Pauli and fermionic polynomials.

**Not yet available:** D1, sparse state/trace routes (the dense routes remain);
D7, adaptive cone cuts (named strengthening cones remain).

### Known performance gaps, to repair

Accepted for this release; future repairs. Measurements compare old `main`
([ddef615](https://github.com/QuantumSOS/NCTSSoS.jl/commit/ddef615)) with the redesign
([6ec1e38](https://github.com/QuantumSOS/NCTSSoS.jl/commit/6ec1e38)), using Mosek 11.2.8:

- Chain40 (contiguous, order 4): solve time 21.5 → 25.6 s.
- Chain40: cones 44 → 114 (variables 19,142 → 6,607; largest side 30 → 16).
- Torus 3×2, first split level: 9.06 times fewer build allocations (budget: 10 times).
- H6: 2.42 times fewer build allocations.
- H8: 2.47 times fewer build allocations.
- H6: cones 12 → 14; largest side +2.
- H8: cones 12 → 14; largest side +2.
