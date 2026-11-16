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
