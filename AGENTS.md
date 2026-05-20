# Repository Guidelines

## 1. Current task

See [`TASK.md`](TASK.md). Read it first; it sets the immediate goal and any
rules for how to iterate on it.

## 2. Project layout (minimal)

- `src/` — library code.
  - `src/types/` — algebras, registries, monomials, polynomials.
  - `src/simplification/` — algebra-specific rewrite rules (one file per algebra).
  - `src/optimization/` — sparsity + moment/SOS relaxations + JuMP model build.
  - `src/states/` — state polynomials.
- `test/` — curated suites (entry: `test/runtests.jl`).
- `docs/` — Documenter site (Literate sources under
  `docs/src/examples/literate/`, generated pages under
  `docs/src/examples/generated/`).

Type hierarchy:

```
AlgebraType
├── MonoidAlgebra        (NonCommutative / Projector / Unipotent)
├── TwistedGroupAlgebra  (Pauli)
└── PBWAlgebra           (Fermionic / Bosonic)
```

Optimization flow: `polyopt()` → `cs_nctssos()` → `compute_sparsity()` →
moment / SOS relaxation → JuMP model.

## 3. Build and test

CI baseline: Julia 1.11; solver: COSMO.

- `make init` — precompile root environment.
- `make test` — full test suite (COSMO).
- `julia --project -e 'using Pkg; Pkg.test()'` — direct test run.
- `make init-docs` then `make servedocs` — set up and preview the docs site.
- `make examples` — regenerate Literate examples (Mosek-only; CI verifies the
  stamp via `julia docs/examples_stamp.jl`). Run after any change under
  `docs/` or any docs-facing content.
