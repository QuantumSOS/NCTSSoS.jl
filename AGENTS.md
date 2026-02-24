# Repository Guidelines

## Project Layout
- `src/` — library code
  - `src/types/` — algebras, registries, monomials, polynomials
  - `src/simplification/` — algebra-specific rewrite rules (one file per algebra)
  - `src/optimization/` — sparsity + moment/SOS relaxations + JuMP model build
  - `src/states/` — state polynomials (quantum information)
- `test/` — curated suites (entry: `test/runtests.jl`)
- `docs/` — Documenter site (`docs/make.jl`, sources in `docs/src/`)
  - Literate sources: `docs/src/examples/literate/`
  - Generated (committed): `docs/src/examples/generated/`
- `test/oracles/` — optional reference outputs (requires external legacy `NCTSSOS` repo)

## Build, Test, and Docs Commands
CI baseline: Julia 1.11; solver: COSMO.
- `make init` — precompile root environment
- `make test` — full test suite (COSMO)
- `make coverage-ci` — CI-style coverage (`lcov.info`)
- `julia --project -e 'using Pkg; Pkg.test()'` — direct test run
- `make init-docs` — set up docs environment (`docs/Project.toml`)
- `make servedocs` — live docs server (LiveServer)
- `make examples` — execute Literate examples, regenerate `docs/src/examples/generated/`, and update `docs/examples_stamp.toml` (requires a Mosek license; CI only verifies the stamp via `julia docs/examples_stamp.jl`)
- `NCTSSOS_PATH=/path/to/NCTSSOS make oracle-chsh` — regenerate oracle values (see `Makefile` `oracle-%`)

## Architecture Overview
Type hierarchy:
```
AlgebraType
├── MonoidAlgebra        (NonCommutative/Projector/Unipotent)
├── TwistedGroupAlgebra  (Pauli)
└── PBWAlgebra           (Fermionic/Bosonic)
```
Core types: `NormalMonomial{A,T}` (immutable word), `Polynomial{A,T,C}` (mutable map), `VariableRegistry{A,T}` (symbol ↔ index).
Optimization flow: `polyopt()` → `cs_nctssos()` → `compute_sparsity()` → moment/SOS relaxation → JuMP model.

## Coding Style & Naming Conventions
- Indent: 4 spaces; no tabs.
- Naming: `CamelCase` types, `snake_case` functions, `UPPER_SNAKE_CASE` constants.
- Type params: keep algebra type `A` first for dispatch (e.g. `Polynomial{A,T,C}`).
- Prefer surgical changes; add regression tests for bug fixes.

## Testing Guidelines
- Canonical instructions: `TESTING.md`.
- Shared infra: `test/TestUtils.jl` defines `SOLVER` (COSMO) and helpers.
- Suites: `test/polynomials/` (no solver), `test/quality/` (Aqua/ExplicitImports/doctests), `test/relaxations/`, `test/correlated_sparsity/`, `test/problems/`.
- Keep tests deterministic and solver-stable (COSMO in CI).

## Commit & Pull Request Guidelines
- Commit messages (recent): `feat: ...`, `fix(scope): ...`, `docs: ...`, `ci: ...` (often with `(#PR)` on merge/squash).
- PRs: small diffs, clear description, and a short “How to test” section (e.g., `make test`). Note any solver requirements (e.g., Mosek-only example regeneration).
