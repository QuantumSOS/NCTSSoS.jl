# Repository Guidelines

If `TASK.md` exists at the repo root, read it after `AGENTS.md` for the current task scope and links into `plan/`.

## Project Layout
- `src/` — library code
  - `src/types/` — algebras, registries, monomials, polynomials
  - `src/simplification/` — algebra-specific rewrite rules
  - `src/optimization/` — sparsity + moment/SOS relaxations + JuMP model build
  - `src/states/` — state polynomials (quantum information)
- `test/` — curated suites (entry: `test/runtests.jl`)
  - `test/data/` — reviewed expectation fixtures
  - `test/TestUtils.jl` — shared infra; defines `SOLVER` (COSMO)
- `docs/` — Documenter site (`docs/make.jl`, sources in `docs/src/`)
  - `docs/src/examples/literate/` — Literate sources
  - `docs/src/examples/generated/` — generated (committed)

## Build, Test, and Docs Commands
CI baseline: Julia 1.11; solver: COSMO.

- `make init` — precompile root environment
- `make test` — full test suite (COSMO)
- `make coverage-ci` — CI-style coverage (`lcov.info`)
- `julia --project -e 'using Pkg; Pkg.test()'` — direct test run
- `make init-docs` — set up docs environment (`docs/Project.toml`)
- `make servedocs` — live docs server (LiveServer)
- `make examples` — execute Literate examples, regenerate `docs/src/examples/generated/`, update `docs/examples_stamp.toml` (requires Mosek)

After any change under `docs/`, run `make examples` then `make servedocs` to preview before handoff.
