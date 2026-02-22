# Testing

Single source of truth for running the test suite.

## Commands

- `make test` — full test suite (COSMO)
- `make coverage-ci` — CI coverage (`lcov.info`)

## Direct Julia (`Pkg.test`)

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Test Layout

- Entry point: `test/runtests.jl` (fixed curated suite, no flags)
- Shared infra: `test/TestUtils.jl` defines `SOLVER` (COSMO) and `flatten_sizes`
- Groups:
  - `test/polynomials/` — polynomial algebra (no solver, no JuMP)
  - `test/quality/` — Aqua, ExplicitImports, doctests
  - `test/relaxations/` — relaxation components (needs `SOLVER`)
  - `test/problems/` — curated problem tests (needs `SOLVER`)

## Solver

All tests use COSMO (open-source).
