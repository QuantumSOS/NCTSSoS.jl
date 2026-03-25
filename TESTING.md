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

The curated `Pkg.test()` / CI suite uses COSMO (open-source).

## Standalone Local-Only Benchmark Scripts

Some heavier literature reproductions live under `test/problems/` but are not
included from `test/runtests.jl`.

- Run all registered local-only scripts:

```bash
make test-local
```

- Run one local-only script directly:

```bash
make test-local-one SCRIPT=path/to/script.jl
```

- To add another local-only script later, register its path in
  `LOCAL_ONLY_TEST_SCRIPTS` in the `Makefile`.

Currently registered:

- `test/problems/trace_polynomial/t1_broyden_banded_trace.jl`
  - local-only Mosek benchmark
  - runs the heavier T1 trace family at `n = 20, 40`
  - requires a local Mosek license
  - equivalent direct command:

```bash
julia --project test/problems/trace_polynomial/t1_broyden_banded_trace.jl
```
