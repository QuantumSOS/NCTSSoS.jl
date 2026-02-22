# Testing

Single source of truth for running the test suite.

## Quick start (preferred: Makefile)

- `make init` — precompile (first run)
- `make test-ci` — CI default: `--polynomials` + `--minimal` (COSMO; fast)
- `make coverage-ci` — CI/Codecov-style coverage (`lcov.info`)
- `make test-minimal` — smoke suite (5 critical paths)
- `make test-polynomials` — algebra only (no solver)
- `make test-relaxations` — relaxation components
- `make test-quality` — Aqua + ExplicitImports + doctests
- `make test` — full suite with `--local` (Mosek; slow)

## Direct Julia (`Pkg.test`)

- Minimal (fast): `julia --project -e 'using Pkg; Pkg.test(test_args=["--minimal"])'`
- CI default: `julia --project -e 'using Pkg; Pkg.test(test_args=["--polynomials","--minimal"])'`
- Full suite (Mosek): `julia --project -e 'using Pkg; Pkg.test(test_args=["--local"])'`

## Flags

Passed via `Pkg.test(test_args=[...])`:

- `--minimal` — fast correctness check (5 paths)
- `--polynomials` — core algebra (no solver)
- `--quality` — quality checks
- `--relaxations` — relaxation components
- `--problems` — problem suites
- `--local` — use Mosek instead of COSMO

## Defaults (no flags)

- Local `Pkg.test()` (no flags): run all test groups; some heavy problem suites only run when `--local` is enabled.
- CI `Pkg.test()` (no flags, `ENV["CI"]=true`): `--polynomials` + `--minimal`.

## Solvers

- Default (open-source): COSMO
- `--local`: Mosek (requires license)

## Test Layout (and `SOLVER` / `USE_LOCAL`)

- Entry point: `test/runtests.jl` (parses flags; runs groups).
- Shared infra: `test/TestUtils.jl` defines `SOLVER`, `SOLVER_NAME`, `USE_LOCAL`.
- Groups:
  - `test/polynomials/` (no solver; no JuMP)
  - `test/relaxations/` (needs `SOLVER`)
  - `test/problems/` (needs `SOLVER`; some suites are heavy)
  - `test/quality/` (Aqua, ExplicitImports, doctests)
- Local-only suites are gated by `USE_LOCAL` in the runners; heavy suites may still run on COSMO but expect it to be slower.

### Local-Only / Recommended `--local`

- `test/problems/condensed_matter/` (precision/perf)
- `test/problems/quantum_networks/` (bilocal network needs `order=3`)
- `test/problems/fermionic/` (precision/perf)
- Some large sweeps (e.g. Bell inequalities, large-scale NC polynomial) are much faster with Mosek.

## Single-file runs

- `make test-file FILE=path/to/file.jl`
