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
  - `test/state_poly/` — state polynomial suite (needs `SOLVER`)
  - `test/correlated_sparsity/` — correlated/term sparsity suite (needs `SOLVER`)
  - `test/trace_poly/` — trace polynomial suite (needs `SOLVER`)
  - `test/v2rdm_structured/` — structured V2RDM benchmark (needs `SOLVER`)
  - `test/problems/` — curated problem tests (needs `SOLVER`)

## Solver

The curated `Pkg.test()` / CI suite uses COSMO (open-source).

The Next functional tests keep the S6 survey witness and program sizes in the
routine runner. Its bound check requires Mosek: COSMO reaches its iteration
limit with the original case settings. S7 still solves with COSMO. After
checking the Mosek licence with a trivial solve, run the licensed check in an
environment with the test dependencies installed:

```julia
using Mosek, NCTSSoS.Next
include("test/next/trace_poly/runtests.jl")
PublicFunctionalChecks.check(;sos_solver=MosekBackend(form=:sos))
```

This uses the same explicit `sos_solver=nothing` skip convention as the Next
problem pilot. It does not run the failed COSMO S6 solve or assert its raw
objective as a bound.

The Next E1–E5 Newton-chip benchmarks retain their old solver exception:
E2–E5 use Clarabel with `verbose=false`. E1 keeps its SOS bound check behind
the same licensed gate; its row and block-size checks run without a licence.
E5 at n=60 uses the reviewed interval `-1e-5 ≤ bound ≤ 1e-6`; all other
benchmark tolerances are unchanged. Run the licensed problem checks with:

```julia
using Mosek, NCTSSoS.Next
include("test/next/problems/runtests.jl")
NextProblemChecks.check(;sos_solver=MosekBackend(form=:sos))
```

The Next Heisenberg star checks (8-site dense, 10-site dense and MF grouping)
take more than 20 seconds per case. The routine runner reports explicit skips;
run them with the original COSMO settings through `references=true`:

```julia
include("test/next/problems/runtests.jl")
NextProblemChecks.check(;references=true)
```

After checking the Mosek licence, add `sos_solver=MosekBackend(form=:sos)`
to include the 8-site SOS check. The Majumdar–Ghosh bound, small spin and
singlet checks, and 16-site structure checks remain in the routine runner.

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

- `test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl`
  - local-only Mosek benchmark for the NC Motzkin polynomial (T4)
  - requires a local Mosek license
  - equivalent direct command:

```bash
julia --project test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl
```
