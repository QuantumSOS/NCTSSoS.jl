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
  - `test/next/core/` — packed words and polynomial algebra
  - `test/next/model/` — public modeling API and symmetry detectors
  - `test/next/pipeline/` — relaxation construction and sparsity
  - `test/next/conic/` — standard forms and solver backends
  - `test/next/state_poly/` and `test/next/trace_poly/` — functional polynomials
  - `test/next/problems/` — reviewed problems and opt-in reference solves
  - `test/next/reference_data/` — frozen standard forms and detector signatures
  - `test/next/v2rdm_structured/` — ported V2RDM reference checks
  - `test/quality/` — Aqua, ExplicitImports, doctests

The `test/next` path is retained; these tests use the promoted `NCTSSoS` API.
`test/data` and `test/oracles` retain numeric fixtures and historical oracles.

## Solver

The curated `Pkg.test()` / CI suite uses COSMO (open-source).

The NCTSSoS functional tests keep the S6 survey witness and program sizes in the
routine runner. Its bound check requires Mosek: COSMO reaches its iteration
limit with the original case settings. S7 still solves with COSMO. After
checking the Mosek licence with a trivial solve, run the licensed check in an
environment with the test dependencies installed:

```julia
using Mosek, NCTSSoS
include("test/next/trace_poly/runtests.jl")
PublicFunctionalChecks.check(;sos_solver=MosekBackend(form=:sos))
```

This uses the same explicit `sos_solver=nothing` skip convention as the NCTSSoS
problem pilot. It does not run the failed COSMO S6 solve or assert its raw
objective as a bound.

The NCTSSoS E1–E5 Newton-chip benchmarks retain their old solver exception:
E2–E5 use Clarabel with `verbose=false`. E1 keeps its SOS bound check behind
the same licensed gate; its row and block-size checks run without a licence.
E5 at n=60 uses the reviewed interval `-1e-5 ≤ bound ≤ 1e-6`; all other
benchmark tolerances are unchanged. Run the licensed problem checks with:

```julia
using Mosek, NCTSSoS
include("test/next/problems/runtests.jl")
NextProblemChecks.check(;sos_solver=MosekBackend(form=:sos))
```

The NCTSSoS Heisenberg star checks (8-site dense, 10-site dense and MF grouping)
take more than 20 seconds per case. The routine runner reports explicit skips;
run them with the original COSMO settings through `references=true`:

```julia
include("test/next/problems/runtests.jl")
NextProblemChecks.check(;references=true)
```

After checking the Mosek licence, add `sos_solver=MosekBackend(form=:sos)`
to include the 8-site SOS check. The Majumdar–Ghosh bound, small spin and
singlet checks, and 16-site structure checks remain in the routine runner.

## Optional local reference checks

The old standalone files and their `test-local` / `test-local-one` Make targets
were removed. Their retained reference checks live in
`test/next/problems/local_references.jl`; the literal Bell table remains beside
them. Routine tests check Bell equations and builds, but skip the licensed solves.

After verifying the Mosek licence with a trivial solve, run the optional checks
in an environment with the test dependencies installed:

```julia
using Mosek
include("test/next/problems/local_references.jl")
LocalReferenceChecks.check(;references=true,mosek=true)
```

Both flags are required. This runs the retained Motzkin, NC polynomial,
Heisenberg, constrained Broyden and selected Bell reference solves. It is not
part of the routine CI budget. Sparse state/trace routes remain deferred;
the deleted Broyden trace script is not a runnable substitute.
