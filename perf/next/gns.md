# GNS reconstruction from a solved result

`gns(result)` selects a covered adjacent pair of Hankel degrees and returns one
`Next.GNSResult`. Use `g[x]`, `evaluate(g, p)` and `expectation(g, p)` with the
variables held by the caller. No registry, moment map or degree is repeated.
Neither word contract nor conic contract changed.

## Three numerical verdicts

- `:moments_reproduced`: flat, positive Hankel extension; the returned model
  reproduces every covered program moment and passes the feasibility checks.
- `:optimizer_certified`: the matrices pass algebra/adjoint, state norm,
  constraint, sector and objective checks, but do not reproduce all moments.
  `expectation(g, p)` is the optimizer's expectation and can differ from
  `expectation(result, p)`. This is a numerical check, not a rigorous certificate.
- `:withheld`: neither test passed. Operators and state are `nothing`; indexing
  and evaluation throw `ArgumentError` with a remedy. A flatness failure is never
  silently changed into a flat extension.

The display reports the verdict, largest program-moment error, ranks, selected
degrees and whether the pair exceeds the relaxation order. The result also
stores Hankel-only error, constraint minimum eigenvalues, equality norms, sector
residuals, objective error, distance to flatness, spectra and conditioning.

## Nine solver runs, with unchanged GNS defaults

Measured on 2026-11-15 in a Linux x86-64 orb, `icelake-server`, Julia 1.12.7,
SCS.jl 2.7.1 (SCS 3.3.1), Mosek.jl 11.2.1. Julia and BLAS use one thread.
No CPU affinity was set. SCS uses the native direct writer, `max_iters=200000`,
and equal absolute/relative tolerances. Mosek uses `MosekBackend()` defaults.
A trivial Mosek optimization checks the licence before measurement.

Reproduce in an environment with NCTSSoS, SCS, Mosek and a valid Mosek licence:

```sh
julia --project=<solver-environment> perf/next/gns.jl
```

The script prints complete spectra and residuals. It tests CHSH at order 4,
the free-algebra problem `X*Y*X` with `1-X^2 >= 0`, `1-Y^2 >= 0` at order 3,
and the four-site XXZ ring with anisotropy 0.5 at order 2. All nine solves
returned `OPTIMAL`; all nine reconstructions returned operators without keyword
overrides. Every selected pair was (full degree 3, principal degree 2).

| Case | Solver tolerance | Ranks principal/full | Max Hankel error | Max program error | Objective error | Verdict |
|---|---|---:|---:|---:|---:|---|
| CHSH | SCS 1e-8 | 4/4 | 2.70e-9 | 7.43e-9 | 8.45e-10 | Moments reproduced |
| CHSH | SCS 1e-6 | 4/4 | 2.70e-9 | 7.43e-9 | 8.45e-10 | Moments reproduced |
| CHSH | Mosek default | 4/4 | 3.27e-10 | 1.00e-9 | 2.11e-9 | Moments reproduced |
| Free | SCS 1e-8 | 5/9 | 1.93e-2 | 1.93e-2 | 1.49e-9 | Optimizer certified |
| Free | SCS 1e-6 | 5/9 | 1.93e-2 | 1.93e-2 | 1.49e-9 | Optimizer certified |
| Free | Mosek default | 5/9 | 1.47e-1 | 1.47e-1 | 1.25e-9 | Optimizer certified |
| XXZ | SCS 1e-8 | 16/16 | 2.82e-12 | 2.82e-12 | 9.62e-13 | Moments reproduced |
| XXZ | SCS 1e-6 | 16/16 | 1.61e-8 | 1.61e-8 | 9.17e-9 | Moments reproduced |
| XXZ | Mosek default | 16/16 | 9.41e-10 | 9.41e-10 | 1.79e-8 | Moments reproduced |

For Free, the smallest constraint eigenvalues are `(7.91e-10, 1.48e-9)` with
SCS and `(6.31e-10, 6.98e-12)` with Mosek. Thus these are feasible numerical
optimizers, not full moment reconstructions. No tolerance was loosened for Mosek.
The two SCS settings stop at the same iterate on CHSH and Free, hence equal rows.

The XXZ program-moment error is relative to the **solved moments**, not exact
diagonalization. Small reconstruction error does not improve inaccurate solver
moments. In particular, the Mosek moment-accuracy caveat in `mosek.md` still
applies. Native SCS at 1e-8 is the acceptance reference for the Pauli relations.

### Rank gaps and positivity

Each gap below is `last retained / first discarded` (two singular values, not
their ratio). Default `atol=nothing` selects the largest ratio above 1000 only
when its discarded tail is at most `check_atol * max(1, σmax)`. The cutoff is ten
times the discarded tail, floored at roundoff and capped below the retained
value; `rtol * σmax` also applies. Otherwise the absolute cutoff is 1e-8.
`check_atol=1e-6`, `rtol=0` and `psd_rtol=1e-8` remain fixed for all cases.

| Case / solver | Principal gap | Full gap | Principal / full cutoff | Minimum full-Hankel eigenvalue |
|---|---|---|---|---:|
| CHSH / SCS 1e-8 | 2.00 / 4.64e-9 | 3.00 / 4.88e-9 | 4.64e-8 / 4.88e-8 | −2.20e-9 |
| CHSH / SCS 1e-6 | 2.00 / 4.64e-9 | 3.00 / 4.88e-9 | 4.64e-8 / 4.88e-8 | −2.20e-9 |
| CHSH / Mosek | 2.00 / 1.37e-9 | 3.00 / 1.82e-9 | 1.37e-8 / 1.82e-8 | 2.07e-10 |
| Free / SCS 1e-8 | 3.20e-2 / 1.80e-9 | 1.92e-5 / 2.42e-9 | 1.80e-8 / 2.42e-8 | −3.86e-10 |
| Free / SCS 1e-6 | 3.20e-2 / 1.80e-9 | 1.92e-5 / 2.42e-9 | 1.80e-8 / 2.42e-8 | −3.86e-10 |
| Free / Mosek | 1.33e-1 / 8.06e-10 | 2.70e-2 / 1.05e-9 | 8.06e-9 / 1.05e-8 | −1.56e-10 |
| XXZ / SCS 1e-8 | 4.00 / 5.46e-12 | 7.00 / 1.76e-11 | 5.46e-11 / 1.76e-10 | −1.76e-11 |
| XXZ / SCS 1e-6 | 4.00 / 3.85e-8 | 7.00 / 1.15e-7 | 3.85e-7 / 1.15e-6 | −1.15e-7 |
| XXZ / Mosek | 4.00 / 3.11e-9 | 7.00 / 5.90e-9 | 3.11e-8 / 5.90e-8 | 6.28e-11 |

XXZ uses a pair above order 2 because all needed moments are covered. Its full
matrix must independently pass the PSD test: for SCS 1e-6 the allowed negativity
is 1.200000008e-7, just above the measured 1.146523946e-7. Tests also supply an
order-one-positive functional whose higher Hankel is indefinite: that pair is
rejected and its operators withheld, even with an `OPTIMAL` status.

## Packed fill and complete reconstruction cost

Same XXZ SCS 1e-8 moments for both implementations, `ComplexF64`, principal
degree 2/full degree 3, full Hankel side 175. Legacy `gns_reconstruct` gets the
same moment map and `hankel_deg=2, atol=1e-6`; new `gns` uses automatic defaults.
Both ranks are 16, principal singular values differ by at most 7.11e-15, and
the old model's largest error over all moments is 2.83e-12.

Five warmed samples per method, explicit GC before each sample, compilation and
solver work outside the timer. Reported bytes are Julia allocated bytes, not
peak RSS. The new fill reuses an output matrix and a prefilled packed-word value
cache; its preflight, cache construction and matrix allocation are outside that
kernel timer. Legacy `hankel_matrix` allocates its output and temporary words.
The complete reconstruction measurements include these setup costs, automatic
degree search and all diagnostics in the new API.

| Region | Median ms | IQR ms | Min–max ms | Allocated bytes |
|---|---:|---:|---:|---:|
| New complete `gns` | 20.134 | 19.742–20.168 | 16.675–22.344 | 20,518,080 |
| Old `gns_reconstruct` | 28.980 | 28.524–29.077 | 26.714–34.500 | 38,440,152 |
| Packed fill, reusable buffers | 1.375 | 1.339–1.376 | 1.194–1.419 | 0 |
| Old `hankel_matrix` | 5.548 | 5.383–5.565 | 5.304–5.919 | 12,208,848 |

Complete reconstruction allocated 46.6% fewer bytes on this case. Timing varied
across runs in the shared orb, with overlapping whole-reconstruction ranges in
an earlier run; this is not a general speedup claim. The separate two-mode
fermionic fill (side 15, including contracting products) also allocated 0 bytes
inside a typed kernel call. An untyped measurement at the harness boundary
allocated 368 bytes; that dispatch cost is not a per-entry allocation.

The final benchmark, including 652 exact-moment assertions and nine fresh solves,
finished in 74.89 s (estimated 60–120 s). Two earlier attempts took 95.78 s and
83.37 s: the first exposed missing adjoint aliases in the added fermionic harness,
and the second measured its untyped call boundary. Both were harness issues;
the final script uses the same coverage preflight as reconstruction.

## Coverage and old API mapping

Values and coordinate coverage are read only through `moment_terms` and
`moment_covered`. A complete row set with discrete and continuous symmetry is
tested. Partially covered real/imaginary components are checked separately;
coordinates may cancel exactly in a complete product. Missing required moments
produce an `ArgumentError` naming split/custom rows/row rules and the remedy.
Sparse grouping and custom families are not supported in this delivery.

| Old entry point | New form |
|---|---|
| `gns_reconstruct`, `reconstruct` | `gns(result)`; one result, no matrices-only wrapper |
| `hankel_matrix`, `hankel_deg` | Internal packed fill and automatic covered degree selection |
| `verify_gns`, `robustness_report` | Always-computed `g.checks`, `g.flatness` and display |
| `test_flatness` | `g.flatness.is_flat` with cutoffs, spectra and residuals |
| `flat_extend` | Dropped as a public step; do not silently alter solved moments |
| `method=:svd` | SVD of the principal block; Hermitian eigenvalues for full spectrum/PSD |
| `method=:cholesky` | Dropped; one rank-revealing reconstruction path |
| Sparse disjoint cliques | Deferred; named `ArgumentError`, not partial matrices |

The old implementation remains unchanged as the oracle. Tests compare ranks,
principal singular values and all Pauli words on identical moments, not matrix
entries (which may differ by a unitary). They also check two-mode CAR with
complex moments, CHSH squares and party commutation, Pauli relations, invalid
constraints, false objective bounds, tolerance remedies, and high-degree
program moments outside the selected Hankel matrix.

## Local verification

The combined driver included `test/TestUtils.jl` (COSMO),
`test/next/runtests.jl`, both native extension test files in separate modules,
`test/quality/Aqua.jl`, `test/quality/ExplicitImports.jl`, and the docs stamp
verification. It passed 43,650 assertions, including the exact fixtures in each
test module, with Aqua 11/11 and ExplicitImports clean. The nine native GNS
cases contributed 333 assertions; coverage/state restrictions contributed 8.
Elapsed time was 857.39 s (estimated 12–16 minutes), including cache rebuilds.

Commands used in the prepared solver environment:

```sh
julia --project=.amp/in/orb-backup/solverenv .amp/in/orb-backup/letters-workload-gate.jl
julia --compiled-modules=existing --project=.amp/in/orb-backup/solverenv perf/next/gns.jl
julia --project=.amp/in/orb-backup/solverenv docs/examples_stamp.jl
git diff --check
```

The local driver is only a wrapper for the committed test entry points above.
No legacy implementation, solver dependency, documentation example or contract
was changed. Sparse reconstruction and the later algebra work packages remain
outside this delivery.
