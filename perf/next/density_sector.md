# DensityMatrix on its sector face

The four-site Ising control now solves to zero with Mosek and SCS. Its density
cone changes from side 16 to side 6. The original formulation returned
`SLOW_PROGRESS` without a bound with Mosek.

## Exact reduction and its boundary

For a positive density matrix, `tr((Q-q)²ρ)=0` implies `(Q-q)ρ=0`.
The pipeline finds the common kernel of the contained sector operators by
exact complex-rational elimination, with no numerical rank tolerance. It
adds the complete off-face equations and replaces the density cone by
`B'ρB ⪰ 0`, where B has full column rank and spans that kernel. These two
conditions are equivalent to the original density positivity and sector
variance. Dropping the off-face equations would be incorrect.

Containment is tested against the density's operator basis, including phases,
not its original site list. This also handles Clifford-transformed algebras.
A global sector does not restrict a proper subsystem unless its operator
belongs to that subsystem's algebra. For example, either one-site density
of the zero-magnetization Bell state is full rank. Multiple contained sectors
use their common kernel. An empty common kernel raises a named `ArgumentError`.
The strengthening report exposes `density_faces` with source, before/after
sides and sector indices; the size display prints the reduction.

Production keeps the nonorthogonal exact basis. The physical-spectrum check
orthogonalizes it only for comparison; congruence, not identical eigenvalues,
is the PSD invariant for a nonorthogonal basis. The test also compares the
emitted density cone's spectrum with the nonorthogonal compression directly.

## Solver and physical-state evidence

Command: `julia --compiled-modules=existing --project perf/next/density_sector.jl check`
with Mosek and SCS available in a separate environment through JULIA_LOAD_PATH.
Julia, BLAS, OpenMP and Mosek each use one thread. No dependency was added.
The Julia wrappers are Mosek 11.2.1 and SCS 2.7.1.

Initial isolated checks (before packaging the reusable acceptance script):

| Case | Old bound | Face bound | Solver |
| --- | ---: | ---: | --- |
| Four-site Ising, order 1, Correlative + TermSplit, antiunitary moments, Q=ΣZ, q=0 | no bound (`SLOW_PROGRESS`) | −5.38196e−10 | Mosek |
| Same Ising control | not measured | −4.94623e−10 | SCS |
| Two-site complex sector Q=3Y₁+4Z₁, q=5 | −1.0071742653917521 | −1.0071742656776286 | Mosek |
| Four-site Heisenberg, S²=0 | −8.000000000000012 | −8.000000000455767 | Mosek |
| CNOT-transformed one-site density algebra, Z₁=1 | not measured | −1.74128e−9 | Mosek |

Both old/new comparisons differ by less than 4.6e−10. The spectra at five
random states per sector agree after orthonormal compression within 1e−12;
all emitted cones and equalities accept those states. The kept regression
also rejects an off-face perturbation with unchanged compressed density,
checks a nonlocal-sector exclusion, trivial/joint/empty sectors, and the
bosonic error's recommendation to use Sector.

The first Ising check took 59.27 s (estimate 30–60 s). The separate complex,
singlet and CNOT checks took 52.55 s (estimate 1–2 min). The sector test file
passed 98 assertions in 141.83 s (estimate 1–2 min), before adding the bosonic
error assertion. The packaged acceptance script passed 55 assertions in
72.61 s on the rebased tree.

## Combination rows and unaffected fingerprints

The refresh checks all 4,864 density+sector rows: 1,900 programs and 2,964
named errors, with no defects. Exactly 1,900 rows change; the 79,651 rows
outside this pair are byte-identical. The selected Ising row now reports
`OPTIMAL/OPTIMAL` with SCS instead of `ITERATION_LIMIT/ITERATION_LIMIT`.
Its reduced/reference bound difference is 1.12e−10; all 257 queries succeed.

The 24-case pilot took 123.70 s. The full refresh took 682.98 s, below the
30–45 minute pilot estimate. The selected solve and queries took 24.33 s,
below the 1–2 minute estimate. Reproduce this subset with:

```sh
julia --project test/next/pipeline/run_combination_matrix.jl --full --strengthening=density+sector --solve density-sector.tsv
```

The existing 95 smoke fingerprints were recomputed: the 94 programs without
this pair are bit-identical; only the expected XXZ density+sector program
changes. This includes 78 unchanged programs from the older 79-case set and
16 bosonic programs. Fingerprint comparison includes all standard-form
dimensions, CSC arrays, right-hand sides, objective, offset and cone labels.

## Combined suite on the base containing #442

The Next suite, native Mosek/SCS checks, Aqua and ExplicitImports passed:
55,887 printed assertions, including 147 Mosek construction/solution checks,
101 SCS writer checks, and Aqua 11/11. ExplicitImports was clean. The suite
took approximately 1,105 s; acceptance, fingerprints and suite together stayed
within the estimated 20–25 minutes. The base was `504e3a2` on `redesign`.

```sh
julia --compiled-modules=existing --project -e 'include("test/next/conic/mosek.jl"); include("test/next/runtests.jl"); include("test/next/conic/extensions.jl"); include("test/quality/Aqua.jl"); include("test/quality/ExplicitImports.jl")'
```

## Exact S² kernel cost through eight sites

Command: `julia --compiled-modules=existing --project perf/next/density_sector.jl kernel`.
Julia 1.12.7, `icelake-server`, one Julia/BLAS thread. Input matrices are built
outside the timer. Each exact kernel call is warmed once, then sampled three
times with GC outside the timer. The exact matrix product S²B is zero and the
kernel dimension equals `binomial(n,n/2)-binomial(n,n/2-1)` at every size.

| Sites | Matrix side | Kernel dimension | Median seconds | Range seconds | Median allocated bytes |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 16 | 2 | 0.000688 | 0.000385–0.001303 | 459,192 |
| 6 | 64 | 5 | 0.008323 | 0.007940–0.009310 | 12,092,720 |
| 8 | 256 | 14 | 0.551286 | 0.484027–1.048805 | 427,786,456 |

The isolated measurement took about 23 s (estimate 1–2 min). These are kernel
costs, not full density-assembly costs or peak memory. There is no existing
DensityMatrix size limit. This measurement does not justify a new kernel cap
through eight sites; larger sizes and arbitrary coefficient bit lengths have
not been measured. The existing full density basis can itself be expensive.

## Forced kernels that remain outside this change

- A moment block containing the linear combination Q−q has that vector in
  its kernel for a state in the sector, since its quadratic form is the
  sector variance. This requires the relevant rows to be present.
- The fermionic G condition includes identity and particle-hole rows. A
  number sector on its modes gives the kernel vector Σa†ᵢaᵢ−N·1.

Those cones are unchanged. Removing density singularity does not claim to
make every relaxation strictly feasible. Sector allocation optimization and
continuous-reduction costs are separate work.
