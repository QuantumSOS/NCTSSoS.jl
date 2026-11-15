# Symmetry-first sparsity experiment

This is a standalone prototype using this branch's detector and
generator-only decomposition. It does not change the production sparse API.

## Branch handoff

This branch preserves the research state for discussion. It is not a production
integration into `redesign`. The current question is why B's global moment
matrix excludes the low-energy assignments admitted by E5. Do not start E7
without a new decision. See the B structure diagnostic at the end of this file.

The scripts and evidence now live in `perf/symmetry-sparsity/`. Historical
commands below use their original `.amp/in/orb-backup/symmetry-sparsity/` path;
replace that prefix with `perf/symmetry-sparsity/` when running this branch.
Relative Julia includes continue to work. Runs overwrite their result files,
so preserve the checked-in evidence or use a separate worktree.

`environment-root.toml` and `environment-docs.toml` record the two manifests.
To restore the measured environment, copy them to `Manifest.toml` and
`docs/Manifest.toml`, respectively, then instantiate both Julia projects.
Their local package paths refer to those original destinations. A valid Mosek
license is required; no license or other credentials are included.

Earlier full-suite evidence is `optimized-full-tests.log`: 12,253 passed,
two broken, 1,307.20 seconds. It predates the latest standalone experiments.
The later prototype and projection checks are recorded below. The handoff
does not claim a new full-suite or benchmark run.

`experiment.jl` supplies signed row actions, symmetry closure of local blocks,
block-orbit transporters and Schreier stabilizers, signed moment coordinates,
level-1 term graphs with invariant chordal completion, and SDP assembly.
No group elements are enumerated by this module. `validate.jl` enumerates the
two-qubit group only as an independent small-system reference.

Run from the repository root with the existing docs environment and a valid
Mosek license. Validation first performs a trivial license-backed solve.

```sh
JULIA_LOAD_PATH="@:$PWD/docs:$PWD:@stdlib" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=4 \
  python3 .amp/in/orb-backup/symmetry-sparsity/limited_run.py \
  julia --project=docs .amp/in/orb-backup/symmetry-sparsity/validate.jl
```

The process guard enforces 600 seconds and 16 GiB sampled process RSS. Internal
budgets cap item counts (including distinct moment words) at 100,000 and stored
group entries/graph edges at 5,000,000. They throw rather than return a partially
closed relaxation. These limits also apply to the initial 20-site feasibility
comparison. A cap is not proof that a route is mathematically infeasible.

## Validation result

The final run passed 722 assertions in 129.184 seconds including startup,
with 2.111 GiB peak sampled RSS. See `validation-vector.log` and `validation.toml`.
Reduced and unreduced D/E/F solutions agreed within 1.71e-8 across nine cases.
Checks include exact Pauli phases against matrices, all 11,520 two-qubit actions,
the 192-action Heisenberg stabilizer, block stabilizers, sign-forced zeros,
an asymmetric complex objective, shared moments, a symmetry-breaking chordal
completion, incomplete-input rejection, and resource caps.

The four-site degree-2/two-site-window diagnostic is deliberately sparse:
E and F give about -3, while D gives -2, equal to the direct ground energy.
With degree 4 and four-site windows, all three give about -2. This is not a
20-site result or evidence that their 20-site bounds coincide.

## Primal and dual formulations

The affine-matrix primal is available with `dual_form=false`. The large
unreduced four-site reference exceeded the memory guard in Mosek's factorization.
`dual_form=true` builds the explicit conic dual with one coefficient equation
per independent moment. Both reduced and unreduced comparisons use this same
form. It recovers primal moments from equality duals, verifies every original
retained PSD block, and checks primal-dual gaps. Two-site tests also compare
against the affine-matrix primal. No solver tolerance was loosened.

The realification doubles the numerical PSD side: a complex Hermitian block of
side m is represented by a real symmetric block of side 2m. Report both when
benchmarking; `sizes` holds the Hermitian side.

## Stage-separated 20-site comparison

`benchmark.jl` takes a route letter A-F, site count, and repetition count.
Use one fresh process per route. Four repetitions mean one first execution and
three warm executions. No concurrent benchmarks were used in the final batch.

```sh
export JULIA_LOAD_PATH="@:$PWD/docs:$PWD:@stdlib"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=4
python3 .amp/in/orb-backup/symmetry-sparsity/limited_run.py \
  julia --project=docs .amp/in/orb-backup/symmetry-sparsity/benchmark.jl E 20 4
```

The driver records stage seconds, allocations, block sizes, energy bounds and
verification residuals in `benchmark-<route>-20.toml`; the matching log records
whole-process time and peak sampled RSS on normal exit. A guard termination
leaves the last active stage in TOML and `LIMIT_EXCEEDED` in the log. Stage sums
exclude loading, the license check, between-sample GC, and report writes.
Allocated bytes are cumulative allocation, not peak resident memory.

All routes use the periodic chain with coupling 1/4 and contiguous degree-4
seed rows. A uses no symmetry; B uses four supplied lattice/axis generators;
C-F use all detected generators. D/F use term-sparsity level 1 with a fixed
active support. E/F start from four-site windows. Mosek uses four threads,
1e-8 feasibility/gap tolerances and a 300-second solver limit. The realification
and explicit dual formulation are the same for all routes.

The final batch uses sorted vector block IDs. Earlier tuple IDs encoded the
number of rows in the Julia type and caused excessive compilation. Files named
`before-vector` are diagnostic history, not the final comparison.

## Final feasibility results

Environment: Julia 1.12.7, Intel Xeon @ 2.60 GHz, one Julia thread and four
BLAS/Mosek threads. The included manifests fix package versions. Local source
is based on 7dc71e3fc73927d6208d1527a284ef7a76689380 with uncommitted changes;
the evidence archive includes that source, not only its base commit.

| Route | Outcome | Whole-process time | Peak sampled RSS |
| --- | --- | ---: | ---: |
| A | 1,910,941 moment words exceeded the 100,000 cap | 50.581 s | 1.245 GiB |
| B | 1,910,941 moment words exceeded the 100,000 cap | 59.060 s | 1.161 GiB |
| C | 7,304,104 moment words exceeded the 100,000 cap | 173.691 s | 3.248 GiB |
| D | Process guard stopped sparsity construction | 600 s guard | Not recorded on guard exit |
| E | Four verified executions | 192.705 s | 1.266 GiB |
| F | One verified first execution | 403.465 s | 1.519 GiB |

Moment counts are before symmetry identifies equivalent words. The item cap
is checked after each parent block's products are generated, so a count can
exceed the cap. A-D have no solved bound in this batch. Do not interpret their
stopping times as times to solution or their caps as inherent route limits.

E's first measured-stage sum was 69.823 s; three warm sums were 25.821, 28.726,
26.110 s (median 26.110 s). F's first measured-stage sum was 356.526 s, including
299.998 s in sparsity construction. Its whole process took 6.72 minutes, within
the 6-8 minute estimate; E's process took 3.21 minutes, within the 3-4 estimate.
Stage sums include verification but exclude startup and license-check overhead.

E keeps 160 parent blocks of side 175, with one block orbit, 12,481 moment words,
nine independent moments, and Hermitian blocks [2,2,3,5,8,8,9,11,11,14]. F has
300 blocks before reduction, four block orbits, 961 moment words, nine independent
moments, and Hermitian blocks [1,1,1,1,2,2,2,2,3,3,4,4,5,5,6,7,8,11,11,14].
Both retain all 67 detected generators; neither enumerates the full group.

| Route | Numerical lower bound | Primal-dual gap | Minimum original-block eigenvalue | Symmetry residual |
| --- | ---: | ---: | ---: | ---: |
| E | -9.999999988364584 | 1.2881e-11 | 5.9902e-10 | 0 |
| F | -9.999999993535125 | 6.9616e-12 | 1.5635e-10 | 0 |

These are floating-point solver results, not rigorous interval certificates.
E was the practical candidate in that capped batch. F's bound differs by only
5.18e-9 and has no observed improvement. The follow-up below resolves B's bound.

## B accuracy follow-up

B's rows, four generators, normalization, objective and PSD constraints are
unchanged. Raise the moment cap to 2,000,000. Project sparse numerical
coefficient matrices K[k] in M = sum(K[k]*y[k]), rather than multiply a dense
matrix of symbolic expressions. The old assembly remains the default reference.
Projection discards real/imaginary coefficients below 1e-12 and checks that the
sum discarded in any reduced entry is below 1e-8. The recovered moments are
checked against the original full 2,401-row matrix, not only the reduced blocks.

An initial attempt was interrupted after 588.094 s: site permutations used
Int64 while the rows used UInt8, causing a full Clifford conversion and
validation on each action. `index-type-oracle.log` verifies all 4,802 translation
and reflection images agree after constructing these permutations with UInt8;
the oracle took 12.58 s. This fixes the harness, not the relaxation.

The 722 existing checks passed after the assembly addition. The seven added
projection checks passed in 67.576 s; they compare every coefficient with the
old complex symbolic projection and compare four-site solves and full PSD
matrices. The four-site bounds differed by about 5e-14. Initial test-harness
errors and the interrupted run remain in separate diagnostic logs.

```sh
export JULIA_LOAD_PATH="@:$PWD/docs:$PWD:@stdlib"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=4
MOMENT_LIMIT=2000000 NUMERIC_PROJECTION=true REPORT_SUFFIX=-solved \
RUN_SECONDS=1800 ALLOW_FEASIBLE=true \
  python3 .amp/in/orb-backup/symmetry-sparsity/limited_run.py \
  julia --project=docs .amp/in/orb-backup/symmetry-sparsity/benchmark.jl B 20 1
```

The solver tolerances and 300 s solver limit are unchanged. ALLOW_FEASIBLE does
not change the solver: it permits post-solve checks when both points are
reported feasible, while preserving the non-optimal termination status. It
does not relax the existing gap, normalization, symmetry or PSD checks; an
extra conic-dual feasibility check is also applied.

| Quantity | B follow-up | E earlier run |
| --- | ---: | ---: |
| Numerical lower bound | -8.909590866404827 | -9.999999988364584 |
| Bound per site | -0.4454795433202413 | -0.4999999994182292 |
| Independent moments | 2,205 | 9 |
| Hermitian blocks | 64 | 10 |
| Largest Hermitian side | 16 | 14 |
| First measured-stage sum | 386.294 s | 69.823 s |

B's whole process took 426.392 s (7.11 minutes), within the final 6-8 minute
estimate, with 6.441 GiB peak sampled RSS. Moment construction took 23.609 s,
decomposition 18.652 s, symbolic assembly 60.655 s, solve 251.324 s, and
verification 29.511 s. No warm B timing was taken. The earlier 398.471 s attempt
stopped on SLOW_PROGRESS before extracting residuals; it is not the final run.

Mosek again reported SLOW_PROGRESS, with both points FEASIBLE_POINT. Recorded
status is `residual_checks_passed_not_optimal`, not `verified` or OPTIMAL.
The moment energy was -8.909590866412945; the primal-dual gap was 8.118e-12.
The original matrix's minimum eigenvalue was -1.763e-10; signed symmetry
residual was zero; maximum conic-dual constraint violation was 5.016e-7.
These satisfy the harness's numerical residual limits, but are not a rigorous
lower-bound certificate. A tiny objective gap alone is not proof of accuracy.

B improves the numerical lower bound over E by 1.09040912196 for the chain,
or 0.05452045610 per site. This compares bound strength, not errors against an
independently computed exact ground energy. E's current sparse constraints are
substantially weaker here; its speed is not evidence of equivalent accuracy.
Symmetry reduction itself does not cause this loss: B and E impose different
matrix constraints. A and C remain unmeasured past their original caps.

## Wider windows: E4, E5 and E6

Vary WINDOW_WIDTH only; retain the degree-4 seed, its 7,801-row closure and
all 67 detected generators. All three use numerical coefficient projection,
the 2-million-moment cap and the same solver settings as the B follow-up.
The harness checks seed-block nesting and coverage of all closed rows;
block-orbit traversal checks invariance of each complete block collection.
No old constraint is lost: each smaller seed block is contained in a larger
one, and applying any group action preserves that containment.

Run sequentially, from the repository root:

```sh
export JULIA_LOAD_PATH="@:$PWD/docs:$PWD:@stdlib"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=4
export MOMENT_LIMIT=2000000 NUMERIC_PROJECTION=true RUN_SECONDS=1800 ALLOW_FEASIBLE=true
for width in 4 5 6; do
  WINDOW_WIDTH=$width REPORT_SUFFIX=-window$width \
    python3 .amp/in/orb-backup/symmetry-sparsity/limited_run.py \
    julia --project=docs .amp/in/orb-backup/symmetry-sparsity/benchmark.jl E 20 4
done
```

Each process has one first execution and three warm executions, with no input
caching between samples. The license check precedes the samples. Measured
stage sums include verification; process times also include loading, license
checks, report writes and between-sample GC. B's prior first run is reused;
there is no warm B sample. E4 is rerun with the same code as E5/E6, rather than
using the earlier symbolic-assembly timing.

| Route | Numerical bound | First stage sum (s) | Warm median (s) | Warm range (s) | Peak sampled RSS (GiB) |
| --- | ---: | ---: | ---: | --- | ---: |
| E4 | -9.999999988364882 | 57.424 | 24.532 | 23.709–26.344 | 1.349 |
| E5 | -9.342585454939290 | 127.321 | 75.430 | 74.620–76.962 | 1.503 |
| E6 | -9.342585453151760 | 281.368 | 226.320 | 211.453–239.272 | 1.778 |
| B, prior run | -8.909590866404827 | 386.294 | Not measured | Not measured | 6.441 |

All twelve E samples reported OPTIMAL and passed every residual check. Their
maximum objective gap was 1.289e-11, minimum original-block eigenvalue was
2.216e-11, maximum conic-dual violation was 2.205e-8, and signed symmetry
residuals were zero. Bounds were identical across all four samples at each
width. B retains the SLOW_PROGRESS qualification above; none of these are
rigorous interval certificates.

All E widths have 160 parent blocks, one block orbit, 42 restricted stabilizer
generators, and ten retained Hermitian blocks. Seed-block sides are 175, 295
and 415. The distinct moment-word counts are 12,481, 49,921 and 180,421; the
independent moment counts are 9, 23 and 74. Full block-size lists are in TOML.

E5 improves E4 by 0.6574145334 for the chain, recovering about 60.29% of the
observed E4-to-B bound difference. E6 improves E5 by only 1.788e-9, below the
solver tolerance, while costing more in all warm samples. E5 remains weaker
than B by 0.4329945885. E5 is the better measured sparse choice at these widths;
do not infer that every further width must plateau or that either bound equals
the exact ground energy.

The main added cost is moment handling, not solving. Warm per-stage medians
for E5/E6 were 37.415/122.199 s for moment coordinates, 27.354/95.233 s for
verification, and 0.170/0.610 s for solving. Per-stage medians need not sum to
the median total. E4/E5/E6 whole-process times were 178.672, 409.976 and
1,014.163 s, totaling 26.714 minutes. This exceeded the initial 10–20 minute
estimate and matched the revised 26–28 minute estimate. No resource cap was hit.

Raw results are `benchmark-E-20-window{4,5,6}.toml` and matching logs. Derived
tables are `window-comparison.csv` and `window-comparison.json`.

## Why E5 and E6 plateau: a six-site witness

`plateau_check.jl` constructs the periodic spin-1/2 Heisenberg matrix directly
in the computational basis, independently of the SDP code. Diagonal bonds
contribute +/-1/4; antiparallel neighbors flip with amplitude 1/2. Its six-site
64-by-64 matrix has ground energy -2.802775637731988, eigenvector residual
6.627e-15. Scaling its energy per bond to 20 bonds gives -9.342585459106628,
matching E5/E6 within 6e-9. The four-site ring similarly gives -10 when scaled.
The four/six/eight-site diagonalizations took 1.64 s (estimate: under 20 s).

A numerical match alone would not explain the plateau. `plateau_witness.jl`
constructs an explicit feasible moment assignment from the six-site ground
state, without solving an SDP or reading its optimal moments:

1. Evaluate the Pauli products of the 415 seed rows in sites 1 through 6.
2. Use the exact signed moment-orbit coordinates of the 20-site E6 relaxation
   to assign all 74 independent variables. Conflicting assignments differ by
   at most 1.684e-15; sign-forced-zero moments satisfy the same tolerance.
3. Check all 67 generator equations on every retained moment word: residual 0.
4. Check the original representative 415-row matrix: minimum eigenvalue
   -9.791e-15. Verified block transporters cover all 160 blocks in one orbit.
5. Check identity normalization and the actual 20-site objective:
   -9.342585459106655.

The witness check passed in 299.900 s, within the 3-5 minute estimate, with
1.235 GiB peak sampled RSS. This is floating-point evidence, not an exact
certificate. It is a moment assignment that passes E6, not a constructed
20-site quantum state. Because E6 implies E5's retained constraints, the
assignment also explains why E5 can attain the same energy. The new E6
constraints do not eliminate this low-energy assignment.

```sh
julia --startup-file=no .amp/in/orb-backup/symmetry-sparsity/plateau_check.jl
JULIA_LOAD_PATH="@:$PWD/docs:$PWD:@stdlib" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=4 \
  python3 .amp/in/orb-backup/symmetry-sparsity/limited_run.py \
  julia --project=docs .amp/in/orb-backup/symmetry-sparsity/plateau_witness.jl
```

Proposed next step, not yet run: construct E7 (535 seed rows per window) with
the same degree-4 seed and all generators. First ask whether this six-site
assignment extends to its new moments and PSD blocks; leave unassigned moments
as variables, never zero-fill them. A feasible extension shows that this
strengthening still admits the plateau. Infeasibility eliminates this witness,
but not necessarily every other optimum at the same energy, so then minimize
the full E7 objective. Take warm timings only after a useful bound improvement.
No claim is made that E7 must improve or that an odd/even rule holds generally.

## Remaining work

The prototype uses straightforward symbolic matrix assembly. Capped routes need
separate follow-up; a partial construction time is not a time to solution.
Production integration belongs in the agreed redesign architecture, not a PR
into main. The owner authorized a branch snapshot for discussion after these
experiments; no merge or production integration is included.

## B structure diagnostic: analyze before expanding

This supersedes the earlier proposed E7 experiment. No E7 run was made.

B imposes one positive-semidefinite matrix on all 2,401 seed rows. Its four
lattice/axis generators preserve that basis. Exact symmetry reduction retains
joint positivity for global linear combinations of these rows; the resulting
64 blocks (largest Hermitian side 16) are not local principal submatrices.
E5 closes rows under all 67 detected generators, then imposes a collection of
160 principal blocks of side 295. Its exact block reduction preserves those
selected constraints, not positivity of the full 7,801-row matrix.

`B_structure.jl` found that 565,201 of B's 2,883,601 upper-triangle row pairs
occur together in E5 blocks, restricted to the original B rows. These are raw
row-pair counts before moment identifications, not independent constraint
counts. No E5 block contains both X1 and X11. This does not identify the single
constraint or smallest block that explains the bound gap.

The script also tests completion of the translation-invariant 20-by-20
X-correlation matrix using the six-site witness's short-range correlations.
Unfixed correlations stay free. With distances 1 through 4 fixed, a positive
completion exists: the best minimum eigenvalue is +0.10361307477029293.
Thus this scalar short-range test does not explain E5's gap. It does not test
all higher-order moments or their joint compatibility with the full symmetry.
With distances 1 through 5 fixed, the best minimum eigenvalue is
-0.029678542412147735. The numerical separating bound, including residual
free-coefficient errors, is -0.029678542410158992. This rejects that extension
of the E6 witness, not every possible E6 optimum.

The orbit of the candidate block {I, X1, ..., X20} has 12 blocks of side 21.
It has not been added to a relaxation or solved. A useful next investigation
is to extract the global combinations in B's positivity certificate and test
which ones E5 lacks, before choosing a strengthening.

```sh
export JULIA_LOAD_PATH="@:$PWD/docs:$PWD:@stdlib"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=4
python3 perf/symmetry-sparsity/limited_run.py \
  julia --project=docs perf/symmetry-sparsity/B_structure.jl
```

Final result: STRUCTURE_PASS, 85.590 seconds and 1.083 GiB peak sampled RSS;
the estimate was about one minute. `B-structure.toml` and `B-structure.log`
contain the final evidence. The initial and rounding logs retain failed
diagnostic attempts. Results are floating-point evidence, not exact proofs.
