# Next scale probe (builds, not an optimization campaign)

[Measurements, limits and profile findings](results.md) are committed with the
raw samples in `results.json` and eight build profiles in `profiles.txt`.
[Precompile scope and first-call measurements](precompile.md) cover the separate
common-path workload decision, with raw records in `precompile.json`.
[TermSplit construction cost](termsplit.md) records the finite-group optimization
and exact before/after program comparisons, with raw records in `termsplit.json`.

Run from the repository root with Python 3.11+ and Julia. Use a disposable Julia
environment containing this checkout, JuMP, MosekTools and SCS. Do not change the
root dependencies. The benchmark checkout is read-only; no runner is imported.

For example, create `/tmp/next-scale-env` with
`julia --project=/tmp/next-scale-env -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.add(["JuMP", "MosekTools", "SCS"])'`.
Use that path for `PATH_TO_ENV` below. The measured package versions are in
`results.json`; a newly resolved environment can select different versions.
Mosek solve checks require a working license. Builds need no license.

```sh
python3 perf/next/run.py chain next 16 4 contiguous --project=PATH_TO_ENV
python3 perf/next/run.py chain legacy 16 4 contiguous --project=PATH_TO_ENV
python3 perf/next/run.py torus next 6 3 ts --project=PATH_TO_ENV
python3 perf/next/run.py hydrogen next 4 1 pqg --project=PATH_TO_ENV --checkout=PATH_TO_BENCHMARK
python3 perf/next/run.py hydrogen next 8 1 t12_number --project=PATH_TO_ENV --checkout=PATH_TO_BENCHMARK
python3 perf/next/run.py bell next 3 2 dense --mode=solve --project=PATH_TO_ENV
```

Each invocation prints a table and creates a worker log, checkpoint TOML and
summary JSON under `--output` (default `.amp/in/orb-backup/scale`). Existing
checkpoints require `--replace`. `--mode=profile` replaces lowering with a fresh
build under Julia's sampling profiler, writing flat (inclusive-count sorted) and
tree profiles. This also permits profiling a build whose lowering exceeds the
budget. Run one process at a time.

For `next`, `--save-program` serializes the first program beside the checkpoint
as `*.toml.program.jls`, outside the build timer. This permits complete
standard-form comparisons without repeating a slow build. The file requires a
compatible Julia/package environment to deserialize; only load trusted files.

The watchdog checks RSS every 100 ms and kills a worker at 600 wall seconds or
20,000,000,000 bytes RSS. These limits include loading, compilation, input
construction, repeated builds and lowering. A limit during warmup is **not** a
warm-build measurement. Completed stages survive a later timeout. RSS is the
whole worker's high-water mark, not retained program size or allocated bytes.
The sampler can overshoot the memory limit between observations.

## Measurement contract

- Julia, BLAS, OpenMP and Mosek use one thread. No concurrent numerical runs.
- Model input construction is measured separately, outside the build timer.
- Warm the exact build once, then measure three fresh builds, collecting garbage
  outside each timer. Report median and range; preserve all samples in JSON.
- Warm lowering separately, then take three lowering samples. Next lowers to
  `moment_standard_form`; legacy lowers to its native `build_jump_model`, without
  an optimizer. **These representations differ**, so compare their time as
  native lowering costs, not identical matrix construction.
- Next lowering sizes are `size(A)` and `nnz(A)`. Legacy sizes count scalar JuMP
  constraint dimensions, variables and affine coefficients (including explicit
  normalization). Next substitutes the identity and omits that variable.
- `rows` means the sum of selected row-family sizes (Next, before cone reduction)
  or the sum of emitted PSD block sides (legacy). Those counts are not generally
  equal after symmetry reduction. Blocks and largest side describe the emitted
  symbolic cones; legacy real embedding can double a complex cone at lowering.
- Profiling is a separate additional build after measurement. No code is tuned.

## Roster and mathematical differences

**Chain:** N=16,32,40,64,100, periodic XX+YY+ZZ, contiguous order-4 rows.
The Next front door uses translation, reflection, two global axis maps and
moment-only complex conjugation, matching `model_chain_case` in
`test/next/pipeline/benchmark_family.jl`. The comparator is the legacy specialized
translation-invariant route with axis permutation symmetry.

**Recorded roster:** 3×2 (`n=6`) and 3×3 (`n=9`) tori, order 3, `dense` and `ts`.
Next includes both translations, coordinate reflections, coordinate exchange for
square tori, global axis maps and moment-only conjugation. Legacy uses dense or
term-sparse construction without those symmetry reductions. Legacy `ts` uses
the roster adapter's MMD; Next uses the default `TermSplit()` with
`MaximalElimination()`. Both measure one sparsity step, not the adapter's solve-driven refinement.
These are not asserted to be identical SDPs. Bonds are unique
undirected nearest-neighbor bonds, each with XX+YY+ZZ and no factor 1/4.

**Scale point, not a QMBCertify reproduction:** 6×6 (`n=36`) torus with the same
symmetries, increasing order. Options `regions` and `regions_ts` use
`PhysicsInformed(level=1)` site groups (graph-radius-one regions plus bond
supports), without and with `TermSplit()`. Order eventually saturates these
finite Pauli regions; it does not become the complete global hierarchy.
Also probe global `dense` and `ts` to expose the cost of the full row basis. No QMBCertify example
or relaxation is specified in either checkout.

**Hydrogen:** H4/H6/H8 STO3G at the committed R=1 Å geometry. Convert only numeric
JSON integral data to TOML using Python's standard library. Use spin-blocked
orbitals (all up, then all down), the benchmark's integral convention and
1e-14 coefficient cutoff. Next uses order-1 rows, P/Q/G (`pqg`) or
P/Q/G/T1/T2prime (`t12`), fixed total number and fixed N_up−N_down=0. Legacy uses
the structured v2RDM builder with spin blocking, spin-resolved traces, explicit
one-D blocks, and no S² condition. Next's `Sector` equations and unsplit condition
blocks are not identical to legacy's specialized spin-blocked representation.
The additional `pqg_number` and `t12_number` variants add separate up/down
`NumberRotation` generators to the same sectors and rows. Their measurements
are appended to `results.md`; they do not replace the original baseline.
The legacy build omits the constant nuclear offset; cone and moment counts do
not depend on it. No legacy hydrogen solve is run here.
The sole hydrogen solve is Next H4/PQG, native SCS, eps_abs=eps_rel=1e-5.
Its total energy includes nuclear repulsion. Compare to the committed FCI value
as a lower bound, not as an assertion that the PQG relaxation is exact.

**Bell:** inn22 n=5,6/order3 dense; n=5/order5 `cs_ts`. Both implementations use
the benchmark file's coefficients (its special I3322 matrix for n=3). Next
projectors have two commuting sites, one per party, with n noncommuting letters
per site. Its correlative graph therefore has two vertices and cannot separate
measurement settings inside a party. Legacy's correlative graph is on generators.
This distinction can change the sparse relaxation; report bounds and sizes
rather than assuming equivalence. Dense n=3/order2 and the largest feasible
roster row are solved with Mosek on both sides; report −minimum as the Bell bound.
No literature value is claimed.
