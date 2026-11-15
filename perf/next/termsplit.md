# TermSplit construction cost

The expensive step was finite-group reduction after splitting, not the support
graph. The 6×6 order-1 graph split took 0.000694 s and 139,896 bytes and produced
row sides `[1,36,36,36]`. The identity-only block then generated a stabilizer of
order 6,912 from seven generators. Its trivial representation did not need a
group-algebra calculation. On 3×3/order 3, 5,255 of 10,913 profile snapshots
included group multiplication, versus about 2% in support splitting.

Two changes retain the original programs:

- Group products follow generator edges and the original enumeration's parent
  links instead of composing complete Clifford maps. Element and generator
  order are unchanged. Storage is O(number of supplied generators × group order);
  a product takes one lookup per generator in the right factor's stored word.
- An identity-only row block uses its exact normalized row directly. It still
  follows the antiunitary real-form check. Other blocks use the existing
  projectors, orbit reduction and splitting algorithm.

## Warm builds

Julia 1.12.7, Linux x86_64, Intel Xeon 2.60 GHz, shared orb, no CPU pinning.
One Julia/BLAS/OpenMP thread; serial processes. Each process warms the build
once, collects three fresh builds with GC between them, then profiles a fresh
build. Times are median [minimum, maximum]; allocations and RSS use decimal GB.
The 600 s / 20 GB cap includes loading, compilation and all repeats.

Baseline library: [91098fb](https://github.com/QuantumSOS/NCTSSoS.jl/commit/91098fb).
Candidate base: [62d151c](https://github.com/QuantumSOS/NCTSSoS.jl/commit/62d151c),
plus this change. The pipeline is unchanged between those bases; #435 adds core
letter methods and precompile workloads. Compilation is outside warm samples.
[Raw samples, source hashes, versions and fingerprints](termsplit.json).

| Torus / order / option | Before warm s [min,max] | After warm s [min,max] | Allocated GB before → after | Peak RSS GB before → after |
|---|---:|---:|---:|---:|
| 6×6 / 1 / `ts` | —, cap during first build | 9.5364 [9.4853,9.6890] | — → 3.73971 | 2.214 → 1.933 |
| 6×6 / 1 / `regions_ts` | —, cap during first warm sample | 0.09197 [0.08176,0.10225] | — → 0.061887 | 2.489 → 2.118 |
| 3×3 / 3 / `ts` | 17.7401 [17.2845,18.7731] | 8.2143 [8.1961,8.2188] | 9.38311 → 1.50312 | 2.382 → 2.167 |

The 3×3 warm build is 2.16 times faster and allocates 6.24 times fewer bytes.
No warm speedup ratio is claimed for a capped baseline. The optimized processes
completed in 174.89, 46.10 and 86.87 s, respectively. RSS includes compilation
and profiling, not just the retained program.

For the uncapped 6×6 global correctness oracle, the original first build took
622.245 s and 410.616 GB (716.43 s for the process). The regional baseline first
build took 508.445 s and 359.466 GB. Candidate first builds took 13.896 s / 4.138 GB
and 4.033 s / 0.503 GB. These are first calls, not warm timings, and include
method compilation; the precompile workloads differ between the bases.

The multiplication proxy used 4,096 actual stabilizer products, each checked
against direct composition: median 28.2 ms / 25,048,280 bytes before, 0.569 ms /
16 bytes after. This is not the application speedup. Group construction itself
grew from 0.553 s / 297 MB to 1.430 s / 594 MB in that proxy. The complete-build
measurements above include that cost. In the final 3×3 profile, multiplication
occupied 19 of 4,292 snapshots; the previous bottleneck is removed.

## Complete programs are unchanged

All three target fingerprints match exactly, including sparse `A` indices and
coefficients, `b`, `c`, the objective offset, and every cone's kind, side, start,
length and label. No solver tolerance is used in this comparison. The 79 saved
programs from the fixed combination smoke slice also retain their fingerprints.

| Torus / order / option | Rows | Moments | Real variables | Cones | Largest side | SHA-256 before and after |
|---|---:|---:|---:|---:|---:|---|
| 6×6 / 1 / `ts` | 109 | 10 | 9 | 11 | 1 | `1a9745a1449149a9ff1c70d3ff313781947fc0f1706d3803414d461de71c4e55` |
| 6×6 / 1 / `regions_ts` | 973 | 4 | 3 | 6 | 2 | `850cdae1809851e1583641b7e45ffe9aae583689426692748779d87cdd1681b6` |
| 3×3 / 3 / `ts` | 2620 | 75 | 74 | 44 | 52 | `b163a1a311b83842c9f366aea0f0334ffeaa21a89681dd50ed3386e23e02ba56` |

The retained group-law tests compare every product and inverse for noncommuting
permutations and signed Clifford maps, including reordered subgroup elements.
They also cover a product near the 65,535-element limit. The identity-only
shortcut is exercised by the existing family/split cases and the target forms.

Validation: 41,060 Next assertions, 147 native Mosek checks, 101 SCS writer
checks, 29 SCS moment checks, Aqua 11/11, and ExplicitImports passed in 614.52 s.
The native Mosek license check solved a one-site problem at
−0.9999999999876992 (exact −1). The docs stamp check passed in 2.28 s.

```sh
julia --project -e 'include("test/next/runtests.jl"); include("test/next/conic/mosek.jl"); include("test/next/conic/extensions.jl"); include("test/quality/Aqua.jl"); include("test/quality/ExplicitImports.jl")'
julia --project docs/examples_stamp.jl
```

## Reproduction

Use the environment described in [README.md](README.md). Run serially:

```sh
python3 perf/next/run.py torus next 36 1 ts --mode=profile --save-program --output=OUTPUT
python3 perf/next/run.py torus next 36 1 regions_ts --mode=profile --save-program --output=OUTPUT
python3 perf/next/run.py torus next 9 3 ts --mode=profile --save-program --output=OUTPUT
```

For the original library, use a worktree at `91098fb` with the same environment
and these versions of `run.py` and `probe.jl`; program capture is outside the
timed builds. The capped global baseline needs a separate uncapped build for
its correctness oracle. Fingerprint each saved program in a compatible Julia
environment with:

```julia
using NCTSSoS, Serialization, SHA
s = NCTSSoS.Next.moment_standard_form(deserialize(only(ARGS)))
io = IOBuffer()
print(io, size(s.A))
write(io, s.A.colptr, s.A.rowval, s.A.nzval, s.b, s.c, s.offset)
print(io, [(c.kind,c.side,c.first,c.length,c.label) for c in s.cones])
println(bytes2hex(sha256(take!(io))))
```

The default `TermSplit()` remains `MaximalElimination`. The legacy benchmark
uses MMD and has a different program; its time is not an equivalent-program
baseline. No relaxation, solver, tolerance, dependency, or core contract changed.
The global 6×6/order-3 memory failure remains outside this work.
