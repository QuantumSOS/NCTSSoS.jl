# Chain construction cost

The periodic Heisenberg chain with the contiguous degree-4 family now builds
faster than the current specialized implementation at all five measured sizes.
The complete Next standard forms are unchanged. Allocated bytes are lower than
the current package, but **peak RSS is not lower at every size**. This is a
build-time improvement, not a claim that every memory budget is met.

## What changed

The lattice builder already checks the exact normalizer relation
`g*T*inv(g) == T` or `inv(T)` for its point group. Moment closure now completes
each translation orbit first and applies these normalizers to one word per
translation orbit. The relation propagates to the other words by the checked
conjugation law. Arbitrary moment-only actions still visit every word; they do
not need to normalize translation. Star and moment-key closure also remain on
every word. A repeated translation vertex still checks its phase relation.

The queue carries each word's phase and conjugation flag, avoiding repeated
relation lookups. The lattice path computes the eight phase/conjugation
coordinate maps once per orbit. The original dictionary representation and
all cached aliases remain, including symmetry-forced zeros used by continuous
closure. Other routes retain their generator traversal. Neither the full
lattice group nor a full site-space moment matrix is built.

## Warm build measurements

Julia 1.12.7, Linux x86_64, Intel Xeon 2.60 GHz, shared orb, no CPU pinning.
One Julia/BLAS/OpenMP thread; serial processes. The unchanged runner warms the
build once, then measures three fresh builds with GC between them. Times are
median [minimum, maximum]. GB denotes decimal bytes. The 300 s / 20 GB limits
cover the whole process, not just the timed builds.

Before and legacy:
[91330d9](https://github.com/QuantumSOS/NCTSSoS.jl/commit/91330d9).
After:
[694a6fa](https://github.com/QuantumSOS/NCTSSoS.jl/commit/694a6faa34b3629f85252356e180e655bcc90d49),
based on the merged bosonic core at
[79ac93d](https://github.com/QuantumSOS/NCTSSoS.jl/commit/79ac93d).
The intervening core change adds bosonic methods; it does not change the Pauli
chain implementation or pipeline. Compilation is outside warm samples.
[Raw samples, source hashes, dimensions and fingerprints](chain.json).

| Sites | Next before s [min, max] | Current specialized s [min, max] | Next after s [min, max] |
|---:|---:|---:|---:|
| 16 | 0.1722 [0.1646, 0.1728] | 0.1856 [0.1811, 0.1886] | 0.1241 [0.1223, 0.1312] |
| 32 | 0.6935 [0.6865, 0.7389] | 0.6837 [0.6806, 0.7019] | 0.5061 [0.4839, 0.5783] |
| 40 | 1.1543 [1.0872, 1.1921] | 1.3088 [1.0291, 1.3174] | 0.7638 [0.7305, 0.8161] |
| 64 | 3.1592 [2.9932, 3.7568] | 2.8535 [2.7917, 2.9202] | 2.1065 [1.9997, 2.7723] |
| 100 | 11.0396 [10.9063, 11.6919] | 10.1023 [9.8865, 11.4132] | 7.1173 [7.0857, 7.7713] |

All three final warm samples are below all three corresponding legacy samples.
These are shared-orb measurements, not a guarantee for every machine. The final
five runner processes took 236.4 s in total. Including fingerprint checks and
cache/profile diagnostics, the final measurement batch took about 6.3 minutes
against an estimate of 6–8 minutes.

| Sites | Legacy allocated GB | Next allocated GB | Legacy peak RSS GB | Next peak RSS GB |
|---:|---:|---:|---:|---:|
| 16 | 0.201842 | 0.085383 | 1.786 | 2.024 |
| 32 | 0.641147 | 0.279818 | 2.367 | 2.139 |
| 40 | 0.925111 | 0.453985 | 2.325 | 2.219 |
| 64 | 2.041106 | 1.030839 | 2.565 | 2.658 |
| 100 | 4.826579 | 2.663757 | 3.417 | 4.067 |

Allocation bytes are nearly unchanged from Next before this optimization.
The change removes repeated action work, not the stored aliases. RSS includes
loading, compilation, build repeats and lowering. The overall memory target is
therefore not established by these results; RSS is higher at 16, 64 and 100
sites. Allocation bytes are not allocation counts, and no 10× allocation claim
is made.

## Why 64 to 100 sites remains a large step

Measured on the final code:

| Sites | Moment labels, including constant | Cached words | Word bytes | Dictionary capacity | Dictionary bytes |
|---:|---:|---:|---:|---:|---:|
| 64 | 6,195 | 4,132,609 | 16 | 8,388,608 | 243,269,744 |
| 100 | 10,065 | 10,518,001 | 32 | 16,777,216 | 754,974,832 |

There are 2.55 times as many cached words, not just 1.62 times as many moment
labels. Word storage doubles above 64 sites, dictionary capacity doubles, and
dictionary storage grows 3.10 times. These are measured storage changes, not
isolated causal timing estimates. The final profile still has 2,467 of 3,475
snapshots under the Fourier-kernel path. Keeping every translated alias remains
a substantial cost and is not fixed by the faster traversal.

## Correctness evidence

- All five complete chain standard forms have identical SHA-256 fingerprints
  before and after: matrix dimensions and CSC arrays, `b`, `c`, objective offset,
  and cone kind, side, position, length and label. This checks the programs, not
  only their dimensions or solver bounds. The Next and legacy cones need not be
  identical: their largest sides remain 16 and 30 respectively.
- 79 combination smoke-program fingerprints match the saved reference.
- 18,332 temporary closure-law assertions passed against the full generator
  traversal on Pauli, fermionic and bosonic words. Each seed used fresh tables.
  Cases include antiunitary phases, real forms and a CNOT that does not normalize
  translation. Returned coordinates, labels, all cached aliases and moment
  queries agree. Runtime: 16.73 s against a 30–60 s estimate.
- The combined Next suite, native Mosek/SCS checks, Aqua and ExplicitImports
  passed: 45,132 assertions, including 147 Mosek checks, 101 SCS writer checks
  and Aqua 11/11; ExplicitImports is clean. Runtime: 951.67 s against a 12–16
  minute estimate. The licensed native Mosek tests verified the license by
  solving, not just importing the package.

## Reproduction

Use the environment described in [README.md](README.md). For each of 16, 32,
40, 64 and 100 sites, run both routes, with separate output directories for the
baseline and candidate checkouts:

```sh
python perf/next/run.py chain next 100 4 contiguous --save-program --seconds 300 --gb 20 --output /tmp/chain-next
python perf/next/run.py chain legacy 100 4 contiguous --seconds 300 --gb 20 --output /tmp/chain-legacy
```

The saved programs are outside the timer. The raw JSON records the exact
commands, versions, source revision, worker hashes and all warm samples.
The combined validation command is:

```sh
julia --project -e 'include("test/next/conic/mosek.jl"); include("test/next/runtests.jl"); include("test/next/conic/extensions.jl"); include("test/quality/Aqua.jl"); include("test/quality/ExplicitImports.jl")'
```

The scratch laws, full-program comparison and cache/profile scripts are kept
in this work thread under `.amp/in/orb-backup/chain-cyclic-laws.jl`,
`chain-equivalence.jl` and `chain-counts.jl` respectively. No public API,
word contract or conic contract changed.
