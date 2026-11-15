# Packed bosons: correctness, width and cost

Run `julia --project perf/next/bosonic.jl cost`. The `narrow` and `fixed`
arguments run the 16-mode width comparison in separate processes. No solver
or license is needed for these measurements. The fixed variant replaces only
`word_type` in its process; it is not a package option.

Environment: Julia 1.12.7, Linux x86-64 orb, `icelake-server`, one Julia
thread and one BLAS thread. Base: `redesign` at `6551587`, plus this PR.
Each region is warmed, then measured three times with GC before each sample.
Compilation, input setup and GC are outside the timed region. These are small
sample distributions, not a claim about other machines or solver time.

## Word layout and scope

`BosonWord{B,K}` stores creation and annihilation exponents in two
`NTuple{K,UInt64}` packs. Each limb holds `64 ÷ B` lanes; unused high bits
are zero. No lane crosses a limb. `word_type` alone chooses B and K from
the registry size and generated degree. Modeling uses B=8 (capacity 255).
Products check exponent sums before packing and stream normal-ordered terms
to the caller's sink. Machine integer coefficient overflow is checked;
Int128 and BigInt can represent larger coefficients. Capacity failures are
`ArgumentError`s, not wrapped exponents or coefficients. BigInt arithmetic
is not allocation-free; the zero-byte claim covers reserved machine-number
sinks.

The existing two contracts are unchanged. Plain assembly, grouping, splitting
and site permutations work through their existing paths. NumberRotation is
implemented in the core. Pipeline Sector and continuous-symmetry admission
remain for the next PR. Stationarity, StateOptimality, DensityMatrix and
FermionicConditions remain rejected. GNS explains that finite matrices cannot
satisfy the bosonic commutation relation and points to `expectation(result,p)`.
Contact D/G and scalar moment equalities are not part of this PR.

## Width: retain B for allocation savings, not for a claimed time gain

Open 4×4 lattice, hopping t=1, interaction U=2, order 2, no sector, plain
`build` only. Both variants produce 561 rows, 29,529 moment keys, 58,904
real variables and one 561-row block.

| Variant | Word bytes | Times (ms), all samples | Allocated bytes, all samples |
|---|---:|---|---|
| B=3, K=1 | 16 | 65.208, 65.519, 60.852 | 30,039,208; 30,036,656; 30,036,656 |
| B=8, K=2 | 32 | 70.641, 64.131, 59.780 | 36,401,920; 36,397,648; 36,397,648 |

Median allocation decreases 17.5%, above the approximately 10% decision
threshold. The time ranges overlap: there is **no demonstrated runtime gain**.
This retains order-dependent variants and repacking in exchange for smaller
words and less allocated storage. The two processes took 35.5 s total.

## L=8 Bose–Hubbard chain, order 2, no sector

Periodic edges, t=1, U=2: `-Σ(bd[i]*b[j]+bd[j]*b[i])+Σ(n[i]^2-n[i])`.
This is the benchmark Hamiltonian without its N=8 sector equality, as agreed.
Next produces 153 rows, 2,445 moment keys, 4,844 real variables and one block;
the CSC standard form has shape 23,409 × 4,844.

| Stage | Next time, median [min,max] ms | Legacy time, median [min,max] ms | Next bytes (median) | Legacy bytes (median) |
|---|---:|---:|---:|---:|
| Model construction | 0.178 [0.170,0.184] | 0.345 [0.321,0.419] | 91,984 | 159,584 |
| Relaxation assembly | 5.081 [5.007,5.246] | 485.280 [473.432,628.201] | 2,529,920 | 276,940,632 |
| Solver-facing lowering | 1.235 [0.923,1.531] | 101.091 [96.793,131.930] | 2,813,440 | 112,976,296 |

Assembly includes basis/sparsity construction and moment assembly in both
routes. **Lowering outputs differ:** Next emits its CSC standard form; legacy
emits JuMP expressions and constraints. The table does not claim equal work
for that row. Neither route constructs an optimizer or solves in this run.
The entire cost process, including warmup/compilation, took 65.1 s.

## Product kernel

The same input in both implementations is `b1² b8 × bd8 bd1²`; its normal
form has six terms. The script checks term-by-term equality before timing.
Each sample evaluates 10,000 products, reads the result length, and sums it.
The packed sink is reserved outside timing; the legacy public product owns
its output.

| Route | Time per product, median [min,max] µs | Bytes per product |
|---|---:|---:|
| Packed, reserved Float64 sink | 0.325 [0.286,0.326] | 0 |
| Legacy normal-monomial product | 14.172 [13.497,15.564] | 7,520 |

## Correctness checks

- `test/next/core/boson.jl`: contract laws; random term-by-term products,
  adjoints and raw-word normal forms against legacy; repeated modes and limb
  boundaries; exponent/coefficient failures; Int128/BigInt exact factorial
  terms; reserved Int, Int128, Complex{Int}, Float64 and ComplexF64 sinks.
- `test/next/model/bosonic.jl`: displaced oscillator minimum −0.49 and stable
  coupled oscillator minimum 0 at orders 2 and 3; no-sector two-site BH;
  Correlative, TermSplit and SitePermutation; named unsupported workflows.
- Native Mosek acceptance: 21/21 checks. BH bounds at orders 2 and 3 were
  −2.2970154730 and −2.2970154414; legacy gave −2.2970152440 and
  −2.2970154127. Both are below the independent truncated-occupation
  diagonalization value −2.2360679775. The solve tolerance is 1e-6.
- COSMO acceptance: 21/21 checks with residual tolerances 1e-10 and
  decomposition disabled, following the existing dense-reference test pattern.
  At order 3, default test settings differed by 1.54e-6 between lowerings;
  the tighter settings reduced this to 6.5e-9. The acceptance tolerance
  remains 1e-6. These settings apply only to the bosonic comparison tests.
- GNS tests cover the substantial-negativity advice, concise display and
  coverage hooks (reserved coordinates must not become solved moments).

For the pipeline combination matrix, use the two-site BH fixture in
`BosonModelChecks.bh`, first without a sector, then add `Sector(sum(bd.*b),2)`
in the pipeline PR. Include site exchange and NumberRotation, both grouping
and split settings, plus rejection rows for unsupported strengthening/GNS.

## Pipeline sectors and number rotation

Bounds and costs measured on `redesign` at `ac82da4` plus the bosonic pipeline changes, with
the same Julia 1.12.7, icelake-server orb and one thread described above.
The historical core-only measurements above are unchanged.

`Sector(Q,q)` now accepts bosons. It adds state equations on the family's
moment words and a zero-variance equation; it does not impose the operator
identity `Q=q`. Products remain normal-ordered sums through the word contract.
`NumberRotation` uses exact diagonal labels, with letters generated from the
registry's algebra. Translation, reflection and moment-only site permutations
needed no boson-specific pipeline path. Stationarity, StateOptimality,
DensityMatrix and FermionicConditions remain rejected for bosons; GNS is
unchanged. The contact D/G benchmark remains out of scope.

### Bounds against legacy and an exact sector Hamiltonian

Run `julia --project perf/next/bosonic_sector.jl` in an environment with
Mosek and MosekTools. A trivial licensed solve runs first. All runs below use
native Mosek for Next and MosekTools for legacy, without dualization. Native
residual/gap tolerances are 1e-10; legacy uses 1e-8. Legacy at 1e-10 stopped
with `SLOW_PROGRESS` on the two-site order-three case and withheld its bound;
that failed attempt is not used as a numerical answer. The acceptance
tolerance remains 1e-6. The first successful retry took 158.07 s. The final
post-rebase run took 124.10 s; bounds and counts were unchanged.

All models have t=1, U=2, chemical potential zero and N=L. Exact diagonalization
uses **only** occupations summing to N, not a truncated infinite-space ground
state. The two-site sector matrix has eigenvalue `1-sqrt(5)`. For L=3 and L=4,
the sector dimensions are 10 and 35. Raw results, including every symmetry
variant and signed difference, are in [bosonic_sector.tsv](bosonic_sector.tsv).

| Sites / edges | Order | Exact sector minimum | Legacy bound | Next bound | Next − legacy |
|---|---:|---:|---:|---:|---:|
| 2 / open | 2 | −1.236067977500 | −1.236067973271 | −1.236067977513 | −4.24e-9 |
| 2 / open | 3 | −1.236067977500 | −1.236067943571 | −1.236067977512 | −3.39e-8 |
| 3 / open | 2 | −2.601332708930 | −2.608691897690 | −2.608691898790 | −1.10e-9 |
| 3 / periodic | 2 | −4.343100154147 | −4.346996021954 | −4.346996021945 | +9.53e-12 |
| 4 / open | 2 | −3.968701041177 | −3.997752575773 | −3.997752587198 | −1.14e-8 |
| 4 / periodic | 2 | −5.613864939118 | −5.628701286552 | −5.628701288616 | −2.06e-9 |

All 20 symmetry variants are below their exact sector minima within 1e-6.
No measured strengthening over legacy is larger than numerical error. Sector
closes more equations than legacy's one-sided localizer, so equality of these
bounds is not a claim that the two constructions are identical.

| Model | Plain real variables | Number rotation | Translation | Dihedral | Dihedral + number | Moment-only translation |
|---|---:|---:|---:|---:|---:|---:|
| 2-site open, order 2 | 165 | 27 | — | — | — | — |
| 2-site open, order 3 | 434 | 52 | — | — | — | — |
| 3-site ring, order 2 | 675 | 127 | 227 | 142 | 26 | 227 |
| 4-site ring, order 2 | 2,034 | 424 | 515 | 331 | 67 | 515 |

Every listed symmetry preserves the plain bound within 1e-6. The open 3- and
4-site number-rotation runs also preserve their plain bounds.

### Full two-site combination product

`CombinationMatrix.bosonic_cases()` enumerates 608 rows: all 38 existing
family choices, four finite choices (none, translation, dihedral, moment-only
translation), continuous choices none/number, and strengthening none/sector.
The fixture is `BosonModelChecks.bh` (chemical potential 0.5) with N=2.
304 order-two programs passed; 304 order-one cases reject missing quartic
objective coverage with a named error. There were no defects. The full run
took 208.93 s (estimate 2–5 minutes), before rebasing from `79ac93d` onto
the merged chain change. The 16 dense bosonic controls passed again in the
post-rebase full suite.

For each built program, the check compares equality subspaces and distinct
cone eigenvalues at two seeded moment vectors against the unreduced reference.
Maximum errors were 3.67e-15 and 1.60e-14 respectively. It also checks all
equations and cone positivity on the exact N=2 ground-state projector.
The Fock representation contains every occupation with total number at most
eight; a degree-six word acts exactly on the N=2 state. A degree guard refuses
to apply this oracle beyond its cutoff. The same physical state is used for
the cases without Sector. The optional three-site full product was not added;
larger finite groups are covered by the chain comparisons above.

Rows are appended in [combination_matrix.tsv.gz](combination_matrix.tsv.gz).
19,760 existing rejection rows were rerun to refresh the NumberRotation and
algebra-guard texts; their outcomes did not change. All 79 pre-existing
Pauli/fermionic/site-word standard-form fingerprints matched byte for byte.
The old 145-case acceptance slice retains its ordering; 16 dense bosonic
controls are appended to it.

Reproduce the bosonic product (one Julia and BLAS thread):

```sh
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project -e '
include("test/next/pipeline/combination_matrix.jl")
C=CombinationMatrix; f=C.fixture(:bosons,2)
rows=[C.run(c;f) for c in C.bosonic_cases()]
@assert all(r->r.outcome!="DEFECT",rows)
C.write_table(stdout,rows)' > bosonic_matrix.tsv
```

### Frozen L=8, N=8 periodic benchmark

Run `julia --project perf/next/bosonic.jl sector`. This is the dense order-two
`bose_hubbard_1d` instance, t=1, U=2, including its particle-number restriction.
Each stage is warmed once, then measured five times with full GC before each
sample. Compilation and GC are excluded from the measured region. Model
construction includes the sector declaration. Assembly includes basis,
sparsity, moment and sector construction. No optimizer or solve is timed.
The final cost process took 70.62 s. Together with the final bound comparisons,
it took 194.72 s, within their combined 3–5 minute estimate. All samples are in
[bosonic_sector_cost.txt](bosonic_sector_cost.txt).

| Stage | Next median [min,max], ms | Legacy median [min,max], ms | Next median allocated bytes | Legacy median allocated bytes |
|---|---:|---:|---:|---:|
| Model construction | 0.205 [0.191,0.283] | 0.360 [0.317,0.417] | 108,480 | 181,312 |
| Relaxation assembly | 378.784 [362.159,429.152] | 458.006 [425.009,471.624] | 402,752,200 | 299,159,336 |
| Solver-facing lowering | 3.503 [3.212,5.998] | 91.783 [88.023,93.379] | 10,599,832 | 116,721,800 |

**These are different programs and different lowering outputs.** Next emits
CSC standard form; legacy emits a JuMP model. Next adds 8,330 real sector
equations; legacy adds 297 polynomial one-sided equations (not the same unit
as real equations). Next has 153 moment rows, 19,173 moment keys, 38,180 real
variables and a standard-form matrix of size 40,069 × 38,180. Sector introduces
degree-six moments beyond the order-two PSD block. Legacy's shorter equation
domain is retained without modification. The Next assembly median is lower,
but the time ranges overlap and **Next allocates 34.6% more**. This does not establish compliance with the
no-worse memory budget, nor an equal-work speedup. Allocated bytes are not
peak RSS; this run does not measure peak memory.

### Combined verification after rebase

The Next suite, native Mosek/SCS checks, Aqua and ExplicitImports passed on
the combined tree: 45,152 printed passing assertions, including 147 Mosek
checks, 101 SCS writer checks and Aqua 11/11. ExplicitImports was clean.
Elapsed time was 994.93 s (16.6 minutes), within the 15–20 minute estimate.
The docs example stamp check also passed; no docs source or generated page
changed.

```sh
julia --project -e 'include("test/next/conic/mosek.jl"); include("test/next/runtests.jl"); include("test/next/conic/extensions.jl"); include("test/quality/Aqua.jl"); include("test/quality/ExplicitImports.jl")'
julia --project docs/examples_stamp.jl
```
