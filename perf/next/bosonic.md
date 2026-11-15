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
