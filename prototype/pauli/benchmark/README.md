# Pauli Heisenberg benchmark — Tencent only

Use the reviewed integrated commit, its transfer locks, and the dedicated
Tencent 128 GB CVM. Do not run measurements on an implementation orb. This
protocol tests the ground-energy lower bound per site of the periodic XXX model
with J1=+1, J2=0 and one `σᵢ·σⱼ/4` per undirected bond. Both solvers receive
the same normalization. See [WANG_MATCH.md](../WANG_MATCH.md) for mathematical
domains, provenance and the differences from the order-4 upstream examples.

## Comparison matrix

| Cases | Order | Wang `extra` | Graph range | Strengthening |
| --- | ---: | ---: | ---: | --- |
| `chain-4`, `chain-6` | 2 | 1, 2 | 2, 3 | Off |
| `chain-8`, `chain-8-rdm` | 2 | 3 | 4 | Off / LSO `lol=4`, PSO 2, RDM 8 |
| `square-4`, `square-4-rdm` | 2 | 0 | 4 | Off / LSO `lol=4`, PSO off, RDM 8 |
| `chain-60`, `chain-80` (and `-rdm` variants) | 2 | 3 | 4 | Same switch as chain-8 |
| `square-6` (and `-rdm` variant) | 2 | 0 | 4 | Same switch as square-4; **exploratory** |

Chain pin: `0b152dfcd9bcde3403aa4d967642ce5261c80712`. Square pin:
`9a3ae67ff2bb6da1c04b72318ecbfcec71db5495` (contains `resort`, which the
chain pin lacks). Call untouched `GSB` with `H_supp=[[1,4]]`, coefficient `3/4`
for a chain or `3/2` for a square, and `SU2_symmetry=false`. Use `QUIET=false`
for diagnostics and seed the stock LSO filter with `20261118`. `L` means the
side length for square cases. See `reference.jl` for complete arguments.

The candidate gets `H` and an order/row/strengthening policy, not geometry
or a generator list. Before timing a case the harness audits word-set equality
against the stock basis expanded over axes and checks retained translations,
finite point and spin actions. A *matched-support* record is not a proof that
strengthened SDPs are identical. Local-products (depth 2) is exploratory; it
does not carry parity claims even if its numerical bound happens to match.
Full rows are an optional control (`--supports full graph_distance local_products`).
Use them on controls, not on the large staged cases. Square-6 radius 4 has
108 additional pair supports; do not report it as matched.

## Replay

Use Julia **1.12.7**, Mosek **11.2.8**, Python 3.11+, one thread each for
Julia/GC, BLAS and Mosek, and Mosek feasibility/gap tolerances `1e-8`.
Set `MOSEKLM_LICENSE_FILE` or use `~/mosek/mosek.lic`; never archive the license.
Set `TRANSFER_LOCKS` to the `locks/` directory from the verified transfer.
Record setup/precompile separately and do not prewarm measured solves.

```sh
bench=$PWD/prototype/pauli/benchmark
julia --startup-file=no --threads=1,0 "$bench/setup.jl" "$HOME/pauli-wang-env" "$TRANSFER_LOCKS"
sudo systemd-run --scope -p MemoryMax=100G -p MemorySwapMax=0 -- \
  sudo -u "$USER" --preserve-env=PATH,MOSEKLM_LICENSE_FILE \
  python3 "$bench/run.py" --env "$HOME/pauli-wang-env" --out "$HOME/pauli-wang-controls" \
  --reviewed-cvm-commit "$(git rev-parse HEAD)" --repeats 3 --timeout 1800 \
  --cases chain-4 chain-6 chain-8 chain-8-rdm square-4 square-4-rdm
```

The runner checks the checkout and locks, CVM RAM, cgroup v2 cap and zero swap,
and a licensed trivial Mosek SDP in each pinned environment. Adjust the exact
privilege wrapper with the CVM operator if necessary; keep the verified cap and
the same user/environment for both solvers. After reviewing the controls, use
separate guarded runs for `chain-60 chain-80`, `square-6`, and the corresponding
`-rdm` cases where budget permits. Keep timeout, OOM, unsupported cases, and
solver failures in the record. Do not retry one side with looser limits.

Each requested repeat starts a fresh worker. With the default three methods,
the starting order is reference, graph-distance, local-products; each repeat
rotates that order by one position. Numerical rejections (including
`SLOW_PROGRESS`) do not cancel later repeats. Timeout, OOM, worker failure and
structural mismatch still stop that method for the case. A run with any failed
attempt exits nonzero after recording the remaining eligible attempts.

## Evidence

- Report independent ED for the small controls (zero-magnetization sector has
  a representative of every spin multiplet on these even-site problems) and
  the finite-ring Bethe numerical oracle for large chains. No exact 6×6 oracle
  is included; report its error as unavailable, not zero.
- Report actual candidate raw rows, affine term count, conic dimensions, and
  emitted LSO, PSO and RDM counts/sides. Stock returns its basis and moment
  pool sizes, but not the emitted LSO rank. Preserve its phase logs and native
  solver residuals, and leave unavailable stock metrics unset.
- Compare fresh-process wall time and peak RSS only on accepted solves; give
  sample counts, medians and spread. A smaller graph is not faster if bound
  accuracy worsens by more than `1e-6` per site. The harness requires `OPTIMAL`,
  MOI feasibility (candidate) or stock printed `PFEAS`,`DFEAS` ≤ `1e-8`,
  relative gap ≤ `1e-7`, and acceptable signed energy error. This is a numerical
  screen, not a certificate. Never report a performance ratio across lanes.
- `reported_objective_per_site` and final primal/dual objectives, relative gap,
  residuals and statuses are diagnostics, even on a rejected attempt. Only
  accepted `results.jsonl` rows have `energy_per_site` and accuracy metrics.
  Stock primal and dual objectives and residuals come from the last completed
  Mosek iteration; stock does not expose separate dual status. Candidate
  `mosek_*_max_violation` fields describe the native returned solution, not
  the stock iteration residuals. Unavailable fields are omitted, never zero-filled.
- Archive `environment.json`, both lock pairs, `results.jsonl`, worker TOMLs
  and logs, native phase times, CVM identity and guard evidence. Verify archive
  integrity outside the VM before cleanup or release.

`run.py --smoke --out NEW_DIRECTORY` is solver-free and only checks ED/topology
and worker invocation. `../test_wang.jl` runs small correctness controls, not
performance measurements.

## Paper-size extension (not a paper-settings reproduction)

[arXiv:2604.01555v2](https://arxiv.org/html/2604.01555v2), Tables 3 and 10,
report −0.4432378 per spin for the periodic 100-site chain and −0.674580 for
the periodic 16×16 square. These are the **new SDP** entries, not exact
energies. `energy_minus_published_sdp_per_site` compares accepted results to
the appropriate printed table entry; it is not an acceptance threshold or
an error against the ground energy. The independent chain oracle is Bethe;
there is no exact oracle for these large squares.

From the verified strengthened chain-80 and square-6 runs, start separate
guarded runs with `--cases chain-100-rdm --supports graph_distance` and
`--cases square-8-rdm --supports graph_distance`. Use chain-90 as an
intermediate rung if chain-100 fails; its oracle is Bethe and it has no paper
table entry. Advance the square only after an accepted candidate, in order:
8 → 10 → 12 → 14 → 16. Stop at the first failed candidate rung. Keep stock-only
runs, if requested, separate; they do not advance the graph campaign.
An unfrustrated chain uses `extra=3` (candidate radius 4), LSO `lol=4`, PSO 2
and RDM 8. The square uses stock `extra=0` versus candidate radius 4, LSO
`lol=4`, PSO off and RDM 8. Both stock calls are the unchanged pinned order-2
`GSB` calls in `reference.jl`. A square word-pool mismatch is exploratory even
if the symmetries match. Retain failures and do not substitute sizes in an
existing run. Report support and settings lanes separately and do not
compare speed ratios across unmatched lanes.

The original detector's audit (before the pair-generator change) found translation order 100 and point-group
order 48 for chain-100. For square-16 it reaches the 2,000,000-trial graph
search limit without discovering regular translations. The runner therefore
records `required_symmetry_missing` and skips that candidate method, while
leaving stock runnable. The registry supports all 256 sites (four 64-site
limbs); the structural obstacle is not word width, the million-row guard, or
RDM size (a supported eight-site cluster when a layout is found). The square
has 46,849 selected candidate words, versus 28,417 stock words. Without
reduction the candidate triangle would contain 1,097,437,825 entries. The
same audit rejects the 14×14 and 12×12 candidate routes; 10×10 and 8×8 retain
both expected translations and all 192 point-group actions. The largest
structurally verified square candidate at that revision was 10×10, with a different row
pool from stock (18,301 vs. 11,101 rows). Audit success does not show that the
solve fits or returns an accepted bound.

The later pair-generator detector passed the larger-square audits. See the
[completed ladder](RESULTS.md) for the measured runs through 16×16, the
explicit-moment submission, and the separate original-coordinate diagnostics.

The paper reports a 32-core, 1 TB machine and Mosek 11.0 but does not list
per-row options. In the [chain pin](https://github.com/wangjie212/QMBCertify/blob/0b152dfcd9bcde3403aa4d967642ce5261c80712/examples/partition_function.jl#L3-L6),
the table's large chain energies are inputs to a partition-function example,
not calls that produced them. The
[historical `H2D.jl`](https://github.com/wangjie212/QMBCertify/blob/337f97d17854af9c38b4d320f4ba78311da71fe8/examples/H2D.jl#L10-L24)
does contain a 16×16 call, shown with its argument values below, but no result
ties it to Table 10:

```julia
GSB([UInt16[1;4]], [3/2], 16, 2, lattice="square", positivity=10,
    pso=2, extra=0, QUIET=false, correlation=true)
```

That historical revision defaults to `lso=true`, `lol=L÷2`, and `three_type=[1,1]`;
it is not the pinned comparison code. Its trailing 100-spin call selects
`d=4, positivity=10, pso=3, extra=0, correlation=true` but inherits the
square's `3/2` coefficient, not the chain's `3/4`. Neither historical call is
evidence of the precise published settings. The later pinned 10-site chain
example even negates its objective, and the pinned 4×4 positive square example
has `d=4, extra=0, rdm=0, pso=0, lso=0`; do not extrapolate either example.
There is no verified paper-settings lane for the staged sizes. Only add one
if a per-size author command or metadata can be tied to the published result.
