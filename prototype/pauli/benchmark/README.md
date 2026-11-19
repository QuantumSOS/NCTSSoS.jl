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
- Archive `environment.json`, both lock pairs, `results.jsonl`, worker TOMLs
  and logs, native phase times, CVM identity and guard evidence. Verify archive
  integrity outside the VM before cleanup or release.

`run.py --smoke --out NEW_DIRECTORY` is solver-free and only checks ED/topology
and worker invocation. `../test_wang.jl` runs small correctness controls, not
performance measurements.
