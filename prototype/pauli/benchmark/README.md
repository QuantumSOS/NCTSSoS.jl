# Pauli Heisenberg benchmark

Do not collect comparison data until the integrated prototype and this harness have been reviewed. Run measurements serially on the dedicated **Tencent Cloud CVM with 128 GB RAM**. This checkout is based on NCTSSoS `e59863669f0fb30497c84f219240c1cf1bf38973`; the runner records the final reviewed descendant commit, not just the base.

## Physics and cases

The target is the unrestricted antiferromagnetic ground-state energy per site of `H = Σ⟨i,j⟩ (XiXj + YiYj + ZiZj)/4`, with `J1=+1`, `J2=0`, no fields, and each undirected periodic nearest-neighbor bond counted once. There are `L` sites/bonds in the ring and `L²` sites/`2L²` bonds in the square torus. Both optimizers are given a per-site objective. The ED worker forms the **total** spin Hamiltonian independently in the zero-magnetization sector and divides its lowest eigenvalue by the site count. This sector contains a member of every spin multiplet for even site count. The 4×4 ED matrix has side 12,870 (about 1.24 GiB for matrix entries alone); the runner computes it once, outside measured solver runs.

| Case ID | Sites | Reference `extra` | Prototype support | Use |
|---|---:|---:|---|---|
| `chain-4`, `chain-6` | 4, 6 | 1, 2 | full (`:ExactNPA`) | ED controls; default |
| `square-4` | 16 | 0 | full (`:ExactNPA`) | ED control and upstream-failure record; default |
| `chain-40` | 40 | 3 | contiguous, range 4 | medium comparison; default |
| `chain-24`, `chain-64` | 24, 64 | 3 | contiguous, range 4 | opt-in scaling |
| `chain-8-rdm` | 8 | 3 | contiguous, range 4 | strengthened ED control; opt-in |
| `chain-40-rdm` | 40 | 3 | contiguous, range 4 | strengthened medium ring; opt-in |

Unstrengthened rows use order 2 without stationarity or extra RDMs; strengthened rows use order 3 on both sides. All cases omit PSD state optimality, sector restrictions, and post-solve certification. The prototype uses `moi_form=:auto`; contiguous cases explicitly select `:circulant_dft`, `:NoElimination`, `:Maximal`, `ts_step=0` and require the reported lattice route. A reference chain `extra=k` selects two-site distances `1:k+1`; order 3 adds the consecutive triple pattern `three_type=[1,1]`. For a square `extra` adds entries to QMBCertify's displacement table (nine base offsets on 4×4), not to the prototype's pair range. Even matching pair placements or ED values does not establish equivalent moment constraints or symmetry reductions.

**Strengthened pairs compare bound quality, not identical programs:** the prototype uses `commutator_linear_order=3` over selected row words and `rdm_positivity_k=4` on explicitly listed translated contiguous four-site regions. QMBCertify uses `rdm=8` on the eight consecutive sites (translation reduction), and `lso=true, lol=4`; its stationarity pool is `generate_mons(L, 4, 7)`, filtered against the retained support with fixed seed `20261118`. `rdm_psd_blocks` and `stationarity_equalities` must both be positive. On the prototype these count emitted RDM PSDs and surviving scalar stationarity rows. On the reference these count PSD constraints added by `posepsd8!` and filtered stationarity equations represented by free SOS multipliers, not an independent numerical rank. The reference and prototype region sizes and stationarity pools differ; the runner records numerical bound differences, not speedup ratios. Opt-in cases are screening targets, not verified memory-fit claims.

## Reference and settings

Pinned [QMBCertify 0b152df](https://github.com/wangjie212/QMBCertify/tree/0b152dfcd9bcde3403aa4d967642ce5261c80712). [`examples/ground_state.jl`](https://github.com/wangjie212/QMBCertify/blob/0b152dfcd9bcde3403aa4d967642ce5261c80712/examples/ground_state.jl) uses `[[1,4]]` (one representative of the same-axis nearest-neighbor correlator). Spin-rotation symmetry gives three components, and a square adds two directions per site: `coe=[3/4]` for `lattice="chain"` and `[3/2]` for `lattice="square"`. The first chain example passes `-coe`; do **not** copy that sign. `GSB` already returns per-site energy with these coefficients; do not divide again. These are harness settings, not claims to reproduce any paper table or historical run.

Unstrengthened defaults: `GSB(supp, coe, L, 2; H_supp=supp, H_coe=coe, lattice, extra, SU2_symmetry=false, lso=false, lol=L, pso=0, rdm=false, energy=[], three_type=[1,1], J2=0, correlation=false, Gram=false, QUIET=true, writetofile=false, mosek_setting=mosek_para(1e-8,1e-8,1e-8,1))`. For strengthened cases change positional order to 3 and use the three keyword overrides above. `L` is the side length of a square, not its site count. The in-memory hook brackets the same `optimize!` call as the earlier Tencent reference hooks, checks exact source anchors, and skips upstream objective/moment extraction after a rejected solve. It does not change either ring SDP or patch the pinned checkout.

**Stock square-4 is blocked:** pinned `eigen_circmat` calls the undefined `resort`, before model construction completes. The runner records `upstream_failure` and no paired difference for this case; any prototype row is standalone against ED. A campaign including this case reports incomplete (nonzero exit). No repaired reference is included. Square-6 is excluded; the historical order-4 `rdm=9`, `pso=3` configuration on another revision reached about 270 GiB.

## Run on the CVM

Use Linux with cgroup v2, Python 3.11+, Julia 1.12, and a readable local Mosek license (`MOSEKLM_LICENSE_FILE` or `~/mosek/mosek.lic`). Do not place license material in the result directory. On an idle machine, from the reviewed checkout:

```sh
bench=$PWD/prototype/pauli/benchmark
julia --startup-file=no --threads=1,0 "$bench/setup.jl" "$HOME/pauli-bench-env"
systemd-run --user --scope -p MemoryMax=100G -p MemorySwapMax=0 -- \
  python3 "$bench/run.py" --env "$HOME/pauli-bench-env" --out "$HOME/pauli-bench-results" \
    --reviewed-cvm-commit "$(git rev-parse HEAD)" --repeats 5 --timeout 1800
```

The commit flag is the operator's attestation that the revision has been reviewed and this is the specified CVM. The runner checks ancestry, checkout cleanliness, host RAM, cgroup memory/swap limits, an actual trivial Mosek SDP solve, and the unchanged dependency manifest. The scope needs permission to set memory limits (ask the CVM operator if user scopes cannot do so). Both methods share one Manifest, local license, Mosek tolerance triple, and single-thread Julia/GC, BLAS and Mosek configuration. No solver fallback. Keep the CVM idle, record the instance SKU and any affinity/NUMA policy with the results, and stage larger cases separately via `--cases`. A failed worker leaves its log and available partial timings; a timeout kills its process group. Do not compare a completed time with a timeout or retry a single side with looser limits.

Archive `environment.json`, `Project.toml`, `Manifest.toml`, `results.jsonl` and logs. Reuse those dependency lockfiles for repeat campaigns (local package paths in the manifest must point to the same pinned checkouts). The setup command resolves an initial manifest; rerunning setup later is not equivalent to using the archived manifest.

## Read the evidence

- `wall_seconds` starts before launching Julia and ends at reaping it (includes process startup, imports, JIT, audit hooks, output, and teardown; polling resolution about 20 ms). Each run is a fresh process with package caches available, not a warm re-solve. The order within a pair alternates. Report each method's sample count, median and spread before discussing any ratio.
- `peak_rss_kib` is Linux per-child `wait4` high-water RSS, not live polling or whole-host memory. The cgroup cap separately charges the whole campaign scope. `build_seconds` covers prototype Hamiltonian creation + `compile_relaxation`, or reference `GSB` assembly up to `optimize!`. `solve_seconds` covers the public prototype `solve_ground_state` call (including setup) or reference JuMP `optimize!` (including transfer). These boundaries are not equivalent; use `optimizer_reported_seconds` for Mosek-only time and `wall_seconds` for the method-level comparison. Prototype internal timing fields are retained separately.
- The acceptance screen requires `OPTIMAL`, feasible primal and dual statuses, finite energy agreeing with the solver objective, and primal–dual objective gap ≤ `1e-7` (scaled by `max(1, |p|, |d|)`). Numerical objectives are not rigorous certificates, and this screen does not independently verify feasibility residuals. Energy must lie between the sum of bond minima and the ED energy (or the staggered product-state upper bound if ED is absent), to `1e-6` per site. `signed_error_per_site = bound − ED`, `exact_gap_per_site = ED − bound`, and `relative_exact_gap = (ED − bound)/abs(ED)`; never clip the signs. The paired difference is **prototype − reference** (positive means a tighter numerical lower-bound estimate). No paired difference is emitted if either side fails; `strength_equivalence` is always `not_established`.

No-solver smoke: `python3 prototype/pauli/benchmark/run.py --smoke --out NEW_DIRECTORY`. This checks analytic ring ED values, torus edges and child measurement. To check pinned reference assembly without a license after environment setup, run the following. Replace the case arguments with `chain 8 3 4 true true` for strengthened assembly, or `square 4 0 4 true false` to observe the upstream failure:

```sh
julia --startup-file=no --threads=1,0 --project="$HOME/pauli-bench-env" \
  "$bench/worker.jl" reference-build-smoke /path/to/new-dry-reference.toml chain 4 1 2 true false
```

Neither smoke mode measures performance or exercises the integrated prototype solve.
