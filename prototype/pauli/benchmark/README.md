# Pauli Heisenberg benchmark

Do not collect comparison data until the integrated prototype and this harness have been reviewed. Run measurements serially on the dedicated **Tencent Cloud CVM with 128 GB RAM**. This checkout is based on NCTSSoS `e59863669f0fb30497c84f219240c1cf1bf38973`; the runner records the final reviewed descendant commit, not just the base.

## Physics and cases

The target is the unrestricted antiferromagnetic ground-state energy per site of `H = Σ⟨i,j⟩ (XiXj + YiYj + ZiZj)/4`, with `J1=+1`, `J2=0`, no fields, and each undirected periodic nearest-neighbor bond counted once. There are `L` sites/bonds in the ring and `L²` sites/`2L²` bonds in the square torus. Both optimizers are given a per-site objective. The ED worker forms the **total** spin Hamiltonian independently in the zero-magnetization sector and divides its lowest eigenvalue by the site count. This sector contains a member of every spin multiplet for even site count. The 4×4 ED matrix has side 12,870 (about 1.24 GiB for matrix entries alone); the runner computes it once, outside measured solver runs.

The prototype compiler receives only `H`, the order and a single policy shared by all geometries: `:GeneralGraphFast` at order 2 or `:TightGroundState` at order 3. No site order, lattice shape, distance cutoff, spatial group, RDM regions or solver form is supplied. Record the support, achieved reductions, surviving scalar equalities, cones, time, memory and bound actually obtained. Requested constraints need not survive reduction. The previous ring-specific comparison is reproducible from `f8b1479`; its timings do not describe this workflow.

Default controls are `chain-4`, `chain-6`, `square-4` (ED) and `chain-40`. Other rings are opt-in with `--cases`. The reference `extra` is 1 and 2 at chain-4 and chain-6, 3 on larger rings, and 0 on square-4. On `chain-8-rdm` and `chain-40-rdm` the reference requests order 3, eight-site RDMs and local stationarity (`lol=4`, seed `20261118`). Stock does not return emitted strengthening counts: do not infer them from the request or claim an audited strengthened pair.

For unstrengthened large rings, compare fresh-process wall time only for accepted solves. Report both errors against the finite-ring numerical Bethe ground energy (not a rigorous certificate); do not claim a win if absolute error worsens by more than `ENERGY_ATOL = 1e-6` per site. Record peak RSS even when time improves. The Bethe calculation runs outside both timed solvers and is checked against ED on 4, 6, 8 and 10 sites.

## Reference and settings

Pinned [QMBCertify 0b152df](https://github.com/wangjie212/QMBCertify/tree/0b152dfcd9bcde3403aa4d967642ce5261c80712). [`examples/ground_state.jl`](https://github.com/wangjie212/QMBCertify/blob/0b152dfcd9bcde3403aa4d967642ce5261c80712/examples/ground_state.jl) uses `[[1,4]]` (one representative of the same-axis nearest-neighbor correlator). Spin-rotation symmetry gives three components, and a square adds two directions per site: `coe=[3/4]` for `lattice="chain"` and `[3/2]` for `lattice="square"`. The first chain example passes `-coe`; do **not** copy that sign. `GSB` already returns per-site energy with these coefficients; do not divide again. These are harness settings, not claims to reproduce any paper table or historical run.

Call the public `QMBCertify.GSB` directly from the clean checkout, with `QUIET=false` for native output. No source rewriting, optimizer hook, square shim or fallback. `reference.jl` records the API arguments: order 2 controls use `SU2_symmetry=false` (the upstream default), `lso=false`, `pso=0`, `rdm=false`; both solvers use one thread and the same `1e-8` tolerance settings. `L` is the side length of a square, not its site count.

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

Record setup wall time separately. Setup installs the local `PauliWorkflowRedesign` package so fresh workers can use its Julia cache. The cache workload builds a symbolic three-site open chain and anisotropic range-one rings at 5 and 65 sites, covering one- and two-limb Pauli words. It does not run a solver or retain a plan. For environments from the original direct-include comparison, add the local prototype to a **copy** of the locked environment, check the dependency diff, and archive the new manifest. Do not replace the original lockfile or run a measured case to warm the cache.

## Read the evidence

- `wall_seconds` starts before launching Julia and ends at reaping it (includes process startup, imports, JIT, audit hooks, output, and teardown; polling resolution about 20 ms). Each run is a fresh process with package caches available, not a warm re-solve. The order within a pair alternates. Report each method's sample count, median and spread before discussing any ratio.
- `peak_rss_kib` is Linux per-child `wait4` high-water RSS. Prototype internal build/solve times are recorded separately. Only compare end-to-end wall time across methods; inspect the unmodified upstream transcript for its phase timings. Unreported diagnostics remain unavailable, not zero.
- For the prototype, require `OPTIMAL` with both MOI statuses `FEASIBLE_POINT`. Stock `GSB` does not return its model: this pin prints a termination line for every nonoptimal result. Require a normal return, one matching `optimum` line, no nonoptimal status, and a complete native Mosek iteration table with final `PFEAS` and `DFEAS` ≤ `1e-8`. Record those printed residuals, not inferred MOI primal/dual statuses. For both methods require finite energy agreeing with the solver objective and relative primal–dual gap ≤ `1e-7` (using the printed final objectives for stock). This is a numerical screen, not a certificate. Energy must lie between the sum of bond minima and ED (or a staggered product-state upper bound) to `1e-6` per site. Report signed errors without clipping. A positive prototype-minus-reference difference means a tighter numerical estimate. Do not emit a pair if either side fails.

No-solver smoke: `python3 prototype/pauli/benchmark/run.py --smoke --out NEW_DIRECTORY`. This checks analytic ring ED values, torus edges and child measurement, not either solver. To call upstream alone, run `worker.jl reference RESULT.toml chain 4 1 false` in the benchmark environment with `--threads=1,0`.
