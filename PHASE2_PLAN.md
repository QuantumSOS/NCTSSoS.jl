# Phase 2 — singular $\mathcal{A}\mathcal{A}^{*}$ diagnostic on H₂ / Nk=2

> Status: **planning**, not yet started. This file specifies the experiment; it does not run it.

## 0. Goal in one sentence

For each residual symmetry block of the H₂ / Nk=2 PQG V2RDM SDP, decide *empirically* which inner linear-solve strategy the production-scale solver (H₄ / Nk=2 and beyond) should use, by measuring the spectrum, kernel, and sparsity of $\mathcal{A}\mathcal{A}^{*}$.

The output of Phase 2 is **a decision table**, not a solver. No SDP gets solved here.

## 1. Why H₂ / Nk=2 (and not H₄ / Nk=2)

- **Small enough for ground truth.** $M = 2\,\text{nk}\,\text{norb} = 16$ spin-orbital modes; the largest PSD block is $O(M^2) = O(256)$; the dual dimension $m$ is comfortably $O(10^3{-}10^4)$. Dense SVD and dense Cholesky are both seconds-to-minutes in BLAS. Every diagnostic in §3 admits an *exact* answer, not an approximation. That makes it a *calibration* run.
- **Big enough for the structure to appear.** Nk=2 means translation symmetry is nontrivial (two momentum sectors $k=0,\,\pi$), Sz blocking is nontrivial, hermiticity matters. The qualitative pattern of "post-Wedderburn $\mathcal{A}\mathcal{A}^{*}$ per block" is exhibited, just at small numbers.
- **Same algebra as production.** Anticommutation rules, PQG structure, hermiticity, particle conservation are all identical to H₄/Nk=2. Numerical magnitudes will scale; *which redundancies exist* is a structural question whose answer transfers.

If the diagnostic at H₂/Nk=2 already shows "no spectral gap" or "kernel is huge", that is a strong falsification signal for the larger production case — facial reduction (Paper 5) is on the critical path before any solver work.

## 2. Prerequisites (not yet in repo)

The repo currently has H₂/Nk=1 integrals. We need H₂/Nk=2:

- [ ] Generate `test/data/assets/h2_chain_nk2_active_2e4o_integrals.txt` (PySCF dump on HAI, mirroring the H₄/Nk=2 pipeline used to produce `h4_chain_nk2_integrals.txt`).
- [ ] Generate companion `test/data/assets/h2_chain_nk2_reference.toml` with PySCF HF / MP2 / CCSD reference energies.
- [ ] Add **one new** script `demos/h2_periodic_nk2_moment_sos.jl` per the implementation contract in §12. Do **not** clone or modify `demos/h2_periodic_moment_sos.jl` or any other file under `demos/`.

These three are setup, not science. Budget: < 1 day.

## 3. Numerical artifacts to extract

The physical solver-relevant operator is $\mathcal{A}\mathcal{A}^{*}$, where $\mathcal{A}$ is the constraint operator that takes the primal moment matrix (or vector of moments) to the constraint values. In the JuMP build of `_solve_complex_moment_problem`, $\mathcal{A}$ is implicit; for diagnostics we materialize it as a sparse matrix.

Two kinds of $\mathcal{A}$ are worth measuring separately:

- **$\mathcal{A}_{\text{full}}$** — combined (Zero equality + HPSD-coupling) constraint operator. This is what an ADMM/AL outer loop sees.
- **$\mathcal{A}_{\text{eq}}$** — equality constraints only. Closer to what facial reduction targets.

Run §3 on both, label clearly.

For each (residual symmetry block $B$, choice of $\mathcal{A}$):

### 3.1 Spectrum of $\mathcal{A}_{B}\mathcal{A}_{B}^{*}$

Compute by **dense SVD** (we are at sub-$10^4$ dimension; this is gold-standard, no Krylov tricks needed at this scale).

Report:

| Quantity | Symbol | Why |
|---|---|---|
| Largest singular value | $\sigma_{1}$ | Conditioning numerator. |
| Numerical rank at tolerances $\tau \in \{10^{-10}, 10^{-12}, 10^{-14}\}$ | $\hat r$ | Sensitivity of "rank" to threshold. |
| Smallest positive singular value | $\sigma_{\hat r}$ | Conditioning denominator on range. |
| Range condition number | $\kappa_{\text{range}} = \sigma_{1}/\sigma_{\hat r}$ | Determines Mazziotti-CG iteration count. |
| Spectral gap ratio | $g = \sigma_{\hat r}/\sigma_{\hat r + 1}$ | $g \gg 1$ ⇒ clean kernel ↔ Tikhonov / Direct viable. $g \approx 1$ ⇒ continuous tail ↔ Tikhonov damages information. |
| Singular value histogram | $\{\log_{10}\sigma_{i}\}$ plot | Visual sanity check. |

### 3.2 Kernel structure

- **Empirical kernel dim**: $m_{B} - \hat r$.
- **Predicted kernel dim**: count, ahead of solving, the algebraic redundancies enforced by the symmetry block. Sources to count separately:
  - Anticommutation simplification residues.
  - Hermiticity relations (Re/Im pairing).
  - Particle number / Sz / momentum sector consistency (constraints that are tautological inside a fixed-quantum-number block).
- **Excess kernel** = empirical − predicted. **This is the load-bearing number.** Excess > 0 means structural redundancies exist that the current symmetry catalog misses; those are facial-reduction targets that aren't yet being captured.
- Save a basis of the empirical kernel (one matrix per block) for inspection. Cluster basis vectors by support pattern (single-orbital vs cross-orbital) — gives hints about which algebraic identity each missed redundancy corresponds to.

### 3.3 Sparsity & fill-in

For the candidate "factor once, reuse forever" path:

- nnz, density of $\mathcal{A}_{B}\mathcal{A}_{B}^{*}$.
- Sparsity pattern image (small enough to plot for H₂/Nk=2).
- Fill-in under **AMD** reordering: nnz of $L$ in $LL^{\top} = \mathcal{A}_{B}\mathcal{A}_{B}^{*} + \delta I$ for tiny $\delta$ (purely a measurement perturbation; we are not solving).
- Same under **METIS** (or AMF if METIS not available).
- Cholesky memory estimate: $\text{nnz}(L) \times 8$ bytes.

### 3.4 Range conditioning (for Mazziotti-CG viability)

Already computed via $\kappa_{\text{range}}$ in 3.1. Convert to estimated CG iterations to $10^{-8}$ residual using the standard bound:
$$
k_{\text{CG}} \approx \tfrac{1}{2}\sqrt{\kappa_{\text{range}}}\, \log(2/10^{-8}) \approx 9\sqrt{\kappa_{\text{range}}}.
$$
Tabulate per block. If any single block requires $> 10^{4}$ inner CG iterations, Mazziotti-CG is in trouble for that block.

### 3.5 (Optional / methodology check) TT-rank

Reshape $\mathcal{A}_{B}\mathcal{A}_{B}^{*}$ into a tensor with one mode per orbital, run successive-SVD TT decomposition at tolerance $10^{-6}$, report TT-ranks. Not load-bearing for the solver decision at H₂/Nk=2 (blocks are too small for compression to matter). Validates the methodology so it can be applied at H₄/Nk=2 later. Skip if time-budget tight.

## 4. Symmetry decomposition

Apply, in order of decreasing certainty/importance:

1. **U(1) particle number** $N$. Sectors: fixed $N \in \{0, 1, \ldots, 2M\}$. Already partly handled as a constraint, but should be exploited as a *block*-decomposition.
2. **U(1) spin** $S_{z}$. Sectors: $S_{z} \in \{-N/2, \ldots, +N/2\}$.
3. **Translation $\mathbb{Z}_{Nk} = \mathbb{Z}_{2}$**. Sectors: $k = 0, \pi$. Trivial for Nk=1; gives 2 momentum blocks for Nk=2.
4. **Hermiticity (Re/Im split).** Real-imaginary lift converts an $n \times n$ HPSD block to a $2n \times 2n$ real PSD block.
5. **Time reversal.** Worth checking — if the Hamiltonian is real-symmetric in the chosen basis, this halves dimensions further.

For each combination of (1)–(5), compute the projected $\mathcal{A}_{B}\mathcal{A}_{B}^{*}$ and run §3.

Use the existing symmetry-adapted-basis machinery in the vault (`topics/symmetry-adapted-basis/`) and any blocking already in `demos/h4_periodic_moment_sos.jl` (`--blocking=momentum|spin|none`). Don't reinvent.

## 5. Decision rules — per block

Once §3 numbers are in, classify each $(B, \mathcal{A})$ pair:

| Classification | Trigger | Implication for production solver |
|---|---|---|
| **Direct Cholesky** | $g > 10$ AND $\text{nnz}(L) \times 8\,\text{bytes} <$ working RAM | ADMM/AL with prefactored inner solve. The dream case. |
| **Tikhonov + CG** | $g > 10$ AND fill-in too large for direct factor | Set $\epsilon \approx \sigma_{\hat r}/100$. CG iter $\sim \sqrt{\sigma_{1}/\epsilon}$. |
| **Mazziotti-CG (no $\epsilon$)** | $g \lesssim 1$ (continuous tail) AND $\sqrt{\kappa_{\text{range}}} \lesssim 10^{3}$ | Vanilla CG, no regularization, periodic re-projection. BPSDP-style. |
| **Facial reduction needed** | excess kernel $>$ 10% of $m_{B}$, OR $\hat r / m_{B} < 0.5$ | Block-level facial reduction is on the critical path. Paper 5 work. |
| **Punt to BM** | block too dense for any of the above AND likely low-rank optimum | Burer–Monteiro contingent fallback (gated). |

## 6. Aggregation across blocks

The instance-level decision is a function of the per-block table:

- **Mostly Direct Cholesky** → ADMM/AL family wins, prefactored inner solve. **Default working hypothesis.**
- **Mixed Direct + Mazziotti-CG** → ADMM/AL with adaptive inner solver per block.
- **A non-trivial fraction "Facial reduction needed"** → upstream facial reduction is *required* before any solver, regardless of family. Paper 5 unblocks the rest.
- **Anywhere "Punt to BM"** → BM activates as a contingent fallback for those specific blocks (per the gated contingency in `topics/bm-benign-landscape/`).

## 7. Honest priors I want this experiment to test

These are my (Yusheng's) current beliefs going in. Phase 2 should confirm or kill them.

1. **Most blocks at H₂/Nk=2 will be tiny and gap-clean.** The interesting question is whether *any* block shows continuous-tail behavior. If even one does at this small scale, the production case will too.
2. **Excess kernel will be small but nonzero.** Anticommutation + hermiticity + symmetry should account for *most* redundancies. Anything left tells me what algebraic identities Paper 5 needs to add to the catalog.
3. **Fill-in will be modest.** AMD or METIS should give $\text{nnz}(L) / m^{2} < 0.1$ at this scale. If not, that's news.
4. **The momentum-blocking already in the codebase is most of the symmetry win.** Adding Sz and time reversal will halve and quarter blocks more, but the qualitative gap structure won't change.

If any of these priors is wrong at H₂/Nk=2, that's the most valuable Phase 2 output — a falsification.

## 8. Falsification budget

| Item | Budget |
|---|---|
| H₂/Nk=2 integrals + reference TOML | 1 day |
| Demo script (clone + retarget) | 0.5 day |
| Diagnostic driver `probes/h2_nk2_aastar_diagnostic.jl` | 1 day |
| Symmetry decomposition wiring | 1–2 days |
| Run + write up `output/phase2/h2_nk2/summary.md` | 0.5 day |
| **Total** | **~5 days** |

If any line item slips to weeks, that itself is information about Phase 3 scope.

## 9. Artifact layout

```
test/data/assets/
  h2_chain_nk2_active_2e4o_integrals.txt       # PySCF dump, prerequisite
  h2_chain_nk2_reference.toml                  # HF/MP2/CCSD references

demos/
  h2_periodic_nk2_moment_sos.jl                # ONE new file; build + solve + per-block dump (§12)

probes/
  h2_nk2_aastar_diagnostic.jl                  # consumes …/blocks/, runs §3 across all blocks, writes summary.md

output/phase2/h2_nk2/
  summary.md                                   # aggregated decision table (§5–§6); written by the probe
  plots/sigma_histogram_<block>.png            # spectrum visualization
  solve/
    objective.json                             # opt value vs HF/MP2/CCSD, residuals, KKT gap, wall-clock
    multipliers.json                           # labelled Lagrange multipliers per equality
  blocks/<symmetry-label>/
    A_full.npz                                 # sparse 𝒱_full restricted to this block (§12.3 step 2)
    A_eq.npz                                   # sparse 𝒱_eq  restricted to this block
    labels.json                                # row/col labels + predicted-redundancy counts
    spectrum.json                              # σ_1, σ_hat_r, g, kappa, etc. (§3.1)
    kernel_basis.npz                           # for inspection (§3.2)
    sparsity.json                              # nnz, AMD/METIS fill-in (§3.3)
    Y_eigvals.json                             # primal ⊙ eigenvalues + numerical rank at τ (§12.3 step 5)
```

## 10. What is *out of scope* for Phase 2

**In scope, to be explicit**: Phase 2 builds *and* solves the H₂/Nk=2 SDP once with `BPSDP.jl` on `HAI`, and extracts both the constraint-side $\mathcal{A}\mathcal{A}^{*}$ diagnostics (§3) and the solution-side moment-matrix rank / residual diagnostics (§12.3 step 5). See §12.

**Out of scope**:

- Modifying `BPSDP.jl` internals (preconditioners, CG kernel, penalty schedule). We use it as-is at the pinned commit on `HAI`.
- Implementing facial reduction. Decision output, not the algorithm.
- Implementing Mazziotti-CG inside NCTSSoS. Decision output only; integration is later.
- Running the same diagnostic at H₄/Nk=2 production scale. Phase 2 is the calibration; H₄/Nk=2 is Phase 3.
- Solver bake-offs against COSMO / libsdp / Mosek. BPSDP is the chosen Phase-2 solver; comparative runtime is a separate experiment.
- Modifying any existing file under `demos/` (see §12.1).

## 11. First-day task list

1. Verify on HAI that the H₄/Nk=2 PySCF integral generation pipeline (under `replicate/`) is reproducible.
2. Adapt that pipeline for nk=2, atoms_per_cell=2, active_space (2e, 4o) → produce `h2_chain_nk2_active_2e4o_integrals.txt`.
3. Write `demos/h2_periodic_nk2_moment_sos.jl` from scratch in the style of `demos/h4_periodic_moment_sos.jl` per §12 (do **not** clone any existing demo).
4. Sync the repo to `/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark` and run the demo on `HAI` via `easy-ssh run` with `BPSDP.jl` developed from `/home/ubuntu/BPSDP.jl`. Expect `OPTIMAL` from BPSDP; sanity-check the optimal value against `h2_chain_nk2_reference.toml`. Read off the per-block moment-matrix rank at $\tau = 10^{-10}$ from the dumped `Y_eigvals.json`.
5. Stub `probes/h2_nk2_aastar_diagnostic.jl` with just §3.1 (dense SVD) on the unsymmetrized $\mathcal{A}\mathcal{A}^{*}$, consuming the per-block sparse matrices the demo wrote to `output/phase2/h2_nk2/blocks/`. Get one block's spectrum out.

After step 5 we know the basic pipeline works and can plug in symmetry blocking and the rest of §3.

---

## 12. Implementation contract for the H₂/Nk=2 demo script

This section is the binding spec for the new demo. §2, §9, §10, §11 all defer to it.

### 12.1 Hard rules

- **Do not modify any existing script under `demos/`.** Every current file there — `h2_periodic_moment_sos.jl`, `h4_periodic_moment_sos.jl`, `h4_periodic_moment_sos_export_libsdp.jl`, `dump_h4_periodic_sdp.jl`, `SDPALibsdpExport.jl`, `solve_libsdp_dat_c.py`, `Project.toml`, `Manifest.toml` — is frozen for Phase 2.
- **Exactly one new file**: `demos/h2_periodic_nk2_moment_sos.jl`. Phase 2 introduces no other Julia scripts under `demos/`. Reusable diagnostic code lives in `probes/`; the demo and the probe communicate via on-disk artifacts (no cross-imports between `demos/` and `probes/`).

### 12.2 Style

Match `demos/h4_periodic_moment_sos.jl` line-for-line in shape:

- Top-level docstring with a usage line, an `# Usage from the repository root:` block, and a paragraph documenting the integral file format.
- `Options` struct + `parse_options(argv)` mirroring the H₄ flag set: `--integrals=…`, `--blocking=momentum|spin|none`, `--include-1d` / `--no-1d`, `--paper-spin` / `--spin-singlet`, `--spin-resolved-trace` / `--no-spin-resolved-trace`, `--singlet-s2` / `--no-singlet-s2`, plus the Phase-2-specific flags listed in §12.4.
- Same helpers in spirit, with the same names: `load_integrals_txt`, `active_hf_energy`, `build_h2_nk2_hamiltonian` (analogue of `build_h4_nk2_hamiltonian`), `build_spin_orbitals`, `build_rdm_bases`, `block_key`, `grouped_blocks`, `grouped_blocks_with_identity`, `single_full_system_clique`, `total_electron_constraint`, `trace_constraints`, `spin_resolved_d_trace_constraints`, `singlet_s2_constraint`, `build_h2_pqg_moment_problem`, `constraint_stats`, `print_summary`.
- Same `print_summary(...)` style — header bar, fixed-width left column, one trailing newline per section, then a closing line that names the artifact.

Targets:

- `nk = 2`, `norb = 4` (active spatial orbitals per k-point), `nelec_per_cell = 2`, `total_electrons = 4`, `M = 2 · nk · norb = 16` spin-orbital modes — matching the “2e4o per cell” active space named in the integral filename and the $M=16$ figure in §1.
- Trace constraints derived from the same closed-shell formulas the H₄ demo uses; do not hardcode the H₄ numbers (e.g. $\mathrm{Tr}\,D = N(N-1)/2$ etc.).

### 12.3 What the script must do, in order

Solving is **mandatory**, not optional. The demo always builds, solves, and dumps. There is no `--no-solve` short-circuit.

1. **Build** the H₂/Nk=2 PQG MomentProblem. Identical control flow to `build_h4_pqg_moment_problem`, just retargeted. Reproduce HF active energy as a sanity check.
2. **Materialize** $\mathcal{A}_{\text{full}}$ (Zero + HPSD-coupling) and $\mathcal{A}_{\text{eq}}$ (Zero only) as sparse matrices, per symmetry block selected by `--blocking`. Write each block to disk:
   - `<output-dir>/blocks/<symmetry-label>/A_full.npz`
   - `<output-dir>/blocks/<symmetry-label>/A_eq.npz`
   - `<output-dir>/blocks/<symmetry-label>/labels.json` (which moments index which rows/columns; predicted-redundancy counts).
3. **Run §3 diagnostics** locally for the diagnostic block selected by `--diagnostic-block` (full sweep across blocks lives in `probes/h2_nk2_aastar_diagnostic.jl`). The demo emits at least: $\sigma_1$, $\hat r$ at the configured tolerances, $\sigma_{\hat r}$, $\kappa_{\text{range}}$, gap $g$, nnz, AMD fill-in. Per-block JSON goes to `<output-dir>/blocks/<symmetry-label>/spectrum.json` and `…/sparsity.json`, exactly per §9.
4. **Solve with `BPSDP.jl` on `HAI`** — see §12.6 for the execution contract. Take the JuMP model produced by NCTSSoS's standard MomentProblem-to-JuMP build path, attach `BPSDP.Optimizer` with the parameters in §12.4, call `optimize!`, and **assert** `termination_status(model) == OPTIMAL`. Anything else aborts the run with a non-zero exit so the failure surfaces in `easy-ssh run` logs.
5. **Extract solution-side diagnostics** from the optimal primal $Y^\star$ and dual $\lambda^\star$:
   - **Per HPSD block $b$**: eigenvalues of $Y_b^\star$, **numerical rank at $\tau \in \{10^{-10}, 10^{-12}, 10^{-14}\}$**, smallest non-zero eigenvalue, eigenvalue histogram (saved to `…/Y_eigvals.json`).
   - Optimal objective vs `hf_active`, MP2, CCSD references from `test/data/assets/h2_chain_nk2_reference.toml`.
   - Primal feasibility residual $\lVert \mathcal{A} y^\star - b\rVert_\infty$ and $\lVert \cdot \rVert_2$.
   - Dual residual / KKT gap, BPSDP wall-clock, iteration count, CG-iteration histogram, termination status.
   - Lagrange multipliers per equality constraint, each row labelled by which trace / which symmetry it enforces (so a non-zero multiplier on a tautological constraint flags an excess-kernel suspect).
6. **Aggregate** §3 + step-5 numbers into the per-block decision table per §5–§6. The aggregated table is written by `probes/h2_nk2_aastar_diagnostic.jl` to `output/phase2/h2_nk2/summary.md`; the demo writes only its own per-block JSON/NPZ artifacts.

### 12.4 Phase-2-specific CLI flags

In addition to the H₄ flag set, the new script accepts:

- `--output-dir=PATH` — root for per-block JSON / NPZ. Default: `output/phase2/h2_nk2/`.
- `--rank-tol=τ₁,τ₂,τ₃` — comma-separated overrides for the default $\{10^{-10}, 10^{-12}, 10^{-14}\}$, applied uniformly to step 3 and step 5.
- `--diagnostic-block=<label>` — which symmetry-block label to inline-diagnose in `print_summary`. Default: the largest block in $m_b$.
- BPSDP knobs (each maps to a `BPSDP.Optimizer` keyword; defaults match the Phase-1 smoke-test recipe and are documented in the script header):
  - `--bpsdp-max-iter=N` (default `5_000`)
  - `--bpsdp-cg-max-iter=N` (default `100`)
  - `--bpsdp-mu-update-frequency=N` (default `25`)
  - `--bpsdp-penalty=ρ` (default `0.1`)
  - `--bpsdp-cg-tol=ε` (default `1e-12`)
  - `--bpsdp-obj-tol=ε` (default `1e-8`)
  - `--bpsdp-err-tol=ε` (default `1e-8`)
  - `--bpsdp-print-level=N` (default `1` so wall-clock and iteration progress land in the `easy-ssh run` log)

No `--solve` / `--no-solve` flag. No `--solver` flag. BPSDP is the only solver.

### 12.5 Out-of-scope reaffirmed

The new script must not:

- Modify NCTSSoS source under `src/`.
- Modify `BPSDP.jl` internals.
- Introduce a new solver dependency beyond `BPSDP.jl` and what `demos/Project.toml` already pins.
- Contain logic meant to evolve across multiple Phase 2 / Phase 3 instances; that lives in `probes/h2_nk2_aastar_diagnostic.jl`. The demo only orchestrates one instance.
- Touch any of `demos/{h2_periodic_moment_sos.jl, h4_periodic_moment_sos.jl, h4_periodic_moment_sos_export_libsdp.jl, dump_h4_periodic_sdp.jl, SDPALibsdpExport.jl, solve_libsdp_dat_c.py, Project.toml, Manifest.toml}`.

### 12.6 Execution contract — must run on `HAI`

BPSDP solves Phase 2; nothing about this is local-friendly. The new demo is designed to be invoked via `easy-ssh run`, with `BPSDP.jl` developed from the synced server checkout. Do **not** add `BPSDP.jl` to the repo's root or `demos/` `Project.toml` — keep it as a `Pkg.develop(path="/home/ubuntu/BPSDP.jl")` step inside the launcher.

- **Repo path on server**: `/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`.
- **`BPSDP.jl` checkout on server**: `/home/ubuntu/BPSDP.jl`, branch `multithread`, commit pinned in `archive/PHASE1_h4_libsdp_bpsdp_roundtrip.md` (`865000a` at the time of archival; verify before each run with `git -C /home/ubuntu/BPSDP.jl rev-parse HEAD`).
- **Launch shape** (canonical; mirror the Phase-1 recipe):

  ```bash
  easy-ssh run 'cd /home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark
  tmpdir=$(mktemp -d)
  julia --startup-file=no --project="$tmpdir" -e "using Pkg; \
    Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\"); \
    Pkg.develop(path=pwd()); \
    Pkg.add([\"JuMP\",\"NPZ\",\"JSON3\",\"TOML\",\"Printf\",\"LinearAlgebra\",\"SparseArrays\"]); \
    Pkg.instantiate()"
  julia --startup-file=no --project="$tmpdir" demos/h2_periodic_nk2_moment_sos.jl \
      --output-dir=output/phase2/h2_nk2'
  ```

  The demo's docstring must include this exact incantation so future runs are copy-paste reproducible.
- **Outputs** land under `output/phase2/h2_nk2/` on the server. Pull them back with `easy-ssh sync` (or whatever the team's standard sync verb is) before reading; do not assume the client has filesystem access during the run.
- **No local Julia run path is supported.** If the demo is invoked locally and `BPSDP` cannot be loaded, it must abort early with a clear error pointing at §12.6 — not silently fall back to a different solver.
- **Determinism guard**: pin `BPSDP.Optimizer(guess_type = :zero, …)` so reruns are bit-stable up to BLAS thread non-determinism. Record `Threads.nthreads()` and `BLAS.get_num_threads()` in `output/phase2/h2_nk2/solve/objective.json`.

---

## Cross-references

- Vault: `~/MyBrain/YushengBrain/inbox/sdp_learning.md` — context on why this experiment is the highest-leverage move in the 6-week plan.
- Vault: `~/MyBrain/YushengBrain/topics/admm-sdp-solver/admm-sdp-solver.{md,typ}` — the structural arbiter argument that motivates measuring $\mathcal{A}\mathcal{A}^{*}$.
- Vault: `~/MyBrain/YushengBrain/topics/sdp-solver-remedies-manybody/` — full survey of remedies whose viability this experiment selects between.
- Vault: `~/MyBrain/YushengBrain/topics/fermionic-moment-sos-program/` — Paper 5 (facial reduction) is the upstream cure activated if §6 says so.
- Repo: `TASK.md` — current operational context (server paths, BPSDP roundtrip recipe).
