# Overview — H₂/H₄ Solver Performance Sweep

## Goal
Produce solver-performance evidence on `solver-evidence` for the periodic V2RDM SDPs built by NCTSSoS, sweeping system × k-points × solver:
- **H₂ at Nk ∈ {2, 3, 4}**, **H₄ at Nk ∈ {2, 3}** (5 cases) — D3
- **Three solvers per case** — D5:
  - COSMO with `CGIndirectKKTSolver`
  - SCS with `IndirectSolver`
  - BPSDP.jl (from `/home/ubuntu/BPSDP.jl`)
- One shared driver, one shell sweep, one summary table per system

The deliverable is primarily a documentation/reproducibility artifact. It will be ported into a Literate example page once the test-stage code is verified.

## In scope
- Build the V2RDM `MomentProblem` via `build_pqg_moment_data` (same API exercised by `test/v2rdm_structured/runtests.jl`). Drop the explicit total-N row, matching the reference repo.
- Lower with:
  ```julia
  build_jump_model(mp; formulation = :psd_blocks, representation = :complex)
  # orphan_policy defaults to :error — DO NOT pass any value. (D1)
  ```
- Solve with **each** of the three solvers, one after another, each from a *fresh* JuMP lowering of the same `MomentProblem`. Solver-specific construction lives in `plan/decisions.md` D5; the canonical code shapes are inlined in `plan/script-layout.md`.
- Generate integrals on the fly per (system, Nk) before solving (D6).
- Capture per-run (run = one solver on one case): wall time (PySCF + build + solve), termination status, `active_objective_Ha`, `total_per_cell_Ha`, `ehf_Ha`, `gap_to_hf_Ha`, primal/dual residuals, iteration counts, solver-specific iterative-solve counts (COSMO CG matvecs, BPSDP inner CG iterations, SCS iters), peak RSS, and `length(orphan_keys(mp))`.
- Emit one JSON per (case, solver) + one rolling `summary.md` table per system with one row per (Nk, solver).

## Out of scope
- Any other `orphan_policy` than `:error` (D1).
- Mosek, Hypatia, ProxSDP, or any other solver outside the three listed (D5).
- Direct-KKT or other non-iterative inner solves on COSMO / SCS.
- H₂ Nk=5 or H₄ Nk ∈ {4, 5} (D3).
- Variation of active space (D4).
- Plotting, decision tables, facial-reduction work, library refactors.
- Local execution. Everything runs on HAI.

## Reference values
- **HF/cell** is the per-(system, Nk) target. Emitted by the PySCF generator into `<integrals_dir>/_meta.json` as `ehf`, alongside `energy_shift = ehf − active_hf`.
- The summary table reports `gap_to_hf_Ha = total_per_cell − ehf`. PQG V2RDM is a lower bound on the active-Hamiltonian ground state, so `gap_to_hf_Ha ≤ 0` is the expected sign.
- The reference repo's `h4_chain_reference.toml` (digitized HF/MP2/CCSD/V2RDM table) is a sanity cross-check only.

## Staging — two-step landing

1. **Test stage (this PR):** code at `test/v2rdm_structured/solver_evidence/`, excluded from `test/runtests.jl`. Invocable on HAI via `bash test/v2rdm_structured/solver_evidence/sweep.sh`.
2. **Docs port (follow-up PR):** Literate page at `docs/src/examples/literate/v2rdm_periodic_solver_evidence.jl`. Live-executes only H₂/Nk=2; embeds the H₂/Nk=3 + H₄ rows as a static table from the committed `summary.md`.

The `demos/` directory has been removed and must not be recreated.

## Deliverables (test stage)

```
.easy-ssh.conf                 # host=HAI, remote_dir=~/NCTSSoS.jl-main
.easy-ssh-ignore               # excludes .git, integrals/, output/, etc.

test/v2rdm_structured/solver_evidence/
  driver.jl                    # shared CLI, build, lower, solve each of {cosmo,scs,bpsdp}, write JSON+md
  h2_launcher.jl               # thin: defaults for H₂ (norb=4, nelec_per_cell=2)
  h4_launcher.jl               # thin: defaults for H₄ (norb=8, nelec_per_cell=4)
  sweep.sh                     # loop over per-system Nk lists; generates integrals on demand
  setup_pyscf_env.sh           # one-shot uv recipe for the PySCF env on HAI
  setup_julia_env.sh           # Pkg.develop BPSDP.jl + Pkg.add the rest
  extract_h4_nk_integrals.py   # ported verbatim from reference repo probes/
  extract_h2_nk_integrals.py   # cloned-and-trimmed from extract_h4 per D2
  README.md                    # how to run on HAI; pointers to TASK.md and decisions.md
```

Generated artifacts (all `.gitignore`-d and excluded from `easy-ssh` sync):

```
integrals/{h2,h4}_chain_nk{N}/
  {h2,h4}_chain_nk{N}_integrals_drop1e-12.txt
  {h2,h4}_chain_nk{N}_meta.json
  {h2,h4}_chain_nk{N}.npz

output/solver-evidence/{h2,h4}/
  nk{N}/{cosmo,scs,bpsdp}/result.json
  nk{N}/{cosmo,scs,bpsdp}/stdout.log
  nk{N}/{cosmo,scs,bpsdp}/stderr.log
  summary.md
```

## Success criterion
From a fresh HAI checkout: `easy-ssh push` → one-shot env setup → `easy-ssh run "bash test/v2rdm_structured/solver_evidence/sweep.sh"` → `easy-ssh pull output/solver-evidence`. The pulled `summary.md` per system shows:
- H₂: 3 cases (Nk=2,3,4) × 3 solvers = 9 rows.
- H₄: 2 cases (Nk=2,3) × 3 solvers = 6 rows.

Each row carries solver, termination status, iterations / inner CG counts, solve wall time, both objective columns, residuals, and `gap_to_hf_Ha`. A reader compares the three solvers on each case and decides which converges, which doesn't, and at what cost.
