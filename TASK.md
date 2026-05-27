# Current Task

Branch `solver-evidence`: add a concise H₂/H₄ periodic V2RDM solver-performance script that runs on the HAI server and sweeps **H₂ at Nk ∈ {2, 3, 4}** and **H₄ at Nk ∈ {2, 3}** (5 cases total). The artifact is primarily for documentation; the test-stage code is a stepping stone to a Literate example page.

## Where it lives

1. **First stage — `test/v2rdm_structured/solver_evidence/`** (excluded from `runtests.jl`; runnable manually on HAI).
2. **Later port — `docs/src/examples/literate/v2rdm_periodic_solver_evidence.jl`** once the test-stage script produces stable numbers.

The `demos/` directory has been removed; do not recreate it.

## Locked decisions

All seven open decisions are resolved. Detail in [`plan/decisions.md`](plan/decisions.md):

- **D1** Lowering call is `build_jump_model(mp; formulation=:psd_blocks, representation=:complex)` — **no `orphan_policy` kwarg**. `:error` (default) must throw if orphans appear; that throw is the alert.
- **D2** Add `extract_h2_nk_integrals.py` cloned from the H₄ generator (2 atoms, lattice 2 Å, [2e, 4orb] active).
- **D3** Per-system Nk caps: **H₂ Nk ∈ {2, 3, 4}**, **H₄ Nk ∈ {2, 3}**. Five cases total.
- **D4** Active spaces unchanged: H₄ [4e, 8orb]/cell; H₂ [2e, 4orb]/cell.
- **D5** Three solvers, each on the iterative inner-linear-solve path: **COSMO with `CGIndirectKKTSolver`**, **SCS with `IndirectSolver`**, **BPSDP.jl** (developed from `/home/ubuntu/BPSDP.jl`). Each gets a fresh JuMP lowering per case.
- **D6** Integrals **computed on the fly** by the sweep; not pre-staged, not committed. Sizes are small (≤15 MB per case).
- **D7** Summary reports `active_objective_Ha`, `energy_shift_Ha`, `total_per_cell_Ha`, `ehf_Ha`, `gap_to_hf_Ha`.

## Detail

- [`plan/decisions.md`](plan/decisions.md) — **start here**: the locked record of D1–D7
- [`plan/overview.md`](plan/overview.md) — goal, scope, deliverables
- [`plan/script-layout.md`](plan/script-layout.md) — file structure to add under `test/`
- [`plan/reference-mock.md`](plan/reference-mock.md) — what to copy and what to cut from the reference repo
- [`plan/server.md`](plan/server.md) — HAI / `easy-ssh` / `uv` conventions
- [`plan/findings/residual-norms.md`](plan/findings/residual-norms.md) — **must read before interpreting solver-residual numbers**: COSMO/SCS report ‖·‖∞, BPSDP reports ‖·‖₂. Same numeric tolerance does not equal same precision.

## Outstanding before the first run

1. HAI is reachable; `.easy-ssh.conf` and `.easy-ssh-ignore` are in place. The first `easy-ssh push` creates `/home/ubuntu/NCTSSoS.jl-main` and primes the remote checkout.
2. First sweep run on H₂/Nk=2 doubles as the end-to-end validation of `build_pqg_moment_data` against real PySCF integrals.
3. Per D1: if `length(orphan_keys(mp)) > 0`, stop and reconsider the formulation. Do not patch around.
