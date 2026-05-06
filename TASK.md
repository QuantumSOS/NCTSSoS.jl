# Current Task Context

- All testing runs on the server via `easy-ssh`; do not run tests locally unless explicitly asked.
- **Active phase: Phase 2.** Full plan: `PHASE2_PLAN.md`. Phase 1 (H₄ periodic SDP via libsdp / BPSDP.jl with `.dat-c` roundtrip) is **done** and archived to `archive/PHASE1_h4_libsdp_bpsdp_roundtrip.md`.
- **Phase 3 infrastructure (parallel track, not Phase 2 work):** `MomentProblem` enrichment + clean lowering. Spec: `MOMENT_PROBLEM_ENRICHMENT_PLAN.md`. Caches a `MomentLinearData` view on `MomentProblem`, eliminating the BPSDP 1×1-cone explosion at H₂/Nk=2 via `:psd_blocks` and fixing four lowering bugs. Precondition for Phase 3 H₄/Nk=2 production; can land before or alongside the Phase 2 diagnostic.

## Phase 2 — singular $\mathcal{A}\mathcal{A}^{*}$ diagnostic on H₂ / Nk=2

**One sentence**: for each residual symmetry block of the H₂/Nk=2 PQG V2RDM SDP, decide *empirically* which inner linear-solve strategy the production solver should use, by measuring the spectrum, kernel, and sparsity of $\mathcal{A}\mathcal{A}^{*}$. Output is a **decision table**, not a solver. **No SDP gets solved.**

Authoritative spec lives in `PHASE2_PLAN.md`. Quick map:

- **§2 — Prerequisites**: `test/data/assets/h2_chain_nk2_active_2e4o_integrals.txt`, `test/data/assets/h2_chain_nk2_reference.toml`, `demos/h2_periodic_nk2_moment_sos.jl`.
- **§3 — Diagnostics per block**: dense SVD of $\mathcal{A}_B \mathcal{A}_B^{*}$ for both $\mathcal{A}_{\text{full}}$ and $\mathcal{A}_{\text{eq}}$ — $\sigma_1$, $\hat r$ at $\tau \in \{10^{-10}, 10^{-12}, 10^{-14}\}$, $\sigma_{\hat r}$, $\kappa_{\text{range}}$, gap $g$, kernel basis, predicted vs empirical kernel (the **excess kernel** is the load-bearing number), nnz, AMD/METIS fill-in.
- **§4 — Symmetries**: $N$, $S_z$, translation $\mathbb{Z}_2$, Re/Im split, time reversal.
- **§5 — Per-block classifier**: {Direct Cholesky | Tikhonov+CG | Mazziotti-CG | **Facial reduction needed** | Punt to BM}.
- **§9 — Artifact layout**: results land in `output/phase2/h2_nk2/{summary.md,plots/,blocks/}`; driver in `probes/h2_nk2_aastar_diagnostic.jl`.
- **§10 — Out of scope**: solving the SDP, implementing facial reduction, Mazziotti-CG integration, H₄/Nk=2 production, end-to-end COSMO/BPSDP timings.
- **§11 — First-day list**: PySCF integral pipeline reuse → produce h2/nk=2 integrals → smoke-test the cloned demo → stub `probes/h2_nk2_aastar_diagnostic.jl` with §3.1 (dense SVD) on one block.

## Server / local path map (still current)

- Main benchmark repo:
  - local: `/Users/exaclior/QuantumSOS/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`
  - server (`HAI`, `easy-ssh` remote): `/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`
- `BPSDP.jl`: local `/Users/exaclior/QuantumSOS/BPSDP.jl` ↔ server `/home/ubuntu/BPSDP.jl` (`multithread` @ `865000a`).
- `libsdp`: local `/Users/exaclior/QuantumSOS/libsdp` ↔ server `/home/ubuntu/libsdp` (`main` @ `fcd9e34`).
- Replication reference: local `/Users/exaclior/QuantumSOS/replicate` ↔ server `/home/ubuntu/replicate`.

## Python environment policy (still current)

Use `uv` for all Python environment management on `HAI`. Do not use system Python packages, `pip install --user`, or ad-hoc global installs. Phase 2 is mostly Julia; if Python is needed (e.g. PySCF integral generation), follow the same pattern.
