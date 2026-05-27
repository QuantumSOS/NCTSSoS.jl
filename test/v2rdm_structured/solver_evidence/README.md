# Periodic V2RDM solver evidence

Manual HAI-only sweep for `solver-evidence`. See repo-root `TASK.md` and `plan/decisions.md`; those are the contract.

## One-time setup on HAI

```bash
bash test/v2rdm_structured/solver_evidence/setup_pyscf_env.sh
bash test/v2rdm_structured/solver_evidence/setup_julia_env.sh
```

`setup_julia_env.sh` creates an external env at `$HOME/.julia-envs/nctssos-solver-evidence`, develops this checkout plus BPSDP from `/home/ubuntu/BPSDP.jl`, and installs JuMP/COSMO/SCS/JSON3/IterativeSolvers/LinearMaps. The last two are required for COSMO's CG-indirect path. Keeping the env outside the repo prevents later `easy-ssh run` pushes from overwriting solver deps in `Project.toml`.

## First smoke

Run the smallest case first:

```bash
H2_NKS="2" H4_NKS="" SOLVERS="cosmo" SOLVE_SECONDS=1800 \
  bash test/v2rdm_structured/solver_evidence/sweep.sh
```

The driver logs `length(NCTSSoS.orphan_keys(mp))` before lowering. If it is nonzero, stop. Do **not** add an orphan-policy fallback.

## Full sweep

```bash
bash test/v2rdm_structured/solver_evidence/sweep.sh
```

Defaults:

- H2: `Nk = 2 3 4`, active `[2e, 4orb]` per cell
- H4: `Nk = 2 3`, active `[4e, 8orb]` per cell
- solvers: `cosmo,scs,bpsdp`
- integrals: generated on demand under `integrals/`
- output: `output/solver-evidence/{h2,h4}/`

Useful overrides:

```bash
SOLVERS="cosmo,scs" JULIA_THREADS=16 SOLVE_SECONDS=7200 FORCE=1 \
  bash test/v2rdm_structured/solver_evidence/sweep.sh
```

Set `SOLVER_EVIDENCE_JULIA_PROJECT=/path/to/env` if you do not want the default external Julia env.

Pull results from local:

```bash
easy-ssh pull output/solver-evidence
```

Each solver writes `result.json`, `stdout.log`, `stderr.log`; each system gets a rolling `summary.md` with active energy, shift, total per-cell energy, HF, and gap to HF.

## Frozen BPSDP-native SDP export

For BPSDP tuning/preconditioner work, freeze the lowered SDP once and skip the NCTSSoS/JuMP rebuild path in later runs:

```bash
bash test/v2rdm_structured/solver_evidence/export_sdp_instances.sh
```

Defaults export H2/Nk=2 and H4/Nk=2 to:

```text
output/sdp_instances/h2_nk2/bpsdp_native_keep.jls
output/sdp_instances/h4_nk2/bpsdp_native_keep.jls
```

Use `DEPENDENT_ROWS=keep,drop` if you also want BPSDP's QR-reduced row variant materialized at export time. The native payload contains `A`, `b`, `c`, block kinds/dimensions, and metadata; load it with `include("test/v2rdm_structured/solver_evidence/bpsdp_native_io.jl")` and `BPSDPNativeIO.read_native(path)`.
