# HAI Server Conventions

All benchmark runs execute on the HAI server. Local copy is for editing and pushing only.

## Verified server state (confirmed 2026-05-14)

- Host `HAI` reachable via plain SSH (entry already in `~/.ssh/config`).
- Julia 1.12.6 at `/home/ubuntu/.juliaup/bin/julia`.
- `uv` 0.11.11 at `/home/ubuntu/.local/bin/uv`.
- BPSDP.jl checkout at `/home/ubuntu/BPSDP.jl` (branch `main`, head `3edebb9`).
- `libsdp` at `/home/ubuntu/libsdp` (head `fcd9e34`).
- Reference repo at `/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`.
- **`/home/ubuntu/NCTSSoS.jl-main` does not yet exist.** The first `easy-ssh push` creates it.

## easy-ssh wiring (already in repo)

```
.easy-ssh.conf
  host='HAI'
  remote_dir='~/NCTSSoS.jl-main'

.easy-ssh-ignore
  .git, .pi, .DS_Store, .vscode, .julia, *.jl.cov, *.jl.mem,
  docs/build, docs/site, integrals/, output/, tmp/
```

Common commands:

```bash
easy-ssh push                    # additive sync; respects .easy-ssh-ignore
easy-ssh push --clean --force    # destructive sync (drops remote files not in local)
easy-ssh run  "<cmd>"            # push, run synchronously
easy-ssh submit "<cmd>"          # push, launch via nohup, detach
easy-ssh logs                    # tail of last submitted job
easy-ssh monitor                 # live tail
easy-ssh pull <path>             # fetch a path back into the local tree
easy-ssh status                  # config / SSH / job status
```

Solver outputs and PySCF integrals stay on the server until you `easy-ssh pull output/solver-evidence` (or pull a single `result.json`).

## Bootstrap on HAI

Two one-shot setup scripts handle the env. Run each once after the first push.

### `setup_pyscf_env.sh` (uv)

```bash
#!/usr/bin/env bash
set -euo pipefail
VENV="${VENV:-$HOME/.venvs/pyscf-h}"
uv venv "$VENV"
uv pip install --python "$VENV/bin/python" pyscf numpy
echo "PySCF venv ready at $VENV"
echo "Export PYSCF_PYTHON=$VENV/bin/python before running sweep.sh"
```

Rules:
- `uv` only. No system Python, no `pip install --user`, no `/tmp` venvs.
- The venv is created once and reused.

### `setup_julia_env.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
julia --startup-file=no --project=. -e '
    using Pkg;
    Pkg.develop(path="/home/ubuntu/BPSDP.jl");
    Pkg.add(["JuMP", "COSMO", "SCS", "JSON3", "IterativeSolvers", "LinearMaps"]);
    Pkg.instantiate();
    Pkg.precompile()'
```

- `IterativeSolvers` + `LinearMaps` are the Requires.jl trigger for COSMO's `CGIndirectKKTSolver`. Both must be in the env or COSMO silently uses the default direct KKT factorization.
- BPSDP is `Pkg.develop`-ed from the local path on HAI.

## Integrals — generated on the fly (D6)

The sweep script auto-generates missing integrals before each (system, Nk) run:

```bash
intdir="integrals/${system}_chain_nk${nk}"
if [[ ! -f "$intdir/${system}_chain_nk${nk}_meta.json" ]]; then
    "$PYSCF_PYTHON" .../extract_${system}_nk_integrals.py --nk="$nk" --outdir="$intdir"
fi
```

Sizes (max):
- H₂/Nk=4, [2e,4o]: ~3 MB text dump.
- H₄/Nk=3, [4e,8o]: ~15 MB text dump.

Both fit under `integrals/`, which is `.gitignore`-d **and** excluded from `easy-ssh` sync.

## Output discipline

- Per-run JSON + bounded solver logs land under `output/solver-evidence/<system>/nk<N>/<solver>/`.
- Output stays on the server. Pull only what you need: `easy-ssh pull output/solver-evidence/h2/summary.md`.
- The committed artifact, if any, is the rolling `summary.md` per system, and only after the numbers stabilize.

## End-to-end first run (recipe)

From the local checkout:

```bash
# 1. Sync the repo to HAI (creates ~/NCTSSoS.jl-main on first run).
easy-ssh push

# 2. One-shot env setup on HAI.
easy-ssh run "bash test/v2rdm_structured/solver_evidence/setup_pyscf_env.sh"
easy-ssh run "bash test/v2rdm_structured/solver_evidence/setup_julia_env.sh"

# 3. Submit the sweep (long-running; detaches).
easy-ssh submit "PYSCF_PYTHON=\$HOME/.venvs/pyscf-h/bin/python \
                 JULIA_NUM_THREADS=8 \
                 bash test/v2rdm_structured/solver_evidence/sweep.sh"

# 4. Watch.
easy-ssh monitor

# 5. Pull results once it ends.
easy-ssh pull output/solver-evidence
```

Result lives in `output/solver-evidence/h{2,4}/summary.md`.
