# Script Layout

Target: a single shared driver plus thin launchers under `test/v2rdm_structured/solver_evidence/`. **No `demos/` directory.** The Julia driver is the only non-trivial file.

## Files to add

Top-level:

```
.easy-ssh.conf      # host=HAI, remote_dir=~/NCTSSoS.jl-main
.easy-ssh-ignore    # excludes .git, integrals/, output/, .julia/, tmp/, .DS_Store
```

Solver-evidence tree:

```
test/v2rdm_structured/solver_evidence/
  driver.jl                    # CLI, build, lower, solve {cosmo,scs,bpsdp}, write JSON+md
  h2_launcher.jl               # launcher: defaults for H₂ (norb=4, nelec_per_cell=2)
  h4_launcher.jl               # launcher: defaults for H₄ (norb=8, nelec_per_cell=4)
  sweep.sh                     # bash sweep; per-system Nk list; integrals on demand
  setup_pyscf_env.sh           # one-shot uv-managed PySCF env recipe
  setup_julia_env.sh           # Pkg.develop BPSDP.jl + Pkg.add the rest
  extract_h4_nk_integrals.py   # verbatim from reference repo probes/
  extract_h2_nk_integrals.py   # cloned-and-trimmed per D2
  README.md                    # how to run on HAI; pointers to decisions.md
```

## Wiring into the existing repo

- **Do not** include these files from `test/runtests.jl` or `test/v2rdm_structured/runtests.jl`. They are not unit tests.
- Document the invocation in `README.md`.
- `.gitignore` already excludes `integrals/` and `output/`.
- `.easy-ssh-ignore` excludes them from sync too — solver outputs come back via `easy-ssh pull`, not via background push.

## Shared driver — required surface

CLI flags:

| flag | purpose |
|---|---|
| `--system=h2\|h4` | system selector |
| `--nk=N` | k-points |
| `--norb=N` | active spatial orbitals; defaulted by launcher |
| `--nelec-per-cell=N` | active electrons per cell; defaulted by launcher |
| `--integrals=PATH` | required: PySCF text dump |
| `--meta=PATH` | required: companion `_meta.json` (gives `ehf` and `energy_shift`) |
| `--output-dir=PATH` | required |
| `--solvers=cosmo,scs,bpsdp` | comma-separated subset; default `cosmo,scs,bpsdp` |
| `--time-limit=SECS` | wall time cap per solver |
| `--repeats=N` | default 1 |

**Per-solver overrides (kept narrow):**

| flag | applies to | default |
|---|---|---|
| `--cosmo-max-iter=N` | COSMO | 200_000 |
| `--cosmo-eps-abs=E` `--cosmo-eps-rel=E` | COSMO | 1e-5 |
| `--cosmo-tol-constant=X` `--cosmo-tol-exponent=X` | COSMO CG-KKT | 1.0, 1.5 |
| `--scs-max-iters=N` | SCS | 200_000 |
| `--scs-eps-abs=E` `--scs-eps-rel=E` | SCS | 1e-5 |
| `--bpsdp-max-iter=N` | BPSDP | 50_000 |
| `--bpsdp-cg-max-iter=N` | BPSDP | 10_000 |
| `--bpsdp-obj-tol=E` `--bpsdp-err-tol=E` | BPSDP | 1e-5 |
| `--bpsdp-cg-tol=E` | BPSDP | 1e-10 |

**No `--orphan-policy` flag.** Lowering hardcodes `orphan_policy = :error` (D1).
**No `--formulation` / `--representation` flags.** Hardcoded to `(:psd_blocks, :complex)`.

## Build + lower step (canonical, built once per case)

```julia
# Parse integrals text → h1e::Dict{Int,Matrix{ComplexF64}},
#                       eri::Dict{NTuple{4,Int},Array{ComplexF64,4}}
# (mirroring the toy parser in test/v2rdm_structured/runtests.jl)

linear = build_pqg_moment_data(h1e, eri;
    nk             = opts.nk,
    norb           = opts.norb,
    nelec_per_cell = opts.nelec_per_cell,
    blocking            = :momentum,
    spin_resolved_trace = true,
    singlet_s2          = true,
    include_one_d       = false,
)

mp = NCTSSoS.MomentProblem(...)   # assembled exactly as the reference repo's helper does
n_orphans = length(NCTSSoS.orphan_keys(mp))
@info "moment-problem build" n_orphans hpsd_block_sizes total_canonical_moments

# n_orphans must be 0. Lowering throws otherwise (default :error policy, D1).
```

Build the `MomentProblem` once per (system, Nk). Lower it **separately** per solver — each solver call goes through `build_jump_model(...)` fresh, so JuMP model state never leaks across solvers.

## Solve step — three branches

```julia
function solve_cosmo!(model, opts)
    using IterativeSolvers, LinearMaps   # Requires.jl trigger
    import COSMO
    kkt = COSMO.with_options(COSMO.CGIndirectKKTSolver;
        tol_constant = opts.cosmo_tol_constant,
        tol_exponent = opts.cosmo_tol_exponent,
    )
    JuMP.set_optimizer(model, JuMP.optimizer_with_attributes(COSMO.Optimizer,
        "kkt_solver"             => kkt,
        "max_iter"               => opts.cosmo_max_iter,
        "eps_abs"                => opts.cosmo_eps_abs,
        "eps_rel"                => opts.cosmo_eps_rel,
        "time_limit"             => opts.time_limit,
        "verbose"                => false, "verbose_timing" => false,
        "decompose"              => false, "complete_dual" => false,
    ))
    JuMP.set_silent(model)
end

function solve_scs_indirect!(model, opts)
    import SCS
    SCS.is_available(SCS.IndirectSolver) || error("SCS.IndirectSolver missing")
    JuMP.set_optimizer(model, JuMP.optimizer_with_attributes(SCS.Optimizer,
        "linear_solver"    => SCS.IndirectSolver,
        "max_iters"        => opts.scs_max_iters,
        "eps_abs"          => opts.scs_eps_abs,
        "eps_rel"          => opts.scs_eps_rel,
        "time_limit_secs"  => opts.time_limit,
        "verbose"          => 0,
    ))
end

function solve_bpsdp!(model, opts, cg_iters::Vector{Int})
    import BPSDP
    monitor = (lvl, outer, inner, _...) -> (push!(cg_iters, Int(inner)); nothing)
    JuMP.set_optimizer(model, () -> BPSDP.Optimizer(
        max_iter                 = opts.bpsdp_max_iter,
        cg_max_iter              = opts.bpsdp_cg_max_iter,
        mu_update_frequency      = 500,
        penalty_parameter        = 0.1,
        cg_convergence           = opts.bpsdp_cg_tol,
        dynamic_cg_convergence   = true,
        sdp_objective_convergence = opts.bpsdp_obj_tol,
        sdp_error_convergence    = opts.bpsdp_err_tol,
        guess_type               = :zero,
        print_level              = 1,
        progress_monitor         = monitor,
        dependent_rows           = :keep,
    ))
end
```

After `optimize!`, mine the raw optimizer for solver-specific state:
- **COSMO** — `raw.results.info.{r_prim, r_dual}`, `raw.results.iter`, `raw.inner.kkt_solver.indirect_kktsolver.multiplications` (total CG matvecs), `raw.results.times.solver_time`.
- **SCS** — `raw.sol.{iterations, solve_time_sec, objective_value, dual_objective_value, raw_status, ret_val}`.
- **BPSDP** — `raw.state.{outer_iterations, primal_error, dual_error, objective_primal, objective_dual, termination_reason}`; plus the monitor's `cg_iters` vector.

## Per-run output

`output/solver-evidence/<system>/nk<N>/<solver>/result.json`:

- run identity: `system`, `nk`, `solver`, `norb`, `nelec_per_cell`, hostname, julia/blas threads, git rev, started_at, finished_at
- problem stats (same across solvers for the same case): hpsd block sizes, total canonical moments, n_orphans (must be 0)
- pyscf metadata: `ehf_Ha`, `energy_shift_Ha`, basis, density_fit
- solve: `lowering_wall_seconds`, `solve_wall_seconds`, iteration count, solver-specific inner-iter total
- objective (D7): `active_objective_Ha`, `total_per_cell_Ha`, `gap_to_hf_Ha`
- residuals: `r_prim`, `r_dual`
- status: `termination_status`, `raw_status`, optional `error`

`output/solver-evidence/<system>/summary.md` — one row per (Nk, solver):

| Nk | solver | status | iter | inner iter | solve s | active Ha | total Ha | HF Ha | gap to HF Ha | r_prim | r_dual |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|

Solver stdout/stderr go to bounded log files per solver (head + tail, truncate >256 KiB).

## Sweep script — required behavior

`test/v2rdm_structured/solver_evidence/sweep.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Per-system Nk lists (D3).
H2_NKS="${H2_NKS:-2 3 4}"
H4_NKS="${H4_NKS:-2 3}"

OUTROOT="${OUTROOT:-output/solver-evidence}"
INTROOT="${INTROOT:-integrals}"
JULIA_THREADS="${JULIA_THREADS:-8}"
SOLVE_SECONDS="${SOLVE_SECONDS:-3600}"
PYSCF_PYTHON="${PYSCF_PYTHON:-$HOME/.venvs/pyscf-h/bin/python}"
SOLVERS="${SOLVERS:-cosmo,scs,bpsdp}"

run_case() {
    local system=$1 nk=$2
    local intdir="$INTROOT/${system}_chain_nk${nk}"
    local txt="$intdir/${system}_chain_nk${nk}_integrals_drop1e-12.txt"
    local meta="$intdir/${system}_chain_nk${nk}_meta.json"

    if [[ ! -f "$txt" || ! -f "$meta" ]]; then
        mkdir -p "$intdir"
        "$PYSCF_PYTHON" test/v2rdm_structured/solver_evidence/extract_${system}_nk_integrals.py \
            --nk="$nk" --outdir="$intdir" --drop-tol=1e-12
    fi

    JULIA_NUM_THREADS="$JULIA_THREADS" \
    julia --startup-file=no --project=. \
        test/v2rdm_structured/solver_evidence/${system}_launcher.jl \
            --nk="$nk" --integrals="$txt" --meta="$meta" \
            --solvers="$SOLVERS" --time-limit="$SOLVE_SECONDS" \
            --output-dir="$OUTROOT/$system/nk$nk"
}

for nk in $H2_NKS; do run_case h2 "$nk"; done
for nk in $H4_NKS; do run_case h4 "$nk"; done
```

Resume-friendly: a case skips an already-completed solver when its `result.json` is present and `--force` is not set.

## What NOT to add

- No per-Nk Julia launchers.
- No `--orphan-policy`, `--formulation`, `--representation` flags.
- No solver other than COSMO / SCS / BPSDP (D5).
- No `demos/` directory.
- No plotting, comparison plots, or classifier code. Tables only.

## Docs port (later PR)

When the test-stage `summary.md` stabilizes:
- Create `docs/src/examples/literate/v2rdm_periodic_solver_evidence.jl`.
- Live-execute H₂/Nk=2 with **COSMO only** (smallest case, single solver) for the doc-build pass.
- Embed the full 9 + 6 = 15-row summary table as static markdown from the committed `summary.md`.
- Update the example TOC and regenerate via `make examples`.
