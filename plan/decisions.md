# Recorded Decisions

Decisions are **locked**. Anything that would change one of these needs to update this file first, then propagate.

**See also:** [`findings/`](findings/) — empirical observations that may force new decisions.
- [`findings/residual-norms.md`](findings/residual-norms.md) — COSMO/SCS use ∞-norm residuals; BPSDP uses 2-norm. Affects D5 tolerance interpretation; D8 pending.

## D1 — `orphan_policy` → leave at default `:error`, no fallback

The canonical lowering call is:

```julia
build_jump_model(mp; formulation = :psd_blocks, representation = :complex)
# orphan_policy defaults to :error — DO NOT pass any value.
```

If the build emits orphans, the throw is the **signal**: the formulation is wrong (likely some constraint we shouldn't be adding, or a missing pivot in the V2RDM block structure). The driver must not silently fall back to `:free_variables` or `:aux_psd_free`.

The driver should still log `length(NCTSSoS.orphan_keys(mp))` before lowering, so the failure mode is one line of context rather than a raw stack trace.

The reference repo's use of `:aux_psd_free` (even with the explicit total-N row dropped) is a *symptom they patched over*, not a recipe to copy.

## D2 — H₂ at Nk ∈ {3} via cloned PySCF script

Yes. `test/v2rdm_structured/solver_evidence/extract_h2_nk_integrals.py` cloned from the H₄ generator:
- 2 H atoms at x ∈ {0, 1} Å, lattice 2.0 Å.
- Active space [2 e⁻, 4 orb] per cell.
- `pseudo="none"` (matches existing H₂/Nk=2 reference).
- Otherwise identical to `extract_h4_nk_integrals.py`.

## D3 — Per-system Nk caps

Sweep range differs per system:
- **H₂: Nk ∈ {2, 3, 4}** (small system; cheap; 3 cases).
- **H₄: Nk ∈ {2, 3}** (large system; 2 cases).

Five runs total. No H₄/Nk=4, no Nk=5 for either system. The earlier "Nk=2..5" scope is retired.

## D4 — Active space unchanged

- H₄: [4 e⁻, 8 orb] / cell.
- H₂: [2 e⁻, 4 orb] / cell.

## D5 — Three solvers, indirect / CG paths only

All three solvers run on every (system, Nk) case. Each uses its *iterative* inner-linear-solve path; "CG-indirect" is COSMO-specific terminology and does **not** apply to SCS or BPSDP.

### COSMO with `CGIndirectKKTSolver`

```julia
using IterativeSolvers, LinearMaps   # Requires.jl trigger
kkt = COSMO.with_options(COSMO.CGIndirectKKTSolver;
    tol_constant = 1.0, tol_exponent = 1.5)
JuMP.optimizer_with_attributes(COSMO.Optimizer,
    "kkt_solver" => kkt,
    "max_iter" => …, "eps_abs" => …, "eps_rel" => …,
    "verbose" => false, "verbose_timing" => false,
    "decompose" => false, "complete_dual" => false)
```

### SCS with `IndirectSolver`

```julia
import SCS
SCS.is_available(SCS.IndirectSolver) || error("SCS.IndirectSolver missing in this build")
JuMP.optimizer_with_attributes(SCS.Optimizer,
    "linear_solver" => SCS.IndirectSolver,
    "max_iters" => …, "eps_abs" => …, "eps_rel" => …,
    "verbose" => 0)
```

SCS's `IndirectSolver` is its CG-based linear solver. There is no `CGIndirectKKTSolver` knob in SCS — the indirect solver *is* the CG path.

### BPSDP.jl

Developed against `/home/ubuntu/BPSDP.jl` on HAI (branch `main`, head `3edebb9` at confirmation time):

```julia
import BPSDP
BPSDP.Optimizer(
    max_iter = …, cg_max_iter = …,
    mu_update_frequency = 500, penalty_parameter = 0.1,
    cg_convergence = 1e-10, dynamic_cg_convergence = true,
    sdp_objective_convergence = 1e-5, sdp_error_convergence = 1e-5,
    guess_type = :zero, print_level = 1,
    progress_monitor = monitor,   # captures inner CG iterations per outer step
    dependent_rows = :keep)
```

BPSDP is BPSDP — it has its own inner CG; there is no "direct" alternative inside it.

### Package envelope

```bash
julia --project=. -e 'using Pkg;
    Pkg.develop(path="/home/ubuntu/BPSDP.jl");
    Pkg.add(["JuMP", "COSMO", "SCS", "JSON3", "IterativeSolvers", "LinearMaps"]);
    Pkg.instantiate(); Pkg.precompile()'
```

`IterativeSolvers` + `LinearMaps` are the trigger packages for COSMO's optional `CGIndirectKKTSolver`. Both must be in the env or COSMO silently uses its default direct KKT factorization.

### Driver responsibilities

The driver loops over the requested solver set per case. CLI: `--solvers=cosmo,scs,bpsdp` (default all three). Each solver gets a **fresh JuMP lowering** of the same `MomentProblem`; do not reuse a model across solvers.

## D6 — Integrals computed on the fly, not pre-staged

The script's primary purpose is documentation / reproducibility. The workflow is:

1. User clones, sets up the PySCF venv once (`setup_pyscf_env.sh`).
2. Sweep launches; for each `(system, Nk)` the launcher checks whether the integral text + meta JSON exist.
3. If missing, the launcher shells out to the PySCF extractor and waits.
4. Then build, solve, write `result.json`, append summary row.

Sizing sanity (no large files):
- H₂/Nk=3, 4 orb, 1 spin: ~2 MB text dump.
- H₄/Nk=3, 8 orb, 1 spin: ~15 MB text dump.

These are small enough to regenerate per-run; they're still `.gitignore`-d under a project-relative `integrals/` directory so a CI clone stays clean.

## D7 — Summary table reports both objective forms

Each row in `output/solver-evidence/<system>/summary.md` carries:

| column | source |
|---|---|
| `active_objective_Ha` | `JuMP.objective_value(model)` |
| `energy_shift_Ha` | `meta.json["energy_shift"]` from the PySCF run |
| `total_per_cell_Ha` | `active_objective + energy_shift` |
| `ehf_Ha` | `meta.json["ehf"]` |
| `gap_to_hf_Ha` | `total_per_cell - ehf` (V2RDM should give a *lower* bound, so this should be ≤ 0) |

## Eventual docs port (recorded so it's not forgotten)

`docs/src/examples/literate/v2rdm_periodic_solver_evidence.jl`:
- Live-executes H₂/Nk=2 end-to-end (smallest case; cheap; demonstrates the full path).
- Embeds H₂/Nk=3 + H₄/Nk=2,3 as a static markdown table read from a committed `summary.md`. Doc build doesn't re-solve those.
- Pulls integrals on the fly via the same PySCF venv at doc-build time only for the live H₂/Nk=2 case; or, if `make examples` is too brittle for that, falls back to a pre-staged Nk=2 integrals file committed under `docs/src/examples/data/`.

## Outstanding (operational — already in motion)

1. **HAI sync.** `easy-ssh` is configured (`.easy-ssh.conf` host `HAI`, remote `~/NCTSSoS.jl-main`). HAI confirmed reachable; remote directory will be created by the first `easy-ssh push`. Julia 1.12.6 and uv 0.11.11 are present. BPSDP.jl is at `/home/ubuntu/BPSDP.jl` (branch `main`, head `3edebb9`).
2. **End-to-end `build_pqg_moment_data` on real PySCF data.** The toy case in `test/v2rdm_structured/runtests.jl` exercises the API shape; the first real-data invocation is H₂/Nk=2 in the sweep itself.
3. **`length(orphan_keys(mp))` empirical result.** Per D1, this is a hard validation step: if it's nonzero, we stop and re-examine the formulation rather than papering over with `:aux_psd_free`.
