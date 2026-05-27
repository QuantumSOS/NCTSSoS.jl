# Reference Mock — What to Copy, What to Cut

Source: `/Users/exaclior/QuantumSOS/NCTSSoS.jl-h4-periodic-v2rdm-benchmark/`. Destination here is **`test/v2rdm_structured/solver_evidence/`**, not `demos/`. The `demos/` directory has been removed from this repo and must not be recreated.

## Port verbatim

- `probes/extract_h4_nk_integrals_pyscf.py` → `test/v2rdm_structured/solver_evidence/extract_h4_nk_integrals.py`. Drop in unchanged. Already emits `ehf`, `active_hf`, `energy_shift` in `_meta.json` — exactly what our summary table reads.
- Text-format spec from that script: `h1e k p q Re Im; eri k1 k2 k3 k4 p r q s Re Im`, chemist order, unnormalized (Julia builder divides by Nk and Nk²).
- Reference integrals at `test/data/assets/h2_chain_nk2_active_2e4o_integrals.txt` are interesting as a known-good cross-check but **are not used** in this branch — we always regenerate via PySCF (D6).

## Clone-and-edit

- `extract_h4_nk_integrals.py` → `extract_h2_nk_integrals.py` (D2):
  - 4 H atoms at x ∈ {0,1,2,3} Å → **2 H atoms at x ∈ {0,1} Å**.
  - `cell.a` first axis 4.0 Å → **2.0 Å**.
  - `--n-active-orb` default 8 → **4**.
  - `--n-active-elec` default 4 → **2**.
  - Stem `h4_chain` → **`h2_chain`**.
  - Keep `pseudo="none"` default (matches existing H₂ reference).
  - Everything else (GDF, `dimension=1`, `low_dim_ft_type=inf_vacuum`, `gth-tzvp`, ERI extraction loop) unchanged.

## Solver wiring — copy these bits from `demos/h4_periodic_nk2_solver_benchmark.jl`

The reference repo's `h4_periodic_nk2_solver_benchmark.jl` has the canonical three-solver wiring. Port these helpers verbatim (or near-verbatim) per D5:

- **COSMO CG-indirect.** `cosmo_cg_optimizer_factory` shape. Includes `using IterativeSolvers; using LinearMaps` at the top of the file (Requires.jl trigger). Plus `optimize_with_redirected_stdio!` + `truncate_log_file!` helpers — COSMO's MOI wrapper can dump the whole JuMP model on failure; muzzle it. Plus `JuMP.set_silent(model)` after `set_optimizer`.
- **SCS indirect.** `scs_indirect_optimizer_factory` shape. The check `SCS.is_available(SCS.IndirectSolver) || error(…)` is important — some SCS.jl builds don't ship the indirect solver.
- **BPSDP.** `bpsdp_optimizer_factory` + `progress_monitor_factory` shape. The monitor's job is to push per-outer `inner` iteration counts into a `Vector{Int}` we then dump in `result.json`.

**Do not** port `install_cosmo_quiet_processconstraints!` unless we actually observe the noisy failure mode. Leave it commented out and document the trigger.

Also port:
- `set_solver!(model, solver, options, cg_iterations)` dispatcher — keep the same three-arm shape; this is the function the driver calls per solver per case.
- The result-mining helpers: `bpsdp_state_summary`, `scs_state_summary`, `cosmo_state_summary`, `cosmo_kkt_summary`. These dig the relevant per-solver diagnostics out of the raw optimizer for `result.json`.

## Other helpers worth porting

- `common_runtime_json` (julia/blas threads, hostname, max RSS, julia version).
- `json_ready`, `write_json` — the recursive sanitizer for non-JSON-safe types.
- `validate_output_dir!` — small but it prevents a typo from nuking `$HOME`.
- `write_markdown_summary` — keep the shape; reduce the columns to D7.

## Cut entirely

- Phase-1 `.dat-c` roundtrip (`probes/solve_datc_bpsdp_primitive.jl`, `probes/h4_nk_periodic_moment_sos_export_libsdp.jl`).
- Phase-2 `AA*` diagnostic (`probes/h2_nk2_aastar_diagnostic.jl`).
- Sector-pruning + Newton-chip experiments.
- The thirty-something solver-tuning CLI flags. Per D5 we keep only the narrow set listed in `plan/script-layout.md` (a few iter / tol knobs per solver, plus `--time-limit`).
- All Mosek-only / GPU paths.

## Sizing target (test-stage)

- `driver.jl`: ≤ ~400 lines (three solvers, three result-mining branches, JSON + md writers).
- `h2_launcher.jl` / `h4_launcher.jl`: ≤ ~30 lines each.
- `sweep.sh`: ≤ ~80 lines.
- `setup_pyscf_env.sh`: ≤ ~20 lines.
- `setup_julia_env.sh`: ≤ ~20 lines.
- `extract_h2_nk_integrals.py` / `extract_h4_nk_integrals.py`: ~220 lines each.
- `README.md`: ≤ ~60 lines.

If the diff blows past these, stop and re-justify before adding more.

## Eventual docs port

When the test-stage script is verified, the Literate page (`docs/src/examples/literate/v2rdm_periodic_solver_evidence.jl`) should be **smaller still** — only the H₂/Nk=2 case live, larger results embedded as a static markdown table — not a copy of the full driver. The driver stays in `test/`.
