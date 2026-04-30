# Current Task Context

- All testing runs on the server via `easy-ssh`; do not run tests locally unless explicitly asked.

## Server / local path map

- Main benchmark repo:
  - local: `/Users/exaclior/QuantumSOS/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`
  - server (`HAI`, `easy-ssh` remote): `/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark`
- `BPSDP.jl`:
  - local synced checkout: `/Users/exaclior/QuantumSOS/BPSDP.jl`
  - server synced checkout: `/home/ubuntu/BPSDP.jl`
  - observed branch/commit on both: `multithread` @ `865000a`
- `libsdp`:
  - local synced checkout: `/Users/exaclior/QuantumSOS/libsdp`
  - server synced checkout: `/home/ubuntu/libsdp`
  - observed branch/commit on both: `main` @ `fcd9e34`
- Replication reference repo:
  - local: `/Users/exaclior/QuantumSOS/replicate`
  - server: `/home/ubuntu/replicate`

## Current objective

- Before touching the large H4 dump, verify the cheap JuMP/BPSDP file-roundtrip recipe below on `HAI`.
- Solve the SDP produced by `demos/h4_periodic_moment_sos.jl`.
- Backend order:
  1. `libsdp` first, using its complex Hermitian BPSDP route.
  2. `BPSDP.jl` second, as the Julia-side/backend comparison path.
- To avoid paying the sparsity/build cost repeatedly, persist the post-sparsity SDP in libsdp-readable SDPA format (`.dat-c` / complex SDPA-sparse), then load that file from `libsdp`.
- Reference implementation for the export/load pattern lives in `QuantumSOS/replicate`:
  - exporter module: `/Users/exaclior/QuantumSOS/replicate/scripts/nctssos_pqg_regen/SDPALibsdpExport.jl`
  - export driver: `/Users/exaclior/QuantumSOS/replicate/scripts/nctssos_pqg_regen/h4_periodic_full_pqg_export_libsdp.jl`
  - libsdp solve driver: `/Users/exaclior/QuantumSOS/replicate/scripts/nctssos_pqg_solve_bpsdp.py`
  - notes: `/Users/exaclior/QuantumSOS/replicate/notes/nctssos_pqg_libsdp_export.md`
- Export shape to reuse: `MomentProblem` / sparse SDP result -> libsdp-complex SDPA-sparse file plus metadata sidecar:
  - `<outdir>/<basename>.dat-c`
  - `<outdir>/<basename>_meta.txt`

## Python / libsdp environment

Use `uv` for all Python environment management on `HAI`; do not use system Python packages, `pip install --user`, or ad-hoc global installs.

Known-good libsdp environment setup:

```bash
easy-ssh run 'cd /home/ubuntu/libsdp
uv venv .venv
. .venv/bin/activate
CMAKE_BUILD_PARALLEL_LEVEL=2 uv pip install numpy /home/ubuntu/libsdp'
```

Important gotcha: after activating `/home/ubuntu/libsdp/.venv`, run Python from the benchmark repo (or `/tmp`), not from `/home/ubuntu/libsdp`; otherwise Python imports the source-tree `libsdp/` package and misses the built `_libsdp` extension in site-packages.

Verified libsdp complex BPSDP smoke test on `HAI` with this uv environment:

```text
problem: min X subject to X = 1, X ⪰ 0, encoded as one 1x1 complex SDPA-sparse block
procedure: minimize
primal objective: +1.000000000000
||Ax - b||: 2.220e-16
dual objective: +1.000000000000
duality gap: 1.110e-16
```

## Cheap JuMP/BPSDP file-roundtrip recipe

Purpose: build a tiny SDP in JuMP, write it with `write_to_file`, read it back with `read_from_file`, and solve the reloaded model with `BPSDP.jl`. Use MOF JSON (`.mof.json`); LP/MPS do not preserve SDP cones. Do not keep a dedicated repo script for this smoke test.

Problem:

```text
min tr(X)
s.t. X ⪰ 0
     X[1,2] = 1
```

Expected result: BPSDP returns `OPTIMAL` with objective approximately `2.0`.

Run from the repo root:

```bash
easy-ssh run 'tmp=$(mktemp -d)
cat > "$tmp/jump_bpsdp_roundtrip.jl" <<\JULIA
using BPSDP
using JuMP

path = joinpath(mktempdir(), "jump_bpsdp_roundtrip.mof.json")

source = Model()
@variable(source, X[1:2, 1:2], PSD)
@constraint(source, X[1, 2] == 1.0)
@objective(source, Min, X[1, 1] + X[2, 2])
write_to_file(source, path)

roundtrip = read_from_file(path)
set_optimizer(
    roundtrip,
    () -> BPSDP.Optimizer(
        max_iter = 5_000,
        cg_max_iter = 100,
        mu_update_frequency = 25,
        penalty_parameter = 0.1,
        cg_convergence = 1e-12,
        dynamic_cg_convergence = false,
        sdp_objective_convergence = 1e-8,
        sdp_error_convergence = 1e-8,
        guess_type = :zero,
        print_level = 0,
    ),
)
set_silent(roundtrip)
optimize!(roundtrip)

term = termination_status(roundtrip)
obj = objective_value(roundtrip)
term == OPTIMAL || error("BPSDP status was $term, expected OPTIMAL")
isapprox(obj, 2.0; atol = 1e-5, rtol = 0) || error("objective was $obj, expected 2.0")
println("BPSDP status: $term")
println("objective: $obj")
JULIA
julia --startup-file=no --project="$tmp" -e "using Pkg; Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\"); Pkg.add(name=\"JuMP\"); Pkg.instantiate()"
julia --startup-file=no --project="$tmp" "$tmp/jump_bpsdp_roundtrip.jl"'
```
