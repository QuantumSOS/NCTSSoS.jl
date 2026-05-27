#!/usr/bin/env bash
set -euo pipefail

JULIA="${JULIA:-julia}"
BPSDP_PATH="${BPSDP_PATH:-/home/ubuntu/BPSDP.jl}"
ENV_DIR="${SOLVER_EVIDENCE_JULIA_PROJECT:-$HOME/.julia-envs/nctssos-solver-evidence}"
export BPSDP_PATH ENV_DIR

"$JULIA" --startup-file=no -e '
using Pkg
Pkg.activate(ENV["ENV_DIR"])
Pkg.develop(path=pwd())
Pkg.develop(path=ENV["BPSDP_PATH"])
Pkg.add(["JuMP", "COSMO", "SCS", "JSON3", "IterativeSolvers", "LinearMaps"])
Pkg.instantiate()
Pkg.precompile()
'

echo "Julia solver-evidence env ready at $ENV_DIR"
