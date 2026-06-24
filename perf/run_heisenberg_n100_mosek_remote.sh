#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

julia --startup-file=no -e '
using Pkg
Pkg.activate(".mosek-env")
Pkg.develop(path=pwd())
Pkg.add(["JuMP", "MosekTools"])
Pkg.instantiate()
Pkg.precompile()
'

export NCTS_HEISENBERG_N="${NCTS_HEISENBERG_N:-100}"
export NCTS_RESULTS_DIR="${NCTS_RESULTS_DIR:-results}"
export NCTS_OFFBLOCK_CHECK="${NCTS_OFFBLOCK_CHECK:-off}"
export NCTS_CHECK_INVARIANCE="${NCTS_CHECK_INVARIANCE:-true}"
export NCTS_MOSEK_LOG="${NCTS_MOSEK_LOG:-0}"

log_path="results/heisenberg_n${NCTS_HEISENBERG_N}_order2_pauli_symmetry_mosek.log"
julia --startup-file=no --project=.mosek-env --threads=auto \
    perf/heisenberg_n100_mosek.jl 2>&1 | tee "$log_path"
