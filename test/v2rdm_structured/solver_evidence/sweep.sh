#!/usr/bin/env bash
set -euo pipefail

H2_NKS="${H2_NKS-2 3 4}"
H4_NKS="${H4_NKS-2 3}"
OUTROOT="${OUTROOT:-output/solver-evidence}"
INTROOT="${INTROOT:-integrals}"
JULIA="${JULIA:-julia}"
JULIA_PROJECT_PATH="${SOLVER_EVIDENCE_JULIA_PROJECT:-$HOME/.julia-envs/nctssos-solver-evidence}"
JULIA_THREADS="${JULIA_THREADS:-8}"
SOLVER_THREADS="${SOLVER_THREADS:-${JULIA_THREADS}}"
SOLVE_SECONDS="${SOLVE_SECONDS:-3600}"
PYSCF_PYTHON="${PYSCF_PYTHON:-$HOME/.venvs/pyscf-h/bin/python}"
SOLVERS="${SOLVERS:-cosmo,scs,bpsdp}"
FORCE="${FORCE:-0}"

force_arg=()
if [[ "$FORCE" == "1" || "$FORCE" == "true" || "$FORCE" == "yes" ]]; then
  force_arg=(--force)
fi

run_case() {
  local system=$1 nk=$2
  local intdir="$INTROOT/${system}_chain_nk${nk}"
  local txt="$intdir/${system}_chain_nk${nk}_integrals_drop1e-12.txt"
  local meta="$intdir/${system}_chain_nk${nk}_meta.json"

  if [[ ! -f "$txt" || ! -f "$meta" ]]; then
    mkdir -p "$intdir"
    "$PYSCF_PYTHON" "test/v2rdm_structured/solver_evidence/extract_${system}_nk_integrals.py" \
      --nk="$nk" --outdir="$intdir" --drop-tol=1e-12
  fi

  JULIA_NUM_THREADS="$JULIA_THREADS" \
    OPENBLAS_NUM_THREADS="$SOLVER_THREADS" \
    OMP_NUM_THREADS="$SOLVER_THREADS" \
    MKL_NUM_THREADS="$SOLVER_THREADS" \
    BLIS_NUM_THREADS="$SOLVER_THREADS" \
    SOLVER_THREADS="$SOLVER_THREADS" \
    "$JULIA" --startup-file=no --project="$JULIA_PROJECT_PATH" \
    "test/v2rdm_structured/solver_evidence/${system}_launcher.jl" \
    --nk="$nk" --integrals="$txt" --meta="$meta" \
    --solvers="$SOLVERS" --time-limit="$SOLVE_SECONDS" \
    --solver-threads="$SOLVER_THREADS" \
    --output-dir="$OUTROOT/$system/nk$nk" "${force_arg[@]}"
}

for nk in $H2_NKS; do run_case h2 "$nk"; done
for nk in $H4_NKS; do run_case h4 "$nk"; done
