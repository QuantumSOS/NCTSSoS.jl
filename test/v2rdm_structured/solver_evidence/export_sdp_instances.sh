#!/usr/bin/env bash
set -euo pipefail

# Generate frozen BPSDP-native SDP payloads for the small periodic V2RDM cases.
# Defaults export H2/Nk=2 and H4/Nk=2 with dependent_rows=:keep.

SYSTEMS="${SYSTEMS:-h2 h4}"
NK="${NK:-2}"
OUTROOT="${OUTROOT:-output/sdp_instances}"
INTROOT="${INTROOT:-integrals}"
JULIA="${JULIA:-julia}"
JULIA_PROJECT_PATH="${SOLVER_EVIDENCE_JULIA_PROJECT:-$HOME/.julia-envs/nctssos-solver-evidence}"
JULIA_THREADS="${JULIA_THREADS:-8}"
SOLVER_THREADS="${SOLVER_THREADS:-${JULIA_THREADS}}"
PYSCF_PYTHON="${PYSCF_PYTHON:-$HOME/.venvs/pyscf-h/bin/python}"
DEPENDENT_ROWS="${DEPENDENT_ROWS:-keep}"
DEPENDENT_ROW_TOL="${DEPENDENT_ROW_TOL:-1e-12}"
FORCE="${FORCE:-0}"
NCTSSOS_GIT_REV="${NCTSSOS_GIT_REV:-$(git rev-parse --short HEAD 2>/dev/null || true)}"
NCTSSOS_GIT_DIRTY="${NCTSSOS_GIT_DIRTY:-$(test -z "$(git status --porcelain 2>/dev/null)" && echo false || echo true)}"
export NCTSSOS_GIT_REV NCTSSOS_GIT_DIRTY

force_arg=()
if [[ "$FORCE" == "1" || "$FORCE" == "true" || "$FORCE" == "yes" ]]; then
  force_arg=(--force)
fi

case_defaults() {
  local system=$1
  case "$system" in
    h2) echo "4 2" ;;
    h4) echo "8 4" ;;
    *) echo "unknown system: $system" >&2; return 2 ;;
  esac
}

run_case() {
  local system=$1 nk=$2
  read -r norb nelec_per_cell < <(case_defaults "$system")

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
    "test/v2rdm_structured/solver_evidence/export_sdp_instance.jl" \
    --system="$system" --nk="$nk" --norb="$norb" --nelec-per-cell="$nelec_per_cell" \
    --integrals="$txt" --meta="$meta" \
    --dependent-rows="$DEPENDENT_ROWS" --dependent-row-tol="$DEPENDENT_ROW_TOL" \
    --solver-threads="$SOLVER_THREADS" \
    --output-dir="$OUTROOT/${system}_nk${nk}" "${force_arg[@]}"
}

for system in $SYSTEMS; do
  run_case "$system" "$NK"
done
