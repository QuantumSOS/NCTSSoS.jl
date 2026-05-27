#!/usr/bin/env bash
set -euo pipefail

# Sequential BPSDP-only tuning sweep for H2/Nk=2.
#
# Why this exists:
#   plan/findings/residual-norms.md showed COSMO/SCS report infinity-norm residuals
#   while upstream BPSDP reported 2-norm residuals, so cross-solver tolerances were
#   not comparable. BPSDP.jl HEAD now uses infinity-norm residuals natively
#   (see BPSDP.jl/docs/src/decisions/DD-014-infinity-residual-stopping.md), so all
#   three solvers can be run at the same numeric tolerance. We target 1e-5 directly.
#
#   Previous runs of this script used sqrt(m) * 1e-5 ~= 2e-3 to compensate for the
#   2-norm BPSDP convention; that scaling is no longer needed and has been removed.
#
# Grid strategy, deliberately not the 90-case Cartesian product:
#   0. Anchor: run the historical original config once:
#        mu=500, penalty=0.1, dynamic_cg=true, cg=1e-10
#   1. Stage 1: sweep the most suspect interaction, (mu_update_frequency, dynamic_cg),
#        with penalty=0.1 and cg=1e-8:
#        mu in {25,100,500,1000,2000}, dynamic_cg in {true,false}
#   2. Stage 2: select the fastest Stage-1/anchor config meeting the fair BPSDP
#        primal/dual/gap target. Around its (mu,dynamic_cg), sweep:
#        penalty in {0.1,1.0,10.0}, cg in {1e-6,1e-8,1e-10}
#
# This spends runs on knobs with plausible first-order effect and only explores
# penalty/CG after choosing the outer-update behavior. It is not full factorial
# astrology. If interactions matter badly, the ranking will show it and we can pay
# for a wider grid later.
#
# Typical HAI launch from repo root:
#   easy-ssh submit "cd /home/ubuntu/NCTSSoS.jl-main && bash test/v2rdm_structured/solver_evidence/tune_bpsdp.sh"
#
# Resume behavior:
#   Existing OUTROOT/<label>/result.json cases are skipped unless FORCE=1.
#
# Useful overrides:
#   OUTROOT=output/bpsdp-tune FORCE=1 STAGE=stage1 bash test/v2rdm_structured/solver_evidence/tune_bpsdp.sh
#   STAGE=stage2 bash test/v2rdm_structured/solver_evidence/tune_bpsdp.sh
#   STAGE=rank bash test/v2rdm_structured/solver_evidence/tune_bpsdp.sh

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../../.." && pwd)"
cd "$REPO_ROOT"

OUTROOT="${OUTROOT:-output/bpsdp-tune-linf}"
INTROOT="${INTROOT:-integrals}"
JULIA="${JULIA:-julia}"
PYTHON="${PYTHON:-python3}"
PYSCF_PYTHON="${PYSCF_PYTHON:-$HOME/.venvs/pyscf-h/bin/python}"
JULIA_PROJECT_PATH="${SOLVER_EVIDENCE_JULIA_PROJECT:-$HOME/.julia-envs/nctssos-solver-evidence}"
JULIA_THREADS="${JULIA_THREADS:-8}"
SOLVER_THREADS="${SOLVER_THREADS:-$JULIA_THREADS}"
SOLVE_SECONDS="${SOLVE_SECONDS:-3600}"
MAX_ITER="${MAX_ITER:-50000}"
CG_MAX_ITER="${CG_MAX_ITER:-10000}"
TARGET_INF="${TARGET_INF:-1e-5}"
BPSDP_OBJ_TOL_REL="${BPSDP_OBJ_TOL_REL:-0.0}"
BPSDP_ERR_TOL_REL="${BPSDP_ERR_TOL_REL:-0.0}"
STAGE1_MU_VALUES="${STAGE1_MU_VALUES:-25 100 500 1000 2000}"
STAGE2_PENALTY_VALUES="${STAGE2_PENALTY_VALUES:-0.1 1.0 10.0}"
STAGE2_CG_VALUES="${STAGE2_CG_VALUES:-1e-6 1e-8 1e-10}"
STAGE="${STAGE:-all}"   # all | stage1 | stage2 | rank
FORCE="${FORCE:-0}"

mkdir -p "$OUTROOT"
if [[ "${TEE_LOG:-1}" != "0" && -z "${TUNE_BPSDP_LOG_ACTIVE:-}" ]]; then
  export TUNE_BPSDP_LOG_ACTIVE=1
  LOG_PATH="${LOG_PATH:-$OUTROOT/tune_bpsdp_$(date -u +%Y%m%dT%H%M%SZ).log}"
  exec > >(tee -a "$LOG_PATH") 2>&1
fi

is_truthy() {
  case "${1:-}" in
    1|true|TRUE|True|yes|YES|Yes|y|Y) return 0 ;;
    *) return 1 ;;
  esac
}

sanitize_label_part() {
  local s="$1"
  s="${s//./p}"
  s="${s//-/m}"
  s="${s//+/p}"
  printf '%s' "$s"
}

compute_fair_tol() {
  # BPSDP.jl now uses infinity-norm residuals natively. No sqrt(m) scaling.
  # FAIR_TOL is whatever ell-infinity target we want all three solvers to meet.
  # BPSDP_M is kept for the metadata log only.
  if [[ -z "${BPSDP_M:-}" && -n "${M_SOURCE_RESULT:-}" && -f "${M_SOURCE_RESULT:-}" ]]; then
    BPSDP_M="$($PYTHON - "$M_SOURCE_RESULT" <<'PY'
import json, sys
with open(sys.argv[1]) as fh:
    data = json.load(fh)
sz = data.get("sdp_size") or {}
print(sz.get("num_equality_constraints") or sz.get("num_constraints_total") or "")
PY
)"
  fi
  BPSDP_M="${BPSDP_M:-40142}"
  BPSDP_TOL="${BPSDP_TOL:-$TARGET_INF}"
  FAIR_TOL="$BPSDP_TOL"
}

prepare_integrals() {
  INTDIR="$INTROOT/h2_chain_nk2"
  INTEGRALS_TXT="$INTDIR/h2_chain_nk2_integrals_drop1e-12.txt"
  INTEGRALS_META="$INTDIR/h2_chain_nk2_meta.json"
  if [[ -f "$INTEGRALS_TXT" && -f "$INTEGRALS_META" ]]; then
    echo "integrals: using cached $INTEGRALS_TXT"
    return 0
  fi
  echo "integrals: missing cached H2/Nk=2 files; generating with $PYSCF_PYTHON"
  mkdir -p "$INTDIR"
  "$PYSCF_PYTHON" "$SCRIPT_DIR/extract_h2_nk_integrals.py" \
    --nk=2 --outdir="$INTDIR" --drop-tol=1e-12
}

write_case_config() {
  local case_dir="$1" label="$2" stage="$3" mu="$4" pen="$5" dyn="$6" cg="$7"
  cat > "$case_dir/config.json.tmp" <<EOF
{
  "label": "$label",
  "stage": "$stage",
  "system": "h2",
  "nk": 2,
  "mu_update_frequency": $mu,
  "penalty_parameter": $pen,
  "dynamic_cg_convergence": $dyn,
  "cg_convergence": $cg,
  "max_iter": $MAX_ITER,
  "cg_max_iter": $CG_MAX_ITER,
  "target_inf_norm": $TARGET_INF,
  "bpsdp_m_for_scaled_tolerance": $BPSDP_M,
  "bpsdp_err_tol": $FAIR_TOL,
  "bpsdp_obj_tol": $FAIR_TOL,
  "bpsdp_err_tol_rel": $BPSDP_ERR_TOL_REL,
  "bpsdp_obj_tol_rel": $BPSDP_OBJ_TOL_REL,
  "solver_threads": $SOLVER_THREADS,
  "julia_threads": $JULIA_THREADS
}
EOF
  mv "$case_dir/config.json.tmp" "$case_dir/config.json"
}

summarize_result() {
  local result_path="$1" label="$2"
  "$PYTHON" - "$result_path" "$label" "$FAIR_TOL" <<'PY'
import json, math, sys
path, label, tol_s = sys.argv[1], sys.argv[2], sys.argv[3]
tol = float(tol_s)
with open(path) as fh:
    r = json.load(fh)
st = r.get("bpsdp_state") or {}

def num(x):
    if isinstance(x, bool) or x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x) if math.isfinite(float(x)) else None
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def fmt(x):
    v = num(x)
    if v is None:
        return ""
    return f"{v:.6g}"

pr = num(st.get("primal_error"))
du = num(st.get("dual_error"))
op = num(st.get("objective_primal"))
od = num(st.get("objective_dual"))
gap = abs(op - od) if op is not None and od is not None else None
feasible = pr is not None and du is not None and gap is not None and pr <= tol and du <= tol and gap <= tol
status = r.get("termination_status") or st.get("termination_reason") or ""
outer = st.get("outer_iterations", r.get("iterations", ""))
inner = r.get("inner_iterations", "")
wall = r.get("solve_wall_seconds", "")
obj = r.get("active_objective_Ha", "")
print(
    f"case={label} status={status} outer={outer} inner={inner} "
    f"wall_s={fmt(wall)} primal={fmt(pr)} dual={fmt(du)} gap={fmt(gap)} "
    f"objective={fmt(obj)} feasible={int(feasible)}"
)
PY
}

run_case() {
  local stage="$1" mu="$2" pen="$3" dyn="$4" cg="$5"
  local dyn_bit="0"
  [[ "$dyn" == "true" ]] && dyn_bit="1"
  local label="mu${mu}_pen$(sanitize_label_part "$pen")_dyn${dyn_bit}_cg$(sanitize_label_part "$cg")"
  local case_dir="$OUTROOT/$label"
  local result_path="$case_dir/result.json"
  local raw_result_path="$case_dir/bpsdp/result.json"
  mkdir -p "$case_dir"

  if ! is_truthy "$FORCE"; then
    if [[ -f "$result_path" ]]; then
      printf 'skip existing %-34s ' "$label"
      summarize_result "$result_path" "$label"
      return 0
    fi
    if [[ -f "$raw_result_path" ]]; then
      cp "$raw_result_path" "$result_path"
      printf 'skip existing %-34s ' "$label"
      summarize_result "$result_path" "$label"
      return 0
    fi
  fi

  write_case_config "$case_dir" "$label" "$stage" "$mu" "$pen" "$dyn" "$cg"

  echo "== run $label stage=$stage mu=$mu penalty=$pen dynamic_cg=$dyn cg=$cg tol=$FAIR_TOL =="
  local force_arg=()
  if is_truthy "$FORCE"; then
    force_arg=(--force)
  fi

  JULIA_NUM_THREADS="$JULIA_THREADS" \
    OPENBLAS_NUM_THREADS="$SOLVER_THREADS" \
    OMP_NUM_THREADS="$SOLVER_THREADS" \
    MKL_NUM_THREADS="$SOLVER_THREADS" \
    BLIS_NUM_THREADS="$SOLVER_THREADS" \
    SOLVER_THREADS="$SOLVER_THREADS" \
    "$JULIA" --startup-file=no --project="$JULIA_PROJECT_PATH" \
    "$SCRIPT_DIR/h2_launcher.jl" \
    --nk=2 --integrals="$INTEGRALS_TXT" --meta="$INTEGRALS_META" \
    --output-dir="$case_dir" --solvers=bpsdp --time-limit="$SOLVE_SECONDS" \
    --solver-threads="$SOLVER_THREADS" \
    --bpsdp-max-iter="$MAX_ITER" --bpsdp-cg-max-iter="$CG_MAX_ITER" \
    --bpsdp-obj-tol="$FAIR_TOL" --bpsdp-err-tol="$FAIR_TOL" \
    --bpsdp-obj-tol-rel="$BPSDP_OBJ_TOL_REL" --bpsdp-err-tol-rel="$BPSDP_ERR_TOL_REL" \
    --bpsdp-mu-update-frequency="$mu" \
    --bpsdp-penalty-parameter="$pen" \
    --bpsdp-dynamic-cg="$dyn" \
    --bpsdp-cg-tol="$cg" \
    "${force_arg[@]}"

  if [[ ! -f "$raw_result_path" ]]; then
    echo "error: driver finished but did not write $raw_result_path" >&2
    return 1
  fi
  cp "$raw_result_path" "$result_path"
  summarize_result "$result_path" "$label"
}

select_stage1_winner() {
  "$PYTHON" - "$OUTROOT" "$FAIR_TOL" <<'PY'
import glob, json, math, os, sys
root, tol_s = sys.argv[1], sys.argv[2]
tol = float(tol_s)

def num(x):
    if isinstance(x, bool) or x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x) if math.isfinite(float(x)) else None
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

rows = []
for result_path in glob.glob(os.path.join(root, "*", "result.json")):
    case_dir = os.path.dirname(result_path)
    cfg_path = os.path.join(case_dir, "config.json")
    if not os.path.isfile(cfg_path):
        continue
    cfg = load_json(cfg_path)
    if cfg.get("stage") not in ("stage1", "anchor"):
        continue
    r = load_json(result_path)
    st = r.get("bpsdp_state") or {}
    pr = num(st.get("primal_error"))
    du = num(st.get("dual_error"))
    op = num(st.get("objective_primal"))
    od = num(st.get("objective_dual"))
    gap = abs(op - od) if op is not None and od is not None else None
    wall = num(r.get("solve_wall_seconds"))
    vals = [v for v in (pr, du, gap) if v is not None]
    maxerr = max(vals) if vals else math.inf
    feasible = pr is not None and du is not None and gap is not None and pr <= tol and du <= tol and gap <= tol
    rows.append({"cfg": cfg, "feasible": feasible, "wall": wall, "maxerr": maxerr})

if not rows:
    raise SystemExit("no Stage-1/anchor result.json files available; run STAGE=stage1 first")

def key(row):
    wall = row["wall"] if row["wall"] is not None else math.inf
    if row["feasible"]:
        return (0, wall)
    return (1, row["maxerr"], wall)

best = sorted(rows, key=key)[0]
cfg = best["cfg"]
dyn = "true" if cfg.get("dynamic_cg_convergence") else "false"
reason = "feasible_wall" if best["feasible"] else "least_bad_max_residual_or_gap"
print("\t".join([
    str(cfg.get("label")),
    str(cfg.get("mu_update_frequency")),
    dyn,
    reason,
]))
PY
}

rank_results() {
  "$PYTHON" - "$OUTROOT" "$FAIR_TOL" <<'PY'
import glob, json, math, os, sys
root, tol_s = sys.argv[1], sys.argv[2]
tol = float(tol_s)

def num(x):
    if isinstance(x, bool) or x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x) if math.isfinite(float(x)) else None
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def s_num(x):
    v = num(x)
    return "" if v is None else f"{v:.12g}"

def load_json(path):
    with open(path) as fh:
        return json.load(fh)

rows = []
for result_path in glob.glob(os.path.join(root, "*", "result.json")):
    case_dir = os.path.dirname(result_path)
    cfg_path = os.path.join(case_dir, "config.json")
    cfg = load_json(cfg_path) if os.path.isfile(cfg_path) else {}
    r = load_json(result_path)
    st = r.get("bpsdp_state") or {}
    pr = num(st.get("primal_error"))
    du = num(st.get("dual_error"))
    op = num(st.get("objective_primal"))
    od = num(st.get("objective_dual"))
    gap = abs(op - od) if op is not None and od is not None else None
    wall = num(r.get("solve_wall_seconds"))
    vals = [v for v in (pr, du, gap) if v is not None]
    maxerr = max(vals) if vals else math.inf
    feasible = pr is not None and du is not None and gap is not None and pr <= tol and du <= tol and gap <= tol
    rows.append({
        "label": cfg.get("label") or os.path.basename(case_dir),
        "stage": cfg.get("stage", ""),
        "mu": cfg.get("mu_update_frequency", ""),
        "dyn": cfg.get("dynamic_cg_convergence", ""),
        "pen": cfg.get("penalty_parameter", ""),
        "cg": cfg.get("cg_convergence", ""),
        "status": r.get("termination_status") or st.get("termination_reason") or "",
        "outer": st.get("outer_iterations", r.get("iterations", "")),
        "inner": r.get("inner_iterations", ""),
        "wall": wall,
        "pr": pr,
        "du": du,
        "gap": gap,
        "objective": num(r.get("active_objective_Ha")),
        "feasible": feasible,
        "maxerr": maxerr,
    })

def key(row):
    wall = row["wall"] if row["wall"] is not None else math.inf
    if row["feasible"]:
        return (0, wall)
    return (1, row["maxerr"], wall)

rows.sort(key=key)
print("rank\tfeasible\tlabel\tstage\tmu\tdynamic_cg\tpenalty\tcg_tol\tstatus\touter\tinner\twall_s\tprimal_error\tdual_error\tobjective_gap\tactive_objective_Ha")
for i, row in enumerate(rows, 1):
    dyn = row["dyn"]
    if isinstance(dyn, bool):
        dyn = "true" if dyn else "false"
    print("\t".join(map(str, [
        i,
        1 if row["feasible"] else 0,
        row["label"],
        row["stage"],
        row["mu"],
        dyn,
        row["pen"],
        row["cg"],
        row["status"],
        row["outer"],
        row["inner"],
        s_num(row["wall"]),
        s_num(row["pr"]),
        s_num(row["du"]),
        s_num(row["gap"]),
        s_num(row["objective"]),
    ])))
PY
}

run_stage1() {
  echo "== anchor: historical original BPSDP config =="
  run_case anchor 500 0.1 true 1e-10

  echo "== stage1: mu_update_frequency x dynamic_cg at penalty=0.1 cg=1e-8 =="
  local mu dyn
  for mu in $STAGE1_MU_VALUES; do
    for dyn in true false; do
      run_case stage1 "$mu" 0.1 "$dyn" 1e-8
    done
  done
}

run_stage2() {
  local winner_line win_label win_mu win_dyn win_reason
  winner_line="$(select_stage1_winner)"
  IFS=$'\t' read -r win_label win_mu win_dyn win_reason <<< "$winner_line"
  echo "== stage2: winner from stage1/anchor: label=$win_label mu=$win_mu dynamic_cg=$win_dyn reason=$win_reason =="
  echo "== stage2: penalty_parameter x cg_convergence around selected mu/dynamic_cg =="
  local pen cg
  for pen in $STAGE2_PENALTY_VALUES; do
    for cg in $STAGE2_CG_VALUES; do
      run_case stage2 "$win_mu" "$pen" "$win_dyn" "$cg"
    done
  done
}

compute_fair_tol
prepare_integrals

echo "== BPSDP tuning sweep =="
echo "repo=$REPO_ROOT"
echo "outroot=$OUTROOT"
echo "stage=$STAGE force=$FORCE"
echo "stage1_mu_values=$STAGE1_MU_VALUES"
echo "stage2_penalty_values=$STAGE2_PENALTY_VALUES"
echo "stage2_cg_values=$STAGE2_CG_VALUES"
echo "threads: JULIA_THREADS=$JULIA_THREADS SOLVER_THREADS=$SOLVER_THREADS"
echo "fair tolerance (BPSDP now reports ell-infinity): target_inf=$TARGET_INF bpsdp_err_tol=$FAIR_TOL bpsdp_obj_tol=$FAIR_TOL rel_err=$BPSDP_ERR_TOL_REL rel_obj=$BPSDP_OBJ_TOL_REL (m=$BPSDP_M, no sqrt(m) scaling)"
echo "limits: max_iter=$MAX_ITER cg_max_iter=$CG_MAX_ITER nominal_time_limit=$SOLVE_SECONDS"

case "$STAGE" in
  all)
    run_stage1
    run_stage2
    ;;
  stage1)
    run_stage1
    ;;
  stage2)
    run_stage2
    ;;
  rank)
    ;;
  *)
    echo "error: STAGE must be one of: all, stage1, stage2, rank" >&2
    exit 2
    ;;
esac

echo "== ranking: feasible rows meet primal, dual, and objective_gap <= $FAIR_TOL =="
rank_results | tee "$OUTROOT/ranking.tsv"
