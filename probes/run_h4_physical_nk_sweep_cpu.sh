#!/usr/bin/env bash
set -u -o pipefail

# Physical H4 Nk sweep:
#   PySCF physical integrals -> NCTSSoS .dat-c export -> BPSDP.jl primitive CPU solve.
# Designed for HAI.  No GPU path exists here.

NKS_CSV=${NKS_CSV:-"2,3,4,5,6,7,8"}
OUTROOT=${OUTROOT:-"output/h4_nk_sweep_physical_cpu_1h"}
PYTHON=${PYTHON:-"/tmp/h4_pyscf_env/bin/python"}
BPSDP_PROJECT=${BPSDP_PROJECT:-"/tmp/h4_bpsdp_jl_env"}
THREADS=${THREADS:-8}
PYSCF_TIMEOUT=${PYSCF_TIMEOUT:-3600}
EXPORT_TIMEOUT=${EXPORT_TIMEOUT:-3600}
SOLVE_SECONDS=${SOLVE_SECONDS:-3600}
SOLVE_MAX_ITER=${SOLVE_MAX_ITER:-100000}
NORB=${NORB:-8}
NELEC_PER_CELL=${NELEC_PER_CELL:-4}
DROP_TOL=${DROP_TOL:-1e-12}

IFS=',' read -r -a NKS <<< "$NKS_CSV"
mkdir -p "$OUTROOT"
SUMMARY="$OUTROOT/summary.md"

value_from_summary() {
    local file=$1
    local key=$2
    awk -F' = ' -v key="$key" '$1 == key {print $2}' "$file" 2>/dev/null | tail -1
}

json_value() {
    local file=$1
    local key=$2
    "$PYTHON" - "$file" "$key" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
print(data.get(sys.argv[2], ""))
PY
}

append_row() {
    local nk=$1 pyscf_code=$2 export_code=$3 solve_code=$4 nkdir=$5
    local solve_summary="$nkdir/solve.summary"
    local meta="$nkdir/integrals/h4_chain_nk${nk}_meta.json"
    local dat="$nkdir/h4_nk${nk}_physical_pqg.dat-c"
    local status stop outer active recovered perr derr seconds ehf shift constraints primal_dim annz
    status=$(value_from_summary "$solve_summary" status)
    stop=$(value_from_summary "$solve_summary" stop_reason)
    outer=$(value_from_summary "$solve_summary" outer_iterations)
    active=$(value_from_summary "$solve_summary" primal_objective)
    recovered=$(value_from_summary "$solve_summary" recovered_total_primal)
    perr=$(value_from_summary "$solve_summary" primal_error)
    derr=$(value_from_summary "$solve_summary" dual_error)
    seconds=$(value_from_summary "$solve_summary" solve_seconds)
    constraints=$(value_from_summary "$solve_summary" constraints)
    primal_dim=$(value_from_summary "$solve_summary" primal_dim)
    annz=$(value_from_summary "$solve_summary" A_nnz)
    ehf=""
    shift=""
    if [[ -f "$meta" ]]; then
        ehf=$(json_value "$meta" ehf)
        shift=$(json_value "$meta" energy_shift)
    fi
    [[ -n "$status" ]] || status="missing"
    [[ -n "$stop" ]] || stop="missing"
    printf '| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' \
        "$nk" "$pyscf_code" "$export_code" "$solve_code" "$status" "$stop" "$outer" \
        "$active" "$recovered" "$perr" "$derr" "$seconds" "$ehf" "$shift" "$constraints" "$primal_dim" "$dat" >> "$SUMMARY"
}

{
    echo '# Physical H4 Nk sweep: PySCF -> NCTSSoS dat-c -> BPSDP CPU'
    echo
    echo "- started: $(date -Is)"
    echo "- nks: $NKS_CSV"
    echo "- active spatial orbitals per k: $NORB"
    echo "- active electrons per cell: $NELEC_PER_CELL"
    echo "- PySCF: gth-tzvp, pseudo=none, dimension=1, low_dim_ft_type=inf_vacuum, density_fit=gdf"
    echo "- BPSDP backend: cpu; JULIA_NUM_THREADS=$THREADS"
    echo "- solve max seconds per Nk: $SOLVE_SECONDS"
    echo "- solve max iterations per Nk: $SOLVE_MAX_ITER"
    echo
    echo '| Nk | PySCF exit | export exit | solve exit | status | stop | outer | active objective | recovered total Ha/cell | primal err | dual err | solve s | HF Ha/cell | shift | constraints | primal dim | dat-c |'
    echo '|---:|---:|---:|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|'
} > "$SUMMARY"

for nk in "${NKS[@]}"; do
    nk=$(echo "$nk" | xargs)
    [[ -n "$nk" ]] || continue
    nkdir="$OUTROOT/nk$nk"
    intdir="$nkdir/integrals"
    mkdir -p "$intdir"

    txt="$intdir/h4_chain_nk${nk}_integrals_drop1e-12.txt"
    meta="$intdir/h4_chain_nk${nk}_meta.json"
    dat="$nkdir/h4_nk${nk}_physical_pqg.dat-c"
    basename="h4_nk${nk}_physical_pqg"

    pyscf_code=0
    if [[ ! -f "$txt" || ! -f "$meta" ]]; then
        echo "[$(date -Is)] PySCF Nk=$nk"
        timeout "$PYSCF_TIMEOUT" env PYTHONUNBUFFERED=1 "$PYTHON" \
            probes/extract_h4_nk_integrals_pyscf.py \
            --nk="$nk" \
            --n-active-orb="$NORB" \
            --n-active-elec="$NELEC_PER_CELL" \
            --density-fit=gdf \
            --drop-tol="$DROP_TOL" \
            --outdir="$intdir" \
            > "$nkdir/pyscf.log" 2>&1
        pyscf_code=$?
    else
        echo "[$(date -Is)] PySCF exists Nk=$nk; skipping"
    fi

    if [[ "$pyscf_code" -ne 0 || ! -f "$txt" ]]; then
        append_row "$nk" "$pyscf_code" -1 -1 "$nkdir"
        continue
    fi

    export_code=0
    if [[ ! -f "$dat" ]]; then
        echo "[$(date -Is)] export Nk=$nk"
        timeout "$EXPORT_TIMEOUT" julia --startup-file=no --project=. \
            probes/h4_nk_periodic_moment_sos_export_libsdp.jl \
            --nk="$nk" \
            --norb="$NORB" \
            --nelec-per-cell="$NELEC_PER_CELL" \
            --integrals="$txt" \
            --outdir="$nkdir" \
            --basename="$basename" \
            --paper-spin \
            --no-1d \
            > "$nkdir/export.log" 2>&1
        export_code=$?
    else
        echo "[$(date -Is)] export exists Nk=$nk; skipping"
    fi

    if [[ "$export_code" -ne 0 || ! -f "$dat" ]]; then
        append_row "$nk" "$pyscf_code" "$export_code" -1 "$nkdir"
        continue
    fi

    shift=$(json_value "$meta" energy_shift)
    echo "[$(date -Is)] solve Nk=$nk CPU shift=$shift"
    timeout "$((SOLVE_SECONDS + 1800))" env JULIA_NUM_THREADS="$THREADS" julia --startup-file=no --project="$BPSDP_PROJECT" \
        probes/solve_datc_bpsdp_primitive.jl \
        --problem="$dat" \
        --backend=cpu \
        --max-iter="$SOLVE_MAX_ITER" \
        --max-seconds="$SOLVE_SECONDS" \
        --monitor-period=25 \
        --summary="$nkdir/solve.summary" \
        --energy-shift="$shift" \
        > "$nkdir/solve.log" 2>&1
    solve_code=$?

    append_row "$nk" "$pyscf_code" "$export_code" "$solve_code" "$nkdir"
done

{
    echo
    echo "- finished: $(date -Is)"
} >> "$SUMMARY"

echo "summary: $SUMMARY"
