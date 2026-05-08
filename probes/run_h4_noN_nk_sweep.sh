#!/usr/bin/env bash
# Sequential H4/Nk=2..6 BPSDP sweep without the explicit N-hat - N constraint.
# Intended to be launched remotely via:
#   RUN_ID=... nohup bash probes/run_h4_noN_nk_sweep.sh

set -Eeuo pipefail

REPO=${REPO:-/home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark}
BPSDP_PATH=${BPSDP_PATH:-/home/ubuntu/BPSDP.jl}
RUN_ID=${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}
OUTROOT=${OUTROOT:-${REPO}/output/h4_nk_no_N_bpsdp/${RUN_ID}}
ENV_DIR=${ENV_DIR:-/tmp/h4_noN_bpsdp_env}
THREADS=${THREADS:-8}
NKS=${NKS:-2,3,4,5,6}
MAX_ITER=${MAX_ITER:-2000}
CG_MAX_ITER=${CG_MAX_ITER:-2000}
ATOL=${ATOL:-1e-4}
PROGRESS_STRIDE=${PROGRESS_STRIDE:-1}

mkdir -p "${OUTROOT}"
LOG=${LOG:-${OUTROOT}/driver.log}
STATUS=${STATUS:-${OUTROOT}/status.txt}

exec > >(tee -a "${LOG}") 2>&1

write_status() {
    local state=$1
    {
        printf 'state=%s\n' "${state}"
        printf 'timestamp=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'run_id=%s\n' "${RUN_ID}"
        printf 'repo=%s\n' "${REPO}"
        printf 'outroot=%s\n' "${OUTROOT}"
        printf 'threads=%s\n' "${THREADS}"
        printf 'nks=%s\n' "${NKS}"
        printf 'max_iter=%s\n' "${MAX_ITER}"
        printf 'cg_max_iter=%s\n' "${CG_MAX_ITER}"
        printf 'atol=%s\n' "${ATOL}"
    } > "${STATUS}"
}

on_exit() {
    code=$?
    if [[ ${code} -eq 0 ]]; then
        write_status finished
    else
        write_status "failed_exit_${code}"
    fi
    exit ${code}
}
trap on_exit EXIT

write_status running

export JULIA_NUM_THREADS=${THREADS}
export OPENBLAS_NUM_THREADS=${THREADS}
export OMP_NUM_THREADS=${THREADS}
export MKL_NUM_THREADS=${THREADS}
export REPO
export BPSDP_PATH

printf '== H4 no-N Nk sweep driver ==\n'
printf 'started=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'repo=%s\n' "${REPO}"
printf 'bpsdp=%s\n' "${BPSDP_PATH}"
printf 'outroot=%s\n' "${OUTROOT}"
printf 'threads=%s\n' "${THREADS}"
printf 'nks=%s\n' "${NKS}"
printf 'max_iter=%s cg_max_iter=%s atol=%s progress_stride=%s\n' "${MAX_ITER}" "${CG_MAX_ITER}" "${ATOL}" "${PROGRESS_STRIDE}"

cd "${REPO}"

echo '== Julia environment =='
julia --startup-file=no --project="${ENV_DIR}" -e '
    using Pkg
    Pkg.develop(path=ENV["REPO"])
    Pkg.develop(path=ENV["BPSDP_PATH"])
    Pkg.add(["JSON3", "JuMP"])
    Pkg.instantiate()
    Pkg.precompile()
'

echo '== Sweep =='
julia --startup-file=no --project="${ENV_DIR}" probes/h4_nk_jump_bpsdp_sweep.jl \
    --nks="${NKS}" \
    --norb=8 \
    --nelec-per-cell=4 \
    --blocking=momentum \
    --paper-spin \
    --no-1d \
    --outdir="${OUTROOT}" \
    --max-iter="${MAX_ITER}" \
    --cg-max-iter="${CG_MAX_ITER}" \
    --atol="${ATOL}" \
    --progress-stride="${PROGRESS_STRIDE}"

echo '== Done =='
printf 'finished=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
