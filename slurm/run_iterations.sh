#!/bin/bash
# =============================================================================
# slurm/run_iterations.sh
# Submits chained jobs of the reference-purification loop, one job per
# iteration (slurm/iter_job.sh), from START_ITER through N_ITER (config).
# Iteration n+1 starts only after iteration n finished (afterok); if iteration
# n reaches a plateau, later jobs exit at once.
#
# Usage:
#   bash slurm/run_iterations.sh <config.yaml> ['<yaml override>'] [start_iter]
#   bash slurm/run_iterations.sh configs/iterations_VUB04.yaml
#   bash slurm/run_iterations.sh configs/iterations_VUB04.yaml '{n_iter: 6}' 3
#
# start_iter (default 1): resume from this iteration instead of the beginning.
# It needs iteration (start_iter - 1)'s state.rds already on disk. If iter_<n>
# folders for n >= start_iter already exist (an earlier, now-superseded run),
# you are asked before they are deleted and rebuilt — nothing is removed
# without that confirmation. A non-interactive shell must set CONFIRM=yes.
# =============================================================================
BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
CONFIG="$(readlink -f "$1" 2>/dev/null)"; OVERRIDE="$2"; START_ITER="${3:-1}"

if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
  echo "Usage: bash slurm/run_iterations.sh <config.yaml> ['<yaml override>'] [start_iter]"; exit 1
fi
if ! [[ "${START_ITER}" =~ ^[0-9]+$ ]] || (( START_ITER < 1 )); then
  echo "ERROR: start_iter must be a positive integer, got '${START_ITER}'"; exit 1
fi

cd "${BASE}" || exit 1
mkdir -p logs

module purge
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load R-bundle-CRAN/2023.12-foss-2023a

OVERRIDE_ARGS=()
[[ -n "${OVERRIDE}" ]] && OVERRIDE_ARGS=(--override "${OVERRIDE}")

# Read ITER_ROOT / N_ITER first, before the full check, to know what's on disk
ENV_LINES=$(Rscript "${BASE}/scripts/config_to_env.R" --config "${CONFIG}" "${OVERRIDE_ARGS[@]}") \
  || { echo "config error"; exit 1; }
eval "${ENV_LINES}"

for v in ITER_ROOT N_ITER; do
  if [[ -z "${!v}" || "${!v}" == "NULL" ]]; then echo "ERROR: ${v} missing in config"; exit 1; fi
done
if (( START_ITER > N_ITER )); then
  echo "ERROR: start_iter (${START_ITER}) is after n_iter (${N_ITER})"; exit 1
fi

# ── Iterations that would be overwritten by this run ─────────────────────────
TO_DELETE=()
for i in $(seq "${START_ITER}" "${N_ITER}"); do
  d="${ITER_ROOT}/iter_${i}"
  [[ -d "${d}" ]] && TO_DELETE+=("${d}")
done

if (( ${#TO_DELETE[@]} > 0 )); then
  echo "These iteration folders already exist and will be deleted and rebuilt:"
  printf '  %s\n' "${TO_DELETE[@]}"
  if [[ "${CONFIRM:-}" != "yes" ]]; then
    if [[ -t 0 ]]; then
      read -r -p "Delete and rebuild them? [y/N] " ans
      [[ "${ans}" =~ ^[Yy]$ ]] || { echo "Aborted; nothing changed."; exit 1; }
    else
      echo "Non-interactive shell: rerun with CONFIRM=yes to allow this, or pick a different iter_root."
      exit 1
    fi
  fi
  rm -rf "${TO_DELETE[@]}"
  echo "Deleted."
fi

# ── Config check ───────────────────────────────────────────────────────────
# --existing-run when resuming (iter_root is expected to hold earlier iterations),
# or when iter_root still has leftover files (e.g. a config.yaml) after the cleanup above.
EXISTING_ARGS=()
if (( START_ITER > 1 )); then
  EXISTING_ARGS+=(--existing-run)
elif [[ -d "${ITER_ROOT}" ]] && [[ -n "$(ls -A "${ITER_ROOT}" 2>/dev/null)" ]]; then
  EXISTING_ARGS+=(--existing-run)
fi

Rscript "${BASE}/scripts/validate_config.R" --config "${CONFIG}" "${OVERRIDE_ARGS[@]}" \
  --start-iter "${START_ITER}" "${EXISTING_ARGS[@]}" \
  || { echo "Config check failed; nothing submitted."; exit 1; }

TIME_LIMIT="${TIME_LIMIT:-08:00:00}"
MEM="${MEM:-128G}"

# Exported for the jobs (sbatch passes the environment on)
export BASE CONFIG OVERRIDE ITER_ROOT

# Keep the exact settings next to the results. A resume doesn't overwrite the
# original config.yaml, so iteration 1's settings stay on record.
mkdir -p "${ITER_ROOT}"
if (( START_ITER == 1 )); then
  cp "${CONFIG}" "${ITER_ROOT}/config.yaml"
  [[ -n "${OVERRIDE}" ]] && echo "${OVERRIDE}" > "${ITER_ROOT}/override.yaml"
else
  stamp=$(date +%Y%m%d_%H%M%S)
  cp "${CONFIG}" "${ITER_ROOT}/config_resume_from_iter${START_ITER}_${stamp}.yaml"
  [[ -n "${OVERRIDE}" ]] && echo "${OVERRIDE}" > "${ITER_ROOT}/override_resume_from_iter${START_ITER}_${stamp}.yaml"
fi

prev=""
for i in $(seq "${START_ITER}" "${N_ITER}"); do
  dep=()
  [[ -n "${prev}" ]] && dep=(--dependency=afterok:"${prev}")

  prev=$(sbatch --parsable "${dep[@]}" --time="${TIME_LIMIT}" --mem="${MEM}" \
    --job-name="iter_${i}" \
    --output="${BASE}/logs/iter_${i}_%j.out" --error="${BASE}/logs/iter_${i}_%j.err" \
    --export=ALL,ITER="${i}" "${BASE}/slurm/iter_job.sh")
  echo "iteration ${i}: job ${prev}"
done
echo "Metrics: grep -h '^\[METRIC\]' ${BASE}/logs/iter_*.out"
