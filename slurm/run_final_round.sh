#!/bin/bash
# =============================================================================
# slurm/run_final_round.sh
# Submits one job running scripts/run_final_round.R: builds the next reference
# from a chosen grid-search combo, then reruns the full inferCNV pipeline
# TWICE (no-denoise, then denoise) and scores/plots both. This is the same
# cost as running scripts/run_infercnv_solo.R twice, not a plotting job -
# sized accordingly by default (roughly double a single iteration's sizing,
# since both passes run sequentially in one job - see run_iterations.sh for
# the single-pass reference sizing).
#
# Checks the config (--type grid_search, same as run_grid_search.sh) and that
# the chosen combo actually exists in grid_search_summary.csv before
# submitting anything.
#
# Usage:
#   bash slurm/run_final_round.sh <config.yaml> <combo_id> ['<yaml override>']
#   bash slurm/run_final_round.sh configs/grid_search_nusa.yaml \
#     kdis1.50_kfre1.60_floor20.0_ovlp0.80_pctfl30_mindensity1.5
#   TIME_LIMIT=12:00:00 MEM=192G bash slurm/run_final_round.sh <config.yaml> <combo_id>
# =============================================================================
BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
CONFIG="$(readlink -f "$1" 2>/dev/null)"; COMBO="$2"; OVERRIDE="$3"

if [[ -z "${CONFIG}" || ! -f "${CONFIG}" || -z "${COMBO}" ]]; then
  echo "Usage: bash slurm/run_final_round.sh <config.yaml> <combo_id> ['<yaml override>']"
  exit 1
fi

cd "${BASE}" || exit 1
mkdir -p logs

module purge
module load inferCNV/1.18.1-foss-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load R-bundle-CRAN/2023.12-foss-2023a

OVERRIDE_ARGS=()
[[ -n "${OVERRIDE}" ]] && OVERRIDE_ARGS=(--override "${OVERRIDE}")

# Check the config before anything is submitted
Rscript "${BASE}/scripts/validate_config.R" --config "${CONFIG}" "${OVERRIDE_ARGS[@]}" --type grid_search \
  || { echo "Config check failed; nothing submitted."; exit 1; }

# Cheap login-node check that the combo actually exists, before burning compute
ENV_LINES=$(Rscript "${BASE}/scripts/config_to_env.R" --config "${CONFIG}" "${OVERRIDE_ARGS[@]}") \
  || { echo "config error"; exit 1; }
eval "${ENV_LINES}"
SUMMARY_CSV="${OUT_DIR}/grid_search_summary.csv"
if [[ ! -f "${SUMMARY_CSV}" ]]; then
  echo "ERROR: ${SUMMARY_CSV} not found - run the grid search first"; exit 1
fi
if ! grep -q "\"${COMBO}\"" "${SUMMARY_CSV}"; then
  echo "ERROR: combo_id '${COMBO}' not found in ${SUMMARY_CSV}"; exit 1
fi

TIME_LIMIT="${TIME_LIMIT:-02:00:00}"
MEM="${MEM:-128G}"

JID=$(sbatch --parsable --time="${TIME_LIMIT}" --mem="${MEM}" \
  --job-name="final_round" \
  --output="${BASE}/logs/final_round_%j.out" --error="${BASE}/logs/final_round_%j.err" \
  --wrap="
    module purge
    module load inferCNV/1.18.1-foss-2023a
    module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
    module load R-bundle-CRAN/2023.12-foss-2023a
    Rscript ${BASE}/scripts/run_final_round.R --config '${CONFIG}' --combo '${COMBO}' ${OVERRIDE:+--override '${OVERRIDE}'}
  ")

echo "Submitted: job ${JID}"
echo "Monitor: tail -f ${BASE}/logs/final_round_${JID}.out"
