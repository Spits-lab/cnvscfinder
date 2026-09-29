#!/bin/bash
# =============================================================================
# slurm/run_grid_search.sh
# Submits one job running scripts/run_grid_search.R (successor to
# legacy/slurm/round_3.sh + legacy/scripts/run_block_test.R). Checks the
# config with --type grid_search before submitting, mirroring how
# slurm/run_iterations.sh checks its own config first.
#
# Usage:
#   bash slurm/run_grid_search.sh <config.yaml> ['<yaml override>']
#   bash slurm/run_grid_search.sh configs/grid_search_nusa.yaml
#   bash slurm/run_grid_search.sh configs/grid_search_nusa.yaml '{ovlp_values: [0.6, 0.7]}'
# =============================================================================
BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
CONFIG="$(readlink -f "$1" 2>/dev/null)"; OVERRIDE="$2"

if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
  echo "Usage: bash slurm/run_grid_search.sh <config.yaml> ['<yaml override>']"; exit 1
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

TIME_LIMIT="${TIME_LIMIT:-06:00:00}"
MEM="${MEM:-96G}"

JID=$(sbatch --parsable --time="${TIME_LIMIT}" --mem="${MEM}" \
  --job-name="grid_search" \
  --output="${BASE}/logs/grid_search_%j.out" --error="${BASE}/logs/grid_search_%j.err" \
  --wrap="
    module purge
    module load inferCNV/1.18.1-foss-2023a
    module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
    module load R-bundle-CRAN/2023.12-foss-2023a
    Rscript ${BASE}/scripts/run_grid_search.R --config '${CONFIG}' ${OVERRIDE:+--override '${OVERRIDE}'}
  ")

echo "Submitted: job ${JID}"
echo "Monitor: tail -f ${BASE}/logs/grid_search_${JID}.out"
