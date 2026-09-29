#!/bin/bash
# =============================================================================
# slurm/plot_results.sh
# Submits one job running both new plotting scripts against the same config:
#   1. scripts/plot_infercnv_heatmap.R  - the iteration's own inferCNV heatmap
#   2. scripts/plot_grid_search_results.R - all grid-search combos in one PDF
# Both are fast (~1 min each in practice), so one small job covers both.
#
# Usage:
#   bash slurm/plot_results.sh <config.yaml> ['<yaml override>']
#   WITH_DENOISE=TRUE TIME_LIMIT=08:00:00 MEM=128G bash slurm/plot_results.sh <config.yaml>
#     - WITH_DENOISE is captured here and baked directly into the --wrap
#       string as an explicit prefix on the Rscript call (same as CONFIG/
#       OVERRIDE already are), not left to sbatch's implicit environment
#       export - that was tried first and never actually reached the job.
#       It ALSO reruns the full inferCNV pipeline with denoise, the same cost
#       as scripts/run_infercnv_solo.R - the defaults below are sized for the
#       fast, non-denoise path only, so bump TIME_LIMIT/MEM when using this
#       (see run_iterations.sh for reference sizing).
# =============================================================================
BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
CONFIG="$(readlink -f "$1" 2>/dev/null)"; OVERRIDE="$2"

if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then
  echo "Usage: bash slurm/plot_results.sh <config.yaml> ['<yaml override>']"; exit 1
fi

cd "${BASE}" || exit 1
mkdir -p logs

OVERRIDE_ARGS=""
[[ -n "${OVERRIDE}" ]] && OVERRIDE_ARGS="--override '${OVERRIDE}'"

TIME_LIMIT="${TIME_LIMIT:-01:30:00}"
MEM="${MEM:-32G}"
WITH_DENOISE="${WITH_DENOISE:-FALSE}"

JID=$(sbatch --parsable --time="${TIME_LIMIT}" --mem="${MEM}" \
  --job-name="plot_results" \
  --output="${BASE}/logs/plot_results_%j.out" --error="${BASE}/logs/plot_results_%j.err" \
  --wrap="
    module purge
    module load inferCNV/1.18.1-foss-2023a
    module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
    module load R-bundle-CRAN/2023.12-foss-2023a
    WITH_DENOISE='${WITH_DENOISE}' Rscript ${BASE}/scripts/plot_infercnv_heatmap.R --config '${CONFIG}' ${OVERRIDE_ARGS}
    Rscript ${BASE}/scripts/plot_grid_search_results.R --config '${CONFIG}' ${OVERRIDE_ARGS}
  ")

echo "Submitted: job ${JID}"
echo "Monitor: tail -f ${BASE}/logs/plot_results_${JID}.out"
