#!/bin/bash
#SBATCH --job-name=overlap_compare
#SBATCH --time=00:20:00
#SBATCH --mem=16G
#SBATCH --output=/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/logs/overlap_compare_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/logs/overlap_compare_%j.err

module purge
module load inferCNV/1.18.1-foss-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load R-bundle-CRAN/2023.12-foss-2023a

BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
cd "${BASE}" || { echo "ERROR: could not cd to ${BASE}"; exit 1; }

mkdir -p logs

echo "=== Starting overlap method comparison ==="
echo "  Working dir: $(pwd)"
echo "  Started:     $(date)"

Rscript scripts/compare_overlap_methods.R
EXIT_CODE=$?

echo "  Finished:    $(date)"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Comparison completed successfully."
else
    echo "Comparison FAILED with exit code ${EXIT_CODE}"
    exit ${EXIT_CODE}
fi
