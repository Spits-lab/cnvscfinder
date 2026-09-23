#!/bin/bash

BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
SAMPLES="${SAMPLES:-RPE NE}"
CUTOFF="${CUTOFF:-0.1}"
THREADS="${THREADS:-4}"

mkdir -p "${BASE}/logs"

# One job per dataset — each runs inferCNV without a reference
for SAMPLE in ${SAMPLES}; do
  sbatch \
    --time=06:00:00 \
    --mem=96G \
    --cpus-per-task="${THREADS}" \
    --job-name="noref_${SAMPLE}" \
    --output="${BASE}/logs/noref_infercnv_${SAMPLE}_%j.out" \
    --wrap="
      module purge
      module load inferCNV/1.18.1-foss-2023a
      module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
      module load R-bundle-CRAN/2023.12-foss-2023a
      Rscript ${BASE}/scripts/run_infercnv_noref.R ${SAMPLE} ${CUTOFF} ${THREADS}
    "
  echo "Submitted ${SAMPLE} ✅"
done

echo "Monitor: tail -f ${BASE}/logs/noref_infercnv_*.out"
