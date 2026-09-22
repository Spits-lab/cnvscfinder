#!/bin/bash

BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

sbatch \
  --time=05:00:00 \
  --mem=128G \
  --job-name="r2_infercnv" \
  --output="${BASE}/logs/pipeline_%j.out" \
  --error="${BASE}/logs/pipeline_%j.err" \
  --wrap="
    module purge
    module load inferCNV/1.18.1-foss-2023a
    module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
    module load R-bundle-CRAN/2023.12-foss-2023a
    Rscript ${BASE}/R/testing.R
  "
