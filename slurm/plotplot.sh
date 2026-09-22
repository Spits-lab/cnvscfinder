BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

sbatch \
  --time=01:00:00 \
  --mem=128G \
  --job-name="plot" \
  --output="${BASE}/logs/plotplot_%j.out" \
  --wrap="
    module purge
    module load inferCNV/1.18.1-foss-2023a
    module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
    module load R-bundle-CRAN/2023.12-foss-2023a
    Rscript ${BASE}/R/ploting_ploting.R
  "
