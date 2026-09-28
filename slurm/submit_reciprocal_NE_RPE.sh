#!/bin/bash
# Submits the reciprocal-overlap pipeline for NE and RPE (inferCNV run without reference).
BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
SAMPLES="${SAMPLES:-NE RPE}"
GROUP="${GROUP:-cell_identity}"   # annotation column used in run_infercnv_noref.R

for SAMPLE in ${SAMPLES}; do
  sbatch \
    --time=04:00:00 \
    --mem=128G \
    --job-name="cnv_${SAMPLE}_recip" \
    --export=ALL,SCRIPT=${BASE}/scripts/run_pipeline.R,METADATA_PATH="${BASE}/data/metadata_${SAMPLE}.rds",GENE_ORDER="${BASE}/data/hg38_gencode_v27.txt",TOOL_OUTDIR="${BASE}/infercnv_results/${SAMPLE}/no_reference",WORKDIR="${BASE}/cnv_results/${SAMPLE}/round_timer/reciprocal",COUNTS_PATH="${BASE}/data/${SAMPLE}_counts.rds",CHROMOSSOME_PATH="${BASE}/data/hg38_chromosome_arms.rds",GROUP_COLS="${GROUP}",N_SPLITS=3,EXECUTION_MODE=single,TOOL=infercnv,CELL_COL=cell_name,BY_COL="${GROUP}",SAMPLE_COL="${GROUP}",MIN_REQUIRED_CELLS=3,P_ARM_PERMISSION=60,Q_ARM_PERMISSION=60,WHOLE_CHR_PERMISSION=65,MIN_OVERLAP=0.75,MIN_OVERLAP_CONSISTENT=0.75,MIN_OVERLAP_NODES=0.75,MIN_REFERENCES=2,REMOVE_REFERENCE=T,CUTOFF=0.1,START_FROM=block2,K_DIS_VALUE=1.4,K_FRE_VALUE=1.40,SENSITIVITY_FLOOR_MB=20,CELL_GROUP_CLUSTER="${GROUP}",OVERLAP_METHOD=reciprocal,CODING_GENES_PATH="${BASE}/data/hg38_expressing_genes.rds",PCT_MAX=45,PCT_FLOOR=30,MIN_EXPR_DENSITY=1.5,MIN_CODING_DENSITY=1.5,MAX_GAP_MB=10,RANGE=0.05,MAX_MB=100 \
    "${BASE}/slurm/submit_full.sh"
  echo "Submitted ${SAMPLE}"
done
