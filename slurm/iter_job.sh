#!/bin/bash
# One iteration of the reference-purification loop, as a single job:
#   split -> inferCNV (pooled reference) -> CNV calling
# Submitted by run_iterations.sh, which exports CONFIG, OVERRIDE (may be empty),
# ITER (iteration number), ITER_ROOT (parent folder of all iterations) and BASE.
module purge
module load inferCNV/1.18.1-foss-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load R-bundle-CRAN/2023.12-foss-2023a

for v in BASE CONFIG ITER ITER_ROOT; do
  if [[ -z "${!v}" ]]; then echo "ERROR: ${v} not set"; exit 1; fi
done

# A plateau in the previous iteration means nothing is left to do
if (( ITER > 1 )) && [[ -f "${ITER_ROOT}/iter_$((ITER-1))/PLATEAU" ]]; then
  echo "Plateau reached in iteration $((ITER-1)); skipping iteration ${ITER}"
  exit 0
fi

ARGS=(--config "${CONFIG}" --iter "${ITER}")
[[ -n "${OVERRIDE}" ]] && ARGS+=(--override "${OVERRIDE}")

t0=$(date +%s)
Rscript "${BASE}/scripts/iter_split.R"        "${ARGS[@]}" || { echo "split failed";   exit 1; }
t1=$(date +%s)
Rscript "${BASE}/scripts/run_infercnv_solo.R" "${ARGS[@]}" || { echo "inferCNV failed"; exit 1; }
t2=$(date +%s)
Rscript "${BASE}/scripts/iter_process.R"      "${ARGS[@]}" || { echo "process failed";  exit 1; }
t3=$(date +%s)

echo "[METRIC] iter=${ITER} stage=job split_s=$((t1-t0)) infercnv_s=$((t2-t1)) process_s=$((t3-t2)) total_s=$((t3-t0))"
