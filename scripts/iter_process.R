#!/usr/bin/env Rscript
# =============================================================================
# scripts/iter_process.R   (step 3 of an iteration)
# Calls CNVs from the iteration's inferCNV run (call_iteration_cnvs(), the steps
# of R/round2.R), stores the scored events and the aneuploid cells in
# iter_<n>/state.rds and logs the iteration metrics. Writes iter_<n>/PLATEAU
# when no reference cell was called aneuploid.
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})
for (f in c("cnv_processing.R", "cnv_scoring.R", "cnv_annotation.R",
            "pipeline.R", "infercnv.R", "config.R", "iteration.R")) {
  source(file.path(BASE, "R", f))
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL),
  make_option("--iter",     type = "integer")
)))
cfg <- load_config(opt$config, opt$override)
if (is.null(opt$iter)) stop("--iter is required")

t0       <- proc.time()
iter_dir <- iter_dir_of(cfg, opt$iter)
state    <- read_iter_state(iter_dir)

obj             <- readRDS(file.path(iter_dir, "infercnv", "run.final.infercnv_obj"))
metadata        <- readRDS(cfg$metadata_path)
chromosome_arms <- readRDS(cfg$chromosome_arms_path)
coding_genes    <- readRDS(cfg$coding_genes_path)

scored <- call_iteration_cnvs(obj, metadata, cfg, coding_genes, chromosome_arms)

aneuploid   <- unique(as.character(scored[[cfg$cell_col]]))
ref_flagged <- length(intersect(state$ref_ids, aneuploid))

state$aneuploid_ids <- aneuploid
state$scored_events <- scored
write_iter_state(iter_dir, state)

# Crude 1q indicator (cells with a chr1 gain); the exact "1q recovered" rule is
# still to be defined, so this is a placeholder count.
n_chr1_gain <- length(unique(scored[[cfg$cell_col]][scored$chr == "chr1" &
                                                     scored$cnv_state == "gain"]))

log_metric(
  iter = opt$iter, stage = "process",
  ref_size = length(state$ref_ids), ref_called_aneuploid = ref_flagged,
  ref_remaining = length(state$ref_ids) - ref_flagged,
  n_cells_with_event = length(aneuploid), n_cells_chr1_gain = n_chr1_gain,
  n_event_rows = nrow(scored),
  elapsed_s = round((proc.time() - t0)[["elapsed"]], 1)
)

if (ref_flagged == 0L) {
  file.create(file.path(iter_dir, "PLATEAU"))
  cat("Plateau: no reference cell called aneuploid in iteration", opt$iter, "\n")
}
