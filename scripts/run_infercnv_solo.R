#!/usr/bin/env Rscript
# scripts/run_infercnv_solo.R   (step 2 of an iteration)
# One pooled inferCNV run: the iteration's reference cells (`refs`) are the
# reference, every other cell is an observation. The annotation is built in
# memory from iter_<n>/state.rds. Output: <iter_dir>/infercnv/run.final.infercnv_obj
#
# Usage: Rscript run_infercnv_solo.R --config cfg.yaml [--override '{...}'] --iter N

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(infercnv)
  library(optparse)
})
source(file.path(BASE, "R/iteration.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL),
  make_option("--iter",     type = "integer")
)))
cfg  <- load_config(opt$config, opt$override)
if (is.null(opt$iter)) stop("--iter is required")

iter_dir <- iter_dir_of(cfg, opt$iter)
outdir   <- file.path(iter_dir, "infercnv")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Data + annotation ────────────────────────────────────────────────────────
state     <- read_iter_state(iter_dir)
counts_mx <- readRDS(cfg$counts_path)

annot <- make_annotation(colnames(counts_mx), state$ref_ids, min_ref = cfg$min_ref %||% 50L)
n_ref <- sum(annot$split_group == "refs")
n_obs <- sum(annot$split_group == "non_refs")
cat("Reference cells:", n_ref, "| Observation cells:", n_obs, "\n")

# ── Run inferCNV ─────────────────────────────────────────────────────────────
options(scipen = 100)
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix       = counts_mx,
  annotations_file        = annot,
  gene_order_file         = cfg$gene_order,
  chr_exclude             = c("chrMT", "chrY", "chrX"),
  ref_group_names         = "refs",
  min_max_counts_per_cell = c(100, 1e6)
)

t_start <- proc.time()
infercnv::run(
  infercnv_obj        = infercnv_obj,
  out_dir             = outdir,
  cutoff              = cfg$cutoff %||% 0.1,
  cluster_by_groups   = TRUE,
  HMM                 = TRUE,
  denoise             = FALSE,
  analysis_mode       = "subclusters",
  output_format       = NA,
  no_plot             = TRUE,
  no_prelim_plot      = TRUE,
  window_length       = 140,
  plot_probabilities  = FALSE,
  plot_steps          = FALSE,
  diagnostics         = FALSE,
  inspect_subclusters = FALSE
)

log_metric(
  iter = opt$iter, stage = "infercnv", ref_size = n_ref, obs_size = n_obs,
  elapsed_s = round((proc.time() - t_start)[["elapsed"]], 1)
)
