#!/usr/bin/env Rscript
# =============================================================================
# scripts/iter_split.R   (step 1 of an iteration)
# Decides which cells are the reference in this iteration and stores them in
# <iter_root>/iter_<n>/state.rds.
#
# Iteration 1 starts from cfg$start_ref (rds character vector of cell ids) or
# cfg$start_results (pipeline_results_*.rds of a previous run; reference = every
# cell without a Block-4 event). Iteration n > 1 applies cfg$ref_rule to the
# state of iteration n-1.
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages(library(optparse))
source(file.path(BASE, "R/iteration.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL),
  make_option("--iter",     type = "integer")
)))

cfg  <- load_config(opt$config, opt$override)
iter <- opt$iter
if (is.null(iter)) stop("--iter is required")
assert_valid_config(cfg, fresh_root = FALSE, iter = iter)

t0       <- proc.time()
iter_dir <- iter_dir_of(cfg, iter)
metadata <- readRDS(cfg$metadata_path)
all_ids  <- as.character(metadata[[cfg$cell_col]])

if (iter == 1L) {
  if (!is.null(cfg$start_ref)) {
    ref_ids <- as.character(readRDS(cfg$start_ref))
  } else if (!is.null(cfg$start_results)) {
    scored  <- readRDS(cfg$start_results)$all_results$block4$scored_events
    ref_ids <- setdiff(all_ids, unique(as.character(scored[[cfg$cell_col]])))
  } else {
    stop("iteration 1 needs start_ref or start_results in the config")
  }
  prev_size <- length(all_ids)
} else {
  prev <- read_iter_state(iter_dir_of(cfg, iter - 1L))
  if (is.null(prev$aneuploid_ids)) stop("iteration ", iter - 1L, " has no aneuploid_ids yet")
  ref_ids   <- next_reference_ids(prev$ref_ids, prev$aneuploid_ids, all_ids,
                                  rule = cfg$ref_rule %||% "remove")
  prev_size <- length(prev$ref_ids)
}

annotation <- make_annotation(all_ids, ref_ids, min_ref = cfg$min_ref %||% 50L)
write_iter_state(iter_dir, list(iter = iter, ref_ids = ref_ids,
                                aneuploid_ids = NULL, scored_events = NULL))

log_metric(
  iter = iter, stage = "split", ref_rule = cfg$ref_rule %||% "remove",
  n_cells = length(all_ids), ref_size = sum(annotation$split_group == "refs"),
  obs_size = sum(annotation$split_group == "non_refs"),
  ref_removed_vs_prev = prev_size - length(ref_ids),
  elapsed_s = round((proc.time() - t0)[["elapsed"]], 1)
)
