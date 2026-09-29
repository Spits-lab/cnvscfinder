#!/usr/bin/env Rscript
# =============================================================================
# scripts/validate_config.R
# Checks a config (iteration loop OR grid search) before jobs are submitted.
# Exits with status 1 and lists every problem if the config is invalid.
#
# Usage: Rscript scripts/validate_config.R --config cfg.yaml [--override '{...}'] \
#          [--type iteration|grid_search] [--existing-run] [--start-iter N]
#   --type          which validator to run [default: iteration]
#   --existing-run  iteration only: the iter_root may already hold results
#                   (skip the empty-folder check)
#   --start-iter    iteration only: iteration this run begins at; > 1 skips
#                   the start_ref/start_results check (a resume reads the
#                   previous iteration's state instead)
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages(library(optparse))
source(file.path(BASE, "R/config.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",       type = "character"),
  make_option("--override",     type = "character", default = NULL),
  make_option("--type",         type = "character", default = "iteration"),
  make_option("--existing-run", action = "store_true", default = FALSE),
  make_option("--start-iter",   type = "integer", default = 1L)
)))

if (!opt$type %in% c("iteration", "grid_search")) {
  message("ERROR: --type must be 'iteration' or 'grid_search', got '", opt$type, "'")
  quit(status = 1)
}
# --existing-run/--start-iter are iteration-only concepts; catch misuse
# rather than silently ignoring them for a grid-search config.
if (opt$type == "grid_search" && (opt$`existing-run` || opt$`start-iter` != 1L)) {
  message("ERROR: --existing-run/--start-iter only apply to --type iteration")
  quit(status = 1)
}

cfg <- load_config(opt$config, opt$override)

# Overlap methods live in a registry inside compute_overlap(); ask it directly.
# Reused for both types below.
overlap_ok <- tryCatch({
  suppressPackageStartupMessages(source(file.path(BASE, "R/cnv_processing.R")))
  function(m) !inherits(try(compute_overlap(1L, 2L, 1L, 2L, method = m), silent = TRUE), "try-error")
}, error = function(e) NULL)
if (is.null(overlap_ok)) message("WARNING: could not load compute_overlap(); overlap methods not checked")

res <- if (opt$type == "iteration") {
  validate_config(cfg, fresh_root = !opt$`existing-run`, iter = opt$`start-iter`,
                  overlap_ok = overlap_ok)
} else {
  validate_grid_search_config(cfg, overlap_ok = overlap_ok)
}

for (w in res$warnings) message("WARNING: ", w)
if (length(res$errors) > 0L) {
  message(length(res$errors), " config problem(s):")
  for (i in seq_along(res$errors)) message("  ", i, ". ", res$errors[i])
  quit(status = 1)
}
message("Config OK (", opt$type, "): ", opt$config)
