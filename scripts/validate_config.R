#!/usr/bin/env Rscript
# =============================================================================
# scripts/validate_config.R
# Checks an iteration config (keys, types, ranges, files) before jobs are
# submitted. Exits with status 1 and lists every problem if the config is invalid.
#
# Usage: Rscript scripts/validate_config.R --config cfg.yaml [--override '{...}'] \
#          [--existing-run] [--start-iter N]
#   --existing-run  the iter_root may already hold results (skip the empty-folder check)
#   --start-iter    iteration this run begins at; > 1 skips the start_ref/start_results
#                   check (a resume reads the previous iteration's state instead)
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages(library(optparse))
source(file.path(BASE, "R/iteration.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",       type = "character"),
  make_option("--override",     type = "character", default = NULL),
  make_option("--existing-run", action = "store_true", default = FALSE),
  make_option("--start-iter",   type = "integer", default = 1L)
)))

cfg <- load_config(opt$config, opt$override)

# Overlap methods live in a registry inside compute_overlap(); ask it directly
overlap_ok <- tryCatch({
  suppressPackageStartupMessages(source(file.path(BASE, "R/cnv_processing.R")))
  function(m) !inherits(try(compute_overlap(1L, 2L, 1L, 2L, method = m), silent = TRUE), "try-error")
}, error = function(e) NULL)
if (is.null(overlap_ok)) message("WARNING: could not load compute_overlap(); overlap methods not checked")

res <- validate_config(cfg, fresh_root = !opt$`existing-run`, iter = opt$`start-iter`,
                      overlap_ok = overlap_ok)
for (w in res$warnings) message("WARNING: ", w)
if (length(res$errors) > 0L) {
  message(length(res$errors), " config problem(s):")
  for (i in seq_along(res$errors)) message("  ", i, ". ", res$errors[i])
  quit(status = 1)
}
message("Config OK: ", opt$config)
