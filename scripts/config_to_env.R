#!/usr/bin/env Rscript
# =============================================================================
# scripts/config_to_env.R
# Turns a YAML config into `export KEY='value'` lines for bash, so the stable
# settings live in one file instead of a 40-variable sbatch command.
# Keys map to the variable names slurm/submit_full.sh reads, upper-cased:
#   min_overlap: 0.75   ->   export MIN_OVERLAP='0.75'
#
# Usage (bash):
#   eval "$(Rscript scripts/config_to_env.R --config configs/iterations_VUB04.yaml \
#           --override '{min_overlap: 0.65, k_dis_value: 1.2}')"
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL,
              help = "YAML fragment overriding config values, e.g. '{min_overlap: 0.65}'")
)))

if (is.null(opt$config) || !file.exists(opt$config)) {
  stop("--config is required and must exist: ", opt$config)
}

config <- yaml::read_yaml(opt$config)
if (!is.null(opt$override)) {
  config <- utils::modifyList(config, yaml::yaml.load(opt$override))
}

fmt <- function(v) {
  if (is.null(v))       return("NULL")
  if (is.logical(v))    return(if (isTRUE(v)) "TRUE" else "FALSE")
  paste(unlist(v), collapse = ",")
}

for (k in names(config)) {
  if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", k)) stop("invalid config key: ", k)
  val <- gsub("'", "'\\\\''", fmt(config[[k]]))   # escape single quotes for bash
  cat(sprintf("export %s='%s'\n", toupper(k), val))
}
