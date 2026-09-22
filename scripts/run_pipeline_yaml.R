#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline_yaml.R
# Reference-independent CNV karyotyping pipeline — YAML-config-driven variant
#
# Mirrors scripts/run_pipeline.R's logic, but reads ALL parameters from a
# single YAML config file instead of ~40 individual CLI flags. This is a
# separate, standalone script — it does not touch or replace run_pipeline.R.
#
# Usage:
#   Rscript run_pipeline_yaml.R --config configs/run_config_TE.yaml
#   Rscript run_pipeline_yaml.R --config configs/run_config_TE.yaml \
#       --override '{min_overlap: 0.65, range: 0.1}'
#   Rscript run_pipeline_yaml.R --help
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# =============================================================================
# Resolve project root from script location
# =============================================================================

args        <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
project_dir <- dirname(dirname(script_path))

# =============================================================================
# Define options — just --config and --override
# =============================================================================

option_list <- list(
  optparse::make_option(
    "--config",
    type    = "character",
    default = NULL,
    help    = "Path to a YAML config file [required]"
  ),
  optparse::make_option(
    "--override",
    type    = "character",
    default = NULL,
    help    = "YAML fragment to override specific --config values, e.g. '{min_overlap: 0.65}' [default: %default]"
  )
)

opt <- optparse::parse_args(
  optparse::OptionParser(
    option_list = option_list,
    description = paste(
      "Reference-independent single-cell CNV karyotyping pipeline (YAML-config variant).",
      "All parameters come from --config; use --override to tweak individual values."
    )
  )
)

if (is.null(opt$config)) {
  stop("--config is required. Usage: Rscript run_pipeline_yaml.R --config path/to/config.yaml")
}
if (!file.exists(opt$config)) {
  stop("--config file not found: ", opt$config)
}

message("Loading config from: ", opt$config)
config <- yaml::read_yaml(opt$config)

if (!is.null(opt$override)) {
  message("Applying overrides: ", opt$override)
  overrides <- yaml::yaml.load(opt$override)
  config    <- utils::modifyList(config, overrides)
}

# `%||%` — fall back to a default if a config key is NULL/absent
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# Validate execution_mode / start_from
# =============================================================================

execution_mode <- config$execution_mode %||% "single"
start_from     <- config$start_from     %||% "block2"

if (!execution_mode %in% c("single", "array")) {
  stop("execution_mode must be 'single' or 'array'. Got: '", execution_mode, "'")
}

if (is.null(config$workdir)) {
  stop("workdir is required in the config.")
}

if (start_from == "block1") {
  if (execution_mode == "array") {
    stop(
      "start_from: block1 not supported with execution_mode: array.\n",
      "For array mode Block1 use: Rscript scripts/run_block1_array.R"
    )
  }
  if (is.null(config$counts_path))     stop("counts_path required for block1")
  if (is.null(config$metadata_path))   stop("metadata_path required for block1")
  if (is.null(config$gene_order_file)) stop("gene_order_file required for block1")
  if (is.null(config$tool_outdir))     stop("tool_outdir required for block1")
}

if (is.null(config$k_value)) {
  stop(
    "k_value is required.\n",
    "Derive it from a trusted calibration point:\n",
    "  k = known_good_threshold / sqrt(known_good_n_cells)"
  )
}

if (start_from == "block2" && is.null(config$tool_outdir)) {
  stop("tool_outdir required when start_from = block2")
}

# =============================================================================
# Load pipeline functions
# =============================================================================

message("Loading pipeline functions...")
source(file.path(project_dir, "R", "cnv_annotation.R"))
source(file.path(project_dir, "R", "cnv_processing.R"))
source(file.path(project_dir, "R", "cnv_scoring.R"))
source(file.path(project_dir, "R", "infercnv.R"))
source(file.path(project_dir, "R", "pipeline.R"))

# =============================================================================
# Load inputs
# =============================================================================

metadata <- NULL
if (!is.null(config$metadata_path)) {
  if (!file.exists(config$metadata_path)) stop("metadata_path not found: ", config$metadata_path)
  message("Loading metadata from: ", config$metadata_path)
  metadata <- readRDS(config$metadata_path)
}

counts_mx <- NULL
if (!is.null(config$counts_path)) {
  if (!file.exists(config$counts_path)) stop("counts_path not found: ", config$counts_path)
  message("Loading counts matrix from: ", config$counts_path)
  counts_mx <- readRDS(config$counts_path)
}

if (start_from %in% c("block3", "block4")) {

  cnv_annotated <- NULL
  if (!is.null(config$cnv_annotated)) {
    if (!file.exists(config$cnv_annotated)) stop("cnv_annotated not found: ", config$cnv_annotated)
    cnv_annotated <- readRDS(config$cnv_annotated)
    message("Loaded cnv_annotated: ", nrow(cnv_annotated), " rows")
  }

  cell_sizes <- NULL
  if (!is.null(config$cell_sizes)) {
    if (!file.exists(config$cell_sizes)) stop("cell_sizes not found: ", config$cell_sizes)
    cell_sizes <- readRDS(config$cell_sizes)
    message("Loaded cell_sizes: ", nrow(cell_sizes), " rows")
  }

} else {
  cnv_annotated <- NULL
  cell_sizes    <- NULL
}

chromosome_arms <- NULL
if (!is.null(config$chromosome_arms_path)) {
  chromosome_arms <- readRDS(config$chromosome_arms_path)
} else {
  chrom_path <- file.path(project_dir, "data", "hg38_chromosome_arms.rds")
  if (!file.exists(chrom_path)) {
    stop(
      "No chromosome_arms_path given and built-in not found at:\n",
      chrom_path, "\n",
      "Either set chromosome_arms_path in the config or save hg38_chromosome_arms.rds to data/"
    )
  }
  message("Loading built-in hg38 chromosome arms from: ", chrom_path)
  chromosome_arms <- readRDS(chrom_path)
}

if (is.null(config$coding_genes_path)) {
  stop("coding_genes_path is required in the config.")
}
if (!file.exists(config$coding_genes_path)) {
  stop("coding_genes_path not found: ", config$coding_genes_path)
}
message("Loading coding genes from: ", config$coding_genes_path)
coding_genes <- readRDS(config$coding_genes_path)

coding_gr <- GenomicRanges::GRanges(
  seqnames = coding_genes$chr,
  ranges   = IRanges::IRanges(start = coding_genes$start, end = coding_genes$end)
)
coding_gr$gene       <- coding_genes$gene_name
coding_expressed_set <- unique(coding_genes$gene_name)
message("Coding genes loaded: ", length(coding_expressed_set))

# =============================================================================
# Create output directory
# =============================================================================

if (!dir.exists(config$workdir)) {
  message("Creating output directory: ", config$workdir)
  dir.create(config$workdir, recursive = TRUE)
}

# =============================================================================
# Print submission summary
# =============================================================================

message(paste0(
  "\n=== Submission parameters ===\n",
  "  EXECUTION_MODE:          ", execution_mode, "\n",
  "  START_FROM:              ", start_from, "\n",
  "  WORKDIR:                 ", config$workdir, "\n",
  "  TOOL:                    ", config$tool %||% "infercnv", "\n",
  "  N_SPLITS:                ", config$n_splits_within %||% 3, "\n",
  "  CUTOFF:                  ", config$cutoff %||% 0.1, "\n",
  "  REMOVE_REFERENCE:        ", config$remove_reference %||% TRUE, "\n",
  "  OVERLAP_METHOD:          ", config$overlap_method %||% "reciprocal", "\n",
  "  MIN_OVERLAP:             ", config$min_overlap %||% 0.8, "\n",
  "  RANGE:                   ", config$range %||% 0.15, "\n",
  "  MAX_MB:                  ", config$max_mb %||% 120, "\n",
  "  SENSITIVITY_FLOOR_MB:    ", config$sensitivity_floor_mb %||% 20, "\n",
  "  K_INTERVAL:              ", config$k_interval %||% 60, "\n",
  "  K_VALUE:                 ", config$k_value, "\n",
  "  GROUP_COLS:              ", paste(config$group_cols %||% "cell_type", collapse = ", "), "\n",
  "  BY_COL:                  ", paste(config$by_col %||% config$group_cols %||% "cell_type", collapse = ", "), "\n",
  "============================="
))

# =============================================================================
# Run pipeline
# =============================================================================

results <- run_full_cnv_pipeline(

  # ---- General --------------------------------------------------------------
  start_from        = start_from,
  execution_mode     = execution_mode,
  save_intermediate = isTRUE(config$save_intermediate),
  workdir           = config$workdir,

  # ---- Block 1 ----------------------------------------------------------------
  counts_mx              = counts_mx,
  metadata               = metadata,
  group_clusters_col     = config$group_clusters_col %||% "cell_type",
  gene_order_file        = config$gene_order_file,
  chromosomes_to_exclude = config$chr_exclude %||% c("MT", "Y"),
  n_splits_within        = config$n_splits_within %||% 3L,
  tool_outdir            = config$tool_outdir,
  tool                   = config$tool %||% "infercnv",
  cutoff                 = config$cutoff %||% 0.1,
  remove_ref             = isTRUE(config$remove_reference %||% TRUE),
  resume_if_exists       = isTRUE(config$resume_if_exists %||% TRUE),
  clonal_col             = config$clonal_col,
  donor_col              = config$donor_col,

  # ---- Block 2 ----------------------------------------------------------------
  pattern                      = config$pattern %||% "^run\\.final",
  k_discrete                   = config$k_interval %||% 60,
  min_overlap_consistent_calls = config$min_overlap_consistent_calls %||% 0.75,
  min_overlap_multiple_nodes   = config$min_overlap_multiple_nodes %||% 0.6,
  filter_seq_mb_init           = config$filter_seq_mb_init %||% 5.0,
  filter_seq_mb_equiv          = config$filter_seq_mb_equiv %||% 7.0,
  min_references                = config$min_references %||% 2L,
  parallel                     = isTRUE(config$parallel),
  cores                        = config$cores %||% 1L,
  pct_max                   = config$pct_max %||% 45,
  pct_floor                 = config$pct_floor %||% 30,
  min_expr_density          = config$min_expr_density %||% 1.5,
  min_coding_density        = config$min_coding_density %||% 1.0,
  max_gap_mb                = config$max_gap_mb %||% 10,
  coding_gr                 = coding_gr,
  coding_expressed_set      = coding_expressed_set,

  # ---- Adaptive overlap ---------------------------------------------------
  range  = config$range  %||% 0.15,
  max_mb = config$max_mb %||% 120,

  # ---- Block 3 -----------------------------------------------------------------
  chromosome_arms = chromosome_arms,
  group_cols      = config$group_cols %||% "cell_type",

  # ---- Block 4 -----------------------------------------------------------------
  cell_sizes            = cell_sizes,
  cnv_annotated         = cnv_annotated,
  by                    = config$by_col %||% config$group_cols,
  sample_col            = config$sample_col %||% "cell_type",
  cell_col              = config$cell_id_col %||% "cell_name",
  min_required_cells    = config$min_required_cells %||% 5L,
  p_arm_permission      = config$p_arm_permission %||% 70,
  q_arm_permission      = config$q_arm_permission %||% 70,
  whole_chr_permission  = config$whole_chr_permission %||% 60,
  min_overlap           = config$min_overlap %||% 0.8,
  sensitivity_floor_mb  = config$sensitivity_floor_mb %||% 20,
  k_threshold_growth    = config$k_value,
  overlap_method        = config$overlap_method %||% "reciprocal"
)

# =============================================================================
# Save outputs
# =============================================================================

out_path <- file.path(config$workdir, "pipeline_results.rds")
message("Saving full pipeline results to: ", out_path)
saveRDS(results, out_path)

message(sprintf("\nDone. Results saved to: %s", out_path))
