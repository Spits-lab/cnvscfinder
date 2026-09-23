#!/usr/bin/env Rscript
# scripts/run_infercnv_noref.R
# Runs inferCNV on one dataset (RPE or NE) WITHOUT a reference group.
# inferCNV then uses the mean expression of all cells as the baseline, so the
# heatmap shows each cell's deviation from the dataset average.
# Cells are grouped in the heatmap by `cell_identity`.
#
# Usage: Rscript run_infercnv_noref.R <sample> [cutoff] [n_threads]
#   sample:    RPE or NE (reads data/<sample>_counts.rds + data/metadata_<sample>.rds)
#   cutoff:    inferCNV min mean count per gene (default 0.1, 10x-style)
#   n_threads: threads for inferCNV (default 1)

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript run_infercnv_noref.R <sample> [cutoff] [n_threads]")
}
sample    <- args[1]
cutoff    <- if (length(args) >= 2L) as.double(args[2])  else 0.1
n_threads <- if (length(args) >= 3L) as.integer(args[3]) else 1L

if (is.na(cutoff) || cutoff < 0) stop("cutoff must be a number >= 0, got: ", args[2])
if (is.na(n_threads) || n_threads < 1L) stop("n_threads must be >= 1, got: ", args[3])

annotation_col <- "cell_identity"

# ── Paths ─────────────────────────────────────────────────────────────────────
counts_path     <- file.path(BASE, "data", paste0(sample, "_counts.rds"))
metadata_path   <- file.path(BASE, "data", paste0("metadata_", sample, ".rds"))
gene_order_path <- file.path(BASE, "data/hg38_gencode_v27.txt")
outdir          <- file.path(BASE, "infercnv_results", sample, "no_reference")

for (p in c(counts_path, metadata_path, gene_order_path)) {
  if (!file.exists(p)) stop("File not found: ", p)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Sample:", sample, "\n")
cat("Loading counts matrix...\n")
counts_mx <- readRDS(counts_path)
cat("Counts dim:", dim(counts_mx), "\n")

cat("Loading metadata...\n")
metadata <- readRDS(metadata_path)
missing_cols <- setdiff(c("cell_name", annotation_col), colnames(metadata))
if (length(missing_cols) > 0L) {
  stop("Metadata missing columns: ", paste(missing_cols, collapse = ", "))
}
if (anyNA(metadata[[annotation_col]])) {
  stop(sum(is.na(metadata[[annotation_col]])), " cells have NA ", annotation_col)
}

# ── Validate cells match ──────────────────────────────────────────────────────
common_cells <- intersect(colnames(counts_mx), metadata$cell_name)
cat("Common cells:", length(common_cells),
    "(counts:", ncol(counts_mx), "| metadata:", nrow(metadata), ")\n")
if (length(common_cells) == 0L) stop("No cells shared between counts and metadata.")

counts_sub <- counts_mx[, common_cells, drop = FALSE]

# inferCNV annotation: rownames = cell names, one column = group label
annot <- data.frame(
  group     = metadata[[annotation_col]][match(common_cells, metadata$cell_name)],
  row.names = common_cells
)

cat("\nCells per", annotation_col, ":\n")
print(table(annot$group))
saveRDS(annot, file.path(outdir, "annotation_used.rds"))

# ── Create inferCNV object ────────────────────────────────────────────────────
# ref_group_names = NULL → no reference; all cells are observations and the
# baseline is the average of all cells.
cat("\nCreating inferCNV object (no reference)...\n")
options(scipen = 100)
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts_sub,
  annotations_file  = annot,
  gene_order_file   = gene_order_path,
  delim             = "\t",
  ref_group_names   = NULL,
  chr_exclude       = c("chrX", "chrY", "chrM")  # gene order uses chr-prefixed names
)
cat("inferCNV object created ✅\n")

# ── Run inferCNV ──────────────────────────────────────────────────────────────
cat("\nRunning inferCNV...\n")
cat("Output dir:", outdir, "\n")
cat("Cutoff:", cutoff, "| threads:", n_threads, "\n")
t_start <- proc.time()

processed_obj <- infercnv::run(
  infercnv_obj       = infercnv_obj,
  out_dir            = outdir,
  cutoff             = cutoff,
  cluster_by_groups  = TRUE,   # heatmap rows grouped by cell_identity
  denoise            = TRUE,
  HMM                = FALSE,  
  no_prelim_plot     = TRUE,
  no_plot            = TRUE,  
  output_format      = "png",
  plot_steps         = FALSE,
  num_threads        = n_threads,
  resume_mode        = TRUE
)

runtime <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nDone in %.1f min ✅\n", runtime / 60))
cat("Heatmap:", file.path(outdir, "infercnv.png"), "\n")
