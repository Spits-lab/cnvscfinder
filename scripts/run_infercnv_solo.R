#!/usr/bin/env Rscript
# scripts/run_round3_infercnv.R
# Runs inferCNV Round 2 with diploid cells as dedicated reference
# No splits — diploid_ref goes directly as reference group

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(infercnv)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
counts_path     <- file.path(BASE, "data/count_mx_VUB04.rds")
metadata_r2_path <- file.path(BASE, "data/metadata0409_round7.rds")
gene_order_path <- file.path(BASE, "data/hg38_gencode_v27.txt")
outdir          <- file.path(BASE, "infercnv_results/VUB04/round7")


dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading counts matrix...\n")
counts_mx <- readRDS(counts_path)
cat("Counts dim:", dim(counts_mx), "\n")

cat("Loading Round 2 metadata...\n")
metadata_r2 <- readRDS(metadata_r2_path)
cat("Metadata rows:", nrow(metadata_r2), "\n")
cat("Groups:\n")
rownames(metadata_r2) <- metadata_r2$cell_name
metadata_r2 <- metadata_r2[,colnames(metadata_r2) %in% "split_group", drop = F]
print(table(metadata_r2$cell_type))

# ── Validate cells match ──────────────────────────────────────────────────────
common_cells <- intersect(
  colnames(counts_mx),
  rownames(metadata_r2)
)
cat("\nCommon cells:", length(common_cells), "\n")

# Subset counts to common cells
counts_sub <- counts_mx[, common_cells, drop = FALSE]


annot <- metadata_r2

cat("\nAnnotation preview:\n")
print(head(annot))
cat("Reference group: diploid_ref\n")
cat("N reference:    ",
    sum(annot$split_group == "ref"), "\n")
cat("N observation:  ",
    sum(annot$split_group != "ref"), "\n")

# ── Create inferCNV object ────────────────────────────────────────────────────
cat("\nCreating inferCNV object...\n")
options(scipen = 100)
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix       = counts_sub,
  annotations_file        = annot,
  gene_order_file         = gene_order_path,
  chr_exclude             = c("MT", "Y", "X"),
  ref_group_names         = "ref",
  min_max_counts_per_cell = c(100, 1e6)
)

cat("inferCNV object created ✅\n")

# ── Run inferCNV ──────────────────────────────────────────────────────────────
cat("\nRunning inferCNV Round 2...\n")
cat("Output dir:", outdir, "\n")
t_start <- proc.time()

processed_obj <- infercnv::run(
  infercnv_obj        = infercnv_obj,
  out_dir             = outdir,
<<<<<<< HEAD
  cutoff              = 0.5,
  cluster_by_groups   = TRUE,
  HMM                 = TRUE,
  denoise             = FALSE,
=======
  cutoff              = 0.1,
  cluster_by_groups   = TRUE,
  HMM                 = TRUE,
  denoise             = T,
>>>>>>> f7a7a33 (feat: initial commit of CNV pipeline scripts)
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

runtime <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nDone in %.1f min ✅\n",
            runtime / 60))
cat("Output saved to:", outdir, "\n")

