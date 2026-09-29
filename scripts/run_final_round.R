#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_final_round.R
# One more purification round, using ONE chosen grid-search combo's aneuploid
# calls to build the next reference, then reruns inferCNV and the final CNV
# calling with THAT combo's own parameters (not the grid's defaults). Reuses
# the existing purification-loop machinery end to end - nothing new:
#   next_reference_ids() / make_annotation()  (R/iteration.R)
#   the same inferCNV rerun pattern             (scripts/run_infercnv_solo.R)
#   call_iteration_cnvs()                       (R/iteration.R) - the full
#     discretize -> collapse -> gene-density filter -> merge -> arm-annotate
#     -> locus-cluster -> threshold-score -> arm-filter -> post-merge chain,
#     same as everywhere else in the pipeline
#   infercnv::plot_cnv()                        (as in scripts/plot_infercnv_heatmap.R)
#   the karyotype/frequency heatmap chain        (as in scripts/plot_grid_search_results.R)
#
# Runs the whole thing TWICE - once without inferCNV's own denoise step, once
# with it (denoise = TRUE, standard sd_amplifier = 1.5) - both fully scored
# and plotted, for a side-by-side check. Not opt-in: unlike
# plot_infercnv_heatmap.R this script already reruns inferCNV once no matter
# what, so there's no cheap default path to protect.
#
# Reuses the grid-search config as-is (no new YAML): --config is
# configs/grid_search_nusa.yaml (or similar); the iteration's OWN saved
# config.yaml (next to iter_1/iter_2/...) supplies counts_path/gene_order/
# cutoff/min_ref, same trick as WITH_DENOISE in plot_infercnv_heatmap.R.
#
# Does NOT touch iter_1/2/3 - everything goes to a new <iter_root>/final_round/
# folder, since the new reference is derived from iter_3's own ref_ids
# (overwriting it while reading from it would be a hazard, and would erase
# the state.rds history the rest of the loop relies on).
#
# No chr1q callout in the karyotype plot title: that's a VUB04-specific
# ground-truth check, doesn't apply here.
#
# Usage: Rscript run_final_round.R --config cfg.yaml --combo <combo_id> [--override '{...}']
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(infercnv)
  library(optparse)
  library(dplyr)
})
for (f in c("cnv_processing.R", "cnv_scoring.R", "cnv_annotation.R",
            "pipeline.R", "infercnv.R", "ploting_functions.R", "config.R", "iteration.R")) {
  source(file.path(BASE, "R", f))
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL),
  make_option("--combo",    type = "character", help = "combo_id from grid_search_summary.csv [required]")
)))
if (is.null(opt$combo)) stop("--combo is required (a combo_id from grid_search_summary.csv)")

cfg <- load_config(opt$config, opt$override)
overlap_ok <- function(m) !inherits(try(compute_overlap(1L, 2L, 1L, 2L, method = m), silent = TRUE), "try-error")
assert_valid_grid_search_config(cfg, overlap_ok = overlap_ok)

resolved <- resolve_iter_infercnv_obj(cfg$iter_path)
if (!is.null(resolved$error)) stop("iter_path: ", resolved$error)
iter_root <- dirname(resolved$iter_dir)
cat("Base iteration: ", resolved$iter_dir, "\n", sep = "")

# ── The chosen combo's row (parameters) and cells it flagged ────────────────
summary_path <- file.path(cfg$out_dir, "grid_search_summary.csv")
if (!file.exists(summary_path)) stop("grid_search_summary.csv not found: ", summary_path)
summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)
combo_row  <- summary_df[summary_df$combo_id == opt$combo, ]
if (nrow(combo_row) != 1L) stop("combo_id '", opt$combo, "' matches ", nrow(combo_row), " row(s) in ", summary_path)

combo_file <- file.path(cfg$out_dir, paste0(opt$combo, "_scored.rds"))
if (!file.exists(combo_file)) stop("scored file for combo not found: ", combo_file)
combo_scored        <- readRDS(combo_file)
combo_aneuploid_ids <- unique(as.character(combo_scored[[cfg$cell_col]]))
cat("Combo '", opt$combo, "': ", length(combo_aneuploid_ids), " cell(s) with a call\n", sep = "")

# ── Next reference: drop cells the chosen combo flagged from iter_3's ref ───
state    <- read_iter_state(resolved$iter_dir)
metadata <- readRDS(cfg$metadata_path)
all_ids  <- as.character(metadata[[cfg$cell_col]])

new_ref_ids <- next_reference_ids(state$ref_ids, combo_aneuploid_ids, all_ids, rule = "remove")
cat("Reference: ", length(state$ref_ids), " -> ", length(new_ref_ids),
    " (", length(state$ref_ids) - length(new_ref_ids), " removed)\n", sep = "")

# ── Iteration's own saved config for the inferCNV plumbing ──────────────────
iter_cfg_path <- file.path(iter_root, "config.yaml")
if (!file.exists(iter_cfg_path)) stop("iteration config not found: ", iter_cfg_path)
iter_cfg <- load_config(iter_cfg_path)
for (k in c("counts_path", "gene_order")) {
  if (is.null(iter_cfg[[k]])) stop(iter_cfg_path, " has no ", k)
}

final_dir <- file.path(iter_root, "final_round")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

counts_mx <- readRDS(iter_cfg$counts_path)
annot     <- make_annotation(colnames(counts_mx), new_ref_ids, min_ref = iter_cfg$min_ref %||% 50L)
cat("Reference cells:", sum(annot$split_group == "refs"),
    "| Observation cells:", sum(annot$split_group == "non_refs"), "\n")

chromosome_arms <- readRDS(cfg$chromosome_arms_path)
coding_genes    <- readRDS(cfg$coding_genes_path)

final_cfg <- modifyList(cfg, list(
  k_discrete           = combo_row$k_dis,
  k_threshold          = combo_row$k_fre,
  sensitivity_floor_mb = combo_row$sens_floor,
  min_overlap          = combo_row$ovlp,
  pct_floor            = combo_row$pct_floor,
  min_expr_density     = combo_row$min_expr_density
))

# ── One pass: rerun inferCNV -> call -> save state -> plot (HMM + karyotype) ──
run_and_score <- function(denoise_flag, out_subdir, label) {
  cat("\n=== Pass: ", label, " (denoise = ", denoise_flag, ") ===\n", sep = "")
  options(scipen = 100)
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix       = counts_mx,
    annotations_file        = annot,
    gene_order_file         = iter_cfg$gene_order,
    chr_exclude             = c("chrMT", "chrY", "chrX"),
    ref_group_names         = "refs",
    min_max_counts_per_cell = c(100, 1e6)
  )

  out_dir <- file.path(final_dir, out_subdir)
  t_start <- proc.time()
  infercnv::run(
    infercnv_obj        = infercnv_obj,
    out_dir             = out_dir,
    cutoff              = iter_cfg$cutoff %||% 0.1,
    cluster_by_groups   = TRUE,
    HMM                 = TRUE,
    denoise             = denoise_flag,
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
  log_metric(stage = paste0("final_round_infercnv_", label), combo = opt$combo,
             denoise = denoise_flag, ref_size = sum(annot$split_group == "refs"),
             elapsed_s = round((proc.time() - t_start)[["elapsed"]], 1))

  final_obj <- readRDS(file.path(out_dir, "run.final.infercnv_obj"))

  t_call <- proc.time()
  scored <- call_iteration_cnvs(final_obj, metadata, final_cfg, coding_genes, chromosome_arms)
  log_metric(stage = paste0("final_round_calling_", label), combo = opt$combo,
             n_rows = nrow(scored), elapsed_s = round((proc.time() - t_call)[["elapsed"]], 1))

  aneuploid <- unique(as.character(scored[[cfg$cell_col]]))
  write_iter_state(out_dir, list(
    iter = "final", combo = opt$combo, denoise = denoise_flag,
    ref_ids = new_ref_ids, aneuploid_ids = aneuploid, scored_events = scored
  ))
  cat(label, ": ", length(new_ref_ids), " reference cells, ", length(aneuploid),
      " with a call, ", nrow(scored), " event rows\n", sep = "")
  cat("State saved to: ", file.path(out_dir, "state.rds"), "\n", sep = "")

  # Plot 1: standard inferCNV HMM heatmap
  plots_dir <- file.path(out_dir, "plots")
  infercnv::plot_cnv(
    infercnv_obj    = final_obj,
    out_dir         = plots_dir,
    output_filename = paste0("infercnv_hmm_plot_final_", label),
    output_format   = "pdf",
    x.range         = "auto",
    x.center        = 1,
    title           = paste0("inferCNV HMM (final round, ", label, ", combo=", opt$combo, ")"),
    color_safe_pal  = FALSE
  )
  cat("HMM plot saved to: ", plots_dir, "\n", sep = "")

  # Plot 2: karyotype/frequency heatmap (same chain as plot_grid_search_results.R)
  if (is.null(scored) || nrow(scored) == 0) {
    cat("No scored events — skipping karyotype plot\n")
    return(invisible(NULL))
  }
  n_total_cells <- dplyr::n_distinct(scored[[cfg$cell_col]])
  title <- sprintf(
    "final round (%s, combo=%s)\ncells=%d | gain=%d | loss=%d",
    label, opt$combo, n_total_cells,
    sum(scored$cnv_state == "gain", na.rm = TRUE),
    sum(scored$cnv_state == "loss", na.rm = TRUE)
  )

  genome_structure <- prepare_genome_structure(chromosome_arms)
  cnv_mapped <- map_cnv_to_genome(
    scored, genome_structure, threshold = 85, arrange_df_cols = cfg$cell_group_cluster
  )

  if (is.null(cnv_mapped) || nrow(cnv_mapped) == 0) {
    cat("Empty cnv_mapped — skipping karyotype plot\n")
    return(invisible(NULL))
  }

  heatmap_plot <- prepare_cnv_plot(
    cnv_mapped, genome_structure, grouping_cols = cfg$cell_group_cluster,
    state_colors = c("gain" = "#E64B35", "loss" = "royalblue4")
  )
  p <- plot_cnv_karyotype(
    heatmap_plot, genome_structure, ideogram_ratio = 0.05,
    arm_colors = c("p" = "paleturquoise4", "cen" = "black", "q" = "red4"),
    show_legend = TRUE
  ) + ggplot2::ggtitle(title)

  karyotype_pdf <- file.path(final_dir, paste0("final_round_karyotype_", label, ".pdf"))
  pdf(karyotype_pdf, width = 20, height = 12)
  print(p)
  dev.off()
  cat("Karyotype plot saved to: ", karyotype_pdf, "\n", sep = "")
  invisible(NULL)
}

run_and_score(denoise_flag = FALSE, out_subdir = "infercnv",          label = "no_denoise")
run_and_score(denoise_flag = TRUE,  out_subdir = "infercnv_denoised", label = "denoised")

cat("\nDone. See: ", final_dir, "\n", sep = "")
