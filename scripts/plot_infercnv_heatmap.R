#!/usr/bin/env Rscript
# =============================================================================
# scripts/plot_infercnv_heatmap.R
# Generalised successor to scripts/plotplotplot.R: draws inferCNV's own
# chromosome x cell heatmap (infercnv::plot_cnv) for one iteration's pooled
# object, resolved from `iter_path` the same way scripts/run_grid_search.R
# does (see R/config.R's resolve_iter_infercnv_obj()) instead of a hand-typed
# path - the exact thing that broke the legacy script when the results
# folder moved. Useful on its own to eyeball an iteration for potential
# losses/coverage gaps before or alongside a grid search on it.
#
# WITH_DENOISE=TRUE (env var, opt-in, off by default): ALSO reruns inferCNV
# with its own standard denoise step (denoise = TRUE, default sd_amplifier =
# 1.5) on the SAME reference split as the resolved iteration, and plots that
# too, for a quick before/after check. This is NOT cheap - it's the full
# inferCNV pipeline again (the same cost as scripts/run_infercnv_solo.R), not
# just a plot - so it needs a much larger --time/--mem than a plain plot.
# Counts/gene order/cutoff/min_ref are read from the resolved iteration's own
# saved config.yaml (written by slurm/run_iterations.sh next to
# iter_1/iter_2/...), not from configs/grid_search_nusa.yaml - no new config
# keys needed there. Output goes to <iter_dir>/infercnv_denoised_check/ - the
# real <iter_dir>/infercnv/ (used by everything else) is never touched.
#
# Usage:
#   Rscript plot_infercnv_heatmap.R --config cfg.yaml [--override '{...}']
#   WITH_DENOISE=TRUE Rscript plot_infercnv_heatmap.R --config cfg.yaml
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(infercnv)
  library(optparse)
})
source(file.path(BASE, "R/config.R"))
source(file.path(BASE, "R/iteration.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL)
)))
cfg <- load_config(opt$config, opt$override)

# base R's as.logical() only recognises "TRUE"/"T"/"true"/... - accept "1" too,
# a common shell-boolean habit that would otherwise silently parse to NA/FALSE.
with_denoise_raw <- Sys.getenv("WITH_DENOISE", "<unset>")
with_denoise <- with_denoise_raw == "1" || isTRUE(as.logical(with_denoise_raw))
cat("WITH_DENOISE env var: '", with_denoise_raw, "' -> with_denoise=", with_denoise, "\n", sep = "")

if (is.null(cfg$iter_path)) stop("config needs iter_path (see scripts/run_grid_search.R's config)")
resolved <- resolve_iter_infercnv_obj(cfg$iter_path)
if (!is.null(resolved$error)) stop("iter_path: ", resolved$error)
cat("Using ", resolved$mode, " mode: ", resolved$iter_dir, "\n", sep = "")

# ── Default plot: the object as-is (unchanged from before) ──────────────────
infercnv_res <- readRDS(resolved$obj_path)
out_dir      <- file.path(resolved$iter_dir, "infercnv", "plots")

infercnv::plot_cnv(
  infercnv_obj    = infercnv_res,
  out_dir         = out_dir,
  output_filename = "infercnv_hmm_plot",
  output_format   = "pdf",
  x.range         = "auto",
  x.center        = 1,
  title           = paste("inferCNV HMM —", basename(resolved$iter_dir)),
  color_safe_pal  = FALSE
)
cat("Saved to:", out_dir, "\n")

# ── WITH_DENOISE=TRUE: rerun with inferCNV's own denoise step, then plot that ──
if (with_denoise) {
  message("WITH_DENOISE=TRUE: rerunning the full inferCNV pipeline (denoise = TRUE). ",
          "This is the same cost as scripts/run_infercnv_solo.R, not a plot - ",
          "make sure this job's --time/--mem are sized for a full inferCNV run, ",
          "not a plain plot (see slurm/run_iterations.sh for reference sizing).")

  iter_root      <- dirname(resolved$iter_dir)
  iter_cfg_path  <- file.path(iter_root, "config.yaml")
  if (!file.exists(iter_cfg_path)) {
    stop("WITH_DENOISE=TRUE needs the iteration's own saved config, not found: ", iter_cfg_path)
  }
  iter_cfg <- load_config(iter_cfg_path)
  for (k in c("counts_path", "gene_order")) {
    if (is.null(iter_cfg[[k]])) stop("WITH_DENOISE=TRUE: ", iter_cfg_path, " has no ", k)
  }

  state     <- read_iter_state(resolved$iter_dir)
  counts_mx <- readRDS(iter_cfg$counts_path)
  annot     <- make_annotation(colnames(counts_mx), state$ref_ids, min_ref = iter_cfg$min_ref %||% 50L)
  cat("Reference cells:", sum(annot$split_group == "refs"),
      "| Observation cells:", sum(annot$split_group == "non_refs"), "\n")

  denoise_dir <- file.path(resolved$iter_dir, "infercnv_denoised_check")
  dir.create(denoise_dir, recursive = TRUE, showWarnings = FALSE)

  options(scipen = 100)
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix       = counts_mx,
    annotations_file        = annot,
    gene_order_file         = iter_cfg$gene_order,
    chr_exclude             = c("chrMT", "chrY", "chrX"),
    ref_group_names         = "refs",
    min_max_counts_per_cell = c(100, 1e6)
  )

  t_start <- proc.time()
  infercnv::run(
    infercnv_obj        = infercnv_obj,
    out_dir             = denoise_dir,
    cutoff              = iter_cfg$cutoff %||% 0.1,
    cluster_by_groups   = TRUE,
    HMM                 = TRUE,
    denoise             = TRUE,
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
  log_metric(stage = "infercnv_denoise_check", iter_dir = resolved$iter_dir,
             elapsed_s = round((proc.time() - t_start)[["elapsed"]], 1))

  denoised_obj  <- readRDS(file.path(denoise_dir, "run.final.infercnv_obj"))
  denoised_plot <- file.path(denoise_dir, "plots")

  infercnv::plot_cnv(
    infercnv_obj    = denoised_obj,
    out_dir         = denoised_plot,
    output_filename = "infercnv_hmm_plot_denoised",
    output_format   = "pdf",
    x.range         = "auto",
    x.center        = 1,
    title           = paste("inferCNV HMM (denoised) —", basename(resolved$iter_dir)),
    color_safe_pal  = FALSE
  )
  cat("Denoised rerun saved to:", denoise_dir, "\n")
  cat("Denoised plot saved to:", denoised_plot, "\n")
}
