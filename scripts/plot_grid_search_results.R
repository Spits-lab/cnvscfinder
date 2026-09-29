#!/usr/bin/env Rscript
# =============================================================================
# scripts/plot_grid_search_results.R
# Generalised successor to R/ploting_ploting.R: loops over every combo's
# *_scored.rds in a grid search's out_dir and plots them all into one PDF, so
# they can be flipped through side by side. Reuses the SAME config as
# scripts/run_grid_search.R (out_dir, chromosome_arms_path, cell_col,
# cell_group_cluster) - no separate plotting config.
#
# Fix in passing (R/ploting_ploting.R): the empty-file skip used return(NULL)
# inside a bare for loop at the script's top level - not inside any function,
# so if it ever actually fired it would raise "no function to return from"
# and kill the whole script instead of skipping to the next file. Uses next
# here instead.
#
# Usage: Rscript plot_grid_search_results.R --config cfg.yaml [--override '{...}']
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})
for (f in c("cnv_annotation.R", "cnv_processing.R", "cnv_scoring.R",
            "ploting_functions.R", "config.R")) {
  source(file.path(BASE, "R", f))
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL)
)))
cfg <- load_config(opt$config, opt$override)

overlap_ok <- function(m) !inherits(try(compute_overlap(1L, 2L, 1L, 2L, method = m), silent = TRUE), "try-error")
assert_valid_grid_search_config(cfg, overlap_ok = overlap_ok)

test_dir <- cfg$out_dir

# ── Find the combos ──────────────────────────────────────────────────────────
cat("Files in ", test_dir, ":\n", sep = "")
scored_files <- list.files(test_dir, pattern = "_scored\\.rds$", full.names = TRUE)
cat("Files found:", length(scored_files), "\n")
if (length(scored_files) == 0L) stop("no *_scored.rds files in ", test_dir, " - run the grid search first")

# ── Genome structure ─────────────────────────────────────────────────────────
chromosome_arms  <- readRDS(cfg$chromosome_arms_path)
genome_structure <- prepare_genome_structure(chromosome_arms)

# ── Open PDF ──────────────────────────────────────────────────────────────────
pdf_path <- file.path(test_dir, "all_plots.pdf")

tryCatch({

  pdf(pdf_path, width = 20, height = 12)

  for (f in scored_files) {

    combo_id <- gsub("_scored\\.rds$", "", basename(f))
    cat("Plotting:", combo_id, "\n")

    tryCatch({

      scored <- readRDS(f)

      if (is.null(scored) || nrow(scored) == 0) {
        cat("  Empty — skipping\n")
        next
      }

      n_total_cells <- dplyr::n_distinct(scored[[cfg$cell_col]])

      # No chr1q callout: that was a VUB04-specific ground-truth check
      # (VUB04's known true event was a 1q gain) and doesn't apply here.
      title <- sprintf(
        "%s\ncells=%d | gain=%d | loss=%d",
        combo_id, n_total_cells,
        sum(scored$cnv_state == "gain", na.rm = TRUE),
        sum(scored$cnv_state == "loss", na.rm = TRUE)
      )

      cnv_mapped <- map_cnv_to_genome(
        scored, genome_structure, threshold = 85, arrange_df_cols = cfg$cell_group_cluster
      )

      if (is.null(cnv_mapped) || nrow(cnv_mapped) == 0) {
        cat("  Empty cnv_mapped — skipping\n")
        next
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

      print(p)
      cat("  Done ✅\n")

    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      cat("  Skipping — continuing loop\n")
    })
  }

}, finally = {
  if (dev.cur() != 1) {
    dev.off()
    cat("\nPDF device closed\n")
  }
})

cat("Saved:", pdf_path, "\n")
