source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_annotation.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_processing.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_scoring.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/ploting_functions.R")


BASE    <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
TEST_DIR <- file.path(BASE, "cnv_results/VUB04/round7/test_run_0907_ref")

# ── List available files ──────────────────────────────────────────────────────
cat("Files in test_run:\n")
list.files(TEST_DIR, pattern = "_scored\\.rds$")

# ── Load file ─────────────────────────────────────────────────────────────────
scored_files <- list.files(
  TEST_DIR,
  pattern    = "_scored\\.rds$",
  full.names = TRUE
)
# ── Genome structure ──────────────────────────────────────────────────────────
chromosome_arms <- readRDS(file.path(BASE,
  "data/hg38_chromosome_arms.rds"))
  
genome_structure <- prepare_genome_structure(chromosome_arms)

cat("Files found:", length(scored_files), "\n")

# ── Open PDF ──────────────────────────────────────────────────────────────────
pdf_path <- file.path(TEST_DIR, "all_plots.pdf")

tryCatch({
  
  pdf(pdf_path, width = 20, height = 12)
  
  for (f in scored_files) {
    
    combo_id <- gsub("_scored\\.rds$", "",
                     basename(f))
    cat("Plotting:", combo_id, "\n")
    
    tryCatch({
      
      scored <- readRDS(f)
      
      if (is.null(scored) || nrow(scored) == 0) {
        cat("  Empty — skipping\n")
        return(NULL)
      }
      
      chr1q <- scored %>%
        dplyr::filter(chr       == "chr1",
                      cnv_state == "gain")
      
      n_chr1q_cells <- dplyr::n_distinct(
        chr1q$cell_name)
      n_total_cells <- dplyr::n_distinct(
        scored$cell_name)
      
      title <- sprintf(
        "%s\ncells=%d | gain=%d | loss=%d | chr1q=%s | chr1q_cells=%d (%.1f%%)",
        combo_id,
        n_total_cells,
        sum(scored$cnv_state == "gain",
            na.rm = TRUE),
        sum(scored$cnv_state == "loss",
            na.rm = TRUE),
        if (nrow(chr1q) > 0)
          paste0(round(max(chr1q$cnv_length_mb,
                           na.rm = TRUE), 1), "Mb")
        else "not detected",
        n_chr1q_cells,
        100 * n_chr1q_cells /
          max(n_total_cells, 1)
      )
      
      cnv_mapped <- map_cnv_to_genome(
        scored,
        genome_structure,
        threshold       = 85,
        arrange_df_cols = "cell_type"
      )
      
      if (is.null(cnv_mapped) ||
          nrow(cnv_mapped) == 0) {
        cat("  Empty cnv_mapped — skipping\n")
        return(NULL)
      }
      
      heatmap_plot <- prepare_cnv_plot(
        cnv_mapped,
        genome_structure,
        grouping_cols = "cell_type",
        state_colors  = c(
          "gain" = "#E64B35",
          "loss" = "royalblue4"
        )
      )
      
      p <- plot_cnv_karyotype(
        heatmap_plot,
        genome_structure,
        ideogram_ratio = 0.05,
        arm_colors     = c(
          "p"   = "paleturquoise4",
          "cen" = "black",
          "q"   = "red4"
        ),
        show_legend = TRUE
      ) + ggplot2::ggtitle(title)
      
      print(p)
      cat("  Done ✅\n")
      
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      cat("  Skipping —  continuing loop\n")
    })
  }

}, finally = {
  # ── Always close device even if error ──────────────────────────────────────
  if (dev.cur() != 1) {
    dev.off()
    cat("\nPDF device closed\n")
  }
})

cat("Saved:", pdf_path, "\n")