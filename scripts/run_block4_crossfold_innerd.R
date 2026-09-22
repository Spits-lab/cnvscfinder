#!/usr/bin/env Rscript
# scripts/run_block4_crossfold_inner.R
# Runs ALL block4 combinations for ONE block2 task
# Clustering hoisted out of inner loop — only reruns per ovlp value

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

option_list <- list(
  optparse::make_option("--task-dir", type = "character", default = NULL),
  optparse::make_option("--b3-path",  type = "character", default = NULL),
  optparse::make_option("--cs-path",  type = "character", default = NULL)
)

opt <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list))

if (is.null(opt$`task-dir`)) stop("--task-dir required")
if (is.null(opt$`b3-path`))  stop("--b3-path required")
if (is.null(opt$`cs-path`))  stop("--cs-path required")

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
source(file.path(BASE, "R/cnv_annotation.R"))
source(file.path(BASE, "R/cnv_processing.R"))
source(file.path(BASE, "R/cnv_scoring.R"))
source(file.path(BASE, "R/pipeline.R"))

cat("Loading block3 from:", opt$`b3-path`, "\n")
cnv_annotated <- readRDS(opt$`b3-path`)
cnv_annotated$cell_type <- "RPE"

cat("Loading cell_sizes from:", opt$`cs-path`, "\n")
cell_sizes <- readRDS(opt$`cs-path`)

chromosome_arms <- readRDS(
  file.path(BASE, "data/hg38_chromosome_arms.rds"))

cat("cnv_annotated rows:", nrow(cnv_annotated), "\n")
cat("cell_sizes:\n")
print(cell_sizes)

# ── Parameter grid ────────────────────────────────────────────────────────────
k_fre_values      <- c(1.5, 1.75, 2.0, 2.25, 2.5)
sens_floor_values <- c(17.5, 20, 22.5, 25)
ovlp_values       <- c(0.75, 0.80, 0.825, 0.85, 0.90)

n_total <- length(k_fre_values) * length(sens_floor_values) * length(ovlp_values)
cat("Total combinations:", n_total, "\n\n")

t_total <- proc.time()
results_summary <- list()

# ── Outer loop: ovlp — clustering runs once per value ─────────────────────────
for (ovlp in ovlp_values) {

  cat(sprintf("\n=== Clustering for ovlp=%.3f ===\n", ovlp))
  t_cluster <- proc.time()

  clustered <- run_cnv_locus_analysis(
    cnv_annotated,
    by             = "cell_type",
    overlap_method = "reciprocal",
    min_ovelap     = ovlp,
    sample_col     = "cell_type",
    cell_col       = "cell_name"
  )

  cat(sprintf("  Clustering done in %.1f sec\n",
              (proc.time() - t_cluster)[["elapsed"]]))

  # ── Inner loops: scoring only — fast ─────────────────────────────────────
  for (k_fre in k_fre_values) {
  for (sens_floor in sens_floor_values) {

    out_dir <- file.path(
      opt$`task-dir`, "block4_crossfold",
      sprintf("kfre%.2f_floor%.2f_ovlp%.3f",
              k_fre, sens_floor, ovlp)
    )
    out_file <- file.path(out_dir, "scored_events.rds")

    if (file.exists(out_file)) {
      cat(sprintf("  Skip (done): kfre=%.2f floor=%.1f ovlp=%.3f\n",
                  k_fre, sens_floor, ovlp))
      scored <- tryCatch(readRDS(out_file), error = function(e) NULL)

    } else {

      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      t_start <- proc.time()

      scored <- tryCatch({

        scored_events <- score_cnv_clusters(
          summary_df           = clustered$cnv_locus_summary,
          clustered_events     = clustered$clustered_events,
          cell_sizes           = cell_sizes,
          by_union             = "cell_type",
          chromosome_arms      = chromosome_arms,
          k                    = k_fre,
          sensitivity_floor_mb = sens_floor,
          min_required_cells   = 3L,
          p_arm_permission     = 60,
          q_arm_permission     = 60,
          whole_chr_permission = 65,
          round_fun            = ceiling
        )

        # dedup overlap fixed at 0.75 — not benchmarked
        deduplicate_cnv_cells(
          scored_events,
          sample_col  = "cell_type",
          cell_col    = "cell_name",
          min_overlap = 0.75
        )

      }, error = function(e) {
        cat("    ERROR:", e$message, "\n")
        NULL
      })

      runtime <- (proc.time() - t_start)[["elapsed"]]

      if (!is.null(scored)) {
        saveRDS(scored, out_file)
        cat(sprintf("    Scored in %.2f sec — %d loci\n",
                    runtime, nrow(scored)))
      }
    }

    # ── Metrics ────────────────────────────────────────────────────────────
    if (is.null(scored) || nrow(scored) == 0) {
      results_summary[[length(results_summary) + 1]] <- data.frame(
        k_fre = k_fre, sens_floor = sens_floor, ovlp = ovlp,
        n_loci = 0L, n_gain = 0L, n_loss = 0L,
        chr1q_detected  = FALSE,
        chr1q_n_cells   = NA_integer_,
        chr1q_length_mb = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    chr1q <- scored %>%
      dplyr::filter(chr == "chr1", arm_class == "q_arm",
                    cnv_state == "gain")

    results_summary[[length(results_summary) + 1]] <- data.frame(
      k_fre  = k_fre, sens_floor = sens_floor, ovlp = ovlp,
      n_loci = nrow(scored),
      n_gain = sum(scored$cnv_state == "gain", na.rm = TRUE),
      n_loss = sum(scored$cnv_state == "loss", na.rm = TRUE),
      chr1q_detected  = nrow(chr1q) > 0,
      chr1q_n_cells   = if (nrow(chr1q) > 0) max(chr1q$n_cells, na.rm = TRUE) else NA_integer_,
      chr1q_length_mb = if (nrow(chr1q) > 0) max(chr1q$cnv_length_mb, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  }
}

# ── Save summary ───────────────────────────────────────────────────────────────
summary_df <- dplyr::bind_rows(results_summary)

summary_path <- file.path(opt$`task-dir`, "block4_crossfold_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

# ── Build combined scored events table for plotting ────────────────────────────
cat("\nBuilding combined scored events table...\n")

all_scored <- lapply(seq_len(nrow(summary_df)), function(i) {
  row <- summary_df[i, ]
  out_file <- file.path(
    opt$`task-dir`, "block4_crossfold",
    sprintf("kfre%.2f_floor%.2f_ovlp%.3f",
            row$k_fre, row$sens_floor, row$ovlp),
    "scored_events.rds"
  )
  if (!file.exists(out_file)) return(NULL)

  scored <- tryCatch(readRDS(out_file), error = function(e) NULL)
  if (is.null(scored) || nrow(scored) == 0) return(NULL)

  scored$k_fre      <- row$k_fre
  scored$sens_floor <- row$sens_floor
  scored$ovlp       <- row$ovlp
  scored$combo_id   <- sprintf("kfre%.2f_floor%.2f_ovlp%.3f",
                               row$k_fre, row$sens_floor, row$ovlp)
  scored
})

all_scored_df <- dplyr::bind_rows(Filter(Negate(is.null), all_scored))

scored_path <- file.path(opt$`task-dir`, "block4_crossfold_all_scored.rds")
saveRDS(all_scored_df, scored_path)

total_runtime <- (proc.time() - t_total)[["elapsed"]]

cat("\n================================================\n")
cat(" Task complete:", basename(opt$`task-dir`), "\n")
cat(" Combinations:  ", n_total, "\n")
cat(" Chr1q detected:", sum(summary_df$chr1q_detected, na.rm = TRUE), "\n")
cat(" Total runtime: ", round(total_runtime / 60, 1), "min\n")
cat(" Summary:       ", summary_path, "\n")
cat(" All scored:    ", scored_path, "\n")
cat("================================================\n")