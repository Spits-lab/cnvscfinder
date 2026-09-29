#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_grid_search.R
# Generalised successor to legacy/scripts/run_block_test.R: a nested grid search over
# k_dis/k_fre/sens_floor/ovlp/pct_floor/min_expr_density, driven entirely by a
# YAML config (dataset, paths, grid values and fixed knobs) instead of being
# hard-coded to one dataset. No metrics/ground-truth scoring here (see
# configs/*.yaml comments) - just a summary CSV (n_rows/n_gain/n_loss/runtime
# per combo) and each combo's scored .rds, for later inspection.
#
# process_cnv_connected()/post_score_merge() below duplicate R/iteration.R's
# copies (deferred consolidation - see the plan's "Deferred findings").
#
# Usage: Rscript run_grid_search.R --config cfg.yaml [--override '{...}']
# =============================================================================

BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})
for (f in c("cnv_processing.R", "cnv_scoring.R", "cnv_annotation.R",
            "pipeline.R", "infercnv.R", "config.R")) {
  source(file.path(BASE, "R", f))
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config",   type = "character"),
  make_option("--override", type = "character", default = NULL)
)))
cfg <- load_config(opt$config, opt$override)

overlap_ok <- function(m) !inherits(try(compute_overlap(1L, 2L, 1L, 2L, method = m), silent = TRUE), "try-error")
assert_valid_grid_search_config(cfg, overlap_ok = overlap_ok)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ── process_cnv_connected() / post_score_merge() (from R/round2.R) ──────────

process_cnv_connected <- function(grp, overlap_method, min_overlap) {
  n  <- nrow(grp)
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr", n),
    ranges   = IRanges::IRanges(start = grp$start, end = grp$end),
    strand   = "*"
  )
  hits <- GenomicRanges::findOverlaps(gr, gr, type = "any", select = "all")
  hits <- hits[S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits)]

  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (length(hits) > 0L) {
    q_idx  <- S4Vectors::queryHits(hits)
    s_idx  <- S4Vectors::subjectHits(hits)
    scores <- compute_overlap(
      q_start = grp$start[q_idx], q_end = grp$end[q_idx],
      s_start = grp$start[s_idx], s_end = grp$end[s_idx],
      method  = overlap_method
    )
    passing <- scores >= min_overlap
    if (any(passing)) {
      g <- igraph::add_edges(g, as.vector(rbind(q_idx[passing], s_idx[passing])))
    }
  }
  grp$cnv_equiv_id <- igraph::components(g)$membership
  grp
}

post_score_merge <- function(scored_events, cell_col = "cell_name",
                             min_overlap = 0.75, overlap_method = "adaptive") {
  deduped <- scored_events %>%
    dplyr::distinct(chr, start, end, cnv_state, .data[[cell_col]], .keep_all = TRUE)

  reclustered <- deduped %>%
    dplyr::group_by(chr, cnv_state) %>%
    dplyr::group_modify(~ process_cnv_connected(.x, overlap_method, min_overlap)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cnv_equiv_id = paste(chr, cnv_state, cnv_equiv_id, sep = "|"))

  bounds <- reclustered %>%
    dplyr::group_by(cnv_equiv_id) %>%
    dplyr::summarise(start_std = min(start, na.rm = TRUE),
                     end_std   = max(end,   na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(len_std = end_std - start_std + 1)

  reclustered %>%
    dplyr::left_join(bounds, by = "cnv_equiv_id") %>%
    dplyr::mutate(start = start_std, end = end_std,
                  cnv_length = len_std, cnv_length_mb = len_std / 1e6) %>%
    dplyr::select(-start_std, -end_std, -len_std)
}

# ── Load shared data ─────────────────────────────────────────────────────────
cat("\n=== Loading inferCNV object ===\n")
resolved <- resolve_iter_infercnv_obj(cfg$iter_path)
if (!is.null(resolved$error)) stop("iter_path: ", resolved$error)
cat("Using ", resolved$mode, " mode: ", resolved$iter_dir, "\n", sep = "")

obj  <- readRDS(resolved$obj_path)
objs <- setNames(list(list(obj@expr.data, obj@gene_order)), cfg$group_value)

metadata        <- readRDS(cfg$metadata_path)
chromosome_arms <- readRDS(cfg$chromosome_arms_path)
coding_genes    <- readRDS(cfg$coding_genes_path)

cell_sizes <- compute_cell_sizes(metadata, group_cols = cfg$cell_group_cluster, cell_col = cfg$cell_col)
cat("Cells: ", cell_sizes$n_total_cells, "\n\n")

# ── Gene sets (shared across the whole grid) ─────────────────────────────────
coding_gr <- GenomicRanges::GRanges(
  seqnames = coding_genes$chr,
  ranges   = IRanges::IRanges(start = coding_genes$start, end = coding_genes$end)
)
coding_gr$gene <- coding_genes$gene_name

expressed_table <- as.data.frame(objs[[cfg$group_value]][[2]])
expressed_table$gene <- rownames(expressed_table)
colnames(expressed_table) <- c("chr", "start", "stop", "gene")

expressed_gr <- GenomicRanges::GRanges(
  seqnames = expressed_table$chr,
  ranges   = IRanges::IRanges(start = expressed_table$start, end = expressed_table$stop)
)
expressed_gr$gene <- expressed_table$gene

coding_expressed_set <- intersect(unique(coding_genes$gene_name), unique(expressed_table$gene))
cat("Coding genes (GTF):      ", length(unique(coding_genes$gene_name)), "\n")
cat("Expressed (@gene_order): ", length(unique(expressed_table$gene)), "\n")
cat("Coding + expressed:      ", length(coding_expressed_set), "\n\n")

# ── Parameter grid ────────────────────────────────────────────────────────────
n_total <- length(cfg$k_dis_values) * length(cfg$ovlp_values) * length(cfg$k_fre_values) *
  length(cfg$sens_floor_values) * length(cfg$pct_floor_values) * length(cfg$min_expr_density_values)

cat(sprintf(paste0(
  "=== Parameter grid ===\n",
  "  k_dis:      %s\n", "  k_fre:      %s\n", "  sens_floor: %s\n", "  ovlp:       %s\n",
  "  pct_floor:  %s\n", "  min_dens:   %s\n", "  Total:      %d combos\n\n"
),
paste(cfg$k_dis_values, collapse = ", "), paste(cfg$k_fre_values, collapse = ", "),
paste(cfg$sens_floor_values, collapse = ", "), paste(cfg$ovlp_values, collapse = ", "),
paste(cfg$pct_floor_values, collapse = ", "), paste(cfg$min_expr_density_values, collapse = ", "),
n_total
))

t_total         <- proc.time()
results_summary <- list()
combo_counter   <- 0L

# ── Outermost loop: k_dis ─────────────────────────────────────────────────────
for (k_dis in cfg$k_dis_values) {

  cat("  load_and_prepare_infercnv_reference...\n")
  gene_level_df <- tryCatch(
    load_and_prepare_infercnv_reference(objs, k = k_dis),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )
  if (is.null(gene_level_df)) { cat("  Skipping k_dis=", k_dis, "\n"); next }
  cat("  gene_level_df rows:", nrow(gene_level_df), "\n")

  collapse_df <- collapse_genes_to_cnv_segments(gene_level_df) %>%
    dplyr::mutate(cnv_length = as.numeric(stop) - as.numeric(start) + 1,
                  cnv_length_mb = cnv_length / 1e6)
  cat("  Segments:", nrow(collapse_df), "\n")

  for (pct_floor in cfg$pct_floor_values) {
  for (min_expr_density in cfg$min_expr_density_values) {

    seg_gr <- GenomicRanges::GRanges(
      seqnames = collapse_df$chr,
      ranges   = IRanges::IRanges(start = collapse_df$start, end = collapse_df$stop)
    )
    seg_gr$seg_idx <- seq_len(nrow(collapse_df))

    cat("Counting coding genes per segment...\n")
    hits_coding <- GenomicRanges::findOverlaps(seg_gr, coding_gr, type = "any")

    if (length(hits_coding) == 0L) {
      message("  Gene density filter: no coding genes overlap any segment — skipping this combo")
      next
    }

    segments <- filter_segments_by_gene_density(
      collapse_df = collapse_df, gene_order = expressed_table, coding_gr = coding_gr,
      coding_expressed_set = coding_expressed_set, pct_max = cfg$pct_max, pct_floor = pct_floor,
      min_expr_density = min_expr_density, min_coding_density = cfg$min_coding_density
    )

    cat("  Merging nearby regions...\n")
    merged <- merge_and_density_check(
      df = segments, gene_order = expressed_table, coding_gr = coding_gr,
      coding_expressed_set = coding_expressed_set, max_gap_mb = cfg$max_gap_mb,
      pct_max = cfg$pct_max, pct_floor = pct_floor, min_expr_density = min_expr_density,
      min_coding_density = cfg$min_coding_density
    )

    merged <- merged %>%
      dplyr::mutate(cnv_length = end - start, cnv_length_mb = cnv_length / 1e6) %>%
      dplyr::filter(cnv_length_mb >= cfg$min_segment_mb)

    cat("  Adding chromosome info...\n")
    cnv_annotated_kdis <- add_chromosome_info(
      merged, chromosome_arms, chr_col = "chr", start_col = "start", end_col = "end"
    ) %>%
      dplyr::mutate(!!cfg$cell_group_cluster := cfg$group_value)

    # ── Second loop: ovlp ─────────────────────────────────────────────────────
    for (ovlp in cfg$ovlp_values) {

      cat(sprintf("\n  --- ovlp=%.2f [%s] ---\n", ovlp, format(Sys.time())))
      t_cluster <- proc.time()

      cat("    run_cnv_locus_analysis...\n")
      clustered_events <- tryCatch(
        run_cnv_locus_analysis(
          cnv_annotated_kdis, by = cfg$cell_group_cluster, overlap_method = cfg$overlap_method,
          min_overlap = ovlp, sample_col = cfg$cell_group_cluster, cell_col = cfg$cell_col,
          range = cfg$range
        ),
        error = function(e) { cat("    ERROR clustering:", e$message, "\n"); NULL }
      )
      if (is.null(clustered_events)) { cat("    Skipping ovlp=", ovlp, "\n"); next }

      clustered_events$clustered_events[[cfg$cell_group_cluster]]  <- cfg$group_value
      clustered_events$cnv_locus_summary[[cfg$cell_group_cluster]] <- cfg$group_value

      cat(sprintf("    Done in %.1f sec — %d loci\n",
                  (proc.time() - t_cluster)[["elapsed"]], nrow(clustered_events$cnv_locus_summary)))

      # ── Inner loops: k_fre × sens_floor ──────────────────────────────────
      for (k_fre in cfg$k_fre_values) {
      for (sens_floor in cfg$sens_floor_values) {

        combo_counter <- combo_counter + 1L
        combo_id <- sprintf(
          "kdis%.2f_kfre%.2f_floor%.1f_ovlp%.2f_pctfl%.0f_mindensity%.1f",
          k_dis, k_fre, sens_floor, ovlp, pct_floor, min_expr_density
        )
        cat(sprintf("    [%d/%d] %s [%s]\n", combo_counter, n_total, combo_id, format(Sys.time())))

        out_file <- file.path(cfg$out_dir, paste0(combo_id, "_scored.rds"))
        runtime_score <- NA_real_

        if (file.exists(out_file)) {
          cat("      Skip (done) ✅\n")
          scored <- tryCatch(readRDS(out_file), error = function(e) NULL)
        } else {
          t_score <- proc.time()
          scored <- tryCatch({

            cat("      [1] prepare_cnv_thresholds...\n")
            thresholded_df <- prepare_cnv_thresholds(
              summary_df = clustered_events$cnv_locus_summary, clustered_events = clustered_events$clustered_events,
              by_union = cfg$cell_group_cluster, cell_sizes = cell_sizes, k = k_fre,
              sensitivity_floor_mb = sens_floor, min_required_cells = cfg$min_required_cells, round_fun = ceiling
            )
            cat("      thresholded rows:", nrow(thresholded_df), "\n")

            cat("      [2] add_arm_percentages...\n")
            thresholded_df <- add_arm_percentages(cnv_df = thresholded_df, chromosome_arms = chromosome_arms)

            cat("      [3] filter_cnv_loci...\n")
            filtered_df <- filter_cnv_loci(
              clustered_events = thresholded_df, p_arm_permission = cfg$p_arm_permission,
              q_arm_permission = cfg$q_arm_permission, whole_chr_permission = cfg$whole_chr_permission
            )
            cat("      After filter:", nrow(filtered_df), "rows\n")
            if (nrow(filtered_df) == 0) { cat("      No events after filter — skipping\n"); return(NULL) }

            cat("      [4] distinct rows...\n")
            deduped_df <- filtered_df %>%
              dplyr::distinct(.data[[cfg$cell_col]], chr, cnv_state, start, end, .keep_all = TRUE)
            cat("      After distinct:", nrow(deduped_df), "rows\n")

            cat("      [5] post_score_merge...\n")
            post_score_merge(
              scored_events = deduped_df, cell_col = cfg$cell_col,
              min_overlap = cfg$post_min_overlap, overlap_method = cfg$post_overlap_method
            )

          }, error = function(e) { cat(sprintf("      ERROR: %s\n", e$message)); NULL })

          runtime_score <- (proc.time() - t_score)[["elapsed"]]
          if (!is.null(scored)) {
            saveRDS(scored, out_file)
            cat(sprintf("      Saved in %.1f sec — %d rows ✅\n", runtime_score, nrow(scored)))
          } else {
            cat("      Failed ❌\n")
          }
        }

        # ── Bookkeeping (no ground truth for this dataset - counts + runtime only) ──
        results_summary[[length(results_summary) + 1L]] <- data.frame(
          k_dis = k_dis, k_fre = k_fre, sens_floor = sens_floor, ovlp = ovlp,
          pct_floor = pct_floor, min_expr_density = min_expr_density, combo_id = combo_id,
          n_rows = if (is.null(scored)) 0L else nrow(scored),
          n_gain = if (is.null(scored)) 0L else sum(scored$cnv_state == "gain", na.rm = TRUE),
          n_loss = if (is.null(scored)) 0L else sum(scored$cnv_state == "loss", na.rm = TRUE),
          runtime_s = runtime_score,
          stringsAsFactors = FALSE
        )

      }  # sens_floor
      }  # k_fre
    }    # ovlp
  }      # min_expr_density
  }      # pct_floor
}        # k_dis

# ── Save summary ───────────────────────────────────────────────────────────────
cat("\n=== Saving summary ===\n")
summary_df   <- dplyr::bind_rows(results_summary)
summary_path <- file.path(cfg$out_dir, "grid_search_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

total_runtime <- (proc.time() - t_total)[["elapsed"]]
cat(sprintf(paste0(
  "\n================================================\n",
  " COMPLETE [%s]\n",
  "================================================\n",
  "  Combinations: %d / %d\n",
  "  Total runtime: %.1f min\n",
  "  Summary:      %s\n",
  "================================================\n"
), format(Sys.time()), nrow(summary_df), n_total, total_runtime / 60, summary_path))

print(summary_df)
