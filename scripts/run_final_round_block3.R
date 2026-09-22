#!/usr/bin/env Rscript
# scripts/run_round3_benchmark.R





process_cnv_connected <- function(
    grp,
    overlap_method,
    min_overlap
) {
  
  n <- nrow(grp)
  
    gr <- GenomicRanges::GRanges(
    seqnames = rep("chr", n), 
    ranges   = IRanges::IRanges(
      start = grp$start,
      end   = grp$end),
    strand = "*"
  )
  
  hits <- GenomicRanges::findOverlaps(
    gr, gr,
    type   = "any",
    select = "all"
  )
  hits <- hits[
    S4Vectors::queryHits(hits) !=
    S4Vectors::subjectHits(hits)
  ]
  
  g <- igraph::make_empty_graph(
    n        = n,
    directed = FALSE
  )
  
  if (length(hits) > 0L) {
    
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    
    scores <- compute_overlap(
      q_start = grp$start[q_idx],
      q_end   = grp$end[q_idx],
      s_start = grp$start[s_idx],
      s_end   = grp$end[s_idx],
      method  = overlap_method
    )
    
    passing <- scores >= min_overlap
    q_pass  <- q_idx[passing]
    s_pass  <- s_idx[passing]
    
    if (length(q_pass) > 0L) {
      g <- igraph::add_edges(
        g,
        as.vector(rbind(q_pass, s_pass))
      )
    }
  }
  
  components     <- igraph::components(g)
  grp$cnv_equiv_id <- components$membership
  
  return(grp)
}



post_score_merge <- function(
    scored_events,
    cell_col       = "cell_name",
    min_overlap    = 0.75,
    overlap_method = "reciprocal"
) {
  
  cat("=== Post-score merge ===\n")
  cat("  Input rows:", nrow(scored_events), "\n")
  
  # ── Step 1: Deduplicate ───────────────────────────────────────────────────
  cat("\n[1] Deduplicating...\n")
  
  deduped <- scored_events %>%
    dplyr::distinct(
      chr, start, end, cnv_state,
      .data[[cell_col]],
      .keep_all = TRUE
    )
  
  cat("  Before:", nrow(scored_events), "\n")
  cat("  After: ", nrow(deduped), "\n")
  
  # ── Step 2: Connected components ─────────────────────────────────────────
  cat("\n[2] Connected components (overlap=",
      min_overlap, ")...\n")
  
  reclustered <- deduped %>%
    dplyr::group_by(
      chr, cnv_state
    ) %>%
    dplyr::group_modify(~ {
      process_cnv_connected(
        grp            = .x,
        overlap_method = overlap_method,
        min_overlap    = min_overlap
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cnv_equiv_id = paste(
        chr,
        cnv_state, cnv_equiv_id,
        sep = "|"
      )
    )
  
  cat("  Clusters:", dplyr::n_distinct(
    reclustered$cnv_equiv_id), "\n")
  
  # ── Step 3: Standardise boundaries ───────────────────────────────────────
  cat("\n[3] Standardising boundaries...\n")
  
  cluster_bounds <- reclustered %>%
    dplyr::group_by(cnv_equiv_id) %>%
    dplyr::summarise(
      start_std         = min(start, na.rm = TRUE),
      end_std           = max(end,   na.rm = TRUE),
      cnv_length_std    = max(end) - min(start) + 1,
      cnv_length_mb_std = (max(end) -
                           min(start) + 1) / 1e6,
      .groups           = "drop"
    )
  
  result <- reclustered %>%
    dplyr::left_join(
      cluster_bounds,
      by = "cnv_equiv_id"
    ) %>%
    dplyr::mutate(
      start         = start_std,
      end           = end_std,
      cnv_length    = cnv_length_std,
      cnv_length_mb = cnv_length_mb_std
    ) %>%
    dplyr::select(
      -start_std, -end_std,
      -cnv_length_std,
      -cnv_length_mb_std
    )
  
  cat("  Rows:", nrow(result), "\n")
  
  # ── Quick diagnostic ──────────────────────────────────────────────────────
  chr1q <- result %>%
    dplyr::filter(
      chr       == "chr1",
      cnv_state == "gain"
    )

  return(result)
}


#!/usr/bin/env Rscript
# scripts/run_round3_benchmark.R

cat("=== Script started:", format(Sys.time()), "===\n")

source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/round2_pipeline.R")
cat("  round2_pipeline.R ✅\n")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_processing.R")
cat("  cnv_processing.R ✅\n")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_scoring.R")
cat("  cnv_scoring.R ✅\n")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_annotation.R")
cat("  cnv_annotation.R ✅\n")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/pipeline.R")
cat("  pipeline.R ✅\n")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/infercnv.R")
cat("  infercnv.R ✅\n\n")

BASE    <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
OUT_DIR <- file.path(BASE, "cnv_results/VUB04/round3_benchmark_total_cell")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load shared data ──────────────────────────────────────────────────────────
cat("=== Loading shared data ===\n")

cat("  Loading inferCNV object...\n")
obj <- readRDS(file.path(BASE,
  "infercnv_results/VUB04/round3",
  "run.final.infercnv_obj"))
cat("  expr.data dim:", dim(obj@expr.data), "\n")

objs <- list(RPE = list(
  obj@expr.data,
  obj@gene_order
))

cat("  objs built ✅\n")

cat("  Loading metadata...\n")
metadata <- readRDS(file.path(BASE,
  "data/metadata_3round.rds"))
cat("  metadata rows:", nrow(metadata), "\n")
cat("  split_group:\n")
print(table(metadata$split_group, useNA = "always"))

cat("  Loading chromosome arms...\n")
chromosome_arms <- readRDS(file.path(BASE,
  "data/hg38_chromosome_arms.rds"))
cat("  chromosome_arms rows:",
    nrow(chromosome_arms), "\n")

# ── Pre-compute ONCE outside all loops ───────────────────────────────────────
cat("\n=== Pre-computing reference split ===\n")

ref_cell_names <- metadata %>%
  dplyr::filter(split_group == "ref") %>%
  dplyr::pull(cell_name)

metadata_obs <- metadata %>%
  dplyr::filter(split_group != "ref") %>%
  dplyr::mutate(cell_type = "RPE")

cell_sizes <- compute_cell_sizes(
  metadata   = metadata,
  group_cols = "cell_type",
  cell_col   = "cell_name"
)

cat("Reference cells:   ", length(ref_cell_names), "\n")
cat("Observation cells: ", cell_sizes$n_total_cells, "\n\n")

# ── Parameter grid ────────────────────────────────────────────────────────────
k_dis_values      <- c(1.5, 2.0)
k_fre_values      <- c(1.5, 1.75, 2.0, 2.5)
sens_floor_values <- c(17.5, 20, 22.5)
ovlp_values       <- c(0.75, 0.80, 0.85)

n_total <- length(k_dis_values)    *
           length(ovlp_values)     *
           length(k_fre_values)    *
           length(sens_floor_values)

cat(sprintf(paste0(
  "=== Parameter grid ===\n",
  "  k_dis:      %s\n",
  "  k_fre:      %s\n",
  "  sens_floor: %s\n",
  "  ovlp:       %s\n",
  "  Total:      %d combos\n\n"
),
paste(k_dis_values,      collapse = ", "),
paste(k_fre_values,      collapse = ", "),
paste(sens_floor_values, collapse = ", "),
paste(ovlp_values,       collapse = ", "),
n_total
))

t_total       <- proc.time()
results_summary <- list()
combo_counter   <- 0L

# ── Outermost loop: k_dis ─────────────────────────────────────────────────────
for (k_dis in k_dis_values) {
  
  
  cat(sprintf(
    " k_dis=%.2f [%s]\n",
    k_dis, format(Sys.time())
  ))
  cat("========================================\n")
  
  t_dis <- proc.time()
  
  # ── Gene level df ─────────────────────────────────────────────────────────
  cat("  load_and_prepare_infercnv_reference...\n")
  gene_level_df <- tryCatch(
    load_and_prepare_infercnv_reference(
      objs, k = k_dis),
    error = function(e) {
      cat("  ERROR:", e$message, "\n")
      NULL
    }
  )
  if (is.null(gene_level_df)) {
    cat("  Skipping k_dis=", k_dis, "\n")
    next
  }
  cat("  gene_level_df rows:", nrow(gene_level_df), "\n")
  
  # ── Collapse ──────────────────────────────────────────────────────────────
  cat("  Collapsing genes to segments...\n")
  collapse_df <- collapse_genes_to_cnv_segments(
    gene_level_df) %>%
    dplyr::mutate(
      cnv_length    = as.numeric(stop) -
                      as.numeric(start) + 1,
      cnv_length_mb = cnv_length / 1e6
    ) %>%
    dplyr::filter(cnv_length_mb >= 15)
  cat("  Segments (>=15Mb):", nrow(collapse_df), "\n")
  
  # ── Merge ─────────────────────────────────────────────────────────────────
  cat("  Merging nearby regions...\n")
  merged <- merge_nearby_regions(
    df      = collapse_df,
    max_gap = 100000L
  ) %>%
    dplyr::mutate(
      cnv_length    = as.numeric(end) -
                      as.numeric(start) + 1,
      cnv_length_mb = cnv_length / 1e6
    )
  cat("  After merge:", nrow(merged), "\n")
  
  # ── Annotate + filter ref cells ───────────────────────────────────────────
  cat("  Adding chromosome info...\n")
  cnv_annotated_kdis <- add_chromosome_info(
    merged,
    chromosome_arms,
    chr_col   = "chr",
    start_col = "start",
    end_col   = "end"
  ) %>%
    dplyr::filter(
      !cell_name %in% ref_cell_names
    ) %>%
    dplyr::mutate(cell_type = "RPE")
  
  cat(sprintf(
    "  Ready: %d segments, %d cells, %.1f sec\n",
    nrow(cnv_annotated_kdis),
    dplyr::n_distinct(cnv_annotated_kdis$cell_name),
    (proc.time() - t_dis)[["elapsed"]]
  ))
  
  # ── Second loop: ovlp ─────────────────────────────────────────────────────
  for (ovlp in ovlp_values) {
    
    cat(sprintf(
      "\n  --- ovlp=%.2f [%s] ---\n",
      ovlp, format(Sys.time())
    ))
    t_cluster <- proc.time()
    
    cat("    run_cnv_locus_analysis...\n")
    clustered_events <- tryCatch(
      run_cnv_locus_analysis(
        cnv_annotated_kdis,
        by             = "cell_type",
        overlap_method = "reciprocal",
        min_ovelap     = ovlp,
        sample_col     = "cell_type",
        cell_col       = "cell_name"
      ),
      error = function(e) {
        cat("    ERROR clustering:", e$message, "\n")
        NULL
      }
    )
    if (is.null(clustered_events)) {
      cat("    Skipping ovlp=", ovlp, "\n")
      next
    }
    
    clustered_events$clustered_events$cell_type <- "RPE"
    clustered_events$cnv_locus_summary$cell_type <- "RPE"
    
    cat(sprintf(
      "    Done in %.1f sec — %d loci\n",
      (proc.time() - t_cluster)[["elapsed"]],
      nrow(clustered_events$cnv_locus_summary)
    ))
    
    # ── Inner loops: k_fre × sens_floor ──────────────────────────────────
    for (k_fre in k_fre_values) {
    for (sens_floor in sens_floor_values) {
      
      combo_counter <- combo_counter + 1L
      
      combo_id <- sprintf(
        "kdis%.2f_kfre%.2f_floor%.1f_ovlp%.2f",
        k_dis, k_fre, sens_floor, ovlp)
      
      cat(sprintf(
        "    [%d/%d] %s [%s]\n",
        combo_counter, n_total,
        combo_id, format(Sys.time())
      ))
      
      out_file <- file.path(OUT_DIR,
        paste0(combo_id, "_scored.rds"))
      
      if (file.exists(out_file)) {
        cat("      Skip (done) ✅\n")
        scored <- tryCatch(
          readRDS(out_file),
          error = function(e) NULL)
      } else {
        
        t_score <- proc.time()
        
        scored <- tryCatch({
          
          # ── Step 1: Thresholds ────────────────────────────────────────
          cat("      [1] prepare_cnv_thresholds...\n")
          thresholded_df <- prepare_cnv_thresholds(
            summary_df           = clustered_events$cnv_locus_summary,
            clustered_events     = clustered_events$clustered_events,
            by_union             = "cell_type",
            cell_sizes           = cell_sizes,
            k                    = k_fre,
            sensitivity_floor_mb = sens_floor,
            min_required_cells   = 3,
            round_fun            = ceiling
          )
          cat("      thresholded rows:", nrow(thresholded_df), "\n")
          
          # ── Step 2: Arm percentages ───────────────────────────────────
          cat("      [2] add_arm_percentages...\n")
          thresholded_df <- add_arm_percentages(
            cnv_df          = thresholded_df,
            chromosome_arms = chromosome_arms
          )
          
          # ── Step 3: Filter ────────────────────────────────────────────
          cat("      [3] filter_cnv_loci...\n")
          filtered_df <- filter_cnv_loci(
            clustered_events     = thresholded_df,
            p_arm_permission     = 65,
            q_arm_permission     = 65,
            whole_chr_permission = 60
          )
          cat("      After filter:", nrow(filtered_df), "rows\n")
          
          if (nrow(filtered_df) == 0) {
            cat("      No events after filter — skipping\n")
            return(NULL)
          }
          
          # ── Step 4: Distinct rows ─────────────────────────────────────
          cat("      [4] distinct rows...\n")
          deduped_df <- filtered_df %>%
            dplyr::distinct(
              cell_name, chr, cnv_state,
              start, end,
              .keep_all = TRUE
            )
          cat("      After distinct:", nrow(deduped_df), "rows\n")
          
          # ── Step 5: Post-score merge ──────────────────────────────────
          # Recluster survivors + standardise boundaries
          cat("      [5] post_score_merge...\n")
          post_score_merge(
            scored_events   = deduped_df,
            cell_col        = "cell_name",
            min_overlap     = 0.75,
            overlap_method  = "reciprocal"
          )
          
        }, error = function(e) {
          cat(sprintf("      ERROR: %s\n", e$message))
          NULL
        })
        
        runtime_score <- (proc.time() -
                          t_score)[["elapsed"]]
        
        if (!is.null(scored)) {
          saveRDS(scored, out_file)
          cat(sprintf(
            "      Saved in %.1f sec — %d rows ✅\n",
            runtime_score, nrow(scored)))
        } else {
          cat("      Failed ❌\n")
        }
      }
      
      # ── Metrics ───────────────────────────────────────────────────────
      if (is.null(scored) || nrow(scored) == 0) {
        results_summary[[
          length(results_summary)+1]] <-
          data.frame(
            k_dis          = k_dis,
            k_fre          = k_fre,
            sens_floor     = sens_floor,
            ovlp           = ovlp,
            combo_id       = combo_id,
            n_rows         = 0L,
            n_gain         = 0L,
            n_loss         = 0L,
            chr1q_detected = FALSE,
            chr1q_n_cells  = NA_integer_,
            chr1q_length_mb = NA_real_,
            stringsAsFactors = FALSE
          )
        next
      }
      
      chr1q <- scored %>%
        dplyr::filter(
          chr       == "chr1",
          cnv_state == "gain"
        )
      
      chr1q_detected <- nrow(chr1q) > 0
      if (chr1q_detected) {
        cat(sprintf(
          "      ✅ Chr1q: %d cells, %.1f Mb\n",
          dplyr::n_distinct(chr1q$cell_name),
          max(chr1q$cnv_length_mb, na.rm = TRUE)
        ))
      }
      
      results_summary[[
        length(results_summary)+1]] <-
        data.frame(
          k_dis           = k_dis,
          k_fre           = k_fre,
          sens_floor      = sens_floor,
          ovlp            = ovlp,
          combo_id        = combo_id,
          n_rows          = nrow(scored),
          n_gain          = sum(scored$cnv_state == "gain",
                                na.rm = TRUE),
          n_loss          = sum(scored$cnv_state == "loss",
                                na.rm = TRUE),
          chr1q_detected  = chr1q_detected,
          chr1q_n_cells   = if (chr1q_detected)
            dplyr::n_distinct(chr1q$cell_name)
            else NA_integer_,
          chr1q_length_mb = if (chr1q_detected)
            max(chr1q$cnv_length_mb, na.rm = TRUE)
            else NA_real_,
          stringsAsFactors = FALSE
        )
    }
    }
    
    cat(sprintf(
      "  ovlp=%.2f complete [%s]\n",
      ovlp, format(Sys.time())
    ))
  }
  
  cat(sprintf(
    "k_dis=%.2f complete [%s]\n",
    k_dis, format(Sys.time())
  ))
}

# ── Save summary ───────────────────────────────────────────────────────────────
cat("\n=== Saving summary ===\n")

summary_df <- dplyr::bind_rows(results_summary)

summary_path <- file.path(OUT_DIR,
  "round3_benchmark_summary.csv")
write.csv(summary_df, summary_path,
          row.names = FALSE)

total_runtime <- (proc.time() - t_total)[["elapsed"]]

cat(sprintf(paste0(
  "\n================================================\n",
  " COMPLETE [%s]\n",
  "================================================\n",
  "  Combinations:    %d / %d\n",
  "  Chr1q detected:  %d\n",
  "  Total runtime:   %.1f min\n",
  "  Summary:         %s\n",
  "================================================\n"
),
format(Sys.time()),
nrow(summary_df),
n_total,
sum(summary_df$chr1q_detected, na.rm = TRUE),
total_runtime / 60,
summary_path
))

cat("\n=== Best combinations ===\n")
summary_df %>%
  dplyr::filter(chr1q_detected) %>%
  dplyr::arrange(dplyr::desc(chr1q_n_cells)) %>%
  dplyr::select(combo_id, chr1q_n_cells,
                chr1q_length_mb, n_rows,
                n_gain, n_loss) %>%
  head(20) %>%
  print()