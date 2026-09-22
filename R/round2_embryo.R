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


source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_processing.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_scoring.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_annotation.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/pipeline.R")
source("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/infercnv.R")


BASE    <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"
OUT_DIR <- file.path(BASE, "cnv_results/TE/single_run")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/infercnv_results/TE/single_run/run.final.infercnv_obj")
objs <- list(TE = list(obj@expr.data, obj@gene_order))

metadata <- readRDS("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/data/metadata1208_TE.rds")
chromosome_arms <- readRDS("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/data/hg19_chromosome_arms.rds")
coding_genes <- readRDS("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/data/hg19_expressing_genes.rds")


gene_level_df <- load_and_prepare_infercnv_reference(objs, k = 1.5)

metadata$cell_type <- "TE"  

collapse_df <- collapse_genes_to_cnv_segments(gene_level_df)

collapse_path <- file.path(OUT_DIR,"collapsed.rds")
saveRDS(collapse_df, file = collapse_path)

expression_df <- obj@gene_order

# ── Build GRanges ─────────────────────────────────────────────────────────────
seg_gr <- GenomicRanges::GRanges(
  seqnames = collapse_df$chr,
  ranges   = IRanges::IRanges(
    start = collapse_df$start,
    end   = collapse_df$stop
  )
)
seg_gr$seg_idx <- seq_len(nrow(collapse_df))

coding_gr <- GenomicRanges::GRanges(
  seqnames = coding_genes$chr,
  ranges   = IRanges::IRanges(
    start = coding_genes$start,
    end   = coding_genes$end
  )
)

coding_gr$gene <- coding_genes$gene_name


cat("coding_genes columns:", colnames(coding_genes), "\n")


expressed_table <- as.data.frame(objs$TE[[2]])
expressed_table$gene <- rownames(expressed_table)
colnames(expressed_table) <- c("chr", "start", "stop", "gene")

expressed_gr <- GenomicRanges::GRanges(
  seqnames = expressed_table$chr,
  ranges   = IRanges::IRanges(
    start = expressed_table$start,
    end   = expressed_table$stop
  )
)

# ── Gene sets ─────────────────────────────────────────────────────────────────
coding_gene_set    <- unique(coding_genes$gene_name)
expressed_gene_set <- unique(expressed_table$gene)

# Overlap set: coding AND expressed
coding_expressed_set <- intersect(
  coding_gene_set,
  expressed_gene_set
)

cat("Coding genes (GTF):      ", length(coding_gene_set), "\n")
cat("Expressed (@gene_order): ", length(expressed_gene_set), "\n")
cat("Coding + expressed:      ", length(coding_expressed_set), "\n")

# ── Count by coordinates (GRanges) ────────────────────────────────────────────
cat("Counting coding genes per segment...\n")
hits_coding <- GenomicRanges::findOverlaps(
  seg_gr, coding_gr, type = "any")


cat("Counting expressed genes per segment...\n")
hits_expressed <- GenomicRanges::findOverlaps(
  seg_gr, expressed_gr, type = "any")

expressed_counts <- data.frame(
  seg_idx = S4Vectors::queryHits(hits_expressed)
) %>%
  dplyr::count(seg_idx,
               name = "n_expressing_genes")



segments <- filter_segments_by_gene_density(
    collapse_df          = collapse_df,
    gene_order           = expressed_table,
    coding_gr            = coding_gr,
    coding_expressed_set = coding_expressed_set,
    pct_max              = 45,
    pct_floor            = 30,
    min_expr_density     = 0.5,
    min_coding_density   = 0.5
  )


merged <- merge_and_density_check(df = segments, 
    gene_order = expressed_table,
    coding_gr = coding_gr,
    coding_expressed_set = coding_expressed_set,
    max_gap_mb         = 10,
    pct_max            = 45,
    pct_floor          = 30,
    min_expr_density   = 0.5,
    min_coding_density = 0.5
)

merged <- merged %>%
  mutate(
    cnv_length = end - start,
    cnv_length_mb = cnv_length/1e6
  ) %>%
  filter(cnv_length_mb >=15)
  

print(colnames(merged))

cnv_annotated <- add_chromosome_info(
  merged,
  chromosome_arms,
  chr_col   = "chr",
  start_col = "start",
  end_col   = "end"
)

#ref_cells <- metadata %>%
#  filter(split_group == "ref")

#metadata <- metadata %>%
#  filter(split_group != "ref")

cell_sizes <- compute_cell_sizes(
  metadata   = metadata,
  group_cols = c("cell_type","embryo"),
  cell_col   = "cell_name"
)


cnv_annotated$cell_type <- "TE"


clustered_events <- run_cnv_locus_analysis(
  cnv_annotated,
  by             = "embryo",
  overlap_method = "reciprocal",
  min_ovelap     = 0.75,
  sample_col     = "reference",
  cell_col       = "cell_name"
)

out_file <- file.path(OUT_DIR,"clustered_events.rds")

  
clustered_events$clustered_events$cell_type <- "TE"
clustered_events$cnv_locus_summary$cell_type <- "TE"


thresholded_df <- prepare_cnv_thresholds(
            summary_df           = clustered_events$cnv_locus_summary,
            clustered_events     = clustered_events$clustered_events,
            by_union             = "embryo",
            cell_sizes           = cell_sizes,
            k                    = 1.5,
            sensitivity_floor_mb = 20,
            min_required_cells   = 3,
            round_fun            = ceiling)

cat("      thresholded rows:", nrow(thresholded_df), "\n")
          

thresholded_df <- add_arm_percentages(cnv_df = thresholded_df, chromosome_arms = chromosome_arms)
          

filtered_df <- filter_cnv_loci(
            clustered_events     = thresholded_df,
            p_arm_permission     = 65,
            q_arm_permission     = 65,
            whole_chr_permission = 60
          )
          
deduped_df <- filtered_df %>%
            dplyr::distinct(
              cell_name, chr, cnv_state,
              start, end,
              .keep_all = TRUE
            )

          
final_score <- post_score_merge(
            scored_events   = deduped_df,
            cell_col        = "cell_name",
            min_overlap     = 0.75,
            overlap_method  = "reciprocal"
          )
          




out_file <- file.path(OUT_DIR,"scored.rds")
        
saveRDS(final_score, out_file)
          
          

