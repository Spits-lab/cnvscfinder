# =============================================================================
# iteration.R
# Domain logic for the iterative euploid-reference purification loop: state
# I/O, the reference-update rule, the annotation builder, and the CNV-calling
# steps run against one pooled inferCNV object. The scripts/iter_*.R entry
# points do the I/O and call these functions.
#
# Config loading/validation is NOT here - see R/config.R (load_config(),
# validate_config()/assert_valid_config()). This file only depends on that
# one (cfg is passed in as plain data); it never sources it.
#
# One folder per iteration: <iter_root>/iter_<n>/
#   state.rds       list(iter, ref_ids, aneuploid_ids, scored_events)
#   infercnv/       inferCNV output (run.final.infercnv_obj)
#   PLATEAU         marker, written when no reference cell was called aneuploid
# =============================================================================

# Fall back to a default when NULL (base R only has this from 4.4)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Logging ──────────────────────────────────────────────────────────────────

#' Emit one greppable metric line: `[METRIC] key=value key=value ...`
#'
#' To strip metrics later: delete this function and its calls, or
#' `grep -v '^\[METRIC\]'` on the logs.
log_metric <- function(...) {
  vals <- list(...)
  if (length(vals) == 0L || is.null(names(vals)) || any(!nzchar(names(vals)))) {
    stop("log_metric() needs named arguments")
  }
  kv <- vapply(names(vals), function(k) {
    v <- vals[[k]]
    paste0(k, "=", if (is.numeric(v)) format(v, digits = 6, scientific = FALSE, trim = TRUE) else as.character(v))
  }, character(1))
  cat("[METRIC] ", paste(kv, collapse = " "), "\n", sep = "")
  invisible(NULL)
}

# ── Iteration folders ────────────────────────────────────────────────────────

iter_dir_of <- function(cfg, iter) {
  if (is.null(cfg$iter_root)) stop("iter_root missing in config")
  file.path(cfg$iter_root, sprintf("iter_%d", iter))
}

# ── Iteration state (one object per iteration) ───────────────────────────────

read_iter_state <- function(iter_dir) {
  path <- file.path(iter_dir, "state.rds")
  if (!file.exists(path)) stop("no state.rds in ", iter_dir)
  readRDS(path)
}

write_iter_state <- function(iter_dir, state) {
  dir.create(iter_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(state, file.path(iter_dir, "state.rds"))
  invisible(state)
}

# ── Reference bookkeeping ────────────────────────────────────────────────────

#' Reference set for the next iteration
#'
#' @param prev_ref_ids character, reference of the previous iteration
#' @param aneuploid_ids character, cells with an event in the previous iteration
#' @param all_cell_ids  character, every cell in the metadata
#' @param rule "remove"    -> previous reference minus aneuploid cells (can only shrink)
#'             "recompute" -> all cells minus aneuploid cells (cells may re-enter)
next_reference_ids <- function(prev_ref_ids, aneuploid_ids, all_cell_ids,
                               rule = c("remove", "recompute")) {
  rule <- match.arg(rule)
  switch(rule,
    remove    = setdiff(prev_ref_ids, aneuploid_ids),
    recompute = setdiff(all_cell_ids, aneuploid_ids)
  )
}

#' inferCNV annotation: reference cells -> `refs`, all others -> `non_refs`
#'
#' @return data.frame, rownames = cell ids, one column `split_group`
make_annotation <- function(cell_ids, ref_ids, min_ref = 50L) {
  cell_ids <- as.character(cell_ids)
  if (anyDuplicated(cell_ids)) stop("duplicated cell ids")
  unknown <- setdiff(ref_ids, cell_ids)
  if (length(unknown) > 0L) stop(length(unknown), " reference cell(s) not in the cell set")

  is_ref <- cell_ids %in% ref_ids
  if (sum(is_ref) < min_ref) stop("reference has ", sum(is_ref), " cells (< min_ref = ", min_ref, ")")
  if (all(is_ref))           stop("no observation cells left outside the reference")

  data.frame(split_group = ifelse(is_ref, "refs", "non_refs"),
             row.names = cell_ids, stringsAsFactors = FALSE)
}

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
                             min_overlap = 0.75, overlap_method = "reciprocal") {
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

#' Call CNVs from one pooled inferCNV object (the steps of R/round2.R)
#'
#' All thresholds come from `cfg` (the YAML config), so each can be overridden.
#'
#' @param obj inferCNV object (run.final.infercnv_obj)
#' @param metadata data.frame with cfg$cell_col
#' @param cfg config list
#' @param coding_genes data.frame (chr, start, end, gene_name)
#' @param chromosome_arms chromosome arm table
#' @return data.frame of scored events (one row per cell and event)
call_iteration_cnvs <- function(obj, metadata, cfg, coding_genes, chromosome_arms) {
  grp_col  <- cfg$cell_group_cluster
  grp_val  <- cfg$group_value
  cell_col <- cfg$cell_col

  gene_level_df <- load_and_prepare_infercnv_reference(
    list(refs = list(obj@expr.data, obj@gene_order)), k = cfg$k_discrete)
  collapse_df <- collapse_genes_to_cnv_segments(gene_level_df)

  expressed_table <- as.data.frame(obj@gene_order)
  expressed_table$gene <- rownames(expressed_table)
  colnames(expressed_table) <- c("chr", "start", "stop", "gene")

  coding_gr <- GenomicRanges::GRanges(
    seqnames = coding_genes$chr,
    ranges   = IRanges::IRanges(start = coding_genes$start, end = coding_genes$end))
  coding_gr$gene <- coding_genes$gene_name
  coding_expressed_set <- intersect(unique(coding_genes$gene_name),
                                    unique(expressed_table$gene))

  segments <- filter_segments_by_gene_density(
    collapse_df = collapse_df, gene_order = expressed_table, coding_gr = coding_gr,
    coding_expressed_set = coding_expressed_set,
    pct_max = cfg$pct_max, pct_floor = cfg$pct_floor,
    min_expr_density = cfg$min_expr_density, min_coding_density = cfg$min_coding_density)

  merged <- merge_and_density_check(
    df = segments, gene_order = expressed_table, coding_gr = coding_gr,
    coding_expressed_set = coding_expressed_set, max_gap_mb = cfg$max_gap_mb,
    pct_max = cfg$pct_max, pct_floor = cfg$pct_floor,
    min_expr_density = cfg$min_expr_density, min_coding_density = cfg$min_coding_density)

  merged <- merged %>%
    dplyr::mutate(cnv_length = end - start, cnv_length_mb = cnv_length / 1e6) %>%
    dplyr::filter(cnv_length_mb >= cfg$min_segment_mb)

  cnv_annotated <- add_chromosome_info(merged, chromosome_arms,
                                       chr_col = "chr", start_col = "start", end_col = "end")

  # The whole sample is analysed as one group (as in round2.R)
  metadata[[grp_col]]      <- grp_val
  cnv_annotated[[grp_col]] <- grp_val
  cell_sizes <- compute_cell_sizes(metadata, group_cols = grp_col, cell_col = cell_col)

  clustered <- run_cnv_locus_analysis(
    cnv_annotated, overlap_method = cfg$overlap_method, min_overlap = cfg$min_overlap,
    range = cfg$range, sensitivity_floor_mb = cfg$sensitivity_floor_mb, max_mb = cfg$max_mb,
    sample_col = "reference", cell_col = cell_col)
  clustered$clustered_events[[grp_col]]   <- grp_val
  clustered$cnv_locus_summary[[grp_col]]  <- grp_val

  thresholded <- prepare_cnv_thresholds(
    summary_df = clustered$cnv_locus_summary, clustered_events = clustered$clustered_events,
    by_union = grp_col, cell_sizes = cell_sizes, k = cfg$k_threshold,
    sensitivity_floor_mb = cfg$sensitivity_floor_mb,
    min_required_cells = cfg$min_required_cells, round_fun = ceiling)
  thresholded <- add_arm_percentages(cnv_df = thresholded, chromosome_arms = chromosome_arms)

  filtered <- filter_cnv_loci(
    clustered_events = thresholded, p_arm_permission = cfg$p_arm_permission,
    q_arm_permission = cfg$q_arm_permission, whole_chr_permission = cfg$whole_chr_permission)
  filtered <- dplyr::distinct(filtered, .data[[cell_col]], chr, cnv_state, start, end,
                              .keep_all = TRUE)

  post_score_merge(filtered, cell_col = cell_col, min_overlap = cfg$post_min_overlap,
                   overlap_method = cfg$post_overlap_method)
}
