# =============================================================================
# iteration.R
# Helpers for the iterative euploid-reference purification loop.
# The scripts/iter_*.R entry points do the I/O and call these functions.
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

# ── Config ───────────────────────────────────────────────────────────────────

#' Read a YAML config; `override` is a YAML fragment, e.g. "{min_overlap: 0.65}"
load_config <- function(config_path, override = NULL) {
  if (is.null(config_path) || !file.exists(config_path)) {
    stop("config file not found: ", config_path)
  }
  cfg <- yaml::read_yaml(config_path)
  if (!is.null(override) && nzchar(override)) {
    cfg <- utils::modifyList(cfg, yaml::yaml.load(override))
  }
  cfg
}

#' Check a config before any job is submitted
#'
#' Collects every problem instead of stopping at the first. Cheap on purpose: it
#' never loads the counts matrix or previous pipeline results (those checks
#' happen in the scripts that load them anyway).
#'
#' @param cfg config list from load_config()
#' @param fresh_root TRUE when launching a new run: iter_root must not hold earlier results
#' @param iter iteration about to run; start_ref/start_results only matter for 1 (NULL = launch check)
#' @param overlap_ok function(method) -> TRUE/FALSE, or NULL to skip the method check
#' @return list(errors = character, warnings = character)
validate_config <- function(cfg, fresh_root = TRUE, iter = NULL, overlap_ok = NULL) {
  err <- character(0); wrn <- character(0)
  add_err <- function(...) err <<- c(err, paste0(...))
  add_wrn <- function(...) wrn <<- c(wrn, paste0(...))

  known <- c(
    "n_iter", "iter_root", "ref_rule", "min_ref", "start_ref", "start_results",
    "time_limit", "mem", "counts_path", "metadata_path", "gene_order",
    "chromosome_arms_path", "coding_genes_path", "cell_col", "cell_group_cluster",
    "group_value", "cutoff", "k_discrete", "pct_max", "pct_floor", "min_expr_density",
    "min_coding_density", "max_gap_mb", "min_segment_mb", "overlap_method",
    "min_overlap", "range", "max_mb", "sensitivity_floor_mb", "k_threshold",
    "min_required_cells", "p_arm_permission", "q_arm_permission",
    "whole_chr_permission", "post_min_overlap", "post_overlap_method")
  optional <- c("start_ref", "start_results")

  # ── Keys ───────────────────────────────────────────────────────────────────
  unknown <- setdiff(names(cfg), known)
  if (length(unknown) > 0L) add_err("unknown key(s), typo? ", paste(unknown, collapse = ", "))
  missing <- setdiff(setdiff(known, optional), names(cfg))
  if (length(missing) > 0L) add_err("missing key(s): ", paste(missing, collapse = ", "))
  nulls <- setdiff(intersect(setdiff(known, optional), names(cfg)),
                   names(Filter(Negate(is.null), cfg)))
  if (length(nulls) > 0L) add_err("key(s) set to null: ", paste(nulls, collapse = ", "))

  has <- function(k) !is.null(cfg[[k]])
  num <- function(k, lo = -Inf, hi = Inf, lo_open = FALSE, integer = FALSE) {
    if (!has(k)) return(invisible())
    v <- cfg[[k]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) return(add_err(k, " must be a single number, got: ", format(v)))
    if (integer && v != round(v)) add_err(k, " must be a whole number, got ", v)
    if (v < lo || v > hi || (lo_open && v <= lo)) {
      add_err(k, " = ", v, " outside the allowed range ",
              if (lo_open) "(" else "[", lo, ", ", hi, "]")
    }
  }

  # ── Numbers ────────────────────────────────────────────────────────────────
  num("n_iter", 1, 50, integer = TRUE)
  num("min_ref", 1, integer = TRUE)
  num("cutoff", 0, lo_open = FALSE)
  for (k in c("min_overlap", "post_min_overlap", "range")) num(k, 0, 1, lo_open = TRUE)
  for (k in c("p_arm_permission", "q_arm_permission", "whole_chr_permission",
              "pct_max", "pct_floor")) num(k, 0, 100, lo_open = TRUE)
  for (k in c("k_discrete", "k_threshold", "sensitivity_floor_mb", "max_mb",
              "min_expr_density", "min_coding_density")) num(k, 0, lo_open = TRUE)
  num("min_required_cells", 1, integer = TRUE)
  num("max_gap_mb", 0)
  num("min_segment_mb", 0)
  if (is.numeric(cfg$pct_floor) && is.numeric(cfg$pct_max) && cfg$pct_floor > cfg$pct_max) {
    add_err("pct_floor (", cfg$pct_floor, ") must be <= pct_max (", cfg$pct_max, ")")
  }
  if (is.numeric(cfg$min_segment_mb) && is.numeric(cfg$max_mb) && cfg$min_segment_mb > cfg$max_mb) {
    add_err("min_segment_mb (", cfg$min_segment_mb, ") must be <= max_mb (", cfg$max_mb, ")")
  }
  if (is.numeric(cfg$min_segment_mb) && is.numeric(cfg$sensitivity_floor_mb) &&
      cfg$min_segment_mb > cfg$sensitivity_floor_mb) {
    add_wrn("min_segment_mb (", cfg$min_segment_mb, ") is above sensitivity_floor_mb (",
            cfg$sensitivity_floor_mb, "): segments the floor would accept are dropped first")
  }

  # ── Choices and formats ────────────────────────────────────────────────────
  if (has("ref_rule") && !cfg$ref_rule %in% c("remove", "recompute")) {
    add_err("ref_rule must be 'remove' or 'recompute', got '", cfg$ref_rule, "'")
  }
  if (!is.null(overlap_ok)) {
    for (k in c("overlap_method", "post_overlap_method")) {
      if (has(k) && !isTRUE(overlap_ok(cfg[[k]]))) add_err(k, " '", cfg[[k]], "' is not a registered overlap method")
    }
  }
  if (has("time_limit") && !grepl("^([0-9]+-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}$", cfg$time_limit)) {
    add_err("time_limit must look like HH:MM:SS, got '", cfg$time_limit, "'")
  }
  if (has("mem") && !grepl("^[0-9]+[MGT]$", cfg$mem)) add_err("mem must look like 128G, got '", cfg$mem, "'")
  for (k in c("cell_col", "cell_group_cluster", "group_value")) {
    if (has(k) && (!is.character(cfg[[k]]) || length(cfg[[k]]) != 1L || !nzchar(cfg[[k]]))) {
      add_err(k, " must be a non-empty string")
    }
  }

  # ── Files ──────────────────────────────────────────────────────────────────
  for (k in c("counts_path", "metadata_path", "gene_order", "chromosome_arms_path", "coding_genes_path")) {
    if (has(k) && !file.exists(cfg[[k]])) add_err(k, ": file not found: ", cfg[[k]])
  }
  if (fresh_root && has("iter_root") && dir.exists(cfg$iter_root) &&
      length(list.files(cfg$iter_root, all.files = TRUE, no.. = TRUE)) > 0L) {
    add_err("iter_root is not empty (would mix with an earlier run): ", cfg$iter_root)
  }

  # ── Starting reference (iteration 1 / launch only) ─────────────────────────
  if (is.null(iter) || iter == 1L) {
    if (has("start_ref") == has("start_results")) {
      add_err("set exactly one of start_ref / start_results (",
              if (has("start_ref")) "both are set" else "neither is set", ")")
    }
    if (has("start_results") && !file.exists(cfg$start_results)) {
      add_err("start_results: file not found: ", cfg$start_results)
    }
    if (has("start_ref")) {
      if (!file.exists(cfg$start_ref)) {
        add_err("start_ref: file not found: ", cfg$start_ref)
      } else {
        ids <- tryCatch(readRDS(cfg$start_ref), error = function(e) NULL)
        if (!is.character(ids) || is.data.frame(ids)) {
          add_err("start_ref must be an rds character vector of cell ids, not ",
                  paste(class(ids), collapse = "/"))
        } else if (has("metadata_path") && file.exists(cfg$metadata_path) && has("cell_col")) {
          md <- readRDS(cfg$metadata_path)
          if (!cfg$cell_col %in% colnames(md)) {
            add_err("cell_col '", cfg$cell_col, "' is not a column of the metadata")
          } else if (length(setdiff(ids, as.character(md[[cfg$cell_col]]))) > 0L) {
            add_err(length(setdiff(ids, as.character(md[[cfg$cell_col]]))),
                    " start_ref cell id(s) are not in the metadata")
          } else if (length(ids) < 0.05 * nrow(md)) {
            add_wrn("start_ref holds only ", length(ids), " of ", nrow(md), " cells")
          }
        }
      }
    }
  }

  # ── Metadata columns (cheap: metadata only) ────────────────────────────────
  if (has("metadata_path") && file.exists(cfg$metadata_path) && has("cell_col")) {
    md <- tryCatch(readRDS(cfg$metadata_path), error = function(e) NULL)
    if (is.null(md) || !is.data.frame(md)) {
      add_err("metadata_path must hold a data.frame")
    } else if (!cfg$cell_col %in% colnames(md)) {
      if (!any(grepl("is not a column of the metadata", err))) {
        add_err("cell_col '", cfg$cell_col, "' is not a column of the metadata")
      }
    } else if (anyDuplicated(md[[cfg$cell_col]])) {
      add_err("metadata has duplicated values in '", cfg$cell_col, "'")
    }
  }

  # ── Resuming (start-iter > 1): the previous iteration must have finished ──
  if (!is.null(iter) && iter > 1L && has("iter_root")) {
    prev_dir  <- iter_dir_of(cfg, iter - 1L)
    prev_path <- file.path(prev_dir, "state.rds")
    if (!file.exists(prev_path)) {
      add_err("resuming at iteration ", iter, " needs iteration ", iter - 1L,
              "'s state.rds, not found in ", prev_dir)
    } else {
      prev <- tryCatch(readRDS(prev_path), error = function(e) NULL)
      if (is.null(prev$ref_ids) || is.null(prev$aneuploid_ids)) {
        add_err("iteration ", iter - 1L, "'s state.rds is incomplete (missing ref_ids or ",
                "aneuploid_ids) — it did not finish the process step")
      }
    }
  }

  list(errors = err, warnings = wrn)
}

#' Stop with every problem listed; print warnings
assert_valid_config <- function(cfg, ...) {
  res <- validate_config(cfg, ...)
  for (w in res$warnings) message("WARNING: ", w)
  if (length(res$errors) > 0L) {
    stop(length(res$errors), " config problem(s):\n",
         paste0("  ", seq_along(res$errors), ". ", res$errors, collapse = "\n"), call. = FALSE)
  }
  invisible(cfg)
}

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
