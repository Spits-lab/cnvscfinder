# =============================================================================
# config.R
# User parameter processing: load a YAML config (with --override merging),
# and validate it. One home for every config-checking concern in this repo:
#   - the iteration/purification loop's config      -> validate_config()
#   - the grid-search benchmark's config             -> validate_grid_search_config()
#   - the shared comparison primitive both build on  -> in_range()/range_desc()
# Also reused (primitives only) by scripts/run_pipeline.R's own check_range(),
# which keeps its own fail-fast style and message wording.
#
# What does NOT belong here: domain logic (CNV calling, clustering, iteration
# bookkeeping) - that stays in R/iteration.R and friends. This file is only
# about reading and checking *parameters*.
# =============================================================================

# Fall back to a default when NULL (base R only has this from 4.4)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Primitives ───────────────────────────────────────────────────────────────

#' Is `val` within [lower, upper] (bounds inclusive/exclusive as given)?
#' Vectorised over `val`.
in_range <- function(val, lower, upper, lower_inclusive = TRUE, upper_inclusive = TRUE) {
  ok_lower <- if (lower_inclusive) val >= lower else val > lower
  ok_upper <- if (upper_inclusive) val <= upper else val < upper
  ok_lower & ok_upper
}

#' Human-readable interval, e.g. "[0, 1]" or "(0, 100]"
range_desc <- function(lower, upper, lower_inclusive = TRUE, upper_inclusive = TRUE) {
  paste0(if (lower_inclusive) "[" else "(", lower, ", ", upper,
         if (upper_inclusive) "]" else ")")
}

# ── Loading ──────────────────────────────────────────────────────────────────

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

# ── Iteration-loop config ────────────────────────────────────────────────────

#' Check the iteration/purification-loop config before any job is submitted
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
    if (!in_range(v, lo, hi, lower_inclusive = !lo_open)) {
      add_err(k, " = ", v, " outside the allowed range ", range_desc(lo, hi, lower_inclusive = !lo_open))
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
  # (path built inline, matching iter_dir_of() in R/iteration.R, to keep this
  # file free of any dependency on that one - see file header)
  if (!is.null(iter) && iter > 1L && has("iter_root")) {
    prev_dir  <- file.path(cfg$iter_root, sprintf("iter_%d", iter - 1L))
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

# ── Grid-search config ───────────────────────────────────────────────────────

#' Check a grid-search config before any combo runs
#'
#' Collects every problem instead of stopping at the first. Cheap: never loads
#' the inferCNV object or metadata beyond checking they exist.
#'
#' @param cfg config list from load_config()
#' @param overlap_ok function(method) -> TRUE/FALSE, or NULL to skip the method check
#' @return list(errors = character, warnings = character)
validate_grid_search_config <- function(cfg, overlap_ok = NULL) {
  err <- character(0); wrn <- character(0)
  add_err <- function(...) err <<- c(err, paste0(...))
  add_wrn <- function(...) wrn <<- c(wrn, paste0(...))

  known <- c(
    "infercnv_obj_path", "metadata_path", "chromosome_arms_path",
    "coding_genes_path", "out_dir", "cell_col", "group_value",
    "cell_group_cluster", "k_dis_values", "k_fre_values", "sens_floor_values",
    "ovlp_values", "pct_floor_values", "min_expr_density_values", "pct_max",
    "min_coding_density", "max_gap_mb", "min_segment_mb", "p_arm_permission",
    "q_arm_permission", "whole_chr_permission", "min_required_cells",
    "overlap_method", "range", "post_min_overlap", "post_overlap_method")

  # ── Keys ───────────────────────────────────────────────────────────────────
  unknown <- setdiff(names(cfg), known)
  if (length(unknown) > 0L) add_err("unknown key(s), typo? ", paste(unknown, collapse = ", "))
  missing <- setdiff(known, names(cfg))
  if (length(missing) > 0L) add_err("missing key(s): ", paste(missing, collapse = ", "))
  nulls <- setdiff(intersect(known, names(cfg)), names(Filter(Negate(is.null), cfg)))
  if (length(nulls) > 0L) add_err("key(s) set to null: ", paste(nulls, collapse = ", "))

  has <- function(k) !is.null(cfg[[k]])
  num_vec <- function(k, lo = -Inf, hi = Inf, lo_open = FALSE) {
    if (!has(k)) return(invisible())
    v <- cfg[[k]]
    if (!is.numeric(v) || length(v) < 1L || anyNA(v)) {
      return(add_err(k, " must be one or more numbers, got: ", paste(v, collapse = ", ")))
    }
    bad <- !in_range(v, lo, hi, lower_inclusive = !lo_open)
    if (any(bad)) {
      add_err(k, " has value(s) outside the allowed range ",
              range_desc(lo, hi, lower_inclusive = !lo_open), ": ",
              paste(v[bad], collapse = ", "))
    }
  }
  num1 <- function(k, lo = -Inf, hi = Inf, lo_open = FALSE) {
    if (!has(k)) return(invisible())
    v <- cfg[[k]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      return(add_err(k, " must be a single number, got: ", format(v)))
    }
    if (!in_range(v, lo, hi, lower_inclusive = !lo_open)) {
      add_err(k, " = ", v, " outside the allowed range ", range_desc(lo, hi, lower_inclusive = !lo_open))
    }
  }

  # ── Grid axes (each one or more values) ───────────────────────────────────
  for (k in c("k_dis_values", "k_fre_values", "sens_floor_values")) num_vec(k, 0, lo_open = TRUE)
  for (k in c("ovlp_values", "post_min_overlap")) num_vec(k, 0, 1, lo_open = TRUE)
  for (k in c("pct_floor_values", "min_expr_density_values")) num_vec(k, 0, lo_open = TRUE)

  # ── Fixed knobs ────────────────────────────────────────────────────────────
  num1("pct_max", 0, 100, lo_open = TRUE)
  num1("min_coding_density", 0, lo_open = TRUE)
  num1("max_gap_mb", 0)
  num1("min_segment_mb", 0)
  for (k in c("p_arm_permission", "q_arm_permission", "whole_chr_permission")) num1(k, 0, 100, lo_open = TRUE)
  num1("min_required_cells", 1)
  num1("range", 0, 1, lo_open = TRUE)

  if (has("pct_floor_values") && has("pct_max") && is.numeric(cfg$pct_floor_values) &&
      is.numeric(cfg$pct_max) && any(cfg$pct_floor_values > cfg$pct_max)) {
    add_err("pct_floor_values has value(s) above pct_max (", cfg$pct_max, ")")
  }

  for (k in c("cell_col", "group_value", "cell_group_cluster")) {
    if (has(k) && (!is.character(cfg[[k]]) || length(cfg[[k]]) != 1L || !nzchar(cfg[[k]]))) {
      add_err(k, " must be a non-empty string")
    }
  }
  if (!is.null(overlap_ok)) {
    for (k in c("overlap_method", "post_overlap_method")) {
      if (has(k) && !isTRUE(overlap_ok(cfg[[k]]))) add_err(k, " '", cfg[[k]], "' is not a registered overlap method")
    }
  } else {
    for (k in c("overlap_method", "post_overlap_method")) {
      if (has(k) && (!is.character(cfg[[k]]) || length(cfg[[k]]) != 1L || !nzchar(cfg[[k]]))) {
        add_err(k, " must be a non-empty string")
      }
    }
  }

  # ── Files ──────────────────────────────────────────────────────────────────
  for (k in c("infercnv_obj_path", "metadata_path", "chromosome_arms_path", "coding_genes_path")) {
    if (has(k) && !file.exists(cfg[[k]])) add_err(k, ": file not found: ", cfg[[k]])
  }

  # ── Grid size sanity ───────────────────────────────────────────────────────
  if (all(vapply(c("k_dis_values", "k_fre_values", "sens_floor_values", "ovlp_values",
                    "pct_floor_values", "min_expr_density_values"), has, logical(1)))) {
    n_total <- length(cfg$k_dis_values) * length(cfg$k_fre_values) *
      length(cfg$sens_floor_values) * length(cfg$ovlp_values) *
      length(cfg$pct_floor_values) * length(cfg$min_expr_density_values)
    if (n_total > 200L) {
      add_wrn("grid has ", n_total, " combinations — confirm this is meant to be that large")
    }
  }

  list(errors = err, warnings = wrn)
}

#' Stop with every problem listed; print warnings
assert_valid_grid_search_config <- function(cfg, ...) {
  res <- validate_grid_search_config(cfg, ...)
  for (w in res$warnings) message("WARNING: ", w)
  if (length(res$errors) > 0L) {
    stop(length(res$errors), " config problem(s):\n",
         paste0("  ", seq_along(res$errors), ". ", res$errors, collapse = "\n"), call. = FALSE)
  }
  invisible(cfg)
}
