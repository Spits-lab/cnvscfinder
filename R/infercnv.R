#' Validate metadata for inferCNV pipeline
#'
#' @param metadata     data.frame with at least cell_name and cell_type_col
#' @param counts_mx    raw counts matrix (genes x cells)
#' @param cell_type_col string, name of the column containing cell type labels
#'
#' @return invisibly returns TRUE if all checks pass; stops with message if not
validate_metadata <- function(metadata,
                              counts_mx,
                              cell_type_col,
                              min_cells = 90) {
  
  # ── 1. Required columns present ───────────────────────────────────────────
  required_cols <- c("cell_name", cell_type_col)
  missing_cols  <- setdiff(required_cols,
                           colnames(metadata))
  
  if (length(missing_cols) > 0) {
    stop(
      "metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nRequired: 'cell_name' and '",
      cell_type_col, "'"
    )
  }
  
  # ── 2. cell_name must be unique ───────────────────────────────────────────
  if (any(duplicated(metadata$cell_name))) {
    n_dup <- sum(duplicated(metadata$cell_name))
    stop(
      n_dup, " duplicated cell_name(s) found. ",
      "Each cell must appear exactly once."
    )
  }
  
  # ── 3. cell_name must match colnames of counts_mx ─────────────────────────
  cells_in_meta   <- metadata$cell_name
  cells_in_counts <- colnames(counts_mx)
  
  only_in_meta   <- setdiff(cells_in_meta,
                            cells_in_counts)
  only_in_counts <- setdiff(cells_in_counts,
                            cells_in_meta)
  
  if (length(only_in_meta) > 0) {
    stop(
      length(only_in_meta),
      " cell(s) in metadata not found in counts_mx. ",
      "First few: ",
      paste(head(only_in_meta, 5), collapse = ", ")
    )
  }
  
  if (length(only_in_counts) > 0) {
    warning(
      length(only_in_counts),
      " cell(s) in counts_mx not in metadata ",
      "— will be dropped. First few: ",
      paste(head(only_in_counts, 5), collapse = ", ")
    )
  }
  
  # ── 4. cell_type_col must not have NA ─────────────────────────────────────
  n_na <- sum(is.na(metadata[[cell_type_col]]))
  if (n_na > 0) {
    stop(
      n_na, " NA value(s) in '", cell_type_col, "'. ",
      "All cells must have a cell type label."
    )
  }
  
  # ── 5. Report cell type sizes and filter ──────────────────────────────────
  type_counts <- table(metadata[[cell_type_col]])
  
  # Identify cell types below threshold
  below_min <- names(type_counts)[
    type_counts < min_cells]
  above_min <- names(type_counts)[
    type_counts >= min_cells]
  
  message("Cell type composition:")
  for (ct in names(type_counts)) {
    n    <- type_counts[[ct]]
    flag <- if (n < min_cells)
      sprintf(" [REMOVED: %d < %d minimum]",
              n, min_cells)
    else ""
    message("  ", ct, ": ", n, " cells", flag)
  }
  
  # Remove cell types below threshold
  if (length(below_min) > 0) {
    message(sprintf(
      "\nRemoving %d cell type(s) with < %d cells: %s",
      length(below_min),
      min_cells,
      paste(below_min, collapse = ", ")
    ))
    
    metadata <- metadata %>%
      dplyr::filter(
        !(!!sym(cell_type_col) %in% below_min)
      )
    
    message(sprintf(
      "Remaining: %d cells across %d cell type(s)",
      nrow(metadata),
      length(above_min)
    ))
  }
  
  # ── 6. At least 2 cell types remaining ───────────────────────────────────
  n_types <- length(unique(
    metadata[[cell_type_col]]))
  
  if (n_types < 2) {
    warning(
      "Only ", n_types,
      " cell type(s) remaining after filtering. ",
      "Across-cell-type comparisons not possible."
    )
  }
  
  if (n_types == 0) {
    stop("No cell types remaining after filtering. ",
         "Lower min_cells threshold.")
  }
  
  # Return filtered metadata
  return(metadata)
}


#' Add A/B/C random split column to metadata for a single cell type
#'
#' Splits cells of one cell type into three roughly equal groups.
#'
#' @param metadata      data.frame (full metadata, not pre-subsetted)
#' @param cell_type_col string, column name for cell type
#' @param cell_type_val string, which cell type to split
#'
#' @return data.frame subset for that cell type with added 'split_group' column
make_splits <- function(metadata, cell_type_col, cell_type_val, n_splits = 3) {
  sub <- metadata[metadata[[cell_type_col]] == cell_type_val, , drop = FALSE]
  
  if (!is.numeric(n_splits) || n_splits < 2 || n_splits != round(n_splits)) {
    stop("n_splits must be an integer >= 2. Got: ", n_splits)
  }
  n_splits <- as.integer(n_splits)
  
  if (n_splits > 26) {
    stop("n_splits cannot exceed 26 (only 26 capital letters available). ",
         "Got: ", n_splits)
  }
  
  sub <- metadata[metadata[[cell_type_col]] == cell_type_val, , drop = FALSE]
  n   <- nrow(sub)
  
  if (n < n_splits) {
    stop(
      "Cell type '", cell_type_val, "' has only ", n, " cell(s) but ",
      "n_splits = ", n_splits, ". Need at least ", n_splits, " cells."
    )
  }
  
  # ── Generate letter labels ─────────────────────────────────────────────────
  labels <- LETTERS[seq_len(n_splits)]  
  
  # ── Compute group sizes ───────────────────────────────────────────────────
  # Each of the first (n_splits - 1) groups gets floor(n / n_splits) cells.
  # The last group gets whatever remains so the total always equals n.
  base_size  <- floor(n / n_splits)
  group_sizes <- rep(base_size, n_splits)
  group_sizes[n_splits] <- n - base_size * (n_splits - 1L)
  
  # ── Build label vector and shuffle ────────────────────────────────────────
  split_labels <- c()
  for (i in seq_len(n_splits)) {
    split_labels <- c(split_labels, rep(labels[i], group_sizes[i]))
  }
  
  sub$split_group <- split_labels
  
  size_summary <- paste(
    mapply(function(lbl, sz) sprintf("%s=%d", lbl, sz), labels, group_sizes),
    collapse = " | "
  )
  message(sprintf("  Split '%s' (n=%d, %d groups): %s",
                  cell_type_val, n, n_splits, size_summary))
  return(sub)
}


#' Build annotations data.frame for inferCNV
#'
#' inferCNV requires a single-column data.frame where:
#'   - rownames = cell barcodes
#'   - single column = group label
#'
#' @param cell_names character vector of cell barcodes
#' @param group_labels character vector of group labels (same length)
#'
#' @return data.frame suitable for inferCNV annotations_file argument
build_annotations_df <- function(cell_names, group_labels) {
  
  if (length(cell_names) != length(group_labels)) {
    stop("cell_names and group_labels must have the same length.")
  }
  
  annot <- data.frame(group = group_labels, row.names = cell_names,
                      stringsAsFactors = FALSE)
  return(annot)
}





# =============================================================================
# INTERNAL: single object builders
# =============================================================================

#' Build one within-celltype inferCNV object
#'
#' Reference group = one of A/B/C
#' Query          = the other two groups
#'
#' @param counts_mx       genes x cells raw count matrix
#' @param split_metadata  metadata subset for one cell type, with split_group col
#' @param ref_group       "A", "B", or "C"
#' @param gene_order_file path to gene order file (hg38/mm10 etc.)
#' @param chr_exclude     chromosomes to exclude (default c("MT","Y"))
#' @param min_max_counts  c(min, max) counts per cell filter
#'
#' @return inferCNV object
.build_within_object <- function(counts_mx,
                                 split_metadata,
                                 ref_group,
                                 gene_order_file,
                                 chr_exclude    = c("MT", "Y"),
                                 min_max_counts = c(100, 1e6)) {
  
  cells     <- split_metadata$cell_name
  sub_counts <- counts_mx[, cells, drop = FALSE]
  
  annot <- build_annotations_df(
    cell_names   = cells,
    group_labels = split_metadata$split_group
  )
  
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix         = sub_counts,
    annotations_file          = annot,
    gene_order_file           = gene_order_file,
    chr_exclude               = chr_exclude,
    ref_group_names           = ref_group,
    min_max_counts_per_cell   = min_max_counts
  )
  
  return(obj)
}


#' Build one across-celltype inferCNV object
#'
#' Query  = all cells of query_type
#' Reference = all cells of ref_type (full cell type, no splitting)
#'
#' @param counts_mx       genes x cells raw count matrix
#' @param metadata        full metadata data.frame
#' @param cell_type_col   column name for cell type labels
#' @param query_type      cell type string for the query
#' @param ref_type        cell type string for the reference
#' @param gene_order_file path to gene order file
#' @param chr_exclude     chromosomes to exclude
#' @param min_max_counts  c(min, max) counts per cell filter
#'
#' @return inferCNV object
.build_across_object <- function(counts_mx,
                                 metadata,
                                 cell_type_col,
                                 query_type,
                                 ref_type,
                                 gene_order_file,
                                 chr_exclude    = c("MT", "Y"),
                                 min_max_counts = c(100, 1e6)) {
  
  # Subset metadata to only query + reference cell type
  sub_meta <- metadata[metadata[[cell_type_col]] %in% c(as.character(query_type), as.character(ref_type)), ,
                       drop = FALSE]
  
  cells      <- sub_meta$cell_name
  sub_counts <- counts_mx[, cells, drop = FALSE]
  
  # Labels are just the cell type names directly
  annot <- build_annotations_df(
    cell_names   = cells,
    group_labels = sub_meta[[cell_type_col]]
  )
  
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix         = sub_counts,
    annotations_file          = annot,
    gene_order_file           = gene_order_file,
    chr_exclude               = chr_exclude,
    ref_group_names           = ref_type,
    min_max_counts_per_cell   = min_max_counts
  )
  
  return(obj)
}


# =============================================================================
# INTERNAL: mode-level builders
# =============================================================================

#' Build all within-celltype objects across all cell types
#'
#' Returns nested list: list[[cell_type]][[ref_group]] = inferCNV object
#'
.build_all_within <- function(counts_mx,
                              metadata,
                              cell_type_col,
                              gene_order_file,
                              chr_exclude,
                              min_max_counts,
                              n_splits_within) {
  
  cell_types         <- unique(metadata[[cell_type_col]])
  all_split_metadata <- list()
  
  message("\n── Building WITHIN objects ──────────────────────────────────────")
  message(sprintf("   n_splits_within = %d  →  refs: %s",
                  n_splits_within,
                  paste(LETTERS[seq_len(n_splits_within)],
                        collapse = ", ")))
  
  objects <- setNames(
    lapply(cell_types, function(ct) {
      
      number_of_cells <- nrow(
        metadata[metadata[[cell_type_col]] == ct, ])
      
      if (number_of_cells >= 100) {
        
        message("\nCell type: ", ct)
        
        # Split cells for this cell type
        split_meta <- make_splits(
          metadata      = metadata,
          cell_type_col = cell_type_col,
          cell_type_val = ct,
          n_splits      = n_splits_within
        )
        
        # Store split registry
        all_split_metadata[[ct]] <<- split_meta
        
        # One inferCNV object per reference group
        refs <- LETTERS[seq_len(n_splits_within)]
        
        objs <- setNames(
          lapply(refs, function(ref) {
            
            others <- setdiff(refs, ref)
            message(sprintf("  Building ref=%s vs (%s)",
                            ref,
                            paste(others,
                                  collapse = "+")))
            
            tryCatch(
              .build_within_object(
                counts_mx       = counts_mx,
                split_metadata  = split_meta,
                ref_group       = ref,
                gene_order_file = gene_order_file,
                chr_exclude     = chr_exclude,
                min_max_counts  = min_max_counts
              ),
              error = function(e) {
                warning(sprintf(
                  "Failed for %s ref=%s: %s",
                  ct, ref, conditionMessage(e)))
                NULL
              }
            )
          }), refs)
        
        Filter(Negate(is.null), objs)  # ← return value for if block
        
      } else {
        
        message(sprintf(
          "\nCell type: %s skipped — low cells (%d)",
          ct, number_of_cells))
        
        cell_types <<- cell_types[
          !(cell_types == ct)]
      }
      
    }), cell_types)
  
  # Combine split metadata
  split_metadata_combined <- do.call(
    rbind, all_split_metadata)
  rownames(split_metadata_combined) <- NULL
  
  return(list(
    objects        = objects,
    split_metadata = split_metadata_combined
  ))
}


#' Build all across-celltype objects across all cell types
#'
#' Returns nested list: list[[query_type]][[ref_type]] = inferCNV object
#'
.build_all_across <- function(counts_mx,
                              metadata,
                              cell_type_col,
                              gene_order_file,
                              chr_exclude,
                              min_max_counts) {
  
  cell_types <- unique(metadata[[cell_type_col]])
  message("\n── Building ACROSS objects ──────────────────────────────────────")
  
  result <- setNames(lapply(cell_types, function(query) {
    
    # Reference = every other cell type
    ref_types <- setdiff(cell_types, query)
    message("\nQuery: ", query, " | References: ",
            paste(ref_types, collapse = ", "))
    
    objs <- setNames(lapply(ref_types, function(ref) {
      
      message("  Building query=", query, " vs ref=", ref)
      
      tryCatch(
        .build_across_object(
          counts_mx       = counts_mx,
          metadata        = metadata,
          cell_type_col   = cell_type_col,
          query_type      = query,
          ref_type        = ref,
          gene_order_file = gene_order_file,
          chr_exclude     = chr_exclude,
          min_max_counts  = min_max_counts
        ),
        error = function(e) {
          warning("Failed to build across object for query=", query,
                  " ref=", ref, ": ", conditionMessage(e))
          NULL
        }
      )
      
    }), ref_types)
    
    Filter(Negate(is.null), objs)
    
  }), cell_types)
  
  return(result)
}


# =============================================================================
# MAIN: exported function
# =============================================================================

#' Create all inferCNV objects for within and/or across comparisons
#'
#' @param counts_mx       genes x cells raw count matrix
#'                        (e.g. from GetAssayData(seurat_obj, layer = "counts"))
#' @param metadata        data.frame with required columns:
#'                          - 'cell_name': must match colnames(counts_mx)
#'                          - cell_type_col: cell type labels
#' @param cell_type_col   string, name of column in metadata with cell types
#'                        (default "cell_type")
#' @param gene_order_file path to inferCNV gene order file
#'                        (e.g. hg38_gencode_v27.txt)
#' @param mode            one of "within", "across", or "both" (default "both")
#' @param chr_exclude     chromosomes to exclude (default c("MT","Y"))
#' @param min_max_counts  c(min, max) counts per cell (default c(100, 1e6))
#'
#' @return named list with elements 'within_cell_type' and/or 'across_cell_type'
#'         each containing nested lists of inferCNV objects:
#'           within_cell_type[[cell_type]][[ref_group]]   (ref_group: A/B/C)
#'           across_cell_type[[query_type]][[ref_type]]
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(infercnv)
#'
#' counts_mx <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
#'
#' metadata <- data.frame(
#'   cell_name = colnames(seurat_obj),
#'   cell_type = seurat_obj$cell_type
#' )
#'
#' obj_list <- make_infercnv_objects(
#'   counts_mx       = counts_mx,
#'   metadata        = metadata,
#'   cell_type_col   = "cell_type",
#'   gene_order_file = "/path/to/hg38_gencode_v27.txt",
#'   mode            = "both"
#' )
#'
#' saveRDS(obj_list, "infercnv_objcomp.rds")
#' }
make_infercnv_objects <- function(counts_mx,
                                  metadata,
                                  cell_type_col   = "cell_type",
                                  gene_order_file,
                                  mode            = "both",
                                  chr_exclude     = c("MT", "Y"),
                                  min_max_counts  = c(100, 1e6),
                                  n_splits_within) {
  
  # ── Input checks ────────────────────────────────────────────────────────
  if (!mode %in% c("within", "across", "both")) {
    stop("mode must be one of: 'within', 'across', 'both'. Got: '", mode, "'")
  }
  
  if (!file.exists(gene_order_file)) {
    stop("gene_order_file not found: ", gene_order_file)
  }
  
  if (!is.matrix(counts_mx) && !inherits(counts_mx, "dgCMatrix")) {
    stop("counts_mx must be a matrix or dgCMatrix (sparse matrix).")
  }
  
  # ── Validate metadata ────────────────────────────────────────────────────
  message("Validating metadata...")
  metadata   <- validate_metadata(metadata, counts_mx, cell_type_col)
  
  # Align counts to metadata (drop cells not in metadata)
  keep_cells <- intersect(colnames(counts_mx), metadata$cell_name)
  counts_mx  <- counts_mx[, keep_cells, drop = FALSE]
  metadata   <- metadata[metadata$cell_name %in% keep_cells, , drop = FALSE]
  
  # ── Build objects ────────────────────────────────────────────────────────
  result <- list()
  
  if (mode %in% c("within", "both")) {
    result$within_cell_type <- .build_all_within(
      counts_mx       = counts_mx,
      metadata        = metadata,
      cell_type_col   = cell_type_col,
      gene_order_file = gene_order_file,
      chr_exclude     = chr_exclude,
      min_max_counts  = min_max_counts,
      n_splits_within
    )
  }
  
  if (mode %in% c("across", "both")) {
    result$across_cell_type <- .build_all_across(
      counts_mx       = counts_mx,
      metadata        = metadata,
      cell_type_col   = cell_type_col,
      gene_order_file = gene_order_file,
      chr_exclude     = chr_exclude,
      min_max_counts  = min_max_counts
    )
  }
  
  # ── Summary ──────────────────────────────────────────────────────────────
  message("\n── Summary ───────────────────────────────────────────────────────")
  
  if (!is.null(result$within_cell_type)) {
    total_within <- sum(sapply(result$within_cell_type$objects, length))
    message("Within objects built: ", total_within,
            " (", length(result$within_cell_type$objects), " cell types x ",length(result$within_cell_type$objects[[1]]),  " refs)")
  }
  
  if (!is.null(result$across_cell_type)) {
    total_across <- sum(sapply(result$across_cell_type, length))
    message("Across objects built: ", total_across)
  }
  
  message("\nDone. Save with: saveRDS(obj_list, 'infercnv_objcomp.rds')")
  
  return(result)
}




# =============================================================================
# infercnv_run.R
# Function to run inferCNV on all objects produced by make_infercnv_objects()
# =============================================================================


#' Run inferCNV on all objects in a nested list
#'
#' Expects the nested list structure produced by make_infercnv_objects():
#'   list$within_cell_type[[cell_type]][[ref_group]]
#'   list$across_cell_type[[query_type]][[ref_type]]
#'
#' Output directories are created automatically:
#'   {base_outdir}/within/{cell_type}/ref_{ref_group}/
#'   {base_outdir}/across/{query_type}/ref_{ref_type}/
#'
#' @param infercnv_obj_list nested list from make_infercnv_objects()
#' @param base_outdir       root output directory (will be created if needed)
#' @param cutoff            minimum average counts per gene for reference cells
#'                          (default 0.1 — correct for RNA-seq, use 1 for 10x)
#' @param cluster_by_groups logical, cluster cells within groups (default TRUE)
#' @param HMM               logical, run HMM CNV prediction (default FALSE)
#' @param denoise           logical, apply denoising (default TRUE)
#' @param analysis_mode     one of "subclusters", "samples", "cells"
#'                          (default "subclusters")
#' @param window_length     smoothing window length (default 140)
#' @param plot_steps        logical, save intermediate plots (default FALSE)
#' @param no_plot           logical, skip final heatmap (default TRUE —
#'                          set FALSE if you want plots)
#' @param no_prelim_plot    logical, skip preliminary plot (default TRUE)
#' @param plot_probabilities logical (default FALSE)
#' @param diagnostics       logical (default FALSE)
#' @param inspect_subclusters logical (default FALSE)
#' @param resume_if_exists  logical, skip runs where out_dir already has
#'                          run.final.infercnv_obj (default TRUE)
#'
#' @return invisibly returns a data.frame log of all runs with status
#'         (success / failed / skipped)
#'
#' @examples
#' \dontrun{
#' obj_list <- readRDS("infercnv_objcomp.rds")
#'
#' run_log <- run_infercnv_objects(
#'   infercnv_obj_list = obj_list,
#'   base_outdir       = "/path/to/output/",
#'   no_plot           = FALSE   # set TRUE on HPC to skip heavy plotting
#' )
#'
#' # Check what failed
#' run_log[run_log$status == "failed", ]
#' }
run_infercnv_objects <- function(infercnv_obj_list,
                                 base_outdir,
                                 cutoff               = 0.1,
                                 cluster_by_groups    = TRUE,
                                 HMM                  = FALSE,
                                 denoise              = TRUE,
                                 analysis_mode        = "subclusters",
                                 window_length        = 140,
                                 plot_steps           = FALSE,
                                 no_plot              = TRUE,
                                 no_prelim_plot       = TRUE,
                                 plot_probabilities   = FALSE,
                                 diagnostics          = FALSE,
                                 inspect_subclusters  = FALSE,
                                 resume_if_exists     = TRUE) {
  
  # ── Input checks ────────────────────────────────────────────────────────
  valid_modes <- c("within_cell_type", "across_cell_type")
  found_modes <- intersect(names(infercnv_obj_list), valid_modes)
  
  if (length(found_modes) == 0) {
    stop(
      "infercnv_obj_list does not contain expected modes. ",
      "Expected names: 'within_cell_type' and/or 'across_cell_type'. ",
      "Got: ", paste(names(infercnv_obj_list), collapse = ", ")
    )
  }
  
  dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Run log ──────────────────────────────────────────────────────────────
  run_log <- data.frame(
    mode      = character(),
    cell_type = character(),
    comp      = character(),
    out_dir   = character(),
    status    = character(),
    message   = character(),
    stringsAsFactors = FALSE
  )
  
  # ── Map mode names to output folder names ────────────────────────────────
  mode_folder_map <- c(
    within_cell_type = "within",
    across_cell_type = "across"
  )
  
  # ── Loop ─────────────────────────────────────────────────────────────────
  for (mode in found_modes) {
    
    mode_folder  <- mode_folder_map[[mode]]
    mode_objects <- infercnv_obj_list[[mode]]
    
    for (cell_type in names(mode_objects)) {
      
      type_objects <- mode_objects[[cell_type]]
      
      for (comp in names(type_objects)) {
        
        infer_obj <- type_objects[[comp]]
        
        # NULL guard — object may have failed during creation
        if (is.null(infer_obj)) {
          message("Skipping NULL object: ", mode, " / ", cell_type, " / ", comp)
          run_log <- rbind(run_log, data.frame(
            mode = mode, cell_type = cell_type, comp = comp,
            out_dir = NA, status = "skipped_null", message = "Object is NULL",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        # Build output directory
        outdir <- file.path(base_outdir, mode_folder, cell_type, comp)
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        
        # Resume check — skip if already completed
        final_obj_path <- file.path(outdir, "run.final.infercnv_obj")
        if (resume_if_exists && file.exists(final_obj_path)) {
          message("Skipping (already complete): ",
                  mode_folder, "/", cell_type, "/", comp)
          run_log <- rbind(run_log, data.frame(
            mode = mode, cell_type = cell_type, comp = comp,
            out_dir = outdir, status = "skipped_exists",
            message = "run.final.infercnv_obj already present",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        message("\n── Running: ", mode_folder, " / ", cell_type, " / ", comp)
        message("   Output: ", outdir)
        
        status  <- "success"
        err_msg <- ""
        
        tryCatch({
          options(scipen = 100)
          infercnv::run(
            infercnv_obj          = infer_obj,
            out_dir               = outdir,
            cutoff                = cutoff,
            cluster_by_groups     = cluster_by_groups,
            HMM                   = HMM,
            denoise               = denoise,
            analysis_mode         = analysis_mode,
            output_format         = NA,
            no_plot               = no_plot,
            no_prelim_plot        = no_prelim_plot,
            window_length         = window_length,
            plot_probabilities    = plot_probabilities,
            plot_steps            = plot_steps,
            diagnostics           = diagnostics,
            inspect_subclusters   = inspect_subclusters
          )
          
          message("   Done: ", mode_folder, " / ", cell_type, " / ", comp)
          
        }, error = function(e) {
          status  <<- "failed"
          err_msg <<- conditionMessage(e)
          warning("FAILED: ", mode_folder, "/", cell_type, "/", comp,
                  "\n  Error: ", err_msg)
        })
        
        run_log <- rbind(run_log, data.frame(
          mode      = mode,
          cell_type = cell_type,
          comp      = comp,
          out_dir   = outdir,
          status    = status,
          message   = err_msg,
          stringsAsFactors = FALSE
        ))
        
        # Clean up memory between runs
        rm(infer_obj)
        gc()
        
      } # comp loop
    } # cell_type loop
  } # mode loop
  
  # ── Final summary ────────────────────────────────────────────────────────
  message("\n── Run Summary ───────────────────────────────────────────────────")
  print(table(run_log$status))
  
  failed <- run_log[run_log$status == "failed", ]
  if (nrow(failed) > 0) {
    message("\nFailed runs:")
    print(failed[, c("mode", "cell_type", "comp", "message")])
  }
  
  # Save log to base_outdir
  log_path <- file.path(base_outdir, "run_log.tsv")
  write.table(run_log, log_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("\nRun log saved to: ", log_path)
  
  invisible(run_log)
}



run_infercnv_pipeline <- function(
    
  # ---- make_infercnv_objects parameters ----------------------------------
  counts_mx,
  metadata,
  cell_type_col   = "cell_type",
  gene_order_file,
  mode            = c("within", "across"),
  chr_exclude     = c("MT", "Y"),
  min_max_counts  = c(100, 1e6),
  n_splits_within = 3,
  
  # ---- run_infercnv_objects parameters -----------------------------------
  base_outdir,
  cutoff          = 0.1,
  cluster_by_groups = TRUE,
  HMM             = FALSE,
  denoise         = TRUE,
  analysis_mode   = "subclusters",
  window_length   = 140,
  no_plot         = TRUE,
  resume_if_exists = TRUE
  
) {
  
  mode <- match.arg(mode)
  
  # ---- Input validation ---------------------------------------------------
  if (!is.matrix(counts_mx) && !inherits(counts_mx, "dgCMatrix")) {
    stop("counts_mx must be a matrix or sparse matrix (dgCMatrix).")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame.")
  }
  if (!cell_type_col %in% colnames(metadata)) {
    stop("cell_type_col '", cell_type_col, "' not found in metadata.")
  }
  if (!file.exists(gene_order_file)) {
    stop("gene_order_file not found: ", gene_order_file)
  }
  if (!dir.exists(base_outdir)) {
    message("base_outdir does not exist — creating: ", base_outdir)
    dir.create(base_outdir, recursive = TRUE)
  }
  if (length(min_max_counts) != 2L) {
    stop("min_max_counts must be a numeric vector of length 2: c(min, max).")
  }
  
  message(sprintf(paste0(
    "inferCNV pipeline starting:\n",
    "  Mode:           %s\n",
    "  Cells:          %d\n",
    "  Genes:          %d\n",
    "  Cell types:     %s\n",
    "  Output dir:     %s"
  ),
  mode,
  ncol(counts_mx),
  nrow(counts_mx),
  paste(unique(metadata[[cell_type_col]]), collapse = ", "),
  base_outdir
  ))
  
  t_start <- proc.time()
  
  # ---- Step 1: make inferCNV objects --------------------------------------
  message("\n[1/2] Creating inferCNV objects...")
  t_make_start <- proc.time()
  
  obj_list <- make_infercnv_objects(
    counts_mx       = counts_mx,
    metadata        = metadata,
    cell_type_col   = cell_type_col,
    gene_order_file = gene_order_file,
    mode            = mode,
    chr_exclude     = chr_exclude,
    min_max_counts  = min_max_counts,
    n_splits_within = n_splits_within
  )
  
  t_make_end <- proc.time()
  
  # ---- Lightweight sanity check between steps -----------------------------
  if (is.null(obj_list)) {
    stop("make_infercnv_objects() returned NULL — check inputs.")
  }
  if (!"objects" %in% names(obj_list[["within_cell_type"]]) &&
      mode %in% c("within")) {
    stop(
      "Expected 'objects' element in obj_list$within_cell_type. ",
      "make_infercnv_objects() may have failed silently."
    )
  }
  
  message(sprintf(
    "Objects created in %.1f seconds.",
    (t_make_end - t_make_start)[["elapsed"]]
  ))
  
  split_metadata <- NULL  # default — only populated for within/both
  
  if (mode %in% c("within")) {
    
    if (!"split_metadata" %in% names(obj_list[["within_cell_type"]])) {
      stop(
        "Expected 'split_metadata' in obj_list$within_cell_type. ",
        "Check make_infercnv_objects() output structure."
      )
    }
    
    split_metadata <- obj_list[["within_cell_type"]][["split_metadata"]]
    obj_list[["within_cell_type"]] <- obj_list[["within_cell_type"]][["objects"]]
    
    message(sprintf(
      "Within mode: split_metadata extracted (%d cell type splits)",
      length(split_metadata)
    ))
  }
  
  
  # ---- Step 2: run inferCNV -----------------------------------------------
  t_run_start <- proc.time()
  
  run_log <- run_infercnv_objects(
    infercnv_obj_list = obj_list,
    base_outdir       = base_outdir,
    cutoff            = cutoff,
    cluster_by_groups = cluster_by_groups,
    HMM               = HMM,
    denoise           = denoise,
    analysis_mode     = analysis_mode,
    window_length     = window_length,
    no_plot           = no_plot,
    resume_if_exists  = resume_if_exists
  )
  
  t_run_end <- proc.time()
  
  # ---- Lightweight sanity check on run_log --------------------------------
  if (is.null(run_log)) {
    warning("run_infercnv_objects() returned NULL run_log — check output directory.")
  }
  
  t_end     <- proc.time()
  runtime   <- (t_end - t_start)[["elapsed"]]
  
  message(sprintf(paste0(
    "\nPipeline complete:\n",
    "  Make time:   %.1f seconds\n",
    "  Run time:    %.1f seconds\n",
    "  Total time:  %.1f seconds (%.1f minutes)"
  ),
  (t_make_end - t_make_start)[["elapsed"]],
  (t_run_end  - t_run_start)[["elapsed"]],
  runtime,
  runtime / 60
  ))
  
  
  list(
    obj_list       = obj_list,
    run_log        = run_log,
    metadata       = split_metadata,
    runtime        = runtime
  )
}



#' Convert a wide expression matrix to long format
#'
#' Reshapes a gene-by-cell expression table into a long-format data frame with
#' one row per gene-cell combination.
#'
#' @param expr_df A data frame or matrix with genes as rows and cells as columns.
melt_expr_to_long <- function(expr_df, cell_prefix = "cell_") {
  
  if (is.null(rownames(expr_df)) || 
      all(rownames(expr_df) == as.character(seq_len(nrow(expr_df))))) {
    stop("expr_df has no meaningful rownames — gene names expected as rownames.")
  }
  
  expr_df |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(
      cols      = -gene,
      names_to  = "cell_name",
      values_to = "state_raw"
    ) |>
    # Strip prefix from cell_name after melting if present — silently skips
    # cells that don't carry the prefix, so mixed naming is handled safely
    dplyr::mutate(
      cell_name = stringr::str_remove(cell_name, paste0("^", cell_prefix))
    )
}




#' Attach genomic coordinates to long-format gene data
#'
#' Merges long-format gene-level data with a gene order table containing genomic
#' coordinates.
#'
#' @param long_df A long-format data frame containing a \code{gene} column.
#' @param gene_order A data frame containing gene coordinates. Row names are
#'   assumed to store gene names.
#'
#' @return A merged data frame containing gene-level values and genomic
#'   coordinate columns.
attach_gene_order <- function(long_df, gene_order) {
  
  gene_order["gene"] <- rownames(gene_order)
  
  required <- c("gene", "chr", "start", "stop")
  if (!all(required %in% colnames(gene_order))) {
    stop("gene_order missing required columns")
  }
  
  merged <- long_df |>
    inner_join(by = "gene",
               y = gene_order)
  
  if (nrow(merged) == 0) {
    stop("Merge produced empty table — gene names may not match between expr and gene_order.")
  }
  if (nrow(merged) < nrow(long_df)) {
    warning(sprintf(
      "Inner join dropped %d rows — some genes have no coordinate annotation.",
      nrow(long_df) - nrow(merged)
    ))
  }
  
  merged
}




#' Discretize raw CNV scores into gain, loss, or neutral states
#' 
#' Converts a continuous inferCNV score into a categorical CNV state using
#' mean +/- k * SD thresholds computed globally across the input data frame.
#' 
#' 
#' Global thresholding is appropriate here because inferCNV output values are
#' already centred by the tool's internal smoothing relative to the reference
#' signal. The neutral state dominates the distribution and anchors the global
#' mean near the neutral baseline, so gains and losses appear as symmetric
#' tail deviations. Per-cell thresholding would artificially force every cell
#' to have calls even in largely euploid cells.
#'
#' @param df A data frame containing a numeric \code{state_raw} column.
#' @param k Numeric multiplier for the standard deviation cutoff. Default is 1.5.
discretize_cnv_state_infer_cnv <- function(df, k =1.5) {
  if (!"state_raw" %in% colnames(df)) {
    stop("Missing required column: state_raw")
  }
  if (!is.numeric(df$state_raw)) {
    stop("state_raw must be numeric")
  }
  if (all(is.na(df$state_raw))) {
    stop("state_raw is entirely NA — check upstream expression extraction")
  }
  
  mu <- mean(df$state_raw, na.rm = TRUE)
  sigma <- sd(df$state_raw, na.rm = TRUE)
  
  upper <- mu + k * sigma
  lower <- mu - k * sigma
  df |>
    mutate(
      state = case_when(
        state_raw < lower ~ "loss",
        state_raw > upper ~ "gain",
        TRUE ~ "neutral"
      )
    )
}




#' Load and prepare inferCNV reference outputs
#'
#' Converts inferCNV outputs into a unified long-format table with genomic
#' coordinates, discrete CNV states, and reference labels.
#'
#' @param infercnv_list A named list of inferCNV objects. Each object is loaded and the 1 and 2 are taken
#'
#' @return A tibble containing gene-level CNV calls across all references.
load_and_prepare_infercnv_reference <- function(infercnv_list) {
  
  refs <- names(infercnv_list)
  
  res <- lapply(seq_along(infercnv_list), function(i) {
    
    infercnv_obj <- infercnv_list[[i]]
    reference <- refs[i]
    melt_expr_to_long(as.data.frame(infercnv_obj[[1]])) |>
      attach_gene_order(infercnv_obj[[2]]) |>
      discretize_cnv_state_infer_cnv() |>
      dplyr::mutate(reference = reference)
  })
  
  dplyr::bind_rows(res)
}



#' Discover and load inferCNV run outputs
#'
#' Finds inferCNV result files matching a pattern inside one or more reference
#' directories and loads them as R objects with \code{readRDS()}.
#'
#' @param base_dir Optional base directory containing the reference-specific
#'   subdirectories.
#' @param ref_dirs Character vector of reference directory names or paths.
#' @param pattern Regular expression used to match inferCNV result files.
#' 
#' @return A named list of loaded inferCNV objects, one per reference directory.
#'
#' @details
#' For each entry in \code{ref_dirs}, the function searches for files matching
#' \code{pattern}. Matching files are loaded with \code{readRDS()} and stored in
#' a list using the reference name as the list element name.
discover_infercnv_runs <- function(base_dir = NULL,
                                   ref_dirs,
                                   pattern = "^run\\.final") {
  ref_paths <- if (!is.null(base_dir)) file.path(base_dir, ref_dirs) else ref_dirs
  names(ref_paths) <- ref_dirs
  
  runs <- lapply(ref_paths, function(ref_path){
    files <- list.files(ref_path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) stop(sprintf("No files matching '%s' found in %s", pattern, ref_path))
    if (length(files)  > 1) stop(sprintf("Multiple files matching '%s' found in %s", pattern, ref_path))
    inferobj <- readRDS(files[[1]])
  })
  
}


load_infercnv_data <- function(base_dir, ref_dirs, pattern = "^run\\.final") {
  
  # Discover and load runs — already written as discover_infercnv_runs()
  infer_objs <- discover_infercnv_runs(
    base_dir = base_dir,
    ref_dirs = ref_dirs,
    pattern  = pattern
  )
  
  # Convert to standard schema — already written
  infer_objs_1 <- lapply(infer_objs, function(x) {
    list(x@expr.data, x@gene_order)
  })
  
  load_and_prepare_infercnv_reference(infer_objs_1)
}





