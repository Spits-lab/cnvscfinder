<<<<<<< HEAD
library(ggplot2)
library(data.table)
library(dplyr)

#' Plot the density distribution of CNV lengths
#'
#' Creates a density plot of CNV lengths in megabases, optionally restricted
#' to a single CNV state such as gain or loss.
#'
#' @param dt A data frame containing at least the columns cnv_length_mb
#'   and optionally cnv_state.
#' @param state Optional character string specifying a CNV state to subset,
#'   such as "gain" or "loss". If NULL, all rows are used.
#' @param thresholds Numeric vector of CNV length thresholds to display as
#'   vertical dashed lines.
#' @param fill_color Fill color for the density polygon. Default is
#' @param title Plot title. Default is "CNV length distribution".
#'
#' @return A ggplot2 object.
plot_cnv_density <- function(
    dt,
    state = NULL,
    thresholds = c(5, 25, 50),
    fill_color = "grey40",
    title = "CNV length distribution"
) {
  
  # Subset efficiently (no copy)
  if (!is.null(state)) {
    dt_sub <- dt[dt$cnv_state == state,]
  } else {
    dt_sub <- dt
  }
  
  ggplot(dt_sub, aes(x = cnv_length_mb)) +
    geom_density(fill = fill_color, alpha = 0.4) +
    geom_vline(
      xintercept = thresholds,
      linetype = "dashed",
      alpha = 0.6
    ) +
    labs(
      title = title,
      x = "CNV length (Mb)",
      y = "Density"
    ) +
    theme_minimal()
}



#' Plot overall and state-specific CNV length distributions
#'
#' Builds three density plots showing CNV length distributions for all CNVs,
#' gain events, and loss events, then stacks them vertically.
#'
#' @param dt A data frame containing CNV length information, including
#'   cnv_length_mb and cnv_state.
#' @param thresholds Numeric vector of thresholds to display as dashed vertical
#'   lines in each panel. Default is c(5, 25, 50).
#'
#' @return A combined patchwork plot object.
plot_all_cnv_distributions <- function(
    dt,
    thresholds = c(5, 25, 50)
) {
  
  p_overall <- plot_cnv_density(
    dt = dt,
    state = NULL,
    thresholds = thresholds,
    fill_color = "grey40",
    title = "Overall CNV length distribution"
  )
  
  p_gain <- plot_cnv_density(
    dt = dt,
    state = "gain",
    thresholds = thresholds,
    fill_color = "steelblue",
    title = "Gain CNV length distribution"
  )
  
  p_loss <- plot_cnv_density(
    dt = dt,
    state = "loss",
    thresholds = thresholds,
    fill_color = "firebrick",
    title = "Loss CNV length distribution"
  )
  
  p_overall / p_gain / p_loss
}




#' Plot density distributions by level and type
#'
#' Creates a density plot for percentage values within a selected level,
#' grouped by type.
#'
#' @param level_name Character string specifying which level in
#'   plot_long$level to plot.
#' @param plot_long A long-format data frame containing at least the columns
#'   level, type, and percentage.
#' @param threshold Numeric threshold shown as a dashed vertical line.
#'
#' @return A ggplot2 object.
#'
#' @details
#' The function also computes the mean percentage per type, although the
#' current vertical reference line uses the supplied threshold value.
make_density_plot <- function(level_name,plot_long, threshold) {
  
  df_sub <- filter(plot_long, level == level_name)
  
  # compute means per type
  mean_df <- df_sub %>%
    group_by(type) %>%
    summarise(mean_value = mean(percentage, na.rm = TRUE),
              .groups = "drop")
  
  ggplot(df_sub,
         aes(x = percentage, fill = type, color = type)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    
    # vertical mean lines (like abline)
    geom_vline(data = mean_df,
               aes(xintercept = threshold),
               linetype = "dashed",
               linewidth = 1.2,
               show.legend = FALSE) +
    
    labs(
      title = level_name,
      x = "Percentage",
      y = "Density"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")
}





################################################################
## Functions to Process Information for conjoined heatmap#######
#################################################################
#' Prepare genome-wide chromosome and arm coordinates
#'
#' Builds chromosome-level cumulative coordinates and lookup tables for
#' chromosome arms, enabling mapping of chromosome-local intervals into a
#' genome-wide coordinate system.
#'
#' @param chromosome_arms A data frame containing chromosome arm annotation.
#'   It should include at least \code{chr}, \code{arm}, \code{arm_start},
#'   \code{arm_end}, and \code{arm_length}.
#'
#' @return A named list containing information the different fraction from the chromossome
#'
prepare_genome_structure <- function(chromosome_arms) {
  
  chromosome_lengths <- chromosome_arms %>%
    group_by(chr) %>%
    summarise(chr_length = sum(arm_length), .groups = "drop") %>%
    mutate(
      chr_length = as.numeric(chr_length),
      chr_num = suppressWarnings(as.numeric(gsub("chr", "", chr)))
    ) %>%
    arrange(chr_num, chr) %>%
    dplyr::select(-chr_num) %>%
    mutate(
      
      chr_start = lag(cumsum(chr_length), default = 0),
      chr_end   = chr_start + chr_length
    )
  
  # Precompute arm lookup tables (vectorized, no match needed later)
  arm_lookup <- chromosome_arms %>%
    left_join(chromosome_lengths, by = "chr") %>%
    mutate(
      genome_arm_start = chr_start + arm_start,
      genome_arm_end   = chr_start + arm_end
    )
  
  p_lookup <- arm_lookup %>%
    filter(arm == "p") %>%
    select(chr,
           p_genome_start = genome_arm_start,
           p_genome_end   = genome_arm_end)
  
  q_lookup <- arm_lookup %>%
    filter(arm == "q") %>%
    select(chr,
           q_genome_start = genome_arm_start,
           q_genome_end   = genome_arm_end)
  
  list(
    chromosome_lengths = chromosome_lengths,
    p_lookup = p_lookup,
    q_lookup = q_lookup,
    arm_lookup = arm_lookup
  )
}

#' Map CNV intervals to genome-wide coordinates
#'
#' Converts chromosome-level CNV coordinates into cumulative genome-wide
#' coordinates for visualization. The function joins chromosome and chromosome-arm
#' lookup tables, computes genome-wide start and end positions, optionally
#' expands large events to full chromosome or arm boundaries for plotting,
#' and encodes CNV state numerically.
#'
#' @param cnv_filtered A data frame of CNV events
#' @param genome_structure
#' @param threshold Optional numeric threshold used to expand plotted CNV
#'   intervals to whole-chromosome or chromosome-arm boundaries.
#' @return A data frame with additional columns
map_cnv_to_genome <- function(cnv_filtered,
                              genome_structure,
                              threshold = NULL,
                              arrange_df_cols) {
  
  chromosome_lengths <- genome_structure$chromosome_lengths
  p_lookup <- genome_structure$p_lookup
  q_lookup <- genome_structure$q_lookup
  
  
  cnv_mapped <- cnv_filtered %>%
    left_join(chromosome_lengths, by = "chr") %>%
    left_join(p_lookup, by = "chr") %>%
    left_join(q_lookup, by = "chr") %>%
    mutate(
      genome_start = chr_start + start,
      genome_end   = chr_start + end
    )
  
  # If threshold is NULL → no override logic
  if (is.null(threshold)) {
    
    cnv_mapped <- cnv_mapped %>%
      mutate(
        genome_start_plot = genome_start,
        genome_end_plot   = genome_end
      )
    
  } else {
    
    cnv_mapped <- cnv_mapped %>%
      mutate(
        genome_start_plot = case_when(
          whole_chromosome_gain > threshold |
            whole_chromosome_loss > threshold ~ chr_start,
          
          p_arm_gain > threshold |
            p_arm_loss > threshold ~ p_genome_start,
          
          q_arm_gain > threshold |
            q_arm_loss > threshold ~ q_genome_start,
          
          TRUE ~ genome_start
        ),
        genome_end_plot = case_when(
          whole_chromosome_gain > threshold |
            whole_chromosome_loss > threshold ~ chr_end,
          
          p_arm_gain > threshold |
            p_arm_loss > threshold ~ p_genome_end,
          
          q_arm_gain > threshold |
            q_arm_loss > threshold ~ q_genome_end,
          
          TRUE ~ genome_end
        )
      )
  }
 
  distinct_dfs_cols <- c("cell_name", arrange_df_cols)
  cell_order <- cnv_filtered %>%
    distinct(across(all_of(distinct_dfs_cols))) %>%
    arrange(all_of(arrange_df_cols))
  
  cnv_mapped <- cnv_mapped %>%
    mutate(
      cell_name = factor(cell_name, levels = cell_order$cell_name),
      cell_id = as.numeric(cell_name),
      cnv_state_numeric = case_when(
        cnv_state == "gain" ~ 1,
        cnv_state == "loss" ~ -1,
        TRUE ~ 0
      )
    ) %>%
    arrange(arrange_df_cols)
}

#' Plot a chromosome ideogram
#'
#' Creates a simple ideogram-style plot of chromosome arms aligned to genome-wide
#' coordinates.
#'
#' @param genome_structure Output from prepare_genome_structure().
#' @param max_cell_id Maximum cell index used to place the ideogram above the
#'   heatmap.
#' @param arm_colors Named vector of colors for chromosome arms.
#'
#' @return A ggplot2 object.
plot_ideogram <- function(genome_structure,
                          max_cell_id,
                          arm_colors = c("p"="#4DBBD5",
                                         "cen"="black",
                                         "q"="#E64B35")) {
  
  chromosome_lengths <- genome_structure$chromosome_lengths
  arm_lookup <- genome_structure$arm_lookup
  
  arm_plot <- arm_lookup %>%
    mutate(
      ymin = max_cell_id + 1,
      ymax = ymin + 1
    )
  
  ggplot(arm_plot) +
    geom_rect(aes(
      xmin = genome_arm_start,
      xmax = genome_arm_end,
      ymin = ymin,
      ymax = ymax,
      fill = arm
    ), color = NA) +
    scale_fill_manual(values = arm_colors, name = "Arm") +
    theme_void() +
    scale_x_continuous(
      breaks = chromosome_lengths$chr_start,
      labels = chromosome_lengths$chr,
      expand = c(0,0),
      limits = c(0, max(chromosome_lengths$chr_end) +1e6)
    )
}



#' Prepare CNV karyotype plot data
#'
#' Assigns plot_idx per cell ordered by grouping columns then CNV burden,
#' validates genome-wide coordinates, and prepares boundary line positions.
#'
#' @param cnv_mapped Data frame with one row per cell per CNV segment.
#'   Must contain: cell_name, genome_start_plot, genome_end_plot,
#'   cnv_state, cnv_length_mb, and all columns in grouping_cols.
#' @param genome_structure Output of prepare_genome_structure().
#' @param grouping_cols Character vector of columns defining cell groups
#'   for ordering and boundary lines. Default c("cell_type").
#'   Can include "embryo_id", "cell_type", or both.
#' @param state_colors Named character vector mapping cnv_state to colours.
#'
#' @return Named list with plot_data, chr_lengths, genome_size,
#'   boundary_lines, cell_order, grouping_cols.
prepare_cnv_plot <- function(
    cnv_mapped,
    genome_structure,
    grouping_cols = c("cell_type"),
    state_colors  = c(
      "gain" = "#E64B35",
      "loss" = "#4DBBD5"
    )
) {
  
  # ---- Input validation ---------------------------------------------------
  required_cols <- c(
    "cell_name", "chr",
    "genome_start_plot", "genome_end_plot",
    "cnv_state", "cnv_length_mb"
  )
  missing_cols <- setdiff(required_cols, colnames(cnv_mapped))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  missing_group <- setdiff(grouping_cols, colnames(cnv_mapped))
  if (length(missing_group) > 0L) {
    stop("grouping_cols not found in cnv_mapped: ",
         paste(missing_group, collapse = ", "))
  }
  
  unknown_states <- setdiff(
    unique(cnv_mapped$cnv_state), names(state_colors)
  )
  if (length(unknown_states) > 0L) {
    warning(
      "cnv_state values not in state_colors — will appear grey: ",
      paste(unknown_states, collapse = ", ")
    )
  }
  
  # ---- Extract genome structure -------------------------------------------
  chr_lengths <- genome_structure$chromosome_lengths
  genome_size <- sum(chr_lengths$chr_length)
  
  # ---- Compute cell ordering ----------------------------------------------
  # Order: grouping_cols → total CNV burden descending → cell_name
  # One row per cell for ordering purposes
  cell_order <- cnv_mapped |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(grouping_cols, "cell_name")
    ))) |>
    dplyr::summarise(
      total_cnv_mb = sum(cnv_length_mb, na.rm = TRUE),
      .groups      = "drop"
    ) |>
    dplyr::arrange(
      dplyr::across(dplyr::all_of(grouping_cols)),
      dplyr::desc(total_cnv_mb),
      cell_name
    ) |>
    dplyr::mutate(plot_idx = dplyr::row_number())
  
  message(sprintf(
    "Cell ordering: %d unique cells across %d group(s)",
    nrow(cell_order),
    dplyr::n_distinct(cell_order[, grouping_cols])
  ))
  
  # ---- Assign plot_idx to cnv_mapped --------------------------------------
  # Join back — every segment row for a cell gets same plot_idx
  plot_data <- cnv_mapped |>
    dplyr::select(-dplyr::any_of("plot_idx")) |>  # remove if exists
    dplyr::left_join(
      cell_order |> dplyr::select(cell_name, plot_idx),
      by = "cell_name"
    )
  
  # Sanity check — every cell should have exactly one plot_idx
  idx_check <- plot_data |>
    dplyr::group_by(cell_name) |>
    dplyr::summarise(
      n_idx = dplyr::n_distinct(plot_idx),
      .groups = "drop"
    ) |>
    dplyr::filter(n_idx > 1L)
  
  if (nrow(idx_check) > 0L) {
    stop(sprintf(
      "%d cell(s) have multiple plot_idx values — check cell_name uniqueness: %s",
      nrow(idx_check),
      paste(idx_check$cell_name[1:min(5, nrow(idx_check))], collapse = ", ")
    ))
  }
  
  # ---- Compute boundary lines ---------------------------------------------
  # Boundaries sit between groups defined by grouping_cols
  # Position = last plot_idx of each group + 0.5
  boundary_lines <- cell_order |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) |>
    dplyr::summarise(
      last_idx = max(plot_idx),
      .groups  = "drop"
    ) |>
    # Don't add boundary after last group
    dplyr::filter(last_idx < max(cell_order$plot_idx)) |>
    dplyr::pull(last_idx) + 0.5
  
  message(sprintf(
    "%d boundary line(s) computed between groups",
    length(boundary_lines)
  ))
  
  # ---- Coordinate validation ----------------------------------------------
  n_na_coords <- sum(
    is.na(plot_data$genome_start_plot) |
      is.na(plot_data$genome_end_plot)
  )
  if (n_na_coords > 0L) {
    warning(sprintf(
      "%d segments have NA genome coordinates — will not be plotted.",
      n_na_coords
    ))
  }
  
  out_of_range <- plot_data |>
    dplyr::filter(
      genome_start_plot < 0 |
        genome_end_plot > genome_size |
        genome_start_plot > genome_end_plot
    )
  
  if (nrow(out_of_range) > 0L) {
    warning(sprintf(
      "%d segments have out-of-range coordinates — check genome_structure.",
      nrow(out_of_range)
    ))
  }
  
  message(sprintf(
    "Plot data ready: %d segments across %d cells",
    nrow(plot_data),
    dplyr::n_distinct(plot_data$cell_name)
  ))
  
  list(
    plot_data      = plot_data,
    chr_lengths    = chr_lengths,
    genome_size    = genome_size,
    boundary_lines = boundary_lines,
    cell_order     = cell_order,
    grouping_cols  = grouping_cols,
    state_colors   = state_colors
  )
}


#' Plot CNV karyotype with ideogram
#'
#' Builds the main CNV heatmap and combines it with a chromosome
#' ideogram using patchwork. Ideogram height is fixed at ideogram_ratio
#' of the total plot height.
#'
#' @param panel_data Output of prepare_cnv_plot().
#' @param genome_structure Output of prepare_genome_structure().
#' @param title Plot title.
#' @param boundary_lines Numeric vector of y positions for group
#'   boundaries. If NULL uses boundary_lines from panel_data.
#' @param ideogram_ratio Fraction of total height for ideogram. Default 0.08.
#' @param arm_colors Named vector of colours for chromosome arms.
#' @param show_legend Logical. Default TRUE.
#'
#' @return A patchwork object combining ideogram + main plot.
plot_cnv_karyotype <- function(
    panel_data,
    genome_structure,
    boundary_lines  = NULL,
    ideogram_ratio  = 0.08,
    arm_colors      = c(
      "p"   = "#4DBBD5",
      "cen" = "black",
      "q"   = "#E64B35"
    ),
    show_legend     = TRUE,
    cell_sizes      = NULL  
) {
  
  plot_data      <- panel_data$plot_data
  chr_lengths    <- panel_data$chr_lengths
  genome_size    <- panel_data$genome_size
  state_colors   <- panel_data$state_colors
  cell_order     <- panel_data$cell_order
  
  # Use provided boundary_lines or fall back to prepared ones
  boundaries <- boundary_lines %||% panel_data$boundary_lines
  
  # ---- Chromosome labels --------------------------------------------------
  chr_labels <- chr_lengths |>
    dplyr::mutate(
      label_pos = chr_start + chr_length / 2,
      label     = gsub("chr", "", chr),
      label     = dplyr::case_when(
        label == "23" ~ "X",
        label == "24" ~ "Y",
        TRUE          ~ label
      )
    )
  
  chr_boundaries <- chr_lengths |>
    dplyr::filter(chr_start > 0)
  
  # ---- Y axis labels ------------------------------------------------------
  # One label per cell — use cell_order for correct positioning
  y_label_data <- panel_data$cell_order |>
    dplyr::group_by(dplyr::across(dplyr::all_of(panel_data$grouping_cols))) |>
    dplyr::summarise(
      mid_idx = mean(plot_idx),   # center of group for label position
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Combine grouping cols into one label if multiple
      group_label = apply(
        dplyr::pick(dplyr::all_of(panel_data$grouping_cols)),
        1,
        paste, collapse = " — "
      )
    )
  
  
  if (!is.null(cell_sizes)) {
    
    # Count CNV positive per group
    cnv_positive <- plot_data %>%
      dplyr::filter(cnv_state != "none") %>%
      dplyr::distinct(
        cell_name,
        dplyr::across(
          dplyr::all_of(panel_data$grouping_cols)
        )
      ) %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(panel_data$grouping_cols)
        )
      ) %>%
      dplyr::summarise(
        n_cnv = dplyr::n(),
        .groups = "drop"
      )
    
    # Total per group from cell_sizes
    total_per_group <- cell_sizes %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(panel_data$grouping_cols)
        )
      ) %>%
      dplyr::summarise(
        n_total = dplyr::n(),
        .groups = "drop"
      )
    
    y_label_data <- y_label_data %>%
      dplyr::left_join(cnv_positive,
                       by = panel_data$grouping_cols) %>%
      dplyr::left_join(total_per_group,
                       by = panel_data$grouping_cols) %>%
      dplyr::mutate(
        pct_cnv   = round(100 * n_cnv / n_total, 1),
        pct_empty = round(100 - pct_cnv, 1),
        group_label = paste0(
          group_label, "\n",
          "CNV: ", pct_cnv, "% | ",
          "Empty: ", pct_empty, "%"
        )
      )
  }
  
  
  # ---- Main plot ----------------------------------------------------------
  p_main <- ggplot2::ggplot(plot_data) +
    
    # White background row per cell
    ggplot2::geom_rect(
      data = plot_data |>
        dplyr::distinct(cell_name, plot_idx),
      ggplot2::aes(
        xmin = 0,
        xmax = genome_size,
        ymin = plot_idx - 0.5,
        ymax = plot_idx + 0.5
      ),
      fill      = "white",
      color     = "grey92",
      linewidth = 0.1
    ) +
    
    # CNV segments
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = genome_start_plot,
        xmax = genome_end_plot,
        ymin = plot_idx - 0.5,
        ymax = plot_idx + 0.5,
        fill = cnv_state
      ),
      color = NA
    ) +
    
    # Chromosome boundaries
    ggplot2::geom_vline(
      data      = chr_boundaries,
      ggplot2::aes(xintercept = chr_start),
      color     = "black",
      linewidth = 0.2,
      alpha     = 0.4
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = chr_labels$label_pos,
      labels = chr_labels$label,
      expand = c(0, 0),
      limits = c(0, genome_size)
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = y_label_data$mid_idx,
      labels = y_label_data$group_label,
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_fill_manual(
      values = state_colors,
      breaks = c("gain", "loss"), 
      name   = "CNV State"
    ) +
    
    ggplot2::labs(
      x     = "Chromosome",
      y     = NULL
    ) +
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_text(size = 7),
      axis.text.x     = ggplot2::element_text(size = 9),
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_rect(
        fill      = NA,
        color     = "black",
        linewidth = 0.3
      ),
      legend.position = if (show_legend) "bottom" else "none"
    )
  
  # ---- Boundary lines -----------------------------------------------------
  n_cells <- nrow(cell_order)
  
  if (!is.null(boundaries) && length(boundaries) > 0L &&
      length(boundaries) < n_cells / 2 + 50) {
    p_main <- p_main +
      ggplot2::geom_hline(
        yintercept = boundaries,
        color      = "black",
        linewidth  = 0.3,
        alpha      = 0.7
      )
  }
  
  # ---- Ideogram -----------------------------------------------------------
  p_ideo <- plot_ideogram(
    genome_structure = genome_structure,
    max_cell_id      = n_cells,
    arm_colors       = arm_colors
  ) +
    ggplot2::scale_x_continuous(
      breaks = chr_labels$label_pos,
      labels = chr_labels$label,
      expand = c(0, 0),
      limits = c(0, genome_size + 1e6)
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_blank()
    )
  
  # ---- Combine with patchwork ---------------------------------------------
  # ideogram_ratio controls relative height
  # e.g. 0.08 = ideogram takes 8% of total height
  patchwork::wrap_plots(
    p_ideo,
    p_main,
    ncol         = 1,
    heights      = c(ideogram_ratio, 1 - ideogram_ratio)
  )
}




# =============================================================================
# add_empty_cells
# Adds rows for cells with no CNV
# so they appear as empty lines in the plot
# Requires: cell_sizes with cell_name + cell_type
# =============================================================================

add_empty_cells <- function(
    cnv_mapped,
    cell_sizes,
    cell_name_col = "cell_name",
    cell_type_col = "cell_type") {
  
  cat("Adding empty cells...\n")
  
  # Cells already in cnv_mapped
  cells_with_cnv <- unique(
    cnv_mapped[[cell_type_coll]])
  
  cat("Cells with CNV:    ",
      length(cells_with_cnv), "\n")
  
  # All cells from cell_sizes
  all_cells <- cell_sizes[[cell_type_coll]]
  
  browser()
  cat("Total cells:       ",
      length(all_cells), "\n")
  
  # Cells with NO CNV
  empty_cells <- cell_sizes %>%
    dplyr::filter(
      !.data[[cell_name_col]] %in% cells_with_cnv
    )
  
  cat("Cells without CNV: ",
      nrow(empty_cells), "\n")
  
  if (nrow(empty_cells) == 0) {
    cat("No empty cells to add\n")
    return(cnv_mapped)
  }
  
  # Build empty rows
  # same columns as cnv_mapped
  # but NA for genomic coordinates
  # cnv_state = "none" for colour mapping
  empty_rows <- empty_cells %>%
    dplyr::select(
      !!cell_name_col := !!cell_name_col,
      !!cell_type_col := !!cell_type_col
    ) %>%
    dplyr::mutate(
      chr                = NA_character_,
      genome_start_plot  = NA_real_,
      genome_end_plot    = NA_real_,
      cnv_state          = "none",
      cnv_length_mb      = 0,
      start              = NA_real_,
      end                = NA_real_,
      cnv_length         = NA_real_,
      arm_class          = NA_character_
    )
  
  # Add any extra columns present in cnv_mapped
  # fill with NA
  missing_cols <- setdiff(
    colnames(cnv_mapped),
    colnames(empty_rows)
  )
  
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      empty_rows[[col]] <- NA
    }
  }
  
  # Bind — empty cells appended
  # prepare_cnv_plot will order them
  # within their cell_type group
  result <- dplyr::bind_rows(
    cnv_mapped,
    empty_rows %>%
      dplyr::select(
        dplyr::all_of(colnames(cnv_mapped))
      )
  )
  
  cat("Total rows after adding empty:",
      nrow(result), "\n")
  cat("Empty cell rows added:",
      nrow(empty_rows), "\n")
  
  return(result)
}



=======
library(ggplot2)
library(data.table)
library(dplyr)

#' Plot the density distribution of CNV lengths
#'
#' Creates a density plot of CNV lengths in megabases, optionally restricted
#' to a single CNV state such as gain or loss.
#'
#' @param dt A data frame containing at least the columns cnv_length_mb
#'   and optionally cnv_state.
#' @param state Optional character string specifying a CNV state to subset,
#'   such as "gain" or "loss". If NULL, all rows are used.
#' @param thresholds Numeric vector of CNV length thresholds to display as
#'   vertical dashed lines.
#' @param fill_color Fill color for the density polygon. Default is
#' @param title Plot title. Default is "CNV length distribution".
#'
#' @return A ggplot2 object.
plot_cnv_density <- function(
    dt,
    state = NULL,
    thresholds = c(5, 25, 50),
    fill_color = "grey40",
    title = "CNV length distribution"
) {
  
  # Subset efficiently (no copy)
  if (!is.null(state)) {
    dt_sub <- dt[dt$cnv_state == state,]
  } else {
    dt_sub <- dt
  }
  
  ggplot(dt_sub, aes(x = cnv_length_mb)) +
    geom_density(fill = fill_color, alpha = 0.4) +
    geom_vline(
      xintercept = thresholds,
      linetype = "dashed",
      alpha = 0.6
    ) +
    labs(
      title = title,
      x = "CNV length (Mb)",
      y = "Density"
    ) +
    theme_minimal()
}



#' Plot overall and state-specific CNV length distributions
#'
#' Builds three density plots showing CNV length distributions for all CNVs,
#' gain events, and loss events, then stacks them vertically.
#'
#' @param dt A data frame containing CNV length information, including
#'   cnv_length_mb and cnv_state.
#' @param thresholds Numeric vector of thresholds to display as dashed vertical
#'   lines in each panel. Default is c(5, 25, 50).
#'
#' @return A combined patchwork plot object.
plot_all_cnv_distributions <- function(
    dt,
    thresholds = c(5, 25, 50)
) {
  
  p_overall <- plot_cnv_density(
    dt = dt,
    state = NULL,
    thresholds = thresholds,
    fill_color = "grey40",
    title = "Overall CNV length distribution"
  )
  
  p_gain <- plot_cnv_density(
    dt = dt,
    state = "gain",
    thresholds = thresholds,
    fill_color = "steelblue",
    title = "Gain CNV length distribution"
  )
  
  p_loss <- plot_cnv_density(
    dt = dt,
    state = "loss",
    thresholds = thresholds,
    fill_color = "firebrick",
    title = "Loss CNV length distribution"
  )
  
  p_overall / p_gain / p_loss
}




#' Plot density distributions by level and type
#'
#' Creates a density plot for percentage values within a selected level,
#' grouped by type.
#'
#' @param level_name Character string specifying which level in
#'   plot_long$level to plot.
#' @param plot_long A long-format data frame containing at least the columns
#'   level, type, and percentage.
#' @param threshold Numeric threshold shown as a dashed vertical line.
#'
#' @return A ggplot2 object.
#'
#' @details
#' The function also computes the mean percentage per type, although the
#' current vertical reference line uses the supplied threshold value.
make_density_plot <- function(level_name,plot_long, threshold) {
  
  df_sub <- filter(plot_long, level == level_name)
  
  # compute means per type
  mean_df <- df_sub %>%
    group_by(type) %>%
    summarise(mean_value = mean(percentage, na.rm = TRUE),
              .groups = "drop")
  
  ggplot(df_sub,
         aes(x = percentage, fill = type, color = type)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    
    # vertical mean lines (like abline)
    geom_vline(data = mean_df,
               aes(xintercept = threshold),
               linetype = "dashed",
               linewidth = 1.2,
               show.legend = FALSE) +
    
    labs(
      title = level_name,
      x = "Percentage",
      y = "Density"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")
}





################################################################
## Functions to Process Information for conjoined heatmap#######
#################################################################
#' Prepare genome-wide chromosome and arm coordinates
#'
#' Builds chromosome-level cumulative coordinates and lookup tables for
#' chromosome arms, enabling mapping of chromosome-local intervals into a
#' genome-wide coordinate system.
#'
#' @param chromosome_arms A data frame containing chromosome arm annotation.
#'   It should include at least \code{chr}, \code{arm}, \code{arm_start},
#'   \code{arm_end}, and \code{arm_length}.
#'
#' @return A named list containing information the different fraction from the chromossome
#'
prepare_genome_structure <- function(chromosome_arms) {
  
  chromosome_lengths <- chromosome_arms %>%
    group_by(chr) %>%
    summarise(chr_length = sum(arm_length), .groups = "drop") %>%
    mutate(
      chr_length = as.numeric(chr_length),
      chr_num = suppressWarnings(as.numeric(gsub("chr", "", chr)))
    ) %>%
    arrange(chr_num, chr) %>%
    dplyr::select(-chr_num) %>%
    mutate(
      
      chr_start = lag(cumsum(chr_length), default = 0),
      chr_end   = chr_start + chr_length
    )
  
  # Precompute arm lookup tables (vectorized, no match needed later)
  arm_lookup <- chromosome_arms %>%
    left_join(chromosome_lengths, by = "chr") %>%
    mutate(
      genome_arm_start = chr_start + arm_start,
      genome_arm_end   = chr_start + arm_end
    )
  
  p_lookup <- arm_lookup %>%
    filter(arm == "p") %>%
    select(chr,
           p_genome_start = genome_arm_start,
           p_genome_end   = genome_arm_end)
  
  q_lookup <- arm_lookup %>%
    filter(arm == "q") %>%
    select(chr,
           q_genome_start = genome_arm_start,
           q_genome_end   = genome_arm_end)
  
  list(
    chromosome_lengths = chromosome_lengths,
    p_lookup = p_lookup,
    q_lookup = q_lookup,
    arm_lookup = arm_lookup
  )
}

#' Map CNV intervals to genome-wide coordinates
#'
#' Converts chromosome-level CNV coordinates into cumulative genome-wide
#' coordinates for visualization. The function joins chromosome and chromosome-arm
#' lookup tables, computes genome-wide start and end positions, optionally
#' expands large events to full chromosome or arm boundaries for plotting,
#' and encodes CNV state numerically.
#'
#' @param cnv_filtered A data frame of CNV events
#' @param genome_structure
#' @param threshold Optional numeric threshold used to expand plotted CNV
#'   intervals to whole-chromosome or chromosome-arm boundaries.
#' @return A data frame with additional columns
map_cnv_to_genome <- function(cnv_filtered,
                              genome_structure,
                              threshold = NULL,
                              arrange_df_cols) {
  
  chromosome_lengths <- genome_structure$chromosome_lengths
  p_lookup <- genome_structure$p_lookup
  q_lookup <- genome_structure$q_lookup
  
  
  cnv_mapped <- cnv_filtered %>%
    left_join(chromosome_lengths, by = "chr") %>%
    left_join(p_lookup, by = "chr") %>%
    left_join(q_lookup, by = "chr") %>%
    mutate(
      genome_start = chr_start + start,
      genome_end   = chr_start + end
    )
  
  # If threshold is NULL → no override logic
  if (is.null(threshold)) {
    
    cnv_mapped <- cnv_mapped %>%
      mutate(
        genome_start_plot = genome_start,
        genome_end_plot   = genome_end
      )
    
  } else {
    
    cnv_mapped <- cnv_mapped %>%
      mutate(
        genome_start_plot = case_when(
          whole_chromosome_gain > threshold |
            whole_chromosome_loss > threshold ~ chr_start,
          
          p_arm_gain > threshold |
            p_arm_loss > threshold ~ p_genome_start,
          
          q_arm_gain > threshold |
            q_arm_loss > threshold ~ q_genome_start,
          
          TRUE ~ genome_start
        ),
        genome_end_plot = case_when(
          whole_chromosome_gain > threshold |
            whole_chromosome_loss > threshold ~ chr_end,
          
          p_arm_gain > threshold |
            p_arm_loss > threshold ~ p_genome_end,
          
          q_arm_gain > threshold |
            q_arm_loss > threshold ~ q_genome_end,
          
          TRUE ~ genome_end
        )
      )
  }
 
  distinct_dfs_cols <- c("cell_name", arrange_df_cols)
  cell_order <- cnv_filtered %>%
    distinct(across(all_of(distinct_dfs_cols))) %>%
    arrange(all_of(arrange_df_cols))
  
  cnv_mapped <- cnv_mapped %>%
    mutate(
      cell_name = factor(cell_name, levels = cell_order$cell_name),
      cell_id = as.numeric(cell_name),
      cnv_state_numeric = case_when(
        cnv_state == "gain" ~ 1,
        cnv_state == "loss" ~ -1,
        TRUE ~ 0
      )
    ) %>%
    arrange(arrange_df_cols)
}

#' Plot a chromosome ideogram
#'
#' Creates a simple ideogram-style plot of chromosome arms aligned to genome-wide
#' coordinates.
#'
#' @param genome_structure Output from prepare_genome_structure().
#' @param max_cell_id Maximum cell index used to place the ideogram above the
#'   heatmap.
#' @param arm_colors Named vector of colors for chromosome arms.
#'
#' @return A ggplot2 object.
plot_ideogram <- function(genome_structure,
                          max_cell_id,
                          arm_colors = c("p"="#4DBBD5",
                                         "cen"="black",
                                         "q"="#E64B35")) {
  
  chromosome_lengths <- genome_structure$chromosome_lengths
  arm_lookup <- genome_structure$arm_lookup
  
  arm_plot <- arm_lookup %>%
    mutate(
      ymin = max_cell_id + 1,
      ymax = ymin + 1
    )
  
  ggplot(arm_plot) +
    geom_rect(aes(
      xmin = genome_arm_start,
      xmax = genome_arm_end,
      ymin = ymin,
      ymax = ymax,
      fill = arm
    ), color = NA) +
    scale_fill_manual(values = arm_colors, name = "Arm") +
    theme_void() +
    scale_x_continuous(
      breaks = chromosome_lengths$chr_start,
      labels = chromosome_lengths$chr,
      expand = c(0,0),
      limits = c(0, max(chromosome_lengths$chr_end) +1e6)
    )
}



#' Prepare CNV karyotype plot data
#'
#' Assigns plot_idx per cell ordered by grouping columns then CNV burden,
#' validates genome-wide coordinates, and prepares boundary line positions.
#'
#' @param cnv_mapped Data frame with one row per cell per CNV segment.
#'   Must contain: cell_name, genome_start_plot, genome_end_plot,
#'   cnv_state, cnv_length_mb, and all columns in grouping_cols.
#' @param genome_structure Output of prepare_genome_structure().
#' @param grouping_cols Character vector of columns defining cell groups
#'   for ordering and boundary lines. Default c("cell_type").
#'   Can include "embryo_id", "cell_type", or both.
#' @param state_colors Named character vector mapping cnv_state to colours.
#'
#' @return Named list with plot_data, chr_lengths, genome_size,
#'   boundary_lines, cell_order, grouping_cols.
prepare_cnv_plot <- function(
    cnv_mapped,
    genome_structure,
    grouping_cols = c("cell_type"),
    state_colors  = c(
      "gain" = "#E64B35",
      "loss" = "#4DBBD5"
    )
) {
  
  # ---- Input validation ---------------------------------------------------
  required_cols <- c(
    "cell_name", "chr",
    "genome_start_plot", "genome_end_plot",
    "cnv_state", "cnv_length_mb"
  )
  missing_cols <- setdiff(required_cols, colnames(cnv_mapped))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  missing_group <- setdiff(grouping_cols, colnames(cnv_mapped))
  if (length(missing_group) > 0L) {
    stop("grouping_cols not found in cnv_mapped: ",
         paste(missing_group, collapse = ", "))
  }
  
  unknown_states <- setdiff(
    unique(cnv_mapped$cnv_state), names(state_colors)
  )
  if (length(unknown_states) > 0L) {
    warning(
      "cnv_state values not in state_colors — will appear grey: ",
      paste(unknown_states, collapse = ", ")
    )
  }
  
  # ---- Extract genome structure -------------------------------------------
  chr_lengths <- genome_structure$chromosome_lengths
  genome_size <- sum(chr_lengths$chr_length)
  
  # ---- Compute cell ordering ----------------------------------------------
  # Order: grouping_cols → total CNV burden descending → cell_name
  # One row per cell for ordering purposes
  cell_order <- cnv_mapped |>
    dplyr::group_by(dplyr::across(dplyr::all_of(
      c(grouping_cols, "cell_name")
    ))) |>
    dplyr::summarise(
      total_cnv_mb = sum(cnv_length_mb, na.rm = TRUE),
      .groups      = "drop"
    ) |>
    dplyr::arrange(
      dplyr::across(dplyr::all_of(grouping_cols)),
      dplyr::desc(total_cnv_mb),
      cell_name
    ) |>
    dplyr::mutate(plot_idx = dplyr::row_number())
  
  message(sprintf(
    "Cell ordering: %d unique cells across %d group(s)",
    nrow(cell_order),
    dplyr::n_distinct(cell_order[, grouping_cols])
  ))
  
  # ---- Assign plot_idx to cnv_mapped --------------------------------------
  # Join back — every segment row for a cell gets same plot_idx
  plot_data <- cnv_mapped |>
    dplyr::select(-dplyr::any_of("plot_idx")) |>  # remove if exists
    dplyr::left_join(
      cell_order |> dplyr::select(cell_name, plot_idx),
      by = "cell_name"
    )
  
  # Sanity check — every cell should have exactly one plot_idx
  idx_check <- plot_data |>
    dplyr::group_by(cell_name) |>
    dplyr::summarise(
      n_idx = dplyr::n_distinct(plot_idx),
      .groups = "drop"
    ) |>
    dplyr::filter(n_idx > 1L)
  
  if (nrow(idx_check) > 0L) {
    stop(sprintf(
      "%d cell(s) have multiple plot_idx values — check cell_name uniqueness: %s",
      nrow(idx_check),
      paste(idx_check$cell_name[1:min(5, nrow(idx_check))], collapse = ", ")
    ))
  }
  
  # ---- Compute boundary lines ---------------------------------------------
  # Boundaries sit between groups defined by grouping_cols
  # Position = last plot_idx of each group + 0.5
  boundary_lines <- cell_order |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) |>
    dplyr::summarise(
      last_idx = max(plot_idx),
      .groups  = "drop"
    ) |>
    # Don't add boundary after last group
    dplyr::filter(last_idx < max(cell_order$plot_idx)) |>
    dplyr::pull(last_idx) + 0.5
  
  message(sprintf(
    "%d boundary line(s) computed between groups",
    length(boundary_lines)
  ))
  
  # ---- Coordinate validation ----------------------------------------------
  n_na_coords <- sum(
    is.na(plot_data$genome_start_plot) |
      is.na(plot_data$genome_end_plot)
  )
  if (n_na_coords > 0L) {
    warning(sprintf(
      "%d segments have NA genome coordinates — will not be plotted.",
      n_na_coords
    ))
  }
  
  out_of_range <- plot_data |>
    dplyr::filter(
      genome_start_plot < 0 |
        genome_end_plot > genome_size |
        genome_start_plot > genome_end_plot
    )
  
  if (nrow(out_of_range) > 0L) {
    warning(sprintf(
      "%d segments have out-of-range coordinates — check genome_structure.",
      nrow(out_of_range)
    ))
  }
  
  message(sprintf(
    "Plot data ready: %d segments across %d cells",
    nrow(plot_data),
    dplyr::n_distinct(plot_data$cell_name)
  ))
  
  list(
    plot_data      = plot_data,
    chr_lengths    = chr_lengths,
    genome_size    = genome_size,
    boundary_lines = boundary_lines,
    cell_order     = cell_order,
    grouping_cols  = grouping_cols,
    state_colors   = state_colors
  )
}


#' Plot CNV karyotype with ideogram
#'
#' Builds the main CNV heatmap and combines it with a chromosome
#' ideogram using patchwork. Ideogram height is fixed at ideogram_ratio
#' of the total plot height.
#'
#' @param panel_data Output of prepare_cnv_plot().
#' @param genome_structure Output of prepare_genome_structure().
#' @param title Plot title.
#' @param boundary_lines Numeric vector of y positions for group
#'   boundaries. If NULL uses boundary_lines from panel_data.
#' @param ideogram_ratio Fraction of total height for ideogram. Default 0.08.
#' @param arm_colors Named vector of colours for chromosome arms.
#' @param show_legend Logical. Default TRUE.
#'
#' @return A patchwork object combining ideogram + main plot.
plot_cnv_karyotype <- function(
    panel_data,
    genome_structure,
    boundary_lines  = NULL,
    ideogram_ratio  = 0.08,
    arm_colors      = c(
      "p"   = "#4DBBD5",
      "cen" = "black",
      "q"   = "#E64B35"
    ),
    show_legend     = TRUE,
    cell_sizes      = NULL  
) {
  
  plot_data      <- panel_data$plot_data
  chr_lengths    <- panel_data$chr_lengths
  genome_size    <- panel_data$genome_size
  state_colors   <- panel_data$state_colors
  cell_order     <- panel_data$cell_order
  
  # Use provided boundary_lines or fall back to prepared ones
  boundaries <- boundary_lines %||% panel_data$boundary_lines
  
  # ---- Chromosome labels --------------------------------------------------
  chr_labels <- chr_lengths |>
    dplyr::mutate(
      label_pos = chr_start + chr_length / 2,
      label     = gsub("chr", "", chr),
      label     = dplyr::case_when(
        label == "23" ~ "X",
        label == "24" ~ "Y",
        TRUE          ~ label
      )
    )
  
  chr_boundaries <- chr_lengths |>
    dplyr::filter(chr_start > 0)
  
  # ---- Y axis labels ------------------------------------------------------
  # One label per cell — use cell_order for correct positioning
  y_label_data <- panel_data$cell_order |>
    dplyr::group_by(dplyr::across(dplyr::all_of(panel_data$grouping_cols))) |>
    dplyr::summarise(
      mid_idx = mean(plot_idx),   # center of group for label position
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Combine grouping cols into one label if multiple
      group_label = apply(
        dplyr::pick(dplyr::all_of(panel_data$grouping_cols)),
        1,
        paste, collapse = " — "
      )
    )
  
  
  if (!is.null(cell_sizes)) {
    
    # Count CNV positive per group
    cnv_positive <- plot_data %>%
      dplyr::filter(cnv_state != "none") %>%
      dplyr::distinct(
        cell_name,
        dplyr::across(
          dplyr::all_of(panel_data$grouping_cols)
        )
      ) %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(panel_data$grouping_cols)
        )
      ) %>%
      dplyr::summarise(
        n_cnv = dplyr::n(),
        .groups = "drop"
      )
    
    # Total per group from cell_sizes
    total_per_group <- cell_sizes %>%
      dplyr::group_by(
        dplyr::across(
          dplyr::all_of(panel_data$grouping_cols)
        )
      ) %>%
      dplyr::summarise(
        n_total = dplyr::n(),
        .groups = "drop"
      )
    
    y_label_data <- y_label_data %>%
      dplyr::left_join(cnv_positive,
                       by = panel_data$grouping_cols) %>%
      dplyr::left_join(total_per_group,
                       by = panel_data$grouping_cols) %>%
      dplyr::mutate(
        pct_cnv   = round(100 * n_cnv / n_total, 1),
        pct_empty = round(100 - pct_cnv, 1),
        group_label = paste0(
          group_label, "\n",
          "CNV: ", pct_cnv, "% | ",
          "Empty: ", pct_empty, "%"
        )
      )
  }
  
  
  # ---- Main plot ----------------------------------------------------------
  p_main <- ggplot2::ggplot(plot_data) +
    
    # White background row per cell
    ggplot2::geom_rect(
      data = plot_data |>
        dplyr::distinct(cell_name, plot_idx),
      ggplot2::aes(
        xmin = 0,
        xmax = genome_size,
        ymin = plot_idx - 0.5,
        ymax = plot_idx + 0.5
      ),
      fill      = "white",
      color     = "grey92",
      linewidth = 0.1
    ) +
    
    # CNV segments
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = genome_start_plot,
        xmax = genome_end_plot,
        ymin = plot_idx - 0.5,
        ymax = plot_idx + 0.5,
        fill = cnv_state
      ),
      color = NA
    ) +
    
    # Chromosome boundaries
    ggplot2::geom_vline(
      data      = chr_boundaries,
      ggplot2::aes(xintercept = chr_start),
      color     = "black",
      linewidth = 0.2,
      alpha     = 0.4
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = chr_labels$label_pos,
      labels = chr_labels$label,
      expand = c(0, 0),
      limits = c(0, genome_size)
    ) +
    
    ggplot2::scale_y_continuous(
      breaks = y_label_data$mid_idx,
      labels = y_label_data$group_label,
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_fill_manual(
      values = state_colors,
      breaks = c("gain", "loss"), 
      name   = "CNV State"
    ) +
    
    ggplot2::labs(
      x     = "Chromosome",
      y     = NULL
    ) +
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_text(size = 7),
      axis.text.x     = ggplot2::element_text(size = 9),
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_rect(
        fill      = NA,
        color     = "black",
        linewidth = 0.3
      ),
      legend.position = if (show_legend) "bottom" else "none"
    )
  
  # ---- Boundary lines -----------------------------------------------------
  n_cells <- nrow(cell_order)
  
  if (!is.null(boundaries) && length(boundaries) > 0L &&
      length(boundaries) < n_cells / 2 + 50) {
    p_main <- p_main +
      ggplot2::geom_hline(
        yintercept = boundaries,
        color      = "black",
        linewidth  = 0.3,
        alpha      = 0.7
      )
  }
  
  # ---- Ideogram -----------------------------------------------------------
  p_ideo <- plot_ideogram(
    genome_structure = genome_structure,
    max_cell_id      = n_cells,
    arm_colors       = arm_colors
  ) +
    ggplot2::scale_x_continuous(
      breaks = chr_labels$label_pos,
      labels = chr_labels$label,
      expand = c(0, 0),
      limits = c(0, genome_size + 1e6)
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_blank()
    )
  
  # ---- Combine with patchwork ---------------------------------------------
  # ideogram_ratio controls relative height
  # e.g. 0.08 = ideogram takes 8% of total height
  patchwork::wrap_plots(
    p_ideo,
    p_main,
    ncol         = 1,
    heights      = c(ideogram_ratio, 1 - ideogram_ratio)
  )
}




# =============================================================================
# add_empty_cells
# Adds rows for cells with no CNV
# so they appear as empty lines in the plot
# Requires: cell_sizes with cell_name + cell_type
# =============================================================================

add_empty_cells <- function(
    cnv_mapped,
    cell_sizes,
    cell_name_col = "cell_name",
    cell_type_col = "cell_type") {
  
  cat("Adding empty cells...\n")
  
  # Cells already in cnv_mapped
  cells_with_cnv <- unique(
    cnv_mapped[[cell_type_coll]])
  
  cat("Cells with CNV:    ",
      length(cells_with_cnv), "\n")
  
  # All cells from cell_sizes
  all_cells <- cell_sizes[[cell_type_coll]]
  
  browser()
  cat("Total cells:       ",
      length(all_cells), "\n")
  
  # Cells with NO CNV
  empty_cells <- cell_sizes %>%
    dplyr::filter(
      !.data[[cell_name_col]] %in% cells_with_cnv
    )
  
  cat("Cells without CNV: ",
      nrow(empty_cells), "\n")
  
  if (nrow(empty_cells) == 0) {
    cat("No empty cells to add\n")
    return(cnv_mapped)
  }
  
  # Build empty rows
  # same columns as cnv_mapped
  # but NA for genomic coordinates
  # cnv_state = "none" for colour mapping
  empty_rows <- empty_cells %>%
    dplyr::select(
      !!cell_name_col := !!cell_name_col,
      !!cell_type_col := !!cell_type_col
    ) %>%
    dplyr::mutate(
      chr                = NA_character_,
      genome_start_plot  = NA_real_,
      genome_end_plot    = NA_real_,
      cnv_state          = "none",
      cnv_length_mb      = 0,
      start              = NA_real_,
      end                = NA_real_,
      cnv_length         = NA_real_,
      arm_class          = NA_character_
    )
  
  # Add any extra columns present in cnv_mapped
  # fill with NA
  missing_cols <- setdiff(
    colnames(cnv_mapped),
    colnames(empty_rows)
  )
  
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      empty_rows[[col]] <- NA
    }
  }
  
  # Bind — empty cells appended
  # prepare_cnv_plot will order them
  # within their cell_type group
  result <- dplyr::bind_rows(
    cnv_mapped,
    empty_rows %>%
      dplyr::select(
        dplyr::all_of(colnames(cnv_mapped))
      )
  )
  
  cat("Total rows after adding empty:",
      nrow(result), "\n")
  cat("Empty cell rows added:",
      nrow(empty_rows), "\n")
  
  return(result)
}



>>>>>>> f7a7a33 (feat: initial commit of CNV pipeline scripts)
