# compare_overlap_methods.R
# Re-runs Block4 (clustering + scoring) only, from an already-saved
# pipeline_results.rds, comparing overlap_method = "adaptive_floor" and
# "reciprocal" against the original saved "adaptive" result — to check
# whether the overlap method explains the smaller-than-expected result
# count. Does not touch Block1/Block2 (no inferCNV re-run).

source("R/cnv_processing.R")
source("R/cnv_annotation.R")
source("R/cnv_scoring.R")

pipeline_results <- readRDS("cnv_results/VUB04/round1_3/pipeline_results.rds")
cnv_annotated    <- pipeline_results[["all_results"]][["block3"]][["cnv_annotated"]]
cell_sizes       <- pipeline_results[["all_results"]][["block3"]][["cell_sizes"]]
original_scored  <- pipeline_results[["all_results"]][["block4"]][["scored_events"]]
chromosome_arms  <- readRDS("data/hg38_chromosome_arms.rds")

message("=== Original (adaptive, from saved run) ===")
message("n rows: ", nrow(original_scored))

# ---- adaptive_floor -------------------------------------------------------
message("\n=== Running adaptive_floor ===")
clustered_af <- run_cnv_locus_analysis(
  cnv_annotated,
  by                    = "cell_type",
  overlap_method        = "adaptive_floor",
  min_overlap           = 0.70,
  sample_col            = "cell_type",
  cell_col              = "cell_name",
  range                 = 0.2,
  sensitivity_floor_mb  = 20,
  max_mb                = 100
)

scored_af <- score_cnv_clusters(
  summary_df            = clustered_af$cnv_locus_summary,
  clustered_events      = clustered_af$clustered_events,
  cell_sizes            = cell_sizes,
  by_union              = "cell_type",
  chromosome_arms       = chromosome_arms,
  min_required_cells    = 3,
  round_fun             = ceiling,
  k                     = 1.4,
  sensitivity_floor_mb  = 20,
  p_arm_permission      = 60,
  q_arm_permission      = 60,
  whole_chr_permission  = 65
)

scored_af <- deduplicate_cnv_cells(scored_af, sample_col = "cell_type", cell_col = "cell_name")
message("n rows: ", nrow(scored_af))

# ---- reciprocal -------------------------------------------------------
message("\n=== Running reciprocal ===")
clustered_rec <- run_cnv_locus_analysis(
  cnv_annotated,
  by             = "cell_type",
  overlap_method = "reciprocal",
  min_overlap    = 0.70,
  sample_col     = "cell_type",
  cell_col       = "cell_name"
)

scored_rec <- score_cnv_clusters(
  summary_df            = clustered_rec$cnv_locus_summary,
  clustered_events      = clustered_rec$clustered_events,
  cell_sizes            = cell_sizes,
  by_union              = "cell_type",
  chromosome_arms       = chromosome_arms,
  min_required_cells    = 3,
  round_fun             = ceiling,
  k                     = 1.4,
  sensitivity_floor_mb  = 20,
  p_arm_permission      = 60,
  q_arm_permission      = 60,
  whole_chr_permission  = 65
)

scored_rec <- deduplicate_cnv_cells(scored_rec, sample_col = "cell_type", cell_col = "cell_name")
message("n rows: ", nrow(scored_rec))

# ---- Summary -------------------------------------------------------
message("\n=== Summary ===")
message(sprintf("adaptive (original saved run): %d rows", nrow(original_scored)))
message(sprintf("adaptive_floor:                 %d rows", nrow(scored_af)))
message(sprintf("reciprocal:                     %d rows", nrow(scored_rec)))

out_path <- "cnv_results/VUB04/round1_3/overlap_method_comparison.rds"
saveRDS(list(adaptive_floor = scored_af, reciprocal = scored_rec), out_path)
message("\nSaved comparison to: ", out_path)
