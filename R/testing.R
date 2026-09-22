source("R/cnv_processing.R")
source("R/cnv_annotation.R")
source("R/cnv_scoring.R")
source("R/pipeline.R")
source("R/infercnv.R")

pipeline_results <- readRDS("cnv_results/VUB04/round1_3/pipeline_results.rds")
cnv_annotated <- pipeline_results[["all_results"]][["block3"]][["cnv_annotated"]]
cell_sizes  <- pipeline_results[["all_results"]][["block3"]][["cell_sizes"]]
chromosome_arms <- readRDS("data/hg38_chromosome_arms.rds")

clustered_events <- run_cnv_locus_analysis(
      cnv_annotated,
      by             = "cell_type",
      overlap_method = "reciprocal",
      min_overlap     = 0.75,
      sample_col     = "cell_type",
      cell_col       = "cell_name",
      range                = 0.2,
      sensitivity_floor_mb = 20,
      max_mb               = 100
    )
    
    scored_events <- score_cnv_clusters(
      summary_df                  = clustered_events$cnv_locus_summary,
      clustered_events            = clustered_events$clustered_events,
      cell_sizes                  = cell_sizes,
      by_union                    = "cell_type",
      chromosome_arms             = chromosome_arms,
      min_required_cells           = 3,
      round_fun                   = ceiling,
      k = 1.5,
      sensitivity_floor_mb = 20,
      p_arm_permission     = 60,
      q_arm_permission     = 60,  
      whole_chr_permission = 65  
    )

    scored_events <-  deduplicate_cnv_cells(scored_events,
                               sample_col   = "cell_type",
                               cell_col     = "cell_name")
    

    saveRDS(list(clustered_events,scored_events), "cnv_results/VUB04/round1_3/scored_events_reciprocal.rds")