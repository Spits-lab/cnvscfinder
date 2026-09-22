source("/rhea/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/round2_pipeline.R")
source("/rhea/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_processing.R")
source("/rhea/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/cnv_scoring.R")
source("/rhea/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/pipeline.R")
source("/rhea/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/R/infercnv.R")


BASE <- "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"


results_r3 <- run_round2_cnv_pipeline(
  infercnv_obj_path = file.path(
    BASE,
    "infercnv_results/VUB04/round3/run.final.infercnv_obj"
  ),
  metadata          =metadata_3round.rds,
  chromosome_arms   = readRDS(file.path(
    BASE, "data/hg38_chromosome_arms.rds")),
  cell_col          = "cell_name",
  group_col         = "cell_type",
  ref_label         = "diploid_ref",
  k_discrete        = 1.0,
  min_segment_mb    = 15,
  max_gap           = 100000L,
  k_fre             = 1.5,
  sensitivity_floor = 20,
  min_overlap       = 0.80,
  workdir           = file.path(
    BASE, "cnv_results/VUB04/round3")
)
