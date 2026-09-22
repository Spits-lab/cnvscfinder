library(infercnv)

infercnv_res <- readRDS("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/infercnv_results/VUB04/round7/run.final.infercnv_obj")

infercnv::plot_cnv( infercnv_obj    = infercnv_res, out_dir         = "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/infercnv_results/VUB04/plots/round7",output_filename = "infer
cnv_hmm_plot",output_format   = "pdf",x.range         = "auto",x.center        = 1, title           = "inferCNV HMM",color_safe_pal  = FALSE)