library(infercnv)


infercnv_res <- readRDS("/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/cnv_results/nusa_gastruloids/iterations/run1/iter_3/infercnv/run.final.infercnv_obj")

infercnv::plot_cnv( infercnv_obj    = infercnv_res, out_dir         = "/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/cnv_results/nusa_gastruloids/iterations/run1/iter_3/infercnv/plots",output_filename = "infercnv_hmm_plot",output_format   = "pdf",x.range         = "auto",x.center        = 1, title           = "inferCNV HMM",color_safe_pal  = FALSE)