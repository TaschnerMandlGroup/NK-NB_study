# =============================================================================
# 09_additional_TF_cluster_stats.R
# Cluster-stratified patient-level tests for the five additional TFs
# perturbed in 08_top_TF_perturbations.py (MEF2C, KLF2, PRDM1, FOS, POU2F2),
# mirroring 04_pseudobulk_stats.R's "Part V" section (written there for
# TBX21/RUNX3 only) -- one lm(cyto_delta ~ group) per NK subtype cluster
# instead of pooling all subtypes together, so ATRXdel-vs-other-group
# contrasts can be read out per cluster.
#
# Uses the same run_perturb_analysis_by_cluster() helper as 04_'s Part V
# (sourced from pseudobulk_helpers.R, not redefined here).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(emmeans); library(ggplot2); library(patchwork)
})

source("/home/rstudio/GRN/pseudobulk_helpers.R")

id_cols <- c("V1", "cell", "group", "cluster", "patient")
FIG_DIR <- "/home/rstudio/mnt_out/NK_project/figures/"

# ---- module-score results, per cluster, for all five TFs x KO/OE ----------
res_mef2c_ko_cl  <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_MEF2C_KO.csv",  "MEF2C KO",  id_cols)
res_mef2c_oe_cl  <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_MEF2C_OE.csv",  "MEF2C OE",  id_cols)
res_klf2_ko_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_KLF2_KO.csv",   "KLF2 KO",   id_cols)
res_klf2_oe_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_KLF2_OE.csv",   "KLF2 OE",   id_cols)
res_prdm1_ko_cl  <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_PRDM1_KO.csv",  "PRDM1 KO",  id_cols)
res_prdm1_oe_cl  <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_PRDM1_OE.csv",  "PRDM1 OE",  id_cols)
res_fos_ko_cl    <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_FOS_KO.csv",    "FOS KO",    id_cols)
res_fos_oe_cl    <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_FOS_OE.csv",    "FOS OE",    id_cols)
res_pou2f2_ko_cl <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_POU2F2_KO.csv", "POU2F2 KO", id_cols)
res_pou2f2_oe_cl <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_POU2F2_OE.csv", "POU2F2 OE", id_cols)

all_results <- list(
  MEF2C  = list(KO = res_mef2c_ko_cl,  OE = res_mef2c_oe_cl),
  KLF2   = list(KO = res_klf2_ko_cl,   OE = res_klf2_oe_cl),
  PRDM1  = list(KO = res_prdm1_ko_cl,  OE = res_prdm1_oe_cl),
  FOS    = list(KO = res_fos_ko_cl,    OE = res_fos_oe_cl),
  POU2F2 = list(KO = res_pou2f2_ko_cl, OE = res_pou2f2_oe_cl)
)

clusters_seen <- unique(unlist(lapply(all_results, function(tf) unlist(lapply(tf, names)))))

# ---- one combined PDF per cluster: 5 TFs x KO/OE = up to 10 panels ---------
for (cl in clusters_seen) {
  panels <- Filter(Negate(is.null), unlist(lapply(names(all_results), function(tf) {
    list(all_results[[tf]]$KO[[cl]]$plot, all_results[[tf]]$OE[[cl]]$plot)
  }), recursive = FALSE))
  if (length(panels) == 0) next

  combined_cl <- wrap_plots(panels, ncol = 4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  out_path <- paste0(FIG_DIR, "Additional_TF_KO_OE_",
                      gsub("[^A-Za-z0-9_]", "-", cl), ".pdf")
  ggsave(out_path, combined_cl, width = 16, height = 14)  # up to 10 panels @ ncol=4 -> 3 rows
  message("Wrote ", out_path)
}

message("[09_additional_TF_cluster_stats] done")
