# =============================================================================
# 11_new_TF_cluster_stats.R
# Cluster-stratified patient-level tests for the five TFs perturbed in
# 10_additional_TF_perturbations.py (JUND, KLF6, STAT1, ZNF281, TAL1),
# mirroring 09_additional_TF_cluster_stats.R's pattern (one lm(cyto_delta ~
# group) per NK subtype cluster, not pooled) for this second batch of TFs.
#
# These five were chosen from 08_'s output specifically for an ATRXdel-vs-
# other-groups signal into RUNX3 (JUND/KLF6/ZNF281), general cytotoxicity-
# module strength (STAT1), or broad direct cytotoxicity-gene connectivity
# that bypasses RUNX3/TBX21 entirely (TAL1) -- see
# 10_additional_TF_perturbations.py's docstring for the specific numbers.
#
# Uses the same helpers as 04_/09_ (sourced from pseudobulk_helpers.R, not
# redefined here): run_perturb_analysis_by_cluster() for the module-score
# panels below, run_single_gene_analysis() (pooled across clusters -- no
# cluster-stratified single-gene helper exists yet) for the RUNX3/TBX21
# readout chains.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(emmeans); library(ggplot2); library(patchwork)
})

source("/home/rstudio/GRN/pseudobulk_helpers.R")

id_cols <- c("V1", "cell", "group", "cluster", "patient")
FIG_DIR <- "/home/rstudio/mnt_out/NK_project/figures/"

# =============================================================================
# (1) module-score results, per cluster, for all five TFs x KO/OE
# =============================================================================
res_jund_ko_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_JUND_KO.csv",   "JUND KO",   id_cols)
res_jund_oe_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_JUND_OE.csv",   "JUND OE",   id_cols)
res_klf6_ko_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_KLF6_KO.csv",   "KLF6 KO",   id_cols)
res_klf6_oe_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_KLF6_OE.csv",   "KLF6 OE",   id_cols)
res_stat1_ko_cl  <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_STAT1_KO.csv",  "STAT1 KO",  id_cols)
res_stat1_oe_cl  <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_STAT1_OE.csv",  "STAT1 OE",  id_cols)
res_znf281_ko_cl <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_ZNF281_KO.csv", "ZNF281 KO", id_cols)
res_znf281_oe_cl <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_ZNF281_OE.csv", "ZNF281 OE", id_cols)
res_tal1_ko_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_TAL1_KO.csv",   "TAL1 KO",   id_cols)
res_tal1_oe_cl   <- run_perturb_analysis_by_cluster("/home/rstudio/GRN/perturb_deltas_TAL1_OE.csv",   "TAL1 OE",   id_cols)

all_results <- list(
  JUND   = list(KO = res_jund_ko_cl,   OE = res_jund_oe_cl),
  KLF6   = list(KO = res_klf6_ko_cl,   OE = res_klf6_oe_cl),
  STAT1  = list(KO = res_stat1_ko_cl,  OE = res_stat1_oe_cl),
  ZNF281 = list(KO = res_znf281_ko_cl, OE = res_znf281_oe_cl),
  TAL1   = list(KO = res_tal1_ko_cl,   OE = res_tal1_oe_cl)
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

  out_path <- paste0(FIG_DIR, "New_TF_KO_OE_",
                      gsub("[^A-Za-z0-9_]", "-", cl), ".pdf")
  ggsave(out_path, combined_cl, width = 16, height = 14)  # up to 10 panels @ ncol=4 -> 3 rows
  message("Wrote ", out_path)
}

# =============================================================================
# (2) readout-chain results (pooled across clusters -- JUND/STAT1 -> RUNX3 +
#     TBX21, KLF6/ZNF281 -> RUNX3 only, TAL1 has no gene-level readout file)
# =============================================================================
res_jund_ko_runx3   <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_JUND_KO_RUNX3.csv",   "JUND KO→RUNX3",   "RUNX3", id_cols)
res_jund_oe_runx3   <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_JUND_OE_RUNX3.csv",   "JUND OE→RUNX3",   "RUNX3", id_cols)
res_jund_ko_tbx21   <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_JUND_KO_TBX21.csv",   "JUND KO→TBX21",   "TBX21", id_cols)
res_jund_oe_tbx21   <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_JUND_OE_TBX21.csv",   "JUND OE→TBX21",   "TBX21", id_cols)
res_klf6_ko_runx3   <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_KLF6_KO_RUNX3.csv",   "KLF6 KO→RUNX3",   "RUNX3", id_cols)
res_klf6_oe_runx3   <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_KLF6_OE_RUNX3.csv",   "KLF6 OE→RUNX3",   "RUNX3", id_cols)
res_stat1_ko_runx3  <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_STAT1_KO_RUNX3.csv",  "STAT1 KO→RUNX3",  "RUNX3", id_cols)
res_stat1_oe_runx3  <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_STAT1_OE_RUNX3.csv",  "STAT1 OE→RUNX3",  "RUNX3", id_cols)
res_stat1_ko_tbx21  <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_STAT1_KO_TBX21.csv",  "STAT1 KO→TBX21",  "TBX21", id_cols)
res_stat1_oe_tbx21  <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_STAT1_OE_TBX21.csv",  "STAT1 OE→TBX21",  "TBX21", id_cols)
res_znf281_ko_runx3 <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_ZNF281_KO_RUNX3.csv", "ZNF281 KO→RUNX3", "RUNX3", id_cols)
res_znf281_oe_runx3 <- run_single_gene_analysis("/home/rstudio/GRN/perturb_deltas_ZNF281_OE_RUNX3.csv", "ZNF281 OE→RUNX3", "RUNX3", id_cols)

readout_block <- wrap_plots(
  res_jund_ko_runx3$plot,   res_jund_oe_runx3$plot,
  res_jund_ko_tbx21$plot,   res_jund_oe_tbx21$plot,
  res_klf6_ko_runx3$plot,   res_klf6_oe_runx3$plot,
  res_stat1_ko_runx3$plot,  res_stat1_oe_runx3$plot,
  res_stat1_ko_tbx21$plot,  res_stat1_oe_tbx21$plot,
  res_znf281_ko_runx3$plot, res_znf281_oe_runx3$plot,
  ncol = 4
) + plot_annotation(title = "Readout-chain response (JUND/KLF6/STAT1/ZNF281 → RUNX3/TBX21)") &
  theme(legend.position = "bottom")

ggsave(paste0(FIG_DIR, "New_TF_readout_chains.pdf"), readout_block, width = 16, height = 12)
message("Wrote ", paste0(FIG_DIR, "New_TF_readout_chains.pdf"))

message("[11_new_TF_cluster_stats] done")
