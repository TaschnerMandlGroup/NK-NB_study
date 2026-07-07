# =============================================================================
# 04_pseudobulk_stats.R
# Cross-group test of in silico TF-perturbation effects on the cytotoxicity
# module, aggregated to ONE value per patient before testing.
#
# Covers four perturbations: TBX21 KO/OE and RUNX3 KO/OE, using one shared
# helper so all four are processed identically, then combined into a single
# 4-panel comparison figure.
#
# Why patient-level: cells from the same donor are not independent. Testing at
# the cell level (pseudoreplication) inflates the FDR and can manufacture
# "significant" differences where none exist -- pseudobulk aggregation to the
# biological replicate is the accepted fix.
#   Squair et al., Nat Commun 2021 (PMID 34584091)
#   Zimmerman et al., Nat Commun 2021; sc-best-practices (DGE chapter)
#
# With 3-5 patients/group (~12-20 donors) a patient-level linear model is
# appropriate and honest for the rebuttal.
#
# patchwork (panel composition): Pedersen, T.L. -- https://patchwork.data-imaging.dev/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(emmeans); library(ggplot2); library(patchwork)
})

id_cols <- c("V1", "cell", "group", "patient")

# ---- helper: pseudobulk -> lm -> emmeans -> pairs -> ggplot panel ----------
run_perturb_analysis <- function(csv_path, label, id_cols) {
  d <- fread(csv_path)
  cyto_cols <- setdiff(colnames(d), id_cols)
  cyto_cols <- cyto_cols[!cyto_cols %in% c("group", "patient")]
  message(label, " cytotoxicity genes: ", paste(cyto_cols, collapse = ", "))
  
  d[, cyto_delta := rowMeans(.SD, na.rm = TRUE), .SDcols = cyto_cols]
  
  # ---- pseudobulk: one value per patient -----------------------------------
  pb <- d[, .(cyto_delta = mean(cyto_delta, na.rm = TRUE),
              n_cells    = .N),
          by = .(patient, group)]
  pb[, group := relevel(factor(group), ref = "control")]
  print(pb[order(group, patient)])
  
  # ---- patient-level model --------------------------------------------------
  fit <- lm(cyto_delta ~ group, data = pb)
  
  cat("\n=== ", label, " ===\n")
  cat("--- lm(cyto_delta ~ group) ---\n"); print(summary(fit))
  cat("--- per-group estimated marginal means (vs 0) ---\n"); print(emmeans(fit, ~ group))
  cat("--- pairwise genotype contrasts (Tukey) ---\n"); print(pairs(emmeans(fit, ~ group)))
  
  # ---- panel for the combined figure ---------------------------------------
  p <- ggplot(pb, aes(group, cyto_delta)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(size = n_cells), width = 0.12, alpha = 0.7) +
    labs(title = label, x = NULL,
         y = "\u0394 cytotoxicity module (patient mean)",
         size = "cells/patient") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  list(pb = pb, fit = fit, plot = p)
}

# ---- run all four perturbations --------------------------------------------
res_tbx21_ko <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas.csv",          "TBX21 KO", id_cols)
res_tbx21_oe <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_OE.csv",       "TBX21 OE", id_cols)
res_runx3_ko <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_KO.csv", "RUNX3 KO", id_cols)
res_runx3_oe <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_OE.csv", "RUNX3 OE", id_cols)

# ---- combined 4-panel figure -------------------------------------------------
combined <- (res_tbx21_ko$plot | res_tbx21_oe$plot) /
  (res_runx3_ko$plot | res_runx3_oe$plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("/home/rstudio/mnt_out/NK_project/figures/TBX21_RUNX3_KO_OE_combined.pdf", combined, width = 9, height = 8)
message("Wrote TBX21_RUNX3_KO_OE_combined.pdf")


# =============================================================================
# RUNX3 -> TBX21: does RUNX3 perturbation shift TBX21 expression itself?
# Note: the base GRN showed NO direct TBX21<->RUNX3 edge in any group (raw
# links_dict check), so any signal here reflects an indirect, multi-step
# propagated effect (n_propagation=3), not a direct regulatory link.
# Single-gene readout (TBX21), not a multi-gene module, so no rowMeans step.
# =============================================================================

run_single_gene_analysis <- function(csv_path, label, gene_col, id_cols) {
  d <- fread(csv_path)
  setnames(d, gene_col, "gene_delta")
  
  pb <- d[, .(gene_delta = mean(gene_delta, na.rm = TRUE),
              n_cells    = .N),
          by = .(patient, group)]
  pb[, group := relevel(factor(group), ref = "control")]
  print(pb[order(group, patient)])
  
  fit <- lm(gene_delta ~ group, data = pb)
  
  cat("\n=== ", label, " ===\n")
  cat("--- lm(gene_delta ~ group) ---\n"); print(summary(fit))
  cat("--- per-group estimated marginal means (vs 0) ---\n"); print(emmeans(fit, ~ group))
  cat("--- pairwise genotype contrasts (Tukey) ---\n"); print(pairs(emmeans(fit, ~ group)))
  
  p <- ggplot(pb, aes(group, gene_delta)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(size = n_cells), width = 0.12, alpha = 0.7) +
    labs(title = label, x = NULL,
         y = "\u0394 TBX21 expression (patient mean)",
         size = "cells/patient") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  list(pb = pb, fit = fit, plot = p)
}

res_runx3_ko_tbx21 <- run_single_gene_analysis(
  "/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_KO_TBX21.csv", "RUNX3 KO \u2192 TBX21", "TBX21", id_cols)
res_runx3_oe_tbx21 <- run_single_gene_analysis(
  "/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_OE_TBX21.csv", "RUNX3 OE \u2192 TBX21", "TBX21", id_cols)

combined_tbx21_link <- (res_runx3_ko_tbx21$plot | res_runx3_oe_tbx21$plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("/home/rstudio/mnt_out/NK_project/figures/RUNX3_to_TBX21_link.pdf", combined_tbx21_link, width = 7, height = 4)
message("Wrote RUNX3_to_TBX21_link.pdf")

# =============================================================================
# part III: top TFs
# =============================================================================

# =============================================================================
# MEF2C KO/OE (edges to both TBX21 and RUNX3) and KLF2 KO/OE (edge to RUNX3)
# =============================================================================

res_mef2c_ko <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_KO.csv", "MEF2C KO", id_cols)
res_mef2c_oe <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_OE.csv", "MEF2C OE", id_cols)
res_klf2_ko  <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_KO.csv",  "KLF2 KO",  id_cols)
res_klf2_oe  <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_OE.csv",  "KLF2 OE",  id_cols)

res_mef2c_ko_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_KO_TBX21.csv", "MEF2C KO -> TBX21", "TBX21", id_cols)
res_mef2c_oe_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_OE_TBX21.csv", "MEF2C OE -> TBX21", "TBX21", id_cols)
res_mef2c_ko_runx3 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_KO_RUNX3.csv", "MEF2C KO -> RUNX3", "RUNX3", id_cols)
res_mef2c_oe_runx3 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_OE_RUNX3.csv", "MEF2C OE -> RUNX3", "RUNX3", id_cols)
res_klf2_ko_runx3  <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_KO_RUNX3.csv",  "KLF2 KO -> RUNX3",  "RUNX3", id_cols)
res_klf2_oe_runx3  <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_OE_RUNX3.csv",  "KLF2 OE -> RUNX3",  "RUNX3", id_cols)

# ---- combined figure: module scores, MEF2C + KLF2, KO vs OE ----------------
combined_mef2c_klf2 <- (res_mef2c_ko$plot | res_mef2c_oe$plot) /
  (res_klf2_ko$plot  | res_klf2_oe$plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("/home/rstudio/mnt_out/NK_project/figures/MEF2C_KLF2_KO_OE_combined.pdf", combined_mef2c_klf2, width = 9, height = 8)
message("Wrote MEF2C_KLF2_KO_OE_combined.pdf")

# =============================================================================
# ===== Part IV: top TF in ATRXdel ============================================
# =============================================================================

# =============================================================================
# MEGA combined figure: all 6 TFs (module scores) + all readout chains
# Module scores: TBX21, RUNX3, MEF2C, KLF2, POU2F2, PRDM1 x KO/OE (12 panels)
# Readout chains: RUNX3->TBX21, MEF2C->TBX21, MEF2C->RUNX3, KLF2->RUNX3,
#                 PRDM1->RUNX3 x KO/OE (10 panels)
# Assembled as one PDF using patchwork::wrap_plots(), grouped into two
# labeled blocks (module scores on top, readout chains below).
#   patchwork: Pedersen, T.L. -- https://patchwork.data-imaging.dev/
#   wrap_plots() docs: https://patchwork.data-imaging.dev/reference/wrap_plots.html
#   plot_annotation() docs: https://patchwork.data-imaging.dev/reference/plot_annotation.html
# =============================================================================

# ---- module-score results for all 6 TFs ------------------------------------
res_tbx21_ko   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas.csv",           "TBX21 KO",  id_cols)
res_tbx21_oe   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_OE.csv",        "TBX21 OE",  id_cols)
res_runx3_ko   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_KO.csv",  "RUNX3 KO",  id_cols)
res_runx3_oe   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_OE.csv",  "RUNX3 OE",  id_cols)
res_mef2c_ko   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_KO.csv",  "MEF2C KO",  id_cols)
res_mef2c_oe   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_OE.csv",  "MEF2C OE",  id_cols)
res_klf2_ko    <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_KO.csv",   "KLF2 KO",   id_cols)
res_klf2_oe    <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_OE.csv",   "KLF2 OE",   id_cols)
res_pou2f2_ko  <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_POU2F2_KO.csv", "POU2F2 KO", id_cols)
res_pou2f2_oe  <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_POU2F2_OE.csv", "POU2F2 OE", id_cols)
res_prdm1_ko   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_PRDM1_KO.csv",  "PRDM1 KO",  id_cols)
res_prdm1_oe   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_PRDM1_OE.csv",  "PRDM1 OE",  id_cols)
res_fos_ko   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_FOS_KO.csv",  "FOS KO",  id_cols)
res_fos_oe   <- run_perturb_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_FOS_OE.csv",  "FOS OE",  id_cols)

# ---- readout-chain results (5 relationships x KO/OE) -----------------------
res_runx3_ko_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_KO_TBX21.csv", "RUNX3 KO\u2192TBX21", "TBX21", id_cols)
res_runx3_oe_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_RUNX3_OE_TBX21.csv", "RUNX3 OE\u2192TBX21", "TBX21", id_cols)
res_mef2c_ko_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_KO_TBX21.csv", "MEF2C KO\u2192TBX21", "TBX21", id_cols)
res_mef2c_oe_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_OE_TBX21.csv", "MEF2C OE\u2192TBX21", "TBX21", id_cols)
res_mef2c_ko_runx3 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_KO_RUNX3.csv", "MEF2C KO\u2192RUNX3", "RUNX3", id_cols)
res_mef2c_oe_runx3 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_MEF2C_OE_RUNX3.csv", "MEF2C OE\u2192RUNX3", "RUNX3", id_cols)
res_klf2_ko_runx3  <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_KO_RUNX3.csv",  "KLF2 KO\u2192RUNX3",  "RUNX3", id_cols)
res_klf2_oe_runx3  <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_KLF2_OE_RUNX3.csv",  "KLF2 OE\u2192RUNX3",  "RUNX3", id_cols)
res_prdm1_ko_runx3 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_PRDM1_KO_RUNX3.csv", "PRDM1 KO\u2192RUNX3", "RUNX3", id_cols)
res_prdm1_oe_runx3 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_PRDM1_OE_RUNX3.csv", "PRDM1 OE\u2192RUNX3", "RUNX3", id_cols)
res_fos_ko_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_FOS_KO_TBX21.csv", "FOS KO\u2192TBX21", "TBX21", id_cols)
res_fos_oe_tbx21 <- run_single_gene_analysis("/home/rstudio/NK-NB_study/code/GRN/perturb_deltas_FOS_OE_TBX21.csv", "FOS OE\u2192TBX21", "TBX21", id_cols)


# ---- assemble module-score block (12 panels, 6 TFs x KO/OE, 4 cols x 3 rows)
module_block <- wrap_plots(
  res_tbx21_ko$plot,  res_tbx21_oe$plot,  res_runx3_ko$plot,  res_runx3_oe$plot,
  res_mef2c_ko$plot,  res_mef2c_oe$plot,  res_klf2_ko$plot,   res_klf2_oe$plot,
  res_pou2f2_ko$plot, res_pou2f2_oe$plot, res_prdm1_ko$plot,  res_prdm1_oe$plot,
  res_fos_ko$plot,res_fos_oe$plot,
  ncol = 4
) + plot_annotation(title = "Module-score response (all 6 TFs)") &
  theme(legend.position = "bottom")

# ---- assemble readout-chain block (10 panels, 5 relationships x KO/OE) -----
readout_block <- wrap_plots(
  res_runx3_ko_tbx21$plot, res_runx3_oe_tbx21$plot,
  res_mef2c_ko_tbx21$plot, res_mef2c_oe_tbx21$plot,
  res_mef2c_ko_runx3$plot, res_mef2c_oe_runx3$plot,
  res_klf2_ko_runx3$plot,  res_klf2_oe_runx3$plot,
  res_prdm1_ko_runx3$plot, res_prdm1_oe_runx3$plot,
  res_fos_ko_tbx21$plot, res_fos_oe_tbx21$plot,
  ncol = 4
) + plot_annotation(title = "Readout-chain response (TF \u2192 TBX21/RUNX3)") &
  theme(legend.position = "bottom")

# ---- stack both blocks into one PDF, collecting all legends into one -------
mega_combined <- (module_block / readout_block) +
  plot_layout(guides = "collect", heights = c(3, 2.5)) &
  theme(legend.position = "bottom")

ggsave("/home/rstudio/mnt_out/NK_project/figures/ALL_TF_perturbation_mega_figure.pdf", mega_combined, width = 16, height = 22)
message("Wrote ALL_TF_perturbation_mega_figure.pdf")


# =============================================================================
# TF-TF regulatory edges: network diagram + bar chart of edge strengths
# Source data: tf_tf_edges_all_groups.csv (all 30 directed TF pairs x 4 groups,
# queried from links.filtered_links -- the network's top-2000/p<0.001 edges).
#   igraph: Csardi & Nepusz (2006), InterJournal Complex Systems, 1695
#   ggraph: Pedersen, T.L. -- https://ggraph.data-imaging.dev/
# =============================================================================

suppressPackageStartupMessages({
  library(tidygraph); library(ggraph); library(igraph)
})

edges <- fread("/home/rstudio/NK-NB_study/code/GRN/tf_tf_edges_all_groups.csv")
found_edges <- edges[found == TRUE]

# ---- network diagram, faceted by group --------------------------------------
graph_list <- list()
for (grp in unique(edges$group)) {
  grp_edges <- found_edges[group == grp]
  if (nrow(grp_edges) == 0) next
  g <- tbl_graph(edges = grp_edges[, .(from = source, to = target, coef_abs, coef_mean)],
                 nodes = data.frame(name = c("TBX21","RUNX3","MEF2C","KLF2","POU2F2","PRDM1")),
                 directed = TRUE)
  
  graph_list[[grp]] <- ggraph(g, layout = "circle") +
    geom_edge_fan(aes(width = coef_abs, color = coef_mean > 0),
                  arrow = arrow(length = unit(3, "mm"), type = "closed"),
                  end_cap = circle(4, "mm"), show.legend = TRUE) +
    scale_edge_width(range = c(0.5, 3), name = "|coef|") +
    scale_edge_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                            labels = c("TRUE" = "activating", "FALSE" = "repressive"),
                            name = "direction") +
    geom_node_point(size = 10, color = "grey85") +
    geom_node_text(aes(label = name), size = 3, fontface = "bold") +
    labs(title = grp) +
    theme_void()
}

network_combined <- wrap_plots(graph_list, ncol = 2) +
  plot_layout(guides = "collect")
ggsave("/home/rstudio/mnt_out/NK_project/figures/TF_TF_network_by_group.pdf", network_combined, width = 10, height = 9)
message("Wrote TF_TF_network_by_group.pdf")

# ---- bar chart of edge strengths, all pairs x all groups --------------------
edges[, pair := paste0(source, "\u2192", target)]
edges[, group := factor(group, levels = c("control","ATRXdel","ATRXwtMYCNwt","MYCNamp"))]

bar_plot <- ggplot(edges[found == TRUE], aes(x = pair, y = coef_mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  coord_flip() +
  labs(x = NULL, y = "coef_mean (signed edge strength)",
       title = "TF\u2192TF edges found in at least one group's fitted network") +
  theme_bw(base_size = 11)

ggsave("/home/rstudio/mnt_out/NK_project/figures/TF_TF_edge_strength_barplot.pdf", bar_plot, width = 8, height = 5)
message("Wrote TF_TF_edge_strength_barplot.pdf")