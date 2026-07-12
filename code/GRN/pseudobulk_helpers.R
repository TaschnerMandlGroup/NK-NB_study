# =============================================================================
# pseudobulk_helpers.R
# Shared pseudobulk -> lm -> emmeans -> pairs -> ggplot-panel helpers, used by
# 04_pseudobulk_stats.R (TBX21/RUNX3 + top-6-TF mega figure) and
# 09_additional_TF_cluster_stats.R (cluster-stratified MEF2C/KLF2/PRDM1/FOS/
# POU2F2). Extracted so the two scripts don't carry duplicate copies of the
# same three functions.
#
# Why patient-level: cells from the same donor are not independent. Testing at
# the cell level (pseudoreplication) inflates the FDR and can manufacture
# "significant" differences where none exist -- pseudobulk aggregation to the
# biological replicate is the accepted fix.
#   Squair et al., Nat Commun 2021 (PMID 34584091)
#   Zimmerman et al., Nat Commun 2021; sc-best-practices (DGE chapter)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(emmeans); library(ggplot2)
})

# ---- shared: Dunnett-style "each group vs control" p-values, positioned
#      above each group's box. Uses trt.vs.ctrl (not all-pairs Tukey) since
#      every comparison here has one fixed reference (control) -- this is
#      the correct multiplicity correction for that design, and it's a
#      single label per non-reference group instead of every pairwise combo.
#      Returns NULL (no labels added) if fewer than 2 groups are present.
add_pvalue_labels <- function(p, pb, fit, value_col) {
  if (uniqueN(pb$group) < 2) return(p)

  dunnett <- as.data.frame(contrast(emmeans(fit, ~ group), method = "trt.vs.ctrl", ref = 1))
  dunnett$group <- sub(" - .*$", "", dunnett$contrast)
  dunnett$label <- ifelse(dunnett$p.value < 0.001, "p<0.001",
                           paste0("p=", formatC(dunnett$p.value, digits = 3, format = "f")))

  ymax <- pb[, .(y = max(get(value_col), na.rm = TRUE)), by = group]
  ymax[, y := y + 0.08 * diff(range(pb[[value_col]], na.rm = TRUE))]
  lbl <- merge(dunnett, ymax, by = "group")

  p + geom_text(data = lbl, aes(x = group, y = y, label = label),
                inherit.aes = FALSE, size = 2.8)
}

# ---- pseudobulk -> lm -> emmeans -> pairs -> ggplot panel, pooled across
#      NK subtype clusters (feeds the combined figures in 04_) --------------
run_perturb_analysis <- function(csv_path, label, id_cols) {
  d <- fread(csv_path)
  cyto_cols <- setdiff(colnames(d), id_cols)
  cyto_cols <- cyto_cols[!cyto_cols %in% c("group", "cluster", "patient")]
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
         y = "Δ cytotoxicity module (patient mean)",
         size = "cells/patient") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p <- add_pvalue_labels(p, pb, fit, "cyto_delta")

  list(pb = pb, fit = fit, plot = p)
}

# ---- same as above, but one model PER NK SUBTYPE CLUSTER -------------------
# Returns a named list keyed by cluster label.
run_perturb_analysis_by_cluster <- function(csv_path, label, id_cols) {
  d <- fread(csv_path)
  cyto_cols <- setdiff(colnames(d), id_cols)
  cyto_cols <- cyto_cols[!cyto_cols %in% c("group", "cluster", "patient")]
  d[, cyto_delta := rowMeans(.SD, na.rm = TRUE), .SDcols = cyto_cols]

  out <- list()
  for (cl in sort(unique(d$cluster))) {
    dt <- d[cluster == cl]
    pb <- dt[, .(cyto_delta = mean(cyto_delta, na.rm = TRUE),
                 n_cells    = .N),
             by = .(patient, group)]
    if (uniqueN(pb$group) < 2 || nrow(pb) < 3) {
      message(label, " [", cl, "]: too few group/patient combinations to fit -- skipped")
      next
    }
    pb[, group := relevel(factor(group), ref = "control")]

    fit <- lm(cyto_delta ~ group, data = pb)
    sub_label <- paste0(label, " [", cl, "]")
    cat("\n=== ", sub_label, " ===\n")
    cat("--- lm(cyto_delta ~ group) ---\n"); print(summary(fit))
    cat("--- pairwise genotype contrasts vs control (Tukey) ---\n"); print(pairs(emmeans(fit, ~ group)))

    p <- ggplot(pb, aes(group, cyto_delta)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
      geom_boxplot(outlier.shape = NA, width = 0.6) +
      geom_jitter(aes(size = n_cells), width = 0.12, alpha = 0.7) +
      labs(title = sub_label, x = NULL,
           y = "Δ cytotoxicity module (patient mean)",
           size = "cells/patient") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    p <- add_pvalue_labels(p, pb, fit, "cyto_delta")

    out[[cl]] <- list(pb = pb, fit = fit, plot = p)
  }
  out
}

# ---- single-gene readout-chain version of run_perturb_analysis (pooled
#      across clusters) -- e.g. RUNX3 KO -> effect on TBX21 expression -------
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
         y = paste0("Δ ", gene_col, " expression (patient mean)"),
         size = "cells/patient") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p <- add_pvalue_labels(p, pb, fit, "gene_delta")

  list(pb = pb, fit = fit, plot = p)
}
