# =============================================================================
# NFAT-target vs. GSK3beta-pathway (Wnt/beta-catenin proxy) module score
# correlation across genetic groups, using Seurat AddModuleScore on the same
# scRNA object used to build the base GRN (script 01/03).
#
# IMPORTANT DIRECTION NOTE (read before interpreting results):
#   - NFAT targets: genes NFAT directly activates. Higher score = more active
#     nuclear NFAT transcriptional output.
#   - GSK3B acts primarily post-translationally, not transcriptionally on
#     itself, so there is no direct "GSK3B target gene" set to score. The
#     literature-standard transcriptional proxy is the canonical Wnt/beta-
#     catenin (TCF/LEF) target gene set, because GSK3B is the primary kinase
#     in the destruction complex that marks beta-catenin for degradation
#     (Jho et al., Mol Cell Biol 2002 for AXIN2 as a universal Wnt-target
#     reporter; He et al., Science 1998 for MYC; Tetsu & McCormick, Nature
#     1999 for CCND1). HIGHER GSK3B activity -> MORE beta-catenin degradation
#     -> LOWER Wnt-target module score. So the Wnt module is an INVERSE proxy;
#     we flip its sign below to get an intuitive "GSK3B activity" score.
#   - Prediction under the GSK3B->NFAT->TBX21 hypothesis (Beals et al.,
#     Genes Dev 1997a/b; PNAS 2002; Frontiers Oncol 2020): higher inferred
#     GSK3B activity should correlate with LOWER NFAT-target module score,
#     since GSK3B phosphorylates NFAT and drives its nuclear export/
#     inactivation.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(data.table)
})

# ---- gene sets ---------------------------------------------------------------
# NFAT target genes (Th1/cytotoxic-lymphocyte-relevant, literature-established):
#   IL2    - Rao lab, classic NFAT target (Rao et al., Annu Rev Immunol 1997)
#   IFNG   - NFAT-responsive in T/NK cells (Sica et al., J Biol Chem 1997)
#   GZMB   - NFAT-responsive in cytotoxic lymphocytes
#   FASLG  - NFAT/AP-1 composite element target (Latinis et al., MCB 1997)
#   TNF    - NFAT target (Tsai et al., MCB 1996)
#   CD69   - early NFAT-responsive activation marker (Castellanos et al., 1997)
#   PTGS2  - NFAT target (Iniguez et al., J Exp Med 2000)
#   NFATC1 - NFAT auto-amplification loop (Serfling et al., Biochim Biophys Acta 2006)
NFAT_TARGETS <- c("IL2", "IFNG", "GZMB", "FASLG", "TNF", "CD69", "PTGS2", "NFATC1")

# Canonical Wnt/beta-catenin (TCF/LEF) target genes - GSK3B-suppressible proxy
WNT_TARGETS <- c("AXIN2", "MYC", "CCND1", "LEF1", "TCF7")

# ---- assumes a Seurat object `nk` with metadata columns matching script 03's
#      GROUP/PATIENT conventions ("group", "rados_et_al_2025"). Adjust names
#      below if your object uses different column names. -----------------------
rds <- "/home/rstudio/mnt_out/NK_project/Rds/ALL_NK_cluster_withNewMetadata.Rds" 
nk <- readRDS(rds)

nk <- AddModuleScore(nk, features = list(NFAT_TARGETS), name = "NFAT_score")
nk <- AddModuleScore(nk, features = list(WNT_TARGETS),  name = "WNT_score")

# Seurat appends "1" to single-list AddModuleScore output columns
nk$NFAT_module <- nk$NFAT_score1
nk$WNT_module  <- nk$WNT_score1
nk$GSK3B_activity_proxy <- -1 * scale(nk$WNT_module)[, 1]  # flip sign: higher = more inferred GSK3B activity

# ---- pseudobulk to patient level, matching the pseudobulk-aggregation logic
#      used throughout the rest of this analysis (Squair et al., Nat Commun 2021)
meta <- nk@meta.data %>%
  as.data.table() %>%
  .[, .(NFAT_module = mean(NFAT_module, na.rm = TRUE),
        GSK3B_activity_proxy = mean(GSK3B_activity_proxy, na.rm = TRUE),
        n_cells = .N),
    by = .(rados_et_al_2025, group)]
meta[, group := relevel(factor(group), ref = "control")]

print(meta[order(group, rados_et_al_2025)])

# ---- correlation, overall and per group --------------------------------------
cat("\n--- Overall correlation (all patients pooled) ---\n")
print(cor.test(meta$GSK3B_activity_proxy, meta$NFAT_module, method = "pearson"))

cat("\n--- Per-group correlation ---\n")
for (grp in unique(meta$group)) {
  sub <- meta[group == grp]
  if (nrow(sub) >= 3) {
    cat("\n", as.character(grp), ":\n")
    print(cor.test(sub$GSK3B_activity_proxy, sub$NFAT_module, method = "pearson"))
  } else {
    cat("\n", as.character(grp), ": n too small for correlation (n =", nrow(sub), ")\n")
  }
}

# ---- does the relationship differ by genotype? (interaction test) -----------
fit <- lm(NFAT_module ~ GSK3B_activity_proxy * group, data = meta)
cat("\n--- lm(NFAT_module ~ GSK3B_activity_proxy * group) ---\n")
print(summary(fit))

# ---- scatter plot -------------------------------------------------------------
p <- ggplot(meta, aes(GSK3B_activity_proxy, NFAT_module, color = group, size = n_cells)) +
  geom_point(alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = 2, linewidth = 0.6) +
  labs(x = "Inferred GSK3\u03b2 activity (\u2212Wnt/\u03b2-catenin target module, patient mean)",
       y = "NFAT-target module score (patient mean)",
       title = "NFAT-target vs. GSK3\u03b2-activity proxy, by genotype",
       size = "cells/patient", color = "group") +
  theme_bw(base_size = 12)

ggsave("/home/rstudio/mnt_out/NK_project/figures/NFAT_vs_GSK3B_correlation.pdf", p, width = 6.5, height = 5)
message("Wrote NFAT_vs_GSK3B_correlation.pdf")



#### ADDITIVE MODEL ######

# ---- additive model: shared slope, group-specific intercept -----------------
fit_additive <- lm(NFAT_module ~ GSK3B_activity_proxy + group, data = meta)
cat("--- lm(NFAT_module ~ GSK3B_activity_proxy + group) ---\n")
print(summary(fit_additive))

# ---- group-level mean comparison: GSK3B activity proxy ----------------------
fit_gsk3b_by_group <- lm(GSK3B_activity_proxy ~ group, data = meta)
cat("\n--- lm(GSK3B_activity_proxy ~ group) ---\n")
print(summary(fit_gsk3b_by_group))
cat("\n--- emmeans (GSK3B activity proxy by group) ---\n")
print(emmeans::emmeans(fit_gsk3b_by_group, ~ group))

# ---- group-level mean comparison: NFAT module --------------------------------
fit_nfat_by_group <- lm(NFAT_module ~ group, data = meta)
cat("\n--- lm(NFAT_module ~ group) ---\n")
print(summary(fit_nfat_by_group))
cat("\n--- emmeans (NFAT module by group) ---\n")
print(emmeans::emmeans(fit_nfat_by_group, ~ group))