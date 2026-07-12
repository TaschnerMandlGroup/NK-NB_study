# =============================================================================
# 01_export_inputs.R
# Export CellOracle inputs from your R objects:
#   (a) scRNA as .h5ad  (raw counts + group + patient in obs, + a UMAP)
#   (b) all_peaks.csv           (chr_start_end, one per line)
#   (c) cicero_connections.csv  (Peak1, Peak2, coaccess)
#
# You need (b)+(c) only to build the BASE GRN (candidate TF->target edges from
# accessible cis-elements). The scRNA drives the per-group GRN fitting.
# ATAC and RNA need NOT be from the same cells.
#
# Refs:
#   CellOracle base-GRN construction .... Kamimoto et al., Nature 2023 (PMID 36755098)
#   Co-accessibility (Cicero) ........... Pliner et al., Mol Cell 2018 (PMID 30078726)
# =============================================================================

rna_rds     <- "/home/rstudio/mnt_out/NK_project/Rds/ALL_NK_cluster_withNewMetadata.Rds"     # scRNA Seurat object, NK cells only
group_col   <- "group"                 # obs col: ATRXmut / ATRXwtMYCNwt / control / MYCNamp
patient_col <- "rados_et_al_2025"                # obs col: biological replicate (3-5 per group)
genome      <- "hg38"                  # CD6 coords in your tracks are hg38
out_dir     <- "/home/rstudio/mnt_out/NK_project/h5ad"

archr_proj_path <- "/home/rstudio/mnt_out/NK_project/ArchR_v20260203/" # loadArchRProject() path
coaccess_cutoff <- 0.5                 # ArchR getCoAccessibility corCutOff

dir.create(out_dir, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(Seurat); library(sceasy); library(reticulate)
})

# scRNA -> h5ad
rna <- readRDS(rna_rds)
stopifnot(all(c(group_col, patient_col) %in% colnames(rna@meta.data)))
DefaultAssay(rna) <- "RNA"

# Keep it lean but preserve raw counts + the two grouping columns + embeddings.
# NOTE: patients within a group MUST be integrated/batch-corrected upstream, so
# the KNN imputation CellOracle runs later reflects state, not donor batch.
rna_export <- DietSeurat(rna, assays = "RNA",
                         dimreducs = intersect(c("UMAP","PCA"), names(rna@reductions)))
# make sure obs carries clean factor labels
rna_export$genetic_group <- as.character(rna_export[[group_col]][,1])
rna_export$patient_id    <- as.character(rna_export[[patient_col]][,1])

sceasy::convertFormat(
  rna_export, from = "seurat", to = "anndata",
  outFile = file.path(out_dir, "nk_rna.h5ad"),
  main_layer = "counts",           # CellOracle imports RAW counts
  drop_single_values = FALSE
)
message("Wrote ", file.path(out_dir, "nk_rna.h5ad"))

# peaks + co-accessibility
if (atac_source == "archr") {
  suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges) })
  proj <- loadArchRProject(archr_proj_path)
  proj <- addCoAccessibility(proj, reducedDims = "IterativeLSI")
  ps  <- getPeakSet(proj)
  pstr <- paste(as.character(seqnames(ps)), start(ps), end(ps), sep = "_")
  cA  <- getCoAccessibility(proj, corCutOff = coaccess_cutoff, returnLoops = FALSE)
  conns <- data.frame(Peak1 = pstr[cA$queryHits],
                      Peak2 = pstr[cA$subjectHits],
                      coaccess = cA$correlation)
  all_peaks <- data.frame(peak = pstr)
} else stop("atac_source must be 'signac' or 'archr'")

write.csv(all_peaks, file.path(out_dir, "all_peaks.csv"), row.names = FALSE)
write.csv(conns,     file.path(out_dir, "cicero_connections.csv"), row.names = FALSE)
message("Wrote all_peaks.csv (", nrow(all_peaks), ") and cicero_connections.csv (", nrow(conns), ")")

##################################
# PART II: NK2+NK3 clusters only #
##################################

rna_rds     <- "/home/rstudio/mnt_out/NK_project/Rds/ALL_NK_NK2+NK2.Rds"  

idxSample <- BiocGenerics::which(proj$predicted_NK_group %in% "immuno")
immunoNK_cells <- proj$cellNames[idxSample]

proj <- subsetArchRProject(
  ArchRProj = proj,
  cells = immunoNK_cells,
  outputDirectory = "ArchRSubset_NK2NK3",
  dropCells = TRUE,
  logFile = NULL,
  threads = 20,
  force = FALSE
)

proj <- addCoAccessibility(proj, reducedDims = "IterativeLSI")
ps  <- getPeakSet(proj)
pstr <- paste(as.character(seqnames(ps)), start(ps), end(ps), sep = "_")
cA  <- getCoAccessibility(proj, corCutOff = coaccess_cutoff, returnLoops = FALSE)
conns <- data.frame(Peak1 = pstr[cA$queryHits],
                    Peak2 = pstr[cA$subjectHits],
                    coaccess = cA$correlation)
all_peaks <- data.frame(peak = pstr)

saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRSubset_NK2NK3", load = FALSE)

write.csv(all_peaks, file.path(out_dir, "all_peaks_NK2NK3.csv"), row.names = FALSE)
write.csv(conns,     file.path(out_dir, "cicero_connections_NK2NK3.csv"), row.names = FALSE)
message("Wrote all_peaks_NK2NK3.csv (", nrow(all_peaks), ") and cicero_connections_NK2NK3.csv (", nrow(conns), ")")


# scRNA -> h5ad
rna <- readRDS(rna_rds)
stopifnot(all(c(group_col, patient_col) %in% colnames(rna@meta.data)))
DefaultAssay(rna) <- "RNA"

# Keep it lean but preserve raw counts + the two grouping columns + embeddings.
# NOTE: patients within a group MUST be integrated/batch-corrected upstream, so
# the KNN imputation CellOracle runs later reflects state, not donor batch.
rna_export <- DietSeurat(rna, assays = "RNA",
                         dimreducs = intersect(c("UMAP","PCA"), names(rna@reductions)))
# make sure obs carries clean factor labels
rna_export$genetic_group <- as.character(rna_export[[group_col]][,1])
rna_export$patient_id    <- as.character(rna_export[[patient_col]][,1])

sceasy::convertFormat(
  rna_export, from = "seurat", to = "anndata",
  outFile = file.path(out_dir, "nk_NK2NK3_rna.h5ad"),
  main_layer = "counts",           # CellOracle imports RAW counts
  drop_single_values = FALSE
)
message("Wrote ", file.path(out_dir, "nk_NK2NK3_rna.h5ad"))
