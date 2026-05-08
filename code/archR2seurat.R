#!/usr/bin/env Rscript

library(ArchR)
library(Seurat)

# Convert ArchR object to Seurat Object
# For now, converts the Peak Matrix and the GeneExpressionMatrix
# Run like this: ./archR2seurat.R <proj_filt>

addArchRThreads(threads = 6)
addArchRGenome("hg38")

### Load ArchR Object ###

args = commandArgs(trailingOnly=TRUE)
proj_filt_path = args[1]

proj_filt <- loadArchRProject(proj_filt_path)

## Matrices ##

# get peaks matrix

raw_peak_mat <- getMatrixFromProject(
  ArchRProj = proj_filt,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

pmat <- raw_peak_mat@assays@data$PeakMatrix
colnames(pmat) <- colnames(raw_peak_mat)

chr <- raw_peak_mat@rowRanges@seqnames
start <- data.frame(raw_peak_mat@rowRanges@ranges)$start
end <- data.frame(raw_peak_mat@rowRanges@ranges)$end

rownames(pmat) <- paste(chr, paste(start, end, sep = "-"), sep = ":")

# get RNA matrix

raw_gex_mat <- getMatrixFromProject(
  ArchRProj = proj_filt,
  useMatrix = "GeneExpressionMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

gmat <- raw_gex_mat@assays@data$GeneExpressionMatrix
colnames(gmat) <- colnames(raw_gex_mat)
rownames(gmat) <- raw_gex_mat@elementMetadata$name

## Metadata ##

metadata <- getCellColData(
  ArchRProj = proj_filt
 )

metadata$Barcodes <- rownames(metadata)

sorted_metadata <- metadata[match(colnames(gmat), metadata$Barcodes),]

### Create a Seurat Object ###

RNA_assay <- CreateAssayObject(counts = gmat)

seurat_obj <- CreateSeuratObject(
        counts = RNA_assay,
        assay = "RNA"
)

Peaks_assay <- CreateAssayObject(counts = pmat)

seurat_obj[["Peaks"]] <- Peaks_assay

seurat_obj <- AddMetaData(object = seurat_obj, metadata = sorted_metadata)

# add UMAP coords

umap_coords <- getEmbedding(ArchRProj = proj_filt, embedding = "UMAP", returnDF = TRUE)
umap_coords <- umap_coords[match(colnames(gmat), rownames(umap_coords)),]

mat_umap_coords <- as.matrix(umap_coords)
umap <- CreateDimReducObject(embeddings = mat_umap_coords)
seurat_obj[["UMAP"]] <- umap

# save seurat object
saveRDS(seurat_obj, file = "seurat_obj.rds")


