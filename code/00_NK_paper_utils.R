#!/usr/bin/env Rscript 

# from: https://www.bioconductor.org/packages/release/bioc/vignettes/aggregateBioVar/inst/doc/multi-subject-scRNA-seq.html
scPHeatmap <- function(sumExp, subjectVar, gtVar, logSample = 1:100, ...) {
  orderSumExp <- sumExp[, order(sumExp[[subjectVar]])]
  sumExpCounts <- as.matrix(
    SummarizedExperiment::assay(orderSumExp, "counts"))
  logcpm <- log2(
    1e6*t(t(sumExpCounts) / colSums(sumExpCounts)) + 1)
  annotations <- data.frame(
    Sex = orderSumExp[[gtVar]],
    Donor = orderSumExp[[subjectVar]])
  rownames(annotations) <- colnames(orderSumExp)
  
  singleCellpHeatmap <- pheatmap::pheatmap(
    mat = logcpm[logSample, ], annotation_col = annotations,
    cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
    scale = "none", ...)
  return(singleCellpHeatmap)
}


# Function to add theme for ggplot
deseq_themes <- function() {
  list(
    theme_classic(),
    lims(x = c(-4, 5), y = c(0, 80)),
    labs(
      x = "log<sub>2</sub> (fold change)",
      y = "-log<sub>10</sub> (p<sub>adj</sub>)"
    ),
    ggplot2::theme(
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown())
  )
}
