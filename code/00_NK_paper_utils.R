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

# Function to map gene names to Entrez IDs
map_gene_to_entrez <- function(gene_names) {
  # Query the org.Hs.eg.db package for the Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_names,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Return the result as a named vector
  return(entrez_ids)
}


# Plot for MiloR DA analysis

plotDAbeeswarm <- function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
           group.by, ")?")
    }
    if (is.numeric(da.res[, group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  }
  else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  if (!is.factor(da.res[, "group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, 
                                               levels = unique(group_by)))
  }
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }
  beeswarm_pos <- ggplot_build(
    da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
      arrange(group_by) %>%
      ggplot(aes(group_by, logFC)) +
      geom_quasirandom()
  )
  
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  n_groups <- unique(da.res$group_by) %>% length()
  
  da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
    mutate(pos_x = pos_x, pos_y = pos_y) %>%
    ggplot(aes(pos_x, pos_y, color = logFC_color)) +
    scale_color_gradient2(
      low = "darkred", mid = "lightgray", high = "#0B5394", midpoint = 0,
      na.value = "gray", name = "logFC (signif only)"
    ) +  # <- updated here
    xlab(group.by) +
    ylab("Log Fold Change") +
    scale_x_continuous(
      breaks = seq(1, n_groups),
      labels = setNames(levels(da.res$group_by), seq(1, n_groups))
    ) +
    geom_point() +
    coord_flip() +
    theme_bw(base_size = 22) +
    theme(strip.text.y = element_text(angle = 0))
  
}
