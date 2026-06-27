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

# ─── Publication labels (HTML/ggtext) ────────────────────────────────────────
# Named vector: internal key → display label with HTML superscripts.
# Used as x_labels / legend_labels in pub_violin().
PUB_GROUP_LABELS <- c(
  "control"      = "CTRL",
  "ATRXdel"      = "ATRX<sup>del</sup>",
  "MYCNamp"      = "MYCN<sup>amp</sup>",
  "ATRXwtMYCNwt" = "ATRX<sup>wt</sup>MYCN<sup>wt</sup>"
)

# ─── pub_violin() ─────────────────────────────────────────────────────────────
# Publication-ready violin plot with optional grouped/dodged layout, stats,
# internal boxplot, and individual data points.
#
# Parameters:
#   data          – data.frame containing x_col and y_col
#   x_col         – column for x-axis categories
#   y_col         – column for numeric y values
#   fill_col      – column for fill/colour grouping; NULL → fill = x_col
#   colors        – named colour vector (names = fill levels)
#   x_order       – display order for x_col levels
#   fill_order    – display order for fill_col levels (defaults to x_order)
#   x_labels      – named HTML character vector for x-axis tick labels
#   legend_labels – named HTML character vector for legend labels
#   y_label       – y-axis title (plain text)
#   y_label_color – colour of y-axis title
#   y_label_face  – font face of y-axis title ("italic", "bold", "bold.italic")
#   title         – optional plot title
#   ref_group     – Wilcoxon vs this reference group (ungrouped violins only);
#                   set NULL to suppress stats
#   comparisons   – explicit list of 2-element vectors for pairwise brackets;
#                   overrides ref_group
#   stat_label    – "p.signif" (stars) or "p.format" (numeric)
#   hide_ns       – hide non-significant brackets (default FALSE)
#   show_points   – overlay jittered individual data points
#   point_size    – size of data points
#   point_alpha   – alpha of data points
#   show_boxplot  – draw thin boxplot inside each violin
#   violin_scale  – "width" or "count"
#   dodge_width   – dodge width for grouped violins
#   base_size     – base font size in pt (10–12 recommended for publication)
#   font_family   – registered font family name (default "Arial")
#
# Returns: a ggplot object
pub_violin <- function(
  data,
  x_col         = "group",
  y_col         = "value",
  fill_col      = NULL,
  colors        = COLOR_CODE_GROUP,
  x_order       = GROUP_ORDER,
  fill_order    = NULL,
  x_labels      = PUB_GROUP_LABELS,
  legend_labels = PUB_GROUP_LABELS,
  y_label       = "",
  y_label_color = "black",
  y_label_face  = "italic",
  title         = NULL,
  ref_group     = "control",
  comparisons   = NULL,
  stat_label    = "p.signif",
  hide_ns       = FALSE,
  show_points   = FALSE,
  point_size    = 0.4,
  point_alpha   = 1,
  show_boxplot  = NULL,   # NULL = auto: TRUE when show_points=FALSE, FALSE otherwise
  violin_scale  = "width",
  dodge_width   = 0.85,
  base_size     = 10,
  font_family   = "Arial",
  x_angle       = 45
) {
  grouped <- !is.null(fill_col)
  if (!grouped) fill_col <- x_col

  # Factor ordering — only keep levels present in the data
  x_present    <- unique(as.character(data[[x_col]]))
  fill_present <- unique(as.character(data[[fill_col]]))
  x_levels    <- c(intersect(x_order, x_present),
                   setdiff(x_present, x_order))
  fo           <- if (!is.null(fill_order)) fill_order else x_order
  fill_levels  <- c(intersect(fo, fill_present),
                    setdiff(fill_present, fo))

  data[[x_col]]    <- factor(as.character(data[[x_col]]),    levels = x_levels)
  data[[fill_col]] <- factor(as.character(data[[fill_col]]), levels = fill_levels)
  data <- data[!is.na(data[[x_col]]) & !is.na(data[[y_col]]), ]

  dodge_pos <- if (grouped) ggplot2::position_dodge(dodge_width) else "identity"

  # Legend labels for fill scale
  leg_lbl <- if (!is.null(legend_labels)) {
    lv <- legend_labels[fill_levels]
    lv[is.na(lv)] <- fill_levels[is.na(lv)]
    lv
  } else {
    ggplot2::waiver()
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x    = .data[[x_col]],
                 y    = .data[[y_col]],
                 fill = .data[[fill_col]])
  ) +
    ggplot2::geom_violin(
      scale     = violin_scale,
      trim      = TRUE,
      linewidth = 0,
      color     = NA,
      position  = dodge_pos,
      alpha     = 1
    ) +
    ggplot2::scale_fill_manual(values = colors, labels = leg_lbl, name = NULL) +
    ggplot2::labs(x = NULL, y = y_label, title = title) +
    ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      axis.title.y      = ggplot2::element_text(
                            face  = y_label_face,
                            color = y_label_color,
                            size  = base_size),
      axis.text.x       = ggtext::element_markdown(size  = base_size - 1,
                                                    color = "black",
                                                    angle = x_angle,
                                                    hjust = if (x_angle != 0) 1 else 0.5,
                                                    vjust = if (x_angle != 0) 1 else 0.5),
      axis.text.y       = ggplot2::element_text(size = base_size - 1,
                                                color = "black"),
      axis.line         = ggplot2::element_line(linewidth = 0.4),
      axis.ticks        = ggplot2::element_line(linewidth = 0.4),
      axis.title.x      = ggplot2::element_blank(),
      legend.text       = ggtext::element_markdown(size = base_size - 1),
      legend.key.size   = ggplot2::unit(3.5, "mm"),
      legend.key        = ggplot2::element_rect(colour = NA, fill = NA),
      legend.background = ggplot2::element_blank(),
      legend.title      = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(size = base_size, face = "bold",
                                                hjust = 0),
      plot.margin       = ggplot2::margin(t = 14, r = 5, b = 3, l = 5),
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(size = base_size, face = "italic")
    )

  # Boxplot: auto-show when no points are drawn; always black lines
  draw_box <- if (is.null(show_boxplot)) FALSE else show_boxplot
  if (draw_box) {
    p <- p + ggplot2::geom_boxplot(
      width         = if (grouped) 0.08 else 0.14,
      outlier.shape = NA,
      position      = dodge_pos,
      fill          = NA,
      color         = "black",
      linewidth     = 0.35
    )
  }

  # Optional individual data points
  if (show_points) {
    p <- p + ggplot2::geom_point(
      position = if (grouped)
        ggplot2::position_jitterdodge(jitter.width = 0.1,
                                      dodge.width  = dodge_width,
                                      seed         = 42)
      else
        ggplot2::position_jitter(width = 0.08, seed = 42),
      size  = point_size,
      alpha = point_alpha,
      color = "black",
      shape = 16
    )
  }

  # X-axis tick labels (HTML markdown)
  if (!is.null(x_labels)) {
    lbl <- x_labels[x_levels]
    lbl[is.na(lbl)] <- x_levels[is.na(lbl)]
    p <- p + ggplot2::scale_x_discrete(labels = lbl)
  }

  # Statistical annotations (ungrouped violins only — for grouped, use
  # rstatix::wilcox_test() + ggpubr::stat_pvalue_manual() externally)
  if (!grouped) {
    comps_to_use <- if (!is.null(comparisons)) {
      comparisons
    } else if (!is.null(ref_group)) {
      # build explicit ref vs each other level → bracket-line style
      lapply(setdiff(x_levels, ref_group), function(g) c(ref_group, g))
    } else {
      NULL
    }
    if (!is.null(comps_to_use)) {
      p <- p + ggpubr::stat_compare_means(
        comparisons  = comps_to_use,
        method       = "wilcox.test",
        label        = stat_label,
        size         = 2.8,
        bracket.size = 0.3,
        tip.length   = 0.01,
        hide.ns      = hide_ns,
        symnum.args  = list(
          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
          symbols   = c("****", "***", "**", "*", "ns")
        )
      )
    }
  }

  p
}

# ─── Publication theme ────────────────────────────────────────────────────────
# Single authoritative pub_theme used by all manuscript figure scripts.
# Source this file and call theme_set(pub_theme) at the top of each script.
pub_theme <- theme_classic(base_size = 9) +
  theme(
    axis.line        = element_line(linewidth = 0.4, color = "black"),
    axis.ticks       = element_line(linewidth = 0.4, color = "black"),
    axis.text        = element_text(size = 8, color = "black"),
    axis.title       = element_text(size = 9, color = "black"),
    plot.title       = element_text(size = 10, face = "bold"),
    plot.subtitle    = element_text(size = 8),
    legend.key.size  = unit(3, "mm"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9),
    strip.background = element_blank(),
    strip.text       = element_text(size = 9)
  )

# ─── Enrichment helpers ───────────────────────────────────────────────────────
# run_enrichment / show_enrichment / run_treeplot are shared by scripts 36 and 37.
# Requires: clusterProfiler, enrichplot, org.Hs.eg.db, msigdbr loaded upstream.

run_enrichment <- function(gene_symbols,
                           type = c("BP", "MF", "CC", "KEGG", "HALLMARK"),
                           term2gene = NULL,
                           label = "") {
  type <- match.arg(type)
  if (length(gene_symbols) == 0) {
    message(label, ": no input genes - skipping"); return(invisible(NULL))
  }
  entrez_ids <- na.omit(map_gene_to_entrez(gene_symbols))
  if (length(entrez_ids) == 0) {
    message(label, ": no valid Entrez IDs - skipping"); return(invisible(NULL))
  }
  tryCatch(
    {
      if (type %in% c("BP", "MF", "CC")) {
        enrichGO(gene = entrez_ids, ont = type, OrgDb = "org.Hs.eg.db",
                 readable = TRUE, pAdjustMethod = "bonferroni",
                 pvalueCutoff = 0.05)
      } else if (type == "KEGG") {
        enrichKEGG(gene = entrez_ids, organism = "hsa", keyType = "kegg",
                   pAdjustMethod = "bonferroni", pvalueCutoff = 0.05)
      } else if (type == "HALLMARK") {
        enricher(gene = entrez_ids, TERM2GENE = term2gene,
                 pAdjustMethod = "bonferroni", pvalueCutoff = 0.05)
      }
    },
    error = function(e) {
      message(label, " (", type, ") error: ", conditionMessage(e)); NULL
    }
  )
}

# Show enrichment dotplot inline; save PDF when out_path provided
show_enrichment <- function(result, label = "", out_path = NULL,
                            out_w = 5, out_h = 4.5) {
  if (is.null(result) || nrow(result@result) == 0) {
    print(paste0("0 enriched terms - ", label)); return(invisible(NULL))
  }
  signif_rows <- result@result[result@result$p.adjust < 0.05, ]
  if (nrow(signif_rows) == 0) {
    print(paste0("0 terms with bonferroni p.adjust < 0.05 - ", label))
    return(invisible(NULL))
  }
  p <- dotplot(result, showCategory = 20) +
    ggtitle(label) + pub_theme +
    theme(axis.text.y = element_text(size = 7))
  print(p)
  if (!is.null(out_path)) {
    pdf(out_path, height = out_h, width = out_w, useDingbats = FALSE)
    print(p); dev.off()
  }
}

run_treeplot <- function(result, label = "", n_cluster = 5, n_words = 1,
                         show_cat = 15, out_path = NULL, out_w = 9, out_h = 4) {
  if (is.null(result) || nrow(result@result) == 0) {
    message("treeplot skipped (no enriched terms): ", label); return(invisible(NULL))
  }
  signif_rows <- result@result[result@result$p.adjust < 0.05, ]
  if (nrow(signif_rows) < 2) {
    message("treeplot skipped (< 2 significant terms): ", label); return(invisible(NULL))
  }
  p <- tryCatch(
    {
      enrichres2 <- pairwise_termsim(result)
      treeplot(enrichres2, showCategory = show_cat, color = "p.adjust",
               offset = 10, offset_tiplab = 0.5,
               nCluster = min(n_cluster, nrow(signif_rows)),
               nWords = n_words, split = "subcategory")
    },
    error = function(e) { message("treeplot error: ", conditionMessage(e)); NULL }
  )
  if (!is.null(p)) {
    print(p)
    if (!is.null(out_path)) {
      pdf(out_path, height = out_h, width = out_w, useDingbats = FALSE)
      print(p); dev.off()
    }
  }
}
