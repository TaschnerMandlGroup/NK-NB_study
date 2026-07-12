# =============================================================================
# 12_TF_network_diagram.R
# Force-directed TF-target network, ATRXdel vs control, built directly from
# grn_edge_diff_ATRXdel_vs_other_groups.csv (already computed in
# 03_grn_perturb.py). Scaled-up version of 04_'s small fixed 6-node circular
# TF-TF diagram: this one includes every gene (TF or target) touched by the
# top differential edges in each cluster, with:
#   - node colour = mean delta_coef_abs of that gene's edges
#                   (red = stronger in ATRXdel, blue = stronger in control)
#   - node size   = degree (number of edges touching that gene)
#   - layout      = force-directed ("fr"), not the small fixed circle layout
#
# N_TOP_EDGES_PER_CLUSTER keeps the plot readable (the reference example had
# ~40 nodes / ~60 edges) -- raise/lower if the result looks too sparse/dense.
#
# Requires ggrepel installed (used internally by geom_node_text(repel=TRUE)).
#   igraph: Csardi & Nepusz (2006), InterJournal Complex Systems, 1695
#   ggraph: Pedersen, T.L. -- https://ggraph.data-imaging.dev/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(igraph); library(ggraph); library(ggplot2)
})

FIG_DIR <- "/home/rstudio/mnt_out/NK_project/figures/"
N_TOP_EDGES_PER_CLUSTER <- 60

edge_diff <- fread("/home/rstudio/GRN/grn_edge_diff_ATRXdel_vs_other_groups.csv")
edge_diff <- edge_diff[group_compared == "control"]  # ATRXdel vs control specifically

for (cl in unique(edge_diff$cluster)) {
  dt <- edge_diff[cluster == cl]
  dt <- dt[order(-abs(delta_coef_abs))][seq_len(min(N_TOP_EDGES_PER_CLUSTER, .N))]
  if (nrow(dt) == 0) next

  g <- graph_from_data_frame(dt[, .(source, target, delta_coef_abs)], directed = TRUE)

  V(g)$deg <- degree(g, mode = "all")
  # per-node mean delta_coef_abs across its incident edges (+ = net stronger
  # in ATRXdel, - = net stronger in control)
  V(g)$mean_delta <- sapply(V(g)$name, function(n) {
    mean(dt$delta_coef_abs[dt$source == n | dt$target == n])
  })

  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.4, colour = "grey60",
                    arrow = arrow(length = unit(2, "mm"), type = "closed"),
                    end_cap = circle(2, "mm"), start_cap = circle(2, "mm")) +
    geom_node_point(aes(size = deg, colour = mean_delta)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) +
    scale_colour_gradient2(low = "steelblue", mid = "grey90", high = "firebrick",
                           midpoint = 0, name = "ATRXdel − control\n(Δ coef_abs)") +
    scale_size(range = c(3, 12), name = "degree") +
    labs(title = paste0("TF-target network: ATRXdel vs control [", cl, "]")) +
    theme_void()

  out_path <- paste0(FIG_DIR, "TF_network_ATRXdel_vs_control_",
                      gsub("[^A-Za-z0-9_]", "-", cl), ".pdf")
  ggsave(out_path, p, width = 10, height = 9)
  message("Wrote ", out_path)
}

message("[12_TF_network_diagram] done")
