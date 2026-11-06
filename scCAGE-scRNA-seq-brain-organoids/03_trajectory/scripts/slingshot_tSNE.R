# =============================================================================
# Slingshot pseudotime on the t-SNE, T0-T2
# Author: Gaia Cocconcelli
# No new computation: pseudotime per cell and curves from slingshot_pseudotime.R
# (T<t>_pseudotime.csv, T<t>_curves.csv). The curve was fitted on the 50 PCs of the
# joint object and is only drawn on the t-SNE (embedCurves), as on UMAP_A.
# Cells outside the trajectory (Unassigned, organoid-only types, stem cells 3) in grey.
# Output: results/slingshot/joint_T<t>/T<t>_pseudotime_tSNE.png
# =============================================================================

suppressMessages(library(ggplot2))

source("config.R")

for (tp in c(0, 1, 2)) {
  tp_path <- file.path(sling_path, paste0("joint_T", tp))
  pt <- read.csv(file.path(tp_path, paste0("T", tp, "_pseudotime.csv")))
  curves <- read.csv(file.path(tp_path, paste0("T", tp, "_curves.csv")))
  cu <- curves[curves$basis == "tsne", ]

  set.seed(101010)
  df <- pt[sample(nrow(pt)), ]
  df <- df[order(!is.na(df$pseudotime)), ] # grey cells first, so they stay below

  p <- ggplot(df, aes(tsne1, tsne2, colour = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.8) + facet_wrap(~ sample) +
    geom_path(data = cu, aes(x, y, group = lineage), inherit.aes = FALSE, linewidth = 1.2) +
    scale_colour_viridis_c(name = "Pseudotime", na.value = "grey85") +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = paste0("T", tp, " - pseudotime")) +
    theme_classic()
  ggsave(file.path(tp_path, paste0("T", tp, "_pseudotime_tSNE.png")),
         p, width = 14, height = 7, dpi = 300, bg = "white")
  message("=== T", tp, ": pseudotime on the t-SNE saved ===")
}
