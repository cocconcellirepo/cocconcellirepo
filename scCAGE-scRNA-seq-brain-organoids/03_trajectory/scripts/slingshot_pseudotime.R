# =============================================================================
# Pseudotime with slingshot (principal curves), CTT and PAT together, one timepoint at a time
# Author: Gaia Cocconcelli
# Idea from La Manno et al. 2018 (Fig. 4a): principal curve through the cells,
# pseudotime = position of each cell along the curve. Code as in OSCA.advanced,
# "Trajectory Analysis", slingshot section:
# the curve is fitted on the 50 PCs of the joint object (same PCA as UMAP_A,
# gene HVGs, CC.Difference regressed) and then only DRAWN on the maps with
# embedCurves(): on UMAP_A and on the t-SNE of 02_celltypist/scripts/05_tsne.R.
# OSCA uses the PCs because UMAP and t-SNE distort the distances between groups of cells.
# Cells in the trajectory (CellTypist, stem threshold 0.1):
# stem cells 1 + stem cells 2 merged in one group "stem cells" (root, start.clus),
# neuroectodermal-like cells, NSC. Unassigned, organoid-only types and
# stem cells 3 (3-6 cells) are left out and drawn in grey.
# If slingshot finds more than one lineage, the pseudotime of a cell is the mean
# over its lineages (rowMeans, as in OSCA).
# T0-T2 only: at T3-T4 these early types are almost gone (see the monocle3 maps for T0-T4).
# Output: results/slingshot/joint_T<t>/
#   T<t>_pseudotime.csv (one row per cell: pseudotime, UMAP_A and t-SNE coordinates)
#   T<t>_curves.csv (curves on UMAP_A and t-SNE), T<t>_slingshot.rds (SlingshotDataSet)
#   T<t>_pseudotime_UMAP_A.png
# The joint .rds files are only read.
# =============================================================================

## ---- Step 0: setup ----
suppressMessages({library(SingleCellExperiment); library(slingshot); library(ggplot2)})

source("config.R")
thr <- 0.1 # stem threshold of the chosen labels

# CellTypist label -> group used in the trajectory (the other labels are left out)
traj_groups <- c("stem cells 1" = "stem cells", "stem cells 2" = "stem cells",
                 "neuroectodermal-like cells" = "neuroectodermal-like cells", "NSC" = "NSC")
root <- "stem cells"

## ---- Step 1: slingshot, one timepoint at a time ----
for (tp in c(0, 1, 2)) {
  tp_out <- file.path(sling_path, paste0("joint_T", tp))
  dir.create(tp_out, showWarnings = FALSE, recursive = TRUE)
  sce <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
  ct <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_stemthr", thr, ".csv")))
  # the export wrote one row per cell, in the order of the joint object
  stopifnot(nrow(ct) == ncol(sce), all(ct$cell_index == seq_len(ncol(sce))),
            all(ct$barcode == colnames(sce)))
  ts <- readRDS(file.path(tsne_path, paste0("joint_T", tp), paste0("T", tp, "_tSNE.rds")))
  stopifnot(identical(ts$cell, colnames(sce)))
  reducedDim(sce, "TSNE") <- as.matrix(ts[, c("tsne1", "tsne2")])
  sce$celltypist <- ct$type
  sce$traj <- unname(traj_groups[ct$type]) # NA = not in the trajectory
  print(table(group = sce$traj, sample = sce$sample, useNA = "ifany"))

  sub <- sce[, !is.na(sce$traj)]
  set.seed(101010)
  sub <- slingshot(sub, clusterLabels = sub$traj, reducedDim = "PCA", start.clus = root)
  print(slingLineages(SlingshotDataSet(sub))) # path(s) through the groups, from the root
  pseudo_paths <- slingPseudotime(sub)
  shared_pseudo <- rowMeans(pseudo_paths, na.rm = TRUE)

  # one row per cell of the joint object, NA outside the trajectory
  pt <- data.frame(cell = colnames(sce), sample = sce$sample, celltypist = sce$celltypist,
                   traj = sce$traj, pseudotime = NA_real_,
                   u1 = reducedDim(sce, "UMAP_A")[, 1], u2 = reducedDim(sce, "UMAP_A")[, 2],
                   tsne1 = ts$tsne1, tsne2 = ts$tsne2)
  pt$pseudotime[match(colnames(sub), pt$cell)] <- shared_pseudo
  print(aggregate(pseudotime ~ sample + traj, pt, median))
  write.csv(pt, file.path(tp_out, paste0("T", tp, "_pseudotime.csv")), row.names = FALSE)
  saveRDS(SlingshotDataSet(sub), file.path(tp_out, paste0("T", tp, "_slingshot.rds")))

  # curves drawn on UMAP_A and on the t-SNE (embedCurves, as in OSCA)
  curve_df <- function(dr) {
    emb <- slingCurves(embedCurves(sub, dr))
    do.call(rbind, lapply(seq_along(emb), function(i) {
      s <- emb[[i]]$s[emb[[i]]$ord, ]
      data.frame(lineage = paste0("Lineage", i), x = s[, 1], y = s[, 2])
    }))
  }
  curves <- rbind(cbind(basis = "umap", curve_df("UMAP_A")), cbind(basis = "tsne", curve_df("TSNE")))
  write.csv(curves, file.path(tp_out, paste0("T", tp, "_curves.csv")), row.names = FALSE)

  # pseudotime and curve on UMAP_A, CTT and PAT side by side (cells out of the trajectory in grey)
  set.seed(101010)
  df <- pt[sample(nrow(pt)), ]
  df <- df[order(!is.na(df$pseudotime)), ] # grey cells first, so they stay below
  cu <- curves[curves$basis == "umap", ]
  p <- ggplot(df, aes(u1, u2, colour = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.8) + facet_wrap(~ sample) +
    geom_path(data = cu, aes(x, y, group = lineage), inherit.aes = FALSE, linewidth = 1.2) +
    scale_colour_viridis_c(name = "Pseudotime", na.value = "grey85") +
    labs(x = "UMAP_A 1", y = "UMAP_A 2",
         title = paste0("T", tp, " - slingshot pseudotime (stem cells -> NE-like -> NSC)")) +
    theme_classic()
  ggsave(file.path(tp_out, paste0("T", tp, "_pseudotime_UMAP_A.png")),
         p, width = 14, height = 7, dpi = 300, bg = "white")
  message("=== T", tp, ": pseudotime saved ===")
  rm(sce, sub); invisible(gc())
}
