# =============================================================================
# t-SNE of the shared embedding (CTT | PAT), coloured by CellTypist cell types, T0-T4
# Author: Gaia Cocconcelli
# Same input as UMAP_A: the joint PCA (50 PCs) of T<t>_joint.rds, built on the
# gene HVGs after the CC.Difference regression, no batch correction. Only the
# drawing algorithm changes (t-SNE instead of UMAP), so the two maps can be put
# side by side. The joint .rds files are only read.
# Labels and colours as 03_celltypist_plots.R:
# T0-T2: stem threshold 0.1, organoid types grouped
#        (T<t>_CellTypist_Cluster_labels_stemthr0.1.csv, column "type")
# T3-T4: best match, every Kanton type with its own colour, Unassigned < 0.4
#        (T<t>_CellTypist_Cluster_labels_nothr_p0.4.csv, column "type")
# t-SNE: scater::runTSNE on all 50 PCs, perplexity 30 (scater default).
# In a t-SNE the local neighbours are kept, but the distances between far
# islands and the size of the islands are not meaningful.
# Output: results/tSNE/joint_T<t>/
#         T<t>_shared_tSNE_CellTypist_<tag>.png
#         T<t>_tSNE.rds, T<t>_tSNE.csv (t-SNE coordinates per cell)
# =============================================================================

suppressMessages({library(SingleCellExperiment); library(scater); library(ggplot2)})

source("config.R")

# palettes of 03_celltypist_plots.R
early_types <- c("stem cells 1", "stem cells 2", "stem cells 3", "neuroectodermal-like cells", "NSC")
pal_early <- c("stem cells 1" = "#1f77b4", "stem cells 2" = "#17becf", "stem cells 3" = "#9467bd",
               "neuroectodermal-like cells" = "#2ca02c", "NSC" = "#d62728")
kanton_stage <- read.csv(file.path(ct_path, "Kanton_Cluster_by_stage.csv"), row.names = 1, check.names = FALSE)
org_types <- setdiff(rownames(kanton_stage), early_types)
org_cols <- c("NSC/radial glia" = "#F3C300", "G2M/S NPCs" = "#654522",
              "G2M/S dorsal progenitors 1" = "#C51B7D", "mesenchymal-like cells" = "#8DB600",
              "midbrain/hindbrain" = "#E68FAC", "radial glia 2" = "#0000CD",
              "ventral progenitors and neurons 1" = "#4B0082", "ventral progenitors and neurons 4" = "#C2B280",
              "G2M/S ventral progenitors and IPs" = "#F99379", "choroid plexus" = "#FF00FF",
              "cortical neurons 1" = "#00FF00", "cortical neurons 2" = "#00FFCC",
              "G2M/S dorsal progenitors 2" = "#0066FF", "gliogenic/outer radial glia" = "#336633",
              "IPs and early cortical neurons" = "#CC66FF", "radial glia 1" = "#110022",
              "ventral progenitors and neurons 2" = "#993333", "ventral progenitors and neurons 3" = "#FF0066")
stopifnot(all(org_types %in% names(org_cols)))
pal_1 <- c(pal_early, "organoid-only types" = "#ff7f0e", "Unassigned" = "grey80") # T0-T2
pal_2 <- c(pal_early, org_cols[org_types], "Unassigned" = "grey80")              # T3-T4

for (tp in 0:4) {
  tag <- if (tp <= 2) "stemthr0.1" else "nothr_p0.4"
  pal <- if (tp <= 2) pal_1 else pal_2
  tp_out <- file.path(tsne_path, paste0("joint_T", tp))
  dir.create(tp_out, showWarnings = FALSE, recursive = TRUE)

  sce <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
  ct <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_", tag, ".csv")))
  # the export wrote one row per cell, in the order of the joint object
  stopifnot(nrow(ct) == ncol(sce), all(ct$cell_index == seq_len(ncol(sce))),
            all(ct$barcode == colnames(sce)), all(ct$type %in% names(pal)))
  lab <- factor(ct$type, levels = names(pal)[names(pal) %in% ct$type])
  print(table(lab, sample = sce$sample))

  set.seed(101010)
  sce <- runTSNE(sce, dimred = "PCA", perplexity = 30, name = "TSNE")
  tsne <- reducedDim(sce, "TSNE")
  ts <- data.frame(cell = colnames(sce), sample = sce$sample, type = ct$type,
                   tsne1 = tsne[, 1], tsne2 = tsne[, 2])
  saveRDS(ts, file.path(tp_out, paste0("T", tp, "_tSNE.rds")))
  write.csv(ts, file.path(tp_out, paste0("T", tp, "_tSNE.csv")), row.names = FALSE)

  df <- data.frame(t1 = tsne[, 1], t2 = tsne[, 2], sample = sce$sample, type = lab)
  set.seed(101010)
  df <- df[sample(nrow(df)), ] # random plotting order, so no group of cells is always on top

  p <- ggplot(df, aes(t1, t2, colour = type)) +
    geom_point(size = 0.4, alpha = 0.8) + facet_wrap(~ sample) +
    scale_colour_manual(values = pal, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "t-SNE 1", y = "t-SNE 2",
         title = paste0("T", tp, " - shared t-SNE, CellTypist", if (tp <= 2) ", stem threshold 0.1" else "")) +
    theme_classic()
  ggsave(file.path(tp_out, paste0("T", tp, "_shared_tSNE_CellTypist_", tag, ".png")),
         p, width = 14, height = 7, dpi = 300, bg = "white")
  message("=== T", tp, ": t-SNE saved ===")
  rm(sce); invisible(gc())
}
