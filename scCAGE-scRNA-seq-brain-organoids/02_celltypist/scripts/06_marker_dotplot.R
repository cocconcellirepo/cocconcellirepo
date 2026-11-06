# =============================================================================
# Dot plot of the marker genes of each CellTypist identity, CTT and PAT side by side
# Author: Gaia Cocconcelli
# Same idea as Figure 3A of Swamy et al. 2021 (GigaScience, scEiaD): the 5 most specific
# genes of each cell type, dot size = % of cells where the gene is detected, colour =
# expression. Here one plot per timepoint, with one column per identity x genotype.
# Marker genes: scran::scoreMarkers (OSCA, chapter "Marker gene detection"), pairwise
# comparisons between identities, blocking on the sample so that CTT and PAT count the
# same. Genes ranked by mean.AUC = mean of the AUCs of the identity against each of the
# other identities, so a gene is on top only if it is higher than in most of the others
# (specificity). I first used findMarkers + summary.AUC, but with the default
# pval.type = "any" summary.AUC is the AUC of only the comparison with the lowest p-value
# (scran docs), so generic genes high against just one identity went on top.
# Top 5 per identity, a gene is shown only once (in the first identity where it is in
# the top 5). logcounts of the joint object (multiBatchNorm, not regressed).
# Identities with at least 20 cells in the timepoint; Unassigned is not an identity
# and is left out. In the plot an identity x genotype column is shown only if it has at
# least 20 cells, and the genes of an identity without any column are not shown.
# Labels: T0-T2 stem threshold 0.1, T3-T4 best match (03_celltypist_plots.R).
# The joint .rds files are only read.
# Run from the project folder, with the timepoints as arguments:
#   Rscript 02_celltypist/scripts/06_marker_dotplot.R 0 1 2 3 4
# Output (results/marker_dotplot/): T<t>_marker_dotplot.png, T<t>_markers_top50.csv (top 50 per identity)
# =============================================================================

suppressMessages({library(SingleCellExperiment); library(scran); library(scuttle); library(ggplot2)})

source("config.R")
dir.create(marker_path, showWarnings = FALSE, recursive = TRUE)
min_cells <- 20
n_top <- 5

# colours and order of the identities as in 03_celltypist_plots.R
pal <- c("stem cells 1" = "#1f77b4", "stem cells 2" = "#17becf", "stem cells 3" = "#9467bd",
         "neuroectodermal-like cells" = "#2ca02c", "NSC" = "#d62728", "organoid-only types" = "#ff7f0e",
         "NSC/radial glia" = "#F3C300", "G2M/S NPCs" = "#654522", "G2M/S dorsal progenitors 1" = "#C51B7D",
         "mesenchymal-like cells" = "#8DB600", "midbrain/hindbrain" = "#E68FAC", "radial glia 2" = "#0000CD",
         "ventral progenitors and neurons 1" = "#4B0082", "ventral progenitors and neurons 4" = "#C2B280",
         "G2M/S ventral progenitors and IPs" = "#F99379", "choroid plexus" = "#FF00FF",
         "cortical neurons 1" = "#00FF00", "cortical neurons 2" = "#00FFCC", "G2M/S dorsal progenitors 2" = "#0066FF",
         "gliogenic/outer radial glia" = "#336633", "IPs and early cortical neurons" = "#CC66FF",
         "radial glia 1" = "#110022", "ventral progenitors and neurons 2" = "#993333",
         "ventral progenitors and neurons 3" = "#FF0066")

for (tp in as.integer(commandArgs(trailingOnly = TRUE))) {
  tag <- if (tp <= 2) "stemthr0.1" else "nothr_p0.4"
  sce <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
  ct <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_", tag, ".csv")))
  stopifnot(nrow(ct) == ncol(sce), all(ct$barcode == colnames(sce)), all(ct$type %in% c(names(pal), "Unassigned")))
  sce$celltype <- ct$type
  sce$geno <- substr(as.character(sce$sample), 1, 3)

  n_type <- table(sce$celltype)
  types <- names(pal)[names(pal) %in% names(n_type)[n_type >= min_cells]]
  sub <- sce[, sce$celltype %in% types]
  rm(sce); invisible(gc())

  mk <- scoreMarkers(sub, groups = sub$celltype, block = sub$sample, assay.type = "logcounts")
  top <- do.call(rbind, lapply(types, function(ty) {
    m <- mk[[ty]]
    m <- m[order(m$mean.AUC, decreasing = TRUE), ]
    k <- seq_len(min(50, nrow(m)))
    data.frame(timepoint = paste0("T", tp), celltype = ty, gene = rownames(m)[k],
               mean.AUC = round(m$mean.AUC[k], 3), min.AUC = round(m$min.AUC[k], 3),
               mean.logFC.cohen = round(m$mean.logFC.cohen[k], 2))
  }))
  write.csv(top, file.path(marker_path, paste0("T", tp, "_markers_top50.csv")), row.names = FALSE)

  # columns: identity x genotype with at least min_cells cells
  grp <- paste(sub$celltype, sub$geno, sep = " | ")
  n_grp <- table(grp)
  grp_levels <- unlist(lapply(types, function(ty) paste(ty, c("CTT", "PAT"), sep = " | ")))
  grp_levels <- grp_levels[grp_levels %in% names(n_grp)[n_grp >= min_cells]]
  plot_types <- types[types %in% sub(" \\| (CTT|PAT)$", "", grp_levels)]

  # top 5 per identity, each gene only once
  used <- character(0); sel <- list()
  for (ty in plot_types) {
    g <- setdiff(top$gene[top$celltype == ty], used)
    g <- head(g, n_top)
    used <- c(used, g)
    sel[[ty]] <- data.frame(celltype = rep(ty, length(g)), gene = g)
  }
  sel <- do.call(rbind, sel)

  keep <- grp %in% grp_levels
  st <- summarizeAssayByGroup(sub[sel$gene, keep], ids = grp[keep],
                              statistics = c("mean", "prop.detected"), assay.type = "logcounts")
  df <- data.frame(gene = rep(rownames(st), ncol(st)), group = rep(colnames(st), each = nrow(st)),
                   mean = as.vector(assay(st, "mean")), detected = 100 * as.vector(assay(st, "prop.detected")))
  df$gene <- factor(df$gene, levels = rev(sel$gene))
  df$group <- factor(df$group, levels = grp_levels)
  gene_cols <- rev(unname(pal[sel$celltype]))

  p <- ggplot(df, aes(group, gene)) +
    geom_point(aes(size = detected, colour = mean)) +
    scale_size(range = c(0, 4), limits = c(0, 100), name = "% detected") +
    scale_colour_viridis_c(option = "magma", begin = 0.1, end = 0.9, name = "Mean logcounts") +
    labs(x = NULL, y = NULL, title = paste0("T", tp, " markers per identity")) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = suppressWarnings(element_text(colour = gene_cols)),
          plot.margin = margin(5.5, 5.5, 5.5, 60))  # room for the first rotated column label
  ggsave(file.path(marker_path, paste0("T", tp, "_marker_dotplot.png")), p,
         width = 4 + 0.3 * length(grp_levels), height = 2 + 0.14 * nrow(sel), dpi = 200,
         bg = "white", limitsize = FALSE)
  message("T", tp, ": ", length(types), " identities (", length(plot_types), " in the plot), ",
          nrow(sel), " genes, ", length(grp_levels), " columns, saved")
  rm(sub, mk); invisible(gc())
}
