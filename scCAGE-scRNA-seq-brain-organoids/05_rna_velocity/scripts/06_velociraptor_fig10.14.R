# =============================================================================
# OSCA.advanced Figure 10.14 with velociraptor, on our scVelo models, T0-T4
# Author: Gaia Cocconcelli
# Steps as in OSCA.advanced (chapter "RNA velocity") and in the velociraptor vignette:
#   embedVelocity(reducedDim(sce, "TSNE"), velo.out)
#   gridVectors(sce, embedded, use.dimred = "TSNE", resolution = 30)
#   plotTSNE(sce, colour_by = "pseudotime", point_alpha = 0.3) + geom_segment(...)
# Differences from the tutorial:
#   - velo.out is not computed with velociraptor::scvelo(): it is our scVelo model
#     (results/velocity/scvelo/<sample>_scvelo.h5ad), read with zellkonverter::readH5AD.
#   - the t-SNE is the t-SNE of the joint PCA (02_celltypist/scripts/05_tsne.R).
#   - pseudotime = scv.tl.velocity_pseudotime on the same model (the function that
#     velociraptor::scvelo() calls), from 05_velocity_pseudotime.py.
# The pseudotime goes from 0 to 1 within each sample. Colour scale limited to the 5th and
# 95th percentile of each sample (values outside get the colour of the limit), as perc=[5, 95]
# in the scVelo tutorial (VelocityBasics, scv.pl.scatter).
# Run from the project folder, with the timepoints as arguments:
#   Rscript 05_rna_velocity/scripts/06_velociraptor_fig10.14.R 3 4
# With "stem_root" as first argument the pseudotime with root = stem cells and end = NSC is used:
#   Rscript 05_rna_velocity/scripts/06_velociraptor_fig10.14.R stem_root 0 1 2
# Output (results/velocity/velociraptor/): <sample>_velocity_pseudotime[_stem_root].png
# =============================================================================

suppressMessages({library(SingleCellExperiment); library(zellkonverter); library(velociraptor)
  library(scater); library(ggplot2)})

source("config.R")
models   <- file.path(velo_path, "scvelo")
out_path <- file.path(velo_path, "velociraptor")

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) > 0 && !grepl("^[0-9]+$", args[1])) args[1] else ""
tps <- as.integer(if (tag == "") args else args[-1])
suffix <- if (tag == "") "" else paste0("_", tag)

for (tp in tps) {
  ts <- read.csv(file.path(tsne_path, paste0("joint_T", tp), paste0("T", tp, "_tSNE.csv")))
  for (s in paste0(c("CTT", "PAT"), tp)) {
    sce <- readH5AD(file.path(models, paste0(s, "_scvelo.h5ad")))
    stopifnot(all(colnames(sce) %in% ts$cell))
    tsne <- as.matrix(ts[match(colnames(sce), ts$cell), c("tsne1", "tsne2")])
    colnames(tsne) <- c("TSNE1", "TSNE2")
    rownames(tsne) <- colnames(sce)
    reducedDim(sce, "TSNE") <- tsne

    pt <- read.csv(file.path(out_path, paste0(s, "_velocity_pseudotime", suffix, ".csv")))
    stopifnot(all(colnames(sce) %in% pt$cell))
    sce$pseudotime <- pt$velocity_pseudotime[match(colnames(sce), pt$cell)]

    embedded <- embedVelocity(reducedDim(sce, "TSNE"), sce)
    grid.df <- gridVectors(sce, embedded, use.dimred = "TSNE", resolution = 30)

    lims <- quantile(sce$pseudotime, c(0.05, 0.95), na.rm = TRUE)   # perc = [5, 95]
    p <- plotTSNE(sce, colour_by = "pseudotime", point_alpha = 0.3) +
      geom_segment(data = grid.df, mapping = aes(x = start.1, y = start.2, xend = end.1, yend = end.2),
                   arrow = arrow(length = unit(0.05, "inches"), type = "closed")) +
      ggtitle(paste(s, "velocity pseudotime", if (tag == "stem_root") "(root stem cells, end NSC)" else "")) +
      scale_colour_viridis_c(name = "pseudotime", limits = lims, oob = scales::squish)
    ggsave(file.path(out_path, paste0(s, "_velocity_pseudotime", suffix, ".png")), p,
           width = 9, height = 6, dpi = 200, bg = "white")
    message(s, " saved, colour limits ", round(lims[1], 3), " - ", round(lims[2], 3))
  }
}
