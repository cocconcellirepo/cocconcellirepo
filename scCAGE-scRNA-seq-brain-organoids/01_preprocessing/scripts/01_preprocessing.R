# =============================================================================
# Pre-processing pipeline - all samples (T0-T4), one function
# Author: Gaia Cocconcelli
# The same function for every library, so that all the samples go through the
# same steps and the only sample-specific parameters are written in samples.csv.
# Steps: Cell Ranger filtered matrix -> QC with fixed thresholds (cells with
# sum < 1500, detected < 2000 or mito > 10% are discarded) -> scDblFinder with the
# doublet rate of the sample (dbr, dbr.sd) -> removal of mito / ribo /
# nuclear-OXPHOS genes -> size factors (quickCluster + computeSumFactors) and
# logcounts -> cell cycle scores (Seurat; CC.Difference = S.Score - G2M.Score).
# HVGs, cell cycle regression, PCA and UMAP are done on the two samples of a
# timepoint together in 02_shared_embedding.R, so they are not repeated here.
# Doublets: dbr = 0.15, dbr.sd = 0.01 for every sample. The score of every
# QC-passed cell is saved, so that the rate can be checked and changed sample by sample.
# Output (results/preprocessing/):
# <LAB>_filtered.rds QC-passed singlets: counts, logcounts, size factors, CC scores
# <LAB>_01_QC_3panel.png, <LAB>_02_QC_UMIvsMito_scatter.png
# <LAB>_scDblFinder_summary.txt singlets / doublets
# <LAB>_scDblFinder_cells.csv score and class of every QC-passed cell
# preprocessing_summary.csv cells kept at each step, per sample
# =============================================================================

suppressMessages({library(Seurat); library(scater); library(scuttle); library(scran)
  library(scDblFinder); library(ggplot2); library(gridExtra)})

source("config.R")
out_path <- pp_path
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
stopifnot(all(dir.exists(samples$data_dir)))
samples[, c("label", "dbr", "dbr.sd")]

process_sample <- function(row) {
  label <- row$label
  message("=== ", label, " ===")

  ## QC (fixed thresholds)
  counts <- Read10X(data.dir = row$data_dir)
  sce <- as.SingleCellExperiment(CreateSeuratObject(counts = counts, project = label,
                                                    min.cells = 0, min.features = 0))
  sce <- addPerCellQCMetrics(sce, subsets = list(Mito = grepl("^MT-", rownames(sce))))
  sce$discard <- sce$sum < 1500 | sce$detected < 2000 | sce$subsets_Mito_percent > 10
  n_barcodes <- ncol(sce)
  message(label, ": QC discards ", sum(sce$discard), " / ", n_barcodes, " cells")

  ggsave(file.path(out_path, paste0(label, "_01_QC_3panel.png")),
         arrangeGrob(
           plotColData(sce, x = "orig.ident", y = "sum", colour_by = "discard") + scale_y_log10() + ggtitle(paste(label, "- Total UMI")),
           plotColData(sce, x = "orig.ident", y = "detected", colour_by = "discard") + scale_y_log10() + ggtitle(paste(label, "- Detected features")),
           plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "discard") + ggtitle(paste(label, "- Mito percent")),
           ncol = 1),
         width = 12, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(out_path, paste0(label, "_02_QC_UMIvsMito_scatter.png")),
         plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "discard") + scale_x_log10() +
           labs(title = paste(label, "- Total UMI vs Mito percent"), x = "Total UMI (log10)", y = "Mito percent"),
         width = 12, height = 10, dpi = 300, bg = "white")
  sce <- sce[, !sce$discard]
  n_qc <- ncol(sce)

  ## doublets (scDblFinder, rate of the sample)
  set.seed(101010)
  sce <- scDblFinder(sce, clusters = TRUE, dbr = row$dbr, dbr.sd = row$dbr.sd, iter = 3)
  write.csv(data.frame(cell = colnames(sce), score = sce$scDblFinder.score, class = sce$scDblFinder.class,
                       sum = sce$sum, detected = sce$detected, mito = sce$subsets_Mito_percent),
            file.path(out_path, paste0(label, "_scDblFinder_cells.csv")), row.names = FALSE)
  dbl <- as.data.frame(table(class = sce$scDblFinder.class))
  dbl$percent <- paste0(round(100 * dbl$Freq / sum(dbl$Freq), 2), "%")
  write.table(dbl, file.path(out_path, paste0(label, "_scDblFinder_summary.txt")),
              quote = FALSE, row.names = FALSE, sep = "\t")
  sce <- sce[, sce$scDblFinder.class == "singlet"]

  ## gene filtering: mito / ribo / nuclear-OXPHOS
  drop_genes <- grepl("^MT-", rownames(sce)) | grepl("^RP[SL]", rownames(sce)) |
    grepl("^ATP5|^COX[0-9]|^NDUF|^UQCR|^SDHA|^SDHB|^SDHC|^SDHD|^ROMO1$", rownames(sce))
  sce <- sce[!drop_genes, ]

  ## normalization (deconvolution size factors)
  set.seed(101010)
  clust <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = clust)
  sce <- logNormCounts(sce)

  ## cell cycle scores (the regression is done later, in the shared embedding)
  seu_cc <- as.Seurat(sce, counts = "counts", data = "logcounts")
  seu_cc <- CellCycleScoring(seu_cc, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  sce$S.Score <- seu_cc$S.Score
  sce$G2M.Score <- seu_cc$G2M.Score
  sce$Phase <- seu_cc$Phase
  sce$CC.Difference <- sce$S.Score - sce$G2M.Score
  rm(seu_cc)

  saveRDS(sce, file.path(out_path, paste0(label, "_filtered.rds")))
  message("=== ", label, ": ", ncol(sce), " cells saved ===")
  data.frame(label = label, barcodes = n_barcodes, after_QC = n_qc, singlets = ncol(sce),
             doublets_pct = round(100 * (n_qc - ncol(sce)) / n_qc, 2))
}

## --- run: all the samples (or only some, e.g. run_labels <- c("CTT3", "PAT3")) ---
run_labels <- samples$label
pp_summary <- do.call(rbind, lapply(run_labels, function(l) {
  res <- process_sample(samples[samples$label == l, ])
  invisible(gc())
  res
}))
print(pp_summary, row.names = FALSE)
write.csv(pp_summary, file.path(out_path, "preprocessing_summary.csv"), row.names = FALSE)
