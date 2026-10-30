# =============================================================================
# Export for CellTypist: Kanton timecourse (cell types) and joint objects T0-T4
# Author: Gaia Cocconcelli
# CellTypist runs in Python, so the data leave R as MatrixMarket files.
# Reference: the whole Kanton timecourse (all stages, h409B2 + H9) with label =
# Cluster, the cell types of Kanton et al. 2019.
# No gene filter: the model can use every gene shared with our data.
# ref_kanton_cluster : raw counts + Cluster of every reference cell
# joint_T<t> : raw counts of CTT<t> + PAT<t>, all genes (T<t>_joint.rds)
# genes_common.txt : genes present in the reference and in all five joint objects
# Kanton_Cluster_by_stage.csv : cells of each cell type at each Kanton stage
# Only raw counts are exported: the normalization CellTypist wants (counts
# scaled to 10,000 per cell, then log1p) is done in Python.
# Output: results/celltypist/
# =============================================================================

suppressMessages({library(SingleCellExperiment); library(Seurat); library(Matrix)})

source("config.R")
dir.create(ct_path, showWarnings = FALSE, recursive = TRUE)

# one folder per dataset: counts (genes x cells), gene names, one row per cell
export_counts <- function(m, cells, dir) {
  dir.create(dir, showWarnings = FALSE)
  m <- Matrix(m, sparse = TRUE)
  if (!all(m@x == round(m@x))) stop(basename(dir), ": the matrix does not contain raw counts")
  writeMM(m, file.path(dir, "counts.mtx"))
  writeLines(rownames(m), file.path(dir, "genes.txt"))
  write.csv(cells, file.path(dir, "cells.csv"), row.names = FALSE)
  message(basename(dir), ": ", nrow(m), " genes x ", ncol(m), " cells")
}

## --- reference: Kanton timecourse, all stages, both lines ---
ref_kanton <- readRDS(kanton_rds)
sce_kanton <- as.SingleCellExperiment(ref_kanton)
rm(ref_kanton); invisible(gc())

# at which stages each cell type is found
stage_lev <- c("iPSCs", "EB", "Neuroectoderm", "Neuroepithelium",
               "Organoid-1M", "Organoid-2M", "Organoid-4M")
tab_cl <- table(Cluster = as.character(sce_kanton$Cluster),
                stage_group = factor(sce_kanton$stage_group, levels = stage_lev))
print(tab_cl)
write.csv(as.data.frame.matrix(tab_cl), file.path(ct_path, "Kanton_Cluster_by_stage.csv"))

export_counts(counts(sce_kanton),
              data.frame(cell = colnames(sce_kanton),
                         label = as.character(sce_kanton$Cluster),
                         stage_group = as.character(sce_kanton$stage_group),
                         line = as.character(sce_kanton$line)),
              file.path(ct_path, "ref_kanton_cluster"))
genes_common <- rownames(sce_kanton)
rm(sce_kanton); invisible(gc())

## --- joint objects T0-T4 ---
for (tp in 0:4) {
  sce_joint <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
  genes_common <- intersect(genes_common, rownames(sce_joint))

  # cell_index = position in the joint object, used to bring the labels back
  cells <- data.frame(cell_index = seq_len(ncol(sce_joint)),
                      barcode = colnames(sce_joint),
                      sample = as.character(sce_joint$sample))
  export_counts(counts(sce_joint), cells, file.path(ct_path, paste0("joint_T", tp)))
  rm(sce_joint); invisible(gc())
}
writeLines(genes_common, file.path(ct_path, "genes_common.txt"))
message(length(genes_common), " genes shared by the reference and T0-T4")
