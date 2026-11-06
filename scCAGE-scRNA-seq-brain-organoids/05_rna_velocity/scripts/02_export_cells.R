# =============================================================================
# Export of our cells for scVelo, T0-T4: barcode, sample, UMAP_A, CellTypist label
# Author: Gaia Cocconcelli
# The .loom files contain every cell called by Cell Ranger, here we write the cells
# kept after our QC (step 1) for each timepoint, with:
#   - barcode as in the joint objects (e.g. CTT0_AAACCTGAGATGTAAC-1)
#   - sample
#   - u1, u2 = coordinates on the shared UMAP_A of that timepoint
#   - celltype = column "type" of the CellTypist labels (step 2): T0-T2 stem threshold 0.1,
#     T3-T4 best match with Unassigned below 0.4
# The joint .rds files are only read, not changed.
# Output: results/velocity/cells/T<t>_cells_for_scvelo.csv
# =============================================================================

suppressMessages(library(SingleCellExperiment))

source("config.R")
out_path <- file.path(velo_path, "cells")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

for (tp in 0:4) {
  tag <- if (tp <= 2) "stemthr0.1" else "nothr_p0.4"
  sce <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
  ct  <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_", tag, ".csv")))
  # same cells, same order as the joint object
  stopifnot(nrow(ct) == ncol(sce), all(ct$cell_index == seq_len(ncol(sce))),
            all(ct$sample == as.character(sce$sample)), !anyNA(ct$barcode),
            all(ct$barcode == colnames(sce)))

  um <- reducedDim(sce, "UMAP_A")
  cells <- data.frame(barcode = ct$barcode, sample = ct$sample,
                      u1 = um[, 1], u2 = um[, 2], celltype = ct$type)
  write.csv(cells, file.path(out_path, paste0("T", tp, "_cells_for_scvelo.csv")), row.names = FALSE)
  print(table(celltype = cells$celltype, sample = cells$sample))
  rm(sce); invisible(gc())
}
