# =============================================================================
# CellTypist labels: are they supported by marker genes and by the probability?
# Author: Gaia Cocconcelli
# For every sample and every CellTypist label (best match, no stem threshold):
# - median probability of the best match, and % of cells just above the
# threshold (probability < 0.6)
# - % of cells with at least one count for pluripotency markers
# (POU5F1, NANOG, LIN28A, DPPA4, TDGF1) and neuroectoderm markers
# (PAX6, SOX1, OTX2, ZIC2, LHX2)
# Detection depends on sequencing depth, so compare labels within a sample.
# Output: results/celltypist/CellTypist_marker_check.csv
# =============================================================================

suppressMessages(library(SingleCellExperiment))

source("config.R")

pluri <- c("POU5F1", "NANOG", "LIN28A", "DPPA4", "TDGF1")
neuro <- c("PAX6", "SOX1", "OTX2", "ZIC2", "LHX2")
early_types <- c("stem cells 1", "stem cells 2", "stem cells 3", "neuroectodermal-like cells", "NSC")

res <- do.call(rbind, lapply(0:4, function(tp) {
  sce <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
  ct <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels.csv")))
  stopifnot(nrow(ct) == ncol(sce), all(ct$cell_index == seq_len(ncol(sce))),
            all(ct$sample == as.character(sce$sample)))
  lab <- ifelse(ct$CellTypist_probmatch == "Unassigned", "Unassigned",
                ifelse(ct$CellTypist %in% early_types, ct$CellTypist, "organoid-only types"))
  mk <- intersect(c(pluri, neuro), rownames(sce))
  if (length(mk) < 10) message("T", tp, ": markers not found: ", paste(setdiff(c(pluri, neuro), mk), collapse = ", "))
  det <- counts(sce)[mk, , drop = FALSE] > 0
  grp <- split(seq_len(ncol(sce)), paste(ct$sample, lab, sep = "|"))
  out <- do.call(rbind, lapply(grp, function(i) {
    data.frame(timepoint = paste0("T", tp), sample = ct$sample[i[1]], label = lab[i[1]],
               n_cells = length(i),
               median_prob = round(median(ct$CellTypist_prob[i]), 2),
               pct_prob_below_0.6 = round(100 * mean(ct$CellTypist_prob[i] < 0.6), 1),
               t(round(100 * rowMeans(det[, i, drop = FALSE]))))
  }))
  rm(sce); invisible(gc())
  out
}))
rownames(res) <- NULL
write.csv(res, file.path(ct_path, "CellTypist_marker_check.csv"), row.names = FALSE)

big <- res[res$n_cells >= 20, ] # tiny groups are not readable
print(big[, c("timepoint", "sample", "label", "n_cells", "median_prob", "pct_prob_below_0.6")], row.names = FALSE)
print(big[, c("sample", "label", intersect(pluri, names(big)))], row.names = FALSE)
print(big[, c("sample", "label", intersect(neuro, names(big)))], row.names = FALSE)
