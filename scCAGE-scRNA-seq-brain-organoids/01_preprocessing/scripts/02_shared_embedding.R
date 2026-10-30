# =============================================================================
# Shared embedding CTT vs PAT (one timepoint at a time), T0-T4 -> UMAP_A
# Author: Gaia Cocconcelli
# CTT and PAT of the same timepoint are put in the same space, without batch
# correction, so that the differences between the two samples are kept.
# UMAP_A is the map used by all the following analyses (cell type labels,
# regulatory element maps, RNA velocity, pseudotime).
# Steps:
# - multiBatchNorm: rescales the size factors of the two libraries;
# - HVGs: modelGeneVarByPoisson blocked by sample, top 10%;
# - cell cycle: CC.Difference regressed on both samples together with
# batchelor::regressBatches (design ~ CC.Difference, keep = 1): the intercept is
# kept, so the corrected values stay on the log-expression scale (they are not
# z-scored). A per-sample regression would remove exactly the CTT vs PAT
# differences we want to look at. The corrected values are used ONLY for the
# PCA: logcounts and counts are not changed (cell type annotation and RNA
# velocity need the counts);
# - PCA (50 PCs) on the corrected values of the HVGs, no scaling;
# - UMAP_A: n_neighbors = 80, min_dist = 0.3.
# Input : results/preprocessing/<LAB>_filtered.rds (01_preprocessing.R)
# Output: results/shared_embedding/joint_T<t>/T<t>_joint.rds
# results/shared_embedding/joint_T<t>/T<t>_shared_UMAP_A_by_sample.png
# =============================================================================

suppressMessages({library(scater); library(scran); library(batchelor); library(ggplot2)})

source("config.R")
keep_cols <- c("sum", "detected", "subsets_Mito_percent", "S.Score", "G2M.Score", "Phase", "CC.Difference")

for (tp in 0:4) {
  labs <- paste0(c("CTT", "PAT"), tp)
  tp_path <- file.path(joint_path, paste0("joint_T", tp))
  dir.create(tp_path, showWarnings = FALSE, recursive = TRUE)

  ## --- load the two samples, keep only counts + the metadata we need ---
  sce_list <- lapply(labs, function(lab) {
    s <- readRDS(file.path(pp_path, paste0(lab, "_filtered.rds")))
    out <- SingleCellExperiment(assays = list(counts = counts(s)), colData = colData(s)[, keep_cols])
    sizeFactors(out) <- sizeFactors(s)
    colnames(out) <- paste0(lab, "_", colnames(s))
    out$sample <- lab
    out$genotype <- sub("[0-9]$", "", lab)
    out
  })
  common_genes <- Reduce(intersect, lapply(sce_list, rownames))
  sce_list <- lapply(sce_list, function(s) s[common_genes, ])

  ## --- joint normalization (rescales size factors between the two libraries) ---
  sce_list <- multiBatchNorm(sce_list[[1]], sce_list[[2]])
  sce_joint <- cbind(sce_list[[1]], sce_list[[2]])
  rm(sce_list)
  sce_joint$sample <- factor(sce_joint$sample, levels = labs)

  ## --- HVGs, blocked by sample ---
  set.seed(101010)
  var_joint <- modelGeneVarByPoisson(sce_joint, block = sce_joint$sample)
  hvg_joint <- getTopHVGs(var_joint, prop = 0.1)

  ## --- CC.Difference regressed on both samples together, intercept kept ---
  design <- model.matrix(~ sce_joint$CC.Difference)
  cc_corrected <- regressBatches(sce_joint, design = design, keep = 1, subset.row = hvg_joint)

  ## --- PCA on the corrected HVGs, then UMAP_A ---
  set.seed(101010)
  reducedDim(sce_joint, "PCA") <- calculatePCA(assay(cc_corrected, "corrected"), ncomponents = 50,
                                               ntop = length(hvg_joint), scale = FALSE)
  rm(cc_corrected)
  set.seed(101010)
  sce_joint <- runUMAP(sce_joint, dimred = "PCA", n_neighbors = 80, min_dist = 0.3, name = "UMAP_A")

  p_overlay <- plotReducedDim(sce_joint, dimred = "UMAP_A", colour_by = "sample", point_size = 0.5) +
    ggtitle(paste0("T", tp, " - shared UMAP_A"))
  ggsave(file.path(tp_path, paste0("T", tp, "_shared_UMAP_A_by_sample.png")),
         p_overlay, width = 12, height = 10, dpi = 300, bg = "white")

  saveRDS(sce_joint, file.path(tp_path, paste0("T", tp, "_joint.rds")))
  message("=== T", tp, ": ", paste(labs, table(sce_joint$sample), collapse = " / "),
          " cells, ", length(hvg_joint), " HVGs, UMAP_A saved ===")
  rm(sce_joint); invisible(gc())
}
