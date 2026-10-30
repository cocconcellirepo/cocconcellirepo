# =============================================================================
# UMAP built on the tCREs (SCAFE), CTT vs PAT, one timepoint at a time (T0-T4)
# Author: Gaia Cocconcelli
# Phase 4 of the tCRE analysis. The coordinates of the cells are computed from the
# tCRE counts only (no genes, no cell identities), in three versions:
# all tCREs -> distal only (D) -> proximal only (P)
# Common set all10: SCAFE aggregate on the pooled signal of the ten libraries
# (CTT and PAT at T0-T4), so the tCRE IDs are the same in every sample.
# Proximal = within 500 nt of an annotated gene TSS (SCAFE default).
# Same steps and parameters as the gene-based shared embedding (02_shared_embedding.R),
# so that only the features change:
# size factors on the tCRE counts of each sample (quickCluster + computeSumFactors)
# -> multiBatchNorm -> HVGs (modelGeneVarByPoisson blocked by sample, prop 0.1)
# -> CC.Difference regressed on both samples together (regressBatches, keep = 1)
# -> PCA (50 PCs, no scaling) -> UMAP (n_neighbors = 80, min_dist = 0.3), no batch correction
# CC.Difference is the cell cycle score of the cell (from the genes): it is regressed
# here too, otherwise the tCRE map and the gene map would differ also for the cell
# cycle and not only for the features.
# The maps are coloured by tCRE values only (UMI, activity, tCREs detected).
# Input : results/scafe/count_all10/<LAB>/<LAB>/matrix/
#         results/scafe/aggregate_all10/annotate/all10/bed/all10.CRE.annot.bed.gz
#         results/shared_embedding/joint_T<t>/T<t>_joint.rds (QC cells, sample, CC.Difference
#         and gene size factors only)
# Output: results/tCRE_embedding/ (overwritten at every run)
#         T<t>_tCRE_embedding.rds (PCA, UMAP and activity of each set)
#         T<t>_tCRE_UMAP_<set>_depth.png
#         T<t>_tCRE_UMAP_<set>_activity.png
#         T<t>_tCRE_UMAP_<set>_activity_overlay.png
# Set tp in Step 0, then run the steps one by one.
# =============================================================================

## ---- Step 0: setup ----
suppressMessages({library(Matrix); library(SingleCellExperiment); library(scater)
  library(scran); library(batchelor); library(ggplot2); library(patchwork)})

source("config.R")
out_path <- tcre_path
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

tp <- 0 # timepoint: 0, 1, 2, 3 or 4
labs <- paste0(c("CTT", "PAT"), tp)

set_title <- c(all = "all tCREs", distal = "distal tCREs", proximal = "proximal tCREs")

## ---- Step 1: QC cells of the joint object and SCAFE matrices ----
# from the joint object only what is needed, then it is removed (it is large)
sce_joint <- readRDS(file.path(joint_path, paste0("joint_T", tp), paste0("T", tp, "_joint.rds")))
stopifnot("CC.Difference" %in% names(colData(sce_joint)))
cells <- colnames(sce_joint)
smp <- as.character(sce_joint$sample)
cc_diff <- sce_joint$CC.Difference
sf_genes <- sizeFactors(sce_joint)
rm(sce_joint); invisible(gc())
table(smp)

# class of each tCRE (P or D) from the annotation of the common set;
# the rows of the count matrices are named chr_start_end_strand
annot <- read.delim(file.path(scafe_path, "aggregate_all10", "annotate", "all10", "bed",
                              "all10.CRE.annot.bed.gz"), header = FALSE)
cre_class <- setNames(sub(":.*", "", annot$V4), paste(annot$V1, annot$V2, annot$V3, annot$V6, sep = "_"))
rm(annot)
table(cre_class)

# one SCAFE count matrix (10x format), cells named <LAB>_<barcode> as in the joint object
read_scafe <- function(lab) {
  d <- file.path(scafe_path, "count_all10", lab, lab, "matrix")
  m <- as(readMM(file.path(d, "matrix.mtx")), "CsparseMatrix")
  rownames(m) <- read.delim(file.path(d, "genes.tsv"), header = FALSE)[, 1]
  colnames(m) <- paste0(lab, "_", readLines(file.path(d, "barcodes.tsv")))
  m
}
cre_list <- lapply(labs, read_scafe)
stopifnot(identical(rownames(cre_list[[1]]), rownames(cre_list[[2]]))) # same tCREs, same order
cre <- do.call(cbind, cre_list)
rm(cre_list)
stopifnot(all(rownames(cre) %in% names(cre_class)), all(cells %in% colnames(cre)))
cre <- cre[, cells] # only QC cells, same order as the joint

cls <- cre_class[rownames(cre)]
sets <- list(all = rep(TRUE, nrow(cre)), distal = cls == "D", proximal = cls == "P")
sapply(sets, sum)
# UMI per cell of each set: a cell with 0 UMI cannot get a size factor
sapply(sets, function(s) summary(colSums(cre[s, ])))

## ---- Step 2: one embedding per tCRE set, and the tCRE values per cell ----
tcre_embedding <- function(m) {
  # size factors on the tCRE counts of each sample, then rescaled between the two libraries
  sce_list <- lapply(labs, function(lab) {
    s <- SingleCellExperiment(assays = list(counts = m[, smp == lab]))
    s$sample <- lab
    set.seed(101010)
    s <- computeSumFactors(s, clusters = quickCluster(s))
    s
  })
  sce_list <- multiBatchNorm(sce_list[[1]], sce_list[[2]])
  s <- cbind(sce_list[[1]], sce_list[[2]])
  rm(sce_list)
  stopifnot(identical(colnames(s), cells))

  set.seed(101010)
  var_t <- modelGeneVarByPoisson(s, block = s$sample)
  hvg_t <- getTopHVGs(var_t, prop = 0.1)

  # CC.Difference regressed on both samples together, intercept kept (as 02_shared_embedding.R)
  design <- model.matrix(~ cc_diff)
  cc_corrected <- regressBatches(s, design = design, keep = 1, subset.row = hvg_t)

  set.seed(101010)
  reducedDim(s, "PCA") <- calculatePCA(assay(cc_corrected, "corrected"), ncomponents = 50,
                                       ntop = length(hvg_t), scale = FALSE)
  rm(cc_corrected)
  set.seed(101010)
  s <- runUMAP(s, dimred = "PCA", n_neighbors = 80, min_dist = 0.3, name = "UMAP")
  list(hvg = hvg_t, sf = sizeFactors(s), umi = colSums(counts(s)),
       pca = reducedDim(s, "PCA"), umap = reducedDim(s, "UMAP"))
}

# activity = log2(UMI of the set / size factor + 1), with the size factors of the joint
# object (from the genes). The tCRE size factors of
# tcre_embedding() are not used here: they come from the same counts, so
# UMI / size factor would be almost the same in every cell.
act <- data.frame(cell = cells, sample = smp, geno = substr(smp, 1, 3))
emb <- list()
for (st in names(sets)) {
  m <- cre[sets[[st]], ]
  act[[paste0(st, "_umi")]] <- colSums(m)
  act[[paste0(st, "_activity")]] <- log2(colSums(m) / sf_genes + 1)
  act[[paste0(st, "_detected")]] <- colSums(m > 0)
  emb[[st]] <- tcre_embedding(m)
  message("T", tp, " ", set_title[st], ": ", sum(sets[[st]]), " tCREs, ",
          length(emb[[st]]$hvg), " HVGs")
  rm(m); invisible(gc())
}
rm(cre); invisible(gc())

## ---- Step 3: save ----
saveRDS(list(cells = cells, sample = smp, activity = act, emb = emb),
        file.path(out_path, paste0("T", tp, "_tCRE_embedding.rds")))
message("=== T", tp, ": tCRE embeddings saved ===")

## ---- Step 4: depth check ----
# correlation of the first 5 PCs with log10(UMI of the set): a high value means
# that the PC follows the depth of the cell more than its state
sapply(emb, function(e) round(cor(e$pca[, 1:5], log10(e$umi)), 2))
# share of variance of the first 5 PCs
sapply(emb, function(e) round(attr(e$pca, "percentVar")[1:5], 2))

## ---- Step 5: depth map ----
# The tCRE UMAP coloured by the UMI of the set (depth check, see Step 4)
for (st in names(emb)) {
  df <- data.frame(u1 = emb[[st]]$umap[, 1], u2 = emb[[st]]$umap[, 2],
                   sample = smp, umi = log10(emb[[st]]$umi))
  set.seed(101010)
  df <- df[sample(nrow(df)), ] # random plotting order, so no group of cells is always on top

  p_umi <- ggplot(df, aes(u1, u2, colour = umi)) +
    geom_point(size = 0.2) + facet_wrap(~ sample) +
    scale_colour_viridis_c(name = "log10 UMI") +
    labs(x = "UMAP 1", y = "UMAP 2", title = paste0("T", tp, " - UMAP on ", set_title[st], ", UMI per cell")) +
    theme_classic(base_size = 10)
  ggsave(file.path(out_path, paste0("T", tp, "_tCRE_UMAP_", st, "_depth.png")),
         p_umi, width = 9, height = 6, dpi = 200, bg = "white")
}

## ---- Step 6: tCRE UMAPs coloured by tCRE activity ----
# Each map is coloured with its own set (distal values on the distal UMAP, ...)

# colour ramps of the overlaid maps (light = low, dark = high)
pal_geno <- list(CTT = colorRampPalette(c("#B5D4F4", "#042C53"))(100),
                 PAT = colorRampPalette(c("#F5C4B3", "#712B13"))(100))
to_colour <- function(v, geno, lim) {
  i <- round((v - lim[1]) / diff(lim) * 99) + 1
  i <- pmin(pmax(i, 1), 100)
  ifelse(geno == "CTT", pal_geno$CTT[i], pal_geno$PAT[i])
}
overlay_plot <- function(df, value, lim, title, legend_title) {
  df$col <- to_colour(df[[value]], df$geno, lim)
  p_map <- ggplot(df, aes(u1, u2, colour = col)) +
    geom_point(size = 0.3) + scale_colour_identity() +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_classic()
  steps <- seq(lim[1], lim[2], length.out = 100)
  bars <- data.frame(value = rep(steps, 2), geno = rep(c("CTT", "PAT"), each = 100),
                     col = c(pal_geno$CTT, pal_geno$PAT))
  p_bar <- ggplot(bars, aes(value, geno, fill = col)) +
    geom_tile() + scale_fill_identity() +
    labs(x = legend_title, y = NULL) +
    theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  p_map / p_bar + plot_layout(heights = c(12, 1))
}

for (st in names(emb)) {
  df <- act
  df$u1 <- emb[[st]]$umap[, 1] # coordinates of the tCRE UMAP
  df$u2 <- emb[[st]]$umap[, 2]
  set.seed(101010)
  df <- df[sample(nrow(df)), ] # random drawing order
  a <- paste0(st, "_activity")
  d <- paste0(st, "_detected")

  # CTT and PAT side by side
  p_act <- ggplot(df, aes(u1, u2, colour = .data[[a]])) +
    geom_point(size = 0.3) + facet_wrap(~ sample) +
    scale_colour_viridis_c(name = "log2(UMI /\nsize factor + 1)") +
    labs(x = "UMAP 1", y = "UMAP 2", title = paste0("T", tp, " - UMAP on ", set_title[st], ": activity")) +
    theme_classic()
  p_det <- ggplot(df, aes(u1, u2, colour = .data[[d]])) +
    geom_point(size = 0.3) + facet_wrap(~ sample) +
    scale_colour_viridis_c(name = "tCREs\ndetected", option = "magma") +
    labs(x = "UMAP 1", y = "UMAP 2", title = paste0("T", tp, " - UMAP on ", set_title[st], ": tCREs detected per cell")) +
    theme_classic()
  ggsave(file.path(out_path, paste0("T", tp, "_tCRE_UMAP_", st, "_activity.png")),
         p_act / p_det, width = 14, height = 12, dpi = 300, bg = "white")

  # CTT and PAT overlaid
  o_act <- overlay_plot(df, a, range(df[[a]]),
                        paste0("T", tp, " - UMAP on ", set_title[st], ": activity, CTT and PAT overlaid"),
                        "log2(UMI / size factor + 1)")
  o_det <- overlay_plot(df, d, range(df[[d]]),
                        paste0("T", tp, " - UMAP on ", set_title[st], ": tCREs detected per cell, CTT and PAT overlaid"),
                        "tCREs detected")
  # o_act and o_det are already patchworks (map + colour bar): wrap_elements() keeps
  # each one as a single block, otherwise patchwork appends o_det to the layout of o_act
  ggsave(file.path(out_path, paste0("T", tp, "_tCRE_UMAP_", st, "_activity_overlay.png")),
         wrap_elements(full = o_act) / wrap_elements(full = o_det),
         width = 8, height = 13, dpi = 300, bg = "white")
}
message("=== T", tp, ": tCRE UMAPs coloured by activity saved ===")
