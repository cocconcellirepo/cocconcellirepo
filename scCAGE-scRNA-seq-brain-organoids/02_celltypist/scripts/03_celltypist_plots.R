# =============================================================================
# CellTypist cell types on the shared UMAP_A (CTT | PAT) and composition, T0-T4
# Author: Gaia Cocconcelli
# Labels from 02_celltypist_annotation.py (model trained on the Kanton timecourse,
# label = Cluster), one label per cell. No new embedding: every cell keeps its
# UMAP_A coordinates (T<t>_joint.rds, only read).
# Block 1, T0-T2: stem threshold thr = 0.1. CellTypist gives every cell an
# independent probability for each Kanton type. A cell is labelled as a stem type
# when its highest probability among stem cells 1/2/3 is >= thr, even if another
# type has a higher probability; otherwise it keeps the CellTypist best match.
# thr was chosen on T0-T2 from the marker analysis.
# Groups: the 5 Kanton types of the early stages one by one, the types found in
# Kanton only in organoids (1-4 months) as "organoid-only types", Unassigned.
# Block 2, T3-T4 (day 13 and day 20): every cell keeps its CellTypist best match
# and every Kanton type has its own colour, in every sample (no grouping of the
# organoid types). No stem threshold (the stem cells are not forced). Unassigned
# = highest probability < p_min_late (0.4): at day 13-20 a part of the CTT cells
# sits between Kanton types (Neuroepithelium -> 1-month organoid), so the CellTypist
# default 0.5 leaves many of them out; 0.4 recovers some of these cells and is not
# lowered further. T0-T2 keep the default 0.5. Every run writes files with its
# own tag (stem threshold + p_min_late).
# Output: results/celltypist/T<t>_CellTypist_Cluster_labels_<tag>.csv
# results/shared_embedding/joint_T<t>/T<t>_shared_UMAP_A_CellTypist_<tag>.png
# results/shared_embedding/joint_T<t>/T<t>_shared_UMAP_A_CellTypist_<tag>_overlay.png (CTT and PAT in one panel)
# results/shared_embedding/Composition_10samples_CellTypist_<tag T3-T4>.png
# results/celltypist/CellTypist_organoid_types_<tag T3-T4>.csv
# =============================================================================

suppressMessages({library(SingleCellExperiment); library(ggplot2)})

source("config.R")
samp_lev <- paste0(c("CTT", "PAT"), rep(0:4, each = 2))

stem_types <- c("stem cells 1", "stem cells 2", "stem cells 3")
early_types <- c(stem_types, "neuroectodermal-like cells", "NSC")
pal_early <- c("stem cells 1" = "#1f77b4", "stem cells 2" = "#17becf", "stem cells 3" = "#9467bd",
               "neuroectodermal-like cells" = "#2ca02c", "NSC" = "#d62728")

# Kanton cells of each type at each stage (written by 01_export_for_celltypist.R)
kanton_stage <- read.csv(file.path(ct_path, "Kanton_Cluster_by_stage.csv"), row.names = 1, check.names = FALSE)
org_types <- setdiff(rownames(kanton_stage), early_types)

# CellTypist output of one timepoint, checked against the joint object
read_ct <- function(tp, sce) {
  ct <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels.csv")))
  pr <- read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_probabilities.csv")),
                 row.names = 1, check.names = FALSE)
  stopifnot(nrow(pr) == nrow(ct), all(rownames(pr) == as.character(ct$cell_index)))
  # the export wrote one row per cell, in the order of the joint object
  stopifnot(nrow(ct) == ncol(sce), all(ct$cell_index == seq_len(ncol(sce))),
            all(ct$sample == as.character(sce$sample)))
  ct$p_stem <- round(apply(pr[, stem_types], 1, max), 3)
  ct$best_stem <- stem_types[apply(pr[, stem_types], 1, which.max)]
  ct
}

# shared UMAP_A, one panel per sample on the same axes; the same map is also saved with
# CTT and PAT in one panel, as T<t>_shared_UMAP_A_by_sample.png (file name + "_overlay")
plot_map <- function(sce, type, pal, title, file) {
  lab <- factor(type, levels = names(pal)[names(pal) %in% type])
  df <- data.frame(u1 = reducedDim(sce, "UMAP_A")[, 1],
                   u2 = reducedDim(sce, "UMAP_A")[, 2],
                   sample = sce$sample, type = lab)
  set.seed(101010)
  df <- df[sample(nrow(df)), ] # random plotting order, so no group of cells is always on top
  p <- ggplot(df, aes(u1, u2, colour = type)) +
    geom_point(size = 0.4, alpha = 0.8) + facet_wrap(~ sample) +
    scale_colour_manual(values = pal, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMAP_A 1", y = "UMAP_A 2", title = title) +
    theme_classic()
  ggsave(file, p, width = 14, height = 7, dpi = 300, bg = "white")
  # one panel: smaller legend in one column, so that the map keeps most of the figure
  p_one <- p + facet_null() +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.5, alpha = 1))) +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
  ggsave(sub("\\.png$", "_overlay.png", file), p_one,
         width = 9, height = 8, dpi = 300, bg = "white")
}

all_labels <- list() # sample + type of every cell, for the composition
org_cells <- list() # cells assigned (probability >= 0.5) to an organoid type

## =====================================================================
## Block 1: T0-T2, stem threshold 0.1 (chosen labels)
## =====================================================================
thr <- 0.1
tag_early <- paste0("stemthr", thr)
pal_1 <- c(pal_early, "organoid-only types" = "#ff7f0e", "Unassigned" = "grey80")

for (tp in 0:2) {
  tp_path <- file.path(joint_path, paste0("joint_T", tp))
  sce <- readRDS(file.path(tp_path, paste0("T", tp, "_joint.rds")))
  ct <- read_ct(tp, sce)
  ct$CellTypist_stem_thr <- ifelse(ct$p_stem >= thr, ct$best_stem, ct$CellTypist)
  ct$type <- ifelse(ct$p_stem >= thr, ct$best_stem,
                    ifelse(ct$CellTypist_probmatch == "Unassigned", "Unassigned",
                           ifelse(ct$CellTypist %in% early_types, ct$CellTypist, "organoid-only types")))
  write.csv(ct, file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_", tag_early, ".csv")),
            row.names = FALSE)
  print(table(type = factor(ct$type, levels = names(pal_1)), sample = sce$sample))

  plot_map(sce, ct$type, pal_1, paste0("T", tp, " - shared UMAP, CellTypist, stem threshold ", thr),
           file.path(tp_path, paste0("T", tp, "_shared_UMAP_A_CellTypist_", tag_early, ".png")))
  all_labels[[paste0("T", tp)]] <- ct[, c("sample", "type")]
  oo <- ct$CellTypist_probmatch != "Unassigned" & ct$CellTypist %in% org_types
  org_cells[[paste0("T", tp)]] <- ct[oo, c("sample", "CellTypist")]
  rm(sce); invisible(gc())
}

## =====================================================================
## Block 2: T3-T4, every Kanton type with its own colour
## =====================================================================
thr_late <- NA # stem threshold; NA = best match only, no stem threshold
p_min_late <- 0.4 # Unassigned if the best match has probability < p_min_late (CellTypist default 0.5)
tag_late <- paste0(if (is.na(thr_late)) "nothr" else paste0("stemthr", thr_late), "_p", p_min_late)

# one colour per Kanton organoid type, fixed, chosen to be as far apart as
# possible from each other and from the early types (CIELAB distance)
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
pal_2 <- c(pal_early, org_cols[org_types], "Unassigned" = "grey80")

for (tp in 3:4) {
  tp_path <- file.path(joint_path, paste0("joint_T", tp))
  sce <- readRDS(file.path(tp_path, paste0("T", tp, "_joint.rds")))
  ct <- read_ct(tp, sce)
  # CellTypist_prob = highest probability of the cell, so prob < p_min_late is the same
  # as "no type at p_min_late" (with 0.5 it is CellTypist_probmatch == "Unassigned")
  ct$type <- ifelse(ct$CellTypist_prob < p_min_late, "Unassigned", ct$CellTypist)
  if (!is.na(thr_late)) ct$type <- ifelse(ct$p_stem >= thr_late, ct$best_stem, ct$type)
  write.csv(ct, file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_", tag_late, ".csv")),
            row.names = FALSE)
  print(table(type = factor(ct$type, levels = names(pal_2)), sample = sce$sample))

  plot_map(sce, ct$type, pal_2,
           paste0("T", tp, " - shared UMAP, CellTypist",
                  if (is.na(thr_late)) "" else paste0(", stem threshold ", thr_late)),
           file.path(tp_path, paste0("T", tp, "_shared_UMAP_A_CellTypist_", tag_late, ".png")))
  all_labels[[paste0("T", tp)]] <- ct[, c("sample", "type")]
  oo <- ct$CellTypist_prob >= p_min_late & ct$CellTypist %in% org_types
  org_cells[[paste0("T", tp)]] <- ct[oo, c("sample", "CellTypist")]
  rm(sce); invisible(gc())
}

## --- organoid types one by one, with their Kanton stages (all 10 samples) ---
oo_all <- do.call(rbind, org_cells)
oo_wide <- as.data.frame.matrix(table(factor(oo_all$CellTypist, levels = org_types),
                                      factor(oo_all$sample, levels = samp_lev)))
oo_out <- data.frame(Kanton_type = org_types,
                     kanton_stage[org_types, c("Organoid-1M", "Organoid-2M", "Organoid-4M")],
                     oo_wide, check.names = FALSE)
print(oo_out, row.names = FALSE)
write.csv(oo_out, file.path(ct_path, paste0("CellTypist_organoid_types_", tag_late, ".csv")), row.names = FALSE)

## =====================================================================
## Composition of the CellTypist cell types in all 10 samples
## =====================================================================
# Labels and groups of the maps above (T0-T2 block 1, T3-T4 block 2). Every
# cell counts once, whatever its number of UMIs.

pal_all <- c(pal_1[names(pal_1) != "Unassigned"], pal_2[setdiff(names(pal_2), names(pal_1))],
             "Unassigned" = "grey80")
comp <- do.call(rbind, all_labels)
stopifnot(all(comp$sample %in% samp_lev), all(comp$type %in% names(pal_all)))

comp_tab <- as.data.frame(table(sample = factor(comp$sample, levels = samp_lev),
                                type = factor(comp$type, levels = names(pal_all))))
comp_tab <- comp_tab[comp_tab$Freq > 0, ]
comp_tab$frac <- comp_tab$Freq / ave(comp_tab$Freq, comp_tab$sample, FUN = sum)
comp_tab$timepoint <- paste0("T", sub("^[A-Z]+", "", comp_tab$sample))
comp_tab$lab <- ifelse(comp_tab$frac >= 0.03, paste0(round(100 * comp_tab$frac), "%"), "")
print(xtabs(Freq ~ sample + type, comp_tab)) # cell counts, to check the bars

n_tab <- aggregate(Freq ~ sample + timepoint, comp_tab, sum)

p_comp <- ggplot(comp_tab, aes(sample, frac, fill = type, group = type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = lab), position = position_stack(vjust = 0.5), size = 3) +
  # n on two lines, so the labels of neighbouring bars do not overlap
  geom_text(data = n_tab, aes(sample, 1.05, label = paste0("n =\n", formatC(Freq, format = "d", big.mark = ","))),
            inherit.aes = FALSE, size = 2.6, lineheight = 0.9) +
  facet_grid(~ timepoint, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_all, name = NULL) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Cells",
       title = paste0("CellTypist cell types in the 10 samples (T0-T2 stem threshold ", thr,
                      if (is.na(thr_late)) ", T3-T4 best match)" else paste0(", T3-T4 stem threshold ", thr_late, ")"))) +
  theme_classic(base_size = 12)
ggsave(file.path(joint_path, paste0("Composition_10samples_CellTypist_", tag_late, ".png")),
       p_comp, width = 12, height = 6, dpi = 300, bg = "white")
