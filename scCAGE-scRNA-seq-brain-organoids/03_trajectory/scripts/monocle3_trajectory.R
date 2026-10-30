# =============================================================================
# Trajectory and pseudotime with monocle3, one genotype at a time (T0 -> T4)
# Author: Gaia Cocconcelli
# Guideline: Sanbomics monocle3 tutorial
# (github.com/mousepixels/sanbomics_scripts/blob/main/monocle3_tutorial.Rmd):
# the samples of all the days are merged in ONE Seurat object (no integration),
# NormalizeData -> FindVariableFeatures (vst, 2,000 genes) -> ScaleData ->
# RunPCA (30 PCs) -> FindNeighbors / RunUMAP (dims 1:30), then monocle3:
# cluster_cells -> learn_graph(use_partition = FALSE) -> order_cells -> plot_cells.
# Here one object per genotype (CTT0 ... CTT4, or PAT0 ... PAT4), so the
# pseudotime of a genotype runs over the five timepoints. CTT and PAT are in two
# different maps: their pseudotime values are not on the same scale.
# Differences from the tutorial:
# - cells and counts from results/preprocessing/<LAB>_filtered.rds (our QC, the same
# cells as the joint objects), so the tutorial QC (nFeature quantiles, mt < 15%)
# is not repeated;
# - JoinLayers() after merge (Seurat 5 keeps one layer per sample; the tutorial
# was written with one matrix);
# - conversion to monocle3 with new_cell_data_set() + the Seurat PCA and UMAP,
# instead of SeuratWrappers::as.cell_data_set() (not installed);
# - root fixed in the code instead of chosen by clicking (order_cells() without
# arguments): the principal graph node closest to most of the T0 stem cells
# (CellTypist stem cells 1 + 2), with the helper of the monocle3 documentation
# (get_earliest_principal_node, "Constructing single-cell trajectories").
# As in monocle3, the graph and the pseudotime are computed on the UMAP.
# Timepoint and library are the same thing here (n = 1 per sample): part of the
# separation between days can be technical.
# Colours: CellTypist labels as in CellTypist_plots.R (T0-T2 stem threshold 0.1,
# organoid types grouped; T3-T4 best match, every Kanton type with its own colour,
# Unassigned < 0.4).
# Step 8 draws the path on the principal graph as numbered steps (see there).
# Output: results/monocle3/<geno>/
# <geno>_monocle3_UMAP_day.png, _CellTypist.png, _pseudotime.png, _pseudotime_con_frecce.png
# <geno>_monocle3_pseudotime.csv (one row per cell), <geno>_monocle3_cds.rds
# Run the steps one by one, first with geno <- "CTT", then with geno <- "PAT".
# =============================================================================

## ---- Step 0: setup ----
suppressMessages({library(SingleCellExperiment); library(Seurat); library(monocle3)
  library(ggplot2); library(patchwork)})

source("config.R")
geno <- "CTT" # "CTT" or "PAT"
out_path <- file.path(pt_path, geno)
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
ct_tag <- c("stemthr0.1", "stemthr0.1", "stemthr0.1", "nothr_p0.4", "nothr_p0.4") # labels of T0 ... T4

# palette of CellTypist_plots.R (T0-T2 groups + T3-T4 Kanton types)
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
pal <- c(pal_early, "organoid-only types" = "#ff7f0e", org_cols[org_types], "Unassigned" = "grey80")
day_pal <- c(T0 = "#c7e9b4", T1 = "#7fcdbb", T2 = "#41b6c4", T3 = "#225ea8", T4 = "#0c2c84") # light to dark, five well separated steps

# root of the trajectory: the principal graph node closest to most of the given cells
# (helper of the monocle3 documentation, used in Steps 5 and 8)
get_earliest_principal_node <- function(cds, cells) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cells, ]))))]
}

## ---- Step 1: one Seurat object per sample (QC cells), CellTypist labels ----
# CellTypist labels of the five timepoints, by cell name (<LAB>_<barcode>)
ct <- do.call(rbind, lapply(0:4, function(tp)
  read.csv(file.path(ct_path, paste0("T", tp, "_CellTypist_Cluster_labels_", ct_tag[tp + 1], ".csv")))[, c("barcode", "type")]))
ct_type <- setNames(ct$type, ct$barcode)

obj_list <- lapply(0:4, function(tp) {
  lab <- paste0(geno, tp)
  s <- readRDS(file.path(pp_path, paste0(lab, "_filtered.rds")))
  m <- counts(s)
  colnames(m) <- paste0(lab, "_", colnames(s)) # same cell names as the joint objects
  o <- CreateSeuratObject(m, project = lab, min.cells = 0, min.features = 200)
  o$sample <- lab
  o$day <- paste0("T", tp)
  o
})
sapply(obj_list, ncol)

## ---- Step 2: merge the five days, Seurat workflow of the tutorial ----
data <- merge(obj_list[[1]], y = obj_list[2:5])
data <- JoinLayers(data)
rm(obj_list); invisible(gc())
stopifnot(all(colnames(data) %in% names(ct_type)))
data$celltypist <- unname(ct_type[colnames(data)])
stopifnot(all(data$celltypist %in% names(pal)))
table(celltypist = data$celltypist, day = data$day)

data <- NormalizeData(object = data, verbose = FALSE)
data <- FindVariableFeatures(object = data, nfeatures = 2000, verbose = FALSE, selection.method = "vst")
data <- ScaleData(data, verbose = FALSE)
data <- RunPCA(data, npcs = 30, verbose = FALSE)
data <- FindNeighbors(data, dims = 1:30)
data <- RunUMAP(data, reduction = "pca", dims = 1:30)

## ---- Step 3: monocle3 object (counts, PCA and UMAP from Seurat) ----
counts_m <- LayerData(data, assay = "RNA", layer = "counts")
cds <- new_cell_data_set(counts_m,
                         cell_metadata = data@meta.data,
                         gene_metadata = data.frame(gene_short_name = rownames(counts_m), row.names = rownames(counts_m)))
reducedDim(cds, "PCA") <- Embeddings(data, "pca")[colnames(cds), ]
reducedDim(cds, "UMAP") <- Embeddings(data, "umap")[colnames(cds), ]

## ---- Step 4: clusters, principal graph ----
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)
# cluster_cells: neighbour graph on the UMAP (k = 20: each cell is linked to its 20 nearest
# cells), then Leiden: groups of cells well connected to each other and little to the rest

## ---- Step 5: pseudotime, root = node closest to most T0 stem cells ----
root_cells <- which(cds$day == "T0" & cds$celltypist %in% c("stem cells 1", "stem cells 2"))
length(root_cells)
root_node <- get_earliest_principal_node(cds, root_cells)
cds <- order_cells(cds, root_pr_nodes = root_node)

pt <- data.frame(cell = colnames(cds), sample = cds$sample, day = cds$day, celltypist = cds$celltypist,
                 pseudotime = pseudotime(cds),
                 umap1 = reducedDim(cds, "UMAP")[, 1], umap2 = reducedDim(cds, "UMAP")[, 2])
# does the pseudotime follow the real days?
print(aggregate(pseudotime ~ day, pt, median))
print(round(tapply(pt$pseudotime, list(pt$celltypist, pt$day), median), 1))
write.csv(pt, file.path(out_path, paste0(geno, "_monocle3_pseudotime.csv")), row.names = FALSE)

## ---- Step 6: maps ----
# Drawn with ggplot from the table of Step 5 (pt) instead of plot_cells(), which draws a
# black point under every cell (outline); small legend in one column, so that the map keeps
# most of the figure (as the one-panel maps of CellTypist_plots.R)
g_full <- principal_graph(cds)[["UMAP"]]
Y_full <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
el_full <- igraph::ends(g_full, igraph::E(g_full))
graph_df <- data.frame(x = Y_full[el_full[, 1], 1], y = Y_full[el_full[, 1], 2],
                       xend = Y_full[el_full[, 2], 1], yend = Y_full[el_full[, 2], 2])
set.seed(101010)
df6 <- pt[sample(nrow(pt)), ] # random plotting order, so no group of cells is always on top
df6$celltypist <- factor(df6$celltypist, levels = names(pal)[names(pal) %in% df6$celltypist])
map6 <- function(colour_by) {
  ggplot(df6, aes(umap1, umap2, colour = .data[[colour_by]])) +
    geom_point(size = 0.3) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic() +
    theme(legend.text = element_text(size = 7), legend.key.size = unit(0.35, "cm"))
}
p_day <- map6("day") + scale_colour_manual(values = day_pal, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) + ggtitle(paste(geno, "- timepoint"))
p_ct <- map6("celltypist") + scale_colour_manual(values = pal, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.5))) + ggtitle(paste(geno, "- CellTypist"))
# trajectory graph drawn twice (white below, dark on top) so it shows on light and dark cells
p_pt <- map6("pseudotime") + scale_colour_viridis_c(option = "plasma", name = "pseudotime") +
  geom_segment(data = graph_df, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE,
               colour = "white", linewidth = 1.1) +
  geom_segment(data = graph_df, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE,
               colour = "grey15", linewidth = 0.45) +
  ggtitle(paste(geno, "- pseudotime"))
for (x in list(list("day", p_day), list("CellTypist", p_ct), list("pseudotime", p_pt))) {
  ggsave(file.path(out_path, paste0(geno, "_monocle3_UMAP_", x[[1]], ".png")), x[[2]],
         width = 9, height = 7, dpi = 300, bg = "white")
}
message("=== ", geno, ": maps saved ===")

## ---- Step 7: save ----
saveRDS(cds, file.path(out_path, paste0(geno, "_monocle3_cds.rds")))
message("=== ", geno, ": monocle3 saved ===")

## ---- Step 8: path on the principal graph as numbered steps ----
# The principal graph has hundreds of nodes. To draw the path as a few numbered steps:
# - pseudotime of each node = distance along the graph from the root (edge length on
#   the UMAP), the same quantity order_cells() gives to the cells;
# - short side branches (shorter than min_branch of the maximum pseudotime) are
#   removed; the root, the branch points and the ends of the remaining branches are kept;
# - along each branch a node is added every `step` of pseudotime (step = maximum / n_steps);
# - nodes are numbered branch by branch, going from the root (depth first, the shorter
#   branch first), arrows go towards higher pseudotime;
# - the boxes name the monocle3 clusters: days with >= 20% of the cells of the cluster
#   and the most frequent CellTypist label (stem cells 1/2/3 = "stem"); clusters with
#   the same name are drawn as one box, only boxes with >= 3% of the cells.
# Step 8 can also be run alone (after Step 0), e.g. after a restart of R: the cds and the
# pseudotime table of this genotype are then read from the files of Steps 5 and 7
if (!exists("cds") || !exists("pt") || !identical(pt$cell, colnames(cds)) || !all(paste0(geno, 0:4) %in% pt$sample)) {
  cds <- readRDS(file.path(out_path, paste0(geno, "_monocle3_cds.rds")))
  pt <- read.csv(file.path(out_path, paste0(geno, "_monocle3_pseudotime.csv")))
  stopifnot(identical(pt$cell, colnames(cds)), all(paste0(geno, 0:4) %in% pt$sample))
  root_node <- get_earliest_principal_node(cds, which(cds$day == "T0" & cds$celltypist %in% c("stem cells 1", "stem cells 2")))
}
n_steps <- 8      # about n_steps numbered nodes along the longest path
min_branch <- 0.2 # side branches shorter than 20% of the maximum pseudotime are not drawn

g <- principal_graph(cds)[["UMAP"]]
Y <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) # node coordinates on the UMAP
stopifnot(all(igraph::V(g)$name %in% rownames(Y)))
el <- igraph::ends(g, igraph::E(g))
igraph::E(g)$len <- sqrt(rowSums((Y[el[, 1], ] - Y[el[, 2], ])^2))
node_pt <- igraph::distances(g, v = root_node, weights = igraph::E(g)$len)[1, ]
# check: pseudotime of the cells vs pseudotime of their closest node (expected close to 1)
cv <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), 1]
print(round(cor(pt$pseudotime, node_pt[igraph::V(g)$name[cv]]), 3))
# monocle3 closes loops in the graph (learn_graph, close_loop = TRUE): only the
# shortest-path tree from the root is kept, so every node is reached by one path
# (the distances from the root do not change); the grey graph on the map stays complete
sp <- igraph::shortest_paths(g, from = root_node, weights = igraph::E(g)$len, output = "epath")$epath
tree_e <- unique(unlist(lapply(sp, as.vector)))
g_tree <- igraph::delete_edges(g, setdiff(seq_len(igraph::ecount(g)), tree_e))
stopifnot(igraph::ecount(g_tree) == igraph::vcount(g_tree) - 1)
message("graph: ", igraph::vcount(g), " nodes, ", igraph::ecount(g) - igraph::ecount(g_tree), " loop edge(s) removed")

# keep the main branches: the tree is oriented from the root and, at every branch point,
# a side branch is drawn only if it goes at least min_branch (fraction of the maximum
# pseudotime) beyond the branch point; the child that goes furthest is always kept, so a
# long branch is never cut into pieces
thr <- min_branch * max(node_pt)
et <- igraph::ends(g_tree, igraph::E(g_tree))
up <- node_pt[et[, 1]] < node_pt[et[, 2]]
dg <- igraph::graph_from_edgelist(cbind(ifelse(up, et[, 1], et[, 2]), ifelse(up, et[, 2], et[, 1])))
# h = highest pseudotime reached by each node and the nodes below it (a child has a higher
# pseudotime than its parent, so going down the pseudotime the children come first)
h <- node_pt[igraph::V(dg)$name]
for (v in names(sort(h, decreasing = TRUE))) {
  ch <- names(igraph::neighbors(dg, v, mode = "out"))
  if (length(ch) > 0) h[v] <- max(h[v], h[ch])
}
drop <- character(0)
for (p in igraph::V(dg)$name) {
  ch <- names(igraph::neighbors(dg, p, mode = "out"))
  if (length(ch) < 2) next
  side <- setdiff(ch, ch[which.max(h[ch])])
  for (cc in side[h[side] - node_pt[p] < thr]) {
    drop <- c(drop, names(igraph::subcomponent(dg, cc, mode = "out")))
  }
}
tr <- igraph::delete_vertices(g_tree, unique(drop))
deg_t <- igraph::degree(tr)
key <- names(deg_t)[deg_t != 2 | names(deg_t) == root_node] # root, branch points, ends
step <- max(node_pt[names(deg_t)]) / n_steps
message("after pruning: ", igraph::vcount(tr), " nodes, ", sum(deg_t >= 3), " branch points, ", sum(deg_t == 1) - (deg_t[root_node] == 1), " branch ends")

# nodes from key node k through its neighbour nb to the next key node
walk <- function(k, nb) {
  path <- c(k, nb)
  while (!(path[length(path)] %in% key)) {
    path <- c(path, setdiff(names(igraph::neighbors(tr, path[length(path)])), path[length(path) - 1]))
  }
  path
}
# highest pseudotime reachable from k through nb (to number the shorter branch first)
down_max <- function(k, nb) {
  comp <- igraph::components(igraph::delete_edges(tr, paste0(k, "|", nb)))$membership
  max(node_pt[names(comp)[comp == comp[nb]]])
}
# numbered nodes of one branch: its first and last node + one every `step`
pick_nodes <- function(path) {
  p <- node_pt[path]
  targets <- if (p[length(p)] - p[1] >= 1.5 * step) seq(p[1] + step, p[length(p)] - step / 2, by = step) else numeric(0)
  unique(c(path[1], path[vapply(targets, function(t) which.min(abs(p - t)), integer(1))], path[length(path)]))
}
num <- root_node     # numbered nodes, in order
segs <- list()       # numbered nodes of each branch
visit <- function(k) {
  nbs <- names(igraph::neighbors(tr, k))
  nbs <- nbs[node_pt[nbs] > node_pt[k]]
  nbs <- nbs[order(vapply(nbs, function(nb) down_max(k, nb), numeric(1)))]
  for (nb in nbs) {
    path <- walk(k, nb)
    sel <- pick_nodes(path)
    segs[[length(segs) + 1]] <<- sel
    num <<- c(num, setdiff(sel, num))
    visit(path[length(path)])
  }
}
visit(root_node)
message("numbered nodes: ", length(num), ", branches: ", length(segs))

seg_cols <- rep(c("#F28E2B", "#1B9E77", "#1F9BD7", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F"),
                length.out = length(segs))
nodes <- data.frame(n = seq_along(num), x = Y[num, 1], y = Y[num, 2])
# arrows between consecutive numbered nodes, shortened so they stop at the circles
r <- 0.03 * diff(range(pt$umap1))
arr <- do.call(rbind, lapply(seq_along(segs), function(i) {
  s <- segs[[i]]
  a <- Y[s[-length(s)], , drop = FALSE]; b <- Y[s[-1], , drop = FALSE]
  d <- sqrt(rowSums((b - a)^2)); ux <- (b[, 1] - a[, 1]) / d; uy <- (b[, 2] - a[, 2]) / d
  data.frame(x = a[, 1] + r * ux, y = a[, 2] + r * uy, xend = b[, 1] - r * ux, yend = b[, 2] - r * uy,
             col = seg_cols[i])
}))
# legend text: the numbers of each branch (on more lines if long) and its pseudotime range
seq_txt <- vapply(segs, function(s) paste(strwrap(paste(match(s, num), collapse = " → "), width = 30),
                                          collapse = "\n"), character(1))
seg_tab <- data.frame(
  text = c("1  root (T0 stem cells)",
           paste0(seq_txt, "\npseudotime ", round(node_pt[vapply(segs, function(s) s[1], character(1))]), "-",
                  round(node_pt[vapply(segs, function(s) s[length(s)], character(1))]))),
  col = c("black", seg_cols))
n_lines <- lengths(strsplit(seg_tab$text, "\n"))
seg_tab$y <- -(cumsum(n_lines) - n_lines + 0.8 * seq_along(n_lines)) # top of each entry
print(seg_tab$text)

# boxes on the monocle3 clusters
short <- c("stem cells 1" = "stem", "stem cells 2" = "stem", "stem cells 3" = "stem",
           "neuroectodermal-like cells" = "NE-like")
pt$lab_short <- ifelse(pt$celltypist %in% names(short), short[pt$celltypist], pt$celltypist)
pt$cluster <- as.character(clusters(cds)[pt$cell])
cl_lab <- vapply(split(seq_len(nrow(pt)), pt$cluster), function(i) {
  day_frac <- table(pt$day[i]) / length(i)
  paste0(paste(names(day_frac)[day_frac >= 0.2], collapse = " + "), " (",
         names(which.max(table(pt$lab_short[i]))), ")")
}, character(1))
pt$box <- unname(cl_lab[pt$cluster])
boxes <- aggregate(cbind(umap1, umap2) ~ box, pt, median)
boxes <- boxes[boxes$box %in% names(which(table(pt$box) >= 0.03 * nrow(pt))), ]
print(table(pt$box))

edges <- data.frame(x = Y[el[, 1], 1], y = Y[el[, 1], 2], xend = Y[el[, 2], 1], yend = Y[el[, 2], 2])
set.seed(101010)
df <- pt[sample(nrow(pt)), ] # random plotting order
p_map <- ggplot() +
  geom_point(data = df, aes(umap1, umap2, colour = pseudotime), size = 0.3) +
  scale_colour_viridis_c(option = "plasma", name = "pseudotime") +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), colour = "grey30", linewidth = 0.5) +
  geom_segment(data = arr, aes(x, y, xend = xend, yend = yend), colour = arr$col, linewidth = 1.4,
               arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
  geom_point(data = nodes, aes(x, y), shape = 21, size = 6, fill = "white", colour = "black", stroke = 0.8) +
  geom_text(data = nodes, aes(x, y, label = n), size = 3, fontface = "bold") +
  # boxes pushed away from the numbered nodes too (nodes passed with an empty label)
  ggrepel::geom_label_repel(data = rbind(data.frame(x = boxes$umap1, y = boxes$umap2, label = boxes$box),
                                         data.frame(x = nodes$x, y = nodes$y, label = "")),
                            aes(x, y, label = label), size = 3.2, fontface = "bold", point.size = 5,
                            label.r = unit(0, "pt"), min.segment.length = Inf, seed = 101010) +
  labs(x = "UMAP 1", y = "UMAP 2", title = paste(geno, "- pseudotime")) +
  theme_classic()
p_leg <- ggplot(seg_tab, aes(0, y, label = text)) +
  geom_text(colour = seg_tab$col, hjust = 0, vjust = 1, size = 3.2, lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(min(seg_tab$y) - max(n_lines) - 15, 0)) +
  theme_void()
ggsave(file.path(out_path, paste0(geno, "_monocle3_UMAP_pseudotime_con_frecce.png")),
       p_map + p_leg + plot_layout(widths = c(3, 1)), width = 13, height = 8, dpi = 300, bg = "white")
message("=== ", geno, ": numbered path saved ===")
