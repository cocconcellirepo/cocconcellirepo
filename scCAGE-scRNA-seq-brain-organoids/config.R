# Paths and sample table used by all the R scripts (edit here only).
# All the scripts are run from this folder.
data_root   <- "data"     # one folder per library, with the Cell Ranger filtered_feature_bc_matrix
kanton_rds  <- file.path(data_root, "timecourse_human_singlecells_GRCh38.rds") # Kanton et al. 2019 timecourse
results_dir <- "results"  # all outputs are written here

pp_path    <- file.path(results_dir, "preprocessing")
joint_path <- file.path(results_dir, "shared_embedding")
ct_path    <- file.path(results_dir, "celltypist")
tsne_path  <- file.path(results_dir, "tSNE")
pt_path    <- file.path(results_dir, "monocle3")
scafe_path <- file.path(results_dir, "scafe")
tcre_path  <- file.path(results_dir, "tCRE_embedding")

samples <- read.csv("samples.csv", stringsAsFactors = FALSE)
samples$data_dir <- file.path(data_root, samples$label, "filtered_feature_bc_matrix")
