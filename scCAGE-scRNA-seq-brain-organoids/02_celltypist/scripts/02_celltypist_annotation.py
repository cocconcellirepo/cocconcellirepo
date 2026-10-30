# =============================================================================
# Cell types with CellTypist, trained on the Kanton timecourse (label = Cluster)
# Author: Gaia Cocconcelli
# Input: the files written by 01_export_for_celltypist.R in results/celltypist/
# One model for all timepoints (T0-T4): Kanton cell types, all stages,
# h409B2 + H9, all genes shared with our data. CellTypist is free to choose any
# of the types, also the organoid ones.
# One label per cell, no majority voting (it would relabel the cells by cluster).
# Output: model_Kanton_Cluster.pkl, T<t>_CellTypist_Cluster_labels.csv,
# T<t>_CellTypist_Cluster_probabilities.csv
# =============================================================================

## ---- Step 0: setup ----
import os
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
import scanpy as sc
import celltypist

ct_path = os.path.join("results", "celltypist")
n_jobs = 8

def read_export(folder):
    # counts.mtx is genes x cells (as written by R), AnnData wants cells x genes
    X = scipy.io.mmread(os.path.join(folder, "counts.mtx")).T.tocsr().astype(np.float32)
    genes = open(os.path.join(folder, "genes.txt")).read().splitlines()
    cells = pd.read_csv(os.path.join(folder, "cells.csv"))
    cells.index = cells.iloc[:, 0].astype(str)
    adata = ad.AnnData(X=X, obs=cells, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()
    # CellTypist input: counts scaled to 10,000 per cell, then log1p, on ALL genes
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

## ---- Step 1: reference ----
ref = read_export(os.path.join(ct_path, "ref_kanton_cluster"))
print(ref.obs["label"].value_counts())

## ---- Step 2: train the model ----
# genes_common.txt = genes present in the reference and in all five joint objects.
# The genes are subset after the normalization, so the cells no longer sum to
# 10,000: check_expression=False skips that check (the normalization is done above).
# feature_selection=True: a first fit on all genes, then a second fit on the
# top 300 genes of each cell type (CellTypist default).
genes = open(os.path.join(ct_path, "genes_common.txt")).read().splitlines()
ref = ref[:, genes].copy()
print(f"training on {ref.n_obs} cells x {ref.n_vars} genes")
model = celltypist.train(ref, labels="label", check_expression=False,
                         feature_selection=True, n_jobs=n_jobs)
model.write(os.path.join(ct_path, "model_Kanton_Cluster.pkl"))
print(f"{len(model.cell_types)} cell types, {len(model.features)} features")

## ---- Step 3: annotate CTT and PAT at each timepoint ----
# "best match": the cell type with the highest probability.
# "prob match", p_thres 0.5: "Unassigned" if no type reaches 0.5, two types
# joined by "|" if both pass.
# The probabilities are one per cell type (one-vs-rest), they do not sum to 1.
for tp in range(5):
    q = read_export(os.path.join(ct_path, f"joint_T{tp}"))
    best = celltypist.annotate(q, model=model, mode="best match")
    prob = celltypist.annotate(q, model=model, mode="prob match", p_thres=0.5)
    out = q.obs[["cell_index", "barcode", "sample"]].copy()
    out["CellTypist"] = best.predicted_labels["predicted_labels"].values
    out["CellTypist_prob"] = best.probability_matrix.max(axis=1).values
    out["CellTypist_probmatch"] = prob.predicted_labels["predicted_labels"].values
    out.to_csv(os.path.join(ct_path, f"T{tp}_CellTypist_Cluster_labels.csv"), index=False)
    best.probability_matrix.to_csv(os.path.join(ct_path, f"T{tp}_CellTypist_Cluster_probabilities.csv"))
    print("")
    print(f"T{tp}: best match, cells per sample")
    print(pd.crosstab(out["CellTypist"], out["sample"]))
    print(f"T{tp}: cells with no type at probability >= 0.5")
    print(out.groupby("sample")["CellTypist_probmatch"].apply(lambda x: (x == "Unassigned").sum()))
