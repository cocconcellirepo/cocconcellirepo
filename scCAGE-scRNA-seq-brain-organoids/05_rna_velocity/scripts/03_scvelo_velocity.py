# =============================================================================
# RNA velocity with scVelo (dynamical model), one sample at a time, T0-T4
# Author: Gaia Cocconcelli
# Steps of the scVelo tutorials (RNA Velocity Basics, Dynamical Modeling), with these changes:
# - only our QC-passed cells, with the shared UMAP_A of the timepoint and the CellTypist
#   labels (02_export_cells.R);
# - genes removed from the loom before the filtering: mitochondrial, ribosomal and
#   OXPHOS genes (same regex as the pre-processing) and the cell cycle genes of
#   GO:0007049 (GO0007049_cell_cycle_genes.txt), so that the arrows are not driven
#   by the cell cycle;
# - instead of the 2,000 most variable genes, the top 110 up-regulated genes of each
#   cell type (Wilcoxon, each type against the rest, logFC > 0, adjusted p < 0.05;
#   types with fewer than 20 cells give no genes, their cells stay). The top 110 per
#   group is the idea of La Manno et al. 2018 (intestine); the test and the labels are
#   ours. The genes depend on the labels, so the arrows are read as movement between
#   these labels;
# - velocity_graph with sqrt_transform = True (variance-stabilizing transform of the
#   velocities before the cosine similarities).
# The script also prints a net flux between pairs of early cell types: the transition
# probabilities summed from type a to type b and back, as a ratio (ad hoc summary,
# not a published metric; CellRank would be the formal alternative).
# Run from the project folder, with the timepoints as arguments:
#   python 05_rna_velocity/scripts/03_scvelo_velocity.py 0 1 2 3 4
# Output (results/velocity/scvelo/): T<t>_stream_UMAP_A.png, <sample>_genes_per_type.csv,
#   <sample>_scvelo.h5ad (read by the next scripts), summary on screen
# =============================================================================

import os, sys
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt, anndata, scanpy as sc, scvelo as scv

scv.settings.verbosity = 1
vel = os.path.join("results", "velocity")
cells_path = os.path.join(vel, "cells")
out = os.path.join(vel, "scvelo")
os.makedirs(out, exist_ok=True)
n_jobs = 8
N_TOP = 110
MIN_CELLS_GROUP = 20
drop_regex = r"^MT-|^RP[SL]|^ATP5|^COX[0-9]|^NDUF|^UQCR|^SDHA|^SDHB|^SDHC|^SDHD|^ROMO1$"
go = set(open("05_rna_velocity/GO0007049_cell_cycle_genes.txt").read().split())
short = {"neuroectodermal-like cells": "NE", "stem cells 1": "stem1", "stem cells 2": "stem2", "NSC": "NSC"}
pairs = [("neuroectodermal-like cells", "stem cells 1"), ("neuroectodermal-like cells", "stem cells 2"),
         ("neuroectodermal-like cells", "NSC")]
pd.set_option("display.width", 250)

# same colours as 03_celltypist_plots.R
pal = {"stem cells 1": "#1f77b4", "stem cells 2": "#17becf", "stem cells 3": "#9467bd",
       "neuroectodermal-like cells": "#2ca02c", "NSC": "#d62728", "organoid-only types": "#ff7f0e",
       "NSC/radial glia": "#F3C300", "G2M/S NPCs": "#654522", "G2M/S dorsal progenitors 1": "#C51B7D",
       "mesenchymal-like cells": "#8DB600", "midbrain/hindbrain": "#E68FAC", "radial glia 2": "#0000CD",
       "ventral progenitors and neurons 1": "#4B0082", "ventral progenitors and neurons 4": "#C2B280",
       "G2M/S ventral progenitors and IPs": "#F99379", "choroid plexus": "#FF00FF",
       "cortical neurons 1": "#00FF00", "cortical neurons 2": "#00FFCC", "G2M/S dorsal progenitors 2": "#0066FF",
       "gliogenic/outer radial glia": "#336633", "IPs and early cortical neurons": "#CC66FF",
       "radial glia 1": "#110022", "ventral progenitors and neurons 2": "#993333",
       "ventral progenitors and neurons 3": "#FF0066", "Unassigned": "#cccccc"}


def net_flux(ad):
    # total transition probability from type a to type b, and back
    T = scv.utils.get_transition_matrix(ad).tocsr()
    ct = ad.obs["celltype"].astype(str).values
    res = {}
    for a, b in pairs:
        ia, ib = np.where(ct == a)[0], np.where(ct == b)[0]
        ab, ba = T[ia][:, ib].sum(), T[ib][:, ia].sum()
        key = f"{short[a]}<->{short[b]}"
        if ab + ba == 0:
            res[key] = "n.a."
        elif ab >= ba:
            res[key] = f"{short[a]}->{short[b]} {ab / max(ba, 1e-12):.2f}"
        else:
            res[key] = f"{short[b]}->{short[a]} {ba / max(ab, 1e-12):.2f}"
    return res


def run_sample(s, cells):
    ad = anndata.read_loom(os.path.join(vel, s, f"{s}.loom"))
    ad.var_names_make_unique()
    # loom barcodes "CTT0:AAAC...x" -> joint barcodes "CTT0_AAAC...-1"
    ad.obs_names = [n.replace(":", "_")[:-1] + "-1" for n in ad.obs_names]
    missing = set(cells["barcode"]) - set(ad.obs_names)
    assert len(missing) == 0, "some of our cells are not in the loom file"
    ad = ad[cells["barcode"].values].copy()
    ad.obsm["X_umap"] = cells[["u1", "u2"]].to_numpy()
    cats = [k for k in pal if k in set(cells["celltype"])]
    ad.obs["celltype"] = pd.Categorical(cells["celltype"].values, categories=cats)
    ad.uns["celltype_colors"] = [pal[k] for k in cats]

    drop = ad.var_names.str.contains(drop_regex)
    in_go = ad.var_names.isin(list(go))
    ad = ad[:, ~(drop | in_go)].copy()
    scv.pp.filter_and_normalize(ad, min_shared_counts=20)
    sc.pp.log1p(ad)
    n_shared = ad.n_vars

    # the top 110 up-regulated genes of each cell type, instead of the HVGs
    counts = ad.obs["celltype"].value_counts()
    groups = [g for g in counts.index if counts[g] >= MIN_CELLS_GROUP]
    sc.tl.rank_genes_groups(ad, "celltype", groups=groups, reference="rest", method="wilcoxon")
    sel = {}
    for g in groups:
        r = sc.get.rank_genes_groups_df(ad, group=g)
        r = r[(r["logfoldchanges"] > 0) & (r["pvals_adj"] < 0.05)]
        sel[g] = list(r["names"][:N_TOP])
    pd.DataFrame([(g, x) for g, v in sel.items() for x in v], columns=["celltype", "gene"]).to_csv(
        os.path.join(out, f"{s}_genes_per_type.csv"), index=False)
    genes = sorted(set().union(*sel.values()))
    print(f"\n===== {s}: {ad.n_obs} cells, genes after min_shared_counts {n_shared}; genes per cell type "
          + ", ".join(f"{g} {len(v)}" for g, v in sel.items()) + f"; union {len(genes)}", flush=True)

    ad = ad[:, genes].copy()
    sc.pp.pca(ad, n_comps=30)
    sc.pp.neighbors(ad, n_pcs=30, n_neighbors=30)
    scv.pp.moments(ad, n_pcs=None, n_neighbors=None)
    scv.tl.recover_dynamics(ad, n_jobs=n_jobs)
    scv.tl.velocity(ad, mode="dynamical")
    scv.tl.velocity_graph(ad, sqrt_transform=True, n_jobs=n_jobs)
    del ad.uns["rank_genes_groups"]   # kept in <sample>_genes_per_type.csv
    ad.write(os.path.join(out, f"{s}_scvelo.h5ad"))
    row = {"sample": s, "genes": ad.n_vars, "n_velocity_genes": int(ad.var["velocity_genes"].sum()), **net_flux(ad)}
    return ad, row


for tp in [int(x) for x in sys.argv[1:]]:
    cells_tp = pd.read_csv(os.path.join(cells_path, f"T{tp}_cells_for_scvelo.csv"))
    adatas, rows = {}, []
    for s in [f"CTT{tp}", f"PAT{tp}"]:
        adatas[s], row = run_sample(s, cells_tp[cells_tp["sample"] == s])
        rows.append(row)
    print("\n" + pd.DataFrame(rows).to_string(index=False), flush=True)

    # CTT and PAT side by side on the same UMAP_A axes
    xy = np.vstack([a.obsm["X_umap"] for a in adatas.values()])
    pad = 0.5
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, (s, a) in zip(axes, adatas.items()):
        scv.pl.velocity_embedding_stream(a, basis="umap", color="celltype", ax=ax, show=False, title=f"{s} velocity",
                                         legend_loc="right margin" if s.startswith("PAT") else "none",
                                         legend_fontsize=7, size=15, alpha=0.6)
        ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
        ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
    fig.savefig(os.path.join(out, f"T{tp}_stream_UMAP_A.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved T{tp}_stream_UMAP_A.png", flush=True)
