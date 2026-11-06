# =============================================================================
# The same velocity streams on the t-SNE, T0-T4
# Author: Gaia Cocconcelli
# Same arrows as the streams on UMAP_A of 03_scvelo_velocity.py: the models are read from
# results/velocity/scvelo/<sample>_scvelo.h5ad (same genes, neighbours, fit and
# velocity_graph), nothing is computed again. The only change is the map: the arrows are
# projected on the t-SNE of the joint PCA (02_celltypist/scripts/05_tsne.R).
# In a t-SNE the local neighbours are kept, but the distances between far islands
# are not meaningful. The .h5ad files are only read.
# Run from the project folder, with the timepoints as arguments:
#   python 05_rna_velocity/scripts/04_velocity_tSNE.py 0 1 2 3 4
# Output: results/velocity/scvelo/T<t>_stream_tSNE.png
# =============================================================================

import os, sys
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt, anndata, scvelo as scv

scv.settings.verbosity = 0
models = os.path.join("results", "velocity", "scvelo")
tsne_dir = os.path.join("results", "tSNE")
pad = 0.5

for tp in [int(x) for x in sys.argv[1:]]:
    ts = pd.read_csv(os.path.join(tsne_dir, f"joint_T{tp}", f"T{tp}_tSNE.csv")).set_index("cell")
    adatas = {}
    for s in [f"CTT{tp}", f"PAT{tp}"]:
        f = os.path.join(models, f"{s}_scvelo.h5ad")
        assert os.path.exists(f), f"{f} not found: run 03_scvelo_velocity.py {tp} first"
        a = anndata.read_h5ad(f)
        assert "velocity_graph" in a.uns, f"{s}: no velocity graph in the model"
        missing = a.obs_names.difference(ts.index)
        assert len(missing) == 0, f"{s}: {len(missing)} cells without t-SNE coordinates"
        a.obsm["X_tsne"] = ts.loc[a.obs_names, ["tsne1", "tsne2"]].values
        if "velocity_tsne" in a.obsm:
            del a.obsm["velocity_tsne"]
        adatas[s] = a
        print(f"read {s}: {a.n_obs} cells", flush=True)

    xy = np.vstack([a.obsm["X_tsne"] for a in adatas.values()])
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, (s, a) in zip(axes, adatas.items()):
        scv.pl.velocity_embedding_stream(a, basis="tsne", color="celltype", ax=ax, show=False, title=f"{s} velocity",
                                         legend_loc="right margin" if s.startswith("PAT") else "none",
                                         legend_fontsize=7, size=15, alpha=0.6)
        ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
        ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
    fig.savefig(os.path.join(models, f"T{tp}_stream_tSNE.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved T{tp}_stream_tSNE.png", flush=True)
