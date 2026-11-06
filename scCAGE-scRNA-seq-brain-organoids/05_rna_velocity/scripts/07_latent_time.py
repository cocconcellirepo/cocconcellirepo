# =============================================================================
# Latent time as in the scVelo tutorial "Dynamical Modeling", T0-T4
# Author: Gaia Cocconcelli
# Tutorial steps:
#   scv.tl.latent_time(adata)
#   scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", size=80)
# on our scVelo models (03_scvelo_velocity.py), default parameters: root cells
# estimated by scVelo, genes with a dynamical fit and fit_likelihood >= 0.1.
# Map: the default basis of scv.pl.scatter, i.e. obsm["X_umap"], which in our models is
# the shared UMAP_A of the timepoint. No arrows.
# Difference from the tutorial: point size left to the scVelo default (set from the
# number of cells), because size = 80 is for the 3,696 cells of the tutorial dataset and
# our samples have 7,000-14,000 cells.
# One sample at a time: the latent time goes from 0 to 1 within each sample.
# The .h5ad files are only read, nothing is saved in them.
# Run from the project folder, with the timepoints as arguments:
#   python 05_rna_velocity/scripts/07_latent_time.py 0 1 2 3 4
# Output (results/velocity/latent_time/): T<t>_latent_time.png (CTT | PAT), <sample>_latent_time.csv
# =============================================================================

import os, sys
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt, anndata, scvelo as scv

scv.settings.verbosity = 1
models = os.path.join("results", "velocity", "scvelo")
out = os.path.join("results", "velocity", "latent_time")
os.makedirs(out, exist_ok=True)

for tp in [int(x) for x in sys.argv[1:]]:
    adatas = {}
    for s in [f"CTT{tp}", f"PAT{tp}"]:
        f = os.path.join(models, f"{s}_scvelo.h5ad")
        assert os.path.exists(f), f"{f} not found"
        a = anndata.read_h5ad(f)
        scv.tl.latent_time(a)
        lt = a.obs["latent_time"].values
        pd.DataFrame({"cell": a.obs_names, "celltype": a.obs["celltype"].astype(str).values,
                      "latent_time": lt}).to_csv(os.path.join(out, f"{s}_latent_time.csv"), index=False)
        print(f"{s}: latent time NA {int(np.isnan(lt).sum())}, < 0.1 {100 * np.nanmean(lt < 0.1):.0f}%, "
              f"> 0.9 {100 * np.nanmean(lt > 0.9):.0f}%, median {np.nanmedian(lt):.2f}", flush=True)
        adatas[s] = a

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, (s, a) in zip(axes, adatas.items()):
        scv.pl.scatter(a, color="latent_time", color_map="gnuplot", ax=ax, show=False, title=f"{s} latent time")
    fig.savefig(os.path.join(out, f"T{tp}_latent_time.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"saved T{tp}_latent_time.png", flush=True)
