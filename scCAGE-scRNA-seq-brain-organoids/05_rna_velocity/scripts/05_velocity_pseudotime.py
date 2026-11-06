# =============================================================================
# Velocity pseudotime (step 1 of OSCA.advanced Figure 10.14), T0-T4
# Author: Gaia Cocconcelli
# velociraptor::scvelo() gets the pseudotime from scv.tl.velocity_pseudotime. Here the
# same scVelo function is run on our models (03_scvelo_velocity.py) instead of a new
# velociraptor fit.
# - T0-T4, default parameters: root and end cells estimated by scVelo.
#   Output: <sample>_velocity_pseudotime.csv
# - T0-T2 only, with "stem_root" as first argument: root = stem cells, end = NSC.
#   root_key = "root_stem" (1 for stem cells 1/2/3, 0 for the others)
#   end_key  = "end_nsc"   (1 for NSC, 0 for the others)
#   In scVelo 0.3 a 0/1 column is turned into one cell: the cell with the most root
#   (end) cells among its neighbours. The distance is a diffusion distance on the
#   velocity graph made symmetric (T + T.T), so the colours follow mostly this choice
#   of root and end point, while the arrows stay those of the velocity. T3-T4 are not
#   included: there are almost no stem cells there.
#   Output: <sample>_velocity_pseudotime_stem_root.csv
# One sample at a time: the pseudotime goes from 0 to 1 within each sample.
# The .h5ad files are only read, nothing is saved in them.
# Run from the project folder:
#   python 05_rna_velocity/scripts/05_velocity_pseudotime.py 0 1 2 3 4
#   python 05_rna_velocity/scripts/05_velocity_pseudotime.py stem_root 0 1 2
# Output folder: results/velocity/velociraptor/
# =============================================================================

import os, sys
import numpy as np, pandas as pd, anndata, scvelo as scv
from scvelo.preprocessing.moments import get_connectivities

scv.settings.verbosity = 1
models = os.path.join("results", "velocity", "scvelo")
out = os.path.join("results", "velocity", "velociraptor")
os.makedirs(out, exist_ok=True)
stem = ["stem cells 1", "stem cells 2", "stem cells 3"]

args = sys.argv[1:]
stem_root = len(args) > 0 and args[0] == "stem_root"
tps = [int(x) for x in (args[1:] if stem_root else args)]


def chosen_cell(a, key):
    # the same choice as VPT.set_iroot in scVelo: neighbours x 0/1 column, highest value
    w = get_connectivities(a).dot(a.obs[key].values)
    return a.obs["celltype"].astype(str).values[np.argmax(w)]


for tp in tps:
    assert not stem_root or tp <= 2, "root = stem cells, only for T0-T2"
    for s in [f"CTT{tp}", f"PAT{tp}"]:
        f = os.path.join(models, f"{s}_scvelo.h5ad")
        assert os.path.exists(f), f"{f} not found"
        a = anndata.read_h5ad(f)
        ct = a.obs["celltype"].astype(str)
        if stem_root:
            a.obs["root_stem"] = ct.isin(stem).astype(float).values
            a.obs["end_nsc"] = (ct == "NSC").astype(float).values
            assert a.obs["root_stem"].sum() > 0 and a.obs["end_nsc"].sum() > 0, f"{s}: no stem cells or no NSC"
            scv.tl.velocity_pseudotime(a, root_key="root_stem", end_key="end_nsc")
            print(f"{s}: root cell = {chosen_cell(a, 'root_stem')}, end cell = {chosen_cell(a, 'end_nsc')}", flush=True)
        else:
            scv.tl.velocity_pseudotime(a)
        pt = a.obs["velocity_pseudotime"].values
        suffix = "_stem_root" if stem_root else ""
        pd.DataFrame({"cell": a.obs_names, "celltype": ct.values, "velocity_pseudotime": pt}).to_csv(
            os.path.join(out, f"{s}_velocity_pseudotime{suffix}.csv"), index=False)
        print(f"{s}: {len(pt)} cells, pseudotime NA {int(np.isnan(pt).sum())}, "
              f"< 0.1 {100 * np.nanmean(pt < 0.1):.0f}%, > 0.9 {100 * np.nanmean(pt > 0.9):.0f}%", flush=True)
