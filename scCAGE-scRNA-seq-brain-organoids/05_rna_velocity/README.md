# 5. RNA velocity with velocyto and scVelo

Spliced and unspliced counts for every library, a dynamical scVelo model per sample (T0-T4), and the
arrows drawn on the same maps as the cell types: the shared UMAP_A of [step 1](../01_preprocessing) and
the t-SNE of [step 2](../02_celltypist). Labels and colours are the CellTypist ones of step 2.

## Pipeline

```
data/<LABEL>/possorted_genome_bam.bam
  └─ 01_run_velocyto.sh          spliced / unspliced counts, one .loom per library
       └─ 02_export_cells.R      QC-passed cells, UMAP_A and labels of each timepoint
            └─ 03_scvelo_velocity.py   one dynamical model per sample, streams on UMAP_A
                 ├─ 04_velocity_tSNE.py          the same arrows on the t-SNE
                 ├─ 05_velocity_pseudotime.py    velocity pseudotime (scVelo)
                 │    └─ 06_velociraptor_fig10.14.R   OSCA Figure 10.14 with velociraptor
                 └─ 07_latent_time.py            latent time of the dynamical model
```

## Choices

| Choice | What | Why |
|---|---|---|
| Cells | only the cells kept by the QC of step 1, with their UMAP_A coordinates | the .loom files contain every barcode called by Cell Ranger; the arrows must be on the same cells as the other maps |
| Genes removed | mitochondrial, ribosomal and OXPHOS genes (regex of the pre-processing) and the GO:0007049 "cell cycle" genes with their child terms (`GO0007049_cell_cycle_genes.txt`, from org.Hs.eg.db) | stem cells and progenitors cycle a lot, and cell cycle genes would make the arrows follow the cycle instead of the differentiation |
| Filtering | `filter_and_normalize(min_shared_counts=20)`, log1p | scVelo tutorial |
| Genes in the model | top 110 up-regulated genes of each cell type (Wilcoxon, type vs rest, logFC > 0, adjusted p < 0.05, types with ≥ 20 cells), instead of the 2,000 most variable genes | the top 110 per group is the approach of La Manno et al. 2018; it keeps genes that separate the identities of our maps |
| Model | PCA 30, 30 neighbours, moments, `recover_dynamics`, dynamical velocity | scVelo Dynamical Modeling tutorial |
| Velocity graph | `sqrt_transform=True` | variance-stabilizing transform of the velocities before the cosine similarities |
| One model per sample | CTT and PAT are fitted separately, then drawn side by side on the same axes | n = 1 library per genotype × timepoint: a joint model would mix library and genotype |
| Figure 10.14 | `embedVelocity` + `gridVectors` of velociraptor on our scVelo model, colour = velocity pseudotime, scale limited to the 5th-95th percentile | the OSCA.advanced figure; the percentile limits are the `perc=[5, 95]` of the scVelo tutorial |
| Root and end (T0-T2) | root = stem cells, end = NSC in `scv.tl.velocity_pseudotime`; T3-T4 with the root estimated by scVelo | stem cells are still there at T0-T2; at T3-T4 they are almost gone |

## Reading the arrows

- The genes in the model depend on the CellTypist labels, so the arrows are read as movement between these labels.
- Velocity pseudotime and latent time go from 0 to 1 inside each sample: values are not compared between samples.
- In a t-SNE the local neighbours are kept, but the distances between far islands are not meaningful.
- With one library per genotype × timepoint, differences between CTT and PAT are descriptive.
- The script also prints a net flux between pairs of early cell types (summed transition probabilities a → b vs b → a). This is an ad hoc summary, not a published metric; CellRank would be the formal way to estimate transitions between groups.

## Output

**Streams on the shared UMAP_A, CTT | PAT**

<img src="figures/T0_stream_UMAP_A.png" width="800">
<img src="figures/T1_stream_UMAP_A.png" width="800">
<img src="figures/T2_stream_UMAP_A.png" width="800">
<img src="figures/T3_stream_UMAP_A.png" width="800">
<img src="figures/T4_stream_UMAP_A.png" width="800">

**The same arrows on the t-SNE**

<img src="figures/T0_stream_tSNE.png" width="800">
<img src="figures/T1_stream_tSNE.png" width="800">
<img src="figures/T2_stream_tSNE.png" width="800">
<img src="figures/T3_stream_tSNE.png" width="800">
<img src="figures/T4_stream_tSNE.png" width="800">

**OSCA Figure 10.14 with velociraptor: velocity pseudotime and grid arrows on the t-SNE.**
T0-T2 with root = stem cells and end = NSC, T3-T4 with the root estimated by scVelo.

| | CTT | PAT |
|---|---|---|
| T0 | <img src="figures/CTT0_velocity_pseudotime_stem_root.png" width="400"> | <img src="figures/PAT0_velocity_pseudotime_stem_root.png" width="400"> |
| T1 | <img src="figures/CTT1_velocity_pseudotime_stem_root.png" width="400"> | <img src="figures/PAT1_velocity_pseudotime_stem_root.png" width="400"> |
| T2 | <img src="figures/CTT2_velocity_pseudotime_stem_root.png" width="400"> | <img src="figures/PAT2_velocity_pseudotime_stem_root.png" width="400"> |
| T3 | <img src="figures/CTT3_velocity_pseudotime.png" width="400"> | <img src="figures/PAT3_velocity_pseudotime.png" width="400"> |
| T4 | <img src="figures/CTT4_velocity_pseudotime.png" width="400"> | <img src="figures/PAT4_velocity_pseudotime.png" width="400"> |

**Latent time of the dynamical model, on UMAP_A**

<img src="figures/T0_latent_time.png" width="800">
<img src="figures/T1_latent_time.png" width="800">
<img src="figures/T2_latent_time.png" width="800">
<img src="figures/T3_latent_time.png" width="800">
<img src="figures/T4_latent_time.png" width="800">

## Run it

From the project folder, after steps 1 and 2 (including `05_tsne.R`). The BAM files and the
`genes.gtf` of the Cell Ranger reference go in `data/` (not included).

```bash
pip install -r requirements.txt
bash 05_rna_velocity/scripts/01_run_velocyto.sh
```
```r
source("05_rna_velocity/scripts/02_export_cells.R")
```
```bash
python 05_rna_velocity/scripts/03_scvelo_velocity.py 0 1 2 3 4
python 05_rna_velocity/scripts/04_velocity_tSNE.py 0 1 2 3 4
python 05_rna_velocity/scripts/05_velocity_pseudotime.py stem_root 0 1 2
python 05_rna_velocity/scripts/05_velocity_pseudotime.py 3 4
Rscript 05_rna_velocity/scripts/06_velociraptor_fig10.14.R stem_root 0 1 2
Rscript 05_rna_velocity/scripts/06_velociraptor_fig10.14.R 3 4
python 05_rna_velocity/scripts/07_latent_time.py 0 1 2 3 4
```

R with SingleCellExperiment, zellkonverter, velociraptor, scater, ggplot2.
