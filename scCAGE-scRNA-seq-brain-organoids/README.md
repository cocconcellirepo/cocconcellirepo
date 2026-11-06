# scCAGE / scRNA-seq of human brain organoids

Single-cell 5′-end RNA-seq (scCAGE-seq, 10x Genomics 5′) of human brain organoids. Two genotypes, a
patient-derived line (PAT) and its isogenic corrected control (CTT), at five timepoints: T0 (iPSCs),
T1, T2, T3 (day 13) and T4 (day 20). Ten libraries, one per genotype × timepoint.

Each step is a folder with its own README: what it does, why each choice was made, and the outputs.

| Step | What it does | Main tools |
|---|---|---|
| [1. Pre-processing](01_preprocessing) | QC, doublet removal, normalization, cell cycle scores; one shared map per timepoint, without batch correction | scater, scran, scDblFinder, batchelor |
| [2. Cell type annotation](02_celltypist) | a CellTypist model trained on the Kanton et al. 2019 organoid timecourse; labels, composition, marker check, t-SNE, marker dot plots | CellTypist, scanpy, scran, ggplot2 |
| [3. Trajectory](03_trajectory) | pseudotime across the five timepoints, one map per genotype, with a numbered path computed from the principal graph; slingshot curves inside T0-T2 | monocle3, slingshot, igraph, ggrepel |
| [4. Regulatory elements](04_tCRE_SCAFE) | transcribed cis-regulatory elements from the 5′ signal, a common set for all samples, and maps built on regulatory elements only | SCAFE, scran, batchelor |
| [5. RNA velocity](05_rna_velocity) | spliced / unspliced counts, one dynamical model per sample, arrows on UMAP_A and t-SNE, velocity pseudotime and latent time | velocyto, scVelo, velociraptor |

Still to add: allele-specific expression.

## Reading the results

There is one library per genotype × timepoint, so genotype and timepoint cannot be separated from
library effects: differences between samples are read as descriptive, not tested.

## Layout

```
config.R        paths and sample table, used by all the R scripts
samples.csv     labels and expected doublet rate of each library
requirements.txt  Python packages for steps 2 and 5
01_preprocessing/  02_celltypist/  03_trajectory/  04_tCRE_SCAFE/  05_rna_velocity/    README, scripts/, figures/
```

Raw data and large objects are not included. All the scripts are run from this folder.
