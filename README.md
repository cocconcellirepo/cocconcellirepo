# Gaia Cocconcelli

Bioinformatician working on single-cell transcriptomics.

I'm doing my MSc thesis in Molecular Biotechnology and Bioinformatics at the University of Milan,
on a project in collaboration with Human Technopole (wet lab) and the University of Milano-Bicocca.
The data are single-cell 5′-end RNA-seq (scCAGE-seq, 10x Genomics 5′) of human brain organoids:
two genotypes followed over five timepoints, from iPSCs to day 20, to see how a disease-associated
mutation changes early neural development.

This profile collects the pipelines I wrote for that project, with their outputs and the reasons
behind each choice. Raw data are not included; every step explains what the input looks like.

## Project

**[scCAGE-scRNA-seq-brain-organoids](scCAGE-scRNA-seq-brain-organoids)**: the whole analysis in one place, step by step.

| Step | What it does |
|---|---|
| [1. Pre-processing](scCAGE-scRNA-seq-brain-organoids/01_preprocessing) | QC, doublet removal, normalization, cell cycle scoring and shared embeddings without batch correction (Bioconductor, OSCA workflow) |
| [2. Cell type annotation](scCAGE-scRNA-seq-brain-organoids/02_celltypist) | A CellTypist model trained on the Kanton et al. 2019 organoid timecourse, with stage-specific thresholds |
| [3. Trajectory](scCAGE-scRNA-seq-brain-organoids/03_trajectory) | Pseudotime across five timepoints with monocle3, and a numbered path on the principal graph computed from the data |

## Coming next

I'm adding the other analyses one at a time, as each one is cleaned up:

- **Transcribed cis-regulatory elements with SCAFE**: promoters and enhancer candidates from 10x 5′ data, a pooled common set for all samples, count matrices, and UMAPs built on regulatory elements only
- **RNA velocity**: velocyto and scVelo
- **Allele-specific expression**: a single locus, from single-cell BAM files

## Tools

R / Bioconductor (SingleCellExperiment, scater, scran, batchelor, scDblFinder), Seurat, monocle3,
slingshot, ggplot2 · Python (CellTypist, scanpy, scVelo) · SCAFE, samtools, pysam · bash on a Linux
server (long jobs with nohup) · git

## Background

BSc in Biological Sciences, Università Politecnica delle Marche · MSc in Molecular Biotechnology and
Bioinformatics, University of Milan · Internship at the European Institute of Oncology (IEO), Milan
