#!/usr/bin/env bash
# Spliced / unspliced counts with velocyto, one library at a time (all ten libraries).
# Only the cells called by Cell Ranger are counted (-b filtered barcodes); the QC of
# step 1 is applied later, in 02_export_cells.R.
# Input : data/<LABEL>/possorted_genome_bam.bam, data/<LABEL>/filtered_feature_bc_matrix/barcodes.tsv.gz,
#         data/genes.gtf (the genes.gtf of the Cell Ranger reference used for the counts)
# Output: results/velocity/<LABEL>/<LABEL>.loom, results/velocity/<LABEL>/velocyto.log
# velocyto writes a cell-sorted copy of the BAM (cellsorted_possorted_genome_bam.bam)
# next to the original one, so data/<LABEL>/ has to be writable during the run.
# Run from the project folder, for example with nohup: it takes several hours per library.
set -uo pipefail

OUT=results/velocity
LIBS="CTT0 PAT0 CTT1 PAT1 CTT2 PAT2 CTT3 PAT3 CTT4 PAT4"

for L in $LIBS; do
  mkdir -p "$OUT/$L"
  zcat "data/$L/filtered_feature_bc_matrix/barcodes.tsv.gz" > "$OUT/$L/barcodes.tsv"
  velocyto run -vv -b "$OUT/$L/barcodes.tsv" -o "$OUT/$L" -e "$L" -@ 8 --samtools-memory 3000 \
    "data/$L/possorted_genome_bam.bam" data/genes.gtf > "$OUT/$L/velocyto.log" 2>&1
  echo "$L exit $? $(date '+%F %T')"
done
echo "ALL DONE $(date '+%F %T')"
