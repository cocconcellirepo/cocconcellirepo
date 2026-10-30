#!/usr/bin/env bash
# SCAFE on all the libraries: workflow.sc.solo on each library, aggregate of their
# pass.ctss files (common set all10), tool.sc.count per library on the common set.
# Commands and options follow the workflows of github.com/chung-lab/SCAFE:
# the count step repeats task 7 of workflow.sc.solo, with the common set as countRegion.
# Each library has its own run_outDir, because --overwrite=yes erases it.
# Input : data/<LABEL>/possorted_genome_bam.bam, data/<LABEL>/filtered_feature_bc_matrix/barcodes.tsv.gz
# Output: results/scafe/sc_solo/<LABEL>/, results/scafe/aggregate_all10/, results/scafe/count_all10/<LABEL>/
# Run from the project folder, for example with nohup: it takes several hours.
set -uo pipefail

GENOME=hg38.gencode_v32
OUT=results/scafe
LIBS="CTT0 PAT0 CTT1 PAT1 CTT2 PAT2 CTT3 PAT3 CTT4 PAT4"
mkdir -p "$OUT"

# 1. sc.solo, one run_outDir per library
for L in $LIBS; do
  echo "=== $(date '+%F %T')  START sc.solo $L ==="
  scafe.workflow.sc.solo \
    --overwrite=yes \
    --run_bam_path="data/$L/possorted_genome_bam.bam" \
    --run_cellbarcode_path="data/$L/filtered_feature_bc_matrix/barcodes.tsv.gz" \
    --genome="$GENOME" \
    --run_tag="$L" \
    --run_outDir="$OUT/sc_solo/$L" \
    --max_thread=10
  echo "=== $(date '+%F %T')  END   sc.solo $L (exit $?) ==="
done

# 2. lib_list with all the libraries, stop if a file is missing
for L in $LIBS; do
  printf "%s\t%s\t%s\n" "$L" \
    "$PWD/$OUT/sc_solo/$L/remove_strand_invader/$L/bed/$L.pass.ctss.bed.gz" \
    "$PWD/$OUT/sc_solo/$L/bam_to_ctss/$L/bed/$L.unencoded_G.collapse.ctss.bed.gz"
done > "$OUT/lib_list_all10.txt"
cut -f2,3 "$OUT/lib_list_all10.txt" | tr '\t' '\n' | while read -r f; do
  [ -s "$f" ] || { echo "MISSING $f"; exit 1; }
done || exit 1

# 3. aggregate on the pooled signal: cluster, filter (logistic regression), annotate
echo "=== $(date '+%F %T')  START aggregate all10 ==="
scafe.workflow.cm.aggregate \
  --overwrite=yes \
  --lib_list_path="$OUT/lib_list_all10.txt" \
  --genome="$GENOME" \
  --max_thread=10 \
  --run_tag=all10 \
  --run_outDir="$OUT/aggregate_all10/"
CRE="$OUT/aggregate_all10/annotate/all10/bed/all10.CRE.coord.bed.gz"
[ -s "$CRE" ] || { echo "MISSING $CRE"; exit 1; }
echo "=== $(date '+%F %T')  END   aggregate all10 ==="

# 4. count per library on the common set
for L in $LIBS; do
  echo "=== $(date '+%F %T')  START count $L ==="
  mkdir -p "$OUT/count_all10/$L"
  scafe.tool.sc.count \
    --overwrite=yes \
    --genome="$GENOME" \
    --countRegion_bed_path="$CRE" \
    --ctss_bed_path="$OUT/sc_solo/$L/bam_to_ctss/$L/bed/$L.CB.ctss.bed.gz" \
    --cellBarcode_list_path="data/$L/filtered_feature_bc_matrix/barcodes.tsv.gz" \
    --outputPrefix="$L" \
    --outDir="$OUT/count_all10/$L/" || exit 1
  echo "=== $(date '+%F %T')  END   count $L ==="
done
echo "=== $(date '+%F %T')  DONE ==="
