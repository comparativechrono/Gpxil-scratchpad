#!/usr/bin/env bash
# HG002's GIAB benchmarks as gpxil files, one per autosome, without the reference plane (gt and flags only), on the
# Gpxil browser's GRCh38 chromosomes:
#
#   ./encode_hg002.sh DATA REF_DIR OUT
#
# DATA holds giab/ (v4.2.1: HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz and _noinconsistent.bed) and giab5/ (v5.0q:
# HG002_GRCh38_v5.0q_smvar.vcf.gz and .benchmark.bed), from GIAB's releases; writes OUT/v4.2.1/ and OUT/v5.0q/.
# About 30 s a chromosome; 30 MB for v4.2.1's 22 files, 38 MB for v5.0q's.
set -u
DATA=$1 REF=$2 OUT=$3
encode() {  # name vcf bed
  mkdir -p "$OUT/$1"
  for i in $(seq 1 22); do
    out="$OUT/$1/HG002_chr$i.gpxil.tif"
    [ -s "$out" ] && continue
    gpxil encode "$REF/chr$i.gpxil.tif" -v "$2" --callable "$3" --contigs "chr$i" --no-ref --planes gt,flags -o "$out" \
      > "$OUT/$1/log_chr$i.txt" 2>&1 || echo "chr$i failed: see $OUT/$1/log_chr$i.txt"
  done
}
encode v4.2.1 "$DATA/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" "$DATA/giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
encode v5.0q "$DATA/giab5/HG002_GRCh38_v5.0q_smvar.vcf.gz" "$DATA/giab5/HG002_GRCh38_v5.0q_smvar.benchmark.bed"
