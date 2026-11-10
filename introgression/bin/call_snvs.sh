#!/bin/bash
# Clair3 on the scaffold set (the SNV filter runs separately: the clair3 env has no bcftools).
# usage: call_snvs.sh <target.fa> <reads.bam> <scaffolds.bed> <sample> <model_dir> <clair3_outdir> <threads>
set -euo pipefail
REF=$1 BAM=$2 BED=$3 SAMPLE=$4 MODEL=$5 OUT=$6 THREADS=$7
[ -d "$MODEL" ] || { echo "Clair3 model missing: $MODEL" >&2; exit 1; }

run_clair3.sh --bam_fn="$BAM" --ref_fn="$REF" --threads=$THREADS --platform=ont \
  --model_path="$MODEL" --output="$OUT" --bed_fn="$BED" --min_coverage=4 --sample_name="$SAMPLE"
[ -s "$OUT/pileup.vcf.gz" ] || { echo "Clair3 wrote no pileup.vcf.gz" >&2; exit 1; }
