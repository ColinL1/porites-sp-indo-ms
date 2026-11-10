#!/bin/bash
# Align one reference assembly to the target and project it onto target coordinates with
# paftools.js call: callable regions (BED) and variants (VCF, haploid calls written as 1/1).
# usage: align_reference.sh <name> <ref.fna> <preset asm10|asm20> <target.fa> <outdir> <workdir> <threads>
set -euo pipefail
NAME=$1 FNA=$2 PRESET=$3 TARGET=$4 OUT=$5 WORK=$6 THREADS=$7

minimap2 -cx "$PRESET" -t $THREADS --cs "$TARGET" "$FNA" 2> "$WORK/$NAME.mm2.log" \
  | sort -k6,6 -k8,8n > "$WORK/$NAME.paf"
paftools.js call -L 10000 -l 5000 "$WORK/$NAME.paf" > "$WORK/$NAME.calls.txt" 2> "$WORK/$NAME.paftools.log"
awk '$1=="R"{print $2"\t"$3"\t"$4}' "$WORK/$NAME.calls.txt" | sort -k1,1 -k2,2n > "$OUT/$NAME.callable.bed"
paftools.js call -s "$NAME" -L 10000 -l 5000 -f "$TARGET" "$WORK/$NAME.paf" 2>> "$WORK/$NAME.paftools.log" \
  | bcftools sort -T "$WORK" -Oz -o "$OUT/$NAME.vcf.gz" -
tabix -f -p vcf "$OUT/$NAME.vcf.gz"
printf "%s\t%s\t%s\t%s\n" "$NAME" "$PRESET" "$(bcftools view -H "$OUT/$NAME.vcf.gz" | wc -l)" \
  "$(awk '{s+=$3-$2} END{print s+0}' "$OUT/$NAME.callable.bed")" > "$WORK/$NAME.stats.tsv"
