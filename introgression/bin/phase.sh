#!/bin/bash
# WhatsHap phasing per scaffold (in parallel), concatenation, block stats, haplotagging.
# usage: phase.sh <target.fa> <reads.bam> <snps.vcf.gz> <scaffolds.bed> <outdir> <haplotagged.bam> <threads>
set -euo pipefail
REF=$1 BAM=$2 SNPS=$3 BED=$4 OUT=$5 HTBAM=$6 THREADS=$7
mkdir -p "$OUT/per_scaffold"

phase_one() {   # one scaffold; whatshap phase is single-threaded
  local c=$1 d=$OUT/per_scaffold
  bcftools view -r "$c" "$SNPS" -Oz -o "$d/$c.in.vcf.gz" 2> "$d/$c.log"
  tabix -f -p vcf "$d/$c.in.vcf.gz"
  if [ "$(bcftools view -H "$d/$c.in.vcf.gz" | wc -l)" -lt 2 ]; then
    cp "$d/$c.in.vcf.gz" "$d/$c.phased.vcf.gz"
  else
    whatshap phase --ignore-read-groups --chromosome "$c" -o "$d/$c.phased.vcf.gz" \
      --reference "$REF" "$d/$c.in.vcf.gz" "$BAM" >> "$d/$c.log" 2>&1
  fi
  tabix -f -p vcf "$d/$c.phased.vcf.gz"
}
export -f phase_one; export REF BAM SNPS OUT
cut -f1 "$BED" | xargs -P $THREADS -I{} bash -c 'phase_one "$1"' _ {}

bcftools concat -Oz -o "$OUT/phased.vcf.gz" $(cut -f1 "$BED" | sed "s#^#$OUT/per_scaffold/#; s#\$#.phased.vcf.gz#")
tabix -f -p vcf "$OUT/phased.vcf.gz"
whatshap stats "$OUT/phased.vcf.gz" --tsv="$OUT/phase_block_stats.tsv" --block-list="$OUT/phase_block_list.tsv" \
  --chr-lengths <(cut -f1,2 "$REF.fai") > "$OUT/phase_stats.txt" 2>&1
whatshap haplotag --ignore-read-groups --skip-missing-contigs -o "$HTBAM" --reference "$REF" \
  --output-haplotag-list "${HTBAM%.bam}.haplotag_list.tsv.gz" "$OUT/phased.vcf.gz" "$BAM" > "${HTBAM%.bam}.haplotag.log" 2>&1
samtools index -@ $THREADS "$HTBAM"
