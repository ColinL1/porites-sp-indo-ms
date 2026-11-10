#!/bin/bash
# Fallback when no --new_mito is given: assemble the target mitogenome from the ONT reads that
# map to the bait (subsampled to 15%, ~200x), Flye --nano-hq. Writes <outdir>/flye/assembly.fasta.
# usage: mito_from_reads.sh <bait.fasta> <reads.fastq.gz> <outdir> <threads> <flye command>
set -euo pipefail
BAIT=$1 READS=$2 OUT=$3 THREADS=$4; shift 4
FLYE=("$@")

minimap2 -ax map-ont -t $THREADS "$BAIT" "$READS" 2>/dev/null \
  | samtools view -b -F 0x904 - | samtools sort -o "$OUT/mito_reads.bam"
samtools index "$OUT/mito_reads.bam"
samtools fastq "$OUT/mito_reads.bam" | gzip > "$OUT/mito_reads.fastq.gz"
echo "mito-mapped reads: $(zcat "$OUT/mito_reads.fastq.gz" | awk 'END{print NR/4}')"
seqkit sample -s11 -p 0.15 "$OUT/mito_reads.fastq.gz" 2>/dev/null | gzip > "$OUT/mito_sub.fastq.gz"
"${FLYE[@]}" --nano-hq "$OUT/mito_sub.fastq.gz" --out-dir "$OUT/flye" --threads $THREADS --min-overlap 1000
