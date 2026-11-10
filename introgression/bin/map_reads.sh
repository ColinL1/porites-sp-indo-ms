#!/bin/bash
# Map the ONT reads to the target assembly; flagstat and mean depth over all positions.
# usage: map_reads.sh <target.fa> <reads.fastq.gz> <sample> <out.bam> <threads>
set -euo pipefail
REF=$1 READS=$2 SAMPLE=$3 BAM=$4 THREADS=$5
SORT_T=4   # piped samtools sort; minimap2 gets the rest
MAP_T=$(( THREADS > SORT_T ? THREADS - SORT_T : THREADS ))

minimap2 -t $MAP_T -ax map-ont --MD -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT" "$REF" "$READS" \
  | samtools sort -@ $SORT_T -m 2G -o "$BAM" -
samtools index -@ $THREADS "$BAM"
samtools flagstat "$BAM" > "${BAM%.bam}.flagstat.txt"
samtools depth -a -@ $THREADS "$BAM" \
  | awk '{s+=$3; n++} END{print "mean_depth_all_positions\t"s/n"\ttotal_pos\t"n}' > "${BAM%.bam}.meandepth.txt"
