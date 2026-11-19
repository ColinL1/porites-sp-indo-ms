#!/usr/bin/env bash
# busco_stage.sh: BUSCO on an assembly stage the pipeline did not assess (e.g. the RagTag scaffolds before gap
# filling), with the command and container of GenomeMaker's BUSCO process (modules/02_multi_assembler.nf:
# quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1, genome mode, offline lineage dataset).
# usage: busco_stage.sh assembly.fasta busco_db_dir lineage name out_dir [threads]
# writes out_dir/<name>.<lineage>.short_summary.txt (and .json, full_table.tsv), as the pipeline names them
set -euo pipefail
ASM=$(readlink -f "$1"); DB=$(readlink -f "$2"); LIN=$3; NAME=$4; OUT=$5; T=${6:-32}
IMG=${BUSCO_IMG:-$HOME/SINGULARITY_CACHE/quay.io-biocontainers-busco-6.1.0--pyhdfd78af_1.img}
mkdir -p "$OUT"; OUT=$(readlink -f "$OUT"); W="$OUT/busco_work_$NAME"; mkdir -p "$W"
singularity exec --containall --home "$W:$W" --pwd "$W" -B "$W" -B "$(dirname "$ASM")" -B "$DB" "$IMG" \
  busco -i "$ASM" -m geno -l "$LIN" -c "$T" -o busco_output --offline --download_path "$DB" -f
cp "$W/busco_output/run_$LIN/full_table.tsv" "$OUT/$NAME.$LIN.full_table.tsv"
cp "$W/busco_output/short_summary.specific.$LIN.busco_output.txt" "$OUT/$NAME.$LIN.short_summary.txt"
cp "$W/busco_output/short_summary.specific.$LIN.busco_output.json" "$OUT/$NAME.$LIN.short_summary.json"
grep -m1 'C:' "$OUT/$NAME.$LIN.short_summary.txt"
