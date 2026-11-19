#!/usr/bin/env bash
# omark_release.sh: OMArk on the release proteins, with the commands and container of GenomeMaker's OMARK
# process (modules/07b_annotation_qc.nf; quay.io/biocontainers/omark:0.5.0--pyhdfd78af_0, OMAmer Metazoa.h5).
# The pipeline ran it at stage 07, before release curation; this repeats it on the release files.
# usage: omark_release.sh proteins.fa annotation.gff3 Metazoa.h5 out_dir [threads]
set -euo pipefail
PROT=$(readlink -f "$1"); GFF=$(readlink -f "$2"); DB=$(readlink -f "$3"); OUT=$4; T=${5:-32}
IMG=${OMARK_IMG:-$HOME/SINGULARITY_CACHE/quay.io-biocontainers-omark-0.5.0--pyhdfd78af_0.img}
mkdir -p "$OUT"; OUT=$(readlink -f "$OUT"); cd "$OUT"
run() { singularity exec --containall --home "$OUT:$OUT" --pwd "$OUT" \
          -B "$OUT" -B "$(dirname "$PROT")" -B "$(dirname "$GFF")" -B "$(dirname "$DB")" "$IMG" "$@"; }
run omamer search --db "$DB" --query "$PROT" --out proteins.omamer --nthreads "$T"
# OMArk's isoform file: one line per gene, transcript IDs separated by ';' (as the pipeline)
awk -F'\t' '!/^#/ && $3=="mRNA" {
    match($9, /ID=[^;]+/);     id = substr($9, RSTART+3, RLENGTH-3)
    match($9, /Parent=[^;]+/); g  = substr($9, RSTART+7, RLENGTH-7)
    iso[g] = iso[g] (iso[g] ? ";" : "") id
} END { for (g in iso) print iso[g] }' "$GFF" > proteins.splice
mkdir -p omark_out
run omark -f proteins.omamer -d "$DB" -i proteins.splice -o omark_out
ls omark_out/*.sum
