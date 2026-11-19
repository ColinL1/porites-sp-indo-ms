#!/usr/bin/env bash
# blob_supplement.sh: supplementary figure, BlobToolKit blob plot of the assembly before filtering, as the
# pipeline renders it (`blobtk plot`, genomehubs/blobtoolkit:4.4.4; phyla of the best BLAST hits, 7 largest + other), printed at
# full page width (7.5 in) so that BlobToolKit's own text is legible (legend ~6.4 pt, ticks ~7.9 pt).
# usage: blob_supplement.sh blobdir out_prefix      (needs singularity and rsvg-convert: conda env r_env_DESeq2)
set -euo pipefail
BLOBDIR=$(readlink -f "$1"); OUT=$2
IMG=${BTK_IMG:-$HOME/SINGULARITY_CACHE/genomehubs-blobtoolkit-4.4.4.img}
OUTDIR=$(readlink -f "$(dirname "$OUT")"); NAME=$(basename "$OUT")
singularity exec --containall -B "$BLOBDIR" -B "$OUTDIR" "$IMG" \
  blobtk plot -v blob -d "$BLOBDIR" -c bestsumorder_phylum --cat-count 8 -o "$OUTDIR/$NAME.svg"
W_IN=7.5
rsvg-convert -f pdf -a -w "$(python3 -c "print(round($W_IN * 96))")" -o "$OUT.pdf" "$OUT.svg"
rsvg-convert -f png -a -w "$(python3 -c "print(round($W_IN * 300))")" -o "$OUT.png" "$OUT.svg"
