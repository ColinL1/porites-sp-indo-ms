#!/usr/bin/env bash
# btk_view.sh: a BlobToolKit plot as BlobToolKit draws it, at print size: `blobtools view` (the BlobToolKit
# viewer, as BlobToolKit's own pipeline makes its static images) with large fonts, in the pipeline's container
# (genomehubs/blobtoolkit:4.4.4); the "Dataset: <blobdir name>" label is removed (source text) and the SVG is
# rasterised by the container's Firefox (it draws the snail's ring labels, which rsvg cannot).
# usage: btk_view.sh blobdir view out_prefix px [--param key=value ...]
#        view: blob | cumulative | snail; px: width of the PNG (e.g. 3.75 in at 600 dpi = 2250)
# writes out_prefix.svg (viewer output, label removed) and out_prefix.png
set -euo pipefail
BLOBDIR=$(readlink -f "$1"); VIEW=$2; OUT=$3; PX=$4; shift 4
IMG=${BTK_IMG:-$HOME/SINGULARITY_CACHE/genomehubs-blobtoolkit-4.4.4.img}
OUTDIR=$(readlink -f "$(dirname "$OUT")"); NAME=$(basename "$OUT")
WORK=$(mktemp -d -p "$OUTDIR" .btk_view.XXXX); trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/home"
PARAMS=(--param largeFonts=true); for p in "$@"; do PARAMS+=("$p"); done
run() { singularity exec --containall --home "$WORK/home:/home/btk" -B "$BLOBDIR" -B "$WORK" "$IMG" "$@"; }
# the viewer's Firefox session can hang without returning; stop it after 5 min and try once more
view() { timeout --kill-after=10 300 singularity exec --containall --home "$WORK/home:/home/btk" -B "$BLOBDIR" \
           -B "$WORK" "$IMG" blobtools view --driver firefox --view "$VIEW" "${PARAMS[@]}" --format svg \
           --out "$WORK" --timeout 240 "$BLOBDIR" > "$WORK/view.log" 2>&1; }
view || { echo "blobtools view failed or timed out; retrying" >&2; rm -f "$WORK"/*.svg; view; } || \
  { cat "$WORK/view.log"; exit 1; }
SVG=$(ls "$WORK"/*.svg | head -1); [ -n "$SVG" ] || { cat "$WORK/view.log"; exit 1; }
python3 - "$SVG" "$OUTDIR/$NAME.svg" <<'EOF'
import re, sys
s = open(sys.argv[1]).read()
s, n = re.subn(r'<text[^>]*>\s*Dataset:[^<]*</text>', '', s)
if n > 1:
    sys.exit(f'more than one "Dataset:" label in {sys.argv[1]} ({n})')
open(sys.argv[2], 'w').write(s)
EOF
# page that draws the SVG at PX wide (height from the viewBox), screenshot by headless Firefox
read -r W H < <(python3 -c "import re,sys; v=re.search(r'viewBox=\"([-0-9. ]+)\"', open(sys.argv[1]).read()).group(1).split(); print(v[2], v[3])" "$OUTDIR/$NAME.svg")
PH=$(python3 -c "print(round($PX * $H / $W))")
cp "$OUTDIR/$NAME.svg" "$WORK/plot.svg"
echo "<html><body style='margin:0;background:#fff'><img src='plot.svg' style='width:${PX}px;height:${PH}px;display:block'></body></html>" > "$WORK/page.html"
run firefox --headless --screenshot "$WORK/shot.png" --window-size="$PX,$PH" "file://$WORK/page.html" > /dev/null 2>&1
mv "$WORK/shot.png" "$OUTDIR/$NAME.png"
