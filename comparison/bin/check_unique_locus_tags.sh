#!/usr/bin/env bash
# check_unique_locus_tags.sh <gbk>...
# funannotate compare pools GO terms, Pfam, InterPro, MEROPS, CAZy and the ortholog groups by locus tag / gene
# ID across all genomes, so a tag used by two genomes silently mixes their annotations. Stop if any is.
set -euo pipefail
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
for g in "$@"; do
    grep -h -o '/locus_tag="[^"]*"' "$g" | sort -u | sed "s|\$|\t$(basename "$g")|" >> "$tmp"
done
dups=$(cut -f1 "$tmp" | sort | uniq -d)
if [ -n "$dups" ]; then
    echo "ERROR: locus tags used by more than one genome ($(echo "$dups" | wc -l) tags), e.g.:" >&2
    { echo "$dups" | head -5 || true; } | while read -r t; do echo "  $t  in: $(grep -F "$t	" "$tmp" | cut -f2 | tr '\n' ' ')" >&2; done
    exit 1
fi
echo "locus tags unique across $# genomes ($(wc -l < "$tmp") tags)"
