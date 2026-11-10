#!/usr/bin/env bash
# funannotate_gene_ids.sh <models.gff3> > out.gff3
# Reshape gene IDs so `funannotate annotate` keeps its annotations. It keeps only annotation lines that start
# with a "base" taken from each protein header "<transcript> <gene>" (library.getGeneBasename, 1.8.17):
#   "_" in gene -> text up to the first "_" of the whole header + "_"; "-" in gene -> up to the first "-"; else the gene.
#  * NCBI style (transcript "rna-XM_1", has a "-"): gene -> "gene-<id>" with "_" -> ".", so the base is "rna".
#  * Other styles (transcript "g1.t2"): gene -> the transcript ID without ".tN" ("g1"), so the base is a prefix of
#    every transcript. Locus tags like ADA42R_000001 are replaced (their "_" would span the header and match nothing).
# Transcript IDs are left alone, so proteins and every search run on them stay valid.
set -euo pipefail
gff=$1
awk -F'\t' -v OFS='\t' '
function fix(id) { gsub(/_/, ".", id); if (id !~ /-/) id = "gene-" id; return id }
function attr(key,   n, a, i) {
    n = split($9, a, ";")
    for (i = 1; i <= n; i++) if (index(a[i], key "=") == 1) return substr(a[i], length(key) + 2)
    return ""
}
/^#/ { if (NR != FNR) print; next }
$3 != "gene" && $3 != "mRNA" && $3 != "tRNA" && $3 != "rRNA" { if (NR != FNR) print; next }
NR == FNR {
    if ($3 == "gene") next
    tid = attr("ID"); p = attr("Parent")
    if (tid !~ /-/) {
        g = tid; sub(/\.t[0-9]+$/, "", g)
        if (g == tid || g ~ /-/ || index(tid, g) != 1) { print "cannot derive a gene ID for transcript " tid > "/dev/stderr"; exit 1 }
        newgene[p] = g
    }
    next
}
$3 == "gene" { id = attr("ID"); $9 = "ID=" ((id in newgene) ? newgene[id] : fix(id)); print; next }
{
    p = attr("Parent"); np = (p in newgene) ? newgene[p] : fix(p)
    n = split($9, a, ";"); $9 = a[1]
    for (i = 2; i <= n; i++) $9 = $9 ";" ((index(a[i], "Parent=") == 1) ? "Parent=" np : a[i])
    print
}
' "$gff" "$gff"
