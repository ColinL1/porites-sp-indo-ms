#!/usr/bin/env python3
"""
Neighbour-joining trees of non-overlapping W-SNP windows for {parent1, parent2, target_a,
target_b, outgroup}, from the quartet geno file. The target is split into two
pseudo-haplotypes, heterozygous sites resolved at random (seed 1). Distances are exact
fractions and NJ follows every tie (common.nj_all): a window with tied joins gets one tree per
tied topology, and twisst_average.py later gives them equal weight. Writes <prefix>_trees.nwk
(one or more trees per window), <prefix>_tree_window.txt (window number of each tree),
<prefix>_windows.tsv and <prefix>_groups.tsv.
"""
import argparse, collections, gzip, random
from common import read_references, by_role, nj_all, p_distances

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--geno", required=True)
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--outgroup", required=True)
ap.add_argument("--window", type=int, default=50)
ap.add_argument("--prefix", required=True)
a = ap.parse_args()
random.seed(1)

refs = read_references(a.references)
P1, P2, T, OG, W = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0], a.target_name, a.outgroup, a.window
hdr = gzip.open(a.geno, "rt").readline().split()
idx = {n: hdr.index(n) for n in (P1, P2, T, OG)}
taxa = [P1, P2, f"{T}_a", f"{T}_b", OG]

rows_by_chrom = collections.defaultdict(list)
with gzip.open(a.geno, "rt") as fh:
    next(fh)
    for line in fh:
        f = line.split()
        row = {tx: f[idx[tx]].split("/")[0] for tx in (P1, P2, OG)}
        if "N" in row.values():
            continue
        g = f[idx[T]].split("/")
        if g[0] == "N":
            continue
        if g[0] != g[1] and random.random() < 0.5:
            g = [g[1], g[0]]
        row[f"{T}_a"], row[f"{T}_b"] = g[0], g[1]
        rows_by_chrom[f[0]].append((int(f[1]), row))

nwin = ntree = 0
with open(f"{a.prefix}_trees.nwk", "w") as trees, open(f"{a.prefix}_windows.tsv", "w") as meta, \
        open(f"{a.prefix}_tree_window.txt", "w") as tw:
    meta.write("chrom\tstart\tend\tn\n")
    for chrom, rows in rows_by_chrom.items():
        rows.sort(key=lambda x: x[0])
        for i in range(0, len(rows) - W + 1, W):
            block = [r for _, r in rows[i:i + W]]
            seqs = {t: "".join(r[t] for r in block) for t in taxa}
            for nwk in nj_all(taxa, p_distances(seqs, taxa)):
                trees.write(nwk + "\n"); tw.write(f"{nwin}\n"); ntree += 1
            meta.write(f"{chrom}\t{rows[i][0]}\t{rows[i+W-1][0]}\t{W}\n")
            nwin += 1
with open(f"{a.prefix}_groups.tsv", "w") as g:
    g.write(f"{P1}\t{P1}\n{P2}\t{P2}\n{T}_a\t{T}\n{T}_b\t{T}\n{OG}\t{OG}\n")
print(f"{OG}: {nwin} windows, {ntree} trees ({ntree - nwin} extra from tied joins)")
