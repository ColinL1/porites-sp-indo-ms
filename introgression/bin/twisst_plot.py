#!/usr/bin/env python3
"""
Twisst summary per outgroup: mean topology weights (twisst_summary.tsv), weights along the
largest scaffolds (twisst_painting.png) and a genome-wide bar chart (twisst_genomewide_bar.png).
T1 = (target, parent1), T2 = (target, parent2), T3 = (parent1, parent2), the species tree.
"""
import argparse, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import read_references, by_role, twisst_weights

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--dir", required=True, help="06_twisst: <outgroup>_weights.csv.gz, <outgroup>_windows.tsv")
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--fai", required=True)
a = ap.parse_args()

refs = read_references(a.references)
P1, P2, T = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0], a.target_name
OGS = by_role(refs, "outgroup")
fai = dict((l.split()[0], int(l.split()[1])) for l in open(a.fai))

data = {}
for og in OGS:
    W = twisst_weights(f"{a.dir}/{og}_weights.csv.gz", P1, P2, T)
    wins = [l.split("\t")[:2] for l in open(f"{a.dir}/{og}_windows.tsv").readlines()[1:]]
    data[og] = [(c, int(s), *w) for (c, s), w in zip(wins, W) if np.isfinite(w).all()]

with open(f"{a.dir}/twisst_summary.tsv", "w") as o:
    o.write(f"outgroup\tn_windows\tT1_{T}_{P1}\tT2_{T}_{P2}\tT3_{P1}_{P2}\tplurality\n")
    for og in OGS:
        m = np.array([r[2:] for r in data[og]]).mean(0)
        o.write(f"{og}\t{len(data[og])}\t{m[0]:.4f}\t{m[1]:.4f}\t{m[2]:.4f}\t{('T1', 'T2', 'T3')[int(m.argmax())]}\n")
print(open(f"{a.dir}/twisst_summary.tsv").read())

present = {r[0] for og in OGS for r in data[og]}
big = [c for c in sorted(fai, key=lambda c: -fai[c]) if fai[c] >= 3_000_000 and c in present][:30]
fig, axes = plt.subplots(len(OGS), 1, figsize=(14, 3.2 * len(OGS)), squeeze=False)
for ax, og in zip(axes[:, 0], OGS):
    by_chrom = collections.defaultdict(list)
    for c, s, t1, t2, t3 in data[og]:
        by_chrom[c].append((s, t1, t2))
    offset, xt, xl = 0, [], []
    for c in big:
        d = np.array(sorted(by_chrom.get(c, [])), float)
        if len(d):
            ax.plot(d[:, 0] / 1e6 + offset, d[:, 1], ".", ms=1.5, color="#1f77b4", alpha=.5, label=f"T1 ({T},{P1})" if offset == 0 else None)
            ax.plot(d[:, 0] / 1e6 + offset, d[:, 2], ".", ms=1.5, color="#d62728", alpha=.5, label=f"T2 ({T},{P2})" if offset == 0 else None)
        xt.append(offset + fai[c] / 2e6); xl.append(c.replace("_RagTag", ""))
        offset += fai[c] / 1e6 + 2
    ax.axhline(1 / 3, color="grey", lw=.6, ls=":")
    ax.set_ylim(-.02, 1.02)
    ax.set_ylabel(f"outgroup {og}\ntopology weight", fontsize=8)
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5, rotation=90)
    ax.legend(fontsize=7, loc="upper right", markerscale=4)
axes[-1, 0].set_xlabel("concatenated scaffold position (Mb)")
fig.tight_layout(); fig.savefig(f"{a.dir}/twisst_painting.png", dpi=130)

fig, ax = plt.subplots(figsize=(6, 4))
x, w = np.arange(len(OGS)), .25
means = [np.array([r[2:] for r in data[og]]).mean(0) for og in OGS]
ax.bar(x - w, [m[2] for m in means], w, label=f"T3 ({P1},{P2})", color="#2ca02c")
ax.bar(x, [m[0] for m in means], w, label=f"T1 ({T},{P1})", color="#1f77b4")
ax.bar(x + w, [m[1] for m in means], w, label=f"T2 ({T},{P2})", color="#d62728")
ax.set_xticks(x); ax.set_xticklabels(OGS, fontsize=8)
ax.axhline(1 / 3, color="grey", lw=.6, ls=":")
ax.set_ylabel("mean topology weight"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{a.dir}/twisst_genomewide_bar.png", dpi=130)
