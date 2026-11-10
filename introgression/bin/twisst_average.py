#!/usr/bin/env python3
"""Average twisst weights over the trees of each window (tied NJ topologies get equal weight).
Each tree's weights are first normalised to sum 1, so every window counts once.
usage: twisst_average.py <per-tree weights.csv.gz> <tree_window.txt> <out weights.csv.gz>"""
import gzip, sys
import numpy as np

src, idx, out = sys.argv[1:4]
head, rows = [], []
with gzip.open(src, "rt") as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("topo"):
            head.append(line)
        elif line.strip():
            rows.append([float(x) for x in line.split()])
win = np.loadtxt(idx, dtype=int, ndmin=1)
W = np.array(rows, float)
if len(W) != len(win):
    raise SystemExit(f"{src}: {len(W)} weight rows but {len(win)} trees")
s = W.sum(1, keepdims=True)
W = np.divide(W, s, out=np.zeros_like(W), where=s > 0)
n = win.max() + 1
acc = np.zeros((n, W.shape[1])); cnt = np.zeros(n)
np.add.at(acc, win, W); np.add.at(cnt, win, (s[:, 0] > 0).astype(float))
with gzip.open(out, "wt") as o:
    o.writelines(head)
    for a, c in zip(acc, cnt):
        o.write("\t".join(f"{x / c:.6g}" if c else "0" for x in a) + "\n")
