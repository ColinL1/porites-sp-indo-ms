#!/usr/bin/env python3
"""Mean T1/T2/T3 weights for every cell of the Twisst robustness grid -> robustness_summary.tsv"""
import argparse, glob, os, re
import numpy as np
from common import read_references, by_role, twisst_weights

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--dir", required=True)
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
a = ap.parse_args()
refs = read_references(a.references)
P1, P2, T = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0], a.target_name

rows = []
for f in glob.glob(f"{a.dir}/*_W*_*_weights.csv.gz"):
    mode, W, method = re.match(r"(\w+?)_W(\d+)_(NJ|ML)_weights", os.path.basename(f)).groups()
    w = twisst_weights(f, P1, P2, T)
    w = w[np.isfinite(w).all(1)]
    rows.append((mode, int(W), method, len(w), *w.mean(0)))
with open(f"{a.dir}/robustness_summary.tsv", "w") as o:
    o.write("representation\twindow_snps\tmethod\tn_windows\tT1\tT2\tT3\tplurality\n")
    for r in sorted(rows):
        o.write("\t".join(map(str, r[:4])) + "".join(f"\t{x:.4f}" for x in r[4:])
                + f"\t{('T1', 'T2', 'T3')[int(np.argmax(r[4:]))]}\n")
print(open(f"{a.dir}/robustness_summary.tsv").read())
