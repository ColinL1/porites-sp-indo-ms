#!/usr/bin/env python3
"""
Does the per-window topology weight follow assembly context rather than biology?

For the rand2 W=200 NJ windows of the robustness grid, per window:
  depth      mean ONT depth (mosdepth 5 kb bins) over the window
  depth_z    |depth - genome-wide median| / median (collapsed repeats give high depth)
  n_joins    RagTag joins (AGP 'U' gaps) inside the window       (only with --agp)
  join_dist  bp from the window midpoint to the nearest join      (only with --agp)
  p1_call    length of the parent1 callable segment at the window midpoint
  p2_call    same for parent2
Correlates each with T1/T2/T3 and compares clean windows (depth within 25% of the median, no
join, both callable segments >= 20 kb) with the rest. With v1.2 (scaffolded on parent1) the
join columns test whether parent1-guided joins push windows towards T1.
Outputs: context.tsv, context_summary.txt
"""
import argparse, collections, csv, gzip
import numpy as np
from common import read_references, by_role, twisst_weights

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", required=True, help="07_twisst_robustness")
ap.add_argument("--depth", required=True, help="mosdepth <prefix>.regions.bed.gz")
ap.add_argument("--aln_dir", required=True)
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--agp", default="")
a = ap.parse_args()
refs = read_references(a.references)
P1, P2, T = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0], a.target_name

Wt = twisst_weights(f"{a.dir}/rand2_W200_NJ_weights.csv.gz", P1, P2, T)
wins = [(c, int(x), int(y)) for c, x, y in (l.split() for l in open(f"{a.dir}/rand2_W200_NJ_windows.tsv"))]

dep = collections.defaultdict(list)
with gzip.open(a.depth, "rt") as fh:
    for line in fh:
        c, x, y, d = line.split()
        dep[c].append((int(x), int(y), float(d)))
alld = np.array([d for v in dep.values() for _, _, d in v])
med = np.median(alld[alld > 0])

joins = collections.defaultdict(list)
if a.agp:
    for line in open(a.agp):
        f = line.rstrip("\n").split("\t")
        if not line.startswith("#") and len(f) > 4 and f[4] == "U":
            joins[f[0]].append((int(f[1]) + int(f[2])) // 2)


def load_segments(p):
    d = collections.defaultdict(list)
    for line in open(p):
        c, s, e = line.split()[:3]
        d[c].append((int(s), int(e)))
    return d
c1 = load_segments(f"{a.aln_dir}/{P1}.callable.bed")
c2 = load_segments(f"{a.aln_dir}/{P2}.callable.bed")


def seg_len(d, c, pos):
    return next((e - s for s, e in d.get(c, []) if s <= pos < e), 0)


rowd = []
for (c, x, y), tw in zip(wins, Wt):
    if not np.isfinite(tw).all():
        continue
    xs = [d for s, e, d in dep.get(c, []) if e > x and s < y]
    dp = np.mean(xs) if xs else np.nan
    mid = (x + y) // 2
    js = joins.get(c, [])
    r = dict(chrom=c, start=x, end=y, T1=tw[0], T2=tw[1], T3=tw[2], depth=dp,
             depth_z=abs(dp - med) / med if dp == dp else np.nan)
    if a.agp:
        r["n_joins"] = sum(x <= j <= y for j in js)
        r["join_dist"] = min(abs(j - mid) for j in js) if js else 1e9
    r["p1_call"] = seg_len(c1, c, mid)
    r["p2_call"] = seg_len(c2, c, mid)
    rowd.append(r)
with open(f"{a.dir}/context.tsv", "w") as o:
    w = csv.DictWriter(o, fieldnames=list(rowd[0]), delimiter="\t")
    w.writeheader(); w.writerows(rowd)

A = {k: np.array([r[k] for r in rowd], float) for k in rowd[0] if k != "chrom"}
ok = np.isfinite(A["depth"])


def corr(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    return np.corrcoef(x[m], y[m])[0, 1] if m.sum() >= 20 else np.nan


ctx = ["depth", "depth_z"] + (["n_joins", "join_dist"] if a.agp else []) + ["p1_call", "p2_call"]
clean = ok & (A["depth_z"] < 0.25) & (A["p1_call"] >= 20000) & (A["p2_call"] >= 20000)
if a.agp:
    clean &= A["n_joins"] == 0
messy = ok & ~clean
with open(f"{a.dir}/context_summary.txt", "w") as o:
    o.write(f"parent1 = {P1}; parent2 = {P2}; AGP: {a.agp or 'none (join columns skipped)'}\n")
    o.write(f"windows = {ok.sum()}; genome-wide median ONT depth = {med:.1f}x\n")
    o.write(f"mean weights: T1={A['T1'][ok].mean():.3f} T2={A['T2'][ok].mean():.3f} T3={A['T3'][ok].mean():.3f}\n\n")
    o.write("Pearson r (context vs topology weight)\ncontext\tr_T1\tr_T2\tr_T3\n")
    for k in ctx:
        o.write(f"{k}\t{corr(A[k], A['T1']):+.3f}\t{corr(A[k], A['T2']):+.3f}\t{corr(A[k], A['T3']):+.3f}\n")
    for lab, m in (("clean", clean), ("messy", messy)):
        o.write(f"\n{lab} windows: {m.sum()}  T1={A['T1'][m].mean():.3f}  T2={A['T2'][m].mean():.3f}  "
                f"T3={A['T3'][m].mean():.3f}  T1-T2={A['T1'][m].mean() - A['T2'][m].mean():+.3f}")
    o.write("\n")
print(open(f"{a.dir}/context_summary.txt").read())
