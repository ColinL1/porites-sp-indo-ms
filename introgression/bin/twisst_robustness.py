#!/usr/bin/env python3
"""
Twisst robustness grid for one outgroup: window size W x tree method x target representation.
  W:       50, 100, 200, 500 SNPs
  method:  NJ (exact distances, every tied topology kept: common.nj_all) and ML (FastTree -nt -gtr)
  target:  rand2   two pseudo-haplotypes, heterozygous sites resolved at random (seed 1)
           phased2 the two WhatsHap haplotypes; windows only inside phase blocks >= 30 kb and
                   with <= 50% unphased heterozygous sites
           cons1   one sequence, heterozygous sites as IUPAC codes
Writes <mode>_W<W>_<NJ|ML>_trees.nwk / _windows.tsv and <mode>_groups.tsv; twisst.py is run
on each afterwards.
"""
import argparse, collections, gzip, random, subprocess, sys
from multiprocessing import Pool
import pysam
from common import read_references, by_role, nj_all, p_distances

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--geno", required=True)
ap.add_argument("--phased_vcf", required=True)
ap.add_argument("--block_list", required=True, help="whatshap stats --block-list output")
ap.add_argument("--references", required=True)
ap.add_argument("--target_name", required=True)
ap.add_argument("--outgroup", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--threads", type=int, default=32)
a = ap.parse_args()

WSIZES = [50, 100, 200, 500]
MINBLK = 30000
IUP = {"AG": "R", "GA": "R", "CT": "Y", "TC": "Y", "AC": "M", "CA": "M", "GT": "K", "TG": "K",
       "AT": "W", "TA": "W", "CG": "S", "GC": "S", "AA": "A", "CC": "C", "GG": "G", "TT": "T"}
refs = read_references(a.references)
P1, P2, T, OG = by_role(refs, "parent1")[0], by_role(refs, "parent2")[0], a.target_name, a.outgroup
TA, TB = f"{T}_a", f"{T}_b"

hdr = gzip.open(a.geno, "rt").readline().split()
ci = {n: i for i, n in enumerate(hdr)}
rows = collections.defaultdict(list)
with gzip.open(a.geno, "rt") as fh:
    next(fh)
    for line in fh:
        f = line.split()
        d = {tx: f[ci[tx]].split("/")[0] for tx in (P1, P2, OG)}
        if "N" in d.values():
            continue
        g = f[ci[T]].split("/")
        if g[0] == "N":
            continue
        d["target"] = (g[0], g[1])
        rows[f[0]].append((int(f[1]), d))
for c in rows:
    rows[c].sort(key=lambda x: x[0])

phased = {}
for r in pysam.VariantFile(a.phased_vcf):
    s = r.samples[0]
    gt = s["GT"]
    if gt is None or None in gt or not s.phased:
        continue
    al = [r.ref.upper()] + [x.upper() for x in r.alts]
    phased[(r.chrom, r.pos)] = (al[gt[0]], al[gt[1]])
blocks = collections.defaultdict(list)
for line in open(a.block_list):
    if line.startswith("#"):
        continue
    _, ch, ps, fr, to, nv = line.split()
    blocks[ch].append((int(fr), int(to)))


def in_long_block(ch, p0, p1):
    return any(s <= p0 and p1 <= e and e - s >= MINBLK for s, e in blocks.get(ch, []))


def nj(seq, taxa):
    """all tied NJ topologies of one window (see common.nj_all)"""
    return nj_all(taxa, p_distances(seq, taxa))


def ml(args):
    seq, taxa = args
    fa = "".join(f">{t}\n{seq[t]}\n" for t in taxa)
    return subprocess.run(["FastTree", "-nt", "-gtr", "-quiet", "-nopr"], input=fa,
                          capture_output=True, text=True).stdout.strip()


def build(mode, W):
    taxa = [P1, P2, T, OG] if mode == "cons1" else [P1, P2, TA, TB, OG]
    aligns, meta = [], []
    for ch, rr in rows.items():
        for i in range(0, len(rr) - W + 1, W):
            seg = rr[i:i + W]
            p0, p1 = seg[0][0], seg[-1][0]
            if mode == "phased2" and not in_long_block(ch, p0, p1):
                continue
            s = {t: [] for t in taxa}
            nN = 0
            for p, d in seg:
                for t in (P1, P2, OG):
                    s[t].append(d[t])
                g = d["target"]
                if mode == "cons1":
                    s[T].append(IUP[g[0] + g[1]])
                elif mode == "rand2":
                    gg = [g[1], g[0]] if (g[0] != g[1] and random.random() < 0.5) else [g[0], g[1]]
                    s[TA].append(gg[0]); s[TB].append(gg[1])
                elif g[0] == g[1]:
                    s[TA].append(g[0]); s[TB].append(g[1])
                elif (ch, p) in phased:
                    ph = phased[(ch, p)]
                    s[TA].append(ph[0]); s[TB].append(ph[1])
                else:
                    s[TA].append("N"); s[TB].append("N"); nN += 1
            if mode == "phased2" and nN > W * 0.5:
                continue
            aligns.append({t: "".join(v) for t, v in s.items()})
            meta.append((ch, p0, p1))
    return taxa, aligns, meta


random.seed(1)
for mode in ("rand2", "phased2", "cons1"):
    for W in WSIZES:
        taxa, aligns, meta = build(mode, W)
        if not aligns:
            print(f"{mode} W={W}: no windows", file=sys.stderr)
            continue
        with Pool(a.threads) as pool:
            njt = pool.starmap(nj, [(x, taxa) for x in aligns], chunksize=64)
            mlt = pool.map(ml, [(x, taxa) for x in aligns], chunksize=16)
        for lab, tr in (("NJ", njt), ("ML", [[t or "();"] for t in mlt])):
            pref = f"{a.outdir}/{mode}_W{W}_{lab}"
            with open(pref + "_trees.nwk", "w") as o, open(pref + "_tree_window.txt", "w") as tw:
                for k, ts in enumerate(tr):
                    for t in ts:
                        o.write(t + "\n"); tw.write(f"{k}\n")
            with open(pref + "_windows.tsv", "w") as o:
                for c, x, y in meta:
                    o.write(f"{c}\t{x}\t{y}\n")
        with open(f"{a.outdir}/{mode}_groups.tsv", "w") as g:
            tips = [T] if mode == "cons1" else [TA, TB]
            g.write(f"{P1}\t{P1}\n{P2}\t{P2}\n" + "".join(f"{x}\t{T}\n" for x in tips) + f"{OG}\t{OG}\n")
        print(f"built {mode} W={W}: {len(aligns)} windows", file=sys.stderr)
